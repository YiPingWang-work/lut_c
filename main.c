#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include "lut/mul_mat_with_lut.h"
#include <time.h>
#include <sys/stat.h>


int M = 512;
int K = 8192;


static inline uint64_t now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ((uint64_t)ts.tv_sec * 1000000000ull + ts.tv_nsec)/1000;
}

static inline float rand_float() {
    return 0.001f + (1.0f - 0.001f) * ((float)random() / (float)RAND_MAX);
}

int load_w_bit(const char *path, block_ifairy *w) {
    FILE *fp = fopen(path, "r");
    if (!fp) {
        perror("fopen");
        return 2;
    }

    int block_idx = 0;
    int bit_cnt = 0;     // 已读 bit 数（0~511）
    int elem_idx = 0;    // 已写 element 数（0~255）
    uint8_t elem = 0;    // 当前 2bit element
    int elem_bits = 0;

    int c;
    while ((c = fgetc(fp)) != EOF && block_idx < (M*K/QK_K)) {
        if (c == '\n' || c == '\r') continue;
        if (c != '0' && c != '1')  continue;

        elem = (elem << 1) | (c - '0');
        elem_bits++;
        bit_cnt++;

        // 收满 2 bit → 写一个 element
        if (elem_bits == 2) {
            int byte_idx = elem_idx >> 2;                 // /4
            int shift    = 6 - ((elem_idx & 3) << 1);     // 高位先放

            w[block_idx].qs[byte_idx] |= (elem & 0x3) << shift;

            elem = 0;
            elem_bits = 0;
            elem_idx++;
        }
        w[block_idx].d_real = rand_float();
        w[block_idx].d_imag = rand_float();

        // 一个 block 完成：512 bit
        if (bit_cnt == QK_K * 2) {
            // printf("block %d: d_real=%f, d_imag=%f\n",
            //        block_idx, w[block_idx].d_real, w[block_idx].d_imag);
            block_idx++;
            bit_cnt = 0;
            elem_idx = 0;
            elem = 0;
            elem_bits = 0;
        }
    }

    fclose(fp);

    if (block_idx != (M*K/QK_K)) {
        fprintf(stderr, "warning: loaded %d / %d blocks\n",
                block_idx, (M*K/QK_K));
    }

    return 0;
}


int load_act_float(const char *path, float *vec, int num_values) {
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    int count = 0;
    float tmp;
    while (count < num_values && fscanf(f, "%f", &tmp) == 1) {
        vec[count++] = tmp;
    }
    fclose(f);
    return (count == num_values) ? 0 : -2;
}


static void print_u8_binary(uint8_t v) {
    for (int i = 7; i >= 0; --i) putchar((v & (1u << i)) ? '1' : '0');
    putchar('\n');
}


/**
 * 创建输出文件夹
 * @param dirname 文件夹名称
 */
static int create_output_dir(const char *dirname) {
    struct stat st = {0};
    if (stat(dirname, &st) == -1) {
        if (mkdir(dirname, 0700) == -1) {
            perror("mkdir");
            return -1;
        }
    }
    return 0;
}

/**
 * 将float数组保存到文件
 * @param folder 文件夹名称
 * @param filename 输出文件名
 * @param data 数据数组
 * @param count 元素个数
 */
static int save_float_array_to_file(const char *folder, const char *filename, const float *data, int count) {
    char filepath[256];
    snprintf(filepath, sizeof(filepath), "%s/%s", folder, filename);
    
    FILE *fp = fopen(filepath, "w");
    if (!fp) {
        perror("fopen");
        return -1;
    }
    
    for (int i = 0; i < count; i++) {
        fprintf(fp, "%.9f\n", data[i]);
    }
    
    fclose(fp);
    return 0;
}


// 平均相对误差: (1/k) * sum |dst - dst_ref| / |dst_ref|
// 这里 k 使用当前输出长度 M*2；当分母为 0 时用 1e-9 防止除零。
static float mean_relative_error_vs_ref(const float *dst, const float *dst_ref) {
    const int count = M * 2;
    double sum = 0.0;
    for (int i = 0; i < count; i++) {
        const float denom = fabsf(dst_ref[i]);
        sum += fabsf(dst[i] - dst_ref[i]) / (denom + 1e-9f);
    }
    return (float)(sum / (double)count);
}


// 最大相对误差: max |dst - dst_ref| / |dst_ref|
// 当分母为 0 时用 1e-9 防止除零。
static float max_relative_error_vs_ref(const float *dst, const float *dst_ref) {
    const int count = M * 2;
    float max_err = 0.0f;
    for (int i = 0; i < count; i++) {
        const float denom = fabsf(dst_ref[i]);
        const float rel = fabsf(dst[i] - dst_ref[i]) / (denom + 1e-9f);
        if (rel > max_err) {
            max_err = rel;
        }
    }
    return max_err;
}


// 余弦相似度: dot(dst, dst_ref) / (||dst|| * ||dst_ref||)
// 当分母接近 0 时返回 0。
static float cosine_similarity_vs_ref(const float *dst, const float *dst_ref) {
    const int count = M * 2;
    double dot = 0.0;
    double norm_dst = 0.0;
    double norm_ref = 0.0;
    for (int i = 0; i < count; i++) {
        const double a = (double)dst[i];
        const double b = (double)dst_ref[i];
        dot += a * b;
        norm_dst += a * a;
        norm_ref += b * b;
    }
    const double denom = sqrt(norm_dst) * sqrt(norm_ref);
    if (denom < 1e-12) {
        return 0.0f;
    }
    return (float)(dot / denom);
}


// 信噪比 SNR(dB): 10 * log10( sum(dst_ref^2) / sum((dst - dst_ref)^2) )
// 当噪声功率接近 0 时返回一个较大值 120 dB。
static float snr_db_vs_ref(const float *dst, const float *dst_ref) {
    const int count = M * 2;
    double signal_power = 0.0;
    double noise_power = 0.0;
    for (int i = 0; i < count; i++) {
        const double ref = (double)dst_ref[i];
        const double err = (double)dst[i] - ref;
        signal_power += ref * ref;
        noise_power += err * err;
    }
    if (noise_power < 1e-24) {
        return 1e10f; // 一致
    }
    return (float)(10.0 * log10(signal_power / noise_power));
}


int write_to_file_int8x16x2_t(const char *filename, int block, int begin, int end, int8x16x2_t v) {

    FILE *fp = fopen(filename, "a");   // 追加写
    if (!fp) {
        perror("fopen");
        return -1;
    }

    int8_t buf0[16];
    int8_t buf1[16];

    vst1q_s8(buf0, v.val[0]);
    vst1q_s8(buf1, v.val[1]);

    fprintf(fp, "%d:%d->%d real: ", block, begin, end);
    for (int i = 0; i < 16; i++) {
        fprintf(fp, "%4d", buf0[i]);
    }

    fprintf(fp, "\n%d:%d->%d imag: ", block, begin, end);
    for (int i = 0; i < 16; i++) {
        fprintf(fp, "%4d", buf1[i]);
    }
    fprintf(fp, "\n");

    fclose(fp);
    return 0;
}


void compare(const block_ifairy *w, const float *act) {
    float *dst1 = calloc(M*2, sizeof(float));
    float *dst2 = calloc(M*2, sizeof(float));
    float *dst3 = calloc(M*2, sizeof(float));
    float *dst4 = calloc(M*2, sizeof(float));
    float *dst5 = calloc(M*2, sizeof(float));

    for(int i = 0; i < M*2; i++) {
        dst1[i] = 0.0f;
        dst2[i] = 0.0f;
        dst3[i] = 0.0f;
        dst4[i] = 0.0f;
        dst5[i] = 0.0f;
    }

    // 验证程序
    long long t0 = now_ns();
    mul_mat_mxk_kx1(K, 0, M-1, w, act, dst1);
    long long t1 = now_ns();
    printf("全精度计算,  耗时: %lld us\n", (t1 - t0));
    

    block_ifairy_q16 *act_block = alloc_act_block_ifairy_q16(K);
    // 16路查表
    act_float_2_block_ifairy_q16(K, act, act_block, 42.6f);
    lut_block *lut = alloc_lut_q8(K);
    block_ifairy_1x3 *_w = alloc_w(M, K);
    t0 = now_ns();
    transpose(M, K, w, _w);
    t1 = now_ns();
    printf("转置权重, 耗时: %f us\n", (float)(t1 - t0)/K);
    t0 = now_ns();
    generate_lut_q8_block_ifairy_q16(K, act_block, lut);
    t1 = now_ns();
    mul_mat_mxk_kx1_with_lut_q8(K, 0, M-1, _w, lut, dst2);
    long long t2 = now_ns();
    printf("本文方法, 耗时: %lld us, 生成LUT耗时: %lld us\n", (t2 - t0), (t1 - t0));
    free_lut_q8(lut);
    

    // 1路查表
    act_float_2_block_ifairy_q16(K, act, act_block, 127.0f);
    int16_t *lut_v_old = alloc_lut_v_q16_old(K);
    float *lut_scale_old = alloc_lut_scale_old(K);
    block_ifairy_1x3_old *_w_old = alloc_w_old(M, K);
    transpose_old(M, K, w, _w_old);
    t0 = now_ns();
    generate_lut_q16_block_ifairy_q16_old(K, act_block, lut_v_old, lut_scale_old);
    t1 = now_ns();
    mul_mat_mxk_kx1_with_lut_q16_old(K, 0, M-1, _w_old, lut_v_old, lut_scale_old, dst3);
    t2 = now_ns();
    printf("传统查表, 耗时: %lld us, 生成LUT耗时: %lld us\n", (t2 - t1), (t1 - t0));

    
    // 传统矩阵乘法
    t0 = now_ns();
    mul_mat_mxk_kx1_q8(K, 0, M-1, w, act_block, dst4);
    t1 = now_ns();
    printf("传统矩阵乘法, 耗时: %lld us\n", (t1 - t0));


    // simd矩阵乘法
    t0 = now_ns();
    mul_mat_mxk_kx1_q8_simd(K, 0, M-1, w, act_block, dst5);
    t1 = now_ns();
    printf("simd矩阵乘法, 耗时: %lld us\n", (t1 - t0));


    
    free_lut_v_q16_old(lut_v_old);
    free_lut_scale_old(lut_scale_old);
    free_w_old(_w_old);
    free_act_block_ifairy_q16(act_block);


    int errors = 0;
    float max_mismatch = 0.0f;
    for (int i = 0; i < M*2; i++) {
        float mismatch = fabsf(dst1[i] - dst2[i])/(fabsf(dst1[i])+1e-9);
        if (mismatch > 0.1) {
            // printf("❌ %d ==> %f, %f, %f\n", i, fabsf(dst1[i] - dst2[i])/(fabsf(dst1[i])+1e-9), dst1[i], dst2[i]);
            if (mismatch > max_mismatch) {
                max_mismatch = mismatch;
            }
            errors++;
        } else {
            // printf("✅ %d ==> %f, %f, %f\n", i, fabsf(dst1[i] - dst2[i])/(fabsf(dst1[i])+1e-9), dst1[i], dst2[i]);
        }
    }
    printf("本文方法: 0.1 accuracy: %f, max_mismatch: %f\n", ((float)M*2 - (float)errors)/((float)M*2), max_mismatch);
    
    errors = 0;
    max_mismatch = 0.0f;
    for (int i = 0; i < M*2; i++) {
        float mismatch = fabsf(dst1[i] - dst3[i])/(fabsf(dst1[i])+1e-9);
        if (mismatch > 0.1) {
            // printf("❌ %d ==> %f, %f, %f\n", i, fabsf(dst1[i] - dst3[i])/(fabsf(dst1[i])+1e-9), dst1[i], dst3[i]);
            if (mismatch > max_mismatch) {
                max_mismatch = mismatch;
            }
            errors++;
        } else {
            // printf("✅ %d ==> %f, %f, %f\n", i, fabsf(dst1[i] - dst3[i])/(fabsf(dst1[i])+1e-9), dst1[i], dst3[i]);
        }
    }
    printf("传统查表: 0.1 accuracy: %f, max_mismatch: %f\n", ((float)M*2 - (float)errors)/((float)M*2), max_mismatch);

    errors = 0;
    max_mismatch = 0.0f;
    for (int i = 0; i < M*2; i++) {
        float mismatch = fabsf(dst1[i] - dst4[i])/(fabsf(dst1[i])+1e-9);
        if (mismatch > 0.1) {
            // printf("❌ %d ==> %f, %f, %f\n", i, fabsf(dst1[i] - dst4[i])/(fabsf(dst1[i])+1e-9), dst1[i], dst4[i]);
            if (mismatch > max_mismatch) {
                max_mismatch = mismatch;
            }
            errors++;
        } else {
            // printf("✅ %d ==> %f, %f, %f\n", i, fabsf(dst1[i] - dst4[i])/(fabsf(dst1[i])+1e-9), dst1[i], dst4[i]);
        }
    }
    printf("传统矩阵乘法: 0.1 accuracy: %f, max_mismatch: %f\n", ((float)M*2 - (float)errors)/((float)M*2), max_mismatch);


    errors = 0;
    max_mismatch = 0.0f;
    for (int i = 0; i < M*2; i++) {
        float mismatch = fabsf(dst1[i] - dst5[i])/(fabsf(dst1[i])+1e-9);
        if (mismatch > 0.1) {
            // printf("❌ %d ==> %f, %f, %f\n", i, fabsf(dst1[i] - dst5[i])/(fabsf(dst1[i])+1e-9), dst1[i], dst5[i]);
            if (mismatch > max_mismatch) {
                max_mismatch = mismatch;
            }
            errors++;
        } else {
            // printf("✅ %d ==> %f, %f, %f\n", i, fabsf(dst1[i] - dst5[i])/(fabsf(dst1[i])+1e-9), dst1[i], dst5[i]);
        }
    }
    printf("simd矩阵乘法: 0.1 accuracy: %f, max_mismatch: %f\n", ((float)M*2 - (float)errors)/((float)M*2), max_mismatch);



    float mre_lut16_vs_ref = mean_relative_error_vs_ref(dst2, dst1);
    float mre_lut1_vs_ref = mean_relative_error_vs_ref(dst3, dst1);
    float xre_lut16_vs_ref = max_relative_error_vs_ref(dst2, dst1);
    float xre_lut1_vs_ref = max_relative_error_vs_ref(dst3, dst1);
    float cos_lut16_vs_ref = cosine_similarity_vs_ref(dst2, dst1);
    float cos_lut1_vs_ref = cosine_similarity_vs_ref(dst3, dst1);
    float snr_lut16_vs_ref = snr_db_vs_ref(dst2, dst1);
    float snr_lut1_vs_ref = snr_db_vs_ref(dst3, dst1);
    printf("平均相对误差(本文方法 vs 基准): %f\n", mre_lut16_vs_ref);
    printf("平均相对误差(传统查表 vs 基准): %f\n", mre_lut1_vs_ref);
    printf("最大相对误差(本文方法 vs 基准): %f\n", xre_lut16_vs_ref);
    printf("最大相对误差(传统查表 vs 基准): %f\n", xre_lut1_vs_ref);
    printf("余弦相似度(本文方法 vs 基准): %f\n", cos_lut16_vs_ref);
    printf("余弦相似度(传统查表 vs 基准): %f\n", cos_lut1_vs_ref);
    printf("信噪比SNR(本文方法 vs 基准): %f dB\n", snr_lut16_vs_ref);
    printf("信噪比SNR(传统查表 vs 基准): %f dB\n", snr_lut1_vs_ref);


    // 保存结果到output文件夹
    create_output_dir("output");
    
    char filename_buf[256];
    snprintf(filename_buf, sizeof(filename_buf), "全精度计算_%d_%d.txt", M, K);
    save_float_array_to_file("output", filename_buf, dst1, M*2);
    
    snprintf(filename_buf, sizeof(filename_buf), "本文方法_%d_%d.txt", M, K);
    save_float_array_to_file("output", filename_buf, dst2, M*2);
    
    snprintf(filename_buf, sizeof(filename_buf), "传统查表_%d_%d.txt", M, K);
    save_float_array_to_file("output", filename_buf, dst3, M*2);


    free(dst1);
    free(dst2);
    free(dst3);
    free(dst4);
    free(dst5);
}


void benchmark_speedup_from_file(const char *matrix_file, const char *act_file) {
    const int test_cases[][2] = {
        {8192*2, 8192*2},
        {16, 16},
        {128, 128},
        {512, 512},
        {1536, 1024},
        {1536, 2048},
        {4096, 2048},
        {4096, 4096},
        {8192, 4096},
        {8192, 8192},
    };
    const int case_n = (int)(sizeof(test_cases) / sizeof(test_cases[0]));

    printf("\n===== Speedup vs 传统矩阵乘法 =====\n");
    printf("LUT计时=生成LUT+计算（不含transpose）\n");
    printf("%-12s %-12s %-12s %-12s %-12s\n",
        "M,K", "全精度", "本文LUT16", "传统LUT1", "SIMD");

    for (int ci = 0; ci < case_n; ci++) {
        M = test_cases[ci][0];
        K = test_cases[ci][1];

        const int blk_n = (K + QK_K - 1) / QK_K;
        const int act_pad_len = blk_n * QK_K * 2;

        block_ifairy *w = calloc((size_t)M * blk_n, sizeof(block_ifairy));
        float *act_raw = calloc((size_t)K * 2, sizeof(float));
        float *act = calloc((size_t)act_pad_len, sizeof(float));

        float *dst_fp = calloc((size_t)M * 2, sizeof(float));
        float *dst_lut16 = calloc((size_t)M * 2, sizeof(float));
        float *dst_lut1 = calloc((size_t)M * 2, sizeof(float));
        float *dst_trad = calloc((size_t)M * 2, sizeof(float));
        float *dst_simd = calloc((size_t)M * 2, sizeof(float));

        if (!w || !act_raw || !act || !dst_fp || !dst_lut16 || !dst_lut1 || !dst_trad || !dst_simd) {
            fprintf(stderr, "OOM at case M=%d K=%d\n", M, K);
            free(w); free(act_raw); free(act);
            free(dst_fp); free(dst_lut16); free(dst_lut1); free(dst_trad); free(dst_simd);
            return;
        }

        int rc = load_w_bit(matrix_file, w);
        if (rc != 0) {
            fprintf(stderr, "load_w_bit failed at M=%d K=%d rc=%d\n", M, K, rc);
            free(w); free(act_raw); free(act);
            free(dst_fp); free(dst_lut16); free(dst_lut1); free(dst_trad); free(dst_simd);
            return;
        }
        rc = load_act_float(act_file, act_raw, K * 2);
        if (rc != 0) {
            fprintf(stderr, "load_act_float failed at M=%d K=%d rc=%d\n", M, K, rc);
            free(w); free(act_raw); free(act);
            free(dst_fp); free(dst_lut16); free(dst_lut1); free(dst_trad); free(dst_simd);
            return;
        }
        memcpy(act, act_raw, (size_t)K * 2 * sizeof(float));

        block_ifairy_q16 *act_block_42 = alloc_act_block_ifairy_q16(K);
        block_ifairy_q16 *act_block_127 = alloc_act_block_ifairy_q16(K);
        lut_block *lut = alloc_lut_q8(K);
        block_ifairy_1x3 *w_t = alloc_w(M, K);

        int16_t *lut_v_old = alloc_lut_v_q16_old(K);
        float *lut_scale_old = alloc_lut_scale_old(K);
        block_ifairy_1x3_old *w_old = alloc_w_old(M, K);

        if (!act_block_42 || !act_block_127 || !lut || !w_t || !lut_v_old || !lut_scale_old || !w_old) {
            fprintf(stderr, "tmp alloc failed at M=%d K=%d\n", M, K);
            free(w); free(act_raw); free(act);
            free(dst_fp); free(dst_lut16); free(dst_lut1); free(dst_trad); free(dst_simd);
            free_act_block_ifairy_q16(act_block_42);
            free_act_block_ifairy_q16(act_block_127);
            free_lut_q8(lut);
            free_w(w_t);
            free_lut_v_q16_old(lut_v_old);
            free_lut_scale_old(lut_scale_old);
            free_w_old(w_old);
            return;
        }

        // transpose 不计入 LUT 时间
        transpose(M, K, w, w_t);
        transpose_old(M, K, w, w_old);

        long long t0, t1;
        double t_fp, t_lut16, t_lut1, t_trad, t_simd;

        t0 = now_ns();
        mul_mat_mxk_kx1(K, 0, M - 1, w, act, dst_fp);
        t1 = now_ns();
        t_fp = (double)(t1 - t0);

        act_float_2_block_ifairy_q16(K, act, act_block_42, 42.6f);
        t0 = now_ns();
        generate_lut_q8_block_ifairy_q16(K, act_block_42, lut);
        mul_mat_mxk_kx1_with_lut_q8(K, 0, M - 1, w_t, lut, dst_lut16);
        t1 = now_ns();
        t_lut16 = (double)(t1 - t0);

        act_float_2_block_ifairy_q16(K, act, act_block_127, 127.0f);
        t0 = now_ns();
        generate_lut_q16_block_ifairy_q16_old(K, act_block_127, lut_v_old, lut_scale_old);
        mul_mat_mxk_kx1_with_lut_q16_old(K, 0, M - 1, w_old, lut_v_old, lut_scale_old, dst_lut1);
        t1 = now_ns();
        t_lut1 = (double)(t1 - t0);

        t0 = now_ns();
        mul_mat_mxk_kx1_q8(K, 0, M - 1, w, act_block_127, dst_trad);
        t1 = now_ns();
        t_trad = (double)(t1 - t0);

        t0 = now_ns();
        mul_mat_mxk_kx1_q8_simd(K, 0, M - 1, w, act_block_127, dst_simd);
        t1 = now_ns();
        t_simd = (double)(t1 - t0);

        char mk[64];
        snprintf(mk, sizeof(mk), "%d*%d", M, K);
        printf("%-12s %-12.4f %-12.4f %-12.4f %-12.4f\n",
                mk,
                t_trad / t_fp,
                t_trad / t_lut16,
                t_trad / t_lut1,
                t_trad / t_simd);

        free_act_block_ifairy_q16(act_block_42);
        free_act_block_ifairy_q16(act_block_127);
        free_lut_q8(lut);
        free_w(w_t);
        free_lut_v_q16_old(lut_v_old);
        free_lut_scale_old(lut_scale_old);
        free_w_old(w_old);

        free(w); free(act_raw); free(act);
        free(dst_fp); free(dst_lut16); free(dst_lut1); free(dst_trad); free(dst_simd);
    }
}


void benchmark_metrics_vs_fp_from_file(const char *matrix_file, const char *act_file) {
    const int test_cases[][2] = {
        {16, 16},
        {128, 128},
        {512, 512},
        {1536, 1024},
        {1536, 2048},
        {4096, 2048},
        {4096, 4096},
        {8192, 4096},
        {8192, 8192},
        {8192*2, 8192*2},
    };
    const int case_n = (int)(sizeof(test_cases) / sizeof(test_cases[0]));

    float mre_lut16[16] = {0}, mre_trad[16] = {0};
    float xre_lut16[16] = {0}, xre_trad[16] = {0};
    float cos_lut16[16] = {0}, cos_trad[16] = {0};
    float snr_lut16[16] = {0}, snr_trad[16] = {0};
    float acc01_lut16[16] = {0}, acc01_trad[16] = {0};
    float maxmis_lut16[16] = {0}, maxmis_trad[16] = {0};
    char dims[16][32];

    printf("\n===== Metrics vs 全精度基准 =====\n");
    printf("方法: 本文加速方法(LUT16) / 传统矩阵乘法(q8)\n");

    for (int ci = 0; ci < case_n; ci++) {
        M = test_cases[ci][0];
        K = test_cases[ci][1];

        const int blk_n = (K + QK_K - 1) / QK_K;
        const int act_pad_len = blk_n * QK_K * 2;

        block_ifairy *w = calloc((size_t)M * blk_n, sizeof(block_ifairy));
        float *act_raw = calloc((size_t)K * 2, sizeof(float));
        float *act = calloc((size_t)act_pad_len, sizeof(float));

        float *dst_fp = calloc((size_t)M * 2, sizeof(float));
        float *dst_lut16 = calloc((size_t)M * 2, sizeof(float));
        float *dst_trad = calloc((size_t)M * 2, sizeof(float));

        if (!w || !act_raw || !act || !dst_fp || !dst_lut16 || !dst_trad) {
            fprintf(stderr, "OOM at case M=%d K=%d\n", M, K);
            free(w); free(act_raw); free(act);
            free(dst_fp); free(dst_lut16); free(dst_trad);
            return;
        }

        int rc = load_w_bit(matrix_file, w);
        if (rc != 0) {
            fprintf(stderr, "load_w_bit failed at M=%d K=%d rc=%d\n", M, K, rc);
            free(w); free(act_raw); free(act);
            free(dst_fp); free(dst_lut16); free(dst_trad);
            return;
        }
        rc = load_act_float(act_file, act_raw, K * 2);
        if (rc != 0) {
            fprintf(stderr, "load_act_float failed at M=%d K=%d rc=%d\n", M, K, rc);
            free(w); free(act_raw); free(act);
            free(dst_fp); free(dst_lut16); free(dst_trad);
            return;
        }
        memcpy(act, act_raw, (size_t)K * 2 * sizeof(float));

        block_ifairy_q16 *act_block_42 = alloc_act_block_ifairy_q16(K);
        block_ifairy_q16 *act_block_127 = alloc_act_block_ifairy_q16(K);
        lut_block *lut = alloc_lut_q8(K);
        block_ifairy_1x3 *w_t = alloc_w(M, K);

        if (!act_block_42 || !act_block_127 || !lut || !w_t) {
            fprintf(stderr, "tmp alloc failed at M=%d K=%d\n", M, K);
            free(w); free(act_raw); free(act);
            free(dst_fp); free(dst_lut16); free(dst_trad);
            free_act_block_ifairy_q16(act_block_42);
            free_act_block_ifairy_q16(act_block_127);
            free_lut_q8(lut);
            free_w(w_t);
            return;
        }

        // 全精度基准
        mul_mat_mxk_kx1(K, 0, M - 1, w, act, dst_fp);

        // 本文加速方法（LUT16）
        transpose(M, K, w, w_t);
        act_float_2_block_ifairy_q16(K, act, act_block_42, 42.6f);
        generate_lut_q8_block_ifairy_q16(K, act_block_42, lut);
        mul_mat_mxk_kx1_with_lut_q8(K, 0, M - 1, w_t, lut, dst_lut16);

        // 传统矩阵乘法（q8）
        act_float_2_block_ifairy_q16(K, act, act_block_127, 127.0f);
        mul_mat_mxk_kx1_q8(K, 0, M - 1, w, act_block_127, dst_trad);

        float mre_lut16_vs_ref = mean_relative_error_vs_ref(dst_lut16, dst_fp);
        float mre_trad_vs_ref = mean_relative_error_vs_ref(dst_trad, dst_fp);
        float xre_lut16_vs_ref = max_relative_error_vs_ref(dst_lut16, dst_fp);
        float xre_trad_vs_ref = max_relative_error_vs_ref(dst_trad, dst_fp);
        float cos_lut16_vs_ref = cosine_similarity_vs_ref(dst_lut16, dst_fp);
        float cos_trad_vs_ref = cosine_similarity_vs_ref(dst_trad, dst_fp);
        float snr_lut16_vs_ref = snr_db_vs_ref(dst_lut16, dst_fp);
        float snr_trad_vs_ref = snr_db_vs_ref(dst_trad, dst_fp);

        int errors_lut16 = 0;
        float max_mismatch_lut16 = 0.0f;
        int errors_trad = 0;
        float max_mismatch_trad = 0.0f;
        for (int i = 0; i < M * 2; i++) {
            float mismatch_lut16 = fabsf(dst_fp[i] - dst_lut16[i]) / (fabsf(dst_fp[i]) + 1e-9f);
            if (mismatch_lut16 > 0.1f) {
                if (mismatch_lut16 > max_mismatch_lut16) {
                    max_mismatch_lut16 = mismatch_lut16;
                }
                errors_lut16++;
            }

            float mismatch_trad = fabsf(dst_fp[i] - dst_trad[i]) / (fabsf(dst_fp[i]) + 1e-9f);
            if (mismatch_trad > 0.1f) {
                if (mismatch_trad > max_mismatch_trad) {
                    max_mismatch_trad = mismatch_trad;
                }
                errors_trad++;
            }
        }
        float acc01_lut16_vs_ref = ((float)M * 2.0f - (float)errors_lut16) / ((float)M * 2.0f);
        float acc01_trad_vs_ref = ((float)M * 2.0f - (float)errors_trad) / ((float)M * 2.0f);

        snprintf(dims[ci], sizeof(dims[ci]), "%d*%d", M, K);
        mre_lut16[ci] = mre_lut16_vs_ref;
        mre_trad[ci] = mre_trad_vs_ref;
        xre_lut16[ci] = xre_lut16_vs_ref;
        xre_trad[ci] = xre_trad_vs_ref;
        cos_lut16[ci] = cos_lut16_vs_ref;
        cos_trad[ci] = cos_trad_vs_ref;
        snr_lut16[ci] = snr_lut16_vs_ref;
        snr_trad[ci] = snr_trad_vs_ref;
        acc01_lut16[ci] = acc01_lut16_vs_ref;
        acc01_trad[ci] = acc01_trad_vs_ref;
        maxmis_lut16[ci] = max_mismatch_lut16;
        maxmis_trad[ci] = max_mismatch_trad;

        free_act_block_ifairy_q16(act_block_42);
        free_act_block_ifairy_q16(act_block_127);
        free_lut_q8(lut);
        free_w(w_t);

        free(w); free(act_raw); free(act);
        free(dst_fp); free(dst_lut16); free(dst_trad);
    }

    printf("\n--- 平均相对误差(MRE) ---\n");
    printf("%-20s", "方法\\维度");
    for (int i = 0; i < case_n; i++) printf(" %-12s", dims[i]);
    printf("\n");
    printf("%-20s", "本文方法");
    for (int i = 0; i < case_n; i++) printf(" %-12.6f", mre_lut16[i]);
    printf("\n");
    printf("%-20s", "传统矩阵方法");
    for (int i = 0; i < case_n; i++) printf(" %-12.6f", mre_trad[i]);
    printf("\n");

    printf("\n--- 最大相对误差(XRE) ---\n");
    printf("%-20s", "方法\\维度");
    for (int i = 0; i < case_n; i++) printf(" %-12s", dims[i]);
    printf("\n");
    printf("%-20s", "本文方法");
    for (int i = 0; i < case_n; i++) printf(" %-12.6f", xre_lut16[i]);
    printf("\n");
    printf("%-20s", "传统矩阵方法");
    for (int i = 0; i < case_n; i++) printf(" %-12.6f", xre_trad[i]);
    printf("\n");

    printf("\n--- 余弦相似度(Cosine) ---\n");
    printf("%-20s", "方法\\维度");
    for (int i = 0; i < case_n; i++) printf(" %-12s", dims[i]);
    printf("\n");
    printf("%-20s", "本文方法");
    for (int i = 0; i < case_n; i++) printf(" %-12.6f", cos_lut16[i]);
    printf("\n");
    printf("%-20s", "传统矩阵方法");
    for (int i = 0; i < case_n; i++) printf(" %-12.6f", cos_trad[i]);
    printf("\n");

    printf("\n--- 信噪比SNR(dB) ---\n");
    printf("%-20s", "方法\\维度");
    for (int i = 0; i < case_n; i++) printf(" %-12s", dims[i]);
    printf("\n");
    printf("%-20s", "本文方法");
    for (int i = 0; i < case_n; i++) printf(" %-12.6f", snr_lut16[i]);
    printf("\n");
    printf("%-20s", "传统矩阵方法");
    for (int i = 0; i < case_n; i++) printf(" %-12.6f", snr_trad[i]);
    printf("\n");

    printf("\n--- 0.1 Accuracy ---\n");
    printf("%-20s", "方法\\维度");
    for (int i = 0; i < case_n; i++) printf(" %-12s", dims[i]);
    printf("\n");
    printf("%-20s", "本文方法");
    for (int i = 0; i < case_n; i++) printf(" %-12.6f", acc01_lut16[i]);
    printf("\n");
    printf("%-20s", "传统矩阵方法");
    for (int i = 0; i < case_n; i++) printf(" %-12.6f", acc01_trad[i]);
    printf("\n");

    printf("\n--- Max Mismatch ---\n");
    printf("%-20s", "方法\\维度");
    for (int i = 0; i < case_n; i++) printf(" %-12s", dims[i]);
    printf("\n");
    printf("%-20s", "本文方法");
    for (int i = 0; i < case_n; i++) printf(" %-12.6f", maxmis_lut16[i]);
    printf("\n");
    printf("%-20s", "传统矩阵方法");
    for (int i = 0; i < case_n; i++) printf(" %-12.6f", maxmis_trad[i]);
    printf("\n");
}

void sample(const block_ifairy *w, const float *act) {
    float *dst2 = calloc(K*2, sizeof(float));
    lut_block *lut = alloc_lut_q8(M);
    block_ifairy_1x3 *_w = alloc_w(M, K);
    transpose(M, K, w, _w);
    while(1) {
        generate_lut_q8(M, act, lut);
        mul_mat_mxk_kx1_with_lut_q8(K, 0, M-1, _w, lut, dst2);
    }
    free_w(_w);
    free_lut_q8(lut);
    free(dst2);
}

int main(int argc, char *argv[]) {
    if (argc == 2 && strcmp(argv[1], "bench") == 0) {
        benchmark_speedup_from_file("./test_data/w.txt", "./test_data/act.txt");
        return 0;
    }
    if (argc == 2 && strcmp(argv[1], "metrics") == 0) {
        benchmark_metrics_vs_fp_from_file("./test_data/w.txt", "./test_data/act.txt");
        return 0;
    }

    if (argc != 3) {
        return -1;
    }
    M = atoi(argv[1]);
    K = atoi(argv[2]);
    const char *matrix_file = "./test_data/w.txt";
    const char *act_file    = "./test_data/act.txt";

    block_ifairy *w = calloc((size_t)M * K/QK_K, sizeof(block_ifairy));
    if (!w) { fprintf(stderr, "OOM w\n"); return 1; }
    float *act = calloc((size_t)K * 2, sizeof(float));
    if (!act) { fprintf(stderr, "OOM act\n"); return 1; }
    int rc = load_w_bit(matrix_file, w);
    if (rc != 0) { fprintf(stderr, "Failed to read matrix (%d)\n", rc); free(w); return 2; }
    rc = load_act_float(act_file, act, K*2);
    if (rc != 0) { fprintf(stderr, "Failed to read act (%d)\n", rc); free(w); return 2; }
    // compare(w, act);
    // sample(w, act);
    free(w);
    free(act);
    // benchmark_speedup_from_file(matrix_file, act_file);
    benchmark_metrics_vs_fp_from_file(matrix_file, act_file);
    return 0;
}
