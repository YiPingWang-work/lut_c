#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include "lut/mul_mat_with_lut.h"
#include <time.h>


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
        w[block_idx].d_real = rand_float() * 2;
        w[block_idx].d_imag = rand_float() * 7;

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

    for(int i = 0; i < M*2; i++) {
        dst1[i] = 0.0f;
        dst2[i] = 0.0f;
    }

    // 验证程序
    long long t0 = now_ns();
    mul_mat_mxk_kx1(K, 0, M-1, w, act, dst1);
    long long t1 = now_ns();
    printf("直接计算,  耗时: %lld us\n", (t1 - t0));
    
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
    printf("查表计算1, 耗时: %lld us, 生成LUT耗时: %lld us\n", (t2 - t1), (t1 - t0));
    free_lut_q8(lut);
    free_w(_w);
    

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
    printf("查表计算2, 耗时: %lld us, 生成LUT耗时: %lld us\n", (t2 - t1), (t1 - t0));
    
    
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
    printf("查表计算1: 0.1 accuracy: %f, max_mismatch: %f\n", ((float)M*2 - (float)errors)/((float)M*2), max_mismatch);
    
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
    printf("查表计算2: 0.1 accuracy: %f, max_mismatch: %f\n", ((float)M*2 - (float)errors)/((float)M*2), max_mismatch);


    free(dst1);
    free(dst2);
    free(dst3);
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
    compare(w, act);
    // sample(w, act);
    free(w);
    free(act);
    return 0;
}
