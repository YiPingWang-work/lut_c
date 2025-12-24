#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include "lut/mul_mat_with_lut.h"
#include <time.h>


#define ROWS 10240
#define COLS 10240

void print_time_ms() {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);

    struct tm tm_now;
    localtime_r(&ts.tv_sec, &tm_now);

    printf("%04d-%02d-%02d %02d:%02d:%02d.%03ld\n",
           tm_now.tm_year + 1900,
           tm_now.tm_mon + 1,
           tm_now.tm_mday,
           tm_now.tm_hour,
           tm_now.tm_min,
           tm_now.tm_sec,
           ts.tv_nsec / 1000000);
}

static inline float rand_float_001_2() {
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
    while ((c = fgetc(fp)) != EOF && block_idx < (ROWS*COLS/QK_K)) {
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

        // 一个 block 完成：512 bit
        if (bit_cnt == QK_K * 2) {
            // w[block_idx].d_real = rand_float_001_2();
            // w[block_idx].d_imag = rand_float_001_2();

            w[block_idx].d_real = 1; // 测试代码
            w[block_idx].d_imag = 1;
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

    if (block_idx != (ROWS*COLS/QK_K)) {
        fprintf(stderr, "warning: loaded %d / %d blocks\n",
                block_idx, (ROWS*COLS/QK_K));
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
    float *dst1 = calloc(COLS*2, sizeof(float));
    float *dst2 = calloc(COLS*2, sizeof(float));

    // printf("查表计算:\n");
    // print_time_ms();
    // lut_block *lut = alloc_lut(ROWS);
    // generate_lut_int8(act, ROWS, lut);
    // mul_mat_nxm_mx1_with_lut(w, COLS, 0, ROWS-1, lut, dst2);
    // print_time_ms();
    // free_lut(lut);
    
    printf("查表计算，12组一量化:\n");
    print_time_ms();
    lut_block_12 *lut_12 = alloc_lut_12(ROWS);
    generate_lut_int8_12(act, ROWS, lut_12);
    mul_mat_nxm_mx1_with_lut_12(w, COLS, 0, ROWS-1, lut_12, dst2);
    print_time_ms();
    free_lut_12(lut_12);


    printf("普通计算:\n");
    print_time_ms();
    mul_mat_nxm_mx1(w, COLS, 0, ROWS-1, act, dst1);
    print_time_ms();

    
    int errors = 0;
    for (int i = 0; i < COLS*2; i++) {
        if ((dst1[i] - dst2[i])/(fabsf(dst1[i])+1e-9) > 3e-1) {
            printf("❌ %d: dst1=%f, dst2=%f\n", i, dst1[i], dst2[i]);
            errors++;
        } else {
            // printf("✅ %d: dst1=%f, dst2=%f\n", i, dst1[i], dst2[i]);
        }
    }
    if (errors == 0) {
        printf("Results match!\n");
    } else {
        printf("Total mismatches: %d\n", errors);
    }

    free(dst1);
    free(dst2);
}

void sample(const block_ifairy *w, const float *act) {
    float *dst2 = calloc(COLS*2, sizeof(float));
    lut_block *lut = alloc_lut(ROWS);
    while(1) {
        generate_lut_int8(act, ROWS, lut);
        mul_mat_nxm_mx1_with_lut(w, COLS, 0, ROWS-1, lut, dst2);
    }
    free_lut(lut);
    free(dst2);
}

int main() {
    const char *matrix_file = "./test_data/w.txt";
    const char *act_file    = "./test_data/act.txt";

    block_ifairy *w = calloc((size_t)ROWS * 4, sizeof(block_ifairy));
    if (!w) { fprintf(stderr, "OOM w\n"); return 1; }
    float *act = calloc((size_t)ROWS * 2, sizeof(float));
    if (!act) { fprintf(stderr, "OOM act\n"); return 1; }
    int rc = load_w_bit(matrix_file, w);
    if (rc != 0) { fprintf(stderr, "Failed to read matrix (%d)\n", rc); free(w); return 2; }
    rc = load_act_float(act_file, act, ROWS*2);
    if (rc != 0) { fprintf(stderr, "Failed to read act (%d)\n", rc); free(w); return 2; }
    compare(w, act);
    // sample(w, act);
    free(w);
    free(act);
    return 0;
}