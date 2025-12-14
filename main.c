#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include "lut/mul_mat_with_lut.h"

static inline void set_block_elem(block_ifairy *b, int idx256, uint8_t code) {
    int byte_idx = idx256 >> 2; /* 4 elems per byte */
    int shift = (3 - (idx256 & 3)) * 2;
    b->qs[byte_idx] &= (uint8_t)~(3u << shift);
    b->qs[byte_idx] |= (uint8_t)((code & 3u) << shift);
}

int read_matrix_01_txt(const char *path, block_ifairy *mat, int rows) {
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    char *line = NULL;
    size_t linecap = 0;
    ssize_t linelen;
    int row = 0;

    while (row < rows && (linelen = getline(&line, &linecap, f)) != -1) {
        /* gather bits '0'/'1' only */
        char bits[1024 + 8];
        int bcnt = 0;
        for (ssize_t i = 0; i < linelen && bcnt < 1024; ++i) {
            if (line[i] == '0' || line[i] == '1') bits[bcnt++] = line[i];
        }
        if (bcnt < 1024) {
            /* allow lines maybe concatenated across file; keep reading until enough */
            while (bcnt < 1024 && (linelen = getline(&line, &linecap, f)) != -1) {
                for (ssize_t i = 0; i < linelen && bcnt < 1024; ++i) {
                    if (line[i] == '0' || line[i] == '1') bits[bcnt++] = line[i];
                }
            }
        }
        if (bcnt != 1024) {
            free(line);
            fclose(f);
            return -2; /* not enough bits for a row */
        }
        /* pack 512 codes (two bits each) into two blocks */
        for (int c = 0; c < 512; ++c) {
            int bit_idx = c * 2;
            uint8_t b0 = bits[bit_idx] - '0';
            uint8_t b1 = bits[bit_idx + 1] - '0';
            uint8_t code = (uint8_t)((b0 << 1) | b1);
            int block = (c >= 256) ? 1 : 0;
            int idx256 = c % 256;
            set_block_elem(&mat[row * 2 + block], idx256, code);
        }
        row++;
    }

    free(line);
    fclose(f);

    if (row != rows) return -3;
    return 0;
}

int read_act_int8_txt(const char *path, int16_t *vec, int num_values) {
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    int count = 0;
    long tmp;
    while (count < num_values && fscanf(f, "%ld", &tmp) == 1) {
        if (tmp < -32768L) tmp = -32768L;
        if (tmp >  32767L) tmp =  32767L;
        vec[count++] = (int16_t)tmp;
    }
    fclose(f);
    return (count == num_values) ? 0 : -2;
}

static void print_u8_binary(uint8_t v) {
    for (int i = 7; i >= 0; --i) putchar((v & (1u << i)) ? '1' : '0');
    putchar('\n');
}

int write_to_file_int8x16x2_t(const char *filename, int index, int8x16x2_t v) {
    FILE *fp = fopen(filename, "a");   // 追加写
    if (!fp) {
        perror("fopen");
        return -1;
    }

    int8_t buf0[16];
    int8_t buf1[16];

    vst1q_s8(buf0, v.val[0]);
    vst1q_s8(buf1, v.val[1]);

    fprintf(fp, "%d vec0: ", index);
    for (int i = 0; i < 16; i++) {
        fprintf(fp, "%4d", buf0[i]);
    }

    fprintf(fp, "\n%d vec1: ", index);
    for (int i = 0; i < 16; i++) {
        fprintf(fp, "%4d", buf1[i]);
    }
    fprintf(fp, "\n");

    fclose(fp);
    return 0;
}

void test_mul_mat_with_lut(const block_ifairy *w, const int16_t *act) { // mxn * nx1
    int8x16x2_t *lut = alloc_lut(1024);
    generate_lut_int8(act, 1024, lut, NULL);
    // 输出lut到文件中
    for (int i = 0; i < (1024+12)/3; i++) {
        write_to_file_int8x16x2_t("./test_data/lut.txt", i, lut[i]);
    }
    int32_t *dst = calloc(2048, sizeof(int32_t));
    mul_mat_nxm_mx1_with_lut(w, 2, 0, 1023, lut, NULL, dst);
    free_lut(lut);
    free(dst);
}

void test_mul_mat(const block_ifairy *w, const int16_t *act) {
    int32_t *dst = calloc(2048, sizeof(int32_t));
    mul_mat_nxm_mx1(w, 2, 0, 1023, act, dst);
    free(dst);
}


int test() {
    const char *matrix_file = "./test_data/matrix_1024x1024.txt";
    const char *act_file    = "./test_data/act_1024.txt";

    const int ROWS = 1024;

    block_ifairy *w = calloc((size_t)ROWS * 2, sizeof(block_ifairy));
    if (!w) { fprintf(stderr, "OOM w\n"); return 1; }
    int16_t *act = calloc((size_t)ROWS * 2, sizeof(int16_t));
    if (!act) { fprintf(stderr, "OOM act\n"); return 1; }
    int rc = read_matrix_01_txt(matrix_file, w, ROWS);
    if (rc != 0) { fprintf(stderr, "Failed to read matrix (%d)\n", rc); free(w); return 2; }
    rc = read_act_int8_txt(act_file, act, ROWS*2);
    if (rc != 0) { fprintf(stderr, "Failed to read act (%d)\n", rc); free(w); return 2; }

    test_mul_mat(w, act);
    printf("\n=====\n");
    test_mul_mat_with_lut(w, act);
    

    free(w);
    free(act);
    return 0;
}


// void test_block_lut() {
//     uint8x16_t iweight_16x3 = {
//         0b00000000, // -1 -1 -1
//         0b00101111, // -i i i = a+bi -> -1 1 1 * (-i) = a+bi *(-i) = -b+ai 
//         0b00011110, // 1 i -i -> -1 -i i
//         0b00000000,
//         0b00000000,
//         0b00000000,
//         0b00000000,
//         0b00000000,
//         0b00000000,
//         0b00000000,
//         0b00000000,
//         0b00000000,
//         0b00000000,
//         0b00000000,
//         0b00000000,
//         0b00000000,
//     };
//     int8x16x2_t ilut = {
//         (int8x16_t){1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16},
//         (int8x16_t){101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116},
//     };
//     int8x16x2_t ans = mul_mat_block_16x3_3x1_with_lut(iweight_16x3, ilut);
//     for (int i = 0; i < 16; i++) {
//         printf("%d: %d, %d\n", i, ans.val[0][i], ans.val[1][i]);
//     }
// }

int main() {
    // test_block_lut();
    return test();
}