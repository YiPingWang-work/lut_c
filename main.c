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

void test_lut(const block_ifairy *w, const int16_t *act, int m, int n) { // mxn * nx1
    int8x16x2_t *lut = alloc_lut(n);
    printf("%d\n", m);
    generate_lut_int8(act, m, lut);
    int32_t *dst = calloc(m, sizeof(int32_t));
    free_lut(lut);
}


void test() {
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

    test_lut(w, act, ROWS, ROWS);
    

    free(w);
    free(act);
    return;
}

int main() {
    uint8_t x = 0b10110011;
    int16_t y[8] = {1, 2, 3, 4, 5, 6, 7, 8};
    int32x2_t z = mul_mat_1x4_4x1(x, y)
    printf("%d %d\n", z[0], z[1]);
    return 0;
}