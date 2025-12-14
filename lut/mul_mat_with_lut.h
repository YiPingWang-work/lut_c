#ifndef MUL_MAT_WITH_LUT
#define MUL_MAT_WITH_LUT

#include <arm_neon.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>

#define QK_K 256
typedef uint16_t ggml_half;


typedef struct {
    uint8_t qs[QK_K/4]; // 2 bits per element
    ggml_half d_real, d_imag;
} block_ifairy;

int8x16x2_t *alloc_lut(int m);
void free_lut(int8x16x2_t *lut);
void generate_lut_int8(const int16_t *act, int m, int8x16x2_t *lut);
void mul_mat_nxm_mx1_with_lut(block_ifairy *weight, int block_n, int row_begin, int row_end, int8x16x2_t *lut, float *lut_scale, int32_t *dst);
void mul_mat_nxm_mx1(block_ifairy *weight, int block_n, int row_begin, int row_end, int16_t *act, int32_t *dst);

#endif