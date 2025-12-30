#ifndef MUL_MAT_WITH_LUT
#define MUL_MAT_WITH_LUT

#include <arm_neon.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <math.h>

#define QK_K 256
// typedef uint16_t ggml_half;


typedef struct {
    uint8_t qs[QK_K/4]; // 2 bits per element
    float d_real, d_imag;
} block_ifairy;

typedef struct {
    int8x16x4_t v[(QK_K+2)/3]; // 每3个复数一个lut条目(ac bd ad bc)，256个复数需要86个条目
    float d_real, d_imag;
} lut_block;

lut_block *alloc_lut(int rows);
void free_lut(lut_block *lut);
void generate_lut_int8(const float *act, int rows, lut_block *lut);
// weight((row_end-row_begin+1)*cols) x act(cols*2，实虚交错) = dst(rows*2，实虚交错)
void mul_mat_nxm_mx1_with_lut(const block_ifairy *weight, int cols, int row_begin, int row_end, const lut_block *lut, float *dst);
void mul_mat_nxm_mx1(const block_ifairy *weight, int block_n, int row_begin, int row_end, const float *act, float *dst);

#endif