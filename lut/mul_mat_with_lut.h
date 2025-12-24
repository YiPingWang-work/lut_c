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
    int8x16x2_t v[(QK_K+2)/3]; // 每3个复数一个lut条目(实部+虚部)，256个复数需要86个条目
    float d_real, d_imag;
} lut_block;


typedef struct {
    int8x16x2_t v[(QK_K+2)/3]; // 每3个复数一个lut条目(实部+虚部)，256个复数需要86个条目
    float16_t d_real[(QK_K+11)/12], d_imag[(QK_K+11)/12]; // 每12个复数一个缩放因子(实部+虚部)，256个复数需要22个条目
} lut_block_12;

lut_block *alloc_lut(int m);
void free_lut(lut_block *lut);
void generate_lut_int8(const float *act, int m, lut_block *lut);
void mul_mat_nxm_mx1_with_lut(const block_ifairy *weight, int cols, int row_begin, int row_end, const lut_block *lut, float *dst);

// 测试
void generate_lut_int8_12(const float *act, int m, lut_block_12 *lut);
void mul_mat_nxm_mx1_with_lut_12(const block_ifairy *weight, int cols, int row_begin, int row_end, const lut_block_12 *lut, float *dst);
lut_block_12 *alloc_lut_12(int m);
void free_lut_12(lut_block_12 *lut);

// 校准代码
void mul_mat_nxm_mx1(const block_ifairy *weight, int block_n, int row_begin, int row_end, const float *act, float *dst);

#endif