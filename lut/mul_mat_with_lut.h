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

// typedef struct {
//     uint8_t qs[(QK_K+2)/3]; // 8 bits 3 elements
//     float d_real, d_imag;
// } block_ifairy_1x3_2;

typedef struct {
    uint8_t qs[(QK_K+2)/3][16]; // 8 bits 3 elements, 16 rows
    float d_real[16];
    float d_imag[16];
} block_ifairy_1x3;

lut_block *alloc_lut(int k);
block_ifairy_1x3 *alloc_new_w(int m, int k);
void free_lut(lut_block *lut);
void free_new_w(block_ifairy_1x3 *w);
void generate_lut_int8(int k, const float *act, lut_block *lut);
void mul_mat_nxm_mx1_with_lut(int k, int row_begin, int row_end, const block_ifairy_1x3 *w, const lut_block *lut, float *dst);
void transpose(int m, int k, const block_ifairy *raw_w, block_ifairy_1x3 *w); 

// === 原程序测试
// int16_t *alloc_lut_2(int n);
// block_ifairy_1x3_2 *alloc_new_w_2(int n, int m);
// float *alloc_lut_scale_2(int n);
// void free_lut_2(int16_t *lut);
// void free_new_w_2(block_ifairy_1x3_2 *w);
// void free_lut_scale_2(float* scale);
// void generate_lut_int8_2(int m, const float *act, int16_t *lut, float *lut_scale);
// void mul_mat_nxm_mx1_with_lut_2(int m, int row_begin, int row_end, const block_ifairy_1x3_2 *w, const int16_t *lut, const float *lut_scale, float *dst);
// void transpose_2(int n, int m, const block_ifairy *raw_w, block_ifairy_1x3_2 *w); 


// === 验证程序
void mul_mat_nxm_mx1(int m, int row_begin, int row_end, const block_ifairy *w, const float *act, float *dst);

#endif