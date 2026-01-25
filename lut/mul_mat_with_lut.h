#ifndef MUL_MAT_WITH_LUT
#define MUL_MAT_WITH_LUT

#include <arm_neon.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <math.h>
#include <stdalign.h>

#define __M__ALIGNED_16 // 要求M必须被16整除
// #define __M__ALIGNED_8  // 要求M必须被8整除但不被16整除
#define QK_K 256


typedef struct {
    uint8_t qs[QK_K/4]; // 2 bits per element
    float d_real, d_imag;
} block_ifairy;

typedef struct __attribute__((aligned(128))) {
    int8x16x4_t v[(QK_K+2)/3]; // 每3个复数一个lut条目(ac bd ad bc)，256个复数需要86个条目，8位
    float d_real, d_imag;
} lut_block;

typedef struct __attribute__((aligned(128))) {
    uint8x16x4_t v[(QK_K+2)/3*2]; // 每3个复数一个lut条目(ac bd ad bc)，256个复数需要86个条目,16位
    float d_real, d_imag;
} lut_block_q16;

typedef struct __attribute__((aligned(128))) {
    uint8_t qs[(QK_K+2)/3]; // 8 bits 3 elements
    float d_real, d_imag;
} block_ifairy_1x3_old;

typedef struct __attribute__((aligned(128))) {
    uint8_t qs[(QK_K+2)/3][16]; // 8 bits 3 elements, 16 rows
    float d_real[16];
    float d_imag[16];
} block_ifairy_1x3;

typedef struct __attribute__((aligned(128))) {
    uint8_t x_real[QK_K], x_imag[QK_K];
    float d_real, d_imag;
} block_ifairy_q16;

void act_float_2_block_ifairy_q16(int k, const float *act_float, block_ifairy_q16 *act_q16, float scale);


// === 16路lut优化
void generate_lut_q8(int k, const float *act, lut_block *lut);
void generate_lut_q8_block_ifairy_q16(int k, const block_ifairy_q16 *act, lut_block *lut);
void transpose(int m, int k, const block_ifairy *raw_w, block_ifairy_1x3 *w);
void mul_mat_mxk_kx1_with_lut_q8(int k, int row_begin, int row_end, const block_ifairy_1x3 *w, const lut_block *lut, float *dst);
void mul_mat_mxk_kxn_with_lut_q8_base(int k, int row_begin, int row_end, int col_begin, int col_end, const block_ifairy_1x3 *w, const lut_block *lut_base, float *dst_base);
void mul_mat_mxk_kxn_with_lut_q8(int k, int n, int row_begin, int row_end, const block_ifairy_1x3 *w, const lut_block *lut, float *dst);
block_ifairy_q16 *alloc_act_block_ifairy_q16(int k);
lut_block *alloc_lut_q8(int k);
block_ifairy_1x3 *alloc_w(int m, int k);
void free_act_block_ifairy_q16(block_ifairy_q16 *act);
void free_lut_q8(lut_block *lut);
void free_w(block_ifairy_1x3 *w);

// === 原程序测试
void generate_lut_q16_block_ifairy_q16_old(int k, const block_ifairy_q16 *act, int16_t *lut_v, float *lut_scale);
void transpose_old(int m, int k, const block_ifairy *raw_w, block_ifairy_1x3_old *w); 
void mul_mat_mxk_kx1_with_lut_q16_old(int k, int row_begin, int row_end, const block_ifairy_1x3_old *w, const int16_t *lut, const float *lut_scale, float *dst);
int16_t *alloc_lut_v_q16_old(int k);
float *alloc_lut_scale_old(int k);
block_ifairy_1x3_old *alloc_w_old(int m, int k);
void free_lut_v_q16_old(int16_t *lut);
void free_lut_scale_old(float* scale);
void free_w_old(block_ifairy_1x3_old *w);

// === 验证程序
void mul_mat_mxk_kx1(int k, int row_begin, int row_end, const block_ifairy *w, const float *act, float *dst);

#endif