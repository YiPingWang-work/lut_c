//
// Created by Yiping Wang on 2025/11/27.
//
#include "mul_mat_with_lut.h"

static inline int8_t set_int8(int16_t v) __attribute__((always_inline));
static inline int8_t set_int8(int16_t v) {
    if (v > 127) return 127;
    if (v < -128) return -128;
    return (int8_t)v;
}


void generate_lut_int8(const float *act, int m, lut_block *lut) { // m个复数
    int block_n = (m+QK_K-1)/QK_K;
    for (int block = 0; block < block_n; block++) {
        int act_begin = block * QK_K * 2, act_end = ((block+1)*QK_K-1)*2 < m*2 ? (block+1)*QK_K*2-1 : m*2;
        float max_real = 0.0f, max_imag = 0.0f;
        // 计算缩放因子
        for (int i = act_begin; i <= act_end; i += 2) {
            float real = fabsf(act[i]);
            float imag = fabsf(act[i+1]);
            if (real > max_real) max_real = real;
            if (imag > max_imag) max_imag = imag;
        }
        float scale_real = max_real / 42.0f;
        float scale_imag = max_imag / 42.0f;
        lut[block].d_real = scale_real;
        lut[block].d_imag = scale_imag;
        float inv_scale_real = 1.0f / scale_real;
        float inv_scale_imag = 1.0f / scale_imag;

        // 生成lut表
        for (int i = act_begin; i <= act_end; i += 6) {
            int8_t r0, r1, r2, i0, i1, i2;
            if (i+6 > act_end) {
                r0 = (int8_t)roundf(act[i]   * inv_scale_real);
                i0 = (int8_t)roundf(-act[i+1] * inv_scale_imag);
                r1 = 0;
                i1 = 0;
                r2 = 0;
                i2 = 0;
            } else {
                r0 = (int8_t)roundf(act[i]   * inv_scale_real);
                i0 = (int8_t)roundf(-act[i+1] * inv_scale_imag);
                r1 = (int8_t)roundf(act[i+2] * inv_scale_real);
                i1 = (int8_t)roundf(-act[i+3] * inv_scale_imag);
                r2 = (int8_t)roundf(act[i+4] * inv_scale_real);
                i2 = (int8_t)roundf(-act[i+5] * inv_scale_imag);
            }
            int8x16_t real, imag;
            real[0]  = -r0 - r1 - r2;   imag[0]  = -i0 - i1 - i2;    // (-1, -1, -1)
            real[1]  = -r0 - r1 + r2;   imag[1]  = -i0 - i1 + i2;    // (-1, -1,  1)
            real[2]  = -r0 - r1 + i2;   imag[2]  = -i0 - i1 - r2;    // (-1, -1, -i)
            real[3]  = -r0 - r1 - i2;   imag[3]  = -i0 - i1 + r2;    // (-1, -1,  i)

            real[4]  = -r0 + r1 - r2;   imag[4]  = -i0 + i1 - i2;    // (-1,  1, -1)
            real[5]  = -r0 + r1 + r2;   imag[5]  = -i0 + i1 + i2;    // (-1,  1,  1)
            real[6]  = -r0 + r1 + i2;   imag[6]  = -i0 + i1 - r2;    // (-1,  1, -i)
            real[7]  = -r0 + r1 - i2;   imag[7]  = -i0 + i1 + r2;    // (-1,  1,  i)

            real[8]  = -r0 + i1 - r2;   imag[8]  = -i0 - r1 - i2;    // (-1, -i, -1)
            real[9]  = -r0 + i1 + r2;   imag[9]  = -i0 - r1 + i2;    // (-1, -i,  1)
            real[10] = -r0 + i1 + i2;   imag[10] = -i0 - r1 - r2;    // (-1, -i, -i)
            real[11] = -r0 + i1 - i2;   imag[11] = -i0 - r1 + r2;    // (-1, -i,  i)

            real[12] = -r0 - i1 - r2;   imag[12] = -i0 + r1 - i2;    // (-1,  i, -1)
            real[13] = -r0 - i1 + r2;   imag[13] = -i0 + r1 + i2;    // (-1,  i,  1)
            real[14] = -r0 - i1 + i2;   imag[14] = -i0 + r1 - r2;    // (-1,  i, -i)
            real[15] = -r0 - i1 - i2;   imag[15] = -i0 + r1 + r2;    // (-1,  i,  i)

            lut[block].v[(i - act_begin)/6] = (int8x16x2_t){.val = {real, imag}};
        }

    }
}


const uint8_t three_vals2index[64] = {
    // -1 * *
    0b000000, // (-1, -1, -1)
    0b000001, // (-1, -1,  1)
    0b000010, // (-1, -1, -i)
    0b000011, // (-1, -1,  i)
    0b000100, // (-1,  1, -1)
    0b000101, // (-1,  1,  1)
    0b000110, // (-1,  1, -i)
    0b000111, // (-1,  1,  i)
    0b001000, // (-1, -i, -1)
    0b001001, // (-1, -i,  1)
    0b001010, // (-1, -i, -i)
    0b001011, // (-1, -i,  i)
    0b001100, // (-1,  i, -1)
    0b001101, // (-1,  i,  1)
    0b001110, // (-1,  i, -i)
    0b001111, // (-1,  i,  i)
    // 1 * *
    0b000101, // ( 1, -1, -1)
    0b000100, // ( 1, -1,  1)
    0b000111, // ( 1, -1, -i)
    0b000110, // ( 1, -1,  i)
    0b000001, // ( 1,  1, -1)
    0b000000, // ( 1,  1,  1)
    0b000011, // ( 1,  1, -i)
    0b000010, // ( 1,  1,  i)
    0b001101, // ( 1, -i, -1)
    0b001100, // ( 1, -i,  1)
    0b001111, // ( 1, -i, -i)
    0b001110, // ( 1, -i,  i)
    0b001001, // ( 1,  i, -1)
    0b001000, // ( 1,  i,  1)
    0b001011, // ( 1,  i, -i)
    0b001010, // ( 1,  i,  i)
    // -i * *
    0b001111, // (-i, -1, -1)
    0b001110, // (-i, -1,  1)
    0b001100, // (-i, -1, -i)
    0b001101, // (-i, -1,  i)
    0b001011, // (-i,  1, -1)
    0b001010, // (-i,  1,  1)
    0b001000, // (-i,  1, -i)
    0b001001, // (-i,  1,  i)
    0b000011, // (-i, -i, -1)
    0b000010, // (-i, -i,  1)
    0b000000, // (-i, -i, -i)
    0b000001, // (-i, -i,  i)
    0b000111, // (-i,  i, -1)
    0b000110, // (-i,  i,  1)
    0b000100, // (-i,  i, -i)
    0b000101, // (-i,  i,  i)
    // i * *
    0b001010, // ( i, -1, -1)
    0b001011, // ( i, -1,  1)
    0b001001, // ( i, -1, -i)
    0b001000, // ( i, -1,  i)
    0b001110, // ( i,  1, -1)
    0b001111, // ( i,  1,  1)
    0b001101, // ( i,  1, -i)
    0b001100, // ( i,  1,  i)
    0b000110, // ( i, -i, -1)
    0b000111, // ( i, -i,  1)
    0b000101, // ( i, -i, -i)
    0b000100, // ( i, -i,  i)
    0b000010, // ( i,  i, -1)
    0b000011, // ( i,  i,  1)
    0b000001, // ( i,  i, -i)
    0b000000, // ( i,  i,  i)
};




static inline int8x16x2_t mul_mat_block_16x3_3x1_with_lut(uint8x16_t iweight_16x3, int8x16x2_t ilut) __attribute__((always_inline));
static inline int8x16x2_t mul_mat_block_16x3_3x1_with_lut(uint8x16_t iweight_16x3, int8x16x2_t ilut) {
    int8x16x2_t ilut_exp[4] = {
        {.val = {ilut.val[0], ilut.val[1]}},
        {.val = {vnegq_s8(ilut.val[0]), vnegq_s8(ilut.val[1])}}, // act * -1
        {.val = {vnegq_s8(ilut.val[1]), ilut.val[0]}}, // act * i
        {.val = {ilut.val[1], vnegq_s8(ilut.val[0])}}, // act * -i
    };

    // 生成flag
    uint8x16_t lut_flag = vshrq_n_u8(vandq_u8(iweight_16x3, vdupq_n_u8(0b00110000)), 4); // 11 00 00
    // index 适配
    uint8x16_t index = {
        three_vals2index[iweight_16x3[0]],
        three_vals2index[iweight_16x3[1]],
        three_vals2index[iweight_16x3[2]],
        three_vals2index[iweight_16x3[3]],
        three_vals2index[iweight_16x3[4]],
        three_vals2index[iweight_16x3[5]],
        three_vals2index[iweight_16x3[6]],
        three_vals2index[iweight_16x3[7]],
        three_vals2index[iweight_16x3[8]],
        three_vals2index[iweight_16x3[9]],
        three_vals2index[iweight_16x3[10]],
        three_vals2index[iweight_16x3[11]],
        three_vals2index[iweight_16x3[12]],
        three_vals2index[iweight_16x3[13]],
        three_vals2index[iweight_16x3[14]],
        three_vals2index[iweight_16x3[15]],
    };

    // 查询lut
    int8x16_t real_exp_00 = vqtbl1q_s8(ilut_exp[0].val[0], index);
    int8x16_t real_exp_01 = vqtbl1q_s8(ilut_exp[1].val[0], index);
    int8x16_t real_exp_10 = vqtbl1q_s8(ilut_exp[2].val[0], index);
    int8x16_t real_exp_11 = vqtbl1q_s8(ilut_exp[3].val[0], index);
    int8x16_t imag_exp_00 = vqtbl1q_s8(ilut_exp[0].val[1], index);
    int8x16_t imag_exp_01 = vqtbl1q_s8(ilut_exp[1].val[1], index);
    int8x16_t imag_exp_10 = vqtbl1q_s8(ilut_exp[2].val[1], index);
    int8x16_t imag_exp_11 = vqtbl1q_s8(ilut_exp[3].val[1], index);

    // 补数据
    int8x16_t real = vdupq_n_s8(0);
    int8x16_t imag = vdupq_n_s8(0);
    real = vbslq_u8(vceqq_u8(lut_flag, vdupq_n_u8(0b00)), real_exp_00, real);
    real = vbslq_u8(vceqq_u8(lut_flag, vdupq_n_u8(0b01)), real_exp_01, real);
    real = vbslq_u8(vceqq_u8(lut_flag, vdupq_n_u8(0b10)), real_exp_10, real);
    real = vbslq_u8(vceqq_u8(lut_flag, vdupq_n_u8(0b11)), real_exp_11, real);
    imag = vbslq_u8(vceqq_u8(lut_flag, vdupq_n_u8(0b00)), imag_exp_00, imag);
    imag = vbslq_u8(vceqq_u8(lut_flag, vdupq_n_u8(0b01)), imag_exp_01, imag);
    imag = vbslq_u8(vceqq_u8(lut_flag, vdupq_n_u8(0b10)), imag_exp_10, imag);
    imag = vbslq_u8(vceqq_u8(lut_flag, vdupq_n_u8(0b11)), imag_exp_11, imag);

    return (int8x16x2_t){.val = {real, imag}};
}


static inline uint32x4x4_t get_iweight_1x12(const block_ifairy *w, int row, int block, int block_n, int i) __attribute__((always_inline));
static inline uint32x4x4_t get_iweight_1x12(const block_ifairy *w, int row, int block, int block_n, int i) {
    uint32x4x4_t iweight_1x12;
    const block_ifairy *w0;
    uint32x4_t v0, v1, v2;
    if (i + 3 >= QK_K / 4) {
        for (int k = 0; k < 4; k++) {
            w0 = &w[(row + 4*k) * block_n + block];
            v0 = (uint32x4_t){
                w0[0 * block_n].qs[i],
                w0[1 * block_n].qs[i],
                w0[2 * block_n].qs[i],
                w0[3 * block_n].qs[i],
            };
            iweight_1x12.val[k] = vshlq_n_u32(v0, 16);
        }
    } else {
        for (int k = 0; k < 4; k++) {
            w0 = &w[(row + 4*k) * block_n + block];
            v0 = (uint32x4_t){
                w0[0 * block_n].qs[i],
                w0[1 * block_n].qs[i],
                w0[2 * block_n].qs[i],
                w0[3 * block_n].qs[i],
            };
            v1 = (uint32x4_t){
                w0[0 * block_n].qs[i + 1],
                w0[1 * block_n].qs[i + 1],
                w0[2 * block_n].qs[i + 1],
                w0[3 * block_n].qs[i + 1],
            };
            v2 = (uint32x4_t){
                w0[0 * block_n].qs[i + 2],
                w0[1 * block_n].qs[i + 2],
                w0[2 * block_n].qs[i + 2],
                w0[3 * block_n].qs[i + 2],
            };
            iweight_1x12.val[k] = vorrq_u32(vorrq_u32(vshlq_n_u32(v0, 16), vshlq_n_u32(v1, 8)), v2);
        }
    }
    return iweight_1x12;
}


static inline uint8x16_t get_iweight_16x3_shift_0(uint32x4x4_t iweight_1x12) __attribute__((always_inline));
static inline uint8x16_t get_iweight_16x3_shift_0(uint32x4x4_t iweight_1x12) {
    const uint32x4_t mask = vdupq_n_u32(0x3F);
    // shift >> 6, mask 6 bits
    uint32x4_t s0 = vandq_u32(vshrq_n_u32(iweight_1x12.val[0], 18), mask);
    uint32x4_t s1 = vandq_u32(vshrq_n_u32(iweight_1x12.val[1], 18), mask);
    uint32x4_t s2 = vandq_u32(vshrq_n_u32(iweight_1x12.val[2], 18), mask);
    uint32x4_t s3 = vandq_u32(vshrq_n_u32(iweight_1x12.val[3], 18), mask);
    // u32 -> u16
    uint16x4_t u0 = vmovn_u32(s0);
    uint16x4_t u1 = vmovn_u32(s1);
    uint16x4_t u2 = vmovn_u32(s2);
    uint16x4_t u3 = vmovn_u32(s3);
    // pack
    uint16x8_t t0 = vcombine_u16(u0, u1);
    uint16x8_t t1 = vcombine_u16(u2, u3);
    // u16 -> u8
    uint8x8_t b0 = vmovn_u16(t0);
    uint8x8_t b1 = vmovn_u16(t1);
    return vcombine_u8(b0, b1);
}


static inline uint8x16_t get_iweight_16x3_shift_1(uint32x4x4_t iweight_1x12) __attribute__((always_inline));
static inline uint8x16_t get_iweight_16x3_shift_1(uint32x4x4_t iweight_1x12) {
    const uint32x4_t mask = vdupq_n_u32(0x3F);
    // shift >> 6, mask 6 bits
    uint32x4_t s0 = vandq_u32(vshrq_n_u32(iweight_1x12.val[0], 12), mask);
    uint32x4_t s1 = vandq_u32(vshrq_n_u32(iweight_1x12.val[1], 12), mask);
    uint32x4_t s2 = vandq_u32(vshrq_n_u32(iweight_1x12.val[2], 12), mask);
    uint32x4_t s3 = vandq_u32(vshrq_n_u32(iweight_1x12.val[3], 12), mask);
    // u32 -> u16
    uint16x4_t u0 = vmovn_u32(s0);
    uint16x4_t u1 = vmovn_u32(s1);
    uint16x4_t u2 = vmovn_u32(s2);
    uint16x4_t u3 = vmovn_u32(s3);
    // pack
    uint16x8_t t0 = vcombine_u16(u0, u1);
    uint16x8_t t1 = vcombine_u16(u2, u3);
    // u16 -> u8
    uint8x8_t b0 = vmovn_u16(t0);
    uint8x8_t b1 = vmovn_u16(t1);
    return vcombine_u8(b0, b1);
}


static inline uint8x16_t get_iweight_16x3_shift_2(uint32x4x4_t iweight_1x12) __attribute__((always_inline));
static inline uint8x16_t get_iweight_16x3_shift_2(uint32x4x4_t iweight_1x12) {
    const uint32x4_t mask = vdupq_n_u32(0x3F);
    // shift >> 6, mask 6 bits
    uint32x4_t s0 = vandq_u32(vshrq_n_u32(iweight_1x12.val[0], 6), mask);
    uint32x4_t s1 = vandq_u32(vshrq_n_u32(iweight_1x12.val[1], 6), mask);
    uint32x4_t s2 = vandq_u32(vshrq_n_u32(iweight_1x12.val[2], 6), mask);
    uint32x4_t s3 = vandq_u32(vshrq_n_u32(iweight_1x12.val[3], 6), mask);
    // u32 -> u16
    uint16x4_t u0 = vmovn_u32(s0);
    uint16x4_t u1 = vmovn_u32(s1);
    uint16x4_t u2 = vmovn_u32(s2);
    uint16x4_t u3 = vmovn_u32(s3);
    // pack
    uint16x8_t t0 = vcombine_u16(u0, u1);
    uint16x8_t t1 = vcombine_u16(u2, u3);
    // u16 -> u8
    uint8x8_t b0 = vmovn_u16(t0);
    uint8x8_t b1 = vmovn_u16(t1);
    return vcombine_u8(b0, b1);
}

static inline uint8x16_t get_iweight_16x3_shift_3(uint32x4x4_t iweight_1x12) __attribute__((always_inline));
static inline uint8x16_t get_iweight_16x3_shift_3(uint32x4x4_t iweight_1x12) {
    const uint32x4_t mask = vdupq_n_u32(0x3F);
    // shift >> 6, mask 6 bits
    uint32x4_t s0 = vandq_u32(iweight_1x12.val[0], mask);
    uint32x4_t s1 = vandq_u32(iweight_1x12.val[1], mask);
    uint32x4_t s2 = vandq_u32(iweight_1x12.val[2], mask);
    uint32x4_t s3 = vandq_u32(iweight_1x12.val[3], mask);
    // u32 -> u16
    uint16x4_t u0 = vmovn_u32(s0);
    uint16x4_t u1 = vmovn_u32(s1);
    uint16x4_t u2 = vmovn_u32(s2);
    uint16x4_t u3 = vmovn_u32(s3);
    // pack
    uint16x8_t t0 = vcombine_u16(u0, u1);
    uint16x8_t t1 = vcombine_u16(u2, u3);
    // u16 -> u8
    uint8x8_t b0 = vmovn_u16(t0);
    uint8x8_t b1 = vmovn_u16(t1);
    return vcombine_u8(b0, b1);
}


void mul_mat_nxm_mx1_with_lut(const block_ifairy *w, int cols, int row_begin, int row_end, const lut_block *lut, float *dst) {
    int block_n = (cols + QK_K - 1) / QK_K;
    for (int row = row_begin; row <= row_end; row+=16) {
        for (int block = 0; block < block_n; block++) {
            int16x8x2_t block_dst_real = {{vdupq_n_s16(0), vdupq_n_s16(0)}};
            int16x8x2_t block_dst_imag = {{vdupq_n_s16(0), vdupq_n_s16(0)}};
            for (int i = 0; i < QK_K/4; i+=3) {
                int max_ii = (i+3 >= QK_K/4) ? 2 : 4;
                uint32x4x4_t iweight_1x12 = get_iweight_1x12(w, row, block, block_n, i);

                for (int ii = 0; ii < max_ii; ii++) {
                    
                    uint8x16_t iweight_16x3 = {};
                    if (ii == 0) {
                        iweight_16x3 = get_iweight_16x3_shift_0(iweight_1x12);
                    } else if (ii == 1) {
                        iweight_16x3 = get_iweight_16x3_shift_1(iweight_1x12);
                    } else if (ii == 2) {
                        iweight_16x3 = get_iweight_16x3_shift_2(iweight_1x12);
                    } else if (ii == 3) {
                        iweight_16x3 = get_iweight_16x3_shift_3(iweight_1x12);
                    }
                    int8x16x2_t iret_ri = mul_mat_block_16x3_3x1_with_lut(iweight_16x3, lut[block].v[4*i/3+ii]);
                    int8x16_t iret_r = iret_ri.val[0];
                    int8x16_t iret_i = iret_ri.val[1];
                    block_dst_real.val[0] = vaddw_s8(block_dst_real.val[0], vget_low_s8(iret_r));
                    block_dst_real.val[1] = vaddw_s8(block_dst_real.val[1], vget_high_s8(iret_r));
                    block_dst_imag.val[0] = vaddw_s8(block_dst_imag.val[0], vget_low_s8(iret_i));
                    block_dst_imag.val[1] = vaddw_s8(block_dst_imag.val[1], vget_high_s8(iret_i));
                }
            }
            // 反量化，写回
            for (int j = 0; j < 16; j++) {
                dst[(row+j)*2  ] += (float)(j < 8 ? block_dst_real.val[0][j] : block_dst_real.val[1][j-8]) * lut[block].d_real * w[(row+j)*block_n+block].d_real;
                dst[(row+j)*2+1] += (float)(j < 8 ? block_dst_imag.val[0][j] : block_dst_imag.val[1][j-8]) * lut[block].d_imag * w[(row+j)*block_n+block].d_imag;
            }
        }
    }
}


lut_block *alloc_lut(int m) {
    int block_n = (m+QK_K-1)/QK_K;
    return calloc(m, sizeof(lut_block));
}


void free_lut(lut_block *lut) {
    free(lut);
}



// ================================ 验证乘法程序 ================================
static inline float32x2_t mul_mat_block_1x4_4x1(uint8_t a, float32_t *b) __attribute__((always_inline));
static inline float32x2_t mul_mat_block_1x4_4x1(uint8_t a, float32_t *b) {
    float32_t acc_r = 0;
    float32_t acc_i = 0;

    for (int k = 0; k < 4; k++) {
        uint8_t code = (a >> (2 * (3-k))) & 0x3;

        float32_t xr = b[2 * k];
        float32_t xi = -b[2 * k + 1]; // 共轭

        switch (code) {
        case 0b00:  // -1
            acc_r += -xr;
            acc_i += -xi;
            break;
        case 0b01:  // +1
            acc_r += xr;
            acc_i += xi;
            break;
        case 0b10:  // -i
            acc_r += xi;
            acc_i += -xr;
            break;
        case 0b11:  // +i
            acc_r += -xi;
            acc_i += xr;
            break;
        }
    }
    return vset_lane_f32(acc_i, vset_lane_f32(acc_r, vdup_n_f32(0), 0), 1);
}


void mul_mat_nxm_mx1(const block_ifairy *w, int cols, int row_begin, int row_end, const float *act, float *dst) {
    int block_n = (cols + QK_K - 1) / QK_K;
    for (int row = row_begin; row <= row_end; row++) {
        for (int block = 0; block < block_n; block++) {
            for (int i = 0; i < QK_K/4; i++) {
                uint8_t w4 = w[row*block_n+block].qs[i];
                float32_t a4[8] = {
                    act[block*QK_K*2+i*8],
                    act[block*QK_K*2+i*8+1],
                    act[block*QK_K*2+i*8+2],
                    act[block*QK_K*2+i*8+3],
                    act[block*QK_K*2+i*8+4],
                    act[block*QK_K*2+i*8+5],
                    act[block*QK_K*2+i*8+6],
                    act[block*QK_K*2+i*8+7],
                };
                float32x2_t ret = mul_mat_block_1x4_4x1(w4, a4);
                dst[2*row]   += (float)ret[0] * w[row*block_n+block].d_real;
                dst[2*row+1] += (float)ret[1] * w[row*block_n+block].d_imag;
            }
        }
    }
}