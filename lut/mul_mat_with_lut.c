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


void generate_lut_int8(const float *act, int rows, lut_block *lut) { // rows个复数
    int block_n = (rows+QK_K-1)/QK_K;
    for (int block = 0; block < block_n; block++) {
        int act_begin = block * QK_K * 2, act_end = ((block+1)*QK_K-1)*2 < rows*2 ? (block+1)*QK_K*2-1 : rows*2;
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
        // printf("block %d: act[%d~%d], scale_real=%f, scale_imag=%f\n", block, act_begin, act_end, scale_real, scale_imag);
        
        // 生成lut表
        for (int i = act_begin; i <= act_end; i += 6) {
            int8_t r0, r1, r2, i0, i1, i2;
            if (i+6 > act_end) {
                r0 = (int8_t)roundf(act[i]    * inv_scale_real);
                i0 = (int8_t)roundf(-act[i+1] * inv_scale_imag);
                r1 = 0;
                i1 = 0;
                r2 = 0;
                i2 = 0;
            } else {
                r0 = (int8_t)roundf(act[i]    * inv_scale_real);
                i0 = (int8_t)roundf(-act[i+1] * inv_scale_imag);
                r1 = (int8_t)roundf(act[i+2]  * inv_scale_real);
                i1 = (int8_t)roundf(-act[i+3] * inv_scale_imag);
                r2 = (int8_t)roundf(act[i+4]  * inv_scale_real);
                i2 = (int8_t)roundf(-act[i+5] * inv_scale_imag);
            }

            int8x16_t ac = {
                -r0 - r1 - r2,
                -r0 - r1 + r2,
                -r0 - r1 + 0,
                -r0 - r1 + 0,

                -r0 + r1 - r2,
                -r0 + r1 + r2,
                -r0 + r1 + 0,
                -r0 + r1 + 0,

                -r0 + 0  - r2,
                -r0 + 0  + r2,
                -r0 + 0  + 0,
                -r0 + 0  + 0,

                -r0 + 0  - r2,
                -r0 + 0  + r2,
                -r0 + 0  + 0,
                -r0 + 0  + 0,
            };

            int8x16_t bd = {
                0 + 0 + 0,
                0 + 0 + 0,
                0 + 0 - i2,
                0 + 0 + i2,

                0 + 0 + 0,
                0 + 0 + 0,
                0 + 0 - i2,
                0 + 0 + i2,

                0 - i1 + 0,
                0 - i1 + 0,
                0 - i1 - i2,
                0 - i1 + i2,

                0 + i1 + 0,
                0 + i1 + 0,
                0 + i1 - i2,
                0 + i1 + i2,
            };

            int8x16_t ad = {
                0 + 0 + 0,
                0 + 0 + 0,
                0 + 0 - r2,
                0 + 0 + r2,

                0 + 0 + 0,
                0 + 0 + 0,
                0 + 0 - r2,
                0 + 0 + r2,

                0 - r1 + 0,
                0 - r1 + 0,
                0 - r1 - r2,
                0 - r1 + r2,

                0 + r1 + 0,
                0 + r1 + 0,
                0 + r1 - r2,
                0 + r1 + r2,
            };
            
            int8x16_t bc = {
                -i0 - i1 - i2,
                -i0 - i1 + i2,
                -i0 - i1 + 0,
                -i0 - i1 + 0,

                -i0 + i1 - i2,
                -i0 + i1 + i2,
                -i0 + i1 + 0,
                -i0 + i1 + 0,

                -i0 + 0  - i2,
                -i0 + 0  + i2,
                -i0 + 0  + 0,
                -i0 + 0  + 0,

                -i0 + 0  - i2,
                -i0 + 0  + i2,
                -i0 + 0  + 0,
                -i0 + 0  + 0,
            };

            lut[block].v[(i - act_begin)/6] = (int8x16x4_t){.val = ac, bd, ad, bc};
        }

    }
}


static const uint8x16x4_t three_vals2index = {
    (uint8x16_t){
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
    },
    (uint8x16_t){
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
    },
    (uint8x16_t){
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
    },
    (uint8x16_t){
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
    }
};


static inline int8x16x4_t mul_mat_block_16x3_3x1_with_lut(uint8x16_t iweight_16x3, int8x16x4_t ilut) __attribute__((always_inline));
static inline int8x16x4_t mul_mat_block_16x3_3x1_with_lut(uint8x16_t iweight_16x3, int8x16x4_t ilut) {

    // index 适配
    uint8x16_t index = vqtbl4q_u8(three_vals2index, iweight_16x3);

    //查询lut
    int8x16_t ac_00 = vqtbl1q_s8(ilut.val[0], index); // ad_10
    int8x16_t ac_11 = vqtbl1q_s8(ilut.val[2], index); // ad_00
    int8x16_t bd_01 = vqtbl1q_s8(ilut.val[1], index); // bc_11
    int8x16_t bd_11 = vqtbl1q_s8(ilut.val[3], index); // bc_00

    int8x16_t ac_01 = vnegq_s8(ac_00);                 // ad_11
    int8x16_t ac_10 = vnegq_s8(ac_11);                 // ad_01
    int8x16_t bd_00 = vnegq_s8(bd_01);                 // bc_10
    int8x16_t bd_10 = vnegq_s8(bd_11);                 // bc_01

    // 拆 2-bit flag
    uint8x16_t fl0 = vtstq_u8(iweight_16x3, vdupq_n_u8(0b00010000));
    uint8x16_t fl1 = vtstq_u8(iweight_16x3, vdupq_n_u8(0b00100000));

    int8x16_t ac_l = vbslq_u8(fl0, ac_01, ac_00);
    int8x16_t ac_h = vbslq_u8(fl0, ac_11, ac_10);
    int8x16_t ad_l = vbslq_u8(fl0, ac_10, ac_11);
    int8x16_t bc_l = vbslq_u8(fl0, bd_10, bd_11);
    int8x16_t bc_h = vbslq_u8(fl0, bd_01, bd_00);
    int8x16_t bd_h = vbslq_u8(fl0, bd_11, bd_10);

    return (int8x16x4_t){{
        vbslq_u8(fl1, ac_h, ac_l),
        vbslq_u8(fl1, ac_l, ad_l),
        vbslq_u8(fl1, bc_h, bc_l),
        vbslq_u8(fl1, bd_h, bc_h)
    }};
}


void mul_mat_nxm_mx1_with_lut(const block_ifairy *w, int cols, int row_begin, int row_end, const lut_block *lut, float *dst) {
    const int block_n = (cols + QK_K - 1) / QK_K;
    for (int row = row_begin; row <= row_end; row+=16) {
        for (int block = 0; block < block_n; block++) {

            int16x8x2_t block_dst_00 = {{vdupq_n_s16(0), vdupq_n_s16(0)}};
            int16x8x2_t block_dst_01 = {{vdupq_n_s16(0), vdupq_n_s16(0)}};
            int16x8x2_t block_dst_10 = {{vdupq_n_s16(0), vdupq_n_s16(0)}};
            int16x8x2_t block_dst_11 = {{vdupq_n_s16(0), vdupq_n_s16(0)}};

            for (int i = 0; i < QK_K/4; i+=3) {
                uint32_t iweight_1x12[16];
                #define MUL_MAT_16X3_3X1_WITH_LUT(ii) do { \
                    int shift = (3-ii)*6; \
                    uint8x16_t iweight_16x3 = { \
                        (iweight_1x12[0]  >> shift), \
                        (iweight_1x12[1]  >> shift), \
                        (iweight_1x12[2]  >> shift), \
                        (iweight_1x12[3]  >> shift), \
                        (iweight_1x12[4]  >> shift), \
                        (iweight_1x12[5]  >> shift), \
                        (iweight_1x12[6]  >> shift), \
                        (iweight_1x12[7]  >> shift), \
                        (iweight_1x12[8]  >> shift), \
                        (iweight_1x12[9]  >> shift), \
                        (iweight_1x12[10] >> shift), \
                        (iweight_1x12[11] >> shift), \
                        (iweight_1x12[12] >> shift), \
                        (iweight_1x12[13] >> shift), \
                        (iweight_1x12[14] >> shift), \
                        (iweight_1x12[15] >> shift), \
                    }; \
                    iweight_16x3 = vandq_u8(iweight_16x3, vdupq_n_u8(0b00111111)); \
\
                    int8x16x4_t iret = mul_mat_block_16x3_3x1_with_lut(iweight_16x3, lut[block].v[4*i/3+ii]); \
\
                    block_dst_00.val[0] = vaddw_s8(block_dst_00.val[0],  vget_low_s8(iret.val[0])); \
                    block_dst_00.val[1] = vaddw_s8(block_dst_00.val[1], vget_high_s8(iret.val[0])); \
                    block_dst_01.val[0] = vaddw_s8(block_dst_01.val[0],  vget_low_s8(iret.val[1])); \
                    block_dst_01.val[1] = vaddw_s8(block_dst_01.val[1], vget_high_s8(iret.val[1])); \
                    block_dst_10.val[0] = vaddw_s8(block_dst_10.val[0],  vget_low_s8(iret.val[2])); \
                    block_dst_10.val[1] = vaddw_s8(block_dst_10.val[1], vget_high_s8(iret.val[2])); \
                    block_dst_11.val[0] = vaddw_s8(block_dst_11.val[0],  vget_low_s8(iret.val[3])); \
                    block_dst_11.val[1] = vaddw_s8(block_dst_11.val[1], vget_high_s8(iret.val[3])); \
                } while(0)

                if (i+3 >= QK_K/4) {
                    for (int j = 0; j < 16 && row+j <= row_end; j++) {
                        const uint8_t *p = &w[(row+j)*block_n + block].qs[i];
                        iweight_1x12[j] = (uint32_t)p[0] << 16;
                    }
                    MUL_MAT_16X3_3X1_WITH_LUT(0);
                    MUL_MAT_16X3_3X1_WITH_LUT(1);
                } else {
                    for (int j = 0; j < 16 && row+j <= row_end; j++) {
                        const uint8_t *p = &w[(row+j)*block_n + block].qs[i];
                        uint32_t v = *(const uint32_t *)p;
                        iweight_1x12[j] = __builtin_bswap32(v) >> 8;
                    }
                    MUL_MAT_16X3_3X1_WITH_LUT(0);
                    MUL_MAT_16X3_3X1_WITH_LUT(1);
                    MUL_MAT_16X3_3X1_WITH_LUT(2);
                    MUL_MAT_16X3_3X1_WITH_LUT(3);
                }
            }

            // 反量化，写回
            for (int j = 0; j < 8 && row+j <= row_end; j++) {
                dst[(row+j)*2  ] += (float)(block_dst_00.val[0][j]) * lut[block].d_real * w[(row+j)*block_n+block].d_real + (float)(block_dst_11.val[0][j]) * lut[block].d_imag * w[(row+j)*block_n+block].d_imag;
                dst[(row+j)*2+1] += (float)(block_dst_01.val[0][j]) * lut[block].d_real * w[(row+j)*block_n+block].d_imag + (float)(block_dst_10.val[0][j]) * lut[block].d_imag * w[(row+j)*block_n+block].d_real;
            }

            for (int j = 8; j < 16 && row+j <= row_end; j++) {
                dst[(row+j)*2  ] += (float)(block_dst_00.val[1][j-8]) * lut[block].d_real * w[(row+j)*block_n+block].d_real + (float)(block_dst_11.val[1][j-8]) * lut[block].d_imag * w[(row+j)*block_n+block].d_imag;
                dst[(row+j)*2+1] += (float)(block_dst_01.val[1][j-8]) * lut[block].d_real * w[(row+j)*block_n+block].d_imag + (float)(block_dst_10.val[1][j-8]) * lut[block].d_imag * w[(row+j)*block_n+block].d_real;
            }
        }
    }
}


lut_block *alloc_lut(int rows) {
    int block_n = (rows+QK_K-1)/QK_K;
    return calloc(block_n, sizeof(lut_block));
}


void free_lut(lut_block *lut) {
    free(lut);
}



// ================================ 验证乘法程序 ================================
static inline float32x2_t mul_mat_block_1x4_4x1(uint8_t a, float32_t *b, float32_t d_real, float32_t d_imag) __attribute__((always_inline));
static inline float32x2_t mul_mat_block_1x4_4x1(uint8_t a, float32_t *b, float32_t d_real, float32_t d_imag) {
    float32_t acc_r = 0;
    float32_t acc_i = 0;

    for (int k = 0; k < 4; k++) {
        uint8_t code = (a >> (2 * (3-k))) & 0x3;

        float32_t xr = b[2 * k];
        float32_t xi = -b[2 * k + 1]; // 共轭

        switch (code) {
        case 0b00:  // -1
            acc_r += -xr * d_real;
            acc_i += -xi * d_real;
            break;
        case 0b01:  // +1
            acc_r += xr * d_real;
            acc_i += xi * d_real;
            break;
        case 0b10:  // -i
            acc_r += xi * d_imag;
            acc_i += -xr * d_imag;
            break;
        case 0b11:  // +i
            acc_r += -xi * d_imag;
            acc_i += xr * d_imag;
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
                float32x2_t ret = mul_mat_block_1x4_4x1(w4, a4, w[row*block_n+block].d_real, w[row*block_n+block].d_imag);
                dst[2*row]   += (float)ret[0];
                dst[2*row+1] += (float)ret[1];
            }
        }
    }
}



// ================================ 验证乘法程序（使用dot） ================================
void mul_mat_nxm_mx1_vecdot(const block_ifairy *w, int cols, int row_begin, int row_end, const int8_t *act, float *dst) { // 使用simd并行乘法计算
    
}