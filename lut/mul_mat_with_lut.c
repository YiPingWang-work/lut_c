//
// Created by Yiping Wang on 2025/11/27.
//
#include "mul_mat_with_lut.h"


void generate_lut_int8(int k, const float *act, lut_block *lut) {
    int blk_n = (k+QK_K-1)/QK_K;
    for (int blk = 0; blk < blk_n; blk++) {
        int act_begin = blk * QK_K * 2, act_end = ((blk+1)*QK_K-1)*2 < k*2 ? (blk+1)*QK_K*2-1 : k*2;
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
        lut[blk].d_real = scale_real;
        lut[blk].d_imag = scale_imag;
        float inv_scale_real = 1.0f / scale_real;
        float inv_scale_imag = 1.0f / scale_imag;
        
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

            lut[blk].v[(i - act_begin)/6] = (int8x16x4_t){.val = ac, bd, ad, bc};
        }

    }
}


static const uint8_t three_vals2index_uint8[64] = {
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


static inline int8x16x4_t mul_mat_block_16x3_3x1_with_lut(uint8x16_t iweight_16x3, int8x16x4_t ilut) __attribute__((always_inline));
static inline int8x16x4_t mul_mat_block_16x3_3x1_with_lut(uint8x16_t iweight_16x3, int8x16x4_t ilut) {

    const uint8x16_t mask_idx = vdupq_n_u8(0b00111111);
    const uint8x16_t mask_b6  = vdupq_n_u8(0b01000000);
    const uint8x16_t mask_b7  = vdupq_n_u8(0b10000000);

    // index 适配
    uint8x16_t index = vandq_u8(iweight_16x3, mask_idx);

    // 设置标识
    uint8x16_t fl0 = vtstq_u8(iweight_16x3, mask_b6);
    uint8x16_t fl1 = vtstq_u8(iweight_16x3, mask_b7);

    //查询lut
    int8x16_t ac_00 = vqtbl1q_s8(ilut.val[0], index); // ad_10
    int8x16_t bd_01 = vqtbl1q_s8(ilut.val[1], index); // bc_11
    int8x16_t ac_11 = vqtbl1q_s8(ilut.val[2], index); // ad_00
    int8x16_t bd_11 = vqtbl1q_s8(ilut.val[3], index); // bc_00

    int8x16_t ac_01 = vnegq_s8(ac_00);                 // ad_11
    int8x16_t ac_10 = vnegq_s8(ac_11);                 // ad_01
    int8x16_t bd_00 = vnegq_s8(bd_01);                 // bc_10
    int8x16_t bd_10 = vnegq_s8(bd_11);                 // bc_01


    int8x16_t ac_l = vbslq_u8(fl0, ac_01, ac_00);
    int8x16_t ad_l = vbslq_u8(fl0, ac_10, ac_11);
    int8x16_t bc_l = vbslq_u8(fl0, bd_10, bd_11);
    int8x16_t ac_h = vbslq_u8(fl0, ac_11, ac_10);
    int8x16_t bc_h = vbslq_u8(fl0, bd_01, bd_00);
    int8x16_t bd_h = vbslq_u8(fl0, bd_11, bd_10);

    return (int8x16x4_t){{
        vbslq_u8(fl1, ac_h, ac_l),
        vbslq_u8(fl1, ac_l, ad_l),
        vbslq_u8(fl1, bc_h, bc_l),
        vbslq_u8(fl1, bd_h, bc_h)
    }};
}


void transpose(int m, int k, const block_ifairy *raw_w, block_ifairy_1x3 *w) {
    const int blk_n = (k+QK_K-1)/QK_K;
    for (int j = 0; j < m; j++) {
        for (int blk = 0; blk < blk_n; blk++) {
            int ii = 0;
            #pragma unroll
            for (int i = 0; i < QK_K/4; i+=3) {
                const uint32_t *p = (uint32_t *)&(raw_w[j*blk_n + blk].qs[i]);
                uint8_t iweight_1x3_0, iweight_1x3_1, iweight_1x3_2, iweight_1x3_3;
                if (i+3 >= QK_K/4) {
                    uint32_t iweight_1x12 = p[0] << 16;
                    iweight_1x3_0 = iweight_1x12 >> 18 & 0b00111111; 
                    iweight_1x3_1 = iweight_1x12 >> 12 & 0b00111111;
                    w[(j/16)*blk_n + blk].qs[ii  ][j%16] = three_vals2index_uint8[iweight_1x3_0] | (iweight_1x3_0 << 2 & 0b11000000);
                    w[(j/16)*blk_n + blk].qs[ii+1][j%16] = three_vals2index_uint8[iweight_1x3_1] | (iweight_1x3_1 << 2 & 0b11000000);
                } else {
                    uint32_t iweight_1x12 = __builtin_bswap32(*p) >> 8;
                    iweight_1x3_0 = iweight_1x12 >> 18 & 0b00111111; 
                    iweight_1x3_1 = iweight_1x12 >> 12 & 0b00111111;
                    iweight_1x3_2 = iweight_1x12 >> 6  & 0b00111111; 
                    iweight_1x3_3 = iweight_1x12       & 0b00111111;
                    w[(j/16)*blk_n + blk].qs[ii  ][j%16] = three_vals2index_uint8[iweight_1x3_0] | (iweight_1x3_0 << 2 & 0b11000000);
                    w[(j/16)*blk_n + blk].qs[ii+1][j%16] = three_vals2index_uint8[iweight_1x3_1] | (iweight_1x3_1 << 2 & 0b11000000);
                    w[(j/16)*blk_n + blk].qs[ii+2][j%16] = three_vals2index_uint8[iweight_1x3_2] | (iweight_1x3_2 << 2 & 0b11000000);
                    w[(j/16)*blk_n + blk].qs[ii+3][j%16] = three_vals2index_uint8[iweight_1x3_3] | (iweight_1x3_3 << 2 & 0b11000000);
                }
                ii += 4;
            }
            w[(j/16)*blk_n + blk].d_real[j%16] = raw_w[j*blk_n + blk].d_real;
            w[(j/16)*blk_n + blk].d_imag[j%16] = raw_w[j*blk_n + blk].d_imag;
        }
    }
}


void mul_mat_mxk_kx1_with_lut(int k, int row_begin, int row_end, const block_ifairy_1x3 *w, const lut_block *lut, float *dst) {
    const int blk_n = (k+QK_K-1)/QK_K;
    const int in_blk_n = (QK_K+2)/3;
    for (int row = row_begin; row <= row_end; row+=16) {
        const block_ifairy_1x3 *w_base = w + (row >> 4)*blk_n;
        for (int blk = 0; blk < blk_n; blk++) {
            
            const block_ifairy_1x3 *w_blk_base = w_base + blk;
            const lut_block *lut_blk_base = lut + blk;

            int16x8x2_t block_dst_00 = {{vdupq_n_s16(0), vdupq_n_s16(0)}};
            int16x8x2_t block_dst_01 = {{vdupq_n_s16(0), vdupq_n_s16(0)}};
            int16x8x2_t block_dst_10 = {{vdupq_n_s16(0), vdupq_n_s16(0)}};
            int16x8x2_t block_dst_11 = {{vdupq_n_s16(0), vdupq_n_s16(0)}};
            
            #pragma unroll
            for (int i = 0; i < in_blk_n; i++) {
                __builtin_prefetch(w_blk_base -> qs[i]);
                uint8x16_t iweight_16x3 = vld1q_u8(w_blk_base -> qs[i]);
                int8x16x4_t iret = mul_mat_block_16x3_3x1_with_lut(iweight_16x3, lut_blk_base -> v[i]);

                block_dst_00.val[0] = vaddw_s8(block_dst_00.val[0],  vget_low_s8(iret.val[0]));
                block_dst_00.val[1] = vaddw_s8(block_dst_00.val[1], vget_high_s8(iret.val[0]));
                block_dst_01.val[0] = vaddw_s8(block_dst_01.val[0],  vget_low_s8(iret.val[1]));
                block_dst_01.val[1] = vaddw_s8(block_dst_01.val[1], vget_high_s8(iret.val[1]));
                block_dst_10.val[0] = vaddw_s8(block_dst_10.val[0],  vget_low_s8(iret.val[2]));
                block_dst_10.val[1] = vaddw_s8(block_dst_10.val[1], vget_high_s8(iret.val[2]));
                block_dst_11.val[0] = vaddw_s8(block_dst_11.val[0],  vget_low_s8(iret.val[3]));
                block_dst_11.val[1] = vaddw_s8(block_dst_11.val[1], vget_high_s8(iret.val[3]));
            }

            // 反量化，写回
            const float lr = lut[blk].d_real;
            const float li = lut[blk].d_imag;
            const int w_row = row >> 4;
            const float *wr = w[w_row * blk_n + blk].d_real;
            const float *wi = w[w_row * blk_n + blk].d_imag;
            float *idst = dst + row * 2;
            for (int j = 0; j <= (row_end - row < 7 ? row_end - row : 7); j++) {
                float _wr = wr[j], _wi = wi[j];
                idst[j*2  ] += (float)block_dst_00.val[0][j] * (lr * _wr) + (float)block_dst_11.val[0][j] * (li * _wi);
                idst[j*2+1] += (float)block_dst_01.val[0][j] * (lr * _wi) + (float)block_dst_10.val[0][j] * (li * _wr);
            }
            if (row_end - row >= 8) {
                for (int j = 8; j <= (row_end - row < 15 ? row_end - row : 15); j++) {
                    int jj = j - 8;
                    float _wr = wr[j], _wi = wi[j];
                    idst[j*2  ] += (float)block_dst_00.val[1][jj] * (lr * _wr) + (float)block_dst_11.val[1][jj] * (li * _wi);
                    idst[j*2+1] += (float)block_dst_01.val[1][jj] * (lr * _wi) + (float)block_dst_10.val[1][jj] * (li * _wr);
                }
            }
        }
    }
}


lut_block *alloc_lut(int k) {
    int blk_n = (k+QK_K-1)/QK_K;
    return calloc(blk_n, sizeof(lut_block));
}


block_ifairy_1x3 *alloc_w(int m, int k) {
    int blk_n = (k+QK_K-1)/QK_K;
    int row_n = (m+15)/16;
    block_ifairy_1x3 *ptr = calloc(blk_n*row_n, sizeof(block_ifairy_1x3));
    return ptr;
}


void free_lut(lut_block *lut) {
    free(lut);
}


void free_w(block_ifairy_1x3 *w) {
    free(w);
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


void mul_mat_mxk_kx1(int k, int row_begin, int row_end, const block_ifairy *w, const float *act, float *dst) {
    int block_n = (k+QK_K-1)/QK_K;
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



// ================================ 原程序测试 ================================

static const int16_t sx[4] = { -1,  1,  0,  0 };
static const int16_t sy[4] = {  0,  0, -1,  1 };

void generate_lut_int8_old(int k, const float *act, int16_t *lut, float *lut_scale) {
    int block_n = (k+QK_K-1)/QK_K;
    int ii = 0;
    for (int blk = 0; blk < block_n; blk++) {
        int act_begin = blk * QK_K * 2, act_end = ((blk+1)*QK_K-1)*2 < k*2 ? (blk+1)*QK_K*2-1 : k*2;
        float max_real = 0.0f, max_imag = 0.0f;
        // 计算缩放因子
        for (int i = act_begin; i <= act_end; i += 2) {
            float real = fabsf(act[i]);
            float imag = fabsf(act[i+1]);
            if (real > max_real) max_real = real;
            if (imag > max_imag) max_imag = imag;
        }
        float scale_real = max_real / 128.0f;
        float scale_imag = max_imag / 128.0f;
        lut_scale[blk*2  ] = scale_real;
        lut_scale[blk*2+1] = scale_imag;
        float inv_scale_real = 1.0f / scale_real;
        float inv_scale_imag = 1.0f / scale_imag;
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

            for (int _i = 0; _i < 4; _i++) {
                for (int _j = 0; _j < 4; _j++) {
                    for (int _k = 0; _k < 4; _k++) {
                        int ac =
                            r0 * sx[_i] +
                            r1 * sx[_j] +
                            r2 * sx[_k];

                        int bd =
                            i0 * sy[_i] +
                            i1 * sy[_j] +
                            i2 * sy[_k];

                        int ad =
                            r0 * sy[_i] +
                            r1 * sy[_j] +
                            r2 * sy[_k];

                        int bc =
                            i0 * sx[_i] +
                            i1 * sx[_j] +
                            i2 * sx[_k];

                        lut[ii++] = ac;
                        lut[ii++] = bd;
                        lut[ii++] = ad;
                        lut[ii++] = bc;
                    }
                }
            }
        }
    }
}


void transpose_old(int m, int k, const block_ifairy *raw_w, block_ifairy_1x3_old *w) {
    const int blk_n = (k+QK_K-1)/QK_K;
    for (int j = 0; j < m; j++) {
        for (int blk = 0; blk < blk_n; blk++) {
            int ii = 0;
            #pragma unroll
            for (int i = 0; i < QK_K/4; i+=3) {
                const uint32_t *p = (uint32_t *)&(raw_w[j*blk_n + blk].qs[i]);
                uint8_t iweight_1x3_0, iweight_1x3_1, iweight_1x3_2, iweight_1x3_3;
                if (i+3 >= QK_K/4) {
                    uint32_t iweight_1x12 = p[0] << 16;
                    iweight_1x3_0 = iweight_1x12 >> 18 & 0b00111111; 
                    iweight_1x3_1 = iweight_1x12 >> 12 & 0b00111111;
                    w[j*blk_n + blk].qs[ii  ] = iweight_1x3_0;
                    w[j*blk_n + blk].qs[ii+1] = iweight_1x3_1;
                } else {
                    uint32_t iweight_1x12 = __builtin_bswap32(*p) >> 8;
                    iweight_1x3_0 = iweight_1x12 >> 18 & 0b00111111; 
                    iweight_1x3_1 = iweight_1x12 >> 12 & 0b00111111;
                    iweight_1x3_2 = iweight_1x12 >> 6  & 0b00111111; 
                    iweight_1x3_3 = iweight_1x12       & 0b00111111;
                    w[j*blk_n + blk].qs[ii  ] = iweight_1x3_0;
                    w[j*blk_n + blk].qs[ii+1] = iweight_1x3_1;
                    w[j*blk_n + blk].qs[ii+2] = iweight_1x3_2;
                    w[j*blk_n + blk].qs[ii+3] = iweight_1x3_3;
                }
                ii += 4;
            }
            w[j*blk_n + blk].d_real = raw_w[j*blk_n + blk].d_real;
            w[j*blk_n + blk].d_imag = raw_w[j*blk_n + blk].d_imag;
        }
    }
}


void mul_mat_mxk_kx1_with_lut_old(int k, int row_begin, int row_end, const block_ifairy_1x3_old *w, const int16_t *lut, const float *lut_scale, float *dst) {
    const int blk_n = (k+QK_K-1)/QK_K;
    const int in_blk_n = (QK_K+2)/3;
    const int next_lut_blk = 256;
    const int lut_blk_n = 4; 
    for (int row = row_begin; row <= row_end; row++) {
        const int w_base = row*blk_n;
        for (int blk = 0; blk < blk_n; blk++) {
            int32x4_t tmp = vdupq_n_s32(0);
            const int16_t *lut_base = lut + (blk*in_blk_n)*next_lut_blk;
            #pragma unroll
            for (int i = 0; i < in_blk_n; i++) {
                const uint8_t iweight_16x3 = w[w_base+blk].qs[i];
                int16x4_t abcd = vld1_s16(lut_base+next_lut_blk*i+iweight_16x3*lut_blk_n);
                tmp = vaddw_s16(tmp, abcd);
            }
            // 反量化，写回
            dst[(row)*2  ] += (float)(tmp[0]) * lut_scale[blk*2] * w[(row)*blk_n+blk].d_real - (float)(tmp[1]) * lut_scale[blk*2+1] * w[(row)*blk_n+blk].d_imag;
            dst[(row)*2+1] += (float)(tmp[2]) * lut_scale[blk*2] * w[(row)*blk_n+blk].d_imag + (float)(tmp[3]) * lut_scale[blk*2+1] * w[(row)*blk_n+blk].d_real;
        }
    }
}


int16_t *alloc_lut_v_old(int k) {
    int blk_n = (k+QK_K-1)/QK_K;
    int in_blk_n = (QK_K+2)/3;
    return calloc(blk_n*in_blk_n*256, sizeof(int16_t));
}


float *alloc_lut_scale_old(int k) {
    int blk_n = (k+QK_K-1)/QK_K;
    return calloc(blk_n*2, sizeof(float));
}


block_ifairy_1x3_old *alloc_w_old(int m, int k) {
    int blk_n = (k+QK_K-1)/QK_K;
    int row_n = (m+15)/16;
    block_ifairy_1x3_old *ptr = calloc(blk_n*row_n*16, sizeof(block_ifairy_1x3_old));
    return ptr;
}


void free_lut_v_old(int16_t *lut) {
    free(lut);
}


void free_lut_scale_old(float *scale) {
    free(scale);
}


void free_w_old(block_ifairy_1x3_old *w) {
    free(w);
}