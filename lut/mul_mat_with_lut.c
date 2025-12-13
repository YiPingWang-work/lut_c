//
// Created by Yiping Wang on 2025/11/27.
//
#include "mul_mat_with_lut.h"

int8_t set_int8(int16_t v) {
    if (v > 127) return 127;
    if (v < -128) return -128;
    return (int8_t)v;
}

void transpose_4x16_to_16x4(int16x4_t src[16], int8x16_t dst[4]) {
    int16_t tmp16x4[16][4];
    for (int i = 0; i < 16; i++) {
        vst1_s16(tmp16x4[i], src[i]);
    }
    int8_t tmp4x16[4][16];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 16; j+=2) {
            tmp4x16[i][j>>1] = set_int8(tmp16x4[j][i]);
            tmp4x16[i][(j>>1)+8] = set_int8(tmp16x4[j+1][i]);
        }
    }
    for (int i = 0; i < 4; i++) {
        dst[i] = vld1q_s8(tmp4x16[i]);
    }
}

int16x4_t get_act(const int16_t *act, int i, int m) {
    m *= 2;
    int16_t ret[4] = {0, 0, 0, 0};
    if (i < m) ret[0] = act[i];
    if (i + 6 < m) ret[1] = act[i + 6];
    if (i + 12 < m) ret[2] = act[i + 12];
    if (i + 18 < m) ret[3] = act[i + 18];
    return vld1_s16(ret);
}


void generate_lut_int8(const int16_t *act, int m, int8x16x2_t *lut) {
    for (int i = 0; i < 2*m; i += 24) { // 每次处理12个复数
        int16x4_t r0 = get_act(act, i,   m);
        int16x4_t i0 = get_act(act, i+1, m);
        int16x4_t r1 = get_act(act, i+2, m);
        int16x4_t i1 = get_act(act, i+3, m);
        int16x4_t r2 = get_act(act, i+4, m);
        int16x4_t i2 = get_act(act, i+5, m);

        int16x4_t val0[16] = {
            // (-1, -1, -1)
            vsub_s16(vsub_s16(vneg_s16(r0), r1), r2),
            vsub_s16(vsub_s16(vneg_s16(i0), i1), i2),

            // (-1, -1,  1)
            vadd_s16(vsub_s16(vneg_s16(r0), r1), r2),
            vadd_s16(vsub_s16(vneg_s16(i0), i1), i2),

            // (-1, -1, -i)  → conj(-i) = +i
            vsub_s16(vsub_s16(vneg_s16(r0), r1), i2),
            vadd_s16(vsub_s16(vneg_s16(i0), i1), r2),

            // (-1, -1,  i)  → conj(i) = -i
            vadd_s16(vsub_s16(vneg_s16(r0), r1), i2),
            vsub_s16(vsub_s16(vneg_s16(i0), i1), r2),

            // (-1,  1, -1)
            vsub_s16(vadd_s16(vneg_s16(r0), r1), r2),
            vsub_s16(vadd_s16(vneg_s16(i0), i1), i2),

            // (-1,  1,  1)
            vadd_s16(vadd_s16(vneg_s16(r0), r1), r2),
            vadd_s16(vadd_s16(vneg_s16(i0), i1), i2),

            // (-1,  1, -i) → conj(-i) = +i
            vsub_s16(vadd_s16(vneg_s16(r0), r1), i2),
            vadd_s16(vadd_s16(vneg_s16(i0), i1), r2),

            // (-1,  1,  i) → conj(i) = -i
            vadd_s16(vadd_s16(vneg_s16(r0), r1), i2),
            vsub_s16(vadd_s16(vneg_s16(i0), i1), r2),
        };

        int16x4_t val1[16] = {
            // (-1, -i, -1) → conj(-i) = +i
            vsub_s16(vadd_s16(vneg_s16(r0), i1), r2),
            vadd_s16(vsub_s16(vneg_s16(i0), r1), i2),

            // (-1, -i,  1)
            vadd_s16(vadd_s16(vneg_s16(r0), i1), r2),
            vsub_s16(vadd_s16(vneg_s16(i0), r1), i2),

            // (-1, -i, -i) → conj(-i) = +i
            vsub_s16(vadd_s16(vneg_s16(r0), i1), i2),
            vadd_s16(vsub_s16(vneg_s16(i0), r1), r2),

            // (-1, -i,  i) → conj(i) = -i
            vadd_s16(vadd_s16(vneg_s16(r0), i1), i2),
            vsub_s16(vsub_s16(vneg_s16(i0), r1), r2),

            // (-1,  i, -1) → conj(i) = -i
            vsub_s16(vsub_s16(vneg_s16(r0), i1), r2),
            vadd_s16(vadd_s16(vneg_s16(i0), r1), i2),

            // (-1,  i,  1)
            vadd_s16(vsub_s16(vneg_s16(r0), i1), r2),
            vsub_s16(vadd_s16(vneg_s16(i0), r1), i2),

            // (-1,  i, -i) → conj(-i) = +i
            vsub_s16(vsub_s16(vneg_s16(r0), i1), i2),
            vadd_s16(vadd_s16(vneg_s16(i0), r1), r2),

            // (-1,  i,  i) → conj(i) = -i
            vadd_s16(vsub_s16(vneg_s16(r0), i1), i2),
            vsub_s16(vsub_s16(vneg_s16(i0), r1), r2),
        };

        int8x16_t t0[4], t1[4];
        transpose_4x16_to_16x4(val0, t0);
        transpose_4x16_to_16x4(val1, t1);

        lut[i/6    ].val[0] = t0[0]; lut[i/6    ].val[1] = t1[0];
        lut[i/6 + 1].val[0] = t0[1]; lut[i/6 + 1].val[1] = t1[1];
        lut[i/6 + 2].val[0] = t0[2]; lut[i/6 + 2].val[1] = t1[2];
        lut[i/6 + 3].val[0] = t0[3]; lut[i/6 + 3].val[1] = t1[3];
    }
    return;
}


int8x16x2_t mul_mat_block_16x3_3x1(uint8x16_t iweight_16x3, int8x16x2_t ilut) {
    // 生成index和flag
    uint8x16_t index = vandq_u8(iweight_16x3, vdupq_n_u8(7)); //  00 01 11
    uint8x16_t flag = vshrq_n_u8(vandq_u8(iweight_16x3, vdupq_n_u8(48)), 4); // 11 00 00

    // 查表获取
    int8x16_t r0 = vqtbl1q_s8(ilut.val[0], index); // 全部的16个数据在ilut0中的实部
    int8x16_t i0 = vqtbl1q_s8(ilut.val[0], vaddq_u8(index, vdupq_n_u8(8))); // 全部的16个数据在ilut0中的虚部
    int8x16_t r1 = vqtbl1q_s8(ilut.val[1], index); // 全部的16个数据在ilut1中的实部
    int8x16_t i1 = vqtbl1q_s8(ilut.val[1], vaddq_u8(index, vdupq_n_u8(8))); // 全部的16个数据在ilut1中的虚部

    // TODO 查表
    uint8x16_t mask = vceqq_u8(vandq_u8(iweight_16x3, vdupq_n_u8(8)), vdupq_n_u8(8)); // 表示该数据是否位于ilut[1]中
    int8x16_t iret_r_tmp0 = vbslq_u8(mask, r1, r0);
    int8x16_t iret_i_tmp0 = vbslq_u8(mask, i1, i0);

    // 如果是 10 or 11 => 先翻转
    mask = vcgeq_u8(flag, vdupq_n_u8(2));
    int8x16_t iret_r_tmp1 = vbslq_u8(mask, iret_i_tmp0, iret_r_tmp0);
    int8x16_t iret_i_tmp1 = vbslq_u8(mask, iret_r_tmp0, iret_i_tmp0);
    int8x16_t iret_r_tmp1_neg = vnegq_s8(iret_r_tmp1);
    int8x16_t iret_i_tmp1_neg = vnegq_s8(iret_i_tmp1);

    // 如果是 01 => ret_r_tmp1 ret_i_tmp1 取负数
    mask = vceqq_u8(flag, vdupq_n_u8(1));
    int8x16_t iret_r_tmp2 = vbslq_u8(mask, iret_r_tmp1_neg, iret_r_tmp1);
    int8x16_t iret_i_tmp2 = vbslq_u8(mask, iret_i_tmp1_neg , iret_i_tmp1);

    // 如果是是10 => ret_r_tmp1 取负数
    mask = vceqq_u8(flag, vdupq_n_u8(2));
    int8x16_t iret_r = vbslq_u8(mask, iret_r_tmp1_neg, iret_r_tmp2);

    // 如果是11 => ret_i_tmp1 取负数
    mask = vceqq_u8(flag, vdupq_n_u8(3));
    int8x16_t iret_i = vbslq_u8(mask, iret_i_tmp1_neg , iret_i_tmp2);

    return (int8x16x2_t){iret_r, iret_i};
}

uint8_t get_3x1(uint32_t iweght_1x12, int i) {
    return (iweght_1x12 >> ((3-i)*6)) & 0x3F;
}

void mul_mat_nxm_mx1_with_lut(block_ifairy *weight, int block_n, int row_begin, int row_end, int8x16x2_t *lut, int32_t *dst) {
    for (int row = row_begin; row <= row_end; row+=16) {
        for (int block = 0; block < block_n; block++) {
            for (int i = 0; i < QK_K; i+=3) {
                uint32_t iweight_1x12[16];
                if (i+3 >= QK_K) {
                    for (int j = 0; j < 16; j++) {
                        iweight_1x12[j] = (weight[(row+j)*block_n+block].qs[i] << 16) | (weight[(row+j)*block_n+block].qs[i+1] << 8);
                    }
                } else {
                    for (int j = 0; j < 16; j++) {
                        iweight_1x12[j] = (weight[(row+j)*block_n+block].qs[i] << 16) | (weight[(row+j)*block_n+block].qs[i+1] << 8) | weight[(row+j)*block_n+block].qs[i+2];
                    }
                }
                for (int ii = 0; ii < 4; ii++) {
                    uint8x16_t iweight_16x3 = {
                        get_3x1(iweight_1x12[0], ii),
                        get_3x1(iweight_1x12[1], ii),
                        get_3x1(iweight_1x12[2], ii),
                        get_3x1(iweight_1x12[3], ii),
                        get_3x1(iweight_1x12[4], ii),
                        get_3x1(iweight_1x12[5], ii),
                        get_3x1(iweight_1x12[6], ii),
                        get_3x1(iweight_1x12[7], ii),
                        get_3x1(iweight_1x12[8], ii),
                        get_3x1(iweight_1x12[9], ii),
                        get_3x1(iweight_1x12[10], ii),
                        get_3x1(iweight_1x12[11], ii),
                        get_3x1(iweight_1x12[12], ii),
                        get_3x1(iweight_1x12[13], ii),
                        get_3x1(iweight_1x12[14], ii),
                        get_3x1(iweight_1x12[15], ii),
                    };
                    int8x16x2_t iret_ri = mul_mat_block_16x3_3x1(iweight_16x3, lut[block*QK_K+i*4+ii]);
                    int8x16_t iret_r = iret_ri.val[0];
                    int8x16_t iret_i = iret_ri.val[1];
                    // TODO 反量化

                    // 回送结果
                    for (int j = 0; j < 16; j++) {
                        dst[2*(row+j)] += (int32_t)iret_r[j];
                        dst[2*(row+j)+1] += (int32_t)iret_i[j];
                    }
                }
            }
        }
    }
}

int8x16x2_t *alloc_lut(int m) {
    return calloc((m+11)/3, sizeof(int8x16x2_t));
}

void free_lut(int8x16x2_t *lut) {
    free(lut);
}

int32x2_t mul_mat_1x4_4x1(uint8_t a, int16_t *b) {
    int32_t acc_r = 0;
    int32_t acc_i = 0;

    for (int k = 0; k < 4; k++) {
        uint8_t code = (a >> (2 * k)) & 0x3;

        int16_t xr = b[2 * k];
        int16_t xi = b[2 * k + 1];

        switch (code) {
        case 0b00:  // -1
            acc_r += -xr;
            acc_i += -xi;
            break;

        case 0b01:  // +1
            acc_r += xr;
            acc_i += xi;
            break;

        case 0b10:  // -i -> conj = +i
            acc_r += -xi;
            acc_i +=  xr;
            break;

        case 0b11:  // +i -> conj = -i
            acc_r +=  xi;
            acc_i += -xr;
            break;
        }
    }

    return vset_lane_s32(acc_i,
           vset_lane_s32(acc_r, vdup_n_s32(0), 0),
           1);
}

void mul_mat_nxm_mx1(block_ifairy *weight, int block_n, int row_begin, int row_end, int16_t *act, int32_t *dst) {
    for (int row = row_begin; row <= row_end; row++) {
        for (int block = 0; block < block_n; block++) {
            for (int i = 0; i < QK_K/4; i++) {
                uint8_t w4 = weight[row*block_n+block].qs[i];
                int32x2_t ret = mul_mat_1x4_4x1(w4, act[2*row]);
                dst[2*row] += (int32_t)ret[0];
                dst[2*row+1] += (int32_t)ret[1];         
            }
        }
    }
}