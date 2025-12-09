//
// Created by Yiping Wang on 2025/11/27.
//
#include <printf.h>
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

int16x4_t get_activation(const int16_t *activation, int i, int m) {
    int16_t ret[4] = {0, 0, 0, 0};
    if (i < m) ret[0] = activation[i];
    if (i + 6 < m) ret[1] = activation[i + 6];
    if (i + 12 < m) ret[2] = activation[i + 12];
    if (i + 18 < m) ret[3] = activation[i + 18];
    return vld1_s16(ret);
}


int generate_lut_int8(const int16_t *activation, int m, int8x16x2_t *lut) {
    int max_lut_m = (m+5)/6;
    for (int i = 0; i < m; i += 24) {
        int16x4_t r0 = get_activation(activation, i, m);
        int16x4_t i0 = get_activation(activation, i+1, m);
        int16x4_t r1 = get_activation(activation, i+2, m);
        int16x4_t i1 = get_activation(activation, i+3, m);
        int16x4_t r2 = get_activation(activation, i+4, m);
        int16x4_t i2 = get_activation(activation, i+5, m);


        int16x4_t val0[16] = {
                // (-1, -1, -1)
                vsub_s16(vsub_s16(vneg_s16(r0), r1), r2),
                vsub_s16(vsub_s16(vneg_s16(i0), i1), i2),
                // (-1, -1, 1)
                vadd_s16(vsub_s16(vneg_s16(r0), r1), r2),
                vadd_s16(vsub_s16(vneg_s16(i0), i1), i2),
                // (-1, -1, -i)
                vadd_s16(vsub_s16(vneg_s16(r0), r1), i2),
                vsub_s16(vsub_s16(vneg_s16(i0), i1), r2),
                // (-1, -1, i)
                vsub_s16(vsub_s16(vneg_s16(r0), r1), i2),
                vadd_s16(vsub_s16(vneg_s16(i0), i1), r2),
                // (-1, 1, -1)
                vsub_s16(vadd_s16(vneg_s16(r0), r1), r2),
                vsub_s16(vadd_s16(vneg_s16(i0), i1), i2),
                // (-1, 1, 1)
                vadd_s16(vadd_s16(vneg_s16(r0), r1), r2),
                vadd_s16(vadd_s16(vneg_s16(i0), i1), i2),
                // (-1, 1, -i)
                vadd_s16(vadd_s16(vneg_s16(r0), r1), i2),
                vsub_s16(vadd_s16(vneg_s16(i0), i1), r2),
                // (-1, 1, i)
                vsub_s16(vadd_s16(vneg_s16(r0), r1), i2),
                vadd_s16(vadd_s16(vneg_s16(i0), i1), r2)
        };

        int16x4_t val1[16] = {
                // (-1, -i, -1)
                vsub_s16(vadd_s16(vneg_s16(r0), i1), r2),
                vsub_s16(vsub_s16(vneg_s16(i0), r1), i2),
                // (-1, -i, 1)
                vadd_s16(vadd_s16(vneg_s16(r0), i1), r2),
                vadd_s16(vsub_s16(vneg_s16(i0), r1), i2),
                // (-1, -i, -i)
                vsub_s16(vadd_s16(vneg_s16(r0), i1), i2),
                vsub_s16(vsub_s16(vneg_s16(i0), r1), r2),
                // (-1, -i, i)
                vadd_s16(vadd_s16(vneg_s16(r0), i1), i2),
                vadd_s16(vsub_s16(vneg_s16(i0), r1), r2),
                // (-1, i, -1)
                vsub_s16(vsub_s16(vneg_s16(r0), i1), r2),
                vsub_s16(vadd_s16(vneg_s16(i0), r1), i2),
                // (-1, i, 1)
                vadd_s16(vsub_s16(vneg_s16(r0), i1), r2),
                vadd_s16(vadd_s16(vneg_s16(i0), r1), i2),
                // (-1, i, -i)
                vsub_s16(vsub_s16(vneg_s16(r0), i1), i2),
                vsub_s16(vadd_s16(vneg_s16(i0), r1), r2),
                // (-1, i, i)
                vadd_s16(vsub_s16(vneg_s16(r0), i1), i2),
                vadd_s16(vadd_s16(vneg_s16(i0), r1), r2)
        };

        int8x16_t t0[4], t1[4];
        transpose_4x16_to_16x4(val0, t0);
        transpose_4x16_to_16x4(val1, t1);
        if (i/6 < max_lut_m) lut[i/6].val[0] = t0[0], lut[i/6].val[1] = t1[0];
        if (i/6+1 < max_lut_m) lut[i/6+1].val[0] = t0[1], lut[i/6+1].val[1] = t1[1];
        if (i/6+2 < max_lut_m) lut[i/6+2].val[0] = t0[2], lut[i/6+2].val[1] = t1[2];
        if (i/6+3 < max_lut_m) lut[i/6+3].val[0] = t0[3], lut[i/6+3].val[1] = t1[3];
    }
    return max_lut_m;
}

uint8_t decode_weight(void * weight, int row, int j) { // TODO
    return 0;
}

void* encode_result() { // TODO
    return NULL;
}

uint8x16_t get_flag(uint8x16_t v) {
    uint8x16_t mask = vdupq_n_u8(0b00110000);
    return vshrq_n_u8(vandq_u8(v, mask), 4);
}

uint8x16_t get_index(uint8x16_t v) {
    uint8x16_t mask = vdupq_n_u8(0b00000111);
    return vandq_u8(v, mask);
}

uint8x16_t use_tlb1(uint8x16_t v) {
    uint8x16_t mask = vdupq_n_u8(0b00001000);
    return vshrq_n_u8(vandq_u8(v, mask), 3);
}

void mul_mat_nxm_mx1(void *weight, int n, int m, int8x16x2_t *lut, void *output) {
    for (int i = 0; i < n; i += 3) { // 每次处理3个复数
        for (int j = 0; j < m; j += 16) { // 16行并行处理

            uint8_t weight_j_i_3x16_arr[16] = { // 获取权重
                decode_weight(weight, j, i),
                decode_weight(weight, j+1, i),
                decode_weight(weight, j+2, i),
                decode_weight(weight, j+3, i),
                decode_weight(weight, j+4, i),
                decode_weight(weight, j+5, i),
                decode_weight(weight, j+6, i),
                decode_weight(weight, j+7, i),
                decode_weight(weight, j+8, i),
                decode_weight(weight, j+9, i),
                decode_weight(weight, j+10, i),
                decode_weight(weight, j+11, i),
                decode_weight(weight, j+12, i),
                decode_weight(weight, j+13, i),
                decode_weight(weight, j+14, i),
                decode_weight(weight, j+15, i),
            };

            uint8x16_t weight_j_i_3x16 = vld1q_u8(weight_j_i_3x16_arr);

            // 生成index和flag
            uint8x16_t index = get_index(weight_j_i_3x16);
            uint8x16_t flag = get_flag(weight_j_i_3x16);


            uint8x16_t tlb1 = use_tlb1(weight_j_i_3x16);

            // 查表获取
            int8x16_t r0 = vqtbl1q_s8(lut[i].val[0], index);
            int8x16_t i0 = vqtbl1q_s8(lut[i].val[1], vshlq_n_u8(index, 1));
            int8x16_t r1 = vqtbl1q_s8(lut[i].val[0], index);
            int8x16_t i1 = vqtbl1q_s8(lut[i].val[1], vshlq_n_u8(index, 1));

            // TODO 查表

            int8x16_t ret_r_tmp0 = {};
            int8x16_t ret_i_tmp0 = {};

            // 如果是 10 or 11 => 先翻转
            uint8x16_t mask = vcgeq_u8(flag, vdupq_n_u8(2));
            int8x16_t ret_r_tmp1 = vbslq_u8(mask, ret_i_tmp0, ret_r_tmp0);
            int8x16_t ret_i_tmp1 = vbslq_u8(mask, ret_r_tmp0, ret_i_tmp0);
            int8x16_t ret_r_tmp1_neg = vnegq_s8(ret_r_tmp1);
            int8x16_t ret_i_tmp1_neg = vnegq_s8(ret_i_tmp1);

            // 如果是 01 => ret_r_tmp1 ret_i_tmp1 取负数
            mask = vceqq_u8(flag, vdupq_n_u8(1));
            int8x16_t ret_r_tmp2 = vbslq_u8(mask, ret_r_tmp1_neg, ret_r_tmp1);
            int8x16_t ret_i_tmp2 = vbslq_u8(mask, ret_i_tmp1_neg , ret_i_tmp1);

            // 如果是是10 => ret_r_tmp1 取负数
            mask = vceqq_u8(flag, vdupq_n_u8(2));
            int8x16_t ret_r = vbslq_u8(mask, ret_r_tmp1_neg, ret_r_tmp2);

            // 如果是11 => ret_i_tmp1 取负数
            mask = vceqq_u8(flag, vdupq_n_u8(3));
            int8x16_t ret_i = vbslq_u8(mask, ret_i_tmp1_neg , ret_i_tmp2);



        }
    }
}
