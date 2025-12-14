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
        for (int j = 0; j < 16; j++) {
            tmp4x16[i][j] = set_int8(tmp16x4[j][i]);
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

void generate_lut_int8(const int16_t *act, int m, int8x16x2_t *lut, float32x2_t *lut_scale) {
    for (int i = 0; i < 2*m; i += 24) { // 每次处理12个复数
        int16x4_t r0 = get_act(act, i,   m);
        int16x4_t i0 = get_act(act, i+1, m);
        int16x4_t r1 = get_act(act, i+2, m);
        int16x4_t i1 = get_act(act, i+3, m);
        int16x4_t r2 = get_act(act, i+4, m);
        int16x4_t i2 = get_act(act, i+5, m);

        int16x4_t real[16] = {
            // 0 (-1, -1, -1)
            vsub_s16(vsub_s16(vneg_s16(r0), r1), r2),
            // 1 (-1, -1,  1)
            vadd_s16(vsub_s16(vneg_s16(r0), r1), r2),
            // 2 (-1, -1, -i)
            vsub_s16(vsub_s16(vneg_s16(r0), r1), i2),
            // 3 (-1, -1,  i) 
            vadd_s16(vsub_s16(vneg_s16(r0), r1), i2),
            // 4 (-1,  1, -1)
            vsub_s16(vadd_s16(vneg_s16(r0), r1), r2),
            // 5 (-1,  1,  1)
            vadd_s16(vadd_s16(vneg_s16(r0), r1), r2),
            // 6 (-1,  1, -i)
            vsub_s16(vadd_s16(vneg_s16(r0), r1), i2),
            // 7 (-1,  1,  i)
            vadd_s16(vadd_s16(vneg_s16(r0), r1), i2),
            // 8 (-1, -i, -1)
            vsub_s16(vsub_s16(vneg_s16(r0), i1), r2),
            // 9 (-1, -i,  1)
            vadd_s16(vsub_s16(vneg_s16(r0), i1), r2),
            // 10 (-1, -i, -i)
            vsub_s16(vsub_s16(vneg_s16(r0), i1), i2),
            // 11 (-1, -i,  i)
            vadd_s16(vsub_s16(vneg_s16(r0), i1), i2),
            // 12 (-1,  i, -1) conj(i)=-i → real = +i1
            vsub_s16(vadd_s16(vneg_s16(r0), i1), r2),
            // 13 (-1,  i,  1)
            vadd_s16(vadd_s16(vneg_s16(r0), i1), r2),
            // 14 (-1,  i, -i)
            vsub_s16(vadd_s16(vneg_s16(r0), i1), i2),
            // 15 (-1,  i,  i)
            vadd_s16(vadd_s16(vneg_s16(r0), i1), i2),
        };

        int16x4_t imag[16] = {
            // 0 (-1, -1, -1)
            vsub_s16(vsub_s16(vneg_s16(i0), i1), i2),
            // 1 (-1, -1,  1)
            vadd_s16(vsub_s16(vneg_s16(i0), i1), i2),
            // 2 (-1, -1, -i)
            vadd_s16(vsub_s16(vneg_s16(i0), i1), r2),
            // 3 (-1, -1,  i)
            vsub_s16(vsub_s16(vneg_s16(i0), i1), r2),
            // 4 (-1,  1, -1)
            vsub_s16(vadd_s16(vneg_s16(i0), i1), i2),
            // 5 (-1,  1,  1)
            vadd_s16(vadd_s16(vneg_s16(i0), i1), i2),
            // 6 (-1,  1, -i)
            vadd_s16(vadd_s16(vneg_s16(i0), i1), r2),
            // 7 (-1,  1,  i)
            vsub_s16(vadd_s16(vneg_s16(i0), i1), r2),
            // 8 (-1, -i, -1)
            vsub_s16(vadd_s16(vneg_s16(i0), r1), i2),
            // 9 (-1, -i,  1)
            vadd_s16(vadd_s16(vneg_s16(i0), r1), i2),
            // 10 (-1, -i, -i)
            vadd_s16(vadd_s16(vneg_s16(i0), r1), r2),
            // 11 (-1, -i,  i)
            vsub_s16(vadd_s16(vneg_s16(i0), r1), r2),
            // 12 (-1,  i, -1)
            vsub_s16(vsub_s16(vneg_s16(i0), r1), i2),
            // 13 (-1,  i,  1)
            vadd_s16(vsub_s16(vneg_s16(i0), r1), i2),
            // 14 (-1,  i, -i)
            vadd_s16(vsub_s16(vneg_s16(i0), r1), r2),
            // 15 (-1,  i,  i)
            vsub_s16(vsub_s16(vneg_s16(i0), r1), r2),
        };

        int8x16_t t0[4], t1[4];
        transpose_4x16_to_16x4(real, t0);
        transpose_4x16_to_16x4(imag, t1);

        lut[i/6    ].val[0] = t0[0]; lut[i/6    ].val[1] = t1[0];
        lut[i/6 + 1].val[0] = t0[1]; lut[i/6 + 1].val[1] = t1[1];
        lut[i/6 + 2].val[0] = t0[2]; lut[i/6 + 2].val[1] = t1[2];
        lut[i/6 + 3].val[0] = t0[3]; lut[i/6 + 3].val[1] = t1[3];
    }
    return;
}

int8x16x2_t mul_mat_block_16x3_3x1_with_lut(uint8x16_t iweight_16x3, int8x16x2_t ilut) {
    int8x16x2_t ilut_exp[4] = {
        {.val = {ilut.val[0], ilut.val[1]}},
        {.val = {vnegq_s8(ilut.val[0]), vnegq_s8(ilut.val[1])}},
        {.val = {vnegq_s8(ilut.val[1]), ilut.val[0]}},
        {.val = {ilut.val[1], vnegq_s8(ilut.val[0])}},
    };

    // 生成flag和index
    uint8x16_t lut_flag = vshrq_n_u8(vandq_u8(iweight_16x3, vdupq_n_u8(0b00110000)), 4); // 11 00 00
    
    // index 适配
    // 如果lut_flag = 00 保持不变
    uint8x16_t index = iweight_16x3;
    // 如果lut_flag = 01 后6位中全部的奇数位01交换
    index = veorq_u8(index, vandq_u8(vceqq_u8(lut_flag, vdupq_n_u8(0b01)) , vdupq_n_u8(0b00010101)));
    // 如果lut_flag = 10 后6位中全部的偶数位01交换
    index = veorq_u8(index, vandq_u8(vceqq_u8(lut_flag, vdupq_n_u8(0b10)) , vdupq_n_u8(0b00101010)));
    // 如果lut_flag = 11 后六位全部01交换
    index = veorq_u8(index, vandq_u8(vceqq_u8(lut_flag, vdupq_n_u8(0b11)) , vdupq_n_u8(0b00111111)));

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

uint8_t get_3x1(uint32_t iweght_1x12, int i) {
    return (iweght_1x12 >> ((3-i)*6)) & 0x3F;
}

void mul_mat_nxm_mx1_with_lut(block_ifairy *weight, int block_n, int row_begin, int row_end, int8x16x2_t *lut, float32x2_t *lut_scale, int32_t *dst) {
    for (int row = row_begin; row <= row_end; row+=16) {
        for (int block = 0; block < block_n; block++) {
            for (int i = 0; i < QK_K/4; i+=3) {
                if (i+3 >= QK_K/4) {
                    // 前3个复数可以做一次查表处理，TODO 但是最后一个复数发生了舍弃
                    uint8x16_t iweight_16x3 = {
                        weight[(row+0)*block_n+block].qs[i] >> 2,
                        weight[(row+1)*block_n+block].qs[i] >> 2,
                        weight[(row+2)*block_n+block].qs[i] >> 2,
                        weight[(row+3)*block_n+block].qs[i] >> 2,
                        weight[(row+4)*block_n+block].qs[i] >> 2,
                        weight[(row+5)*block_n+block].qs[i] >> 2,
                        weight[(row+6)*block_n+block].qs[i] >> 2,
                        weight[(row+7)*block_n+block].qs[i] >> 2,
                        weight[(row+8)*block_n+block].qs[i] >> 2,
                        weight[(row+9)*block_n+block].qs[i] >> 2,
                        weight[(row+10)*block_n+block].qs[i] >> 2,
                        weight[(row+11)*block_n+block].qs[i] >> 2,
                        weight[(row+12)*block_n+block].qs[i] >> 2,
                        weight[(row+13)*block_n+block].qs[i] >> 2,
                        weight[(row+14)*block_n+block].qs[i] >> 2,
                        weight[(row+15)*block_n+block].qs[i] >> 2,
                    };
                    int lut_i = (block*QK_K+4*i)/3;
                    int8x16x2_t iret_ri = mul_mat_block_16x3_3x1_with_lut(iweight_16x3, lut[lut_i]);
                    int8x16_t iret_r = iret_ri.val[0];
                    int8x16_t iret_i = iret_ri.val[1];
                    // TODO 反量化
                    // 回送结果
                    for (int j = 0; j < 16; j++) {
                        dst[2*(row+j)] += (int32_t)iret_r[j];
                        dst[2*(row+j)+1] += (int32_t)iret_i[j];
                    }
                } else {
                    uint32_t iweight_1x12[16];
                    for (int j = 0; j < 16; j++) {
                        iweight_1x12[j] = (weight[(row+j)*block_n+block].qs[i] << 16) | (weight[(row+j)*block_n+block].qs[i+1] << 8) | weight[(row+j)*block_n+block].qs[i+2];
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
                        int lut_i = (block*QK_K+4*i)/3+ii;
                        int8x16x2_t iret_ri = mul_mat_block_16x3_3x1_with_lut(iweight_16x3, lut[lut_i]);
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
    for (int row = row_begin; row <= row_end; row++) {
        if (row < 10) printf("%d:   %d, %d\n", row, dst[2*row], dst[2*row+1]);
        if (row > 1000) printf("%d:   %d, %d\n", row, dst[2*row], dst[2*row+1]);
    }
}

int8x16x2_t *alloc_lut(int m) {
    return calloc((m+11)/3, sizeof(int8x16x2_t));
}

void free_lut(int8x16x2_t *lut) {
    free(lut);
}

int32x2_t mul_mat_block_1x4_4x1(uint8_t a, int16_t *b) {
    int32_t acc_r = 0;
    int32_t acc_i = 0;

    for (int k = 0; k < 4; k++) {
        uint8_t code = (a >> (2 * (3-k))) & 0x3;

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
                int16_t a4[8] = {
                    act[block*QK_K+i*8],
                    act[block*QK_K+i*8+1],
                    act[block*QK_K+i*8+2],
                    act[block*QK_K+i*8+3],
                    act[block*QK_K+i*8+4],
                    act[block*QK_K+i*8+5],
                    act[block*QK_K+i*8+6],
                    act[block*QK_K+i*8+7],
                };
                int32x2_t ret = mul_mat_block_1x4_4x1(w4, a4);
                dst[2*row] += (int32_t)ret[0];
                dst[2*row+1] += (int32_t)ret[1];         
            }
        }
        if (row < 10) printf("%d:   %d, %d\n", row, dst[2*row], dst[2*row+1]);
        if (row > 1000) printf("%d:   %d, %d\n", row, dst[2*row], dst[2*row+1]);
    }
}