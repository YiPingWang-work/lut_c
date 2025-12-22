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

void generate_lut_int8(const int16_t *act, int m, lut_block *lut) { // m个复数
    int block_n = (m+QK_K-1)/QK_K;
    for (int block = 0; block < block_n; block++) {
        int act_begin = block * QK_K * 2, act_end = ((block+1)*QK_K-1)*2 < m*2 ? (block+1)*QK_K*2-1 : m*2;
        for (int i = act_begin; i <= act_end; i += 6) {
            int16_t r0, r1, r2, i0, i1, i2;
            if (i+6 > act_end) {
                r0 = act[i];
                i0 = -act[i+1];
                r1 = 0;
                i1 = 0;
                r2 = 0;
                i2 = 0;
            } else {
                r0 = act[i];
                i0 = -act[i+1];
                r1 = act[i+2];
                i1 = -act[i+3];
                r2 = act[i+4];
                i2 = -act[i+5];
            }
            int8x16_t real, imag;
            real[0]  = set_int8(-r0 - r1 - r2);   imag[0]  = set_int8(-i0 - i1 - i2);    // (-1, -1, -1)
            real[1]  = set_int8(-r0 - r1 + r2);   imag[1]  = set_int8(-i0 - i1 + i2);    // (-1, -1,  1)
            real[2]  = set_int8(-r0 - r1 + i2);   imag[2]  = set_int8(-i0 - i1 - r2);    // (-1, -1, -i)
            real[3]  = set_int8(-r0 - r1 - i2);   imag[3]  = set_int8(-i0 - i1 + r2);    // (-1, -1,  i)

            real[4]  = set_int8(-r0 + r1 - r2);   imag[4]  = set_int8(-i0 + i1 - i2);    // (-1,  1, -1)
            real[5]  = set_int8(-r0 + r1 + r2);   imag[5]  = set_int8(-i0 + i1 + i2);    // (-1,  1,  1)
            real[6]  = set_int8(-r0 + r1 + i2);   imag[6]  = set_int8(-i0 + i1 - r2);    // (-1,  1, -i)
            real[7]  = set_int8(-r0 + r1 - i2);   imag[7]  = set_int8(-i0 + i1 + r2);    // (-1,  1,  i)

            real[8]  = set_int8(-r0 + i1 - r2);   imag[8]  = set_int8(-i0 - r1 - i2);    // (-1, -i, -1)
            real[9]  = set_int8(-r0 + i1 + r2);   imag[9]  = set_int8(-i0 - r1 + i2);    // (-1, -i,  1)
            real[10] = set_int8(-r0 + i1 + i2);   imag[10] = set_int8(-i0 - r1 - r2);    // (-1, -i, -i)
            real[11] = set_int8(-r0 + i1 - i2);   imag[11] = set_int8(-i0 - r1 + r2);    // (-1, -i,  i)

            real[12] = set_int8(-r0 - i1 - r2);   imag[12] = set_int8(-i0 + r1 - i2);    // (-1,  i, -1)
            real[13] = set_int8(-r0 - i1 + r2);   imag[13] = set_int8(-i0 + r1 + i2);    // (-1,  i,  1)
            real[14] = set_int8(-r0 - i1 + i2);   imag[14] = set_int8(-i0 + r1 - r2);    // (-1,  i, -i)
            real[15] = set_int8(-r0 - i1 - i2);   imag[15] = set_int8(-i0 + r1 + r2);    // (-1,  i,  i)
    
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

static inline uint8_t get_3x1(uint32_t iweght_1x12, int i) __attribute__((always_inline));
static inline uint8_t get_3x1(uint32_t iweght_1x12, int i) {
    return (iweght_1x12 >> ((3-i)*6)) & 0x3F;
}

void mul_mat_nxm_mx1_with_lut(block_ifairy *weight, int block_n, int row_begin, int row_end, lut_block *lut, float32x2_t *lut_scale, int32_t *dst) {
    for (int row = row_begin; row <= row_end; row+=16) {
        for (int block = 0; block < block_n; block++) {
            for (int i = 0; i < QK_K/4; i+=3) {
                uint32_t iweight_1x12[16];
                int max_ii;
                if (i+3 >= QK_K/4) {
                    max_ii = 2;
                    for (int j = 0; j < 16; j++) {
                        iweight_1x12[j] = (weight[(row+j)*block_n+block].qs[i] << 16);
                    }
                } else {
                    max_ii = 4;
                    for (int j = 0; j < 16; j++) {
                        iweight_1x12[j] = (weight[(row+j)*block_n+block].qs[i] << 16) | (weight[(row+j)*block_n+block].qs[i+1] << 8) | weight[(row+j)*block_n+block].qs[i+2];
                    }
                }
                for (int ii = 0; ii < max_ii; ii++) {
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
                    int8x16x2_t iret_ri = mul_mat_block_16x3_3x1_with_lut(iweight_16x3, lut[block].v[4*i/3+ii]);
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

lut_block *alloc_lut(int m) {
    int block_n = (m+QK_K-1)/QK_K;
    return calloc(m, sizeof(lut_block));
}

void free_lut(lut_block *lut) {
    free(lut);
}



// ================================ 验证乘法程序 ================================
static inline int32x2_t mul_mat_block_1x4_4x1(uint8_t a, int16_t *b) __attribute__((always_inline));
static inline int32x2_t mul_mat_block_1x4_4x1(uint8_t a, int16_t *b) {
    int32_t acc_r = 0;
    int32_t acc_i = 0;

    for (int k = 0; k < 4; k++) {
        uint8_t code = (a >> (2 * (3-k))) & 0x3;

        int16_t xr = b[2 * k];
        int16_t xi = -b[2 * k + 1]; // 共轭

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
                    act[block*QK_K*2+i*8],
                    act[block*QK_K*2+i*8+1],
                    act[block*QK_K*2+i*8+2],
                    act[block*QK_K*2+i*8+3],
                    act[block*QK_K*2+i*8+4],
                    act[block*QK_K*2+i*8+5],
                    act[block*QK_K*2+i*8+6],
                    act[block*QK_K*2+i*8+7],
                };
                int32x2_t ret = mul_mat_block_1x4_4x1(w4, a4);
                dst[2*row] += (int32_t)ret[0];
                dst[2*row+1] += (int32_t)ret[1];         
            }
        }
    }
}