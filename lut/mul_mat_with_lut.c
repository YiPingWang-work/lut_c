//
// Created by Yiping Wang on 2025/11/27.
//
#include "mul_mat_with_lut.h"


void act_float_2_block_ifairy_q16(int k, const float *act_float, block_ifairy_q16 *act_q16, float scale) {
    const int _blk_n = (k+QK_K-1)/QK_K;
    for (int _blk = 0; _blk < _blk_n; _blk++) {
        int _act_begin = _blk * QK_K * 2, _act_end = ((_blk+1)*QK_K-1)*2 < k*2 ? (_blk+1)*QK_K*2-1 : k*2;
        float _max_real = 0.0f, _max_imag = 0.0f;
        // 计算缩放因子
        for (int _i = _act_begin; _i <= _act_end; _i += 2) {
            float _real = fabsf(act_float[_i  ]);
            float _imag = fabsf(act_float[_i+1]);
            if (_real > _max_real) _max_real = _real;
            if (_imag > _max_imag) _max_imag = _imag;
        }
        float _scale_real = _max_real / scale;
        float _scale_imag = _max_imag / scale;

        act_q16[_blk].d_real = _scale_real;
        act_q16[_blk].d_imag = _scale_imag;
        float _inv_scale_real = 1.0f / _scale_real;
        float _inv_scale_imag = 1.0f / _scale_imag;
        for (int _i = _act_begin; _i <= _act_end; _i += 2) {
            act_q16[_blk].x_real[(_i-_act_begin)/2] = (uint8_t)(roundf(act_float[_i  ] * _inv_scale_real));
            act_q16[_blk].x_imag[(_i-_act_begin)/2] = (uint8_t)(roundf(act_float[_i+1] * _inv_scale_imag));
        }
    }
}



// ================================ 16路LUT ================================
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


void generate_lut_q8(int k, const float *act, lut_block *lut) {
    const int _blk_n = (k+QK_K-1)/QK_K;
    for (int _blk = 0; _blk < _blk_n; _blk++) {
        int _act_begin = _blk * QK_K * 2, _act_end = ((_blk+1)*QK_K-1)*2 < k*2 ? (_blk+1)*QK_K*2-1 : k*2;
        float _max_real = 0.0f, _max_imag = 0.0f;
        // 计算缩放因子
        for (int _i = _act_begin; _i <= _act_end; _i += 2) {
            float _real = fabsf(act[_i  ]);
            float _imag = fabsf(act[_i+1]);
            if (_real > _max_real) _max_real = _real;
            if (_imag > _max_imag) _max_imag = _imag;
        }
        float _scale_real = _max_real / 42.6f;
        float _scale_imag = _max_imag / 42.6f;
        lut[_blk].d_real = _scale_real;
        lut[_blk].d_imag = _scale_imag;
        float _inv_scale_real = 1.0f / _scale_real;
        float _inv_scale_imag = 1.0f / _scale_imag;
        
        // 生成lut表
        for (int _i = _act_begin; _i <= _act_end; _i += 6) {
            int8_t _r0, _r1, _r2, _i0, _i1, _i2;
            if (_i+6 > _act_end) {
                _r0 = (int8_t)roundf(act[_i]    * _inv_scale_real);
                _i0 = (int8_t)roundf(-act[_i+1] * _inv_scale_imag);
                _r1 = 0;
                _i1 = 0;
                _r2 = 0;
                _i2 = 0;
            } else {
                _r0 = (int8_t)roundf(act[_i]    * _inv_scale_real);
                _i0 = (int8_t)roundf(-act[_i+1] * _inv_scale_imag);
                _r1 = (int8_t)roundf(act[_i+2]  * _inv_scale_real);
                _i1 = (int8_t)roundf(-act[_i+3] * _inv_scale_imag);
                _r2 = (int8_t)roundf(act[_i+4]  * _inv_scale_real);
                _i2 = (int8_t)roundf(-act[_i+5] * _inv_scale_imag);
            }

            int8x16_t _ac = {
                -_r0 - _r1 - _r2,
                -_r0 - _r1 + _r2,
                -_r0 - _r1,
                -_r0 - _r1,

                -_r0 + _r1 - _r2,
                -_r0 + _r1 + _r2,
                -_r0 + _r1,
                -_r0 + _r1,

                -_r0 - _r2,
                -_r0 + _r2,
                -_r0,
                -_r0,

                -_r0 - _r2,
                -_r0 + _r2,
                -_r0,
                -_r0,
            };

            int8x16_t _bd = {
                0,
                0,
                -_i2,
                _i2,

                0,
                0,
                -_i2,
                _i2,

                -_i1,
                -_i1,
                -_i1 - _i2,
                -_i1 + _i2,

                _i1,
                _i1,
                _i1 - _i2,
                _i1 + _i2,
            };

            int8x16_t _ad = {
                0,
                0,
                -_r2,
                _r2,

                0,
                0,
                -_r2,
                _r2,

                -_r1,
                -_r1,
                -_r1 - _r2,
                -_r1 + _r2,

                _r1,
                _r1,
                _r1 - _r2,
                _r1 + _r2,
            };
            
            int8x16_t _bc = {
                -_i0 - _i1 - _i2,
                -_i0 - _i1 + _i2,
                -_i0 - _i1,
                -_i0 - _i1,

                -_i0 + _i1 - _i2,
                -_i0 + _i1 + _i2,
                -_i0 + _i1,
                -_i0 + _i1,

                -_i0 - _i2,
                -_i0 + _i2,
                -_i0,
                -_i0,

                -_i0 - _i2,
                -_i0 + _i2,
                -_i0,
                -_i0,
            };

            lut[_blk].v[(_i-_act_begin)/6] = (int8x16x4_t){.val = _ac, _bd, _ad, _bc};
        }

    }
}


void generate_lut_q8_block_ifairy_q16(int k, const block_ifairy_q16 *act, lut_block *lut) {
    const int _blk_n = (k+QK_K-1)/QK_K;
    for (int _blk = 0; _blk < _blk_n; _blk++) {
        lut[_blk].d_real = act[_blk].d_real;
        lut[_blk].d_imag = act[_blk].d_imag;
        int8_t _r0, _i0, _r1, _i1, _r2, _i2;
        int _i = 0;
        for (; _i+3 < QK_K; _i+=3) {
            
            _r0 = (int8_t) act[_blk].x_real[_i  ];
            _i0 = (int8_t)-act[_blk].x_imag[_i  ];
            _r1 = (int8_t) act[_blk].x_real[_i+1];
            _i1 = (int8_t)-act[_blk].x_imag[_i+1];
            _r2 = (int8_t) act[_blk].x_real[_i+2];
            _i2 = (int8_t)-act[_blk].x_imag[_i+2];
            

            int8x16_t _ac = {
                -_r0 - _r1 - _r2,
                -_r0 - _r1 + _r2,
                -_r0 - _r1,
                -_r0 - _r1,

                -_r0 + _r1 - _r2,
                -_r0 + _r1 + _r2,
                -_r0 + _r1,
                -_r0 + _r1,

                -_r0 - _r2,
                -_r0 + _r2,
                -_r0,
                -_r0,

                -_r0 - _r2,
                -_r0 + _r2,
                -_r0,
                -_r0,
            };

            int8x16_t _bd = {
                0,
                0,
                -_i2,
                _i2,

                0,
                0,
                -_i2,
                _i2,

                -_i1,
                -_i1,
                -_i1 - _i2,
                -_i1 + _i2,

                _i1,
                _i1,
                _i1 - _i2,
                _i1 + _i2,
            };

            int8x16_t _ad = {
                0,
                0,
                -_r2,
                _r2,

                0,
                0,
                -_r2,
                _r2,

                -_r1,
                -_r1,
                -_r1 - _r2,
                -_r1 + _r2,

                _r1,
                _r1,
                _r1 - _r2,
                _r1 + _r2,
            };
            
            int8x16_t _bc = {
                -_i0 - _i1 - _i2,
                -_i0 - _i1 + _i2,
                -_i0 - _i1,
                -_i0 - _i1,

                -_i0 + _i1 - _i2,
                -_i0 + _i1 + _i2,
                -_i0 + _i1,
                -_i0 + _i1,

                -_i0 - _i2,
                -_i0 + _i2,
                -_i0,
                -_i0,

                -_i0 - _i2,
                -_i0 + _i2,
                -_i0,
                -_i0,
            };

            lut[_blk].v[_i/3] = (int8x16x4_t){.val= _ac, _bd, _ad, _bc};
        }
        
        // 处理剩余部分
        _r0 = (int8_t) act[_blk].x_real[_i  ];
        _i0 = (int8_t)-act[_blk].x_imag[_i  ];
        _r1 = (int8_t) 0;
        _i1 = (int8_t) 0;
        _r2 = (int8_t) 0;
        _i2 = (int8_t) 0;

        int8x16_t _ac = {
            -_r0 - _r1 - _r2,
            -_r0 - _r1 + _r2,
            -_r0 - _r1,
            -_r0 - _r1,

            -_r0 + _r1 - _r2,
            -_r0 + _r1 + _r2,
            -_r0 + _r1,
            -_r0 + _r1,

            -_r0 - _r2,
            -_r0 + _r2,
            -_r0,
            -_r0,

            -_r0 - _r2,
            -_r0 + _r2,
            -_r0,
            -_r0,
        };

        int8x16_t _bd = {
            0,
            0,
            -_i2,
            _i2,

            0,
            0,
            -_i2,
            _i2,

            -_i1,
            -_i1,
            -_i1 - _i2,
            -_i1 + _i2,

            _i1,
            _i1,
            _i1 - _i2,
            _i1 + _i2,
        };

        int8x16_t _ad = {
            0,
            0,
            -_r2,
            _r2,

            0,
            0,
            -_r2,
            _r2,

            -_r1,
            -_r1,
            -_r1 - _r2,
            -_r1 + _r2,

            _r1,
            _r1,
            _r1 - _r2,
            _r1 + _r2,
        };
        
        int8x16_t _bc = {
            -_i0 - _i1 - _i2,
            -_i0 - _i1 + _i2,
            -_i0 - _i1,
            -_i0 - _i1,

            -_i0 + _i1 - _i2,
            -_i0 + _i1 + _i2,
            -_i0 + _i1,
            -_i0 + _i1,

            -_i0 - _i2,
            -_i0 + _i2,
            -_i0,
            -_i0,

            -_i0 - _i2,
            -_i0 + _i2,
            -_i0,
            -_i0,
        };

        lut[_blk].v[_i/3] = (int8x16x4_t){.val= _ac, _bd, _ad, _bc};
    }
}


void transpose(int m, int k, const block_ifairy *raw_w, block_ifairy_1x3 *w) {
    const int _blk_n = (k+QK_K-1)/QK_K;
    const int i_in_blk_n = (QK_K+2)/3;
    for (int _j = 0; _j < m; _j++) {
        for (int _blk = 0; _blk < _blk_n; _blk++) {
            int _ii = 0;
            #pragma unroll
            for (int _i = 0; _i < QK_K/4; _i+=3) {
                const uint32_t *_p = (uint32_t *)&(raw_w[_j*_blk_n + _blk].qs[_i]);
                uint8_t _iweight_1x3_0, _iweight_1x3_1, _iweight_1x3_2, _iweight_1x3_3;
                if (_i+3 >= QK_K/4) {
                    uint32_t _iweight_1x12 = _p[0] << 16;
                    _iweight_1x3_0 = _iweight_1x12 >> 18 & 0b00111111; 
                    _iweight_1x3_1 = _iweight_1x12 >> 12 & 0b00111111;
                    w[(_j/16)*_blk_n + _blk].qs[_ii  ][_j%16] = three_vals2index_uint8[_iweight_1x3_0] | (_iweight_1x3_0 << 2 & 0b11000000);
                    w[(_j/16)*_blk_n + _blk].qs[_ii+1][_j%16] = three_vals2index_uint8[_iweight_1x3_1] | (_iweight_1x3_1 << 2 & 0b11000000);
                } else {
                    uint32_t _iweight_1x12 = __builtin_bswap32(*_p) >> 8;
                    _iweight_1x3_0 = _iweight_1x12 >> 18 & 0b00111111; 
                    _iweight_1x3_1 = _iweight_1x12 >> 12 & 0b00111111;
                    _iweight_1x3_2 = _iweight_1x12 >> 6  & 0b00111111; 
                    _iweight_1x3_3 = _iweight_1x12       & 0b00111111;
                    w[(_j/16)*_blk_n + _blk].qs[_ii  ][_j%16] = three_vals2index_uint8[_iweight_1x3_0] | (_iweight_1x3_0 << 2 & 0b11000000);
                    w[(_j/16)*_blk_n + _blk].qs[_ii+1][_j%16] = three_vals2index_uint8[_iweight_1x3_1] | (_iweight_1x3_1 << 2 & 0b11000000);
                    w[(_j/16)*_blk_n + _blk].qs[_ii+2][_j%16] = three_vals2index_uint8[_iweight_1x3_2] | (_iweight_1x3_2 << 2 & 0b11000000);
                    w[(_j/16)*_blk_n + _blk].qs[_ii+3][_j%16] = three_vals2index_uint8[_iweight_1x3_3] | (_iweight_1x3_3 << 2 & 0b11000000);
                }
                _ii += 4;
            }
            w[(_j/16)*_blk_n + _blk].d_real[_j%16] = raw_w[_j*_blk_n + _blk].d_real;
            w[(_j/16)*_blk_n + _blk].d_imag[_j%16] = raw_w[_j*_blk_n + _blk].d_imag;
        }
    }
}


void transpose_tmp(int m, int k, const block_ifairy *raw_w, uint8_t *w, float *w_scale_real, float *w_scale_imag) {
    const int _blk_n = (k+QK_K-1)/QK_K;
    const int _in_blk_n = (QK_K+2)/3;
    for (int _j = 0; _j < m; _j++) {
        for (int _blk = 0; _blk < _blk_n; _blk++) {
            int _ii = 0;
            #pragma unroll
            for (int _i = 0; _i < QK_K/4; _i+=3) {
                const uint32_t *_p = (uint32_t *)&(raw_w[_j*_blk_n + _blk].qs[_i]);
                uint8_t _iweight_1x3_0, _iweight_1x3_1, _iweight_1x3_2, _iweight_1x3_3;
                if (_i+3 >= QK_K/4) {
                    uint32_t _iweight_1x12 = _p[0] << 16;
                    _iweight_1x3_0 = _iweight_1x12 >> 18 & 0b00111111; 
                    _iweight_1x3_1 = _iweight_1x12 >> 12 & 0b00111111;
                    w[(((_j/16)*_blk_n + _blk)*_in_blk_n+_ii  )*16+_j%16] = three_vals2index_uint8[_iweight_1x3_0] | (_iweight_1x3_0 << 2 & 0b11000000);
                    w[(((_j/16)*_blk_n + _blk)*_in_blk_n+_ii+1)*16+_j%16] = three_vals2index_uint8[_iweight_1x3_1] | (_iweight_1x3_1 << 2 & 0b11000000);
                } else {
                    uint32_t _iweight_1x12 = __builtin_bswap32(*_p) >> 8;
                    _iweight_1x3_0 = _iweight_1x12 >> 18 & 0b00111111; 
                    _iweight_1x3_1 = _iweight_1x12 >> 12 & 0b00111111;
                    _iweight_1x3_2 = _iweight_1x12 >> 6  & 0b00111111; 
                    _iweight_1x3_3 = _iweight_1x12       & 0b00111111;
                    w[(((_j/16)*_blk_n + _blk)*_in_blk_n+_ii  )*16+(_j%16)] = three_vals2index_uint8[_iweight_1x3_0] | (_iweight_1x3_0 << 2 & 0b11000000);
                    w[(((_j/16)*_blk_n + _blk)*_in_blk_n+_ii+1)*16+(_j%16)] = three_vals2index_uint8[_iweight_1x3_1] | (_iweight_1x3_1 << 2 & 0b11000000);
                    w[(((_j/16)*_blk_n + _blk)*_in_blk_n+_ii+2)*16+(_j%16)] = three_vals2index_uint8[_iweight_1x3_2] | (_iweight_1x3_2 << 2 & 0b11000000);
                    w[(((_j/16)*_blk_n + _blk)*_in_blk_n+_ii+3)*16+(_j%16)] = three_vals2index_uint8[_iweight_1x3_3] | (_iweight_1x3_3 << 2 & 0b11000000);
                }
                _ii += 4;
            }
            w_scale_real[((_j/16)*_blk_n + _blk)*16 +(_j%16)] = raw_w[_j*_blk_n + _blk].d_real;
            w_scale_imag[((_j/16)*_blk_n + _blk)*16 +(_j%16)] = raw_w[_j*_blk_n + _blk].d_imag;
        }
    }
}


static inline int8x16x4_t mul_mat_block_16x3_3x1_with_lut(uint8x16_t iweight_16x3, int8x16x4_t ilut) __attribute__((always_inline));
static inline int8x16x4_t mul_mat_block_16x3_3x1_with_lut(uint8x16_t iweight_16x3, int8x16x4_t ilut) {

    const uint8x16_t _mask_idx = vdupq_n_u8(0b00111111);
    const uint8x16_t _mask_b6  = vdupq_n_u8(0b01000000);
    const uint8x16_t _mask_b7  = vdupq_n_u8(0b10000000);

    // index 适配
    uint8x16_t _index = vandq_u8(iweight_16x3, _mask_idx);

    // 设置标识
    uint8x16_t _fl0 = vtstq_u8(iweight_16x3, _mask_b6);
    uint8x16_t _fl1 = vtstq_u8(iweight_16x3, _mask_b7);

    //查询lut
    int8x16_t _ac_00 = vqtbl1q_s8(ilut.val[0], _index);  // ad_10
    int8x16_t _bd_01 = vqtbl1q_s8(ilut.val[1], _index);  // bc_11
    int8x16_t _ac_11 = vqtbl1q_s8(ilut.val[2], _index);  // ad_00
    int8x16_t _bd_11 = vqtbl1q_s8(ilut.val[3], _index);  // bc_00

    int8x16_t _ac_01 = vnegq_s8(_ac_00);                 // ad_11
    int8x16_t _ac_10 = vnegq_s8(_ac_11);                 // ad_01
    int8x16_t _bd_00 = vnegq_s8(_bd_01);                 // bc_10
    int8x16_t _bd_10 = vnegq_s8(_bd_11);                 // bc_01


    int8x16_t _ac_l = vbslq_u8(_fl0, _ac_01, _ac_00);
    int8x16_t _ad_l = vbslq_u8(_fl0, _ac_10, _ac_11);
    int8x16_t _bc_l = vbslq_u8(_fl0, _bd_10, _bd_11);
    int8x16_t _ac_h = vbslq_u8(_fl0, _ac_11, _ac_10);
    int8x16_t _bc_h = vbslq_u8(_fl0, _bd_01, _bd_00);
    int8x16_t _bd_h = vbslq_u8(_fl0, _bd_11, _bd_10);

    return (int8x16x4_t){
        .val =  vbslq_u8(_fl1, _ac_h, _ac_l),
                vbslq_u8(_fl1, _ac_l, _ad_l),
                vbslq_u8(_fl1, _bc_h, _bc_l),
                vbslq_u8(_fl1, _bd_h, _bc_h)
    };
}


void mul_mat_mxk_kx1_with_lut_q8(int k, int row_begin, int row_end, const block_ifairy_1x3 *w, const lut_block *lut, float *dst) {
    const int _blk_n = (k+QK_K-1)/QK_K;
    const int _in_blk_n = (QK_K+2)/3;
    for (int _row = row_begin; _row <= row_end; _row+=16) {
        const block_ifairy_1x3 *_w_base = w + (_row >> 4)*_blk_n;
        for (int _blk = 0; _blk < _blk_n; _blk++) {
            
            const block_ifairy_1x3 *_w_blk_base = _w_base + _blk;
            const lut_block *_lut_blk_base = lut + _blk;

            int16x8x2_t _block_dst_00 = {.val=vdupq_n_s16(0), vdupq_n_s16(0)};
            int16x8x2_t _block_dst_01 = {.val=vdupq_n_s16(0), vdupq_n_s16(0)};
            int16x8x2_t _block_dst_10 = {.val=vdupq_n_s16(0), vdupq_n_s16(0)};
            int16x8x2_t _block_dst_11 = {.val=vdupq_n_s16(0), vdupq_n_s16(0)};
            
            #pragma unroll
            for (int _i = 0; _i < _in_blk_n; _i++) {
                uint8x16_t _iweight_16x3 = vld1q_u8(_w_blk_base -> qs[_i]);
                int8x16x4_t _iret = mul_mat_block_16x3_3x1_with_lut(_iweight_16x3, _lut_blk_base -> v[_i]);

                _block_dst_00.val[0] = vaddw_s8(_block_dst_00.val[0],  vget_low_s8(_iret.val[0]));
                _block_dst_00.val[1] = vaddw_s8(_block_dst_00.val[1], vget_high_s8(_iret.val[0]));
                _block_dst_01.val[0] = vaddw_s8(_block_dst_01.val[0],  vget_low_s8(_iret.val[1]));
                _block_dst_01.val[1] = vaddw_s8(_block_dst_01.val[1], vget_high_s8(_iret.val[1]));
                _block_dst_10.val[0] = vaddw_s8(_block_dst_10.val[0],  vget_low_s8(_iret.val[2]));
                _block_dst_10.val[1] = vaddw_s8(_block_dst_10.val[1], vget_high_s8(_iret.val[2]));
                _block_dst_11.val[0] = vaddw_s8(_block_dst_11.val[0],  vget_low_s8(_iret.val[3]));
                _block_dst_11.val[1] = vaddw_s8(_block_dst_11.val[1], vget_high_s8(_iret.val[3]));
            }

            // 反量化，写回
            const float _lr = lut[_blk].d_real;
            const float _li = lut[_blk].d_imag;
            const int _w_row = _row >> 4;
            const float *_wr = w[_w_row * _blk_n + _blk].d_real;
            const float *_wi = w[_w_row * _blk_n + _blk].d_imag;
            float *_idst = dst + _row * 2;

#if defined(__M__ALIGNED_16)
                // 加载缩放因子到向量寄存器
                float32x4_t _vlr = vdupq_n_f32(_lr);
                float32x4_t _vli = vdupq_n_f32(_li);
                
                // 处理前8个元素
                float32x4_t _vwr_lo = vld1q_f32(_wr);
                float32x4_t _vwr_hi = vld1q_f32(_wr + 4);
                float32x4_t _vwi_lo = vld1q_f32(_wi);
                float32x4_t _vwi_hi = vld1q_f32(_wi + 4);
                
                float32x4_t _vdst00_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_00.val[0])));
                float32x4_t _vdst00_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_00.val[0])));
                float32x4_t _vdst01_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_01.val[0])));
                float32x4_t _vdst01_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_01.val[0])));
                float32x4_t _vdst10_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_10.val[0])));
                float32x4_t _vdst10_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_10.val[0])));
                float32x4_t _vdst11_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_11.val[0])));
                float32x4_t _vdst11_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_11.val[0])));
                
                float32x4_t _vreal_lo = vmlaq_f32(vmulq_f32(_vdst00_lo, vmulq_f32(_vlr, _vwr_lo)), _vdst11_lo, vmulq_f32(_vli, _vwi_lo));
                float32x4_t _vreal_hi = vmlaq_f32(vmulq_f32(_vdst00_hi, vmulq_f32(_vlr, _vwr_hi)), _vdst11_hi, vmulq_f32(_vli, _vwi_hi));
                float32x4_t _vimag_lo = vmlaq_f32(vmulq_f32(_vdst01_lo, vmulq_f32(_vlr, _vwi_lo)), _vdst10_lo, vmulq_f32(_vli, _vwr_lo));
                float32x4_t _vimag_hi = vmlaq_f32(vmulq_f32(_vdst01_hi, vmulq_f32(_vlr, _vwi_hi)), _vdst10_hi, vmulq_f32(_vli, _vwr_hi));
                
                float32x4x2_t _vout_lo = {_vreal_lo, _vimag_lo};
                float32x4x2_t _vout_hi = {_vreal_hi, _vimag_hi};
                
                float32x4x2_t _vold_lo = vld2q_f32(_idst);
                float32x4x2_t _vold_hi = vld2q_f32(_idst + 8);
                _vout_lo.val[0] = vaddq_f32(_vout_lo.val[0], _vold_lo.val[0]);
                _vout_lo.val[1] = vaddq_f32(_vout_lo.val[1], _vold_lo.val[1]);
                _vout_hi.val[0] = vaddq_f32(_vout_hi.val[0], _vold_hi.val[0]);
                _vout_hi.val[1] = vaddq_f32(_vout_hi.val[1], _vold_hi.val[1]);
                
                vst2q_f32(_idst, _vout_lo);
                vst2q_f32(_idst + 8, _vout_hi);
                
                
                // 处理后8个元素
                _vwr_lo = vld1q_f32(_wr + 8);
                _vwr_hi = vld1q_f32(_wr + 12);
                _vwi_lo = vld1q_f32(_wi + 8);
                _vwi_hi = vld1q_f32(_wi + 12);
                
                _vdst00_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_00.val[1])));
                _vdst00_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_00.val[1])));
                _vdst01_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_01.val[1])));
                _vdst01_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_01.val[1])));
                _vdst10_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_10.val[1])));
                _vdst10_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_10.val[1])));
                _vdst11_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_11.val[1])));
                _vdst11_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_11.val[1])));
                
                _vreal_lo = vmlaq_f32(vmulq_f32(_vdst00_lo, vmulq_f32(_vlr, _vwr_lo)), _vdst11_lo, vmulq_f32(_vli, _vwi_lo));
                _vreal_hi = vmlaq_f32(vmulq_f32(_vdst00_hi, vmulq_f32(_vlr, _vwr_hi)), _vdst11_hi, vmulq_f32(_vli, _vwi_hi));
                _vimag_lo = vmlaq_f32(vmulq_f32(_vdst01_lo, vmulq_f32(_vlr, _vwi_lo)), _vdst10_lo, vmulq_f32(_vli, _vwr_lo));
                _vimag_hi = vmlaq_f32(vmulq_f32(_vdst01_hi, vmulq_f32(_vlr, _vwi_hi)), _vdst10_hi, vmulq_f32(_vli, _vwr_hi));
                
                _vout_lo = (float32x4x2_t){_vreal_lo, _vimag_lo};
                _vout_hi = (float32x4x2_t){_vreal_hi, _vimag_hi};
                
                _vold_lo = vld2q_f32(_idst + 16);
                _vold_hi = vld2q_f32(_idst + 24);
                _vout_lo.val[0] = vaddq_f32(_vout_lo.val[0], _vold_lo.val[0]);
                _vout_lo.val[1] = vaddq_f32(_vout_lo.val[1], _vold_lo.val[1]);
                _vout_hi.val[0] = vaddq_f32(_vout_hi.val[0], _vold_hi.val[0]);
                _vout_hi.val[1] = vaddq_f32(_vout_hi.val[1], _vold_hi.val[1]);
                
                vst2q_f32(_idst + 16, _vout_lo);
                vst2q_f32(_idst + 24, _vout_hi);
                
#elif defined(__M__ALIGNED_8)
                // 加载缩放因子到向量寄存器
                float32x4_t _vlr = vdupq_n_f32(_lr);
                float32x4_t _vli = vdupq_n_f32(_li);
                
                // 处理前8个元素
                float32x4_t _vwr_lo = vld1q_f32(_wr);
                float32x4_t _vwr_hi = vld1q_f32(_wr + 4);
                float32x4_t _vwi_lo = vld1q_f32(_wi);
                float32x4_t _vwi_hi = vld1q_f32(_wi + 4);
                
                float32x4_t _vdst00_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_00.val[0])));
                float32x4_t _vdst00_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_00.val[0])));
                float32x4_t _vdst01_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_01.val[0])));
                float32x4_t _vdst01_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_01.val[0])));
                float32x4_t _vdst10_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_10.val[0])));
                float32x4_t _vdst10_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_10.val[0])));
                float32x4_t _vdst11_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_11.val[0])));
                float32x4_t _vdst11_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_11.val[0])));
                
                float32x4_t _vreal_lo = vmlaq_f32(vmulq_f32(_vdst00_lo, vmulq_f32(_vlr, _vwr_lo)), _vdst11_lo, vmulq_f32(_vli, _vwi_lo));
                float32x4_t _vreal_hi = vmlaq_f32(vmulq_f32(_vdst00_hi, vmulq_f32(_vlr, _vwr_hi)), _vdst11_hi, vmulq_f32(_vli, _vwi_hi));
                float32x4_t _vimag_lo = vmlaq_f32(vmulq_f32(_vdst01_lo, vmulq_f32(_vlr, _vwi_lo)), _vdst10_lo, vmulq_f32(_vli, _vwr_lo));
                float32x4_t _vimag_hi = vmlaq_f32(vmulq_f32(_vdst01_hi, vmulq_f32(_vlr, _vwi_hi)), _vdst10_hi, vmulq_f32(_vli, _vwr_hi));
                
                float32x4x2_t _vout_lo = {_vreal_lo, _vimag_lo};
                float32x4x2_t _vout_hi = {_vreal_hi, _vimag_hi};
                
                float32x4x2_t _vold_lo = vld2q_f32(_idst);
                float32x4x2_t _vold_hi = vld2q_f32(_idst + 8);
                _vout_lo.val[0] = vaddq_f32(_vout_lo.val[0], _vold_lo.val[0]);
                _vout_lo.val[1] = vaddq_f32(_vout_lo.val[1], _vold_lo.val[1]);
                _vout_hi.val[0] = vaddq_f32(_vout_hi.val[0], _vold_hi.val[0]);
                _vout_hi.val[1] = vaddq_f32(_vout_hi.val[1], _vold_hi.val[1]);
                
                vst2q_f32(_idst, _vout_lo);
                vst2q_f32(_idst + 8, _vout_hi);

                if (row_end - _row >= 8) {
                    for (int _j = 8; _j < 16; _j++) {
                        int _jj = _j - 8;
                        float __wr = _wr[_j], __wi = _wi[_j];
                        _idst[_j*2  ] += (float)_block_dst_00.val[1][_jj] * (_lr * __wr) + (float)_block_dst_11.val[1][_jj] * (_li * __wi);
                        _idst[_j*2+1] += (float)_block_dst_01.val[1][_jj] * (_lr * __wi) + (float)_block_dst_10.val[1][_jj] * (_li * __wr);
                    }
                }
#else
                for (int _j = 0; _j <= (row_end - _row < 7 ? row_end - _row : 7); _j++) {
                    float __wr = _wr[_j], __wi = _wi[_j];
                    _idst[_j*2  ] += (float)_block_dst_00.val[0][_j] * (_lr * __wr) + (float)_block_dst_11.val[0][_j] * (_li * __wi);
                    _idst[_j*2+1] += (float)_block_dst_01.val[0][_j] * (_lr * __wi) + (float)_block_dst_10.val[0][_j] * (_li * __wr);
                }
                if (row_end - _row >= 8) {
                    for (int _j = 8; _j <= (row_end - _row < 15 ? row_end - _row : 15); _j++) {
                        int _jj = _j - 8;
                        float __wr = _wr[_j], __wi = _wi[_j];
                        _idst[_j*2  ] += (float)_block_dst_00.val[1][_jj] * (_lr * __wr) + (float)_block_dst_11.val[1][_jj] * (_li * __wi);
                        _idst[_j*2+1] += (float)_block_dst_01.val[1][_jj] * (_lr * __wi) + (float)_block_dst_10.val[1][_jj] * (_li * __wr);
                    }
                }
#endif
        }
    }
}


void generate_lut_q8_block_ifairy_q16_tmp(int k, const block_ifairy_q16 *act, int8_t *lut_v, float *lut_scale) {
    const int _blk_n = (k+QK_K-1)/QK_K;
    int ii = 0;
    for (int _blk = 0; _blk < _blk_n; _blk++) {
        lut_scale[_blk*2  ] = act[_blk].d_real;
        lut_scale[_blk*2+1] = act[_blk].d_imag;
        int8_t _r0, _i0, _r1, _i1, _r2, _i2;
        int _i = 0;
        for (; _i+3 < QK_K; _i+=3) {
            
            _r0 = (int8_t) act[_blk].x_real[_i  ];
            _i0 = (int8_t)-act[_blk].x_imag[_i  ];
            _r1 = (int8_t) act[_blk].x_real[_i+1];
            _i1 = (int8_t)-act[_blk].x_imag[_i+1];
            _r2 = (int8_t) act[_blk].x_real[_i+2];
            _i2 = (int8_t)-act[_blk].x_imag[_i+2];
            

            {
                lut_v[ii++] = -_r0 - _r1 - _r2;
                lut_v[ii++] = -_r0 - _r1 + _r2;
                lut_v[ii++] = -_r0 - _r1;
                lut_v[ii++] = -_r0 - _r1;

                lut_v[ii++] = -_r0 + _r1 - _r2;
                lut_v[ii++] = -_r0 + _r1 + _r2;
                lut_v[ii++] = -_r0 + _r1;
                lut_v[ii++] = -_r0 + _r1;

                lut_v[ii++] = -_r0 - _r2;
                lut_v[ii++] = -_r0 + _r2;
                lut_v[ii++] = -_r0;
                lut_v[ii++] = -_r0;

                lut_v[ii++] = -_r0 - _r2;
                lut_v[ii++] = -_r0 + _r2;
                lut_v[ii++] = -_r0;
                lut_v[ii++] = -_r0;
            }
            {
                lut_v[ii++] =  0;
                lut_v[ii++] =  0;
                lut_v[ii++] = -_i2;
                lut_v[ii++] = _i2;

                lut_v[ii++] = 0;
                lut_v[ii++] = 0;
                lut_v[ii++] = -_i2;
                lut_v[ii++] = _i2;

                lut_v[ii++] = -_i1;
                lut_v[ii++] = -_i1;
                lut_v[ii++] = -_i1 - _i2;
                lut_v[ii++] = -_i1 + _i2;

                lut_v[ii++] = _i1;
                lut_v[ii++] = _i1;
                lut_v[ii++] = _i1 - _i2;
                lut_v[ii++] = _i1 + _i2;
            }
            {
                lut_v[ii++] = 0;
                lut_v[ii++] = 0;
                lut_v[ii++] = -_r2;
                lut_v[ii++] = _r2;

                lut_v[ii++] = 0;
                lut_v[ii++] = 0;
                lut_v[ii++] = -_r2;
                lut_v[ii++] = _r2;

                lut_v[ii++] = -_r1;
                lut_v[ii++] = -_r1;
                lut_v[ii++] = -_r1 - _r2;
                lut_v[ii++] = -_r1 + _r2;

                lut_v[ii++] = _r1;
                lut_v[ii++] = _r1;
                lut_v[ii++] = _r1 - _r2;
                lut_v[ii++] = _r1 + _r2;
            }

            {
                lut_v[ii++] = -_i0 - _i1 - _i2;
                lut_v[ii++] = -_i0 - _i1 + _i2;
                lut_v[ii++] = -_i0 - _i1;
                lut_v[ii++] = -_i0 - _i1;

                lut_v[ii++] = -_i0 + _i1 - _i2;
                lut_v[ii++] = -_i0 + _i1 + _i2;
                lut_v[ii++] = -_i0 + _i1;
                lut_v[ii++] = -_i0 + _i1;

                lut_v[ii++] = -_i0 - _i2;
                lut_v[ii++] = -_i0 + _i2;
                lut_v[ii++] = -_i0;
                lut_v[ii++] = -_i0;

                lut_v[ii++] = -_i0 - _i2;
                lut_v[ii++] = -_i0 + _i2;
                lut_v[ii++] = -_i0;
                lut_v[ii++] = -_i0;
            }

        }
        
        // 处理剩余部分
        _r0 = (int8_t) act[_blk].x_real[_i  ];
        _i0 = (int8_t)-act[_blk].x_imag[_i  ];
        _r1 = (int8_t) 0;
        _i1 = (int8_t) 0;
        _r2 = (int8_t) 0;
        _i2 = (int8_t) 0;

        {
            lut_v[ii++] = -_r0 - _r1 - _r2;
            lut_v[ii++] = -_r0 - _r1 + _r2;
            lut_v[ii++] = -_r0 - _r1;
            lut_v[ii++] = -_r0 - _r1;

            lut_v[ii++] = -_r0 + _r1 - _r2;
            lut_v[ii++] = -_r0 + _r1 + _r2;
            lut_v[ii++] = -_r0 + _r1;
            lut_v[ii++] = -_r0 + _r1;

            lut_v[ii++] = -_r0 - _r2;
            lut_v[ii++] = -_r0 + _r2;
            lut_v[ii++] = -_r0;
            lut_v[ii++] = -_r0;

            lut_v[ii++] = -_r0 - _r2;
            lut_v[ii++] = -_r0 + _r2;
            lut_v[ii++] = -_r0;
            lut_v[ii++] = -_r0;
        }

        {
            lut_v[ii++] = 0;
            lut_v[ii++] = 0;
            lut_v[ii++] = -_i2;
            lut_v[ii++] = _i2;

            lut_v[ii++] = 0;
            lut_v[ii++] = 0;
            lut_v[ii++] = -_i2;
            lut_v[ii++] = _i2;

            lut_v[ii++] = -_i1;
            lut_v[ii++] = -_i1;
            lut_v[ii++] = -_i1 - _i2;
            lut_v[ii++] = -_i1 + _i2;

            lut_v[ii++] = _i1;
            lut_v[ii++] = _i1;
            lut_v[ii++] = _i1 - _i2;
            lut_v[ii++] = _i1 + _i2;
        }

        {
            lut_v[ii++] = 0;
            lut_v[ii++] = 0;
            lut_v[ii++] = -_r2;
            lut_v[ii++] = _r2;

            lut_v[ii++] = 0;
            lut_v[ii++] = 0;
            lut_v[ii++] = -_r2;
            lut_v[ii++] = _r2;

            lut_v[ii++] = -_r1;
            lut_v[ii++] = -_r1;
            lut_v[ii++] = -_r1 - _r2;
            lut_v[ii++] = -_r1 + _r2;

            lut_v[ii++] = _r1;
            lut_v[ii++] = _r1;
            lut_v[ii++] = _r1 - _r2;
            lut_v[ii++] = _r1 + _r2;
        }
        
        {
            lut_v[ii++] = -_i0 - _i1 - _i2;
            lut_v[ii++] = -_i0 - _i1 + _i2;
            lut_v[ii++] = -_i0 - _i1;
            lut_v[ii++] = -_i0 - _i1;

            lut_v[ii++] = -_i0 + _i1 - _i2;
            lut_v[ii++] = -_i0 + _i1 + _i2;
            lut_v[ii++] = -_i0 + _i1;
            lut_v[ii++] = -_i0 + _i1;

            lut_v[ii++] = -_i0 - _i2;
            lut_v[ii++] = -_i0 + _i2;
            lut_v[ii++] = -_i0;
            lut_v[ii++] = -_i0;

            lut_v[ii++] = -_i0 - _i2;
            lut_v[ii++] = -_i0 + _i2;
            lut_v[ii++] = -_i0;
            lut_v[ii++] = -_i0;
        }
    }
}


void mul_mat_mxk_kx1_with_lut_q8_tmp(int k, int row_begin, int row_end, const uint8_t *w, const float *w_scale_real, const float *w_scale_imag, const int8_t *lut_v, const float *lut_scale, float *dst) {
    const int _blk_n = (k+QK_K-1)/QK_K;
    const int _in_blk_n = (QK_K+2)/3;
    for (int _row = row_begin; _row <= row_end; _row+=16) {
        const uint8_t *_w_base = w + (_row >> 4)*_blk_n*_in_blk_n*16;
        for (int _blk = 0; _blk < _blk_n; _blk++) {
            
            const uint8_t *_w_blk_base = _w_base + _blk * _in_blk_n * 16;
            const int8_t *_lut_blk_base = lut_v + _blk * _in_blk_n * 64;

            int16x8x2_t _block_dst_00 = {.val=vdupq_n_s16(0), vdupq_n_s16(0)};
            int16x8x2_t _block_dst_01 = {.val=vdupq_n_s16(0), vdupq_n_s16(0)};
            int16x8x2_t _block_dst_10 = {.val=vdupq_n_s16(0), vdupq_n_s16(0)};
            int16x8x2_t _block_dst_11 = {.val=vdupq_n_s16(0), vdupq_n_s16(0)};
            
            #pragma unroll
            for (int _i = 0; _i < _in_blk_n; _i++) {
                uint8x16_t _iweight_16x3 = vld1q_u8(_w_blk_base + _i * 16);
                int8x16x4_t ilut = vld1q_s8_x4(_lut_blk_base + _i * 64);

                int8x16x4_t _iret = mul_mat_block_16x3_3x1_with_lut(_iweight_16x3, ilut);

                _block_dst_00.val[0] = vaddw_s8(_block_dst_00.val[0],  vget_low_s8(_iret.val[0]));
                _block_dst_00.val[1] = vaddw_s8(_block_dst_00.val[1], vget_high_s8(_iret.val[0]));
                _block_dst_01.val[0] = vaddw_s8(_block_dst_01.val[0],  vget_low_s8(_iret.val[1]));
                _block_dst_01.val[1] = vaddw_s8(_block_dst_01.val[1], vget_high_s8(_iret.val[1]));
                _block_dst_10.val[0] = vaddw_s8(_block_dst_10.val[0],  vget_low_s8(_iret.val[2]));
                _block_dst_10.val[1] = vaddw_s8(_block_dst_10.val[1], vget_high_s8(_iret.val[2]));
                _block_dst_11.val[0] = vaddw_s8(_block_dst_11.val[0],  vget_low_s8(_iret.val[3]));
                _block_dst_11.val[1] = vaddw_s8(_block_dst_11.val[1], vget_high_s8(_iret.val[3]));
            }

            // 反量化，写回
            const float _lr = lut_scale[_blk*2  ];
            const float _li = lut_scale[_blk*2+1];
            const int _w_row = _row >> 4;
            const float *_wr = w_scale_real + (_w_row * _blk_n + _blk)*16;
            const float *_wi = w_scale_imag + (_w_row * _blk_n + _blk)*16;
            float *_idst = dst + _row * 2;

#if defined(__M__ALIGNED_16)
                // 加载缩放因子到向量寄存器
                float32x4_t _vlr = vdupq_n_f32(_lr);
                float32x4_t _vli = vdupq_n_f32(_li);
                
                // 处理前8个元素
                float32x4_t _vwr_lo = vld1q_f32(_wr);
                float32x4_t _vwr_hi = vld1q_f32(_wr + 4);
                float32x4_t _vwi_lo = vld1q_f32(_wi);
                float32x4_t _vwi_hi = vld1q_f32(_wi + 4);
                
                float32x4_t _vdst00_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_00.val[0])));
                float32x4_t _vdst00_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_00.val[0])));
                float32x4_t _vdst01_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_01.val[0])));
                float32x4_t _vdst01_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_01.val[0])));
                float32x4_t _vdst10_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_10.val[0])));
                float32x4_t _vdst10_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_10.val[0])));
                float32x4_t _vdst11_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_11.val[0])));
                float32x4_t _vdst11_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_11.val[0])));
                
                float32x4_t _vreal_lo = vmlaq_f32(vmulq_f32(_vdst00_lo, vmulq_f32(_vlr, _vwr_lo)), _vdst11_lo, vmulq_f32(_vli, _vwi_lo));
                float32x4_t _vreal_hi = vmlaq_f32(vmulq_f32(_vdst00_hi, vmulq_f32(_vlr, _vwr_hi)), _vdst11_hi, vmulq_f32(_vli, _vwi_hi));
                float32x4_t _vimag_lo = vmlaq_f32(vmulq_f32(_vdst01_lo, vmulq_f32(_vlr, _vwi_lo)), _vdst10_lo, vmulq_f32(_vli, _vwr_lo));
                float32x4_t _vimag_hi = vmlaq_f32(vmulq_f32(_vdst01_hi, vmulq_f32(_vlr, _vwi_hi)), _vdst10_hi, vmulq_f32(_vli, _vwr_hi));
                
                float32x4x2_t _vout_lo = {_vreal_lo, _vimag_lo};
                float32x4x2_t _vout_hi = {_vreal_hi, _vimag_hi};
                
                float32x4x2_t _vold_lo = vld2q_f32(_idst);
                float32x4x2_t _vold_hi = vld2q_f32(_idst + 8);
                _vout_lo.val[0] = vaddq_f32(_vout_lo.val[0], _vold_lo.val[0]);
                _vout_lo.val[1] = vaddq_f32(_vout_lo.val[1], _vold_lo.val[1]);
                _vout_hi.val[0] = vaddq_f32(_vout_hi.val[0], _vold_hi.val[0]);
                _vout_hi.val[1] = vaddq_f32(_vout_hi.val[1], _vold_hi.val[1]);
                
                vst2q_f32(_idst, _vout_lo);
                vst2q_f32(_idst + 8, _vout_hi);
                
                
                // 处理后8个元素
                _vwr_lo = vld1q_f32(_wr + 8);
                _vwr_hi = vld1q_f32(_wr + 12);
                _vwi_lo = vld1q_f32(_wi + 8);
                _vwi_hi = vld1q_f32(_wi + 12);
                
                _vdst00_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_00.val[1])));
                _vdst00_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_00.val[1])));
                _vdst01_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_01.val[1])));
                _vdst01_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_01.val[1])));
                _vdst10_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_10.val[1])));
                _vdst10_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_10.val[1])));
                _vdst11_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_11.val[1])));
                _vdst11_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_11.val[1])));
                
                _vreal_lo = vmlaq_f32(vmulq_f32(_vdst00_lo, vmulq_f32(_vlr, _vwr_lo)), _vdst11_lo, vmulq_f32(_vli, _vwi_lo));
                _vreal_hi = vmlaq_f32(vmulq_f32(_vdst00_hi, vmulq_f32(_vlr, _vwr_hi)), _vdst11_hi, vmulq_f32(_vli, _vwi_hi));
                _vimag_lo = vmlaq_f32(vmulq_f32(_vdst01_lo, vmulq_f32(_vlr, _vwi_lo)), _vdst10_lo, vmulq_f32(_vli, _vwr_lo));
                _vimag_hi = vmlaq_f32(vmulq_f32(_vdst01_hi, vmulq_f32(_vlr, _vwi_hi)), _vdst10_hi, vmulq_f32(_vli, _vwr_hi));
                
                _vout_lo = (float32x4x2_t){_vreal_lo, _vimag_lo};
                _vout_hi = (float32x4x2_t){_vreal_hi, _vimag_hi};
                
                _vold_lo = vld2q_f32(_idst + 16);
                _vold_hi = vld2q_f32(_idst + 24);
                _vout_lo.val[0] = vaddq_f32(_vout_lo.val[0], _vold_lo.val[0]);
                _vout_lo.val[1] = vaddq_f32(_vout_lo.val[1], _vold_lo.val[1]);
                _vout_hi.val[0] = vaddq_f32(_vout_hi.val[0], _vold_hi.val[0]);
                _vout_hi.val[1] = vaddq_f32(_vout_hi.val[1], _vold_hi.val[1]);
                
                vst2q_f32(_idst + 16, _vout_lo);
                vst2q_f32(_idst + 24, _vout_hi);
                
#elif defined(__M__ALIGNED_8)
                // 加载缩放因子到向量寄存器
                float32x4_t _vlr = vdupq_n_f32(_lr);
                float32x4_t _vli = vdupq_n_f32(_li);
                
                // 处理前8个元素
                float32x4_t _vwr_lo = vld1q_f32(_wr);
                float32x4_t _vwr_hi = vld1q_f32(_wr + 4);
                float32x4_t _vwi_lo = vld1q_f32(_wi);
                float32x4_t _vwi_hi = vld1q_f32(_wi + 4);
                
                float32x4_t _vdst00_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_00.val[0])));
                float32x4_t _vdst00_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_00.val[0])));
                float32x4_t _vdst01_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_01.val[0])));
                float32x4_t _vdst01_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_01.val[0])));
                float32x4_t _vdst10_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_10.val[0])));
                float32x4_t _vdst10_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_10.val[0])));
                float32x4_t _vdst11_lo = vcvtq_f32_s32(vmovl_s16(vget_low_s16(_block_dst_11.val[0])));
                float32x4_t _vdst11_hi = vcvtq_f32_s32(vmovl_s16(vget_high_s16(_block_dst_11.val[0])));
                
                float32x4_t _vreal_lo = vmlaq_f32(vmulq_f32(_vdst00_lo, vmulq_f32(_vlr, _vwr_lo)), _vdst11_lo, vmulq_f32(_vli, _vwi_lo));
                float32x4_t _vreal_hi = vmlaq_f32(vmulq_f32(_vdst00_hi, vmulq_f32(_vlr, _vwr_hi)), _vdst11_hi, vmulq_f32(_vli, _vwi_hi));
                float32x4_t _vimag_lo = vmlaq_f32(vmulq_f32(_vdst01_lo, vmulq_f32(_vlr, _vwi_lo)), _vdst10_lo, vmulq_f32(_vli, _vwr_lo));
                float32x4_t _vimag_hi = vmlaq_f32(vmulq_f32(_vdst01_hi, vmulq_f32(_vlr, _vwi_hi)), _vdst10_hi, vmulq_f32(_vli, _vwr_hi));
                
                float32x4x2_t _vout_lo = {_vreal_lo, _vimag_lo};
                float32x4x2_t _vout_hi = {_vreal_hi, _vimag_hi};
                
                float32x4x2_t _vold_lo = vld2q_f32(_idst);
                float32x4x2_t _vold_hi = vld2q_f32(_idst + 8);
                _vout_lo.val[0] = vaddq_f32(_vout_lo.val[0], _vold_lo.val[0]);
                _vout_lo.val[1] = vaddq_f32(_vout_lo.val[1], _vold_lo.val[1]);
                _vout_hi.val[0] = vaddq_f32(_vout_hi.val[0], _vold_hi.val[0]);
                _vout_hi.val[1] = vaddq_f32(_vout_hi.val[1], _vold_hi.val[1]);
                
                vst2q_f32(_idst, _vout_lo);
                vst2q_f32(_idst + 8, _vout_hi);

                if (row_end - _row >= 8) {
                    for (int _j = 8; _j < 16; _j++) {
                        int _jj = _j - 8;
                        float __wr = _wr[_j], __wi = _wi[_j];
                        _idst[_j*2  ] += (float)_block_dst_00.val[1][_jj] * (_lr * __wr) + (float)_block_dst_11.val[1][_jj] * (_li * __wi);
                        _idst[_j*2+1] += (float)_block_dst_01.val[1][_jj] * (_lr * __wi) + (float)_block_dst_10.val[1][_jj] * (_li * __wr);
                    }
                }
#else
                for (int _j = 0; _j <= (row_end - _row < 7 ? row_end - _row : 7); _j++) {
                    float __wr = _wr[_j], __wi = _wi[_j];
                    _idst[_j*2  ] += (float)_block_dst_00.val[0][_j] * (_lr * __wr) + (float)_block_dst_11.val[0][_j] * (_li * __wi);
                    _idst[_j*2+1] += (float)_block_dst_01.val[0][_j] * (_lr * __wi) + (float)_block_dst_10.val[0][_j] * (_li * __wr);
                }
                if (row_end - _row >= 8) {
                    for (int _j = 8; _j <= (row_end - _row < 15 ? row_end - _row : 15); _j++) {
                        int _jj = _j - 8;
                        float __wr = _wr[_j], __wi = _wi[_j];
                        _idst[_j*2  ] += (float)_block_dst_00.val[1][_jj] * (_lr * __wr) + (float)_block_dst_11.val[1][_jj] * (_li * __wi);
                        _idst[_j*2+1] += (float)_block_dst_01.val[1][_jj] * (_lr * __wi) + (float)_block_dst_10.val[1][_jj] * (_li * __wr);
                    }
                }
#endif
        }
    }
}


void mul_mat_mxk_kxn_with_lut_q8_base(int k, int row_begin, int row_end, int col_begin, int col_end, const block_ifairy_1x3 *w, const lut_block *lut_base, float *dst_base) {
    const int _blk_n = (k+QK_K-1)/QK_K;
    for (int _col = col_begin; _col <= col_end; _col++) {
        float *_dst = dst_base + _col*k;
        const lut_block *_lut = lut_base + _col*_blk_n;
        mul_mat_mxk_kx1_with_lut_q8(k, row_begin, row_end, w, lut_base, dst_base);
    }
}


void mul_mat_mxk_kxn_with_lut_q8(int k, int n, int row_begin, int row_end, const block_ifairy_1x3 *w, const lut_block *lut, float *dst) {
    mul_mat_mxk_kxn_with_lut_q8_base(k, row_begin, row_end, 0, n-1, w, lut, dst);
}


block_ifairy_q16 *alloc_act_block_ifairy_q16(int k) {
    int _blk_n = (k+QK_K-1)/QK_K;
    return calloc(_blk_n, sizeof(block_ifairy_q16)); 
}


lut_block *alloc_lut_q8(int k) {
    int _blk_n = (k+QK_K-1)/QK_K;
    return calloc(_blk_n, sizeof(lut_block));
}


block_ifairy_1x3 *alloc_w(int m, int k) {
    int _blk_n = (k+QK_K-1)/QK_K;
    int _row_n = (m+15)/16;
    block_ifairy_1x3 *_ptr = calloc(_blk_n*_row_n, sizeof(block_ifairy_1x3));
    return _ptr;
}


uint8_t *alloc_w_tmp(int m, int k) {
    int _blk_n = (k+QK_K-1)/QK_K;
    int _in_blk_n = (QK_K+2)/3;
    int _row_n = (m+15)/16;
    return calloc(_blk_n * _in_blk_n * 16 * _row_n, sizeof(uint8_t));
}


float *alloc_w_scale_real_tmp(int m, int k) {
    int _blk_n = (k+QK_K-1)/QK_K;
    return calloc(_blk_n * m, sizeof(float));
}


float *alloc_w_scale_imag_tmp(int m, int k) {
    int _blk_n = (k+QK_K-1)/QK_K;
    return calloc(_blk_n * m, sizeof(float));
}


int8_t *alloc_lut_v_q8(int k) {
    int _blk_n = (k+QK_K-1)/QK_K;
    int _in_blk_n = (QK_K+2)/3;
    return calloc(_blk_n * _in_blk_n * 64, sizeof(int8_t));
}


float *alloc_lut_scale_q8(int k) {
    int _blk_n = (k+QK_K-1)/QK_K;
    return calloc(_blk_n * 2, sizeof(float));
}


void free_act_block_ifairy_q16(block_ifairy_q16 *act) {
    free(act);
}


void free_lut_q8(lut_block *lut) {
    free(lut);
}


void free_w(block_ifairy_1x3 *w) {
    free(w);
}


void free_w_tmp(uint8_t *w) {
    free(w);
}


void free_w_scale_real_tmp(float *w_scale_real) {
    free(w_scale_real);
}


void free_w_scale_imag_tmp(float *w_scale_imag) {
    free(w_scale_imag);
}


void free_lut_v_q8_tmp(int8_t *lut_v) {
    free(lut_v);
}


void free_lut_scale_q8_tmp(float *lut_scale) {
    free(lut_scale);
}



// ================================ 乘法程序 ================================
static inline float32x2_t mul_mat_block_1x4_4x1(uint8_t a, float32_t *b, float32_t d_real, float32_t d_imag) __attribute__((always_inline));
static inline float32x2_t mul_mat_block_1x4_4x1(uint8_t a, float32_t *b, float32_t d_real, float32_t d_imag) {
    float32_t _acc_r = 0;
    float32_t _acc_i = 0;

    for (int _k = 0; _k < 4; _k++) {
        uint8_t _code = (a >> (2 * (3-_k))) & 0x3;

        float32_t _xr = b[2 * _k];
        float32_t _xi = -b[2 * _k + 1]; // 共轭

        switch (_code) {
        case 0b00:  // -1
            _acc_r += -_xr * d_real;
            _acc_i += -_xi * d_real;
            break;
        case 0b01:  // +1
            _acc_r += _xr * d_real;
            _acc_i += _xi * d_real;
            break;
        case 0b10:  // -i
            _acc_r += _xi * d_imag;
            _acc_i += -_xr * d_imag;
            break;
        case 0b11:  // +i
            _acc_r += -_xi * d_imag;
            _acc_i += _xr * d_imag;
            break;
        }
    }
    return vset_lane_f32(_acc_i, vset_lane_f32(_acc_r, vdup_n_f32(0), 0), 1);
}


void mul_mat_mxk_kx1(int k, int row_begin, int row_end, const block_ifairy *w, const float *act, float *dst) {
    int _block_n = (k+QK_K-1)/QK_K;
    for (int _row = row_begin; _row <= row_end; _row++) {
        for (int _block = 0; _block < _block_n; _block++) {
            for (int _i = 0; _i < QK_K/4; _i++) {
                uint8_t _w4 = w[_row*_block_n+_block].qs[_i];
                float32_t _a4[8] = {
                    act[_block*QK_K*2+_i*8],
                    act[_block*QK_K*2+_i*8+1],
                    act[_block*QK_K*2+_i*8+2],
                    act[_block*QK_K*2+_i*8+3],
                    act[_block*QK_K*2+_i*8+4],
                    act[_block*QK_K*2+_i*8+5],
                    act[_block*QK_K*2+_i*8+6],
                    act[_block*QK_K*2+_i*8+7],
                };
                float32x2_t _ret = mul_mat_block_1x4_4x1(_w4, _a4, w[_row*_block_n+_block].d_real, w[_row*_block_n+_block].d_imag);
                dst[2*_row]   += (float)_ret[0];
                dst[2*_row+1] += (float)_ret[1];
            }
        }
    }
}



// ================================ 1路LUT ================================
static const int8_t sx[4] = { -1,  1,  0,  0 };
static const int8_t sy[4] = {  0,  0, -1,  1 };


void generate_lut_q16_block_ifairy_q16_old(int k, const block_ifairy_q16 *act, int16_t *lut_v, float *lut_scale) {
    int _block_n = (k+QK_K-1)/QK_K;
    int _ii = 0;
    for (int _blk = 0; _blk < _block_n; _blk++) {
        lut_scale[_blk*2  ] = act[_blk].d_real;
        lut_scale[_blk*2+1] = act[_blk].d_imag;
        // 生成lut表
        for (int _i = 0; _i < QK_K; _i+=3) {
            int8_t _r0, _i0, _r1, _i1, _r2, _i2;
            if (_i+3 >= QK_K) {
                _r0 = (int8_t) act[_blk].x_real[_i  ];
                _i0 = (int8_t) act[_blk].x_imag[_i  ];
                _r1 = (int8_t) 0;
                _i1 = (int8_t) 0;
                _r2 = (int8_t) 0;
                _i2 = (int8_t) 0;
            } else {
                _r0 = (int8_t) act[_blk].x_real[_i  ];
                _i0 = (int8_t) act[_blk].x_imag[_i  ];
                _r1 = (int8_t) act[_blk].x_real[_i+1];
                _i1 = (int8_t) act[_blk].x_imag[_i+1];
                _r2 = (int8_t) act[_blk].x_real[_i+2];
                _i2 = (int8_t) act[_blk].x_imag[_i+2];
            } 
            for (int __i = 0; __i < 4; __i++) {
                for (int __j = 0; __j < 4; __j++) {
                    for (int __k = 0; __k < 4; __k++) {
                        int16_t _ac =
                            _r0 * sx[__i] +
                            _r1 * sx[__j] +
                            _r2 * sx[__k];

                        int16_t _bd =
                            _i0 * sy[__i] +
                            _i1 * sy[__j] +
                            _i2 * sy[__k];

                        int16_t _ad =
                            _r0 * sy[__i] +
                            _r1 * sy[__j] +
                            _r2 * sy[__k];

                        int16_t _bc =
                            _i0 * sx[__i] +
                            _i1 * sx[__j] +
                            _i2 * sx[__k];

                        lut_v[_ii++] = _ac;
                        lut_v[_ii++] = _bd;
                        lut_v[_ii++] = _ad;
                        lut_v[_ii++] = _bc;
                    }
                }
            }
        }
    }
}


void transpose_old(int m, int k, const block_ifairy *raw_w, block_ifairy_1x3_old *w) {
    const int _blk_n = (k+QK_K-1)/QK_K;
    for (int _j = 0; _j < m; _j++) {
        for (int _blk = 0; _blk < _blk_n; _blk++) {
            int _ii = 0;
            #pragma unroll
            for (int _i = 0; _i < QK_K/4; _i+=3) {
                const uint32_t *_p = (uint32_t *)&(raw_w[_j*_blk_n + _blk].qs[_i]);
                uint8_t _iweight_1x3_0, _iweight_1x3_1, _iweight_1x3_2, _iweight_1x3_3;
                if (_i+3 >= QK_K/4) {
                    uint32_t _iweight_1x12 = _p[0] << 16;
                    _iweight_1x3_0 = _iweight_1x12 >> 18 & 0b00111111; 
                    _iweight_1x3_1 = _iweight_1x12 >> 12 & 0b00111111;
                    w[_j*_blk_n + _blk].qs[_ii  ] = _iweight_1x3_0;
                    w[_j*_blk_n + _blk].qs[_ii+1] = _iweight_1x3_1;
                } else {
                    uint32_t _iweight_1x12 = __builtin_bswap32(*_p) >> 8;
                    _iweight_1x3_0 = _iweight_1x12 >> 18 & 0b00111111; 
                    _iweight_1x3_1 = _iweight_1x12 >> 12 & 0b00111111;
                    _iweight_1x3_2 = _iweight_1x12 >> 6  & 0b00111111; 
                    _iweight_1x3_3 = _iweight_1x12       & 0b00111111;
                    w[_j*_blk_n + _blk].qs[_ii  ] = _iweight_1x3_0;
                    w[_j*_blk_n + _blk].qs[_ii+1] = _iweight_1x3_1;
                    w[_j*_blk_n + _blk].qs[_ii+2] = _iweight_1x3_2;
                    w[_j*_blk_n + _blk].qs[_ii+3] = _iweight_1x3_3;
                }
                _ii += 4;
            }
            w[_j*_blk_n + _blk].d_real = raw_w[_j*_blk_n + _blk].d_real;
            w[_j*_blk_n + _blk].d_imag = raw_w[_j*_blk_n + _blk].d_imag;
        }
    }
}


void mul_mat_mxk_kx1_with_lut_q16_old(int k, int row_begin, int row_end, const block_ifairy_1x3_old *w, const int16_t *lut, const float *lut_scale, float *dst) {
    const int _blk_n = (k+QK_K-1)/QK_K;
    const int _in_blk_n = (QK_K+2)/3;
    const int _next_lut_blk = 256;
    const int _lut_blk_n = 4; 
    for (int _row = row_begin; _row <= row_end; _row++) {
        const int _w_base = _row*_blk_n;
        for (int _blk = 0; _blk < _blk_n; _blk++) {
            int32x4_t _isum = vdupq_n_s32(0);
            const int16_t *_lut_base = lut + (_blk*_in_blk_n)*_next_lut_blk;
            #pragma unroll
            for (int _i = 0; _i < _in_blk_n; _i++) {
                const uint8_t _iweight_16x3 = w[_w_base+_blk].qs[_i];
                int16x4_t _abcd = vld1_s16(_lut_base+_next_lut_blk*_i+_iweight_16x3*_lut_blk_n);
                _isum = vaddw_s16(_isum, _abcd);
            }
            // 反量化，写回
            dst[(_row)*2  ] += (float)(_isum[0]) * lut_scale[_blk*2] * w[(_row)*_blk_n+_blk].d_real + (float)(_isum[1]) * lut_scale[_blk*2+1] * w[(_row)*_blk_n+_blk].d_imag;
            dst[(_row)*2+1] += (float)(_isum[2]) * lut_scale[_blk*2] * w[(_row)*_blk_n+_blk].d_imag - (float)(_isum[3]) * lut_scale[_blk*2+1] * w[(_row)*_blk_n+_blk].d_real;
        }
    }
}


int16_t *alloc_lut_v_q16_old(int k) {
    int _blk_n = (k+QK_K-1)/QK_K;
    int _in_blk_n = (QK_K+2)/3;
    return calloc(_blk_n*_in_blk_n*256, sizeof(int16_t));
}


float *alloc_lut_scale_old(int k) {
    int _blk_n = (k+QK_K-1)/QK_K;
    return calloc(_blk_n*2, sizeof(float));
}


block_ifairy_1x3_old *alloc_w_old(int m, int k) {
    int _blk_n = (k+QK_K-1)/QK_K;
    int _row_n = (m+15)/16;
    block_ifairy_1x3_old *_ptr = calloc(_blk_n*_row_n*16, sizeof(block_ifairy_1x3_old));
    return _ptr;
}


void free_lut_v_q16_old(int16_t *lut) {
    free(lut);
}


void free_lut_scale_old(float *scale) {
    free(scale);
}


void free_w_old(block_ifairy_1x3_old *w) {
    free(w);
}