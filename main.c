#include <stdio.h>
#include <arm_neon.h>
#include <stdint.h>
#include "lut/mul_mat_with_lut.h"

// void print_int8x16(int8x16_t v) {
//     int8_t tmp[16];
//     vst1q_s8(tmp, v);
//     for (int i = 0; i < 16; i++) {
//         printf("%d ", tmp[i]);
//     }
//     printf("\n");
// }
//
//
// int main() {
//     int16_t act[26];
//     // 示例：生成 0~25
//     for (int i = 0; i < 26; i++) {
//         act[i] = (int16_t)(i);
//     }
//
//     // 打印数组
//     for (int i = 0; i < 26; i++) {
//         if (i % 6 == 0) {
//             printf("[");
//         }
//         printf(" %d ", act[i]);
//         if (i % 6 == 5) {
//             printf("]");
//         }
//     }
//     printf("\n\n");
//
//
//     int8x16_pair lut[4];
//
//     printf("%d\n",generate_lut_int8(act, 26, lut));
//
//     for (int i = 0; i < 5; i++) {
//         printf("%d:\n", i);
//         print_int8x16(lut[i].v[0]);
//         print_int8x16(lut[i].v[1]);
//     }
//
//     return 0;
// }



int main() {
    uint8x16_t flag = {
            3,
            1,
            2,
            3,
            1,
            0,
            2,
            1,

            3,
            3,
            0,
            0,
            1,
            2,
            1,
            3,
    };
    int8x16_t ret_r_tmp0 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    int8x16_t ret_i_tmp0 = vaddq_s8(vdupq_n_u8(100), ret_r_tmp0);

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


    for (int i = 0; i < 16; i++) {
        printf("%d %d\n", ret_r[i], ret_i[i]);
    }

}






/*
        int16x4_t v0[16] =
        {
                // (-1, -1, -1)
                -r0 - r1 - r2,
                -i0 - i1 - i2,

                // (-1, -1, 1)
                -r0 - r1 + r2,
                -i0 - i1 + i2,

                // (-1, -1, -i)
                -r0 - r1 + i2,
                -i0 - i1 - r2,

                // (-1, -1, i)
                -r0 - r1 - i2,
                -i0 - i1 + r2,

                // (-1, 1, -1)
                -r0 + r1 - r2,
                -i0 + i1 - i2,

                // (-1, 1, 1)
                -r0 + r1 + r2,
                -i0 + i1 + i2,

                // (-1, 1, -i)
                -r0 + r1 + i2,
                -i0 + i1 - r2,

                // (-1, 1, i)
                -r0 + r1 - i2,
                -i0 + i1 + r2
        };


        int16x4_t v1[16] = {
                // (-1, -i, -1)
                -r0 + i1 - r2,
                -i0 - r1 - i2,

                // (-1, -i, 1)
                -r0 + i1 + r2,
                -i0 - r1 + i2,

                // (-1, -i, -i)
                -r0 + i1 - i2,
                -i0 - r1 - r2,

                // (-1, -i, i)
                -r0 + i1 + i2,
                -i0 - r1 + r2,

                // (-1, i, -1)
                -r0 - i1 - r2,
                -i0 + r1 - i2,

                // (-1, i, 1)
                -r0 - i1 + r2,
                -i0 + r1 + i2,

                // (-1, i, -i)
                -r0 - i1 - i2,
                -i0 + r1 - r2,

                // (-1, i, i)
                -r0 - i1 + i2,
                -i0 + r1 + r2
        };
*/