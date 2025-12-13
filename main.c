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
    uint8x16_t iweight_3x16 = {0, 15, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

    int8x16_t ret_r_tmp0 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    int8x16_t ret_i_tmp0 = vaddq_s8(vdupq_n_u8(100), ret_r_tmp0);
    uint8x16_t index = vandq_u8(iweight_3x16, vdupq_n_u8(7));
    int8x16x2_t ilut = {
        ret_r_tmp0,
        ret_i_tmp0,
    };


    int8x16_t r0 = vqtbl1q_s8(ilut.val[0], index); // 全部的16个数据在ilut0中的实部
    int8x16_t i0 = vqtbl1q_s8(ilut.val[0], vaddq_u8(index, vdupq_n_u8(8))); // 全部的16个数据在ilut0中的虚部
    int8x16_t r1 = vqtbl1q_s8(ilut.val[1], index); // 全部的16个数据在ilut1中的实部
    int8x16_t i1 = vqtbl1q_s8(ilut.val[1], vaddq_u8(index, vdupq_n_u8(8))); // 全部的16个数据在ilut1中的虚部



    for (int i = 0; i < 16; i++) {
        printf("%d %d\n", r0[i], i0[i]);
    }

    printf(" ===\n");

    for (int i = 0; i < 16; i++) {
        printf("%d %d\n", r1[i], i1[i]);
    }
    
    printf(" ===\n");

    uint8x16_t mask = vceqq_u8(vandq_u8(iweight_3x16, vdupq_n_u8(8)), vdupq_n_u8(8));


    int8x16_t iret_r_tmp0 = vbslq_u8(mask, r1, r0);
    int8x16_t iret_i_tmp0 = vbslq_u8(mask, i1, i0);

    for (int i = 0; i < 16; i++) {
        printf("%d %d %d\n", mask[i], iret_r_tmp0[i], iret_i_tmp0[i]);
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