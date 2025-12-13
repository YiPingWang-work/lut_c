#ifndef MUL_MAT_WITH_LUT
#define MUL_MAT_WITH_LUT

#include <arm_neon.h>

#define QK_K 256
typedef uint16_t ggml_half;


typedef struct {
    uint8_t qs[QK_K/4]; // 2 bits per element
    ggml_half d_real, d_imag;
} block_ifairy;

static const uint8x16_t uint8x16_t_0_15 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
static const uint8x16_t uint8x16_t_0_15_swapped = {1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14};


void generate_lut_int8(const int16_t *activation, int m, int8x16x2_t *lut);
void mul_mat_nxm_mx1(void *weight, int n, int m, int8x16x2_t *lut, int8_t *output);

#endif