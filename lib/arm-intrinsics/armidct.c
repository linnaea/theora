#include "armfrag.h"

#if defined(OC_ARM_ASM_NEON)
typedef int32_t ogg_int32_t;
#include "../dct.h"
#include "neon_a64_compat.h"
#include "neon_transpose8x8.h"

static int16x8_t scale16(int16x8_t a, int32_t b) {
    int32x4_t al = vmovl_s16(vget_low_s16(a));
    int32x4_t ah = vmovl_high_s16(a);
    al = vmulq_n_s32(al, b);
    ah = vmulq_n_s32(ah, b);
    return vuzp2q_s16(vreinterpretq_s16_s32(al), vreinterpretq_s16_s32(ah));
}

static void idct8x8(const int16x8_t x[8], int16x8_t y[8]) {
    int16x8_t t[8], r;

    /*Stage 1:*/
    /*0-1 butterfly.*/
    t[0] = scale16(vaddq_s16(x[0], x[4]), OC_C4S4);
    t[1] = scale16(vsubq_s16(x[0], x[4]), OC_C4S4);
    /*2-3 rotation by 6pi/16.*/
    t[2] = vsubq_s16(scale16(x[2], OC_C6S2), scale16(x[6], OC_C2S6));
    t[3] = vaddq_s16(scale16(x[2], OC_C2S6), scale16(x[6], OC_C6S2));
    /*4-7 rotation by 7pi/16.*/
    t[4] = vsubq_s16(scale16(x[1], OC_C7S1), scale16(x[7], OC_C1S7));
    t[7] = vaddq_s16(scale16(x[1], OC_C1S7), scale16(x[7], OC_C7S1));
    /*5-6 rotation by 3pi/16.*/
    t[5] = vsubq_s16(scale16(x[5], OC_C3S5), scale16(x[3], OC_C5S3));
    t[6] = vaddq_s16(scale16(x[5], OC_C5S3), scale16(x[3], OC_C3S5));

    /*Stage 2:*/
    /*4-5 butterfly.*/
    r = vaddq_s16(t[4], t[5]);
    t[5] = scale16(vsubq_s16(t[4], t[5]), OC_C4S4);
    t[4] = r;
    /*7-6 butterfly.*/
    r = vaddq_s16(t[7], t[6]);
    t[6] = scale16(vsubq_s16(t[7], t[6]), OC_C4S4);
    t[7] = r;

    /*Stage 3:*/
    /*0-3 butterfly.*/
    r = vaddq_s16(t[0], t[3]);
    t[3] = vsubq_s16(t[0], t[3]);
    t[0] = r;
    /*1-2 butterfly.*/
    r = vaddq_s16(t[1], t[2]);
    t[2] = vsubq_s16(t[1], t[2]);
    t[1] = r;
    /*6-5 butterfly.*/
    r = vaddq_s16(t[6], t[5]);
    t[5] = vsubq_s16(t[6], t[5]);
    t[6] = r;

    /*Stage 4:*/
    /*0-7 butterfly.*/
    y[0] = vaddq_s16(t[0], t[7]);
    y[7] = vsubq_s16(t[0], t[7]);
    /*1-6 butterfly.*/
    y[1] = vaddq_s16(t[1], t[6]);
    y[6] = vsubq_s16(t[1], t[6]);
    /*2-5 butterfly.*/
    y[2] = vaddq_s16(t[2], t[5]);
    y[5] = vsubq_s16(t[2], t[5]);
    /*3-4 butterfly.*/
    y[3] = vaddq_s16(t[3], t[4]);
    y[4] = vsubq_s16(t[3], t[4]);
}

static void idct8x4(const int16x8_t x[4], int16x8_t y[8]) {
    int16x8_t t[8], r;

    /*Stage 1:*/
    t[0] = scale16(x[0], OC_C4S4);
    t[1] = scale16(x[0], OC_C4S4);
    t[2] = scale16(x[2], OC_C6S2);
    t[3] = scale16(x[2], OC_C2S6);
    t[4] = scale16(x[1], OC_C7S1);
    t[7] = scale16(x[1], OC_C1S7);
    t[5] = vnegq_s16(scale16(x[3], OC_C5S3));
    t[6] = scale16(x[3], OC_C3S5);

    /*Stage 2:*/
    r = vaddq_s16(t[4], t[5]);
    t[5] = scale16(vsubq_s16(t[4], t[5]), OC_C4S4);
    t[4] = r;
    r = vaddq_s16(t[7], t[6]);
    t[6] = scale16(vsubq_s16(t[7], t[6]), OC_C4S4);
    t[7] = r;

    /*Stage 3:*/
    r = vaddq_s16(t[0], t[3]);
    t[3] = vsubq_s16(t[0], t[3]);
    t[0] = r;
    r = vaddq_s16(t[1], t[2]);
    t[2] = vsubq_s16(t[1], t[2]);
    t[1] = r;
    r = vaddq_s16(t[6], t[5]);
    t[5] = vsubq_s16(t[6], t[5]);
    t[6] = r;

    /*Stage 4:*/
    y[0] = vaddq_s16(t[0], t[7]);
    y[7] = vsubq_s16(t[0], t[7]);
    y[1] = vaddq_s16(t[1], t[6]);
    y[6] = vsubq_s16(t[1], t[6]);
    y[2] = vaddq_s16(t[2], t[5]);
    y[5] = vsubq_s16(t[2], t[5]);
    y[3] = vaddq_s16(t[3], t[4]);
    y[4] = vsubq_s16(t[3], t[4]);
}

static void scale_final(const int16x8_t y[8], int16_t _y[64]) {
    for (int i = 0; i < 8; i++)
        vst1q_s16(&_y[i * 8], vrshrq_n_s16(y[i], 4));
}

static void oc_idct8x8_4(int16_t _y[64], int16_t _x[64]) {
    int16x8_t x[8], y[8];
    for (int i = 0; i < 4; i++)
        x[i] = vld1q_s16(&_x[i * 8]);
    idct8x4(x, y);
    v_transpose8x8(y, x);
    idct8x4(x, y);
    scale_final(y, _y);
    for (int i = 0; i < 4; i++)
        vst1_s16(&_x[i * 8], vcreate_s16(0));
}

static void oc_idct8x8_8(int16_t _y[64], int16_t _x[64]) {
    int16x8_t x[8], y[8];
    for (int i = 0; i < 8; i++)
        x[i] = vld1q_s16(&_x[i * 8]);
    idct8x8(x, y);
    v_transpose8x8(y, x);
    idct8x8(x, y);
    scale_final(y, _y);
    for (int i = 0; i < 8; i++)
        vst1q_s16(&_x[i * 8], veorq_s16(x[0], x[0]));
}

void oc_idct8x8_neon(int16_t _y[64],int16_t _x[64],int _last_zzi) {
    if (_last_zzi <= 10)
        oc_idct8x8_4(_y, _x);
    else
        oc_idct8x8_8(_y, _x);
}
#endif
