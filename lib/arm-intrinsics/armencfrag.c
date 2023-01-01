#include <arm_neon.h>
#include <stdint.h>
#include <stdlib.h>

#ifndef __aarch64__
#  define vuzp1q_s64(a, b) (vcombine_s64(vget_low_s64(a), vget_low_s64(b)))
#  define vuzp2q_s64(a, b) (vcombine_s64(vget_high_s64(a), vget_high_s64(b)))

static inline int32_t vaddlvq_s16(int16x8_t p) {
    int32x4_t q = vpaddlq_s16(p);
    int32x2_t r = vpadd_s32(vget_high_s32(q), vget_low_s32(q));
    return (vget_lane_s32(r, 0) + vget_lane_s32(r, 1));
}

static inline uint32_t vaddlvq_u16(uint16x8_t p) {
    uint32x4_t q = vpaddlq_u16(p);
    uint32x2_t r = vpadd_u32(vget_high_u32(q), vget_low_u32(q));
    return (vget_lane_u32(r, 0) + vget_lane_u32(r, 1));
}

static inline uint32_t vaddvq_u32(uint32x4_t p) {
    uint32x2_t r = vpadd_u32(vget_high_u32(p), vget_low_u32(p));
    return (vget_lane_u32(r, 0) + vget_lane_u32(r, 1));
}
#endif

static inline unsigned oc_hadamard_satd_neon_8x8(int *_dc,const int16x8_t _rows[8]) {
    // in: 9 bit [-255..255]
    // Hadamard stage 1: 4x2
    int16x8_t p0 = vaddq_s16(_rows[0], _rows[4]);
    int16x8_t p1 = vaddq_s16(_rows[1], _rows[5]);
    int16x8_t p2 = vaddq_s16(_rows[2], _rows[6]);
    int16x8_t p3 = vaddq_s16(_rows[3], _rows[7]);
    int16x8_t p4 = vsubq_s16(_rows[0], _rows[4]);
    int16x8_t p5 = vsubq_s16(_rows[1], _rows[5]);
    int16x8_t p6 = vsubq_s16(_rows[2], _rows[6]);
    int16x8_t p7 = vsubq_s16(_rows[3], _rows[7]);
    // 10 bit

    // Hadamard stage 2: 2x4
    int16x8_t q0 = vaddq_s16(p0, p2);
    int16x8_t q1 = vaddq_s16(p1, p3);
    int16x8_t q2 = vsubq_s16(p0, p2);
    int16x8_t q3 = vsubq_s16(p1, p3);
    int16x8_t q4 = vaddq_s16(p4, p6);
    int16x8_t q5 = vaddq_s16(p5, p7);
    int16x8_t q6 = vsubq_s16(p4, p6);
    int16x8_t q7 = vsubq_s16(p5, p7);
    // 11 bit

    // Hadamard stage 3: 1x8
    p0 = vaddq_s16(q0, q1);
    p1 = vsubq_s16(q0, q1);
    p2 = vaddq_s16(q2, q3);
    p3 = vsubq_s16(q2, q3);
    p4 = vaddq_s16(q4, q5);
    p5 = vsubq_s16(q4, q5);
    p6 = vaddq_s16(q6, q7);
    p7 = vsubq_s16(q6, q7);
    // 12 bit

    if (_dc) *_dc = vaddlvq_s16(p0);

    // Transpose 8x8
    int16x8x2_t s0 = vtrnq_s16(p0, p1);
    int16x8x2_t s1 = vtrnq_s16(p2, p3);
    int16x8x2_t s2 = vtrnq_s16(p4, p5);
    int16x8x2_t s3 = vtrnq_s16(p6, p7);
    /*  1  9  3 11  5 13  7 15  s0[0] = VTRN1.1N p0 p1 */
    /*  2 10  4 12  6 14  8 16  s0[1] = VTRN2.1N p0 p1 */
    /* 17 25 19 27 21 29 23 31  s1[0] = VTRN1.1N p2 p3 */
    /* 18 26 20 28 22 30 24 32  s1[1] = VTRN2.1N p2 p3 */
    /* 33 41 35 43 37 45 39 47  s2[0] = VTRN1.1N p4 p5 */
    /* 34 42 36 44 38 46 40 48  s2[1] = VTRN2.1N p4 p5 */
    /* 49 57 51 59 53 61 55 63  s3[0] = VTRN1.1N p6 p7 */
    /* 50 58 52 60 54 62 56 64  s3[1] = VTRN2.1N p6 p7 */

    int32x4x2_t t0 = vtrnq_s32(vreinterpretq_s32_s16(s0.val[0]), vreinterpretq_s32_s16(s1.val[0]));
    int32x4x2_t t1 = vtrnq_s32(vreinterpretq_s32_s16(s0.val[1]), vreinterpretq_s32_s16(s1.val[1]));
    int32x4x2_t t2 = vtrnq_s32(vreinterpretq_s32_s16(s2.val[0]), vreinterpretq_s32_s16(s3.val[0]));
    int32x4x2_t t3 = vtrnq_s32(vreinterpretq_s32_s16(s2.val[1]), vreinterpretq_s32_s16(s3.val[1]));
    /*  1  9 17 25  5 13 21 29  t0[0] = VTRN1.2N s0[0] s1[0] */
    /*  2 10 18 26  6 14 22 30  t1[0] = VTRN1.2N s0[1] s1[1] */
    /*  3 11 19 27  7 15 23 31  t0[1] = VTRN2.2N s0[0] s1[0] */
    /*  4 12 20 28  8 16 24 32  t1[1] = VTRN2.2N s0[1] s1[1] */
    /* 33 41 49 57 37 45 53 61  t2[0] = VTRN1.2N s2[0] s3[0] */
    /* 34 42 50 58 38 46 54 62  t3[0] = VTRN1.2N s2[1] s3[1] */
    /* 35 43 51 59 39 47 55 63  t2[1] = VTRN2.2N s2[0] s3[0] */
    /* 36 44 52 60 40 48 56 64  t3[1] = VTRN2.2N s2[1] s3[1] */

    int64x2_t u0 = vuzp1q_s64(vreinterpretq_s64_s32(t0.val[0]), vreinterpretq_s64_s32(t2.val[0]));
    int64x2_t u1 = vuzp1q_s64(vreinterpretq_s64_s32(t1.val[0]), vreinterpretq_s64_s32(t3.val[0]));
    int64x2_t u2 = vuzp1q_s64(vreinterpretq_s64_s32(t0.val[1]), vreinterpretq_s64_s32(t2.val[1]));
    int64x2_t u3 = vuzp1q_s64(vreinterpretq_s64_s32(t1.val[1]), vreinterpretq_s64_s32(t3.val[1]));
    int64x2_t u4 = vuzp2q_s64(vreinterpretq_s64_s32(t0.val[0]), vreinterpretq_s64_s32(t2.val[0]));
    int64x2_t u5 = vuzp2q_s64(vreinterpretq_s64_s32(t1.val[0]), vreinterpretq_s64_s32(t3.val[0]));
    int64x2_t u6 = vuzp2q_s64(vreinterpretq_s64_s32(t0.val[1]), vreinterpretq_s64_s32(t2.val[1]));
    int64x2_t u7 = vuzp2q_s64(vreinterpretq_s64_s32(t1.val[1]), vreinterpretq_s64_s32(t3.val[1]));
    /*  1  9 17 25 33 41 49 57  u0 = VUZP1 t0[0] t2[0] */
    /*  2 10 18 26 34 42 50 58  u1 = VUZP1 t1[0] t3[0] */
    /*  3 11 19 27 35 43 51 59  u2 = VUZP1 t0[1] t2[1] */
    /*  4 12 20 28 36 44 52 60  u3 = VUZP1 t1[1] t3[1] */
    /*  5 13 21 29 37 45 53 61  u4 = VUZP2 t0[0] t2[0] */
    /*  6 14 22 30 38 46 54 62  u5 = VUZP2 t1[0] t3[0] */
    /*  7 15 23 31 39 47 55 63  u6 = VUZP2 t0[1] t2[1] */
    /*  8 16 24 32 40 48 56 64  u7 = VUZP2 t1[1] t3[1] */

    // Hadamard stage 1: 4x2
    p0 = vaddq_s16(vreinterpretq_s16_s64(u0), vreinterpretq_s16_s64(u4));
    p1 = vaddq_s16(vreinterpretq_s16_s64(u1), vreinterpretq_s16_s64(u5));
    p2 = vaddq_s16(vreinterpretq_s16_s64(u2), vreinterpretq_s16_s64(u6));
    p3 = vaddq_s16(vreinterpretq_s16_s64(u3), vreinterpretq_s16_s64(u7));
    p4 = vsubq_s16(vreinterpretq_s16_s64(u0), vreinterpretq_s16_s64(u4));
    p5 = vsubq_s16(vreinterpretq_s16_s64(u1), vreinterpretq_s16_s64(u5));
    p6 = vsubq_s16(vreinterpretq_s16_s64(u2), vreinterpretq_s16_s64(u6));
    p7 = vsubq_s16(vreinterpretq_s16_s64(u3), vreinterpretq_s16_s64(u7));
    // 13 bit

    // Hadamard stage 2: 2x4
    q0 = vaddq_s16(p0, p2);
    q1 = vaddq_s16(p1, p3);
    q2 = vsubq_s16(p0, p2);
    q3 = vsubq_s16(p1, p3);
    q4 = vaddq_s16(p4, p6);
    q5 = vaddq_s16(p5, p7);
    q6 = vsubq_s16(p4, p6);
    q7 = vsubq_s16(p5, p7);
    // 14 bit

    // Hadamard stage 3: 1x8, ABS
    p0 = vabsq_s16(vaddq_s16(q0, q1));
    p2 = vabsq_s16(vaddq_s16(q2, q3));
    p4 = vabsq_s16(vaddq_s16(q4, q5));
    p6 = vabsq_s16(vaddq_s16(q6, q7));
    // 15->14 bit (ABS)

    // SUM
    uint16x8_t v0 = vreinterpretq_u16_s16(vaddq_s16(p0, p4));
    uint16x8_t v2 = vreinterpretq_u16_s16(vaddq_s16(p2, p6));
    // 15 bit

    v0 = vaddq_u16(v0, v2);
    v2 = vreinterpretq_u16_s16(vabaq_s16(vabaq_s16(vabaq_s16(vabdq_s16(q0, q1), q2, q3), q4, q5), q6, q7));
    // 16 bit full

    return vaddlvq_u16(v0) + vaddlvq_u16(v2);
}

unsigned oc_enc_frag_satd_neon(const unsigned char *_src,
                               const unsigned char *_ref,int _ystride,unsigned _thresh){
    int16x8_t buf[8];
    for (int i = 0; i < 8; i++)
        buf[i] = vreinterpretq_s16_u16(vsubl_u8(
                vld1_u8(&_src[i*_ystride]),
                vld1_u8(&_ref[i*_ystride])));

    return oc_hadamard_satd_neon_8x8(NULL,buf);
}

unsigned oc_enc_frag_satd2_neon(const unsigned char *_src,
                                const unsigned char *_ref1,const unsigned char *_ref2,int _ystride,unsigned _thresh){
    int16x8_t buf[8];
    for (int i = 0; i < 8; i++)
        buf[i] = vreinterpretq_s16_u16(vsubq_u16(
                vmovl_u8(vld1_u8(&_src[i*_ystride])),
                vshrq_n_u16(vaddl_u8(vld1_u8(&_ref1[i*_ystride]), vld1_u8(&_ref2[i*_ystride])), 1)));

    return oc_hadamard_satd_neon_8x8(NULL,buf);
}

unsigned oc_enc_frag_intra_satd_neon(const unsigned char *_src,int _ystride){
    int16x8_t buf[8];
    int dc;
    for (int i = 0; i < 8; i++)
        buf[i] = vreinterpretq_s16_u16(vmovl_u8(vld1_u8(&_src[i*_ystride])));

    unsigned sad = oc_hadamard_satd_neon_8x8(&dc,buf);
    return sad - dc;
}

unsigned oc_enc_frag_sad_neon(const unsigned char *_src,
                              const unsigned char *_ref,int _ystride){
    uint16x8_t a;

    for(int i = 0; i < 8; i++) {
        uint8x8_t src = vld1_u8(&_src[i * _ystride]);
        uint8x8_t ref = vld1_u8(&_ref[i * _ystride]);
        a = i ? vabal_u8(a, src, ref) : vabdl_u8(src, ref);
    }

    return vaddlvq_u16(a);
}

unsigned oc_enc_frag_sad_thresh_neon(const unsigned char *_src,
                                     const unsigned char *_ref,int _ystride,unsigned _thresh){
    return oc_enc_frag_sad_neon(_src, _ref, _ystride);
}

unsigned oc_enc_frag_sad2_thresh_neon(const unsigned char *_src,
                                      const unsigned char *_ref1, const unsigned char *_ref2,
                                      int _ystride, unsigned _thresh){
    uint16x8_t a;

    for(int i = 0; i < 8; i++) {
        uint16x8_t src = vmovl_u8(vld1_u8(&_src[i * _ystride]));
        uint16x8_t ref = vshrq_n_u16(vaddl_u8(vld1_u8(&_ref1[i * _ystride]), vld1_u8(&_ref2[i * _ystride])), 1);
        a = i ? vabaq_u16(a, src, ref) : vabdq_u16(src, ref);
    }

    return vaddlvq_u16(a);
}
