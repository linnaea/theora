#include "wav128frag.h"

#ifdef OC_WASM_SIMD128
#include <wasm_simd128.h>

void oc_frag_recon_intra_wav128(unsigned char *_dst,
                                int _ystride, const int16_t _residue[64]) {
    v128_t ref = wasm_i16x8_const_splat(128);

    for (int i = 0; i < 4; i++) {
        v128_t res0 = wasm_v128_load(&_residue[(i * 2 + 0) * 8]);
        v128_t frm0 = wasm_i16x8_add(res0, ref);
        v128_t res1 = wasm_v128_load(&_residue[(i * 2 + 1) * 8]);
        v128_t frm1 = wasm_i16x8_add(res1, ref);
        v128_t frm = wasm_u8x16_narrow_i16x8(frm0, frm1);
        wasm_v128_store64_lane(&_dst[(i * 2 + 0) * _ystride], frm, 0);
        wasm_v128_store64_lane(&_dst[(i * 2 + 1) * _ystride], frm, 1);
    }
}

void oc_frag_recon_inter_wav128(unsigned char *_dst, const unsigned char *_src,
                                int _ystride,const int16_t _residue[64]) {
    for (int i = 0; i < 4; i++) {
        v128_t res0 = wasm_v128_load(&_residue[(i * 2 + 0) * 8]);
        v128_t res1 = wasm_v128_load(&_residue[(i * 2 + 1) * 8]);
        v128_t ref0 = wasm_u16x8_load8x8(&_src[(i * 2 + 0) * _ystride]);
        v128_t ref1 = wasm_u16x8_load8x8(&_src[(i * 2 + 1) * _ystride]);
        v128_t frm0 = wasm_i16x8_add(res0, ref0);
        v128_t frm1 = wasm_i16x8_add(res1, ref1);
        v128_t frm = wasm_u8x16_narrow_i16x8(frm0, frm1);
        wasm_v128_store64_lane(&_dst[(i * 2 + 0) * _ystride], frm, 0);
        wasm_v128_store64_lane(&_dst[(i * 2 + 1) * _ystride], frm, 1);
    }
}

void oc_frag_recon_inter2_wav128(unsigned char *_dst,const unsigned char *_src1, const unsigned char *_src2,
                                 int _ystride,const int16_t _residue[64]) {
    for (int i = 0; i < 4; i++) {
        v128_t res0 = wasm_v128_load(&_residue[(i * 2 + 0) * 8]);
        v128_t res1 = wasm_v128_load(&_residue[(i * 2 + 1) * 8]);
        v128_t ref0 = wasm_i16x8_add(wasm_u16x8_load8x8(&_src1[(i * 2 + 0) * _ystride]),
                                     wasm_u16x8_load8x8(&_src2[(i * 2 + 0) * _ystride]));
        v128_t ref1 = wasm_i16x8_add(wasm_u16x8_load8x8(&_src1[(i * 2 + 1) * _ystride]),
                                     wasm_u16x8_load8x8(&_src2[(i * 2 + 1) * _ystride]));
        v128_t frm0 = wasm_i16x8_add(res0, wasm_u16x8_shr(ref0, 1));
        v128_t frm1 = wasm_i16x8_add(res1, wasm_u16x8_shr(ref1, 1));
        v128_t frm = wasm_u8x16_narrow_i16x8(frm0, frm1);
        wasm_v128_store64_lane(&_dst[(i * 2 + 0) * _ystride], frm, 0);
        wasm_v128_store64_lane(&_dst[(i * 2 + 1) * _ystride], frm, 1);
    }
}

#endif
