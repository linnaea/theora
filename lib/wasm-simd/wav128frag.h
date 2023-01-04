#ifndef _wasm_simd_wav128frag_H
#define _wasm_simd_wav128frag_H 1
#include <stdint.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef OC_WASM_SIMD128
void oc_frag_recon_intra_wav128(unsigned char *_dst,
                                int _ystride, const int16_t _residue[64]);
void oc_frag_recon_inter_wav128(unsigned char *_dst, const unsigned char *_src,
                                int _ystride,const int16_t _residue[64]);
void oc_frag_recon_inter2_wav128(unsigned char *_dst,const unsigned char *_src1, const unsigned char *_src2,
                                 int _ystride,const int16_t _residue[64]);
#endif

#endif //_wasm_simd_wav128frag_H
