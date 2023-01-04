#ifndef _wasm_simd_wav128int_H
#define _wasm_simd_wav128int_H 1

# ifdef OC_WASM_SIMD128
#  include "wav128frag.h"
#  define oc_frag_recon_intra(_enc,...) oc_frag_recon_intra_wav128(__VA_ARGS__)
#  define oc_frag_recon_inter(_enc,...) oc_frag_recon_inter_wav128(__VA_ARGS__)
#  define oc_frag_recon_inter2(_enc,...) oc_frag_recon_inter2_wav128(__VA_ARGS__)
# endif

#endif //_wasm_simd_wav128int_H
