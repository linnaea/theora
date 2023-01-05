#ifndef _arm_intrinsics_armfrag_H
#define _arm_intrinsics_armfrag_H 1
#include <stdint.h>
#include <stddef.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef OC_ARM_ASM_NEON
void oc_idct8x8_neon(int16_t _y[64],int16_t _x[64],int _last_zzi);
void oc_frag_copy_neon(unsigned char *_dst,const unsigned char *_src,int _ystride);
void oc_frag_copy_list_neon(unsigned char *_dst_frame, const unsigned char *_src_frame,int _ystride,
                            const ptrdiff_t *_fragis,ptrdiff_t _nfragis,const ptrdiff_t *_frag_buf_offs);
void oc_frag_recon_intra_neon(unsigned char *_dst,
                              int _ystride, const int16_t _residue[64]);
void oc_frag_recon_inter_neon(unsigned char *_dst, const unsigned char *_src,
                              int _ystride,const int16_t _residue[64]);
void oc_frag_recon_inter2_neon(unsigned char *_dst,const unsigned char *_src1, const unsigned char *_src2,
                               int _ystride,const int16_t _residue[64]);
#endif

#endif //_arm_intrinsics_armfrag_H
