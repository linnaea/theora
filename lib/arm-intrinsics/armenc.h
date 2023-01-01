/********************************************************************
 *                                                                  *
 * THIS FILE IS PART OF THE OggTheora SOFTWARE CODEC SOURCE CODE.   *
 * USE, DISTRIBUTION AND REPRODUCTION OF THIS LIBRARY SOURCE IS     *
 * GOVERNED BY A BSD-STYLE SOURCE LICENSE INCLUDED WITH THIS SOURCE *
 * IN 'COPYING'. PLEASE READ THESE TERMS BEFORE DISTRIBUTING.       *
 *                                                                  *
 * THE Theora SOURCE CODE IS COPYRIGHT (C) 2002-2010                *
 * by the Xiph.Org Foundation and contributors http://www.xiph.org/ *
 *                                                                  *
 ********************************************************************

  function:
    last mod: $Id: x86int.h 17344 2010-07-21 01:42:18Z tterribe $

 ********************************************************************/
#if !defined(_arm_armenc_H)
# define _arm_armenc_H (1)
# include "armint.h"

# if defined(OC_ARM_ASM)
#  define oc_enc_vtable_init oc_enc_accel_init_arm
#  define OC_ENC_USE_VTABLE (1)
# endif

# include "../encint.h"

# if defined(OC_ARM_ASM)
void oc_enc_accel_init_arm(oc_enc_ctx *_enc);

#  if defined(OC_ARM_ASM_NEON)
unsigned oc_enc_frag_satd_neon(const unsigned char *_src,
                               const unsigned char *_ref,int _ystride,unsigned _thresh);
unsigned oc_enc_frag_satd2_neon(const unsigned char *_src,
                                const unsigned char *_ref1,const unsigned char *_ref2,int _ystride,unsigned _thresh);
unsigned oc_enc_frag_intra_satd_neon(const unsigned char *_src,int _ystride);
unsigned oc_enc_frag_sad_neon(const unsigned char *_src,
                              const unsigned char *_ref,int _ystride);
unsigned oc_enc_frag_sad_thresh_neon(const unsigned char *_src,
                                     const unsigned char *_ref,int _ystride,unsigned _thresh);
unsigned oc_enc_frag_sad2_thresh_neon(const unsigned char *_src,
                                      const unsigned char *_ref1, const unsigned char *_ref2,
                                      int _ystride, unsigned _thresh);
#  endif

# endif

#endif
