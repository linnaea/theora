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
#if !defined(_arm_armint_H)
# define _arm_armint_H (1)
# include "../internal.h"

# if defined(OC_ARM_ASM)
#  define oc_state_accel_init oc_state_accel_init_arm
#  if defined(__aarch64__) || defined(_M_ARM64)
#   define oc_idct8x8(_state,...) oc_idct8x8_neon(__VA_ARGS__)
#   define oc_frag_copy(_state,...) oc_frag_copy_neon(__VA_ARGS__)
#   define oc_frag_copy_list(_enc,...) oc_frag_copy_list_neon(__VA_ARGS__)
#   define oc_frag_recon_intra(_enc,...) oc_frag_recon_intra_neon(__VA_ARGS__)
#   define oc_frag_recon_inter(_enc,...) oc_frag_recon_inter_neon(__VA_ARGS__)
#   define oc_frag_recon_inter2(_enc,...) oc_frag_recon_inter2_neon(__VA_ARGS__)
#  else
#   define OC_STATE_USE_VTABLE (1)
#  endif
# endif

# include "../state.h"
# include "armcpu.h"

# if defined(OC_ARM_ASM)
#  include "armfrag.h"
void oc_state_accel_init_arm(oc_theora_state *_state);
# endif

#endif