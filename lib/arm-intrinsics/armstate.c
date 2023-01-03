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
    last mod: $Id: x86state.c 17344 2010-07-21 01:42:18Z tterribe $

 ********************************************************************/
#include "armint.h"

#if defined(OC_ARM_ASM)

void oc_state_accel_init_arm(oc_theora_state *_state){
  oc_state_accel_init_c(_state);
  _state->cpu_flags=oc_cpu_flags_get();
# if defined(OC_ENC_USE_VTABLE)
#  if defined(OC_ARM_ASM_NEON)
  if(cpu_flags & OC_CPU_ARM_NEON){
    _state->opt_vtable.frag_copy=oc_frag_copy_neon;
    _state->opt_vtable.frag_copy_list=oc_frag_copy_list_neon;
    _state->opt_vtable.frag_recon_intra=oc_frag_recon_intra_neon;
    _state->opt_vtable.frag_recon_inter=oc_frag_recon_inter_neon;
    _state->opt_vtable.frag_recon_inter2=oc_frag_recon_inter2_neon;
  }
#  endif
# endif
}

#endif
