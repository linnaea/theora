#define _GNU_SOURCE
#include "enc_worker.h"
#include <string.h>
#include <stdatomic.h>
#include <stdio.h>
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

static void oc_cost_inter_nomv_dct(oc_enc_ctx *_enc,unsigned _mbi,int _mb_mode, oc_dct_cost_table *_dct){
  oc_cost_inter_dct(_enc,_mbi,_mb_mode,0,_dct);
}

static void oc_cost_inter1mv_dct(oc_enc_ctx *_enc,unsigned _mbi,int _mb_mode,oc_mv _mv, oc_dct_cost_table *_dct){
  oc_cost_inter_dct(_enc,_mbi,_mb_mode,_mv,_dct);
}

/*Estimate the R-D cost of the DCT coefficients given the SATD of a block after
   prediction.*/
static unsigned oc_dct_cost2(oc_enc_ctx *_enc,unsigned *_ssd,
 int _qii,int _pli,int _qti,unsigned _satd,int shift){
  unsigned rmse;
  int      bin;
  int      dx;
  int      y0;
  int      z0;
  int      dy;
  int      dz;
  /*SATD metrics for chroma planes vary much less than luma, so we scale them
     by 4 to distribute them into the mode decision bins more evenly.*/
  _satd<<=_pli+1&2;
  bin=OC_MINI(_satd>>shift,OC_COMP_BINS-2);
  dx=_satd-(bin<<shift);
  y0=_enc->mode_rd[_qii][_pli][_qti][bin].rate;
  z0=_enc->mode_rd[_qii][_pli][_qti][bin].rmse;
  dy=_enc->mode_rd[_qii][_pli][_qti][bin+1].rate-y0;
  dz=_enc->mode_rd[_qii][_pli][_qti][bin+1].rmse-z0;
  rmse=OC_MAXI(z0+(dz*dx>>shift),0);
  *_ssd=rmse*rmse>>2*OC_RMSE_SCALE-OC_BIT_SCALE;
  return OC_MAXI(y0+(dy*dx>>shift),0);
}

void oc_cost_dct_fill(oc_enc_ctx *_enc, int nqis, int qti, int pli, int bi, unsigned satd, int shift, oc_dct_cost_table *_dct) {
  if(!shift) shift=_enc->sp_level<OC_SP_LEVEL_NOSATD?OC_SATD_SHIFT:OC_SAD_SHIFT;
  for(int qii=0;qii<nqis;qii++)
    _dct->block[bi][qii].cost = oc_dct_cost2(_enc, &_dct->block[bi][qii].ssd, qii, pli, qti, satd, shift);
}

static void oc_skip_cost_ssd(oc_enc_ctx *_enc,
 unsigned _mbi,const unsigned _rd_scale[5],unsigned _ssd[12]){
  const unsigned char   *src;
  const unsigned char   *ref;
  int                    ystride;
  const oc_fragment     *frags;
  const ptrdiff_t       *frag_buf_offs;
  const ptrdiff_t       *sb_map;
  const oc_mb_map_plane *mb_map;
  const unsigned char   *map_idxs;
  oc_mv                 *mvs;
  int                    map_nidxs;
  unsigned               uncoded_ssd;
  int                    mapii;
  int                    mapi;
  int                    pli;
  int                    bi;
  ptrdiff_t              fragi;
  ptrdiff_t              frag_offs;
  int                    borderi;
  src=_enc->state.ref_frame_data[OC_FRAME_IO];
  ref=_enc->state.ref_frame_data[OC_FRAME_PREV];
  ystride=_enc->state.ref_ystride[0];
  frags=_enc->state.frags;
  frag_buf_offs=_enc->state.frag_buf_offs;
  sb_map=_enc->state.sb_maps[_mbi>>2][_mbi&3];
  mvs=_enc->mb_info[_mbi].block_mv;
  for(bi=0;bi<4;bi++){
    fragi=sb_map[bi];
    borderi=frags[fragi].borderi;
    frag_offs=frag_buf_offs[fragi];
    if(borderi<0){
      uncoded_ssd=oc_enc_frag_ssd(_enc,src+frag_offs,ref+frag_offs,ystride);
    }
    else{
      uncoded_ssd=oc_enc_frag_border_ssd(_enc,
       src+frag_offs,ref+frag_offs,ystride,_enc->state.borders[borderi].mask);
    }
    /*Scale to match DCT domain and RD.*/
    uncoded_ssd=OC_RD_SKIP_SCALE(uncoded_ssd,_rd_scale[bi]);
    /*Motion is a special case; if there is more than a full-pixel motion
       against the prior frame, penalize skipping.
      TODO: The factor of two here is a kludge, but it tested out better than a
       hard limit.*/
    if(mvs[bi]!=0)uncoded_ssd*=2;
    _ssd[bi]=uncoded_ssd;
  }
  mb_map=(const oc_mb_map_plane *)_enc->state.mb_maps[_mbi];
  map_nidxs=OC_MB_MAP_NIDXS[_enc->state.info.pixel_fmt];
  map_idxs=OC_MB_MAP_IDXS[_enc->state.info.pixel_fmt];
  map_nidxs=(map_nidxs-4>>1)+4;
  mapii=4;
  mvs=_enc->mb_info[_mbi].unref_mv;
  for(pli=1;pli<3;pli++){
    ystride=_enc->state.ref_ystride[pli];
    for(;mapii<map_nidxs;mapii++){
      mapi=map_idxs[mapii];
      bi=mapi&3;
      fragi=mb_map[pli][bi];
      borderi=frags[fragi].borderi;
      frag_offs=frag_buf_offs[fragi];
      if(borderi<0){
        uncoded_ssd=oc_enc_frag_ssd(_enc,src+frag_offs,ref+frag_offs,ystride);
      }
      else{
        uncoded_ssd=oc_enc_frag_border_ssd(_enc,
         src+frag_offs,ref+frag_offs,ystride,_enc->state.borders[borderi].mask);
      }
      /*Scale to match DCT domain and RD.*/
      uncoded_ssd=OC_RD_SKIP_SCALE(uncoded_ssd,_rd_scale[4]);
      /*Motion is a special case; if there is more than a full-pixel motion
         against the prior frame, penalize skipping.
        TODO: The factor of two here is a kludge, but it tested out better than
         a hard limit*/
      if(mvs[OC_FRAME_PREV]!=0)uncoded_ssd*=2;
      _ssd[mapii]=uncoded_ssd;
    }
    map_nidxs=(map_nidxs-4<<1)+4;
  }
}

static unsigned oc_mb_activity(oc_enc_ctx *_enc,unsigned _mbi,
 unsigned _activity[4]){
  const unsigned char *src;
  const ptrdiff_t     *frag_buf_offs;
  const ptrdiff_t     *sb_map;
  unsigned             luma;
  int                  ystride;
  ptrdiff_t            frag_offs;
  ptrdiff_t            fragi;
  int                  bi;
  frag_buf_offs=_enc->state.frag_buf_offs;
  sb_map=_enc->state.sb_maps[_mbi>>2][_mbi&3];
  src=_enc->state.ref_frame_data[OC_FRAME_IO];
  ystride=_enc->state.ref_ystride[0];
  luma=0;
  for(bi=0;bi<4;bi++){
    const unsigned char *s;
    unsigned             x;
    unsigned             x2;
    unsigned             act;
    int                  i;
    int                  j;
    fragi=sb_map[bi];
    frag_offs=frag_buf_offs[fragi];
    /*TODO: This could be replaced with SATD^2, since we already have to
       compute SATD.*/
    x=x2=0;
    s=src+frag_offs;
    for(i=0;i<8;i++){
      for(j=0;j<8;j++){
        unsigned c;
        c=s[j];
        x+=c;
        x2+=c*c;
      }
      s+=ystride;
    }
    luma+=x;
    act=(x2<<6)-x*x;
    if(act<8<<12){
      /*The region is flat.*/
      act=OC_MINI(act,5<<12);
    }
    else{
      unsigned e1;
      unsigned e2;
      unsigned e3;
      unsigned e4;
      /*Test for an edge.
        TODO: There are probably much simpler ways to do this (e.g., it could
         probably be combined with the SATD calculation).
        Alternatively, we could split the block around the mean and compute the
         reduction in variance in each half.
        For a Gaussian source the reduction should be
         (1-2/pi) ~= 0.36338022763241865692446494650994.
        Significantly more reduction is a good indication of a bi-level image.
        This has the advantage of identifying, in addition to straight edges,
         small text regions, which would otherwise be classified as "texture".*/
      e1=e2=e3=e4=0;
      s=src+frag_offs-1;
      for(i=0;i<8;i++){
        for(j=0;j<8;j++){
          e1+=abs((s[j+2]-s[j]<<1)+(s-ystride)[j+2]-(s-ystride)[j]
           +(s+ystride)[j+2]-(s+ystride)[j]);
          e2+=abs(((s+ystride)[j+1]-(s-ystride)[j+1]<<1)
           +(s+ystride)[j]-(s-ystride)[j]+(s+ystride)[j+2]-(s-ystride)[j+2]);
          e3+=abs(((s+ystride)[j+2]-(s-ystride)[j]<<1)
           +(s+ystride)[j+1]-s[j]+s[j+2]-(s-ystride)[j+1]);
          e4+=abs(((s+ystride)[j]-(s-ystride)[j+2]<<1)
           +(s+ystride)[j+1]-s[j+2]+s[j]-(s-ystride)[j+1]);
        }
        s+=ystride;
      }
      /*If the largest component of the edge energy is at least 40% of the
         total, then classify the block as an edge block.*/
      if(5*OC_MAXI(OC_MAXI(e1,e2),OC_MAXI(e3,e4))>2*(e1+e2+e3+e4)){
         /*act=act_th*(act/act_th)**0.7
              =exp(log(act_th)+0.7*(log(act)-log(act_th))).
           Here act_th=5.0 and 0x394A=oc_blog32_q10(5<<12).*/
         act=oc_bexp32_q10(0x394A+(7*(oc_blog32_q10(act)-0x394A+5)/10));
      }
    }
    _activity[bi]=act;
  }
  return luma;
}

static void oc_mb_activity_fast(oc_enc_ctx *_enc,unsigned _mbi,
 unsigned _activity[4],const unsigned _intra_satd[12]){
  int bi;
  for(bi=0;bi<4;bi++){
    unsigned act;
    act=(11*_intra_satd[bi]>>8)*_intra_satd[bi];
    if(act<8<<12){
      /*The region is flat.*/
      act=OC_MINI(act,5<<12);
    }
    _activity[bi]=act;
  }
}

/*Compute the masking scales for the blocks in a macro block.
  All masking is computed from the luma blocks.
  We derive scaling factors for the chroma blocks from these, and use the same
   ones for all chroma blocks, regardless of the subsampling.
  It's possible for luma to be perfectly flat and yet have high chroma energy,
   but this is unlikely in non-artificial images, and not a case that has been
   addressed by any research to my knowledge.
  The output of the masking process is two scale factors, which are fed into
   the various R-D optimizations.
  The first, rd_scale, is applied to D in the equation
    D*rd_scale+lambda*R.
  This is the form that must be used to properly combine scores from multiple
   blocks, and can be interpreted as scaling distortions by their visibility.
  The inverse, rd_iscale, is applied to lambda in the equation
    D+rd_iscale*lambda*R.
  This is equivalent to the first form within a single block, but much faster
   to use when evaluating many possible distortions (e.g., during actual
   quantization, where separate distortions are evaluated for every
   coefficient).
  The two macros OC_RD_SCALE(rd_scale,d) and OC_RD_ISCALE(rd_iscale,lambda) are
   used to perform the multiplications with the proper re-scaling for the range
   of the scaling factors.
  Many researchers apply masking values directly to the quantizers used, and
   not to the R-D cost.
  Since we generally use MSE for D, rd_scale must use the square of their
   values to generate an equivalent effect.*/
static unsigned oc_mb_masking(unsigned _rd_scale[5],unsigned _rd_iscale[5],
 const ogg_uint16_t _chroma_rd_scale[2],const unsigned _activity[4],
 unsigned _activity_avg,unsigned _luma,unsigned _luma_avg){
  unsigned activity_sum;
  unsigned la;
  unsigned lb;
  unsigned d;
  int      bi;
  int      bi_min;
  int      bi_min2;
  /*The ratio lb/la is meant to approximate
     ((((_luma-16)/219)*(255/128))**0.649**0.4**2), which is the
     effective luminance masking from~\cite{LKW06} (including the self-masking
     deflator).
    The following actually turns out to be a pretty good approximation for
     _luma>75 or so.
    For smaller values luminance does not really follow Weber's Law anyway, and
     this approximation gives a much less aggressive bitrate boost in this
     region.
    Though some researchers claim that contrast sensitivity actually decreases
     for very low luminance values, in my experience excessive brightness on
     LCDs or buggy color conversions (e.g., treating Y' as full-range instead
     of the CCIR 601 range) make artifacts in such regions extremely visible.
    We substitute _luma_avg for 128 to allow the strength of the masking to
     vary with the actual average image luminance, within certain limits (the
     caller has clamped _luma_avg to the range [90,160], inclusive).
    @ARTICLE{LKW06,
      author="Zhen Liu and Lina J. Karam and Andrew B. Watson",
      title="{JPEG2000} Encoding With Perceptual Distortion Control",
      journal="{IEEE} Transactions on Image Processing",
      volume=15,
      number=7,
      pages="1763--1778",
      month=Jul,
      year=2006
    }*/
#if 0
  la=_luma+4*_luma_avg;
  lb=4*_luma+_luma_avg;
#else
  /*Disable luminance masking.*/
  la=lb=1;
#endif
  activity_sum=0;
  for(bi=0;bi<4;bi++){
    unsigned a;
    unsigned b;
    activity_sum+=_activity[bi];
    /*Apply activity masking.*/
    a=_activity[bi]+4*_activity_avg;
    b=4*_activity[bi]+_activity_avg;
    d=OC_RD_SCALE(b,1);
    /*And luminance masking.*/
    d=(a+(d>>1))/d;
    _rd_scale[bi]=(d*la+(lb>>1))/lb;
    /*And now the inverse.*/
    d=OC_MAXI(OC_RD_ISCALE(a,1),1);
    d=(b+(d>>1))/d;
    _rd_iscale[bi]=(d*lb+(la>>1))/la;
  }
  /*Now compute scaling factors for chroma blocks.
    We start by finding the two smallest iscales from the luma blocks.*/
  bi_min=_rd_iscale[1]<_rd_iscale[0];
  bi_min2=1-bi_min;
  for(bi=2;bi<4;bi++){
    if(_rd_iscale[bi]<_rd_iscale[bi_min]){
      bi_min2=bi_min;
      bi_min=bi;
    }
    else if(_rd_iscale[bi]<_rd_iscale[bi_min2])bi_min2=bi;
  }
  /*If the minimum iscale is less than 1.0, use the second smallest instead,
     and force the value to at least 1.0 (inflating chroma is a waste).*/
  if(_rd_iscale[bi_min]<(1<<OC_RD_ISCALE_BITS))bi_min=bi_min2;
  d=OC_MINI(_rd_scale[bi_min],1<<OC_RD_SCALE_BITS);
  _rd_scale[4]=OC_RD_SCALE(d,_chroma_rd_scale[0]);
  d=OC_MAXI(_rd_iscale[bi_min],1<<OC_RD_ISCALE_BITS);
  _rd_iscale[4]=OC_RD_ISCALE(d,_chroma_rd_scale[1]);
  return activity_sum;
}

static unsigned oc_mb_intra_dct(oc_enc_ctx *_enc,unsigned _mbi,
                                unsigned _intra_satd[12], oc_dct_cost_table *_dct){
  const unsigned char   *src;
  const ptrdiff_t       *frag_buf_offs;
  const ptrdiff_t       *sb_map;
  const oc_mb_map_plane *mb_map;
  const unsigned char   *map_idxs;
  int                    map_nidxs;
  int                    mapii;
  int                    mapi;
  int                    ystride;
  int                    pli;
  int                    bi;
  ptrdiff_t              fragi;
  ptrdiff_t              frag_offs;
  unsigned               luma;
  int                    dc;
  int                    nqis;
  nqis=_enc->state.nqis;
  frag_buf_offs=_enc->state.frag_buf_offs;
  sb_map=_enc->state.sb_maps[_mbi>>2][_mbi&3];
  src=_enc->state.ref_frame_data[OC_FRAME_IO];
  ystride=_enc->state.ref_ystride[0];
  luma=0;
  for(bi=0;bi<4;bi++){
    fragi=sb_map[bi];
    frag_offs=frag_buf_offs[fragi];
    _intra_satd[bi]=oc_enc_frag_intra_satd(_enc,&dc,src+frag_offs,ystride);
    _intra_satd[bi]+=abs(dc);
    luma+=dc;
    oc_cost_dct_fill(_enc, nqis, 0, 0, bi, _intra_satd[bi], OC_SATD_SHIFT, _dct);
  }
  mb_map=(const oc_mb_map_plane *)_enc->state.mb_maps[_mbi];
  map_idxs=OC_MB_MAP_IDXS[_enc->state.info.pixel_fmt];
  map_nidxs=OC_MB_MAP_NIDXS[_enc->state.info.pixel_fmt];
  /*Note: This assumes ref_ystride[1]==ref_ystride[2].*/
  ystride=_enc->state.ref_ystride[1];
  nqis=1;
  for(mapii=4;mapii<map_nidxs;mapii++){
    mapi=map_idxs[mapii];
    pli=mapi>>2;
    bi=mapi&3;
    fragi=mb_map[pli][bi];
    frag_offs=frag_buf_offs[fragi];
    if(_enc->sp_level<OC_SP_LEVEL_NOSATD){
      _intra_satd[mapii]=oc_enc_frag_intra_satd(_enc,&dc,src+frag_offs,ystride);
      _intra_satd[mapii]+=abs(dc);
    }
    else{
      _intra_satd[mapii]=oc_enc_frag_intra_sad(_enc,src+frag_offs,ystride);
    }
    oc_cost_dct_fill(_enc, nqis, 0, pli, mapii, _intra_satd[mapii], 0, _dct);
  }
  return luma;
}

void oc_cost_inter_dct(oc_enc_ctx *_enc, unsigned _mbi, int _mb_mode, oc_mv _mv, oc_dct_cost_table *_dct){
  const unsigned char   *src;
  const unsigned char   *ref;
  int                    ystride;
  const ptrdiff_t       *frag_buf_offs;
  const ptrdiff_t       *sb_map;
  const oc_mb_map_plane *mb_map;
  const unsigned char   *map_idxs;
  int                    map_nidxs;
  int                    mapii;
  int                    mapi;
  int                    mv_offs[2];
  int                    pli;
  int                    bi;
  ptrdiff_t              fragi;
  ptrdiff_t              frag_offs;
  int                    dc;
  unsigned               satd;
  int                    nqis;
  nqis=_enc->state.nqis;
  src=_enc->state.ref_frame_data[OC_FRAME_IO];
  ref=_enc->state.ref_frame_data[OC_FRAME_FOR_MODE(_mb_mode)];
  ystride=_enc->state.ref_ystride[0];
  frag_buf_offs=_enc->state.frag_buf_offs;
  sb_map=_enc->state.sb_maps[_mbi>>2][_mbi&3];
  if(oc_state_get_mv_offsets(&_enc->state,mv_offs,0,_mv)>1){
    for(bi=0;bi<4;bi++){
      fragi=sb_map[bi];
      frag_offs=frag_buf_offs[fragi];
      if(_enc->sp_level<OC_SP_LEVEL_NOSATD){
        satd=oc_enc_frag_satd2(_enc,&dc,src+frag_offs,
                               ref+frag_offs+mv_offs[0],ref+frag_offs+mv_offs[1],ystride);
        satd+=abs(dc);
      }
      else{
        satd=oc_enc_frag_sad2_thresh(_enc,src+frag_offs,
                                     ref+frag_offs+mv_offs[0],ref+frag_offs+mv_offs[1],ystride,UINT_MAX);
      }
      oc_cost_dct_fill(_enc, nqis, 1, 0, bi, satd, 0, _dct);
    }
  }
  else{
    for(bi=0;bi<4;bi++){
      fragi=sb_map[bi];
      frag_offs=frag_buf_offs[fragi];
      if(_enc->sp_level<OC_SP_LEVEL_NOSATD){
        satd=oc_enc_frag_satd(_enc,&dc,src+frag_offs,
                              ref+frag_offs+mv_offs[0],ystride);
        satd+=abs(dc);
      }
      else{
        satd=oc_enc_frag_sad(_enc,src+frag_offs,
                             ref+frag_offs+mv_offs[0],ystride);
      }
      oc_cost_dct_fill(_enc, nqis, 1, 0, bi, satd, 0, _dct);
    }
  }
  mb_map=(const oc_mb_map_plane *)_enc->state.mb_maps[_mbi];
  map_idxs=OC_MB_MAP_IDXS[_enc->state.info.pixel_fmt];
  map_nidxs=OC_MB_MAP_NIDXS[_enc->state.info.pixel_fmt];
  /*Note: This assumes ref_ystride[1]==ref_ystride[2].*/
  ystride=_enc->state.ref_ystride[1];
  nqis=1;
  if(oc_state_get_mv_offsets(&_enc->state,mv_offs,1,_mv)>1){
    for(mapii=4;mapii<map_nidxs;mapii++){
      mapi=map_idxs[mapii];
      pli=mapi>>2;
      bi=mapi&3;
      fragi=mb_map[pli][bi];
      frag_offs=frag_buf_offs[fragi];
      if(_enc->sp_level<OC_SP_LEVEL_NOSATD){
        satd=oc_enc_frag_satd2(_enc,&dc,src+frag_offs,
                               ref+frag_offs+mv_offs[0],ref+frag_offs+mv_offs[1],ystride);
        satd+=abs(dc);
      }
      else{
        satd=oc_enc_frag_sad2_thresh(_enc,src+frag_offs,
                                     ref+frag_offs+mv_offs[0],ref+frag_offs+mv_offs[1],ystride,UINT_MAX);
      }
      oc_cost_dct_fill(_enc, nqis, 1, pli, mapii, satd, 0, _dct);
    }
  }
  else{
    for(mapii=4;mapii<map_nidxs;mapii++){
      mapi=map_idxs[mapii];
      pli=mapi>>2;
      bi=mapi&3;
      fragi=mb_map[pli][bi];
      frag_offs=frag_buf_offs[fragi];
      if(_enc->sp_level<OC_SP_LEVEL_NOSATD){
        satd=oc_enc_frag_satd(_enc,&dc,src+frag_offs,
                              ref+frag_offs+mv_offs[0],ystride);
        satd+=abs(dc);
      }
      else{
        satd=oc_enc_frag_sad(_enc,src+frag_offs,
                             ref+frag_offs+mv_offs[0],ystride);
      }
      oc_cost_dct_fill(_enc, nqis, 1, pli, mapii, satd, 0, _dct);
    }
  }
}

enum oc_enc_wrk_computed {
  WRK_I1MV_EST = 1 << 0,
  WRK_G1MV_EST = 1 << 1,
  WRK_4MV_EST = 1 << 2,
  WRK_I1MV_REFINED = 1 << 3,
  WRK_G1MV_REFINED = 1 << 4,
  WRK_4MV_REFINED = 1 << 5
};

struct oc_enc_wmb_data {
  unsigned luma;
  unsigned activity;
  unsigned rd_scale[5];
  unsigned rd_iscale[5];
  unsigned skip_ssd[12];
  oc_dct_cost_table intra_dct;
  oc_dct_cost_table i0mv_dct;
  oc_dct_cost_table g0mv_dct;

  enum oc_enc_wrk_computed computed;
  oc_dct_cost_table i1mv_dct;
  oc_dct_cost_table g1mv_dct;
  oc_dct_cost_table i1mvr_dct;
  oc_dct_cost_table g1mvr_dct;
};

struct oc_enc_worker_ctrl {
  int threads;
  pthread_rwlock_t frame_state_lock;
  pthread_mutex_t worker_comp_lock;
  pthread_cond_t worker_comp_signal;
  pthread_barrier_t end_sync;
  pthread_barrier_t start_sync;

  pthread_t *workers;
  struct oc_enc_wmb_data *mbs;
  atomic_flag *started;
  unsigned char *finished;

  unsigned char recode, intra, stop, stall;
  unsigned mbi_start, mbi_end;

  atomic_uint current_mbi;
};

static void inter_analysis(oc_enc_ctx *_enc,int _recode,unsigned mbi,struct oc_enc_wmb_data *d) {
  unsigned             activity_avg;
  unsigned             luma_avg;
  const ogg_uint16_t  *chroma_rd_scale;
  oc_mb_enc_info      *embs;
  int                  sp_level;
  unsigned             sbi;
  unsigned             quadi;
  sp_level=_enc->sp_level;
  activity_avg=_enc->activity_avg;
  luma_avg=OC_CLAMPI(90<<8,_enc->luma_avg,160<<8);
  chroma_rd_scale=_enc->chroma_rd_scale[OC_INTER_FRAME][_enc->state.qis[0]];
  embs=_enc->mb_info;
  sbi=mbi>>2;
  quadi=mbi&3;

  if(_enc->state.sb_flags[sbi].quad_valid&1<<quadi){
    unsigned intra_satd[12];
    unsigned activity[4];
    unsigned luma;
    luma=oc_mb_intra_dct(_enc,mbi,intra_satd,&d->intra_dct);
    /*Activity masking.*/
    if(sp_level<OC_SP_LEVEL_FAST_ANALYSIS)
      oc_mb_activity(_enc,mbi,activity);
    else
      oc_mb_activity_fast(_enc,mbi,activity,intra_satd);

    d->luma = luma;
    d->activity=oc_mb_masking(d->rd_scale,d->rd_iscale,chroma_rd_scale,activity,activity_avg,luma,luma_avg);

    /*Motion estimation:
      We always do a basic 1MV search for all macroblocks, coded or not,
       keyframe or not.*/
    if(!_recode&&sp_level<OC_SP_LEVEL_NOMC)oc_mcenc_search(_enc,mbi);

    /*Find the block choice with the lowest estimated coding cost.
      If a Cb or Cr block is coded but no Y' block from a macro block then
       the mode MUST be OC_MODE_INTER_NOMV.
      This is the default state to which the mode data structure is
       initialised in encoder and decoder at the start of each frame.*/
    /*Block coding cost is estimated from correlated SATD metrics.*/
    /*At this point, all blocks that are in frame are still marked coded.*/
    if(!_recode){
      embs[mbi].unref_mv[OC_FRAME_GOLD]=embs[mbi].analysis_mv[0][OC_FRAME_GOLD];
      embs[mbi].unref_mv[OC_FRAME_PREV]=embs[mbi].analysis_mv[0][OC_FRAME_PREV];
      embs[mbi].refined=0;
    }

    /*Estimate the cost in a delta frame for various modes.*/
    oc_skip_cost_ssd(_enc,mbi,d->rd_scale,d->skip_ssd);
    oc_cost_inter_nomv_dct(_enc,mbi,OC_MODE_INTER_NOMV,&d->i0mv_dct);
    oc_cost_inter_nomv_dct(_enc,mbi,OC_MODE_GOLDEN_NOMV,&d->g0mv_dct);
    d->computed = 0;

    if(sp_level<OC_SP_LEVEL_NOMC){
      if(!(embs[mbi].refined&(1<<OC_MODE_INTER_MV))){
        oc_mcenc_refine1mv(_enc,mbi,OC_FRAME_PREV);
        embs[mbi].refined|=(1<<OC_MODE_INTER_MV);
      }

      oc_cost_inter1mv_dct(_enc,mbi,OC_MODE_INTER_MV,
                            embs[mbi].unref_mv[OC_FRAME_PREV],&d->i1mv_dct);
      oc_cost_inter1mv_dct(_enc,mbi,OC_MODE_INTER_MV,
                            embs[mbi].analysis_mv[0][OC_FRAME_PREV],&d->i1mvr_dct);
      oc_cost_inter1mv_dct(_enc,mbi,OC_MODE_GOLDEN_MV,
                            embs[mbi].unref_mv[OC_FRAME_GOLD],&d->g1mv_dct);
      d->computed = WRK_G1MV_EST | WRK_I1MV_EST | WRK_I1MV_REFINED;

      if(_enc->threads->stall)
        return;

      if(!(embs[mbi].refined&(1<<OC_MODE_GOLDEN_MV))){
        oc_mcenc_refine1mv(_enc,mbi,OC_FRAME_GOLD);
        embs[mbi].refined|=(1<<OC_MODE_GOLDEN_MV);
      }

      oc_cost_inter1mv_dct(_enc,mbi,OC_MODE_GOLDEN_MV,
                            embs[mbi].analysis_mv[0][OC_FRAME_GOLD],&d->g1mvr_dct);

      d->computed |= WRK_G1MV_REFINED;
    }
  }
}

static void intra_analysis(oc_enc_ctx *_enc,int _recode,unsigned mbi,struct oc_enc_wmb_data *d) {
  const oc_sb_map     *sb_maps;
  oc_fragment         *frags;
  unsigned             activity_avg;
  unsigned             luma_avg;
  const ogg_uint16_t  *chroma_rd_scale;
  unsigned             sbi;
  unsigned             quadi;
  sb_maps=(const oc_sb_map *)_enc->state.sb_maps;
  frags=_enc->state.frags;
  activity_avg=_enc->activity_avg;
  luma_avg=OC_CLAMPI(90<<8,_enc->luma_avg,160<<8);
  chroma_rd_scale=_enc->chroma_rd_scale[OC_INTER_FRAME][_enc->state.qis[0]];
  sbi=mbi>>2;
  quadi=mbi&3;

  if(_enc->state.sb_flags[sbi].quad_valid&1<<quadi){
    unsigned  intra_satd[12];
    unsigned  activity[4];
    unsigned  luma;
    int       bi;
    /*Activity masking.*/
    if(_enc->sp_level<OC_SP_LEVEL_FAST_ANALYSIS){
      luma=oc_mb_activity(_enc,mbi,activity);
    }
    else{
      luma=oc_mb_intra_dct(_enc,mbi,intra_satd,&d->intra_dct);
      oc_mb_activity_fast(_enc,mbi,activity,intra_satd);
      for(bi=0;bi<4;bi++)frags[sb_maps[sbi][quadi][bi]].qii=0;
    }
    d->luma = luma;
    d->activity = oc_mb_masking(d->rd_scale,d->rd_iscale,chroma_rd_scale,activity,activity_avg,luma,luma_avg);
    /*Motion estimation:
      We do a basic 1MV search for all macroblocks, coded or not,
       keyframe or not, unless we aren't using motion estimation at all.*/
    if(!_recode&&_enc->state.curframe_num>0&&
     _enc->sp_level<OC_SP_LEVEL_NOMC&&_enc->keyframe_frequency_force>1){
      oc_mcenc_search(_enc,mbi);
    }
  }
}

void oc_enc_worker_start(oc_enc_ctx *enc, unsigned sbi_start, unsigned sbi_end, int recode, int intra) {
  struct oc_enc_worker_ctrl *w = enc->threads;
#ifdef HAVE_PTHREAD
  pthread_barrier_wait(&w->end_sync);
#endif
  unsigned current_capacity = w->mbi_end - w->mbi_start;
  w->mbi_start = sbi_start << 2;
  w->mbi_end = sbi_end << 2;
  unsigned blocks = w->mbi_end - w->mbi_start;
  w->recode = !!recode;
  w->intra = !!intra;
  w->stall = 1;
  atomic_store(&w->current_mbi, w->mbi_start);
  if (!w->workers) blocks = 1;

  if (current_capacity < blocks) {
    free(w->finished);
    free(w->started);
    free(w->mbs);
    w->started = calloc(blocks, sizeof(*w->started));
    w->finished = calloc(blocks, sizeof(*w->finished));
    w->mbs = calloc(blocks, sizeof(*w->mbs));
  } else {
    memset(w->finished, 0, blocks * sizeof(*w->finished));
  }

  for(int i = 0; i < blocks; i++)
    atomic_flag_clear(&w->started[i]);

#ifdef HAVE_PTHREAD
  pthread_barrier_wait(&w->start_sync);
#endif
}

static void do_work(oc_enc_ctx *enc, unsigned wmbi, unsigned mbi) {
  struct oc_enc_worker_ctrl *w = enc->threads;

  if(w->intra)
    intra_analysis(enc, w->recode, mbi, &w->mbs[wmbi]);
  else
    inter_analysis(enc, w->recode, mbi, &w->mbs[wmbi]);

  w->finished[wmbi] = 1;
}

#ifdef HAVE_PTHREAD
static void *worker_thread(void *_enc) {
  oc_enc_ctx *enc = _enc;
  struct oc_enc_worker_ctrl *w = enc->threads;
  do {
    pthread_barrier_wait(&w->start_sync);

    for (int n = 0; n < 8; n++) {
      unsigned mbi = atomic_fetch_add(&w->current_mbi, 1);
      unsigned wmbi = mbi - w->mbi_start;
      if (mbi >= w->mbi_end) break;
      if (!atomic_flag_test_and_set(&w->started[wmbi]))
        do_work(enc, wmbi, mbi);
      pthread_cond_signal(&w->worker_comp_signal);
    }

    while (1) {
      unsigned mbi = atomic_fetch_add(&w->current_mbi, 4);
      unsigned wmbi = mbi - w->mbi_start;

      for (int n = 0; n < 4; n++) {
        if (mbi >= w->mbi_end) break;
        if (!atomic_flag_test_and_set(&w->started[wmbi]))
          do_work(enc, wmbi, mbi);
        mbi++; wmbi++;
      }

      if (mbi >= w->mbi_end) break;
      pthread_cond_signal(&w->worker_comp_signal);
    }

    pthread_barrier_wait(&w->end_sync);
  } while(!w->stop);

  return NULL;
}
#endif

void oc_enc_worker_free(oc_enc_ctx* enc) {
  struct oc_enc_worker_ctrl *w = enc->threads;
  if (!w) return;
  if (w->stop) return;
  enc->threads = NULL;

  w->stop = 1;
#ifdef HAVE_PTHREAD
  pthread_barrier_wait(&w->end_sync);
  for(int i = 0; i < w->threads - 1; i++)
    pthread_join(w->workers[i], NULL);

  pthread_cond_destroy(&w->worker_comp_signal);
  pthread_mutex_destroy(&w->worker_comp_lock);
  pthread_rwlock_destroy(&w->frame_state_lock);
  pthread_barrier_destroy(&w->end_sync);
  pthread_barrier_destroy(&w->start_sync);
  free(w->workers);
#endif
  free(w->finished);
  free(w->started);
  free(w->mbs);

  free(w);
}

int oc_enc_worker_set_threads(oc_enc_ctx* enc, int n) {
  int current_threads = 0;
  if (enc->threads) current_threads = enc->threads->threads;
  if (n < 1) return current_threads;
  if (n == current_threads) return n;
  if (n > 10) n = 10;

  oc_enc_worker_free(enc);
  struct oc_enc_worker_ctrl *w = enc->threads = malloc(sizeof(struct oc_enc_worker_ctrl));
#ifdef HAVE_PTHREAD
  pthread_barrier_init(&w->start_sync, NULL, n);
  pthread_barrier_init(&w->end_sync, NULL, n);
  pthread_rwlock_init(&w->frame_state_lock, NULL);
  pthread_mutex_init(&w->worker_comp_lock, NULL);
  pthread_cond_init(&w->worker_comp_signal, NULL);
  w->threads = n;
#else
  w->threads = n = 1;
#endif
  w->stop = 0;
  w->mbi_start = 0;
  w->mbi_end = 0;
  w->mbs = NULL;
  w->started = NULL;
  w->finished = NULL;

  w->workers = n > 1 ? calloc(n-1, sizeof(pthread_t)) : NULL;
#ifdef HAVE_PTHREAD
  pthread_attr_t worker_attr;
  struct sched_param worker_sch;
  int worker_policy;
  pthread_attr_init(&worker_attr);
#ifdef SCHED_IDLE
  pthread_attr_setschedpolicy(&worker_attr, SCHED_IDLE);
#endif
  pthread_attr_getschedpolicy(&worker_attr, &worker_policy);
  pthread_attr_getschedparam(&worker_attr, &worker_sch);
  worker_sch.sched_priority = sched_get_priority_min(worker_policy);
  pthread_attr_setschedparam(&worker_attr, &worker_sch);

  for(int i = 0; i < n-1; i++) {
    pthread_create(&w->workers[i], &worker_attr, worker_thread, enc);
    char th_name[17];
    snprintf(th_name, 17, "theora E%d", i+1);
    pthread_setname_np(w->workers[i], th_name);
  }

  pthread_attr_destroy(&worker_attr);
  pthread_barrier_wait(&w->start_sync);
#endif
  return n;
}

void oc_enc_worker_wait(oc_enc_ctx *enc, unsigned mbi) {
  struct oc_enc_worker_ctrl *w = enc->threads;
  unsigned wmbi;
  if (!w->workers) {
    atomic_flag_clear(w->started);
    *w->finished = 0;
    wmbi = 0;
  } else {
    wmbi = mbi - w->mbi_start;
  }

  if (!atomic_flag_test_and_set(&w->started[wmbi])) {
    w->stall = 1;
    do_work(enc, wmbi, mbi);
  } else {
#ifdef HAVE_PTHREAD
    w->stall = 0;
    while (1) {
      if(w->finished[wmbi]) break;
      w->stall = 1;
      struct timespec poll_time;
      clock_gettime(CLOCK_REALTIME, &poll_time);
      poll_time.tv_nsec += 4000000;
      if (poll_time.tv_nsec > 999999999) {
        poll_time.tv_nsec -= 1000000000;
        poll_time.tv_sec += 1;
      }

      pthread_mutex_lock(&w->worker_comp_lock);
      pthread_cond_timedwait(&w->worker_comp_signal, &w->worker_comp_lock, &poll_time);
      pthread_mutex_unlock(&w->worker_comp_lock);
    }
#else
    oc_assume(0);
#endif
  }
}

void oc_enc_worker_get_rd_acc(oc_enc_ctx *enc, unsigned mbi,
                              unsigned **rd_scale, unsigned **rd_iscale, oc_dct_cost_table **intra_dct,
                              int64_t *luma_accum, int64_t *activity_accum) {
  struct oc_enc_worker_ctrl *w = enc->threads;
  unsigned wmbi = w->workers ? mbi - w->mbi_start : 0;

  oc_assume(w->finished[wmbi]);
  *rd_scale = w->mbs[wmbi].rd_scale;
  *rd_iscale = w->mbs[wmbi].rd_iscale;
  *intra_dct = &w->mbs[wmbi].intra_dct;
  *luma_accum += w->mbs[wmbi].luma;
  *activity_accum += w->mbs[wmbi].activity;
}

void oc_enc_worker_get_nomv_dct(oc_enc_ctx *enc, unsigned mbi, unsigned **skip_ssd,
                                oc_dct_cost_table **inter_dct, oc_dct_cost_table **golden_dct) {
  struct oc_enc_worker_ctrl *w = enc->threads;
  unsigned wmbi = w->workers ? mbi - w->mbi_start : 0;

  oc_assume(w->finished[wmbi]);
  oc_assume(!w->intra);
  *skip_ssd = w->mbs[wmbi].skip_ssd;
  *inter_dct = &w->mbs[wmbi].i0mv_dct;
  *golden_dct = &w->mbs[wmbi].g0mv_dct;
}

void oc_enc_worker_get_mv_dct(oc_enc_ctx *enc, unsigned mbi, int which, int refined, oc_dct_cost_table **dct) {
  struct oc_enc_worker_ctrl *w = enc->threads;
  unsigned wmbi = w->workers ? mbi - w->mbi_start : 0;

  oc_assume(w->finished[wmbi]);
  oc_assume(!w->intra);

  struct oc_enc_wmb_data *d = &w->mbs[wmbi];
  enum oc_enc_wrk_computed computed;
  int ref_frame;
  switch(which) {
    case OC_MODE_INTER_MV:
      computed = refined ? WRK_I1MV_REFINED : WRK_I1MV_EST;
      *dct = refined ? &d->i1mvr_dct : &d->i1mv_dct;
      ref_frame = OC_FRAME_PREV;
      break;
    case OC_MODE_GOLDEN_MV:
      computed = refined ? WRK_G1MV_REFINED : WRK_G1MV_EST;
      *dct = refined ? &d->g1mvr_dct : &d->g1mv_dct;
      ref_frame = OC_FRAME_GOLD;
      break;
    default: abort();
  }

  if(!(d->computed & computed)) {
    oc_mb_enc_info *embs = enc->mb_info;
    oc_mv          *mvs = refined ? embs[mbi].analysis_mv[0] : embs[mbi].unref_mv;

    if(refined) {
      if(!(embs[mbi].refined&(1<<which))){
        oc_mcenc_refine1mv(enc,mbi,ref_frame);
        embs[mbi].refined|=(1<<which);
      }
    }

    oc_cost_inter1mv_dct(enc,mbi,which,mvs[ref_frame],*dct);
    d->computed |= computed;
  }
}
