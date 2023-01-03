#include <arm_neon.h>
#include <stdint.h>
#include <stddef.h>

void oc_frag_copy_neon(unsigned char *_dst,const unsigned char *_src,int _ystride) {
  for (int i = 0; i < 8; i++)
    vst1_u8(&_dst[i * _ystride], vld1_u8(&_src[i * _ystride]));
}

void oc_frag_copy_list_neon(unsigned char *_dst_frame, const unsigned char *_src_frame,int _ystride,
                            const ptrdiff_t *_fragis,ptrdiff_t _nfragis,const ptrdiff_t *_frag_buf_offs) {
  ptrdiff_t fragii;
  for (fragii = 0; fragii < _nfragis; fragii++) {
    ptrdiff_t frag_buf_off;
    frag_buf_off = _frag_buf_offs[_fragis[fragii]];
    oc_frag_copy_neon(_dst_frame + frag_buf_off,
                      _src_frame + frag_buf_off, _ystride);
  }
}

void oc_frag_recon_intra_neon(unsigned char *_dst,
                              int _ystride, const int16_t _residue[64]) {
  for (int i = 0; i < 8; i++)
    vst1_u8(&_dst[i * _ystride],
            vqmovun_s16(vaddq_s16(vld1q_s16(&_residue[i * 8]),
                                  vdupq_n_s16(128))));
}

void oc_frag_recon_inter_neon(unsigned char *_dst, const unsigned char *_src,
                              int _ystride,const int16_t _residue[64]) {
  for (int i = 0; i < 8; i++)
    vst1_u8(&_dst[i * _ystride],
            vqmovun_s16(vaddq_s16(vld1q_s16(&_residue[i * 8]),
                                  vreinterpretq_s16_u16(vmovl_u8(vld1_u8(&_src[i * _ystride]))))));
}

void oc_frag_recon_inter2_neon(unsigned char *_dst,const unsigned char *_src1, const unsigned char *_src2,
                               int _ystride,const int16_t _residue[64]) {
  for (int i = 0; i < 8; i++)
    vst1_u8(&_dst[i * _ystride],
            vqmovun_s16(vaddq_s16(vld1q_s16(&_residue[i * 8]),
                                  vreinterpretq_s16_u16(vmovl_u8(vhadd_u8(vld1_u8(&_src1[i * _ystride]),
                                                                          vld1_u8(&_src2[i * _ystride])))))));
}
