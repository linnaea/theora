#ifndef THEORA_ENC_WORKER_H
#define THEORA_ENC_WORKER_H
#include "encint.h"

void oc_cost_inter_satd(oc_enc_ctx *_enc, unsigned _mbi, int _mb_mode, oc_mv _mv, unsigned frag_satd[12]);
void oc_enc_worker_start(oc_enc_ctx *enc, unsigned sbi_start, unsigned sbi_end, int recode, int intra);

void oc_enc_worker_wait(oc_enc_ctx *_enc, unsigned mbi);

void oc_enc_worker_get_rd_acc(oc_enc_ctx *_enc, unsigned mbi,
                              unsigned **rd_scale, unsigned **rd_iscale,
                              int64_t *luma_accum, int64_t *activity_accum);

void oc_enc_worker_get_nomv_satd(oc_enc_ctx *_enc, unsigned mbi, unsigned **skip_ssd, unsigned **intra_satd,
                                 unsigned **inter_satd, unsigned **golden_satd);

void oc_enc_worker_get_mv_satd(oc_enc_ctx *_enc, unsigned mbi, int which, int refined, unsigned **satd);

#endif //THEORA_ENC_WORKER_H
