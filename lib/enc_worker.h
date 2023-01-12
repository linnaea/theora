#ifndef THEORA_ENC_WORKER_H
#define THEORA_ENC_WORKER_H
#include "encint.h"

//[pli][bi][qii][cost, ssd]
typedef struct oc_dct_cost_table {
  struct {
    unsigned cost, ssd;
  } block[12][3];
} oc_dct_cost_table;

void oc_cost_dct_fill(oc_enc_ctx *_enc, int nqis, int qti, int pli, int bi, unsigned satd, int shift, oc_dct_cost_table *_dct);
void oc_cost_inter_dct(oc_enc_ctx *_enc, unsigned _mbi, int _mb_mode, oc_mv _mv, oc_dct_cost_table *_dct);
void oc_enc_worker_start(oc_enc_ctx *enc, unsigned sbi_start, unsigned sbi_end, int recode, int intra);

void oc_enc_worker_wait(oc_enc_ctx *_enc, unsigned mbi);

void oc_enc_worker_get_rd_acc(oc_enc_ctx *enc, unsigned mbi,
                              unsigned **rd_scale, unsigned **rd_iscale, oc_dct_cost_table **intra_dct,
                              int64_t *luma_accum, int64_t *activity_accum);

void oc_enc_worker_get_nomv_dct(oc_enc_ctx *_enc, unsigned mbi, unsigned **skip_ssd,
                                 oc_dct_cost_table **inter_dct, oc_dct_cost_table **golden_dct);

void oc_enc_worker_get_mv_dct(oc_enc_ctx *enc, unsigned mbi, int which, int refined, oc_dct_cost_table **dct);

#endif //THEORA_ENC_WORKER_H
