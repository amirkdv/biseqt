#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdio.h>
#include <string.h>

#include "pwlib.h"

/**
 * Given an alignment ::alignment returns the length of its opseq on the
 * "from" or "to" sequence.
 *
 * @param aln The alignment of interest.
 * @param on A single character of either 'S' or 'T'.
 * @return A non-negative integer corresponding to the length of the indicated
 *   sequence that is covered by the alignment and -1 if an error occurs.
 */
int aln_seq_len(alignment* aln, char on) {
  if (on != 'S' && on != 'T') {
    return -1;
  }
  int sum = 0;
  for (int i = 0; i < strlen(aln->opseq); i++) {
    if ((on == 'S' && aln->opseq[i] != 'I') ||
        (on == 'T' && aln->opseq[i] != 'D')) {
      sum ++;
    }
  }
  return sum;
}

int extend_1d_once(alignment* aln, alnframe* frame, seedext_params* params, int forward) {
  alignment* res;
  alnprob prob;
  alnframe ext_frame;
  dptable table;
  intpair opt;
  std_alntype type;
  int res_opseq_len, seg_opseq_len;
  char* opseq;
  int failure = 0;

  int S_len = aln_seq_len(aln, 'S'),
      T_len = aln_seq_len(aln, 'T');
  intpair ext_S_range, ext_T_range;
  if (forward) {
    ext_S_range = (intpair) {aln->S_idx + S_len, aln->S_idx + S_len + params->window};
    ext_T_range = (intpair) {aln->T_idx + T_len, aln->T_idx + T_len + params->window};
  } else {
    ext_S_range = (intpair) {aln->S_idx - params->window, aln->S_idx};
    ext_T_range = (intpair) {aln->T_idx - params->window, aln->T_idx};
  }

  type = forward ? START_ANCHORED_OVERLAP : END_ANCHORED_OVERLAP;
  ext_frame = (alnframe) {
    .S=frame->S, .T=frame->T, .S_range=ext_S_range, .T_range=ext_T_range
  };
  std_alnparams alnparams = (std_alnparams) {.type=type};
  prob = (alnprob) {
    .frame=&ext_frame, .scores=params->scores, .mode=STD_MODE, .std_params=&alnparams
  };
  table = (dptable) {.prob=&prob, .cells=NULL, .row_lens=NULL, .num_rows=-1};
  if (dptable_init(&table) == -1) {
    return -1;
  }
  opt = dptable_solve(&table);
  if (opt.i == -1 || opt.j == -1) {
    failure = 1;
  }

  res = dptable_traceback(&table, opt);
  if (res == NULL) {
    failure = 1;
  }
  dptable_free(&table);
  if (failure) {
    return -1;
  }

  // There is an alignment:
  seg_opseq_len = strlen(aln->opseq);
  res_opseq_len = strlen(res->opseq);
  aln->score += res->score;
  opseq = malloc(seg_opseq_len + res_opseq_len + 1);
  if (forward) {
    strncpy(opseq, aln->opseq, seg_opseq_len);
    strncpy(opseq + seg_opseq_len, aln->opseq, res_opseq_len);
  } else {
    aln->S_idx = res->S_idx;
    aln->T_idx = res->T_idx;
    strncpy(opseq, res->opseq, res_opseq_len);
    strncpy(opseq + res_opseq_len, aln->opseq, seg_opseq_len);
  }
  // strncpy does not null-terminate:
  opseq[seg_opseq_len + res_opseq_len] = '\0';
  aln->opseq = opseq;
  return 0;
}

int extend_1d(alignment* res, alignment* aln, alnframe* frame, seedext_params* params, int forward) {

  alignment cur_aln = *aln;
  double cur_min = aln->score;
  int num_mins = 0, actual_window;
  int S_end, T_end, S_wiggle, T_wiggle, min_wiggle;
  seedext_params* curparams = malloc(sizeof(seedext_params));
  while (1) {
    if (forward) {
      S_end = cur_aln.S_idx + aln_seq_len(&cur_aln, 'S');
      T_end = cur_aln.T_idx + aln_seq_len(&cur_aln, 'T');
      S_wiggle = frame->S_range.j - frame->S_range.i - S_end;
      T_wiggle = frame->T_range.j - frame->T_range.i - T_end;
      min_wiggle = S_wiggle < T_wiggle ? S_wiggle : T_wiggle;
    } else {
      min_wiggle = cur_aln.S_idx < cur_aln.T_idx ? cur_aln.S_idx : cur_aln.T_idx;
    }

    actual_window = params->window < min_wiggle ? params->window : min_wiggle;
    if (actual_window < 0) {
      // hit the end
      *res = cur_aln;
      return 0;
    }

    *curparams = *params;
    curparams->window = actual_window;
    if (extend_1d_once(&cur_aln, frame, curparams, forward) == -1) {
      // No nonempty alignment found:
      return -1;
    }

    if (cur_aln.score < cur_min) {
      num_mins +=1;
      if (num_mins >= params->max_new_mins) {
        return -1;
      }
    }
  }
  return -1;
}

alignment* extend(alignment** alns, int num_alns, alnframe* frame, seedext_params* params) {

  alignment fwd, bwd;
  alignment* aln;
  char* opseq;
  int fwd_aln_len, bwd_aln_len, seg_aln_len;
  double overlap_score;
  for (int i = 0; i < num_alns; i ++) {
    if (extend_1d(&fwd, alns[i], frame, params, 1) == -1) {
      continue;
    }
    if (extend_1d(&bwd, alns[i], frame, params, 0) == -1) {
      continue;
    }
    overlap_score = bwd.score + fwd.score - alns[i]->score;
    if (overlap_score < params->min_score) {
      continue;
    }
    // Found a fully extending segment; return it:
    fwd_aln_len = strlen(fwd.opseq);
    bwd_aln_len = strlen(bwd.opseq);
    seg_aln_len = strlen(alns[i]->opseq);
    opseq = malloc(fwd_aln_len + bwd_aln_len - seg_aln_len + 1);
    strncpy(opseq, bwd.opseq, bwd_aln_len - seg_aln_len);
    strncpy(opseq + bwd_aln_len - seg_aln_len, fwd.opseq, fwd_aln_len);
    // strncpy does not null terminate:
    opseq[fwd_aln_len + bwd_aln_len - seg_aln_len] = '\0';
    // build the overall alignment
    aln = malloc(sizeof(alignment));
    aln->S_idx = bwd.S_idx;
    aln->T_idx = bwd.T_idx;
    aln->score = overlap_score;
    aln->opseq = opseq;
    return aln;
  }
  return NULL;
}
