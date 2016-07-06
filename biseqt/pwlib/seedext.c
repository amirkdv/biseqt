#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdio.h>
#include <string.h>

#include "pwlib.h"

/**
 * Given an alignment ::transcript returns the length of its opseq on the
 * "from" or "to" sequence.
 *
 * @param tx The transcript of interest.
 * @param on A single character of either 'S' or 'T'.
 * @return A non-negative integer corresponding to the length of the indicated
 *   sequence that is covered by the transcript and -1 if an error occurs.
 */
int tx_seq_len(transcript* tx, char on) {
  if (on != 'S' && on != 'T') {
    return -1;
  }
  int sum = 0;
  for (int i = 0; i < strlen(tx->opseq); i++) {
    if ((on == 'S' && tx->opseq[i] != 'I') ||
        (on == 'T' && tx->opseq[i] != 'D')) {
      sum ++;
    }
  }
  return sum;
}

/**
 * Given a segment transcript, extends it in the given direction by one window.
 *
 * @param tx The transcript to be extended.
 * @param frame The frame defining the alignment of the original sequences.
 * @param params Seed extension parameters.
 * @param forward Either of 0 or 1 indicating the direction of extension.
 *
 * @return -1 if an error occurs and 0 otherwise.
 */
int extend_1d_once(transcript* tx, alnframe* frame, seedext_params* params, int forward) {

  transcript* res;
  alnprob prob;
  alnframe ext_frame;
  dptable table;
  intpair opt;
  std_alntype type;
  int res_opseq_len, seg_opseq_len;
  char* opseq;
  int failure = 0;

  int S_len = tx_seq_len(tx, 'S'),
      T_len = tx_seq_len(tx, 'T');
  intpair ext_S_range, ext_T_range;
  if (forward) {
    ext_S_range = (intpair) {tx->S_idx + S_len, tx->S_idx + S_len + params->window};
    ext_T_range = (intpair) {tx->T_idx + T_len, tx->T_idx + T_len + params->window};
  } else {
    ext_S_range = (intpair) {tx->S_idx - params->window, tx->S_idx};
    ext_T_range = (intpair) {tx->T_idx - params->window, tx->T_idx};
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
  seg_opseq_len = strlen(tx->opseq);
  res_opseq_len = strlen(res->opseq);
  tx->score += res->score;
  opseq = malloc(seg_opseq_len + res_opseq_len + 1);
  if (forward) {
    strncpy(opseq, tx->opseq, seg_opseq_len);
    strncpy(opseq + seg_opseq_len, tx->opseq, res_opseq_len);
  } else {
    tx->S_idx = res->S_idx;
    tx->T_idx = res->T_idx;
    strncpy(opseq, res->opseq, res_opseq_len);
    strncpy(opseq + res_opseq_len, tx->opseq, seg_opseq_len);
  }
  // strncpy does not null-terminate:
  opseq[seg_opseq_len + res_opseq_len] = '\0';
  tx->opseq = opseq;
  return 0;
}

/**
 * Given a segment fully extends it in one direction. A fully extended segment
 * (in one direction) is one that hits the boundary (beginning or end) of either
 * of the sequences.
 *
 * @param res Extended segment to be populated here.
 * @param tx The original segment.
 * @param frame The frame defining the alignment of the original sequences.
 * @param params Seed extension parameters.
 * @param forward Either of 0 or 1 indicating the direction of extension.
 *
 * @return 0 if the seed successfully extends to the boundary of either of
 *   the sequences and -1 otherwise.
 */
int extend_1d(transcript* res, transcript* tx, alnframe* frame, seedext_params* params, int forward) {

  transcript cur_tx = *tx;
  double cur_min = tx->score;
  int num_mins = 0, actual_window;
  int S_end, T_end, S_wiggle, T_wiggle, min_wiggle;
  seedext_params* curparams = malloc(sizeof(seedext_params));
  while (1) {
    if (forward) {
      S_end = cur_tx.S_idx + tx_seq_len(&cur_tx, 'S');
      T_end = cur_tx.T_idx + tx_seq_len(&cur_tx, 'T');
      S_wiggle = frame->S_range.j - frame->S_range.i - S_end;
      T_wiggle = frame->T_range.j - frame->T_range.i - T_end;
      min_wiggle = S_wiggle < T_wiggle ? S_wiggle : T_wiggle;
    } else {
      min_wiggle = cur_tx.S_idx < cur_tx.T_idx ? cur_tx.S_idx : cur_tx.T_idx;
    }

    actual_window = params->window < min_wiggle ? params->window : min_wiggle;
    if (actual_window < 0) {
      // hit the end
      *res = cur_tx;
      return 0;
    }

    *curparams = *params;
    curparams->window = actual_window;
    if (extend_1d_once(&cur_tx, frame, curparams, forward) == -1) {
      // No nonempty alignment found:
      return -1;
    }

    if (cur_tx.score < cur_min) {
      num_mins +=1;
      if (num_mins >= params->max_new_mins) {
        return -1;
      }
    }
  }
  return -1;
}

/**
 * Given an array of segments tries to extend all in both directions and returns
 * a fully extended segment as soon as it finds one.
 *
 * @param txs The original segments.
 * @param num_txs The number of provided segments.
 * @param frame The frame defining the alignment of the original sequences.
 * @param params Seed extension parameters.
 *
 * @return 0 if the seed successfully extends to the boundary of either of
 *   the sequences and -1 otherwise.
 */
transcript* extend(transcript** txs, int num_txs, alnframe* frame, seedext_params* params) {

  transcript fwd, bwd;
  transcript* tx;
  char* opseq;
  int fwd_tx_len, bwd_tx_len, seg_tx_len;
  double overlap_score;
  for (int i = 0; i < num_txs; i ++) {
    if (extend_1d(&fwd, txs[i], frame, params, 1) == -1) {
      continue;
    }
    if (extend_1d(&bwd, txs[i], frame, params, 0) == -1) {
      continue;
    }
    overlap_score = bwd.score + fwd.score - txs[i]->score;
    if (overlap_score < params->min_score) {
      continue;
    }
    // Found a fully extending segment; return it:
    fwd_tx_len = strlen(fwd.opseq);
    bwd_tx_len = strlen(bwd.opseq);
    seg_tx_len = strlen(txs[i]->opseq);
    opseq = malloc(fwd_tx_len + bwd_tx_len - seg_tx_len + 1);
    strncpy(opseq, bwd.opseq, bwd_tx_len - seg_tx_len);
    strncpy(opseq + bwd_tx_len - seg_tx_len, fwd.opseq, fwd_tx_len);
    // strncpy does not null terminate:
    opseq[fwd_tx_len + bwd_tx_len - seg_tx_len] = '\0';
    // build the overall transcript
    tx = malloc(sizeof(transcript));
    tx->S_idx = bwd.S_idx;
    tx->T_idx = bwd.T_idx;
    tx->score = overlap_score;
    tx->opseq = opseq;
    return tx;
  }
  return NULL;
}
