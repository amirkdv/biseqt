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
 * Given a segment `seg`, extends it in the given direction by one window.
 *
 * @param res Extended segment to be populated here.
 * @param seg The original segment.
 * @param S The integer array for the "from" sequence containing indices of
 *   its letters in alphabet.
 * @param T The integer array for the "to" sequence containing indices of
 *   its letters in alphabet.
 * @param scores Alignment parameters to be used over the window.
 * @param window The length of the extension alignment window.
 * @param forward Either of 0 or 1 indicating the direction of extension.
 *
 * @return -1 if an error occurs and 0 otherwise.
 */
int extend_1d_once(segment* res, segment* seg,
  int* S, int* T, alnscores* scores,
  int window, int forward) {

  alnprob prob;
  alnframe frame;
  dptable table;
  intpair opt;
  std_alntype type;
  transcript* tx;
  int tx_opseq_len, seg_opseq_len;
  char* opseq;
  int failure = 0;

  int S_len = tx_seq_len(seg->tx, 'S'),
      T_len = tx_seq_len(seg->tx, 'T');
  intpair S_range, T_range;
  if (forward) {
    S_range = (intpair) {seg->tx->S_idx + S_len, seg->tx->S_idx + S_len + window};
    T_range = (intpair) {seg->tx->T_idx + T_len, seg->tx->T_idx + T_len + window};
  } else {
    S_range = (intpair) {seg->tx->S_idx - window, seg->tx->S_idx};
    T_range = (intpair) {seg->tx->T_idx - window, seg->tx->T_idx};
  }

  type = forward ? START_ANCHORED_OVERLAP : END_ANCHORED_OVERLAP;
  frame = (alnframe) {.S=S, .T=T, .S_range=S_range, .T_range=T_range};
  std_alnparams params = (std_alnparams) {.type=type};
  prob = (alnprob) {.frame=&frame, .scores=scores, .mode=STD_MODE, .std_params=&params};
  table = (dptable) {.prob=&prob, .cells=NULL, .row_lens=NULL, .num_rows=-1};
  if (dptable_init(&table) == -1) {
    return -1;
  }
  opt = dptable_solve(&table);
  if (opt.i == -1 || opt.j == -1) {
    failure = 1;
  }

  tx = dptable_traceback(&table, opt);
  if (tx == NULL) {
    failure = 1;
  }
  dptable_free(&table);
  if (failure) {
    return -1;
  }

  // There is an alignment, return the corresponding extended segment:
  seg_opseq_len = strlen(seg->tx->opseq);
  tx_opseq_len = strlen(tx->opseq);
  tx->score += seg->tx->score;
  opseq = malloc(seg_opseq_len + tx_opseq_len + 1);
  if (forward) {
    tx->S_idx = seg->tx->S_idx;
    tx->T_idx = seg->tx->T_idx;
    strncpy(opseq, seg->tx->opseq, seg_opseq_len);
    strncpy(opseq + seg_opseq_len, tx->opseq, tx_opseq_len);
  } else {
    strncpy(opseq, tx->opseq, tx_opseq_len);
    strncpy(opseq + tx_opseq_len, seg->tx->opseq, seg_opseq_len);
  }
  // strncpy does not null-terminate:
  opseq[seg_opseq_len + tx_opseq_len] = '\0';
  tx->opseq = opseq;

  res->S_id = seg->S_id;
  res->T_id = seg->T_id;
  res->tx = tx;
  return 0;
}

/**
 * Given a segment fully extends it in one direction. A fully extended segment
 * (in one direction) is one that hits the boundary (beginning or end) of either
 * of the sequences.
 *
 * @param res Extended segment to be populated here.
 * @param seg The original segment.
 * @param S The integer array for the "from" sequence containing indices of
 *   its letters in alphabet.
 * @param T The integer array for the "to" sequence containing indices of
 *   its letters in alphabet.
 * @param S_len The length of the "from" sequence.
 * @param T_len The length of the "to" sequence.
 * @param scores Alignment parameters to be used over the window.
 * @param window The length of the extension alignment window.
 * @param max_new_mins Maximum number of new minima observed in the score
 *    random walk until the segement is dropped (i.e -1 is returned).
 * @param forward Either of 0 or 1 indicating the direction of extension.
 * @param debug Whether to dump score random walks for all tried extensions. If
 *    truthy, each segment's score random walk is written as a line in a file
 *    "scores.txt" (The file is never truncated; all data is appended to it).
 *
 * @return 0 if the seed successfully extends to the boundary of either of
 *   the sequences and -1 otherwise.
 */
int extend_1d(segment* res, segment* seg, int* S, int* T, int S_len, int T_len,
  alnscores* scores, int window, int max_new_mins, int forward, int debug) {

  FILE* f;
  if (debug) {
    f = fopen("scores.txt", "a");
    if (f == NULL) {
      printf("Failed to open file\n");
      exit(1);
    }
    fprintf(f, "(%d:%d,%d:%d) [", seg->S_id, seg->tx->S_idx, seg->T_id, seg->tx->T_idx);
  }

  segment cur_seg = *seg;
  double cur_min = seg->tx->score;
  int num_mins = 0, retcode, actual_window;
  int S_end, T_end, S_wiggle, T_wiggle, min_wiggle;
  while (1) {
    if (forward) {
      S_end = cur_seg.tx->S_idx + tx_seq_len(cur_seg.tx, 'S');
      T_end = cur_seg.tx->T_idx + tx_seq_len(cur_seg.tx, 'T');
      S_wiggle = S_len - S_end;
      T_wiggle = T_len - T_end;
      min_wiggle = S_wiggle < T_wiggle ? S_wiggle : T_wiggle;
    } else {
      min_wiggle = cur_seg.tx->S_idx < cur_seg.tx->T_idx ?
        cur_seg.tx->S_idx : cur_seg.tx->T_idx;
    }

    actual_window = window < min_wiggle ? window : min_wiggle;
    if (actual_window == 0) {
      // hit the end
      *res = cur_seg;
      if (debug) {
        fprintf(f, "] +\n");
        fclose(f);
      }
      return 0;
    }

    retcode = extend_1d_once(&cur_seg, &cur_seg, S, T, scores, actual_window, forward);
    if (retcode == -1) {
      // No nonempty alignment found:
      if (debug) {
        fprintf(f, "] -\n");
      }
      return -1;
    }
    if (debug) {
      fprintf(f, "%.2f,", cur_seg.tx->score);
    }

    if (cur_seg.tx->score < cur_min) {
      num_mins +=1;
      if (num_mins >= max_new_mins) {
        if (debug) {
          fprintf(f, "] -\n");
          fclose(f);
        }
        return -1;
      }
    }
  }
  if (debug) {
    fprintf(f, "] -\n");
    fclose(f);
  }
  return -1;
}

/**
 * Given an array of segments tries to extend all in both directions and returns
 * a fully extended segment as soon as it finds one.
 *
 * @param segs The original segments.
 * @param num_segs The number of provided segments.
 * @param S The integer array for the "from" sequence containing indices of
 *   its letters in alphabet.
 * @param T The integer array for the "to" sequence containing indices of
 *   its letters in alphabet.
 * @param S_len The length of the "from" sequence.
 * @param T_len The length of the "to" sequence.
 * @param scores Alignment parameters to be used over the window.
 * @param window The length of the extension alignment window.
 * @param max_new_mins Maximum number of new minima observed in the score
 *    random walk until the segement is dropped (i.e -1 is returned).
 * @param min_overlap_score The minimum overall score required for a fully
 *    extended segment to be reported as an overlap alignment.
 * @param debug Whether to dump score random walks for all tried extensions. If
 *    truthy, each segment's score random walk is written as a line in a file
 *    "scores.txt" (The file is never truncated; all data is appended to it).
 *
 * @return 0 if the seed successfully extends to the boundary of either of
 *   the sequences and -1 otherwise.
 */
segment* extend(segment** segs, int num_segs, int* S, int* T, int S_len, int T_len,
  alnscores* scores, int window, int max_new_mins, double min_overlap_score, int debug) {

  segment fwd, bwd, *res;
  transcript* tx;
  char* opseq;
  int fwd_tx_len, bwd_tx_len, seg_tx_len, retcode, overlap_score;
  for (int i = 0; i < num_segs; i ++) {
    retcode = extend_1d(&fwd, segs[i],
        S, T, S_len, T_len, scores,
        window, max_new_mins, 1, debug);
    if (retcode == -1) {
      continue;
    }
    retcode = extend_1d(&bwd, segs[i],
      S, T, S_len, T_len, scores,
      window, max_new_mins, 0, debug);
    if (retcode == -1) {
      continue;
    }
    overlap_score = bwd.tx->score + fwd.tx->score - segs[i]->tx->score;
    if (overlap_score <= min_overlap_score) {
      continue;
    }
    // Found a fully extending segment; return it:
    fwd_tx_len = strlen(fwd.tx->opseq);
    bwd_tx_len = strlen(bwd.tx->opseq);
    seg_tx_len = strlen(segs[i]->tx->opseq);
    opseq = malloc(fwd_tx_len + bwd_tx_len - seg_tx_len + 1);
    strncpy(opseq, bwd.tx->opseq, bwd_tx_len - seg_tx_len);
    strncpy(opseq + bwd_tx_len - seg_tx_len, fwd.tx->opseq, fwd_tx_len);
    // strncpy does not null terminate:
    opseq[fwd_tx_len + bwd_tx_len - seg_tx_len] = '\0';
    // build the overall transcript
    tx = malloc(sizeof(transcript));
    tx->S_idx = bwd.tx->S_idx;
    tx->T_idx = bwd.tx->T_idx;
    tx->score = overlap_score;
    tx->opseq = opseq;

    // build the fully extended segment
    res = malloc(sizeof(segment));
    res->S_id = segs[i]->S_id;
    res->T_id = segs[i]->T_id;
    res->tx = tx;
    return res;
  }
  return NULL;
}
