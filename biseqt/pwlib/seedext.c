#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdio.h>
#include <string.h>

#include "pwlib.h"

/**
 * Given an alignment returns the length of its transcript on the
 * "from" or "to" sequence.
 *
 * @param aln The alignment of interest.
 * @param on A single character of either 'O' (for origin) or 'M' (for mutant).
 * @return A non-negative integer corresponding to the length of the indicated
 *   sequence that is covered by the alignment and -1 if an error occurs.
 */
int aln_seq_len(alignment* aln, char on) {
  if (on != 'O' && on != 'M') {
    return -1;
  }
  int sum = 0;
  for (int i = 0; i < strlen(aln->transcript); i++) {
    if ((on == 'O' && aln->transcript[i] != 'I') ||
        (on == 'M' && aln->transcript[i] != 'D')) {
      sum ++;
    }
  }
  return sum;
}

int extend_1d_once(alignment* aln, alnframe* frame, seedext_params* params,
                   int forward) {
  alignment* res;
  alnprob prob;
  alnframe ext_frame;
  dptable table;
  intpair opt;
  std_alntype type;
  int res_transcript_len, seg_transcript_len;
  char* transcript;
  int failure = 0;

  int origin_len = aln_seq_len(aln, 'O'),
      mutant_len = aln_seq_len(aln, 'M');
  intpair ext_origin_range, ext_mutant_range;
  if (forward) {
    ext_origin_range = (intpair) {
      aln->origin_idx + origin_len,
      aln->origin_idx + origin_len + params->window
    };
    ext_mutant_range = (intpair) {
      aln->mutant_idx + mutant_len,
      aln->mutant_idx + mutant_len + params->window
    };
  } else {
    ext_origin_range = (intpair) {
      aln->origin_idx - params->window,
      aln->origin_idx
    };
    ext_mutant_range = (intpair) {
      aln->mutant_idx - params->window,
      aln->mutant_idx
    };
  }

  type = forward ? START_ANCHORED_OVERLAP : END_ANCHORED_OVERLAP;
  ext_frame = (alnframe) {
    .origin=frame->origin,
    .mutant=frame->mutant,
    .origin_range=ext_origin_range,
    .mutant_range=ext_mutant_range
  };
  std_alnparams alnparams = (std_alnparams) {.type=type};
  prob = (alnprob) {
    .frame=&ext_frame,
    .scores=params->scores,
    .mode=STD_MODE,
    .std_params=&alnparams
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
  seg_transcript_len = strlen(aln->transcript);
  res_transcript_len = strlen(res->transcript);
  aln->score += res->score;
  transcript = malloc(seg_transcript_len + res_transcript_len + 1);
  if (forward) {
    strncpy(transcript, aln->transcript, seg_transcript_len);
    strncpy(transcript + seg_transcript_len, aln->transcript, res_transcript_len);
  } else {
    aln->origin_idx = res->origin_idx;
    aln->mutant_idx = res->mutant_idx;
    strncpy(transcript, res->transcript, res_transcript_len);
    strncpy(transcript + res_transcript_len, aln->transcript, seg_transcript_len);
  }
  // strncpy does not null-terminate:
  transcript[seg_transcript_len + res_transcript_len] = '\0';
  aln->transcript = transcript;
  return 0;
}

int extend_1d(alignment* res, alignment* aln, alnframe* frame,
              seedext_params* params, int forward) {
  alignment cur_aln = *aln;
  double cur_min = aln->score;
  int num_mins = 0, actual_window;
  int origin_end, mutant_end, origin_width, mutant_width, min_width;
  intpair origin_range = frame->origin_range,
          mutant_range = frame->mutant_range;
  seedext_params* curparams = malloc(sizeof(seedext_params));

  while (1) {
    if (forward) {
      origin_end = cur_aln.origin_idx + aln_seq_len(&cur_aln, 'O');
      mutant_end = cur_aln.mutant_idx + aln_seq_len(&cur_aln, 'M');
      origin_width = origin_range.j - origin_range.i - origin_end;
      mutant_width = mutant_range.j - mutant_range.i - mutant_end;
      min_width = origin_width < mutant_width ? origin_width : mutant_width;
    } else {
      if (cur_aln.origin_idx < cur_aln.mutant_idx) {
        min_width = cur_aln.origin_idx;
      } else {
        min_width = cur_aln.mutant_idx;
      }
    }

    actual_window = params->window < min_width ? params->window : min_width;
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

alignment* extend(alignment** alns, int num_alns, alnframe* frame,
                  seedext_params* params) {
  alignment fwd, bwd;
  alignment* aln;
  char* transcript;
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
    fwd_aln_len = strlen(fwd.transcript);
    bwd_aln_len = strlen(bwd.transcript);
    seg_aln_len = strlen(alns[i]->transcript);
    transcript = malloc(fwd_aln_len + bwd_aln_len - seg_aln_len + 1);
    strncpy(transcript, bwd.transcript, bwd_aln_len - seg_aln_len);
    strncpy(transcript + bwd_aln_len - seg_aln_len, fwd.transcript, fwd_aln_len);
    // strncpy does not null terminate:
    transcript[fwd_aln_len + bwd_aln_len - seg_aln_len] = '\0';
    // build the overall alignment
    aln = malloc(sizeof(alignment));
    aln->origin_idx = bwd.origin_idx;
    aln->mutant_idx = bwd.mutant_idx;
    aln->score = overlap_score;
    aln->transcript = transcript;
    return aln;
  }
  return NULL;
}
