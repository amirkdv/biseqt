#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdio.h>
#include <string.h>

#include "libalign.h"

/**
 * Given an alignment problem definition, creates and initializes the dynamic
 * programming table.
 *
 * @param def
 *    Alignment problem definition. The only members that are relevant at this
 *    point is the corresponding frames of the two sequences which determine
 *    the size of the table).
 *
 * @return pointer to the `malloc` ed dynamic programming table.
 */
align_dp_cell** init_dp_table(align_problem* def) {
  int n = def->S_max_idx - def->S_min_idx;
  int m = def->T_max_idx - def->T_min_idx;
  align_dp_cell** P = malloc((n+1) * sizeof(align_dp_cell *));
  if (P == NULL) {
    printf("Failed to allocate memory.\n");
    return NULL;
  }
  // We need an additional row/col in the beginning. Table indices are therefore
  // exactly one ahead of subproblem indices.
  for (int i = 0; i < n+1; i++) {
    P[i] = malloc((m+1) * sizeof(align_dp_cell));
    if (P[i] == NULL) {
      printf("Failed to allocate memory.\n");
      return NULL;
    }
    for (int j = 0; j < m+1; j++) {
      P[i][j] = (align_dp_cell) {i, j, 0, NULL};
    }
  }
  return P;
}

/**
 * Frees the allocated memory for a given alignment problem so that we can reuse
 * the same ::align_problem over and over.
 *
 * @param P the dynamic programming table.
 * @param row_cnt the number of rows in the table.
 * @param col_cnt the number of columns in the table.
 */
void free_dp_table(align_dp_cell** P, int row_cnt, int col_cnt) {
  if (P == NULL) {
    return;
  }
  int i,j;
  for (i = 0; i < row_cnt; i++) {
    for (j = 0; j < col_cnt; j++) {
      if (P[i][j].num_choices > 0) {
        free(P[i][j].choices);
      }
    }
    free(P[i]);
  }
  free(P);
}

/**
 * Given an alignment with allocated DP table, solves the alignment problem (i.e
 * populates the alignment table). The optimal score and transcript can then
 * be obtained by using `traceback`. Half of the constraints imposed by
 * an alignment type (see ::align_type) are implemented here where we decide
 * what positions on the table can be starting positions of alignments. The
 * other half of the constraints concerning the ending position of the alignment
 * on the table is encapsulated in `traceback`.
 *
 * The rules of starting alignments are:
 * - Local and end-anchored alignments (except for overlap ones) can start
 *   anywhere.
 * - Overlap alignments (except for start anchored ones) can start anywhere.
 * - Global and start anchored alignments must start at the top left corner.
 *
 * @param P The dynamic programming table,
 * @param def The alignment problem definition.
 * @return The optimal cell of the DP table for the alignment to *end* at.
 */
align_dp_cell* solve(align_dp_cell** P, align_problem* def) {
  int n = def->S_max_idx - def->S_min_idx;
  int m = def->T_max_idx - def->T_min_idx;
  double max_score, prev_score, max_prev_score, del_score, ins_score;
  int num_choices, max_prev_choice_idx, num_max_scores;
  int i,j,k;
  int *max_score_alts = NULL;
  align_choice *alts = NULL;
  int s,t;

  // Base case
  P[0][0].num_choices = 1;
  P[0][0].choices = malloc(sizeof(align_choice));
  if (P[0][0].choices == NULL) {
    printf("Failed to allocate memory.\n");
    return NULL;
  }
  P[0][0].choices[0].op = 'B';
  P[0][0].choices[0].diversion = 0;
  P[0][0].choices[0].score = 0;
  P[0][0].choices[0].base = NULL;

  // Populate the table
  for (i = 0; i < n+1; i++) {
    for (j = 0; j < m+1; j++) {
      if (i == 0 && j == 0) {
        // Already dealt with
        continue;
      }
      // Sane default so we can "continue" when a cell is not worth pursuing
      P[i][j].num_choices = 0;

      if (alts != NULL) {
        free(alts);
      }
      if (max_score_alts != NULL) {
        free(max_score_alts);
      }
      // Allocate for all 4 possible choices (B,M/S,I,D)
      alts = malloc(4 * sizeof(align_choice));
      if (alts == NULL) {
        printf("Failed to allocate memory.\n");
        return NULL;
      }
      // Build all the alternatives at cell (i,j)
      num_choices = 0;
      if (def->type == LOCAL || // local alignments can start anywhere
          def->type == END_ANCHORED || // end-anchored alignments can ...
          // Overlap alignments and end-anchored alignments can start anywhere
          // on either of the left or top edges:
          (def->type == OVERLAP && (i == 0 || j == 0)) ||
          (def->type == END_ANCHORED_OVERLAP && (i == 0 || j == 0))
        ) {
        alts[num_choices].op = 'B';
        alts[num_choices].score = 0;
        alts[num_choices].diversion = 0;
        alts[num_choices].base = NULL;

        num_choices++;
      }
      // the indices in the table are on ahead of the indices of letters:
      s = def->S[def->S_min_idx + i - 1];
      t = def->T[def->T_min_idx + j - 1];

      if (def->params->content_dependent_gap_scores == NULL) {
        del_score = def->params->gap_extend_score;
        ins_score = def->params->gap_extend_score;
      } else {
        // the indices in the table are on ahead of the indices of letters:
        del_score = def->params->content_dependent_gap_scores[s];
        ins_score = def->params->content_dependent_gap_scores[t];
      }
      // To (i-1,j)
      if (i > 0) {
        // Are there choices for (i-1,j) and are we inside the band?
        if (P[i-1][j].num_choices > 0 && (
              def->params->max_diversion < 0 ||
              // diversion of all choices in the same cell are identical:
              abs(P[i-1][j].choices[0].diversion - 1) <= def->params->max_diversion
            )) {
          max_prev_choice_idx = 0;
          max_prev_score = -INT_MAX;
          for (k = 0; k < P[i-1][j].num_choices; k++) {
            prev_score = P[i-1][j].choices[k].score + del_score;
            if (P[i-1][j].choices[k].op != 'D') {
              prev_score += def->params->gap_open_score;
            }
            if (prev_score > max_prev_score) {
              max_prev_score = prev_score;
              max_prev_choice_idx = k;
            }
          }
          alts[num_choices].op = 'D';
          alts[num_choices].diversion = P[i-1][j].choices[0].diversion - 1;
          alts[num_choices].score = max_prev_score;
          alts[num_choices].base = &(P[i-1][j].choices[max_prev_choice_idx]);

          num_choices++;
        }
      }
      // To (i,j-1)
      if (j > 0) {
        // Are there choices for (i,j-1) and are we inside the band?
        if (P[i][j-1].num_choices > 0 && (
              def->params->max_diversion < 0 ||
              // diversion of all choices in the same cell are identical:
              abs(P[i][j-1].choices[0].diversion + 1) <= def->params->max_diversion
            )) {
          max_prev_choice_idx = 0;
          max_prev_score = - INT_MAX;
          for (k = 0; k < P[i][j-1].num_choices; k++) {
            prev_score = P[i][j-1].choices[k].score + ins_score;
            if (P[i][j-1].choices[k].op != 'I') {
              prev_score += def->params->gap_open_score;
            }
            if (prev_score > max_prev_score) {
              max_prev_score = prev_score;
              max_prev_choice_idx = k;
            }
          }
          alts[num_choices].op = 'I';
          alts[num_choices].diversion = P[i][j-1].choices[0].diversion - 1;
          alts[num_choices].score = max_prev_score;
          alts[num_choices].base = &(P[i][j-1].choices[max_prev_choice_idx]);

          num_choices++;
        }
      }
      // To (i-1,j-1)
      if ((i > 0) && (j > 0)) {
        // Are there choices for (i-1,j-1)? (We must be in band if (i-1,j-1) is:
        if (P[i-1][j-1].num_choices > 0) {
          // All choices to (i-1,j-1) have the same score as there is no
          // gap open distinction, add the choice right away:
          if (s == t) {
            alts[num_choices].op = 'M';
          } else {
            alts[num_choices].op = 'S';
          }
          alts[num_choices].diversion = P[i-1][j-1].choices[0].diversion;
          alts[num_choices].score = P[i-1][j-1].choices[0].diversion;
          alts[num_choices].base = &(P[i-1][j-1].choices[0]);
          alts[num_choices].score = P[i-1][j-1].choices[0].score
            + def->params->subst_scores[s][t];

          num_choices++;
        }
      }

      // ================== Find the best alternatives ==================
      if (num_choices == 0) {
        continue;
      }

      // indices of maximum choices in the `alts' array
      max_score_alts = malloc(num_choices * sizeof(int));
      if (max_score_alts == NULL) {
        printf("Failed to allocate memory.\n");
        return NULL;
      }
      num_max_scores = 0;
      max_score = alts[0].score;
      for (k = 0; k < num_choices; k++){
        if (alts[k].score == max_score) {
          max_score_alts[num_max_scores] = k;
          num_max_scores++;
        }
        else if (alts[k].score > max_score) {
          max_score = alts[k].score;
          max_score_alts[0] = k;
          num_max_scores = 1;
        }
      }
      P[i][j].num_choices = num_max_scores;
      P[i][j].choices = malloc(num_max_scores * sizeof(align_choice));
      if (P[i][j].choices == NULL) {
        printf("Failed to allocate memory.\n");
        return NULL;
      }
      for (k = 0; k < num_max_scores; k++) {
        P[i][j].choices[k] = alts[max_score_alts[k]];
      }
    }
  }
  free(alts);
  free(max_score_alts);
  return find_optimal(P, def);
}

/**
 * Finds the optimal cell to start traceback from given a populated DP table
 * and an alignment problem definition (to know where to look for the optimal
 * cell). This is internally used by `solve` in the final step. The rules are
 * as follows:
 * - Global and end-anchored alignments must end at the bottom right corner.
 * - Overlap alignments (except end-anchored ones) can end anywhere on either
 *   of the bottom or right edges; the best is found.
 * - Local and start-anchored alignments (except for overlap ones) can end
 *   anywhere; the best is found.
 *
 * @param P  The *solved* (populated) dynamic programming table.
 * @param def The alignment problem definition.
 *
 * @return The optimal cell of the table for the alignment to *end* at.
 */
align_dp_cell* find_optimal(align_dp_cell** P, align_problem* def) {
  double max;
  int i,j;
  int row = -1, col = -1;
  if (def->type == GLOBAL ||
      def->type == END_ANCHORED ||
      def->type == END_ANCHORED_OVERLAP) {
    // Global and end-anchored alignments must end at the bottom right corner
    row = def->S_max_idx - def->S_min_idx;
    col = def->T_max_idx - def->T_min_idx;
    if (P[row][col].num_choices == 0) {
      return NULL;
    }
  }
  else if (def->type == OVERLAP || def->type == START_ANCHORED_OVERLAP) {
    // Overlap alignments (except end-anchored ones) can end anywhere on either
    // of the bottom or right edges; find the best:
    max = -INT_MAX;
    for (i = 0; i < def->S_max_idx - def->S_min_idx + 1; i++){
      for (j = 0; j < def->T_max_idx - def->T_min_idx + 1; j++) {
        // Are we on the bottom row or the right column?
        if (i != def->S_max_idx - def->S_min_idx &&
            j != def->T_max_idx - def->T_min_idx) {
          continue;
        }
        if (P[i][j].num_choices == 0) {
          continue;
        }
        if (P[i][j].choices[0].score > max) {
          row = i;
          col = j;
          max = P[i][j].choices[0].score;
        }
      }
    }
  }
  else if (def->type == LOCAL || def->type == START_ANCHORED) {
    // Local and start-anchored alignments (except for overlap ones) can end
    // anywhere; find the best:
    max = P[0][0].choices[0].score;
    for (i = 0; i < def->S_max_idx - def->S_min_idx + 1; i++){
      for (j = 0; j < def->T_max_idx - def->T_min_idx + 1; j++) {
        if (P[i][j].num_choices == 0) {
          continue;
        }
        if (P[i][j].choices[0].score > max) {
          row = i;
          col = j;
          max = P[i][j].choices[0].score;
        }
      }
    }
  }
  if (row == -1 || col == -1 || P[row][col].num_choices == 0) {
    return NULL;
  }
  return &(P[row][col]);
}

/**
 * Traces back *an* alignment (and not all alignments with identical scores)
 * from a given cell all the way back to a cell with a `NULL` base (i.e an
 * alignment start cell).
 *
 * @param P The solved (i.e populated) alignment DP table.
 * @param def The alignment problem definition.
 * @param end The desired ending point of alignment which becomes the starting
 *    point of traceback.
 *
 * @note
 *    Finding more than one optimal alignment is a nontrivial search problem,
 *    and requires some sort of global state keeping to avoid convoluted
 *    recursions. I couldn't get it right in the first go; leave for later.
 */
transcript* traceback(align_dp_cell** P, align_problem* def, align_dp_cell* end) {
  char op, *opseq;
  transcript* tx = malloc(sizeof(transcript));
  int S_idx = end->row,
      T_idx = end->col,
      len = S_idx + T_idx + 1,
      pos = len - 1;
  align_dp_cell cur = P[S_idx][T_idx];
  // We write ops to rev_opseq backwards starting from the end (position `len')
  char rev_opseq[len];
  while (1) {
    op = cur.choices[0].op;
    if (op == 'B') {
      break;
    }
    pos--;
    rev_opseq[pos] = op;
    if (op == 'M' || op == 'S') {
      S_idx --;
      T_idx --;
    }
    if (op == 'I') {
      T_idx --;
    }
    if (op == 'D') {
      S_idx --;
    }
    cur = P[S_idx][T_idx];
  }
  if (pos == len - 1) {
    // empty opseq
    return NULL;
  }

  len = len - pos - 1;
  opseq = malloc(len + 1);
  if (opseq == NULL) {
    printf("Failed to allocate memory.\n");
    return NULL;
  }
  strncpy(opseq, rev_opseq + pos, len);
  // strncpy does not null terminate:
  opseq[len] = '\0';

  tx->S_idx = S_idx + def->S_min_idx;
  tx->T_idx = T_idx + def->T_min_idx;
  tx->score = end->choices[0].score;
  tx->opseq = opseq;
  return tx;
}

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
 * @param params Alignment parameters to be used over the window.
 * @param window The length of the extension alignment window.
 * @param forward Either of 0 or 1 indicating the direction of extension.
 *
 * @return -1 if an error occurs and 0 otherwise.
 */
int extend_1d_once(segment* res, segment* seg,
  int* S, int* T, align_params* params,
  int window, int forward) {

  align_problem def;
  align_dp_cell **P, *opt;
  align_type type;
  transcript* tx;
  int tx_opseq_len, seg_opseq_len;
  char* opseq;
  int failure = 0;

  int S_len = tx_seq_len(seg->tx, 'S'),
      T_len = tx_seq_len(seg->tx, 'T');
  int S_min_idx, S_max_idx, T_min_idx, T_max_idx;
  if (forward) {
    S_min_idx = seg->tx->S_idx + S_len;
    T_min_idx = seg->tx->T_idx + T_len;
    S_max_idx = S_min_idx + window;
    T_max_idx = T_min_idx + window;
  } else {
    S_max_idx = seg->tx->S_idx;
    T_max_idx = seg->tx->T_idx;
    S_min_idx = S_max_idx - window;
    T_min_idx = T_max_idx - window;
  }

  type = forward ? START_ANCHORED_OVERLAP : END_ANCHORED_OVERLAP;
  def = (align_problem) {
    .S=S, .T=T, .type=type, .params=params,
    .S_min_idx=S_min_idx, .S_max_idx=S_max_idx,
    .T_min_idx=T_min_idx, .T_max_idx=T_max_idx,
  };
  P = init_dp_table(&def);
  if (P == NULL) {
    return -1;
  }
  opt = solve(P, &def);
  if (opt == NULL) {
    failure = 1;
  }

  tx = traceback(P, &def, opt);
  if (tx == NULL) {
    failure = 1;
  }
  free_dp_table(P, def.S_max_idx - def.S_min_idx + 1, def.T_max_idx - def.T_min_idx + 1);
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
 * @param params Alignment parameters to be used over the window.
 * @param window The length of the extension alignment window.
 * @param max_succ_drops Maximum number of "drops" until the
 *    segment is dropped (i.e -1 is returned).
 * @param drop_threshold What constitutes a drop in the score when comparing
 *    an (n+1)-th window extension and the n-th window extension.
 * @param forward Either of 0 or 1 indicating the direction of extension.
 *
 * @return 0 if the seed successfully extends to the boundary of either of
 *   the sequences and -1 otherwise.
 */
int extend_1d(segment* res, segment* seg, int* S, int* T, int S_len, int T_len,
  align_params* params, int window, int max_succ_drops, double drop_threshold,
  int forward) {

  segment cur_seg = *seg;
  double prev_score = seg->tx->score;
  int i, retcode, score_history[max_succ_drops];
  for (i = 0; i < max_succ_drops; i++) {
    score_history[i] = drop_threshold;
  }
  int S_end, T_end, S_wiggle, T_wiggle, min_wiggle, actual_window, drops;
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
      return 0;
    }

    retcode = extend_1d_once(&cur_seg, &cur_seg, S, T, params, actual_window, forward);
    if (retcode == -1) {
      // No nonempty alignment found:
      return -1;
    }

    for (i = 0; i < max_succ_drops - 1; i++) {
      score_history[i] = score_history[i+1];
    }
    // TODO is this correct?
    score_history[max_succ_drops-1] = cur_seg.tx->score - prev_score;
    drops = 0;
    for (i = 0; i < max_succ_drops; i ++) {
      if (score_history[i] < drop_threshold) {
        drops += 1;
      }
    }

    if (drops >= max_succ_drops) {
      return -1;
    }
    prev_score = cur_seg.tx->score;
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
 * @param params Alignment parameters to be used over the window.
 * @param window The length of the extension alignment window.
 * @param max_succ_drops Maximum number of "drops" until the
 *    segment is dropped (i.e -1 is returned).
 * @param drop_threshold What constitutes a drop in the score when comparing
 *    an (n+1)-th window extension and the n-th window extension.
 * @param forward Either of 0 or 1 indicating the direction of extension.
 *
 * @return 0 if the seed successfully extends to the boundary of either of
 *   the sequences and -1 otherwise.
 */
segment* extend(segment** segs, int num_segs, int* S, int* T, int S_len, int T_len,
  align_params* params, int window, int max_succ_drops, double drop_threshold) {

  segment fwd, bwd, *res;
  transcript* tx;
  char* opseq;
  int fwd_tx_len, bwd_tx_len, seg_tx_len, retcode;
  for (int i = 0; i < num_segs; i ++) {
    retcode = extend_1d(&fwd, segs[i],
        S, T, S_len, T_len, params,
        window, max_succ_drops, drop_threshold, 1);
    if (retcode == -1 || fwd.tx->score <= drop_threshold) {
      continue;
    }
    retcode = extend_1d(&bwd, segs[i],
      S, T, S_len, T_len, params,
      window, max_succ_drops, drop_threshold, 0);
    if (retcode == -1 || bwd.tx->score <= drop_threshold) {
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
    tx->score = bwd.tx->score + fwd.tx->score - segs[i]->tx->score;
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

/**
 * Translates a given sequence to integer indices of each letter
 * in the corresponding alphabet. This is the array that has to be provided
 * to all other functions, e.g solve().
 *
 * @param alphabet The alphabet over which the sequence is defined.
 * @param sequence The sequence as a single string, potentially having multiple
 *   characters per letter.
 */
int* idxseq_from_charseq(sequence_alphabet* alphabet, char* sequence) {
  int i,j, cur_idx, length = strlen(sequence) / alphabet->letter_length;
  int* idx_seq = malloc(length * sizeof(int));
  char* cur = malloc(alphabet->letter_length);
  for (i = 0; i < length; i++) {
    cur_idx = -1;
    snprintf(cur,
      alphabet->letter_length + 1,
      "%s", sequence + i*(alphabet->letter_length)
    );
    // TODO Would it make things better if the following was cached somehow?
    for (j = 0; j < alphabet->length; j++) {
      if (strcmp(alphabet->letters[j], cur) == 0) {
        cur_idx = j;
        break;
      }
    }
    if (cur_idx == -1) {
      printf("Invalid letter: %s\n", cur);
    }
    idx_seq[i] = cur_idx;
  }
  free(cur);
  return idx_seq;
}

/**
 * Helper method for debugging purposes. Only works on system with procfs.
 */
void _print_mem_usage() {
  int tot[1];
  int res[1];
  FILE *fp = fopen("/proc/self/statm", "r");
  fscanf(fp, "%d %d", tot, res);
  printf("%.2f MB (tot), %.2f MB (res)\n",
    tot[0]/1024.0, res[0]/1024.0);
  fclose(fp);
}
