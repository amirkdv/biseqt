#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdio.h>
#include <string.h>

#include "pwlib.h"

/**
 * Given an alignment problem definition, creates and initializes the dynamic
 * programming table.
 *
 * @param prob
 *    Alignment problem definition. The only members that are relevant at this
 *    point is the corresponding frames of the two sequences which determine
 *    the size of the table).
 *
 * @return pointer to the `malloc` ed dynamic programming table.
 */
dpcell** stdpw_init(std_alnprob* prob) {
  int n = prob->frame->S_max_idx - prob->frame->S_min_idx;
  int m = prob->frame->T_max_idx - prob->frame->T_min_idx;
  dpcell** P = malloc((n+1) * sizeof(dpcell *));
  if (P == NULL) {
    printf("Failed to allocate memory (`init_dp_table()`).\n");
    return NULL;
  }
  // We need an additional row/col in the beginning. Table indices are therefore
  // exactly one ahead of subproblem indices.
  for (int i = 0; i < n+1; i++) {
    P[i] = malloc((m+1) * sizeof(dpcell));
    if (P[i] == NULL) {
      printf("Failed to allocate memory (`init_dp_table()`).\n");
      return NULL;
    }
    for (int j = 0; j < m+1; j++) {
      P[i][j] = (dpcell) {i, j, 0, NULL};
    }
  }
  return P;
}

/**
 * Frees the allocated memory for a given alignment problem so that we can reuse
 * the same ::std_alnprob over and over.
 *
 * @param P the dynamic programming table.
 * @param row_cnt the number of rows in the table.
 * @param col_cnt the number of columns in the table.
 */
void stdpw_free(dpcell** P, int row_cnt, int col_cnt) {
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
 * an alignment type are implemented here where we decide
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
 * @param prob The alignment problem definition.
 * @return The optimal cell of the DP table for the alignment to *end* at.
 */
dpcell* stdpw_solve(dpcell** P, std_alnprob* prob) {
  int n = prob->frame->S_max_idx - prob->frame->S_min_idx;
  int m = prob->frame->T_max_idx - prob->frame->T_min_idx;
  double max_score, prev_score, max_prev_score, del_score, ins_score;
  int num_choices, max_prev_choice_idx, num_max_scores;
  int i,j,k;
  int *max_score_alts = NULL;
  alnchoice *alts = NULL;
  int s,t;

  // Base case
  P[0][0].num_choices = 1;
  P[0][0].choices = malloc(sizeof(alnchoice));
  if (P[0][0].choices == NULL) {
    printf("Failed to allocate memory (`std_solve()`).\n");
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
      alts = malloc(4 * sizeof(alnchoice));
      if (alts == NULL) {
        printf("Failed to allocate memory (`std_solve()`).\n");
        return NULL;
      }
      // Build all the alternatives at cell (i,j)
      num_choices = 0;
      if (prob->type == LOCAL || // local alignments can start anywhere
          prob->type == END_ANCHORED || // end-anchored alignments can ...
          // Overlap alignments and end-anchored alignments can start anywhere
          // on either of the left or top edges:
          (prob->type == OVERLAP && (i == 0 || j == 0)) ||
          (prob->type == END_ANCHORED_OVERLAP && (i == 0 || j == 0))
        ) {
        alts[num_choices].op = 'B';
        alts[num_choices].score = 0;
        alts[num_choices].diversion = 0;
        alts[num_choices].base = NULL;

        num_choices++;
      }
      // the indices in the table are on ahead of the indices of letters:
      s = prob->frame->S[prob->frame->S_min_idx + i - 1];
      t = prob->frame->T[prob->frame->T_min_idx + j - 1];

      if (prob->params->content_dependent_gap_scores == NULL) {
        del_score = prob->params->gap_extend_score;
        ins_score = prob->params->gap_extend_score;
      } else {
        // the indices in the table are on ahead of the indices of letters:
        del_score = prob->params->content_dependent_gap_scores[s];
        ins_score = prob->params->content_dependent_gap_scores[t];
      }
      // To (i-1,j)
      if (i > 0) {
        // Are there choices for (i-1,j) and are we inside the band?
        if (P[i-1][j].num_choices > 0 && (
              prob->bradius < 0 ||
              // diversion of all choices in the same cell are identical:
              abs(P[i-1][j].choices[0].diversion - 1) <= prob->bradius
            )) {
          max_prev_choice_idx = 0;
          max_prev_score = -INT_MAX;
          for (k = 0; k < P[i-1][j].num_choices; k++) {
            prev_score = P[i-1][j].choices[k].score + del_score;
            if (P[i-1][j].choices[k].op != 'D') {
              prev_score += prob->params->gap_open_score;
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
              prob->bradius < 0 ||
              // diversion of all choices in the same cell are identical:
              abs(P[i][j-1].choices[0].diversion + 1) <= prob->bradius
            )) {
          max_prev_choice_idx = 0;
          max_prev_score = - INT_MAX;
          for (k = 0; k < P[i][j-1].num_choices; k++) {
            prev_score = P[i][j-1].choices[k].score + ins_score;
            if (P[i][j-1].choices[k].op != 'I') {
              prev_score += prob->params->gap_open_score;
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
            + prob->params->subst_scores[s][t];

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
        printf("Failed to allocate memory (`std_solve()`).\n");
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
      P[i][j].choices = malloc(num_max_scores * sizeof(alnchoice));
      if (P[i][j].choices == NULL) {
        printf("Failed to allocate memory (`std_solve()`).\n");
        return NULL;
      }
      for (k = 0; k < num_max_scores; k++) {
        P[i][j].choices[k] = alts[max_score_alts[k]];
      }
    }
  }
  free(alts);
  free(max_score_alts);
  return stdpw_find_optimal(P, prob);
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
 * @param prob The alignment problem definition.
 *
 * @return The optimal cell of the table for the alignment to *end* at.
 */
dpcell* stdpw_find_optimal(dpcell** P, std_alnprob* prob) {
  double max;
  int i,j;
  int row = -1, col = -1;
  if (prob->type == GLOBAL ||
      prob->type == END_ANCHORED ||
      prob->type == END_ANCHORED_OVERLAP) {
    // Global and end-anchored alignments must end at the bottom right corner
    row = prob->frame->S_max_idx - prob->frame->S_min_idx;
    col = prob->frame->T_max_idx - prob->frame->T_min_idx;
    if (P[row][col].num_choices == 0) {
      return NULL;
    }
  }
  else if (prob->type == OVERLAP || prob->type == START_ANCHORED_OVERLAP) {
    // Overlap alignments (except end-anchored ones) can end anywhere on either
    // of the bottom or right edges; find the best:
    max = -INT_MAX;
    for (i = 0; i < prob->frame->S_max_idx - prob->frame->S_min_idx + 1; i++){
      for (j = 0; j < prob->frame->T_max_idx - prob->frame->T_min_idx + 1; j++) {
        // Are we on the bottom row or the right column?
        if (i != prob->frame->S_max_idx - prob->frame->S_min_idx &&
            j != prob->frame->T_max_idx - prob->frame->T_min_idx) {
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
  else if (prob->type == LOCAL || prob->type == START_ANCHORED) {
    // Local and start-anchored alignments (except for overlap ones) can end
    // anywhere; find the best:
    max = P[0][0].choices[0].score;
    for (i = 0; i < prob->frame->S_max_idx - prob->frame->S_min_idx + 1; i++){
      for (j = 0; j < prob->frame->T_max_idx - prob->frame->T_min_idx + 1; j++) {
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
 * @param prob The alignment problem definition.
 * @param end The desired ending point of alignment which becomes the starting
 *    point of traceback.
 *
 * @note
 *    Finding more than one optimal alignment is a nontrivial search problem,
 *    and requires some sort of global state keeping to avoid convoluted
 *    recursions. I couldn't get it right in the first go; leave for later.
 */
transcript* stdpw_traceback(dpcell** P, std_alnprob* prob, dpcell* end) {
  char op, *opseq;
  transcript* tx = malloc(sizeof(transcript));
  int S_idx = end->row,
      T_idx = end->col,
      len = S_idx + T_idx + 1,
      pos = len - 1;
  dpcell cur = P[S_idx][T_idx];
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
    printf("Failed to allocate memory (`tracback()`).\n");
    return NULL;
  }
  strncpy(opseq, rev_opseq + pos, len);
  // strncpy does not null terminate:
  opseq[len] = '\0';

  tx->S_idx = S_idx + prob->frame->S_min_idx;
  tx->T_idx = T_idx + prob->frame->T_min_idx;
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
