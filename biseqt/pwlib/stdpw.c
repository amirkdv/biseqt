#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdio.h>
#include <string.h>

#include "pwlib.h"

/**
 * Given an alignment with allocated DP table, solves the alignment problem (i.e
 * populates the alignment table) and returns the optimal ending point for the
 * alignment.  The optimal score and transcript can then be obtained by using
 * `traceback`. Half of the constraints imposed by an alignment type are
 * implemented here where we decide what positions on the table can be starting
 * positions of alignments. The other half of the constraints concerning the
 * ending position of the alignment on the table is encapsulated in `traceback`.
 *
 * @param T The dynamic programming table.
 * @return The optimal cell for the alignment to end at or {-1,-1} if error.
 */
gridcoord stdpw_solve(dptable* T) {
  int num_choices, num_max_scores;
  int i,j,k;
  int *max_score_choices = NULL;
  double max_score;
  alnchoice *choices = NULL;

  // Populate the table
  for (i = 0; i < T->num_rows; i++) {
    for (j = 0; j < T->num_cols; j++) {
      if (choices != NULL) {
        free(choices);
      }
      if (max_score_choices != NULL) {
        free(max_score_choices);
      }
      // Allocate for all 4 possible choices (B,M/S,I,D)
      choices = malloc(4 * sizeof(alnchoice));
      if (choices == NULL) {
        printf("Failed to allocate memory (`stdpw_solve()`).\n");
        return (gridcoord){-1, -1};
      }
      // Build all the alternatives at cell (i,j)
      num_choices = 0;
      num_choices += (_alnchoice_B (T, i, j, &choices[num_choices])      == 0) ? 1 : 0;
      num_choices += (_alnchoice_ID(T, i, j, &choices[num_choices], 'D') == 0) ? 1 : 0;
      num_choices += (_alnchoice_ID(T, i, j, &choices[num_choices], 'I') == 0) ? 1 : 0;
      num_choices += (_alnchoice_MS(T, i, j, &choices[num_choices])      == 0) ? 1 : 0;

      // Find the best alternatives
      if (num_choices == 0) {
        T->cells[i][j].num_choices = 0;
        continue;
      }

      // indices of maximum choices in the `choices' array
      max_score_choices = malloc(num_choices * sizeof(int));
      if (max_score_choices == NULL) {
        printf("Failed to allocate memory (`stdpw_solve()`).\n");
        return (gridcoord){-1, -1};
      }
      num_max_scores = 0;
      max_score = choices[0].score;
      for (k = 0; k < num_choices; k++){
        if (choices[k].score == max_score) {
          max_score_choices[num_max_scores] = k;
          num_max_scores++;
        }
        else if (choices[k].score > max_score) {
          max_score = choices[k].score;
          max_score_choices[0] = k;
          num_max_scores = 1;
        }
      }
      T->cells[i][j].num_choices = num_max_scores;
      T->cells[i][j].choices = malloc(num_max_scores * sizeof(alnchoice));
      if (T->cells[i][j].choices == NULL) {
        printf("Failed to allocate memory (`stdpw_solve()`).\n");
        return (gridcoord){-1, -1};
      }
      for (k = 0; k < num_max_scores; k++) {
        T->cells[i][j].choices[k] = choices[max_score_choices[k]];
      }
    }
  }
  free(choices);
  free(max_score_choices);
  return stdpw_find_optimal(T);
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
 * @param T The *solved* dynamic programming table.
 *
 * @return The optimal cell for the alignment to end at or {-1,-1} if error.
 */
gridcoord stdpw_find_optimal(dptable* T) {
  std_alnprob* prob = T->std_prob;
  dpcell** P = T->cells;
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
      return (gridcoord){-1, -1};
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
    return (gridcoord){-1, -1};
  }
  return (gridcoord) {.row=row, .col=col};
}

/**
 * Traces back *an* alignment (and not all alignments with identical scores)
 * from a given cell all the way back to a cell with a `NULL` base (i.e an
 * alignment start cell).
 *
 * @param T The *solved* dynamic programming table.
 * @param end The desired ending point of alignment which becomes the starting
 *    point of traceback.
 *
 * @note
 *    Finding more than one optimal alignment is a nontrivial search problem,
 *    and requires some sort of global state keeping to avoid convoluted
 *    recursions. I couldn't get it right in the first go; leave for later.
 */
transcript* stdpw_traceback(dptable* T, gridcoord end) {
  std_alnprob* prob = T->std_prob;
  dpcell** P = T->cells;
  char op, *opseq;
  transcript* tx = malloc(sizeof(transcript));
  int S_idx = end.row,
      T_idx = end.col,
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
  tx->score = P[end.row][end.col].choices[0].score;
  tx->opseq = opseq;
  return tx;
}
