#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdio.h>
#include <string.h>

#include "pwlib.h"

//FIXME docs

/**
 * Given a ::dptable containing the alignment problem definition and the mode of
 * alignment, the table's dimensions are calculated and the appropriate amount of
 * memory is allocated with initialized ::dpcell entries.
 *
 * In banded mode the region of interest, which is not rectangular in the
 * default (x,y) coordinate system, is mapped to an alternative coordinate
 * system (refered to as (d,a)-coordinates).
 *
 * @param T
 *    Table to be initialized. The only fields that must be
 *    already populated are the alignment mode and the problem definition.
 *
 * @return 0 if successful and -1 if an error occurs.
 */
// FIXME move this and the next one to pw.c
int dptable_init(dptable* T) {
  int i, j;
  intpair frame_dims = _frame_dims(T->prob->frame);
  switch (T->prob->mode) {
    case STD_MODE:
      T->table_dims = (intpair) {frame_dims.i + 1, frame_dims.j + 1};
      break;
    case BANDED_MODE:
      T->table_dims = (intpair) {
        1 + 2 * T->prob->banded_params->radius,
        1 + (frame_dims.i < frame_dims.j ? frame_dims.i : frame_dims.j)
      };
      T->row_lens = malloc(T->table_dims.i * sizeof(int));
      for (i = 0; i < T->table_dims.i; i++) {
        T->row_lens[i] = _da_row_len(frame_dims, i + T->prob->banded_params->radius);
      }
      break;
  }
  T->cells = malloc(T->table_dims.i * sizeof(dpcell *));
  if (T->cells == NULL) {
    printf("Failed to allocate memory (`dptable_init()`).\n");
    return -1;
  }
  for (i = 0; i < T->table_dims.i; i++) {
    T->cells[i] = malloc(T->table_dims.j * sizeof(dpcell));
    if (T->cells[i] == NULL) {
      printf("Failed to allocate memory (`dptable_init()`).\n");
      return -1;
    }
    for (j = 0; j < T->table_dims.j; j++) {
      T->cells[i][j] = (dpcell) {.num_choices=0, .choices=NULL};
    }
  }
  return 0;
}

/**
 * Frees the allocated memory for the cells of a given ::dptable.
 *
 * @param T the dynamic programming table containing all the info.
 */
void dptable_free(dptable* T) {
  if (T == NULL) {
    return;
  }
  int i,j;
  for (i = 0; i < T->table_dims.i; i++) {
    for (j = 0; j < T->table_dims.j; j++) {
      if (T->cells[i][j].num_choices > 0) {
        free(T->cells[i][j].choices);
      }
    }
    free(T->cells[i]);
  }
  free(T->cells);
}

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
intpair dptable_solve(dptable* T) {
  int num_choices, num_max_scores;
  intpair cellpos; // dpcell position in the table, could be either xy/da
  int i,j,k;
  int *max_score_choices = NULL;
  double max_score;
  alnchoice *choices = NULL;
  intpair frame_dims = _frame_dims(T->prob->frame);

  // It's simpler to always populate the table in the natural order of xy
  // coordinates (note that order of dependencies is obviously statisfied in any
  // coordinate system).
  for (i = 0; i <= frame_dims.i; i++) {
    for (j = 0; j <= frame_dims.j; j++) {
      // the coordinates of our cell in the coordinate system in use: xy or da.
      cellpos = _cellpos_from_xy(T->prob, i, j);
      if (choices != NULL) {
        free(choices);
      }
      if (max_score_choices != NULL) {
        free(max_score_choices);
      }
      // Allocate for all 4 possible choices (B,M/S,I,D)
      choices = malloc(4 * sizeof(alnchoice));
      if (choices == NULL) {
        printf("Failed to allocate memory (`dptable_solve()`).\n");
        return (intpair){-1, -1};
      }
      // Find all possible moves: _alnchoice_X functions populate the choice*
      // they're given and return 0 if move is allowed.
      num_choices = 0;
      num_choices += (_alnchoice_B(T, cellpos, &choices[num_choices]) == 0) ? 1 : 0;
      num_choices += (_alnchoice_D(T, cellpos, &choices[num_choices]) == 0) ? 1 : 0;
      num_choices += (_alnchoice_I(T, cellpos, &choices[num_choices]) == 0) ? 1 : 0;
      num_choices += (_alnchoice_M(T, cellpos, &choices[num_choices]) == 0) ? 1 : 0;

      if (num_choices == 0) {
        T->cells[i][j].num_choices = 0;
        continue;
      }

      // Find the highest scoring alternatives
      max_score_choices = malloc(num_choices * sizeof(int));
      if (max_score_choices == NULL) {
        printf("Failed to allocate memory (`dptable_solve()`).\n");
        return (intpair){-1, -1};
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
      T->cells[cellpos.i][cellpos.j].num_choices = num_max_scores;
      T->cells[cellpos.i][cellpos.j].choices = malloc(num_max_scores * sizeof(alnchoice));
      if (T->cells[cellpos.i][cellpos.j].choices == NULL) {
        printf("Failed to allocate memory (`dptable_solve()`).\n");
        return (intpair){-1, -1};
      }
      for (k = 0; k < num_max_scores; k++) {
        T->cells[cellpos.i][cellpos.j].choices[k] = choices[max_score_choices[k]];
      }
    }
  }
  free(choices);
  free(max_score_choices);
  return (T->prob->mode == STD_MODE ? _std_find_optimal(T) : _banded_find_optimal(T));
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
transcript* dptable_traceback(dptable* T, intpair end) {
  char *opseq;
  transcript* tx = malloc(sizeof(transcript));
  int len = end.i + end.j + 1, // FIXME not unnecessarily too big?
      pos = len - 1;
  alnchoice* cur = &(T->cells[end.i][end.j].choices[0]);
  intpair starts = _xy_from_cellpos(T->prob, end.i, end.j);
  // We write ops to rev_opseq backwards starting from the end (position `len')
  char rev_opseq[len];
  while (cur->base != NULL) {
    pos--;
    rev_opseq[pos] = cur->op;
    starts.i -= (cur->op == 'I' ? 0 : 1);
    starts.j -= (cur->op == 'D' ? 0 : 1);
    cur = cur->base;
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

  tx->S_idx = starts.i + T->prob->frame->S_range.i;
  tx->T_idx = starts.j + T->prob->frame->T_range.i;
  tx->score=T->cells[end.i][end.j].choices[0].score;
  tx->opseq=opseq;
  return tx;
}
