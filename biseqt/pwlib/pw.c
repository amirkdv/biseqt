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
int dptable_init(dptable* T) {
  switch(T->prob->mode) {
    case STD_MODE:
      if (_std_table_init_dims(T) == -1) {
        return -1;
      }
      break;
    case BANDED_MODE:
      if (_banded_table_init_dims(T) == -1) {
        return -1;
      }
      break;
    default:
      _panick("Shouldn't have happened!");
  }
  return _table_init_cells(T);
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
  for (i = 0; i < T->num_rows; i++) {
    for (j = 0; j < T->row_lens[i]; j++) {
      if (T->cells[i][j].num_choices > 0) {
        free(T->cells[i][j].choices);
      }
    }
    free(T->cells[i]);
  }
  if (T->prob->mode == BANDED_MODE) {
    free(T->row_lens);
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
  int x, y, k;
  int *max_score_choices = NULL;
  double max_score;
  alnchoice *choices = NULL;
  intpair ylim, xlim = _xlim(T->prob);
  // It's simpler to always populate the table in the natural order of xy
  // coordinates (note that order of dependencies is obviously statisfied in any
  // coordinate system).
  for (x = xlim.i; x < xlim.j; x++) {
    ylim = _ylim(T->prob, x);
    for (y = ylim.i; y < ylim.j; y++) {
      // the coordinates of our cell in the dynamic programming table.
      if (choices != NULL) {
        free(choices);
        choices = NULL;
      }
      if (max_score_choices != NULL) {
        free(max_score_choices);
        max_score_choices = NULL;
      }
      // Allocate for all 4 possible choices (B,M/S,I,D)
      choices = malloc(4 * sizeof(alnchoice));
      if (choices == NULL) {
        _panick("Failed to allocate memory.");
      }

      // Find all possible moves: _alnchoice_X functions populate the choice*
      // they're given and return 0 if move is allowed.
      num_choices = 0;
      num_choices += (_alnchoice_B(T, x, y, &choices[num_choices]) == 0) ? 1 : 0;
      num_choices += (_alnchoice_D(T, x, y, &choices[num_choices]) == 0) ? 1 : 0;
      num_choices += (_alnchoice_I(T, x, y, &choices[num_choices]) == 0) ? 1 : 0;
      num_choices += (_alnchoice_M(T, x, y, &choices[num_choices]) == 0) ? 1 : 0;

      cellpos = _cellpos_from_xy(T->prob, x, y);
      if (num_choices == 0) {
        T->cells[cellpos.i][cellpos.j].num_choices = 0;
        continue;
      }

      // Find the highest scoring alternatives
      max_score_choices = malloc(num_choices * sizeof(int));
      if (max_score_choices == NULL) {
        _panick("Failed to allocate memory.");
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
        _panick("Failed to allocate memory.");
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
  if (tx == NULL) {
    _panick("Failed to allocate memory.");
  }
  intpair xy = _xy_from_cellpos(T->prob, end.i, end.j);
  int len = xy.i + xy.j + 1, // FIXME not unnecessarily too big?
      pos = len - 1;
  alnchoice* cur = &(T->cells[end.i][end.j].choices[0]);
  // We write ops to rev_opseq backwards starting from the end (position `len')
  char rev_opseq[len];
  while (cur->base != NULL) {
    pos--;
    rev_opseq[pos] = cur->op;
    xy.i -= (cur->op == 'I' ? 0 : 1);
    xy.j -= (cur->op == 'D' ? 0 : 1);
    cur = cur->base;
  }
  if (pos == len - 1) {
    // empty opseq
    return NULL;
  }

  len = len - pos - 1;
  opseq = malloc(len + 1);
  if (opseq == NULL) {
    _panick("Failed to allocate memory.");
  }
  strncpy(opseq, rev_opseq + pos, len);
  // strncpy does not null terminate:
  opseq[len] = '\0';

  tx->S_idx = xy.i + T->prob->frame->S_range.i;
  tx->T_idx = xy.j + T->prob->frame->T_range.i;
  tx->score=T->cells[end.i][end.j].choices[0].score;
  tx->opseq=opseq;
  return tx;
}
