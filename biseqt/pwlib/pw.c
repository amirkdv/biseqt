#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdio.h>
#include <string.h>

#include "pwlib.h"

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

alignment* dptable_traceback(dptable* T, intpair end) {
  char *opseq;
  alignment* aln = malloc(sizeof(alignment));
  if (aln == NULL) {
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

  aln->S_idx = xy.i + T->prob->frame->S_range.i;
  aln->T_idx = xy.j + T->prob->frame->T_range.i;
  aln->score=T->cells[end.i][end.j].choices[0].score;
  aln->opseq=opseq;
  return aln;
}
