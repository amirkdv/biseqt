#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <limits.h>

#include "pwlib.h"

//FIXME docs


intpair _framedims(alnframe* frame) {
  return (intpair) {
    frame->S_max_idx - frame->S_min_idx,
    frame->T_max_idx - frame->T_min_idx
  };
}

bool _cellpos_valid(dptable* T, intpair pos) {
  if (pos.i < 0 || pos.j < 0 || pos.i >= T->num_rows) {
    return false;
  }
  if (T->mode == STD_MODE && pos.j >= T->num_cols) {
    return false;
  }
  if (T->mode == BANDED_MODE && pos.j >= T->row_lens[pos.i]) {
    return false;
  }
  return true;
}


/**
 * Converts the coordinates of a cell in the dynamic programming table from
 * the (d,a) system to the (x,y) system.
 *
 * Note: the first coordinate (d) is expected to be in the correct frame (not
 * adjusted to be positive for memory access).
 */
intpair _da_from_xy(int x, int y) {
  return (intpair) {x - y, x < y ? x : y};
}

/**
 * Converts the coordinates of a cell in the dynamic programming table from
 * the (x,y) system to the (d,a) system.
 *
 * Note: the first coordinate (d) must be adjusted based on band radius before
 * being used to access the dynamic programming table.
 */
intpair _xy_from_da(int d, int a) {
  return (intpair) {a + (d > 0 ? 0 : d), a - (d > 0 ? d : 0)};
}

intpair _cellpos_from_xy(int x, int y, alnmode mode, int T_len) {
  if (mode == STD_MODE) {
    return (intpair) {x, y};
  } else if (mode == BANDED_MODE) {
    intpair pos = _da_from_xy(x, y);
    return (intpair) {pos.i + T_len, pos.i};
  } else {
    _panick("Unknown alignment mode given to _cellpos_from_xy()");
    // for the compiler:
    return (intpair) {-1, -1};
  }
}

/**
 * Gives the maximum possible value of `a` for a given shift in standard (d,a)
 * format.
 */
int _da_row_len(intpair framedims, int d) {
  return 1 + (d > 0 ? 0 : d) + (framedims.i - d > framedims.j ? framedims.j : framedims.i - d);
}

/**
 * If possible populates an alignment 'B' move for a given position of the table
 * (in either coordinate systems depending on alignment mode).
 *
 * @return 0 if choice successfully built, -1 if choice not possible.
 */
int _alnchoice_B(dptable *T, intpair pos, alnchoice* choice) {
  // return -1 if the current cell is not elligible for a B:
  switch (T->mode) {
    case STD_MODE:
      if (pos.i == 0 && pos.j == 0) {
        break;
      }
      if (T->std_prob->type == LOCAL || T->std_prob->type == END_ANCHORED) {
        break;
      }
      if (T->std_prob->type == OVERLAP && (pos.i == 0 || pos.j == 0)) {
        break;
      }
      if (T->std_prob->type == END_ANCHORED_OVERLAP && (pos.i == 0 || pos.j == 0)) {
        break;
      }
      // not elligible:
      return -1;
    case BANDED_MODE:
      if (pos.i == 0 && pos.j == 0) {
        break;
      }
      if (T->banded_prob->type == B_OVERLAP && pos.j == 0) {
        break;
      }
      // not elligible:
      return -1;
    default:
      _panick("Invalid alignment mode");
      // for compiler's sake:
      return -1;
  }
  choice->op = 'B';
  choice->score = 0;
  choice->base = NULL;
  return 0;
}

/**
 * If possible populates an alignment 'M' or 'S' move for a given position of
 * the table (in either coordinate systems depending on alignment mode).
 *
 * @return 0 if choice successfully built, -1 if choice not possible.
 */
int _alnchoice_M(dptable *T, intpair pos, alnchoice* choice) {
  intpair prev,  // the position of base cell (in whichever system `pos` is),
            xy,  // the position of `pos` in xy-system,
         chars; // the indices of characters in question of S and T,
  double score;
  alnframe* frame;
  switch (T->mode) {
    case STD_MODE:
      frame = T->std_prob->frame;
      xy = pos;
      prev = (intpair) {pos.i - 1, pos.j - 1};
      break;
    case BANDED_MODE:
      frame = T->banded_prob->frame;
      prev = (intpair) {pos.i, pos.j-1};
      xy = _xy_from_da(pos.i - T->banded_prob->bradius, pos.j);
      break;
    default:
      _panick("Invalid alignment mode");
      // for compiler's sake:
      return -1;
  }
  if (!_cellpos_valid(T, prev) || T->cells[prev.i][prev.j].num_choices < 1) {
    return -1;
  }
  // pos is guaranteed to be in xy-system now:
  chars = (intpair) {
    frame->S[frame->S_min_idx + xy.i - 1],
    frame->T[frame->T_min_idx + xy.j - 1]
  };

  score = T->std_prob->scores->subst_scores[chars.i][chars.j];
  choice->op = (chars.i == chars.j ? 'M' : 'S');
  // No path-dependence for S/M; all previous bases are the same:
  choice->score = T->cells[prev.i][prev.j].choices[0].score + score;
  choice->base = &(T->cells[prev.i][prev.j].choices[0]);
  return 0;
}

intpair _prev_cell_I(dptable* T, intpair pos) {
  switch (T->mode) {
    case STD_MODE:
      return (intpair) {pos.i, pos.j - 1};
    case BANDED_MODE:
      if (pos.j > T->banded_prob->bradius) {
        return (intpair) {pos.i + 1, pos.j};
      } else {
        return (intpair) {pos.i + 1, pos.j - 1};
      }
    default:
      _panick("Invalid alignment mode");
      // for compiler's sake:
      return (intpair) {-1, -1};
  }
}

intpair _prev_cell_D(dptable* T, intpair pos) {
  switch (T->mode) {
    case STD_MODE:
      return (intpair) {pos.i - 1, pos.j};
    case BANDED_MODE:
      if (pos.j < T->banded_prob->bradius) {
        return (intpair) {pos.i - 1, pos.j};
      } else {
        return (intpair) {pos.i - 1, pos.j - 1};
      }
    default:
      _panick("Invalid alignment mode");
      // for compiler's sake:
      return (intpair) {-1, -1};
  }
}

double _ge_score(dptable *T, intpair pos, char op) {
  alnscores* scores = (T->mode == STD_MODE ? T->std_prob->scores : T->banded_prob->scores);
  if (scores->content_dependent_gap_scores == NULL) {
    return scores->gap_extend_score;
  } else {
    alnframe* frame = (T->mode == STD_MODE ? T->std_prob->frame : T->banded_prob->frame);
    intpair xy = (T->mode == STD_MODE ? pos : _xy_from_da(pos.i - T->banded_prob->bradius, pos.j));
    int content = (op == 'D' ? frame->S_min_idx + xy.i : frame->T_min_idx + xy.j);
    return scores->content_dependent_gap_scores[content];
  }
}

//FIXME this kind of silly thing would be resolved when refactoring alnprob
//happens.
double _go_score(dptable* T) {
  switch (T->mode) {
    case STD_MODE:
      return T->std_prob->scores->gap_open_score;
    case BANDED_MODE:
      return T->banded_prob->scores->gap_open_score;
    default:
      _panick("Invalid alignment mode");
      // for compiler's sake:
      return -1;
  }
}

int _alnchoice_ID(dptable* T, intpair pos, alnchoice* choice, char op) {
  int base_idx;
  double score, max_score;
  double ge_score, go_score;

  if (op != 'I' && op != 'D') {
    _panick("Unknown op given to _alnchoice_ID\n");
  }

  intpair prev = (op == 'I' ? _prev_cell_I(T, pos) : _prev_cell_D(T, pos));
  if (!_cellpos_valid(T, prev)) {
    return -1;
  }
  if (T->cells[prev.i][prev.j].num_choices <= 0) {
    return -1;
  }
  //FIXME
  /*printf("(%d, %d) ---%c--> (%d,%d)\n", pos.i, pos.j, op, prev.i, prev.j);*/
  /*exit(1);*/

  go_score = _go_score(T);
  ge_score = _ge_score(T, pos, op);

  base_idx = 0;
  max_score = -INT_MAX;
  for (int k = 0; k < T->cells[prev.i][prev.j].num_choices; k++) {
    score = T->cells[prev.i][prev.j].choices[k].score + ge_score;
    if (T->cells[prev.i][prev.j].choices[k].op != op) {
      score += go_score;
    }
    if (score > max_score) {
      max_score = score;
      base_idx = k;
    }
  }
  choice->op = op;
  choice->score = max_score;
  choice->base = &(T->cells[prev.i][prev.j].choices[base_idx]);
  return 0;
}

int _alnchoice_D(dptable* T, intpair pos, alnchoice* choice) {
  return _alnchoice_ID(T, pos, choice, 'D');
}

int _alnchoice_I(dptable* T, intpair pos, alnchoice* choice) {
  return _alnchoice_ID(T, pos, choice, 'I');
}

/**
 */
intpair _std_find_optimal(dptable* T) {
  std_alnprob* prob = T->std_prob;
  double max;
  int i,j;
  int row = -1, col = -1;
  if (prob->type == GLOBAL ||
      prob->type == END_ANCHORED ||
      prob->type == END_ANCHORED_OVERLAP) {
    // Global and end-anchored alignments must end at the bottom right corner
    row = prob->frame->S_max_idx - prob->frame->S_min_idx;
    col = prob->frame->T_max_idx - prob->frame->T_min_idx;
    if (T->cells[row][col].num_choices == 0) {
      return (intpair){-1, -1};
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
        if (T->cells[i][j].num_choices == 0) {
          continue;
        }
        if (T->cells[i][j].choices[0].score > max) {
          row = i;
          col = j;
          max = T->cells[i][j].choices[0].score;
        }
      }
    }
  }
  else if (prob->type == LOCAL || prob->type == START_ANCHORED) {
    // Local and start-anchored alignments (except for overlap ones) can end
    // anywhere; find the best:
    max = T->cells[0][0].choices[0].score;
    for (i = 0; i < prob->frame->S_max_idx - prob->frame->S_min_idx + 1; i++){
      for (j = 0; j < prob->frame->T_max_idx - prob->frame->T_min_idx + 1; j++) {
        if (T->cells[i][j].num_choices == 0) {
          continue;
        }
        if (T->cells[i][j].choices[0].score > max) {
          row = i;
          col = j;
          max = T->cells[i][j].choices[0].score;
        }
      }
    }
  }
  if (row == -1 || col == -1 || T->cells[row][col].num_choices == 0) {
    return (intpair){-1, -1};
  }
  return (intpair) {row, col};
}

/**
 */
intpair _banded_find_optimal(dptable* T) {
  alnframe* frame = T->banded_prob->frame;
  double max;
  int i,j;
  int row = -1, col = -1;
  if (T->banded_prob->type == B_GLOBAL) {
    row = frame->S_max_idx - frame->S_min_idx;
    col = frame->T_max_idx - frame->T_min_idx;
  }
  else if (T->banded_prob->type == B_OVERLAP) {
    max = -INT_MAX;
    for (i = 0; i < T->num_rows; i++){
      j = T->row_lens[i];
      if (T->cells[i][j].num_choices == 0) {
        continue;
      }
      if (T->cells[i][j].choices[0].score > max) {
        row = i;
        col = j;
        max = T->cells[i][j].choices[0].score;
      }
    }
  }
  if (row == -1 || col == -1 || T->cells[row][col].num_choices == 0) {
    return (intpair){-1, -1};
  }
  return (intpair) {row, col};
}

/**
 * Helper method for debugging purposes. Only works on systems with procfs.
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

void _panick(char* message) {
  printf("Panick: %s\n", message);
  exit(1);
}
