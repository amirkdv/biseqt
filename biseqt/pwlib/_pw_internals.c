#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <limits.h>

#include "pwlib.h"

//FIXME docs


intpair _frame_dims(alnframe* frame) {
  return (intpair) {
    frame->S_range.j - frame->S_range.i,
    frame->T_range.j - frame->T_range.i
  };
}

bool _cellpos_valid(dptable* T, intpair pos) {
  if (pos.i < 0 || pos.j < 0 || pos.i >= T->table_dims.i) {
    return false;
  }
  if (T->prob->mode == STD_MODE && pos.j >= T->table_dims.j) {
    return false;
  }
  if (T->prob->mode == BANDED_MODE && pos.j >= T->row_lens[pos.i]) {
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

intpair _cellpos_from_xy(alnprob* prob, int x, int y) {
  if (prob->mode == STD_MODE) {
    return (intpair) {x, y};
  } else if (prob->mode == BANDED_MODE) {
    intpair pos = _da_from_xy(x, y);
    return (intpair) {pos.i + prob->banded_params->radius - prob->banded_params->ctrdiag, pos.i};
  } else {
    _panick("Unknown alignment mode given to _cellpos_from_xy()");
  }
  // for the compiler:
  return (intpair) {-1, -1};
}

intpair _xy_from_cellpos(alnprob* prob, int i, int j) {
  if (prob->mode == STD_MODE) {
    return (intpair) {i, j};
  } else if (prob->mode == BANDED_MODE) {
    i += prob->banded_params->ctrdiag - prob->banded_params->radius;
    return _xy_from_da(i, j);
  } else {
    _panick("Unknown alignment mode given to _cellpos_from_xy()");
  }
  // for the compiler:
  return (intpair) {-1, -1};
}

intpair _prev_cell_M(alnprob* prob, int i, int j) {
  switch (prob->mode) {
    case STD_MODE:
      return (intpair) {i - 1, j - 1};
    case BANDED_MODE:
      return (intpair) {i, j - 1};
    default:
      _panick("Invalid alignment mode");
  }
  // for the compiler:
  return (intpair) {-1, -1};
}

intpair _prev_cell_I(alnprob* prob, int i, int j) {
  switch (prob->mode) {
    case STD_MODE:
      return (intpair) {i, j - 1};
    case BANDED_MODE:
      if (j > prob->banded_params->radius) {
        return (intpair) {i + 1, j};
      } else {
        return (intpair) {i + 1, j - 1};
      }
    default:
      _panick("Invalid alignment mode");
  }
  // for the compiler:
  return (intpair) {-1, -1};
}

intpair _prev_cell_D(alnprob* prob, int i, int j) {
  switch (prob->mode) {
    case STD_MODE:
      return (intpair) {i - 1, j};
    case BANDED_MODE:
      if (j < prob->banded_params->radius) {
        return (intpair) {i - 1, j};
      } else {
        return (intpair) {i - 1, j - 1};
      }
    default:
      _panick("Invalid alignment mode");
  }
  // for the compiler:
  return (intpair) {-1, -1};
}

double _ge_score(alnprob* prob, intpair pos, char op) {
  if (prob->scores->content_dependent_gap_scores == NULL) {
    return prob->scores->gap_extend_score;
  }
  intpair xy = _xy_from_cellpos(prob, pos.i, pos.j);
  if (op == 'D') {
    return prob->scores->content_dependent_gap_scores[prob->frame->S_range.i + xy.i];
  }
  return prob->scores->content_dependent_gap_scores[prob->frame->T_range.i + xy.j];
}


/**
 * Gives the maximum possible value of `a` for a given shift in standard (d,a)
 * format.
 */
int _da_row_len(intpair dims, int d) {
  return 1 + (d > 0 ? 0 : d) + (dims.i - d > dims.j ? dims.j : dims.i - d);
}

/**
 * If possible populates an alignment 'B' move for a given position of the table
 * (in either coordinate systems depending on alignment mode).
 *
 * @return 0 if choice successfully built, -1 if choice not possible.
 */
int _alnchoice_B(dptable *T, intpair pos, alnchoice* choice) {
  // return -1 if the current cell is not elligible for a B:
  std_alntype std_type;
  switch (T->prob->mode) {
    case STD_MODE:
      std_type = T->prob->std_params->type;
      if (pos.i == 0 && pos.j == 0) {
        break;
      }
      if (std_type == LOCAL || std_type == END_ANCHORED) {
        break;
      }
      if (std_type == OVERLAP && (pos.i == 0 || pos.j == 0)) {
        break;
      }
      if (std_type == END_ANCHORED_OVERLAP && (pos.i == 0 || pos.j == 0)) {
        break;
      }
      // not elligible:
      return -1;
    case BANDED_MODE:
      if (pos.i == 0 && pos.j == 0) {
        break;
      }
      if (T->prob->banded_params->type == B_OVERLAP && pos.j == 0) {
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
  intpair chars; // the indices of characters in question of S and T,
  double score;
  alnframe* frame = T->prob->frame;
  intpair xy = _xy_from_cellpos(T->prob, pos.i, pos.j);
  intpair prev = _prev_cell_M(T->prob, pos.i, pos.j);

  if (!_cellpos_valid(T, prev) || T->cells[prev.i][prev.j].num_choices < 1) {
    return -1;
  }
  // pos is guaranteed to be in xy-system now:
  chars = (intpair) {
    frame->S[frame->S_range.i + xy.i - 1],
    frame->T[frame->T_range.i + xy.j - 1]
  };

  score = T->prob->scores->subst_scores[chars.i][chars.j];
  choice->op = (chars.i == chars.j ? 'M' : 'S');
  // No path-dependence for S/M; all previous bases are the same:
  choice->score = T->cells[prev.i][prev.j].choices[0].score + score;
  choice->base = &(T->cells[prev.i][prev.j].choices[0]);
  return 0;
}

int _alnchoice_ID(dptable* T, intpair pos, alnchoice* choice, char op) {
  int base_idx;
  double score, max_score, ge_score;
  intpair prev;

  switch (op) {
    case 'I':
      prev = _prev_cell_I(T->prob, pos.i, pos.j);
      break;
    case 'D':
      prev = _prev_cell_D(T->prob, pos.i, pos.j);
      break;
    default:
      _panick("Unknown op given to _alnchoice_ID\n");
  }
  if (!_cellpos_valid(T, prev)) {
    return -1;
  }
  if (T->cells[prev.i][prev.j].num_choices <= 0) {
    return -1;
  }

  ge_score = _ge_score(T->prob, pos, op);

  base_idx = 0;
  max_score = -INT_MAX;
  for (int k = 0; k < T->cells[prev.i][prev.j].num_choices; k++) {
    score = T->cells[prev.i][prev.j].choices[k].score + ge_score;
    if (T->cells[prev.i][prev.j].choices[k].op != op) {
      score += T->prob->scores->gap_open_score;
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
  std_alntype type = T->prob->std_params->type;
  alnframe* frame = T->prob->frame;
  double max;
  int i,j;
  int row = -1, col = -1;
  if (type == GLOBAL || type == END_ANCHORED || type == END_ANCHORED_OVERLAP) {
    // Global and end-anchored alignments must end at the bottom right corner
    row = frame->S_range.j - frame->S_range.i;
    col = frame->T_range.j - frame->T_range.i;
    if (T->cells[row][col].num_choices == 0) {
      return (intpair){-1, -1};
    }
  }
  else if (type == OVERLAP || type == START_ANCHORED_OVERLAP) {
    // Overlap alignments (except end-anchored ones) can end anywhere on either
    // of the bottom or right edges; find the best:
    max = -INT_MAX;
    for (i = 0; i < frame->S_range.j - frame->S_range.i + 1; i++){
      for (j = 0; j < frame->T_range.j - frame->T_range.i + 1; j++) {
        // Are we on the bottom row or the right column?
        if (i != frame->S_range.j - frame->S_range.i &&
            j != frame->T_range.j - frame->T_range.i) {
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
  else if (type == LOCAL || type == START_ANCHORED) {
    // Local and start-anchored alignments (except for overlap ones) can end
    // anywhere; find the best:
    max = T->cells[0][0].choices[0].score;
    for (i = 0; i < frame->S_range.j - frame->S_range.i + 1; i++){
      for (j = 0; j < frame->T_range.j - frame->T_range.i + 1; j++) {
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
  if (row == -1 || col == -1) {
    return (intpair){-1, -1};
  }
  return (intpair) {row, col};
}

/**
 */
intpair _banded_find_optimal(dptable* T) {
  alnframe* frame = T->prob->frame;
  double max;
  int i,j;
  int row = -1, col = -1;
  if (T->prob->banded_params->type == B_GLOBAL) {
    row = frame->S_range.j - frame->S_range.i;
    col = frame->T_range.j - frame->T_range.i;
  }
  else if (T->prob->banded_params->type == B_OVERLAP) {
    max = -INT_MAX;
    for (i = 0; i < T->table_dims.i; i++){
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
