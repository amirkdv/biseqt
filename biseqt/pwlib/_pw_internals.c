#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <limits.h>

#include "pwlib.h"

int _std_table_init_dims(dptable* T) {
  int xmax = T->prob->frame->origin_range.j - T->prob->frame->origin_range.i,
      ymax = T->prob->frame->mutant_range.j - T->prob->frame->mutant_range.i;
  int i;
  T->num_rows = xmax + 1;

  // Caclculate row lengths:
  T->row_lens = malloc(T->num_rows * sizeof(int));
  if (T->row_lens == NULL) {
    _panick("Failed to allocated memory.");
  }
  for (i = 0; i < T->num_rows; i++) {
    T->row_lens[i] = ymax + 1;
  }
  return 0;
}

int _banded_table_init_dims(dptable* T) {
  int xmax = T->prob->frame->origin_range.j - T->prob->frame->origin_range.i,
      ymax = T->prob->frame->mutant_range.j - T->prob->frame->mutant_range.i,
      dend = xmax - ymax,
      dmax = T->prob->banded_params->dmax,
      dmin = T->prob->banded_params->dmin,
      d, i;
  if (dmax > xmax || dmin < -ymax) {
    dmax = (dmax > xmax ? xmax : dmax);
    dmin = (dmin < -ymax ? -ymax: dmin);
    printf("Band [%d, %d] exceeds table limits, reduced it to [%d, %d].\n",
      T->prob->banded_params->dmin, T->prob->banded_params->dmax, dmin, dmax);
    T->prob->banded_params->dmax = dmax;
    T->prob->banded_params->dmin = dmin;
  }
  // in global mode, make sure the end points are in band:
  if (T->prob->banded_params->type == B_GLOBAL &&
      (dend > dmax || dend < dmin || dmax * dmin > 0)
    ){
    printf("End points not within band for global alignment!\n");
    return -1;
  }
  T->num_rows = 1 + dmax - dmin;

  if (T->num_rows < 0) {
    printf("Invalid band: [%d, %d]!\n", dmin, dmax);
    return -1;
  }

  // Caclculate row lengths:
  T->row_lens = malloc(T->num_rows * sizeof(int));
  if (T->row_lens == NULL) {
    _panick("Failed to allocated memory.");
  }
  for (i = 0; i < T->num_rows; i++) {
    // the actual shift is dmin + i since i is the adjusted to [0, dmax-dmin].
    d = dmin + i;
    T->row_lens[i] = 1 + (d > 0 ? 0 : d) + (xmax - d > ymax ? ymax : xmax - d);
    if (T->row_lens[i] <= 0) {
      _panick("This shouldn't have happened: row length is negative!");
    }
  }
  return 0;
}

int _table_init_cells(dptable* T) {
  int i, j;
  T->cells = malloc(T->num_rows * sizeof(dpcell *));
  if (T->cells == NULL) {
    _panick("Failed to allocated memory.");
  }
  for (i = 0; i < T->num_rows; i++) {
    T->cells[i] = malloc(T->row_lens[i] * sizeof(dpcell));
    if (T->cells[i] == NULL) {
      _panick("Failed to allocated memoryx.");
    }
    for (j = 0; j < T->row_lens[i]; j++) {
      T->cells[i][j] = (dpcell) {.num_choices=0, .choices=NULL};
    }
  }
  return 0;
}

bool _cellpos_valid(dptable* T, intpair pos) {
  if (pos.i < 0 ||
      pos.j < 0 ||
      pos.i >= T->num_rows ||
      pos.j >= T->row_lens[pos.i]
    ) {
    return false;
  }
  return true;
}

intpair _cellpos_from_xy(alnprob* prob, int x, int y) {
  if (prob->mode == STD_MODE) {
    return (intpair) {x, y};
  } else if (prob->mode == BANDED_MODE) {
    intpair pos = (intpair) {x - y, x < y ? x : y};
    // pos.i is the actual d, we want it adjusted to [0, dmax-dmin].
    return (intpair) {pos.i - prob->banded_params->dmin, pos.j};
  } else {
    _panick("Unknown alignment mode.");
  }
  return (intpair) {-1, -1}; // for the compiler
}

intpair _xy_from_cellpos(alnprob* prob, int i, int j) {
  switch (prob->mode) {
    case STD_MODE:
      return (intpair) {i, j};
    case BANDED_MODE:
      // i is the adjusted d and j is a:
      i += prob->banded_params->dmin;
      // is is now the actual d and the formula is:
      // (x,y) = (a+max(d,0), a-min(d,0))
      return (intpair) {j + (i > 0 ? i : 0), j - (i > 0 ? 0 : i)};
    default:
      _panick("Unknown alignment mode.");
  }
  return (intpair) {-1, -1}; // for the compiler
}

intpair _xlim(alnprob* prob) {
  int dmin, dmax;
  int origin_len = prob->frame->origin_range.j - prob->frame->origin_range.i,
      mutant_len = prob->frame->mutant_range.j - prob->frame->mutant_range.i;
  switch (prob->mode) {
    case STD_MODE:
      return (intpair) {0, origin_len + 1};
    case BANDED_MODE:
      dmin = prob->banded_params->dmin;
      dmax = prob->banded_params->dmax;
      return (intpair) {
        dmin > 0 ? dmin : 0,
        1 + (origin_len > mutant_len + dmax ? mutant_len + dmax : origin_len)
      };
    default:
      _panick("Unknown alignment mode.");
  }
  return (intpair) {-1, -1}; // for the compiler
}

intpair _ylim(alnprob* prob, int x) {
  int dmin, dmax;
  int mutant_len = prob->frame->mutant_range.j - prob->frame->mutant_range.i;
  switch (prob->mode) {
    case STD_MODE:
      return (intpair) {0, mutant_len + 1};
    case BANDED_MODE:
      dmin = prob->banded_params->dmin;
      dmax = prob->banded_params->dmax;
      return (intpair) {
        x - dmax > 0 ? x - dmax : 0,
        1 + (mutant_len > x - dmin ? x - dmin : mutant_len)
      };
    default:
      _panick("Unknown alignment mode.");
  }
  return (intpair) {-1, -1}; // for the compiler
}

/**
 * If possible populates an alignment 'B' move for a given position of the table
 * (in either coordinate systems depending on alignment mode).
 *
 * @return 0 if choice successfully built, -1 if choice not possible.
 */
int _alnchoice_B(dptable *T, int x, int y, alnchoice* choice) {
  std_alntype std_type;
  banded_alntype banded_type;
  intpair cellpos;
  switch (T->prob->mode) {
    case STD_MODE:
      std_type = T->prob->std_params->type;
      // the origin is always allowed as a start position:
      if (x == 0 && y == 0) {
        break;
      }
      // local and end-anchored alignments can start anywhere:
      if (std_type == LOCAL || std_type == END_ANCHORED) {
        break;
      }
      // overlap and end-anchored overlap alignments must start on the edges:
      if ((std_type == OVERLAP || std_type == END_ANCHORED_OVERLAP) &&
          (x == 0 || y == 0)) {
        break;
      }
      // B not allowed otherwise:
      return -1;
    case BANDED_MODE:
      cellpos = _cellpos_from_xy(T->prob, x, y);
      banded_type = T->prob->banded_params->type;
      // global alignments must start at the origin:
      if (banded_type == B_GLOBAL && x == 0 && y == 0) {
        break;
      }
      // overlap alignments can start anywhere with a = 0
      if (banded_type == B_OVERLAP && cellpos.j == 0) {
        break;
      }
      // B not allowed otherwise:
      return -1;
    default:
      _panick("Invalid alignment mode");
      // for compiler's sake:
      return -1;
  }
  choice->mins_cd = T->prob->max_new_mins;
  choice->cur_min = -INT_MAX;
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
int _alnchoice_M(dptable *T, int x, int y, alnchoice* choice) {
  intpair chars; // the indices of characters in question of origin and T,
  double score;
  alnframe* frame = T->prob->frame;
  intpair prev = _cellpos_from_xy(T->prob, x-1, y-1);

  if (!_cellpos_valid(T, prev) || T->cells[prev.i][prev.j].num_choices < 1) {
    return -1;
  }
  // pos is guaranteed to be in xy-system now:
  chars = (intpair) {
    frame->origin[frame->origin_range.i + x - 1],
    frame->mutant[frame->mutant_range.i + y - 1]
  };

  score = T->prob->scores->subst_scores[chars.i][chars.j];
  choice->op = (chars.i == chars.j ? 'M' : 'S');
  // No path-dependence for S/M; all previous bases are the same:
  choice->score = T->cells[prev.i][prev.j].choices[0].score + score;
  choice->base = &(T->cells[prev.i][prev.j].choices[0]);
  if (T->prob->max_new_mins >0 &&
      choice->score < T->cells[prev.i][prev.j].choices[0].cur_min) {
    choice->mins_cd = T->cells[prev.i][prev.j].choices[0].mins_cd - 1;
    choice->cur_min = choice->score;
    if (choice->mins_cd < 0) {
      return -1;
    }
  }
  return 0;
}

int _alnchoice_ID(dptable* T, int x, int y, alnchoice* choice, char op) {
  int base_idx;
  double score, max_score;
  intpair prev;

  switch (op) {
    case 'I':
      prev = _cellpos_from_xy(T->prob, x, y-1);
      break;
    case 'D':
      prev = _cellpos_from_xy(T->prob, x-1, y);
      break;
    default:
      _panick("Unknown op given to _alnchoice_ID\n");
  }
  if (!_cellpos_valid(T, prev) || T->cells[prev.i][prev.j].num_choices < 1) {
    return -1;
  }

  base_idx = 0;
  max_score = -INT_MAX;
  for (int k = 0; k < T->cells[prev.i][prev.j].num_choices; k++) {
    score = T->cells[prev.i][prev.j].choices[k].score +
            T->prob->scores->gap_extend_score;
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
  if (T->prob->max_new_mins >0 &&
      choice->score < T->cells[prev.i][prev.j].choices[0].cur_min) {
    choice->mins_cd = T->cells[prev.i][prev.j].choices[0].mins_cd - 1;
    choice->cur_min = choice->score;
    if (choice->mins_cd < 0) {
      return -1;
    }
  }
  return 0;
}

int _alnchoice_D(dptable* T, int x, int y, alnchoice* choice) {
  return _alnchoice_ID(T, x, y, choice, 'D');
}

int _alnchoice_I(dptable* T, int x, int y, alnchoice* choice) {
  return _alnchoice_ID(T, x, y, choice, 'I');
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
    row = frame->origin_range.j - frame->origin_range.i;
    col = frame->mutant_range.j - frame->mutant_range.i;
    if (T->cells[row][col].num_choices == 0) {
      return (intpair){-1, -1};
    }
  }
  else if (type == OVERLAP || type == START_ANCHORED_OVERLAP) {
    // Overlap alignments (except end-anchored ones) can end anywhere on either
    // of the bottom or right edges; find the best:
    max = -INT_MAX;
    for (i = 0; i < frame->origin_range.j - frame->origin_range.i + 1; i++){
      for (j = 0; j < frame->mutant_range.j - frame->mutant_range.i + 1; j++) {
        // Are we on the bottom row or the right column?
        if (i != frame->origin_range.j - frame->origin_range.i &&
            j != frame->mutant_range.j - frame->mutant_range.i) {
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
    for (i = 0; i < frame->origin_range.j - frame->origin_range.i + 1; i++){
      for (j = 0; j < frame->mutant_range.j - frame->mutant_range.i + 1; j++) {
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
  intpair opt = (intpair) {-1, -1};
  if (T->prob->banded_params->type == B_GLOBAL) {
    opt = _cellpos_from_xy(T->prob,
      frame->origin_range.j - frame->origin_range.i,
      frame->mutant_range.j - frame->mutant_range.i
    );
    /*printf("(d+,a)=(%d,%d)\n", opt.i, opt.j);*/ // FIXME segfaults, does it?
    if (T->cells[opt.i][opt.j].num_choices < 1) {
      return (intpair) {-1, -1};
    } else {
      return opt;
    }
  }
  if (T->prob->banded_params->type == B_OVERLAP) {
    max = -INT_MAX;
    for (i = 0; i < T->num_rows; i++){
      j = T->row_lens[i] - 1;
      if (T->cells[i][j].num_choices > 1 &&
        T->cells[i][j].choices[0].score > max) {
        max = T->cells[i][j].choices[0].score;
        opt = (intpair) {i, j};
      }
    }
    if (opt.i != -1 && opt.j != -1) {
      return opt;
    } else {
      return (intpair) {-1, -1};
    }
  }
  return (intpair) {-1, -1};
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
