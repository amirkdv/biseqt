#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pwlib.h"


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
  int S_len, T_len;
  alnframe* frame;
  switch (T->mode) {
    case STD_MODE:
      frame = T->std_prob->frame;
      S_len = frame->S_max_idx - frame->S_min_idx;
      T_len = frame->T_max_idx - frame->T_min_idx;
      T->num_rows = S_len + 1;
      T->num_cols = T_len + 1;
      break;
    case BANDED_MODE:
      frame = T->banded_prob->frame;
      S_len = frame->S_max_idx - frame->S_min_idx;
      T_len = frame->T_max_idx - frame->T_min_idx;
      T->num_rows = 1 + 2 * T->banded_prob->bradius;
      T->num_cols = 1 + (S_len < T_len ? S_len : T_len);
      break;
  }
  T->cells = malloc(T->num_rows * sizeof(dpcell *));
  if (T->cells == NULL) {
    printf("Failed to allocate memory (`bandedpw_init()`).\n");
    return -1;
  }
  for (int i = 0; i < T->num_rows; i++) {
    T->cells[i] = malloc(T->num_cols * sizeof(dpcell));
    if (T->cells[i] == NULL) {
      printf("Failed to allocate memory (`bandedpw_init()`).\n");
      return -1;
    }
    for (int j = 0; j < T->num_cols; j++) {
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
  for (i = 0; i < T->num_rows; i++) {
    for (j = 0; j < T->num_cols; j++) {
      if (T->cells[i][j].num_choices > 0) {
        free(T->cells[i][j].choices);
      }
    }
    free(T->cells[i]);
  }
  free(T->cells);
}

/**
 * Given a ::dptable dispatches to the write solver to find the optimal ending
 * position of the alignment.
 *
 * @param T
 *    Initialized alignment table to be solved.
 *
 * @return The optimal cell for the alignment to end at or {-1,-1} if error.
 */
gridcoord dptable_solve(dptable* T) {
  switch (T->mode) {
    case STD_MODE:
      return stdpw_solve(T);
    case BANDED_MODE:
      return (gridcoord) {-1, -1};
      /*return bandedpw_solve(T);*/
  }
  return (gridcoord) {-1, -1};
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
