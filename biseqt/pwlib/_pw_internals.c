#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>

#include "pwlib.h"

/**
 * If possible populates an alignment 'B' move for a given position of the table
 * (in either coordinate systems depending on alignment mode).
 *
 * @return 0 if choice successfully built, -1 if choice not possible.
 */
int _alnchoice_B(dptable *T, int i, int j, alnchoice* choice) {
  // return -1 if the current cell is not elligible for a B:
  switch (T->mode) {
    case STD_MODE:
      if (i == 0 && j == 0) {
        break;
      }
      if (T->std_prob->type == LOCAL || T->std_prob->type == END_ANCHORED) {
        break;
      }
      if (T->std_prob->type == OVERLAP && (i == 0 || j == 0)) {
        break;
      }
      if (T->std_prob->type == END_ANCHORED_OVERLAP && (i == 0 || j == 0)) {
        break;
      }
      // not elligible:
      return -1;
    case BANDED_MODE:
      //FIXME
      break;
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
int _alnchoice_MS(dptable *T, int i, int j, alnchoice* choice) {
  int s,t, iprev, jprev;
  double score;
  switch (T->mode) {
    case STD_MODE:
      if (i <=0 || j <= 0) {
        return -1;
      }
      iprev = i - 1;
      jprev = j - 1;
      s = T->std_prob->frame->S[T->std_prob->frame->S_min_idx + i - 1];
      t = T->std_prob->frame->T[T->std_prob->frame->T_min_idx + j - 1];
      score = T->std_prob->scores->subst_scores[s][t];
      break;
    case BANDED_MODE:
      //FIXME
      break;
  }
  if (T->cells[iprev][jprev].num_choices <= 0) {
    return -1;
  }
  choice->op = (s==t ? 'M' : 'S');
  // All choices to (i-1,j-1) have the same score:
  choice->score = T->cells[iprev][jprev].choices[0].score + score;
  choice->base = &(T->cells[iprev][jprev].choices[0]);
  return 0;
}

/**
 * If possible populates an alignment 'M' or 'S' move for a given position of
 * the table (in either coordinate systems depending on alignment mode).
 *
 * @return 0 if choice successfully built, -1 if choice not possible.
 */
int _alnchoice_ID(dptable* T, int i, int j, alnchoice* choice, char op) {
  int iprev, jprev, base_idx, content;
  double score, max_score, gap_open_score, gap_extend_score;
  switch (T->mode) {
    case STD_MODE:
      if (op == 'D') {
        if (i <= 0) {
          return -1;
        }
        iprev = i-1;
        jprev = j;
        content = T->std_prob->frame->S[T->std_prob->frame->S_min_idx + i - 1];
      } else if (op == 'I') {
        if (j <= 0) {
          return -1;
        }
        iprev = i;
        jprev = j-1;
        content = T->std_prob->frame->T[T->std_prob->frame->T_min_idx + j - 1];
      } else {
        return -1;
      }
      gap_open_score = T->std_prob->scores->gap_open_score;
      if (T->std_prob->scores->content_dependent_gap_scores == NULL) {
        gap_extend_score = T->std_prob->scores->gap_extend_score;
      } else {
        gap_extend_score = T->std_prob->scores->content_dependent_gap_scores[content];
      }
      break;
    case BANDED_MODE:
      //FIXME
      break;
  }
  if (T->cells[iprev][jprev].num_choices <= 0) {
    return -1;
  }
  base_idx = 0;
  max_score = -INT_MAX;
  for (int k = 0; k < T->cells[iprev][jprev].num_choices; k++) {
    score = T->cells[iprev][jprev].choices[k].score + gap_extend_score;
    score += T->cells[iprev][jprev].choices[k].op == op ? 0 : gap_open_score;
    if (score > max_score) {
      max_score = score;
      base_idx = k;
    }
  }
  choice->op = op;
  choice->score = max_score;
  choice->base = &(T->cells[iprev][jprev].choices[base_idx]);
  return 0;
}

