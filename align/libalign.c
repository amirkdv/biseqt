#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdio.h>
#include <string.h>

#include "libalign.h"

/**
 * Given an alignment problem definition, creates and initializes the DP matrix.
 *
 * @param def (align_problem*): pointer to problem definition struct. The only
 *    members that are relevant at this point is the corresponding frames of S
 *    and T (which determine the size of the table).
 *
 * @return align_dp_cell**: pointer to the malloced DP matrix.
 */
align_dp_cell** define(align_problem* def) {
  int n = def->S_max_idx - def->S_min_idx;
  int m = def->T_max_idx - def->T_min_idx;
  align_dp_cell** P = malloc((n+1) * sizeof(align_dp_cell *));
  if (P == NULL) {
    printf("Failed to allocate memory.\n");
    return NULL;
  }
  // We need an additional row/col in the beginning. Table indices are therefore
  // exactly one ahead of subproblem indices.
  for (int i = 0; i < n+1; i++) {
    P[i] = malloc((m+1) * sizeof(align_dp_cell));
    if (P[i] == NULL) {
      printf("Failed to allocate memory.\n");
      return NULL;
    }
    for (int j = 0; j < m+1; j++) {
      P[i][j] = (align_dp_cell) {i, j, 0, NULL};
    }
  }
  return P;
}


/**
 * Frees the allocated memory for a given alignment problem so that we can reuse
 * the same align_problem* over and over.
 */
void free_dp_table(align_dp_cell** P, int size) {
  if (P == NULL) {
    return;
  }
  for (int i = 0; i < size; i++) {
    free(P[i]);
  }
  free(P);
}

/**
 * Given an alignment with allocated DP table, solves the alignment problem (i.e
 * populates the alignment table). The optimal score and transcript can then
 * be obtained by calling traceback on P.
 */
align_dp_cell* solve(align_dp_cell** P, align_problem* def) {
  int n = def->S_max_idx - def->S_min_idx;
  int m = def->T_max_idx - def->T_min_idx;
  double max_score, prev_score, max_prev_score;
  int num_choices, max_prev_choice_idx, num_max_scores;
  int i,j,k;
  int *max_score_alts;
  align_choice *alts;
  int s,t;

  // Base case
  P[0][0].num_choices = 1;
  P[0][0].choices = malloc(sizeof(align_choice));
  if (P[0][0].choices == NULL) {
    printf("Failed to allocate memory.\n");
    return NULL;
  }
  P[0][0].choices[0].op = 'B';
  P[0][0].choices[0].diversion = 0;
  P[0][0].choices[0].score = 0;
  P[0][0].choices[0].base = NULL;

  // Populate the table
  int infocounter = 0;
  for (i = 0; i < n+1; i++) {
    for (j = 0; j < m+1; j++) {
      if (i == 0 && j == 0) {
        // Already dealt with
        continue;
      }
      // Sane default so we can continue when a cell is not worth pursuing
      P[i][j].num_choices = 0;

      // Allocate for all 4 possible choices (B,M/S,I,D)
      alts = malloc(4 * sizeof(align_choice));
      if (alts == NULL) {
        printf("Failed to allocate memory.\n");
        return NULL;
      }
      // Build all the alternatives at cell (i,j)
      num_choices = 0;
      if (def->type == LOCAL ||
          (def->type == OVERLAP && j == 0) ||
          (def->type == OVERLAP && i == 0) ||
          (def->type == END_ANCHORED)
        ) {
        // 1. local alignments can start anywhere,
        // 2. overlap alignments can start anywhere on the left column or the
        //    top row,
        // 3. end-anchored alignments can start anywhere.
        alts[num_choices].op = 'B';
        alts[num_choices].score = 0;
        alts[num_choices].diversion = 0;
        alts[num_choices].base = NULL;

        num_choices++;
      }
      // To (i-1,j)
      if (i > 0) {
        // Are there choices for (i-1,j) and are we inside the band?
        if (P[i-1][j].num_choices > 0 && (
              def->params->max_diversion < 0 ||
              // diversion of all choices in the same cell are identical:
              abs(P[i-1][j].choices[0].diversion - 1) <= def->params->max_diversion
            )) {
          max_prev_choice_idx = 0;
          max_prev_score = -INT_MAX;
          for (k = 0; k < P[i-1][j].num_choices; k++) {
            prev_score = P[i-1][j].choices[k].score
              + def->params->gap_extend_score;
            if (P[i-1][j].choices[k].op != 'D') {
              prev_score += def->params->gap_open_score;
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
              def->params->max_diversion < 0 ||
              // diversion of all choices in the same cell are identical:
              abs(P[i][j-1].choices[0].diversion + 1) <= def->params->max_diversion
            )) {
          max_prev_choice_idx = 0;
          max_prev_score = - INT_MAX;
          for (k = 0; k < P[i][j-1].num_choices; k++) {
            prev_score = P[i][j-1].choices[k].score
              + def->params->gap_extend_score;
            if (P[i][j-1].choices[k].op != 'I') {
              prev_score += def->params->gap_open_score;
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
          // the indices in the table are on ahead of the indices of letters:
          s = def->S[def->S_min_idx + i - 1];
          t = def->T[def->T_min_idx + j - 1];
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
            + def->params->subst_scores[s][t];

          num_choices++;
        }
      }

      // ================== Find the best alternatives ==================
      if (num_choices == 0) {
        continue;
      }
      infocounter ++;
      if (infocounter % 1000 == 0) {
        /*print_mem_usage();*/
      }

      // indices of maximum choices in the `alts' array
      max_score_alts = malloc(num_choices * sizeof(int));
      if (max_score_alts == NULL) {
        printf("Failed to allocate memory.\n");
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
      P[i][j].choices = malloc(num_max_scores * sizeof(align_choice));
      if (P[i][j].choices == NULL) {
        printf("Failed to allocate memory.\n");
        return NULL;
      }
      for (k = 0; k < num_max_scores; k++) {
        P[i][j].choices[k] = alts[max_score_alts[k]];
      }
      free(max_score_alts);
    }
  }
  free(alts);
  return find_optimal(P, def);
}

/**
 * Finds the optimal cell to start traceback from given a populated DP table
 * and an alignment problem definition (to know where to look for the optimal
 * cell).
 * @param P (align_dp_cell**):  pointer to a solved alignment DP table.
 * @param def (align_problem*): pointer to the align_problem struct.
 *
 * @return align_dp_cell*: pointer to the optimal cell of the table.
 */
align_dp_cell* find_optimal(align_dp_cell** P, align_problem* def) {
  double max;
  int i,j;
  int row = -1, col = -1;
  if (def->type == GLOBAL || def->type == END_ANCHORED) {
    // Global and end-anchored alignments must end at the bottom right corner
    row = def->S_max_idx - def->S_min_idx;
    col = def->T_max_idx - def->T_min_idx;
    if (P[row][col].num_choices == 0) {
      return NULL;
    }
  }
  else if (def->type == OVERLAP) {
    // Overlap alignments must end on the bottom row or the right column,
    // find the best:
    max = -INT_MAX;
    for (i = 0; i < def->S_max_idx - def->S_min_idx + 1; i++){
      for (j = 0; j < def->T_max_idx - def->T_min_idx + 1; j++) {
        // Are we on the bottom row or the right column?
        if (i != def->S_max_idx - def->S_min_idx &&
            j != def->T_max_idx - def->T_min_idx) {
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
  else if (def->type == LOCAL || def->type == START_ANCHORED) {
    // Local and start-anchored alignments can end anywhere, find the best:
    max = P[0][0].choices[0].score;
    for (i = 0; i < def->S_max_idx - def->S_min_idx + 1; i++){
      for (j = 0; j < def->T_max_idx - def->T_min_idx + 1; j++) {
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
 * Traces back the calculated optimal solutions backwards starting from a given
 * "end" point:
 * @param P (align_dp_cell**):  pointer to a solved alignment DP table.
 * @param def (align_problem*): pointer to the align_problem struct.
 * @param end (align_dp_cell*): ending point of alignment, starting point of
 *    traceback
 *
 * @return char*: string with the following format:
 *
 *           (<Si,Tj>),<score>:<opseq>
 *
 *     Si and Tj are integers specifying the positiong along each string where
 *     the alignment begins. Score is the score of the transcript to 2 decimal
 *     places. What follows the ':' is a sequence of "ops" defined as follows:
 *       B begin
 *       M match
 *       S substitution
 *       I insert
 *       D delete
 *     All op sequences begin with a B and insertion/deletions are meant to
 *     mean "from S to T".
 *
 * Limitations:
 *  - Finding more than one optimal alignment is a nontrivial search problem,
 *    and requires some sort of global state keeping to avoid convoluted
 *    recursions. I couldn't get it right in the first go; leave for later.
 */
char *traceback(align_dp_cell** P, align_problem* def, align_dp_cell* end) {
  char *infostr, *transcript, *ret;
  char op;
  int idx_S = end->row,
      idx_T = end->col;
  align_dp_cell curr = P[idx_S][idx_T];
  // string manipulation indices for the transcript:
  int len, infolen, pos;
  // allocate more than enough memory for the transcript as we trace back.
  len = idx_S + idx_T + 1;
  // We write ops to rev_transcript backwards starting from the end (position `len')
  char rev_transcript[len];
  pos = len - 1;
  while (1) {
    op = curr.choices[0].op;
    if (op == 'B') {
      break;
    }
    pos--;
    rev_transcript[pos] = op;
    if (op == 'M' || op == 'S') {
      idx_S --;
      idx_T --;
    }
    if (op == 'I') {
      idx_T --;
    }
    if (op == 'D') {
      idx_S --;
    }
    curr = P[idx_S][idx_T];
  }
  if (pos == len - 1) {
    // empty opseq
    return NULL;
  }
  // build the info string: "(<idx_S>,<idx_T>),<score>:"
  infolen = 8; // for "(,),:" and the 2 decimal points
  if (idx_S + def->S_min_idx > 0 ) {
    infolen += (int)floor(log10(idx_S + def->S_min_idx)) + 1;
  }
  if (idx_T + def->T_min_idx > 0 ) {
    infolen += (int)floor(log10(idx_T + def->T_min_idx)) + 1;
  }
  if (end->choices[0].score != 0) {
    infolen += (int)floor(log10(fabs(end->choices[0].score))) + 3;
    if (end->choices[0].score < 0) {
      infolen += 1;
    }
  }
  infostr = malloc(infolen);
  sprintf(infostr, "(%d,%d),%.2f:", (idx_S + def->S_min_idx), (idx_T + def->T_min_idx), end->choices[0].score);
  // the backtraced transcript was written backwords to the end of rev_transcript
  len = len - pos - 1;
  transcript = malloc(len + 1);
  strncpy(transcript, rev_transcript + pos, len);
  transcript[len] = '\0';
  ret = malloc(infolen + len + 1);
  sprintf(ret, "%s%s", infostr, transcript);
  return ret;
}

/**
 * Translates a given sequence to integer indices of each letter
 * in the corresponding alphabet. This is done at the first step of solve()
 * so that we don't have to worry about alphabet pecularities (e.g letters
 * longer than one character each).
 */
int* idxseq_from_charseq(sequence_alphabet* alphabet, char* sequence, int length) {
  int i,j;
  int cur_idx;
  int* idx_seq = malloc(length * sizeof(int));
  char* cur = malloc(alphabet->letter_length);
  for (i = 0; i < length; i++) {
    cur_idx = -1;
    snprintf(cur,
      alphabet->letter_length + 1,
      "%s", sequence + i*(alphabet->letter_length)
    );
    for (j = 0; j < alphabet->length; j++) {
      if (strcmp(alphabet->letters[j], cur) == 0) {
        cur_idx = j;
        break;
      }
    }
    if (cur_idx == -1) {
      printf("Invalid letter: %s\n", cur);
    }
    idx_seq[i] = cur_idx;
  }
  free(cur);
  return idx_seq;
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

