typedef enum {
  GLOBAL, // Needleman-Wunsch
  LOCAL,  // Smith-Waterman
  START_ANCHORED, // Solves for the optimal local alignment starting at the beginning
  END_ANCHORED, // Solves for the optimal local alignment ending at the end
  OVERLAP // Solves for the optimal suffix-prefix alignment
} align_type;

typedef struct {
  int length; // number of letters in the alphabet
  int letter_length; // all letters in the alphabet must have the same length
  char** letters; // each letter can be a string
} sequence_alphabet;

typedef struct {
  sequence_alphabet* alphabet;
  double **subst_scores;
  double gap_open_score;
  double gap_extend_score;
  int max_diversion; // band length (banded alignment only if >= 0)
} align_params;

typedef struct {
  int* S; // "from" sequence, each letter translated to its position in alphabet
  int* T; // "to" sequence, each letter ...
  int S_min_idx; // opening index for the frame of S: inclusive
  int S_max_idx; // closing index for the frame of S: exclusive
  int T_min_idx; // opening ... of T: inclusive
  int T_max_idx; // closing ... of T: exclusive
  align_type type;
  align_params *params;
} align_problem;


typedef struct align_choice { // we need the name to be able to reference ourselves
  char op; // 'M'atch, 'S'ubstitute, 'I'nsert, 'D'elete, 'B'egin
  int diversion; // used to ignore out-of-band positions in banded alignment.
  double score;
  struct align_choice *base; // the choice on which the score depends on.
} align_choice;

typedef struct {
  int row; // vertical position (row number) in the DP table
  int col; // horizontal position (col number) in the DP table
  int num_choices; // number of optimal scoring choices for the subproblem
  struct align_choice *choices; // optimal scoring choices
} align_dp_cell;

align_dp_cell** define(align_problem* def);
align_dp_cell* solve(align_dp_cell** P, align_problem* def);
align_dp_cell* find_optimal(align_dp_cell** P, align_problem* def);
char* traceback(align_dp_cell** P, align_problem* def, align_dp_cell* end);
int* idxseq_from_charseq(sequence_alphabet* alphabet, char* sequence, int length);
void free_dp_table(align_dp_cell** P, int size);
