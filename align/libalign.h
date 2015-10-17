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
  char* S; // "from" sequence, order of S/T only affects I/D's in transcript
  char* T; // "to" sequence
  int S_min_idx; // opening index for the frame of S: inclusive
  int S_max_idx; // closing index for the frame of S: exclusive
  int T_min_idx; // ... of T: inclusive
  int T_max_idx; // ... of T: exclusive
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
  int idx_S; // horizontal position in the DP table
  int idx_T; // vertical position in the DP table
  int num_choices; // number of optimal scoring choices for the subproblem
  struct align_choice *choices; // optimal scoring choices
  double score;
} align_dp_cell;

align_dp_cell** define(align_problem* def);
int solve(align_dp_cell** P, align_problem* def);
char* traceback(align_dp_cell** P, align_problem* def);
int* _to_idx_sequence(sequence_alphabet* alphabet, char* sequence, int length);
