/**
 * Classifies alignment problems into a number of types. Note that Any of these
 * problems can be solved in a banded (and thus naturally non-optimal) fashion,
 * see ::align_params.
*/
typedef enum {
  GLOBAL, /**< Equivalent to the [Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) algorithm.*/
  LOCAL,  /**< Equivalent to the [Smith-Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) algorithm.*/
  START_ANCHORED, /**< Local alignments starting at the beginning of both sequences.*/
  END_ANCHORED, /**< Local alignments ending at the end of both sequences.*/
  OVERLAP, /**< Suffix-prefix alignments of either direction (including substring alignments).*/
  START_ANCHORED_OVERLAP, /**< Suffix-prefix alignments starting at the beginning of both sequences.*/
  END_ANCHORED_OVERLAP, /**< Suffix-prefix alignments ending at the end of both sequences.*/
} align_type;

/**
 * Defines an alphabet for sequences. The alphabet can contain arbitrary
 * letters; potentially letters of longer than one character. The only
 * requirement is that all letters have the same character length.
 */
typedef struct {
  int length; /**< The number of letters in the alphabet.*/
  int letter_length; /**< All letters in the alphabet must have the same length.*/
  char** letters; /**< The letters of the alphabet. */
} sequence_alphabet;

/**
 * Groups together the parameters for a pairwise alignment algorithm used to
 * solve an ::align_problem.
 */
typedef struct {
  sequence_alphabet* alphabet; /**< The alphabet of the aligned sequences. */
  double **subst_scores;
  double gap_open_score;
  double gap_extend_score;
  int max_diversion; // band length (banded alignment only if >= 0)
} align_params;

/**
 * Groups together all the information defining a pairwise alignment problem.
 * An pairwaise alignment problem is defined over a pair of sequences by
 * a "frame" on each sequence (the `(S|T)_(min|max)_idx` members), an alignment
 * type and a set of alignment parameters.
 */
typedef struct {
  int* S; /**< The "from" sequence, each letter translated to its position in alphabet.*/
  int* T; /**< The "to" sequence, each letter translated to its position in alphabet.*/
  int S_min_idx; /**< The opening index for the frame of S, inclusive.*/
  int S_max_idx; /**< The closing index for the frame of S: exclusive.*/
  int T_min_idx; /**< The opening index for the frame of T: inclusive.*/
  int T_max_idx; /**< The closing index for the frame of T: exclusive.*/
  align_type type; /**< The type of alignment problem.*/
  align_params *params; /**< The parameters defining an optimal solution.*/
} align_problem;

/**
 * The fundamental unit of solving pairwise alignment problems. Each cell of the
 * dynamic programming table corrsponds to a subproblem of the overall
 * alignment problem. For each subproblem, a "choice" corresponds to a single
 * edit operation (which applies to the very last positions of the subproblem)
 * and a "base" which is a choice for the cell to which the op points (the
 * "base" is `NULL` for the choice of starting an alignment at a given position).
 *
 * Each choice can thus be traced all the way back via the chain of "base"
 * pointers to the beginning of the alignment it corresponds to. It is the
 * total score of this alignment that is stored at each cell for each choice.
 *
 * This architecture allows us to easily capture path-dependent scoring schemes,
 * for example, the affine gap penalty.
 */
typedef struct align_choice { // we need the name to be able to reference ourselves
  char op; /**< Either of `M` for match, `S` for substitution, `I` for insertion, and `D` for deletion. */
  int diversion; /**< Used to ignore out-of-band positions in banded alignment.*/
  double score; /**< The overall score of the alignment corresponding to the chain of bases leading to here.*/
  struct align_choice *base; /**< The choice for the consequent subproblem on which this choice and its score depends.*/
} align_choice;

/**
 * Each cell of the dynamic programming table.
 * @note The `row` and `col` of each cell corresponds to the position in the
 * DP table. Therefore, these are always one larger than the 0-starting
 * position of the subproblem in the "from" and "to" sequences.
 */
typedef struct {
  int row; /**< Vertical position (row number) in the table.*/
  int col; /**< horizontal position (col number) in the table.*/
  int num_choices; /**< The number of optimal scoring choices for the subproblem.*/
  struct align_choice *choices; /**< The optimal scoring choices. We need to keep all in case of path-dependent scoring schemes, e.g affine gap.*/
} align_dp_cell;

/**
 * Represents a solution to an alignment problem.
 */
typedef struct {
  int S_idx; /**< The starting position of the alignment along the "from" sequence.*/
  int T_idx; /**< The starting position of the alignment along the "to" sequence.*/
  double score; /**< The score of the alignment according to the corresponding ::align_params.*/
  char* opseq; /**< The sequence of edit operations, as in ::align_choice, that defines the alignment.*/
} transcript;

align_dp_cell** init_dp_table(align_problem* def);
void free_dp_table(align_dp_cell** P, int row_cnt, int col_cnt);
align_dp_cell* solve(align_dp_cell** P, align_problem* def);
align_dp_cell* find_optimal(align_dp_cell** P, align_problem* def);
transcript* traceback(align_dp_cell** P, align_problem* def, align_dp_cell* end);
int* idxseq_from_charseq(sequence_alphabet* alphabet, char* sequence, int length);
// So that we can free memory directly from python
void free(void *);
size_t strlen(const char*);
