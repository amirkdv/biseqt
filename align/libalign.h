/**
 * Classifies alignment problems into a number of types. Each type specifies
 * the boundary conditions (where an alignment can start and where it can end)
 * on the dynamic programming table.
 *
 * @note Any of these alignments can be made banded using `max_diversion` of
 * ::align_params.
*/
typedef enum {
  GLOBAL, /**< Solve for optimal global alignments, equivalent to the
    Needleman-Wunsch algorithm.*/
  LOCAL,  /**< Solve for optimal local alignments, equivalent to the
    Smith-Waterman algorithm.*/
  START_ANCHORED, /**< Solve for optimal local alignments starting at the
    start of alignment frame for both sequences.*/
  END_ANCHORED, /**< Solve for optimal local alignments ending at the
    end of alignment frame for both sequences.*/
  OVERLAP, /**< Solve for suffix-prefix alignments of either direction
    (including substring alignments).*/
  START_ANCHORED_OVERLAP, /**< Solve or suffix-prefix alignments starting
    at the start of alignment frame for both sequences.*/
  END_ANCHORED_OVERLAP, /**< Solve or suffix-prefix alignments ending
    at the end of alignment frame for both sequences.*/
} align_type;

/**
 * Defines an alphabet for sequences. The alphabet can contain arbitrary
 * letters; potentially letters of longer than one character. The only
 * requirement is that all letters have the same character length.
 */
typedef struct {
  int length; /**< The number of letters in the alphabet.*/
  int letter_length; /**< The common length of all letters in the alphabet.*/
  char** letters; /**< The letters of the alphabet. */
} sequence_alphabet;

/**
 * Groups together the parameters for a pairwise alignment algorithm used to
 * solve an ::align_problem.
 */
typedef struct {
  sequence_alphabet* alphabet; /**< The alphabet of the aligned sequences. */
  double **subst_scores; /**< The substitution score matrix, ordered the same
    as letters in the alphabet.*/
  double gap_open_score; /**< The gap open score in the affine gap
    penalty scheme, use 0 for a linear gap model. This is in effect regardless
    of content-dependence of gap extension score. */
  double gap_extend_score; /**< The gap extension score, or the linear gap
    score when using a linear gap model. This is only in effect when
    `content_dependent_gap_scores == NULL`.*/
  double* content_dependent_gap_scores; /**< The gap extension scores based on
    what letter is inserted/deleted. If this is `NULL`, `gap_extend_score` will
    be in effect and all gaps of equal length are scored identically.*/
  int max_diversion; /**< band length (banded alignment only if `>= 0`).*/
} align_params;

/**
 * Groups together all the information defining a pairwise alignment problem.
 * An pairwaise alignment problem is defined over a pair of sequences by
 * a "frame" on each sequence (the `(S|T)_(min|max)_idx` members), an alignment
 * type and a set of alignment parameters.
 */
typedef struct {
  int* S; /**< The "from" sequence where each letter is represented by its
    position in alphabet.*/
  int* T; /**< The "to" sequence where each letter is represented by its
    position in alphabet.*/
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
typedef struct align_choice {
  char op; /**< Either of `M` for match, `S` for substitution, `I`
    for insertion, and `D` for deletion. */
  int diversion; /**< Used to ignore out-of-band positions in the dynamic
    programming table when performing banded alignments.*/
  double score; /**< The overall score of the alignment corresponding to the
    chain of bases leading to here.*/
  struct align_choice *base; /**< The choice for the consequent subproblem on
    which this choice and its score depends.*/
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
  int num_choices; /**< The number of optimal scoring choices for the
    corresponding subproblem.*/
  struct align_choice *choices; /**< The optimal scoring choices. We need to
    keep all in case of path-dependent scoring schemes, e.g affine gap.*/
} align_dp_cell;

/**
 * Represents a solution to an alignment problem.
 */
typedef struct {
  int S_idx; /**< The starting position of the alignment along the
    "from" sequence. This is *not* relative to the start of frame.*/
  int T_idx; /**< The starting position of the alignment along the "to"
    sequence. This is *not* relative to the start of frame.*/
  double score; /**< The score of the alignment according to the
    corresponding ::align_params.*/
  char* opseq; /**< The sequence of edit operations, as in ::align_choice,
    that defines the alignment.*/
} transcript;

/**
 * A segment corresponds to a local alignment of two sequences.
 * Segments are the currency of alignment by seed extension. Each segment is
 * uniquely identified by two sequence IDs and a transcript.
 */
typedef struct segment {
  int S_id; /**< The identifier of the "from" sequence. */
  int T_id; /**<The identifier of the "to" sequence. */
  transcript* tx; /**< The transcript of the local alignment. */
} segment;

int tx_seq_len(transcript* tx, char on);
int extend_1d_once(segment* res, segment* seg,
  int* S, int* T, align_params* params,
  int window, int forward);
int extend_1d(segment* res, segment* seg,
  int* S, int* T, int S_len, int T_len, align_params* params,
  int window, int max_new_mins, int forward, int debug);
segment* extend(segment** segs, int num_segs, int* S, int* T, int S_len, int T_len,
  align_params* params, int window, int max_new_mins, double min_overlap_score, int debug);
align_dp_cell** init_dp_table(align_problem* def);
void free_dp_table(align_dp_cell** P, int row_cnt, int col_cnt);
align_dp_cell* solve(align_dp_cell** P, align_problem* def);
align_dp_cell* find_optimal(align_dp_cell** P, align_problem* def);
transcript* traceback(align_dp_cell** P, align_problem* def, align_dp_cell* end);
int* idxseq_from_charseq(sequence_alphabet* alphabet, char* sequence);
// So that we can free memory directly from python
void free(void *);
size_t strlen(const char*);
