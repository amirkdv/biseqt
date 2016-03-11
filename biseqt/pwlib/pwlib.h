/**
 * Indicates which of the pairwise alignment algorithms should be used, cf.
 * ::dptable.  There are two available variations of the standard dynamic
 * programming algorithm.
 *
 * First, the "standard" mode, in which time and space are quadratic in sequence
 * lengths.  Second, the "banded" mode, in which both time and space complexity
 * are linear in sequence lengths (in fact, the shorter of the sequence
 * lengths). This mode is not well-posed for local alignments.
 *
 * @note Any of these alignments can be made banded using `bradius` of
 * ::alnscores. This, however, only reduces the time complexity but not the
 * memory complexity.
*/
typedef enum {
  STD_MODE, /**< Standard alignment algorithm.*/
  BANDED_MODE, /**< Banded alignment algorithm.*/
} alnmode;

/**
 * Defines the boundary conditions (i.e where an alignment can start and
 * where it can end) on the dynamic programming table for standard alignments.
 *
 * @note Any of these alignments can be made banded using `bradius` of
 * ::alnscores. This, however, only reduces the time complexity but not the
 * memory complexity.
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
} std_alntype;

/**
 * Defines the boundary conditions (i.e where an alignment can start and
 * where it can end) on the dynamic programming table for banded alignments.
*/
typedef enum {
  B_GLOBAL, /**< Solve for banded global alignments.*/
  B_OVERLAP, /**< Solve for banded suffix-prefix alignments of either direction
    (including substring alignments, if within band).*/
} banded_alntype;

/**
 * Groups together the parameters for a pairwise alignment algorithm used to
 * solve an alignment problem.
 */
typedef struct {
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
} alnscores;

/**
 * Defines the skeleton of a pairwise alignment problem: two sequences and
 * their start/end of frame positions.
 */
typedef struct {
  int* S; /**< The "from" sequence as an array of integers.*/
  int* T; /**< The "to" sequence as an array of integers.*/
  int S_min_idx; /**< The opening index for the frame of S, inclusive.*/
  int S_max_idx; /**< The closing index for the frame of S, exclusive.*/
  int T_min_idx; /**< The opening index for the frame of T, inclusive.*/
  int T_max_idx; /**< The closing index for the frame of T, exclusive.*/
} alnframe;

/**
 * Groups together all information defining a standard pairwise alignment
 * problem: a frame, alignment type, alignment parameters, and a band radius
 * (optional).
 */
typedef struct {
  alnframe *frame; /**< The skeleton of the DP table. */
  alnscores *scores; /**< The parameters defining an optimal solution.*/
  std_alntype type; /**< The subtype of standard algorithms.*/
} std_alnprob;

/**
 * Groups together all information defining a banded pairwise alignment
 * problem: a frame, alignment type, alignment parameters, starting diagonal,
 * and a band radius.
 */
typedef struct {
  alnframe* frame; /**< The skeleton of the DP table. */
  banded_alntype type; /**< The subtype of standard algorithms.*/
  alnscores* scores; /**< The scores for various alignment moves..*/
  int bradius; /**< The band radius, inclusive. */
  int ctrdiag; /**< The diagonal at the center of band (must be between `-|T|`
    and `+|S|`. */
} banded_alnprob;

/**
 * The fundamental unit of solving standard pairwise alignment problems. Each
 * cell of the dynamic programming table corrsponds to a subproblem of the
 * overall alignment problem. For each subproblem, a "choice" corresponds to a
 * single edit operation (which applies to the very last positions of the
 * subproblem) and a "base" which is a choice for the cell to which the op
 * points (the "base" is `NULL` for the choice of starting an alignment at a
 * given position).
 *
 * Each choice can thus be traced all the way back via the chain of "base"
 * pointers to the beginning of the alignment it corresponds to. It is the total
 * score of this alignment that is stored at each cell for each choice.
 *
 * This architecture allows us to easily capture path-dependent scoring schemes,
 * for example, the affine gap penalty.
 */
typedef struct alnchoice {
  char op; /**< Either of `M` for match, `S` for substitution, `I`
    for insertion, and `D` for deletion. */
  double score; /**< The overall score of the alignment corresponding to the
    chain of bases leading to here.*/
  struct alnchoice *base; /**< The choice for the consequent subproblem on
    which this choice and its score depends.*/
} alnchoice;

/**
 * Any pair of integers.
 */
typedef struct {
  int i;
  int j;
} gridcoord;

/**
 * A single cell of the dynamic programming table. Each cell has an array of
 * optimal choices linking to choices in dependent cells.
 */
typedef struct {
  int num_choices; /**< The number of optimal scoring choices for the
    corresponding subproblem (i.e number of elements in choices).*/
  struct alnchoice *choices; /**< The optimal scoring choices. We need to
    keep all in case of path-dependent scoring schemes, e.g affine gap.*/
} dpcell;

/**
 * The full dynamic programming table for a standard or banded alignment.
 */
typedef struct {
  dpcell** cells; /**< The DP table in (i,j) or (d,a) coordinates. */
  int num_rows; /**< The number of rows in P. */
  int num_cols; /**< The number of rows in P. */
  alnmode mode; /**< Indicates which of standard or overlap alignment this is. */
  union {
    std_alnprob* std_prob; /**< The alignment problem for standard algorithms. */
    banded_alnprob* banded_prob; /**< The alignment problem for banded algorithms. */
  };
} dptable;

/**
 * Represents a solution to an alignment problem.
 */
typedef struct {
  int S_idx; /**< The starting position of the alignment along the
    "from" sequence. This is *not* relative to the start of alignment frame.*/
  int T_idx; /**< The starting position of the alignment along the "to"
    sequence. This is *not* relative to the start of alignment frame.*/
  double score; /**< The score of the alignment according to the
    corresponding ::alnscores.*/
  char* opseq; /**< The sequence of edit operations, as in ::alnchoice
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

int dptable_init(dptable* T);
void dptable_free(dptable* T);
gridcoord dptable_solve(dptable* T);

segment* extend(segment** segs, int num_segs, int* S, int* T, int S_len, int T_len,
  alnscores* scores, int window, int max_new_mins, double min_overlap_score, int debug);

gridcoord stdpw_solve(dptable* T);
gridcoord stdpw_find_optimal(dptable* T);
transcript* stdpw_traceback(dptable* T, gridcoord end);

gridcoord bandedpw_solve(dptable* T);
gridcoord bandedpw_find_optimal(dptable* T);
transcript* bandedpw_traceback(dptable* T, dpcell* end);

// Internals
int _alnchoice_B(dptable *T, int i, int j, alnchoice* choice);
int _alnchoice_MS(dptable *T, int i, int j, alnchoice* choice);
int _alnchoice_ID(dptable* T, int i, int j, alnchoice* choice, char op);
