// TODO The C component could be hacked using "define" to override any of the
// functions in here. document it and expose the user's additional c|h files.
// cf. http://stackoverflow.com/a/11758777
//
// TODO _pw_internals.c requires better documentation
//
// TODO consider exposing the start and end positions of the alignment
// instead of defining new types in code. This can simply be an array of
// intpairs for either of start and end that are looped through in _alnchoice_B
// and _std_find_optimal or _banded_find_optimal.

/**
 * Any pair of integers.
 */
typedef struct {
  int i; /**< The "first" element. **/
  int j; /**< The "second" element. **/
} intpair;

/**
 * Indicates which of the pairwise alignment algorithms should be used, see
 * ::dptable.  There are two available variations of the standard dynamic
 * programming algorithm.
 *
 * First, the "standard" mode, in which time and space are quadratic in sequence
 * lengths.  Second, the "banded" mode, in which both time and space complexity
 * are linear in sequence lengths (in fact, the shorter of the sequence
 * lengths). This mode is not well-posed for local alignments.
*/
typedef enum {
  STD_MODE, /**< Standard alignment algorithm.*/
  BANDED_MODE, /**< Banded alignment algorithm.*/
} alnmode;

/**
 * Defines the boundary conditions (i.e where an alignment can start and
 * where it can end) on the dynamic programming table for standard alignments.
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
    penalty scheme, use 0 for a linear gap model.*/
  double gap_extend_score; /**< The gap extension score, or the linear gap
    score when using a linear gap model.*/
} alnscores;

/**
 * Defines the skeleton of a pairwise alignment problem: two sequences and
 * their start/end of frame positions.
 */
typedef struct {
  int* origin; /**< The "from" sequence as an array of integers.*/
  int* mutant; /**< The "to" sequence as an array of integers.*/
  intpair origin_range; /**< Vertical span of the frame in [min, max) format. */
  intpair mutant_range; /**< Horizontal span of the frame in [min, max) format. */
} alnframe;

/**
 * Groups together all parameters of a standard pairwise alignment problem:
 * for now only an alignment type.
 */
typedef struct {
  std_alntype type; /**< The subtype of standard algorithms.*/
} std_alnparams;

/**
 * Groups together all parameters of a banded pairwise alignment problem:
 * alignment type, scores, and bounding diagonals.
 */
typedef struct {
  banded_alntype type; /**< The subtype of banded algorithms.*/
  int dmin; /**< The upper diagonal bounding the band. */
  int dmax; /**< The lower diagonal bounding the band. */
} banded_alnparams;

/**
 * The definition of a pairwise alignment problem the contents of which
 * is enough to define and solve an alignment.
 */
typedef struct {
  alnframe *frame; /**< The skeleton of the DP table. */
  alnscores *scores; /**< The parameters defining an optimal solution.*/
  int max_new_mins; /**< Maximum number of times a new minimum score is
    tolerated before alignment is aborted (only operative if positive). Not
    to be confused with ::seedext_params.max_new_mins. **/
  alnmode mode; /**< Indicates which of standard or banded alignment this is. */
  union {
    std_alnparams* std_params; /**< The parameters for standard algorithms. */
    banded_alnparams* banded_params; /**< The parameters for banded algorithms. */
  };
} alnprob;

/**
 * The fundamental unit of solving pairwise alignment problems. Each
 * cell of the dynamic programming table corrsponds to a subproblem of the
 * overall alignment problem. For each subproblem, a "choice" corresponds to a
 * single edit operation (which applies to the very last positions of the
 * subproblem) and a "base" which is a choice residing at the cell to which the
 * op points (the "base" is `NULL` for the choice of starting an alignment at a
 * given position).
 *
 * Each choice can thus be traced all the way back via the chain of "base"
 * pointers to the beginning of the alignment to which it belongs. It is the
 * total score of this alignment that is stored at each cell for each choice.
 *
 * This architecture allows us to easily capture path-dependent scoring schemes,
 * for example, the affine gap penalty.
 */
typedef struct alnchoice {
  char op; /**< Either of `M` for match, `origin` for substitution, `I`
    for insertion, and `D` for deletion. */
  double score; /**< The overall score of the alignment corresponding to the
    chain of bases leading to here.*/
  struct alnchoice *base; /**< The choice for the consequent subproblem on
    which this choice and its score depends.*/
  int mins_cd; /**< The countdown of the number of new score minima encountered
    along the path leading to this choice. */
  int cur_min; /**< The current minimum score along the path leading to this
    choice. */
} alnchoice;

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
  int num_rows; /**< The height of the DP table in memory. */
  int* row_lens; /**< Row lengths (width) in the DP table at each row. */
  alnprob* prob; /**< The alignment problem.*/
} dptable;

/**
 * Represents a solution to an alignment problem.
 */
typedef struct alignment {
  int origin_idx; /**< The starting position of the alignment along the
    original sequence. This is *not* relative to the start of alignment frame.*/
  int mutant_idx; /**< The starting position of the alignment along the "to"
    sequence. This is *not* relative to the start of alignment frame.*/
  double score; /**< The score of the alignment according to the
    corresponding ::alnscores.*/
  char* transcript; /**< The sequence of edit operations, as in ::alnchoice
    that defines the alignment.*/
} alignment;

/**
 * Groups together seed extension parameters.
 */
typedef struct {
  alnscores* scores; /**< Substitution and gap scores. **/
  int window; /**< The width of the rolling window of alignment. **/
  double min_score; /**< The minimum required score for an overlap to be reported. **/
  int max_new_mins; /**< Maximum number of times a new minimum score is
    tolerated before alignment is aborted. This differs from
    ::alnprob.max_new_mins in that this parameter applies to the number of
    times the score over an entire rolling frame reaches new minima and not
    single edit operations.
    **/
} seedext_params;

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
int dptable_init(dptable* T);

/**
 * Frees the allocated memory for the cells of a given ::dptable.
 *
 * @param T the dynamic programming table containing all the info.
 */
void dptable_free(dptable* T);

/**
 * Traces back *an* alignment (and not all alignments with identical scores)
 * from a given cell all the way back to a cell with a `NULL` base (i.e an
 * alignment start cell).
 *
 * @param T The *solved* dynamic programming table.
 * @param end The desired ending point of alignment which becomes the starting
 *    point of traceback.
 *
 * @note
 *    Finding more than one optimal alignment is a nontrivial search problem,
 *    and requires some sort of global state keeping to avoid convoluted
 *    recursions. I couldn't get it right in the first go; leave for later.
 */
alignment* dptable_traceback(dptable* T, intpair end);

/**
 * Given an ::alignment with allocated DP table, solves the alignment problem (i.e
 * populates the alignment table) and returns the optimal ending point for the
 * alignment.  The optimal score and edit transcript can then be obtained by
 * using `dptable_traceback`. Half of the constraints imposed by an alignment
 * type are implemented here where we decide what positions on the table can be
 * starting positions of alignments. The other half of the constraints
 * concerning the ending position of the alignment on the table are enforced
 * in `dptable_traceback`.
 *
 * @param T The dynamic programming table.
 * @return The optimal cell for the alignment to end at or {-1,-1} if error.
 */
intpair dptable_solve(dptable* T);

/**
 * Given an array of partial alignments tries to extend all in both directions
 * and returns a fully extended alignment as soon as it finds one.
 *
 * @param alns The original partial alignments.
 * @param num_alns The number of provided partial alignments.
 * @param frame The frame defining the alignment of the original sequences.
 * @param params Seed extension parameters.
 *
 * @return 0 if the seed successfully extends to the boundary of either of
 *   the sequences and -1 otherwise.
 */
alignment* extend(alignment** alns, int num_alns, alnframe* frame, seedext_params* params);

/**
 * Given a partial alignment fully extends it in one direction. A fully
 * extended alignment (in one direction) is one that hits the boundary
 * (beginning or end) of either of the sequences.
 *
 * @param res Extended alignment to be populated here.
 * @param aln The original partial alignment.
 * @param frame The frame defining the alignment of the original sequences.
 * @param params Seed extension parameters.
 * @param forward Either of 0 or 1 indicating the direction of extension.
 *
 * @return 0 if the seed successfully extends to the boundary of either of
 *   the sequences and -1 otherwise.
 */
int extend_1d(alignment* res, alignment* aln, alnframe* frame, seedext_params* params, int forward);

/**
 * Given an alignment, extends it in the given direction by one window.
 *
 * @param aln The alignment to be extended.
 * @param frame The frame defining the alignment of the original sequences.
 * @param params Seed extension parameters.
 * @param forward Either of 0 or 1 indicating the direction of extension.
 *
 * @return -1 if an error occurs and 0 otherwise.
 */
int extend_1d_once(alignment* aln, alnframe* frame, seedext_params* params, int forward);

// --------------------- Internals ---------------------
int _alnchoice_B(dptable *T, int x, int y, alnchoice* choice);
int _alnchoice_M(dptable *T, int x, int y, alnchoice* choice);
int _alnchoice_I(dptable* T, int x, int y, alnchoice* choice);
int _alnchoice_D(dptable* T, int x, int y, alnchoice* choice);

intpair _std_find_optimal(dptable* T);
intpair _banded_find_optimal(dptable* T);

int _std_table_init_dims(dptable* T);
int _banded_table_init_dims(dptable* T);
int _table_init_cells(dptable* T);

intpair _xlim(alnprob* prob);
intpair _ylim(alnprob* prob, int x);

intpair _cellpos_from_xy(alnprob* prob, int x, int y);
intpair _xy_from_cellpos(alnprob* prob, int i, int j);

void _print_mem_usage();
void _panick(char* message);
