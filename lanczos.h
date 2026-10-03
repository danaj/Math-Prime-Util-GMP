#ifndef MPU_LANCZOS_H
#define MPU_LANCZOS_H

#include "ptypes.h"

/* The panel dense solver wins through this measured reduced-column
 * crossover.  Both solvers remain valid on either side of it. */
#define LA_DENSE_CROSSOVER_COLS 1536UL

typedef struct {
  /* 32-bit indices cover SIQS through its 431-bit developer hard cap. */
  uint32_t *data;
  uint32_t weight;
  uint32_t orig;
} la_col_t;

/* Return the bit showing whether a column belongs to a dependency. */
static INLINE uint64_t la_get_null_entry(const uint64_t *nullrows,
                                         unsigned long column,
                                         unsigned long dependency) {
  return nullrows[column] & ((uint64_t)1 << dependency);
}

/* Peel singleton rows and trim excess columns before solving. */
extern void la_reduce_matrix(unsigned long *nrows,
                             unsigned long *ncols,
                             la_col_t *cols);

/* Find up to 64 exact nullspace dependencies by dense elimination. */
extern uint64_t *la_dense_nullspace(unsigned long nrows,
                                    unsigned long ncols,
                                    const la_col_t *cols,
                                    uint64_t *mask);

/* Find nullspace dependencies with the sparse block-Lanczos solver. */
extern uint64_t *la_block_lanczos(unsigned long nrows,
                                  unsigned long dense_rows,
                                  unsigned long ncols,
                                  la_col_t *cols,
                                  uint32_t seed1,
                                  uint32_t seed2,
                                  uint64_t *mask);

/* Find a wider dependency sample by retaining all rows in the iteration. */
extern uint64_t *la_block_lanczos_wide(unsigned long nrows,
                                       unsigned long dense_rows,
                                       unsigned long ncols,
                                       la_col_t *cols,
                                       uint32_t seed1,
                                       uint32_t seed2,
                                       uint64_t *mask);

#ifdef PSIQS
/* Use a per-solve worker pool for large unpacked matrices; otherwise serial. */
extern uint64_t *la_block_lanczos_threaded(unsigned long nrows,
                                          unsigned long dense_rows,
                                          unsigned long ncols,
                                          la_col_t *cols,
                                          uint32_t seed1,
                                          uint32_t seed2,
                                          uint64_t *mask,
                                          uint32_t nthreads,
                                          int retain_all_rows);
#endif

#endif
