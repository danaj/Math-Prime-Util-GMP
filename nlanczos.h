#ifndef MPU_NLANCZOS_H
#define MPU_NLANCZOS_H

#include "ptypes.h"

/*
 * Keep the new solver's symbols distinct while it is tested alongside the
 * established implementation.  The aliases at the end let siqs.c switch
 * solvers by changing only its include.
 */
typedef struct {
  unsigned long *data;
  unsigned long weight;
  unsigned long orig;
} nlanczos_col_t;

/* The reducer establishes its own deterministic light-to-heavy order. */
#define MPU_LANCZOS_REDUCER_SORTS_COLUMNS 1

extern uint64_t nlanczos_get_null_entry(const uint64_t *nullrows,
                                        unsigned long i,
                                        unsigned long l);
extern void nlanczos_reduce_matrix(unsigned long *nrows,
                                   unsigned long *ncols,
                                   nlanczos_col_t *cols);
extern uint64_t *nlanczos_dense_nullspace64(unsigned long nrows,
                                             unsigned long ncols,
                                             const nlanczos_col_t *cols,
                                             uint64_t *mask);
extern uint64_t *nlanczos_block_lanczos(unsigned long nrows,
                                         unsigned long dense_rows,
                                         unsigned long ncols,
                                         nlanczos_col_t *cols,
                                         uint32_t seed1,
                                         uint32_t seed2,
                                         uint64_t *mask);

#ifndef MPU_NLANCZOS_EXPLICIT_NAMES
#define la_col_t              nlanczos_col_t
#define getNullEntry          nlanczos_get_null_entry
#define reduce_matrix         nlanczos_reduce_matrix
#define dense_nullspace64     nlanczos_dense_nullspace64
#define block_lanczos         nlanczos_block_lanczos
#endif

#endif
