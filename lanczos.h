#ifndef MPU_LANCZOS_H
#define MPU_LANCZOS_H

#include "ptypes.h"

typedef struct {
  unsigned long *data;        /* The list of occupied rows in this column */
  unsigned long weight;       /* Number of nonzero entries in this column */
  unsigned long orig;         /* Original relation number */
} la_col_t;

extern uint64_t getNullEntry(const uint64_t *nullrows, unsigned long i,
                             unsigned long l);
extern void reduce_matrix(unsigned long *nrows, unsigned long *ncols,
                          la_col_t *cols);
/* Return up to 64 exact dependencies, using the same per-column packed result
 * and mask convention as block_lanczos().  The input columns are unchanged. */
extern uint64_t *dense_nullspace64(unsigned long nrows,
                                    unsigned long ncols,
                                    const la_col_t *cols,
                                    uint64_t *mask);
/* Returns dependency vectors and a nonzero mask, or NULL after all retries. */
extern uint64_t *block_lanczos(unsigned long nrows,
                               unsigned long dense_rows,
                               unsigned long ncols, la_col_t *cols,
                               uint32_t seed1, uint32_t seed2,
                               uint64_t *mask);

#endif
