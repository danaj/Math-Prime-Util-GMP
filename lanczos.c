/*
 * Independent sparse block-Lanczos solver for Math::Prime::Util::GMP.
 *
 * The recurrence follows Peter Montgomery's block-Lanczos algorithm over
 * GF(2).  The cache layout and post-Lanczos treatment of dense rows were
 * informed by the public-domain msieve implementation, but this code is
 * written for MPU's smaller, single-process SIQS matrices and API.
 *
 * Copyright (c) 2026 Dana Jacobsen.  See LICENSE for redistribution terms.
 */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lanczos.h"
#include "utility.h"

#ifndef UINT32_MAX
#define UINT32_MAX ((uint32_t)-1)
#endif
#ifndef UINT64_MAX
#define UINT64_MAX ((uint64_t)-1)
#endif
#ifndef UINT64_C
#define UINT64_C(value) ((uint64_t)(value))
#endif

#define NLA_EXTRA_COLUMNS 64UL
#define NLA_MAX_ATTEMPTS 100U
#define NLA_RAND_MULT 2131995753U

/* Packing pays for its conversion at about 1024 active rows and is fastest
 * across the matrix sizes SIQS commonly reaches.  By roughly 40K columns the
 * input's column layout is faster.  A round 32K crossover avoids a second
 * packed format.  On the measured SIQS workloads near this crossover,
 * subsecond linear-algebra differences are immaterial next to the sieving
 * time; the upper cutoff also keeps the representation and kernels simple. */
#define NLA_PACK_MIN_ROWS 1024UL
#define NLA_PACK_MAX_COLS 32768UL
#define NLA_POST_ROWS 48U
#define NLA_PACKED_DENSE_ROWS 64U
#define NLA_DENSE_PANEL_BITS 8U
#if NLA_DENSE_PANEL_BITS < 2 || NLA_DENSE_PANEL_BITS > 8
#error "NLA_DENSE_PANEL_BITS must be between 2 and 8"
#endif
#define NLA_DENSE_PANEL_COMBINATIONS (1U << NLA_DENSE_PANEL_BITS)
#define NLA_BIT(i) (UINT64_C(1) << (i))

typedef struct {
  unsigned long row;
  unsigned long count;
} nla_row_info_t;

typedef struct {
  const la_col_t *cols;
  unsigned long input_rows;
  unsigned long input_dense_rows;
  unsigned long ncols;
  unsigned long active_rows;
  unsigned long iteration_rows;
  unsigned long image_rows;
  unsigned int post_rows;
  unsigned int packed_dense_rows;
  unsigned int sparse_rows;
  int packed;
  uint64_t *post_bits;
  uint64_t *dense_bits;
  size_t *row_offsets;
  uint16_t *row_columns;
} nla_matrix_t;

static size_t nla_array_bytes(size_t count, size_t item_size) {
  if (item_size != 0 && count > (size_t)-1 / item_size)
    croak("lanczos: allocation size overflow");
  return count * item_size;
}

static unsigned long nla_ceil_div(unsigned long value,
                                  unsigned long divisor) {
  return value / divisor + (value % divisor != 0);
}

static void *nla_malloc(size_t count, size_t item_size) {
  size_t bytes = nla_array_bytes(count, item_size);
  void *p = malloc(bytes != 0 ? bytes : 1);
  if (p == NULL)
    croak("lanczos: unable to allocate memory");
  return p;
}

static void *nla_calloc(size_t count, size_t item_size) {
  void *p;
  (void)nla_array_bytes(count, item_size);
  p = calloc(count != 0 ? count : 1, item_size);
  if (p == NULL)
    croak("lanczos: unable to allocate memory");
  return p;
}

static unsigned int nla_ctz64(uint64_t value) {
#if defined(__GNUC__) || defined(__clang__)
  return (unsigned int)__builtin_ctzll(value);
#else
  unsigned int bit = 0;
  while ((value & UINT64_C(1)) == 0) {
    value >>= 1;
    bit++;
  }
  return bit;
#endif
}

static uint64_t nla_low_mask(unsigned int bits) {
  return bits >= 64U ? UINT64_MAX : (NLA_BIT(bits) - UINT64_C(1));
}

static uint32_t nla_rand32(uint32_t *low, uint32_t *high) {
  uint64_t product = (uint64_t)(*low) * NLA_RAND_MULT + *high;
  *low = (uint32_t)product;
  *high = (uint32_t)(product >> 32);
  return *low;
}

static int nla_compare_columns(const void *a, const void *b) {
  const la_col_t *x = (const la_col_t *)a;
  const la_col_t *y = (const la_col_t *)b;
  if (x->weight < y->weight) return -1;
  if (x->weight > y->weight) return 1;
  if (x->orig < y->orig) return -1;
  if (x->orig > y->orig) return 1;
  return 0;
}

typedef struct {
  unsigned long nrows;
  unsigned long ncols;
  la_col_t *cols;
  unsigned long *counts;
  size_t *offsets;
  unsigned long *incidence;
  unsigned long *queue;
  unsigned long queue_head;
  unsigned long queue_tail;
  unsigned char *alive;
  unsigned long live_cols;
} nla_prune_t;

static void nla_prune_column(nla_prune_t *p, unsigned long column) {
  unsigned long i;
  la_col_t *c;
  if (!p->alive[column])
    return;
  p->alive[column] = 0;
  p->live_cols--;
  c = p->cols + column;
  for (i = 0; i < c->weight; i++) {
    unsigned long row = c->data[i];
    if (p->counts[row] == 0)
      croak("lanczos: inconsistent matrix row count");
    p->counts[row]--;
    if (p->counts[row] == 1) {
      if (p->queue_tail >= p->nrows)
        croak("lanczos: singleton queue overflow");
      p->queue[p->queue_tail++] = row;
    }
  }
}

static void nla_prune_singletons(nla_prune_t *p) {
  while (p->queue_head < p->queue_tail) {
    unsigned long row = p->queue[p->queue_head++];
    size_t i;
    if (p->counts[row] != 1)
      continue;
    for (i = p->offsets[row]; i < p->offsets[row + 1]; i++) {
      unsigned long column = p->incidence[i];
      if (p->alive[column]) {
        nla_prune_column(p, column);
        break;
      }
    }
  }
}

/*
 * Peel singleton rows, then discard the heaviest excess columns until the
 * remaining 2-core has at most 64 more columns than active rows.  Row numbers
 * are intentionally left unchanged: SIQS debug checks still refer to factor-
 * base row numbers.  The solver builds its own compact row map later.
 */
void la_reduce_matrix(unsigned long *nrows, unsigned long *ncols,
                      la_col_t *cols) {
  nla_prune_t p;
  unsigned long row, column, i;
  unsigned long live_rows;
  size_t entries = 0;
  size_t *next;

  if (*ncols == 0)
    return;
  qsort(cols, (size_t)*ncols, sizeof(*cols), nla_compare_columns);

  memset(&p, 0, sizeof(p));
  p.nrows = *nrows;
  p.ncols = *ncols;
  p.cols = cols;
  p.live_cols = *ncols;
  p.counts = (unsigned long *)nla_calloc((size_t)*nrows,
                                         sizeof(*p.counts));
  p.offsets = (size_t *)nla_calloc((size_t)*nrows + 1,
                                    sizeof(*p.offsets));
  p.queue = (unsigned long *)nla_malloc((size_t)*nrows,
                                        sizeof(*p.queue));
  p.alive = (unsigned char *)nla_malloc((size_t)*ncols,
                                        sizeof(*p.alive));
  memset(p.alive, 1, (size_t)*ncols);

  for (column = 0; column < *ncols; column++) {
    const la_col_t *c = cols + column;
    if ((size_t)c->weight > (size_t)-1 - entries)
      croak("lanczos: matrix weight overflow");
    entries += (size_t)c->weight;
    for (i = 0; i < c->weight; i++) {
      row = c->data[i];
      if (row >= *nrows)
        croak("lanczos: matrix row is out of range");
      p.counts[row]++;
    }
  }

  for (row = 0; row < *nrows; row++)
    p.offsets[row + 1] = p.offsets[row] + (size_t)p.counts[row];
  p.incidence = (unsigned long *)nla_malloc(entries,
                                             sizeof(*p.incidence));
  next = (size_t *)nla_malloc((size_t)*nrows, sizeof(*next));
  if (*nrows != 0)
    memcpy(next, p.offsets, (size_t)*nrows * sizeof(*next));
  for (column = 0; column < *ncols; column++) {
    const la_col_t *c = cols + column;
    for (i = 0; i < c->weight; i++) {
      row = c->data[i];
      p.incidence[next[row]++] = column;
    }
  }
  free(next);

  for (row = 0; row < *nrows; row++)
    if (p.counts[row] == 1)
      p.queue[p.queue_tail++] = row;

  for (;;) {
    unsigned long remove_count;
    nla_prune_singletons(&p);
    live_rows = 0;
    for (row = 0; row < *nrows; row++)
      if (p.counts[row] != 0)
        live_rows++;

    if (p.live_cols <= live_rows ||
        p.live_cols - live_rows <= NLA_EXTRA_COLUMNS)
      break;

    remove_count = p.live_cols - live_rows - NLA_EXTRA_COLUMNS;
    for (column = *ncols; column-- > 0 && remove_count != 0;) {
      if (p.alive[column]) {
        nla_prune_column(&p, column);
        remove_count--;
      }
    }
  }

  for (column = i = 0; column < *ncols; column++) {
    if (!p.alive[column]) {
      free(cols[column].data);
      cols[column].data = NULL;
      continue;
    }
    if (i != column) {
      cols[i] = cols[column];
      cols[column].data = NULL;
    }
    i++;
  }
  *ncols = i;

  if (get_verbose_level() > 3)
    printf("Lanczos reduced to %lu active rows x %lu columns\n",
           live_rows, *ncols);

  free(p.alive);
  free(p.queue);
  free(p.incidence);
  free(p.offsets);
  free(p.counts);
}

static int nla_verify_sparse(unsigned long nrows, unsigned long ncols,
                             const la_col_t *cols,
                             const uint64_t *dependencies) {
  uint64_t *parity = (uint64_t *)calloc(nrows != 0 ? (size_t)nrows : 1,
                                        sizeof(*parity));
  unsigned long column, i;
  int valid = 1;
  if (parity == NULL)
    return 0;
  for (column = 0; column < ncols; column++) {
    uint64_t bits = dependencies[column];
    if (bits == 0)
      continue;
    for (i = 0; i < cols[column].weight; i++)
      parity[cols[column].data[i]] ^= bits;
  }
  for (i = 0; i < nrows; i++) {
    if (parity[i] != 0) {
      valid = 0;
      break;
    }
  }
  free(parity);
  return valid;
}

/*
 * Exact elimination for small matrices.  A row-major echelon form has two
 * useful properties here: elimination can start at the pivot word, and all
 * selected dependencies can be back-substituted together in uint64_t lanes.
 */
uint64_t *la_dense_nullspace(unsigned long nrows,
                             unsigned long ncols,
                             const la_col_t *cols,
                             uint64_t *mask) {
  unsigned char *row_used = NULL;
  unsigned long *row_map = NULL, *pivot_columns = NULL;
  uint64_t *matrix = NULL, *result = NULL, *back_tables = NULL;
  uint64_t *panel_table = NULL;
  unsigned long active_rows = 0, column_words;
  unsigned long rank = 0, column, row, word, other;
  unsigned int dependency_count = 0;
  size_t matrix_words;

  *mask = 0;
  if (ncols == 0)
    return NULL;

  if (nrows > (unsigned long)((size_t)-1 / sizeof(*row_used)) ||
      nrows > (unsigned long)((size_t)-1 / sizeof(*row_map)) ||
      ncols > (unsigned long)((size_t)-1 / sizeof(*pivot_columns)) ||
      ncols > (unsigned long)((size_t)-1 / sizeof(*result)))
    return NULL;
  row_used = (unsigned char *)calloc(nrows != 0 ? (size_t)nrows : 1,
                                     sizeof(*row_used));
  row_map = (unsigned long *)malloc(
      nrows != 0 ? (size_t)nrows * sizeof(*row_map) : 1);
  pivot_columns = (unsigned long *)malloc(
      ncols != 0 ? (size_t)ncols * sizeof(*pivot_columns) : 1);
  result = (uint64_t *)calloc(ncols != 0 ? (size_t)ncols : 1,
                              sizeof(*result));
  if (row_used == NULL || row_map == NULL ||
      pivot_columns == NULL || result == NULL)
    goto allocation_failure;

  for (column = 0; column < ncols; column++) {
    for (row = 0; row < cols[column].weight; row++) {
      unsigned long index = cols[column].data[row];
      if (index >= nrows)
        croak("lanczos: dense solver matrix row is out of range");
      row_used[index] = 1;
    }
  }
  for (row = 0; row < nrows; row++)
    if (row_used[row])
      row_map[row] = active_rows++;

  column_words = nla_ceil_div(ncols, 64UL);
  if (column_words != 0 &&
      (size_t)active_rows > (size_t)-1 / (size_t)column_words)
    goto allocation_failure;
  matrix_words = (size_t)active_rows * (size_t)column_words;
  if (matrix_words > (size_t)-1 / sizeof(*matrix))
    goto allocation_failure;
  matrix = (uint64_t *)calloc(matrix_words != 0 ? matrix_words : 1,
                              sizeof(*matrix));
  if (matrix == NULL)
    goto allocation_failure;
  if ((size_t)column_words >
      (size_t)-1 / NLA_DENSE_PANEL_COMBINATIONS / sizeof(*panel_table))
    goto allocation_failure;
  panel_table = (uint64_t *)calloc(
      (size_t)column_words * NLA_DENSE_PANEL_COMBINATIONS,
      sizeof(*panel_table));
  if (panel_table == NULL)
    goto allocation_failure;

  for (column = 0; column < ncols; column++)
    for (row = 0; row < cols[column].weight; row++) {
      unsigned long mapped = row_map[cols[column].data[row]];
      matrix[(size_t)mapped * column_words + (column >> 6)] ^=
          NLA_BIT(column & 63UL);
    }

  for (column = 0; column < ncols && rank < active_rows;) {
    uint64_t pivot_bits[NLA_DENSE_PANEL_BITS];
    unsigned int panel_count = 0;
    unsigned long panel_word = column >> 6;
    unsigned long word_end = (panel_word + 1U) << 6;
    unsigned int entry;

    if (word_end > ncols)
      word_end = ncols;
    for (; column < word_end && rank + panel_count < active_rows &&
           panel_count < NLA_DENSE_PANEL_BITS; column++) {
      unsigned long selected = active_rows;
      uint64_t pivot_bit = NLA_BIT(column & 63UL);
      for (row = rank + panel_count; row < active_rows; row++) {
        const uint64_t *candidate = matrix + (size_t)row * column_words;
        uint64_t reduced = candidate[panel_word];
        unsigned int pivot;
        for (pivot = 0; pivot < panel_count; pivot++)
          if (reduced & pivot_bits[pivot])
            reduced ^= matrix[(size_t)(rank + pivot) * column_words +
                              panel_word];
        if (reduced & pivot_bit) {
          selected = row;
          break;
        }
      }
      if (selected != active_rows) {
        uint64_t *pivot_row =
            matrix + (size_t)(rank + panel_count) * column_words;
        unsigned int pivot;
        if (selected != rank + panel_count) {
          uint64_t *selected_row =
              matrix + (size_t)selected * column_words;
          for (word = panel_word; word < column_words; word++) {
            uint64_t temporary = pivot_row[word];
            pivot_row[word] = selected_row[word];
            selected_row[word] = temporary;
          }
        }
        for (pivot = 0; pivot < panel_count; pivot++) {
          const uint64_t *earlier =
              matrix + (size_t)(rank + pivot) * column_words;
          if (pivot_row[panel_word] & pivot_bits[pivot])
            for (word = panel_word; word < column_words; word++)
              pivot_row[word] ^= earlier[word];
        }
        pivot_columns[rank + panel_count] = column;
        pivot_bits[panel_count++] = pivot_bit;
      }
    }
    if (panel_count == 0)
      continue;

    /* Make the small pivot block an identity.  The remaining rows can then
     * select a precomputed pivot-row combination directly from their bits. */
    for (entry = panel_count; entry-- > 0;) {
      const uint64_t *later =
          matrix + (size_t)(rank + entry) * column_words;
      unsigned int earlier;
      for (earlier = 0; earlier < entry; earlier++) {
        uint64_t *earlier_row =
            matrix + (size_t)(rank + earlier) * column_words;
        if (earlier_row[panel_word] & pivot_bits[entry])
          for (word = panel_word; word < column_words; word++)
            earlier_row[word] ^= later[word];
      }
    }

    for (entry = 1; entry < (1U << panel_count); entry++) {
      unsigned int pivot = nla_ctz64(entry);
      unsigned int previous = entry & (entry - 1U);
      const uint64_t *pivot_row =
          matrix + (size_t)(rank + pivot) * column_words;
      uint64_t *combination =
          panel_table + (size_t)entry * column_words;
      const uint64_t *previous_combination =
          panel_table + (size_t)previous * column_words;
      for (word = panel_word; word < column_words; word++)
        combination[word] = previous_combination[word] ^ pivot_row[word];
    }

    for (other = rank + panel_count; other < active_rows; other++) {
      uint64_t *candidate = matrix + (size_t)other * column_words;
      unsigned int selected = 0;
      for (entry = 0; entry < panel_count; entry++)
        if (candidate[panel_word] & pivot_bits[entry])
          selected |= 1U << entry;
      if (selected != 0) {
        const uint64_t *combination =
            panel_table + (size_t)selected * column_words;
        for (word = panel_word; word < column_words; word++)
          candidate[word] ^= combination[word];
      }
    }
    rank += panel_count;
  }

  for (column = 0, other = 0;
       column < ncols && dependency_count < 64U; column++) {
    if (other < rank && pivot_columns[other] == column) {
      other++;
      continue;
    }
    result[column] = NLA_BIT(dependency_count);
    dependency_count++;
  }

  /* Each row says pivot + later columns = 0.  Solve all selected free
   * variables together, from the last column toward the first.  Four-column
   * tables replace the many indexed result loads in the dense tail. */
  if (dependency_count != 0) {
    unsigned long group_count = nla_ceil_div(ncols, 4UL);
    if ((size_t)group_count > (size_t)-1 / 16U / sizeof(*back_tables))
      goto allocation_failure;
    back_tables = (uint64_t *)calloc((size_t)group_count * 16U,
                                     sizeof(*back_tables));
    if (back_tables == NULL)
      goto allocation_failure;

    for (column = ncols; column-- > 0;) {
      unsigned long group = column >> 2;
      unsigned long base = group << 2;
      if (rank != 0 && pivot_columns[rank - 1U] == column) {
        const uint64_t *pivot_row;
        uint64_t dependencies = 0;
        unsigned long limit = base + 4UL;
        unsigned long scan_group;
        rank--;
        if (limit > ncols)
          limit = ncols;
        pivot_row = matrix + (size_t)rank * column_words;
        for (other = column + 1U; other < limit; other++)
          if (pivot_row[other >> 6] & NLA_BIT(other & 63UL))
            dependencies ^= result[other];
        for (scan_group = group + 1U;
             scan_group < group_count; scan_group++) {
          unsigned long group_base = scan_group << 2;
          unsigned int nibble = (unsigned int)
              ((pivot_row[group_base >> 6] >> (group_base & 63UL)) & 0xfU);
          dependencies ^= back_tables[scan_group * 16UL + nibble];
        }
        result[column] = dependencies;
      }

      if ((column & 3UL) == 0) {
        uint64_t *one_table = back_tables + group * 16UL;
        unsigned int value;
        for (value = 1; value < 16U; value++) {
          unsigned int bit = nla_ctz64(value);
          unsigned long result_column = base + bit;
          one_table[value] = one_table[value & (value - 1U)];
          if (result_column < ncols)
            one_table[value] ^= result[result_column];
        }
      }
    }
  }

  free(back_tables);
  free(panel_table);
  free(matrix);
  free(pivot_columns);
  free(row_map);
  free(row_used);

  if (dependency_count == 0 ||
      !nla_verify_sparse(nrows, ncols, cols, result)) {
    free(result);
    return NULL;
  }
  *mask = nla_low_mask(dependency_count);
  return result;

allocation_failure:
  free(back_tables);
  free(panel_table);
  free(matrix);
  free(result);
  free(pivot_columns);
  free(row_map);
  free(row_used);
  return NULL;
}

static int nla_compare_rows(const void *a, const void *b) {
  const nla_row_info_t *x = (const nla_row_info_t *)a;
  const nla_row_info_t *y = (const nla_row_info_t *)b;
  if (x->count > y->count) return -1;
  if (x->count < y->count) return 1;
  if (x->row < y->row) return -1;
  if (x->row > y->row) return 1;
  return 0;
}

static void nla_matrix_append_mapped(nla_matrix_t *matrix,
                                     size_t *row_cursor,
                                     unsigned long column,
                                     uint32_t mapped) {
  unsigned int dense_end = matrix->post_rows + matrix->packed_dense_rows;

  if (mapped < matrix->post_rows) {
    matrix->post_bits[column] ^= NLA_BIT(mapped);
  } else if (mapped < dense_end) {
    matrix->dense_bits[column] ^=
        NLA_BIT(mapped - matrix->post_rows);
  } else {
    size_t row = (size_t)(mapped - dense_end);
    size_t position = row_cursor[row]++;
    if (position >= matrix->row_offsets[row + 1U])
      croak("lanczos: packed matrix count mismatch");
    matrix->row_columns[position] = (uint16_t)column;
  }
}

static int nla_input_dense_bit(const la_col_t *column,
                               unsigned long row) {
  const unsigned long *words = column->data + column->weight;
  return (words[row >> 5] & ((unsigned long)1 << (row & 31UL))) != 0;
}

static void nla_matrix_init(nla_matrix_t *matrix,
                            unsigned long nrows,
                            unsigned long dense_rows,
                            unsigned long ncols,
                            const la_col_t *cols) {
  unsigned long *counts;
  uint32_t *row_map;
  nla_row_info_t *row_info;
  size_t *row_cursor = NULL;
  size_t offset = 0;
  unsigned long row, column, i, active = 0;

  memset(matrix, 0, sizeof(*matrix));
  matrix->cols = cols;
  matrix->input_rows = nrows;
  matrix->input_dense_rows = dense_rows;
  matrix->ncols = ncols;
  matrix->iteration_rows = nrows;
  matrix->image_rows = nrows;

  if (dense_rows > nrows)
    croak("lanczos: dense row count exceeds matrix row count");
  if (nrows > UINT32_MAX || ncols > UINT32_MAX)
    croak("lanczos: matrix dimensions exceed internal index range");

  counts = (unsigned long *)nla_calloc((size_t)nrows, sizeof(*counts));
  for (column = 0; column < ncols; column++) {
    const la_col_t *c = cols + column;
    for (i = 0; i < c->weight; i++) {
      row = c->data[i];
      if (row >= nrows)
        croak("lanczos: matrix row is out of range");
      counts[row]++;
    }
    for (row = 0; row < dense_rows; row++)
      if (nla_input_dense_bit(c, row))
        counts[row]++;
  }
  for (row = 0; row < nrows; row++)
    if (counts[row] != 0)
      active++;
  matrix->active_rows = active;

  /* Small matrices avoid conversion overhead and use the input columns.  The
   * second test also guarantees that packed matrix-vector kernels always have
   * at least 64 iteration rows for their fixed-size dense operations. */
  if (active < NLA_PACK_MIN_ROWS || active <= NLA_POST_ROWS + 64U ||
      ncols > NLA_PACK_MAX_COLS) {
    free(counts);
    return;
  }

#if NLA_PACK_MAX_COLS > 65536UL
#error "NLA_PACK_MAX_COLS must fit uint16_t column indexes"
#endif
#if NLA_POST_ROWS > 64 || NLA_PACKED_DENSE_ROWS > 64
#error "packed and post-Lanczos row groups must fit in uint64_t"
#endif

  matrix->packed = 1;
  matrix->post_rows = NLA_POST_ROWS;
  matrix->packed_dense_rows = NLA_PACKED_DENSE_ROWS;
  if (matrix->packed_dense_rows > active - matrix->post_rows)
    matrix->packed_dense_rows = (unsigned int)(active - matrix->post_rows);
  matrix->iteration_rows = active - matrix->post_rows;
  matrix->image_rows = active;
  matrix->sparse_rows = (unsigned int)(matrix->iteration_rows -
                                       matrix->packed_dense_rows);

  row_info = (nla_row_info_t *)nla_malloc((size_t)active,
                                           sizeof(*row_info));
  for (row = i = 0; row < nrows; row++) {
    if (counts[row] != 0) {
      row_info[i].row = row;
      row_info[i].count = counts[row];
      i++;
    }
  }
  qsort(row_info, (size_t)active, sizeof(*row_info), nla_compare_rows);
  row_map = (uint32_t *)nla_malloc((size_t)nrows, sizeof(*row_map));
  for (i = 0; i < active; i++)
    row_map[row_info[i].row] = (uint32_t)i;

  matrix->post_bits = (uint64_t *)nla_calloc((size_t)ncols,
                                              sizeof(*matrix->post_bits));
  matrix->dense_bits = (uint64_t *)nla_calloc((size_t)ncols,
                                               sizeof(*matrix->dense_bits));
  matrix->row_offsets = (size_t *)nla_malloc(
      (size_t)matrix->sparse_rows + 1U, sizeof(*matrix->row_offsets));

  matrix->row_offsets[0] = 0;
  for (row = 0; row < matrix->sparse_rows; row++) {
    size_t count = (size_t)row_info[matrix->post_rows +
                                        matrix->packed_dense_rows + row].count;
    if (count > (size_t)-1 - offset)
      croak("lanczos: packed matrix weight overflow");
    offset += count;
    matrix->row_offsets[row + 1U] = offset;
  }
  matrix->row_columns = (uint16_t *)nla_malloc(
      offset, sizeof(*matrix->row_columns));
  row_cursor = (size_t *)nla_malloc((size_t)matrix->sparse_rows,
                                    sizeof(*row_cursor));
  if (matrix->sparse_rows != 0)
    memcpy(row_cursor, matrix->row_offsets,
           (size_t)matrix->sparse_rows * sizeof(*row_cursor));

  for (column = 0; column < ncols; column++) {
    const la_col_t *c = cols + column;
    for (i = 0; i < c->weight; i++)
      nla_matrix_append_mapped(matrix, row_cursor, column,
                               row_map[c->data[i]]);
    for (row = 0; row < dense_rows; row++)
      if (nla_input_dense_bit(c, row))
        nla_matrix_append_mapped(matrix, row_cursor,
                                 column, row_map[row]);
  }
  for (row = 0; row < matrix->sparse_rows; row++)
    if (row_cursor[row] != matrix->row_offsets[row + 1U])
      croak("lanczos: packed matrix count mismatch");

  if (get_verbose_level() > 3) {
    double mb = ((double)offset * sizeof(*matrix->row_columns) +
                 ((double)matrix->sparse_rows + 1.0) *
                     sizeof(*matrix->row_offsets) +
                 2.0 * (double)ncols * sizeof(uint64_t)) / 1048576.0;
    printf("Lanczos packed %lu x %lu matrix: %u post rows, "
           "%u packed-dense rows, %u sparse rows and %lu entries "
           "(%.1f MB)\n",
           matrix->iteration_rows, ncols, matrix->post_rows,
           matrix->packed_dense_rows, matrix->sparse_rows,
           (unsigned long)offset, mb);
  }

  free(row_cursor);
  free(row_map);
  free(row_info);
  free(counts);
}

static void nla_matrix_clear(nla_matrix_t *matrix) {
  free(matrix->row_columns);
  free(matrix->row_offsets);
  free(matrix->dense_bits);
  free(matrix->post_bits);
  memset(matrix, 0, sizeof(*matrix));
}

/* Build the eight byte-index tables in Gray-code order (255 XORs each). */
static void nla_precompute_small(const uint64_t *matrix, uint64_t *table) {
  unsigned int byte;
  for (byte = 0; byte < 8U; byte++) {
    uint64_t accum = 0;
    unsigned int previous = 0;
    unsigned int i;
    table[byte * 256U] = 0;
    for (i = 1; i < 256U; i++) {
      unsigned int gray = i ^ (i >> 1);
      unsigned int changed = gray ^ previous;
      accum ^= matrix[byte * 8U + nla_ctz64(changed)];
      table[byte * 256U + gray] = accum;
      previous = gray;
    }
  }
}

static uint64_t nla_apply_small(uint64_t row, const uint64_t *table) {
  return table[0U * 256U + (unsigned int)( row        & 0xffU)] ^
         table[1U * 256U + (unsigned int)((row >>  8) & 0xffU)] ^
         table[2U * 256U + (unsigned int)((row >> 16) & 0xffU)] ^
         table[3U * 256U + (unsigned int)((row >> 24) & 0xffU)] ^
         table[4U * 256U + (unsigned int)((row >> 32) & 0xffU)] ^
         table[5U * 256U + (unsigned int)((row >> 40) & 0xffU)] ^
         table[6U * 256U + (unsigned int)((row >> 48) & 0xffU)] ^
         table[7U * 256U + (unsigned int)( row >> 56)];
}

/* Fixed 64x64 products do not amortize the larger byte tables. */
static void nla_precompute_nibbles(const uint64_t *matrix,
                                   uint64_t *table) {
  unsigned int nibble;
  for (nibble = 0; nibble < 16U; nibble++) {
    uint64_t accum = 0;
    unsigned int previous = 0;
    unsigned int i;
    table[nibble * 16U] = 0;
    for (i = 1; i < 16U; i++) {
      unsigned int gray = i ^ (i >> 1);
      unsigned int changed = gray ^ previous;
      accum ^= matrix[nibble * 4U + nla_ctz64(changed)];
      table[nibble * 16U + gray] = accum;
      previous = gray;
    }
  }
}

static uint64_t nla_apply_nibbles(uint64_t row, const uint64_t *table) {
  return table[ 0U * 16U + (unsigned int)( row        & 0xfU)] ^
         table[ 1U * 16U + (unsigned int)((row >>  4) & 0xfU)] ^
         table[ 2U * 16U + (unsigned int)((row >>  8) & 0xfU)] ^
         table[ 3U * 16U + (unsigned int)((row >> 12) & 0xfU)] ^
         table[ 4U * 16U + (unsigned int)((row >> 16) & 0xfU)] ^
         table[ 5U * 16U + (unsigned int)((row >> 20) & 0xfU)] ^
         table[ 6U * 16U + (unsigned int)((row >> 24) & 0xfU)] ^
         table[ 7U * 16U + (unsigned int)((row >> 28) & 0xfU)] ^
         table[ 8U * 16U + (unsigned int)((row >> 32) & 0xfU)] ^
         table[ 9U * 16U + (unsigned int)((row >> 36) & 0xfU)] ^
         table[10U * 16U + (unsigned int)((row >> 40) & 0xfU)] ^
         table[11U * 16U + (unsigned int)((row >> 44) & 0xfU)] ^
         table[12U * 16U + (unsigned int)((row >> 48) & 0xfU)] ^
         table[13U * 16U + (unsigned int)((row >> 52) & 0xfU)] ^
         table[14U * 16U + (unsigned int)((row >> 56) & 0xfU)] ^
         table[15U * 16U + (unsigned int)( row >> 60)];
}

static void nla_small_multiply(const uint64_t *left,
                               const uint64_t *right,
                               uint64_t *product,
                               uint64_t *table) {
  uint64_t temporary[64];
  unsigned int i;
  nla_precompute_nibbles(right, table);
  for (i = 0; i < 64U; i++)
    temporary[i] = nla_apply_nibbles(left[i], table);
  memcpy(product, temporary, sizeof(temporary));
}

/* Input and output must not alias. */
static void nla_small_transpose(const uint64_t *input, uint64_t *output) {
  static const uint64_t masks[6] = {
    UINT64_C(0x00000000ffffffff), UINT64_C(0x0000ffff0000ffff),
    UINT64_C(0x00ff00ff00ff00ff), UINT64_C(0x0f0f0f0f0f0f0f0f),
    UINT64_C(0x3333333333333333), UINT64_C(0x5555555555555555)
  };
  unsigned int stage;
  memcpy(output, input, 64U * sizeof(*output));
  for (stage = 0; stage < 6U; stage++) {
    unsigned int shift = 32U >> stage;
    unsigned int row;
    uint64_t mask = masks[stage];
    for (row = 0; row < 64U; row = (row + shift + 1U) & ~shift) {
      uint64_t low = output[row];
      uint64_t high = output[row + shift];
      uint64_t swap = ((low >> shift) ^ high) & mask;
      output[row] = low ^ (swap << shift);
      output[row + shift] = high ^ swap;
    }
  }
}

static void nla_vector_small_mask_acc(const uint64_t *vector,
                                      const uint64_t *small,
                                      uint64_t *output,
                                      unsigned long length,
                                      uint64_t mask,
                                      uint64_t *table) {
  unsigned long i;
  nla_precompute_small(small, table);
  for (i = 0; i < length; i++)
    output[i] = (output[i] & mask) ^ nla_apply_small(vector[i], table);
}

/* Compute transpose(left) * right for two length-n arrays of 64-bit rows. */
static void nla_inner_product(const uint64_t *left,
                              const uint64_t *right,
                              uint64_t *product,
                              unsigned long length,
                              uint64_t *bins) {
  unsigned long i;
  unsigned int byte;
  memset(bins, 0, 8U * 256U * sizeof(*bins));
  for (i = 0; i < length; i++) {
    uint64_t x = left[i];
    uint64_t y = right[i];
    bins[0U * 256U + (unsigned int)( x        & 0xffU)] ^= y;
    bins[1U * 256U + (unsigned int)((x >>  8) & 0xffU)] ^= y;
    bins[2U * 256U + (unsigned int)((x >> 16) & 0xffU)] ^= y;
    bins[3U * 256U + (unsigned int)((x >> 24) & 0xffU)] ^= y;
    bins[4U * 256U + (unsigned int)((x >> 32) & 0xffU)] ^= y;
    bins[5U * 256U + (unsigned int)((x >> 40) & 0xffU)] ^= y;
    bins[6U * 256U + (unsigned int)((x >> 48) & 0xffU)] ^= y;
    bins[7U * 256U + (unsigned int)( x >> 56)] ^= y;
  }
  for (byte = 0; byte < 8U; byte++) {
    uint64_t *byte_bins = bins + byte * 256U;
    unsigned int bit;
    for (bit = 8U; bit-- > 1U;) {
      uint64_t accum = 0;
      unsigned int half = 1U << bit;
      unsigned int value;
      for (value = 0; value < half; value++) {
        uint64_t high = byte_bins[half + value];
        accum ^= high;
        byte_bins[value] ^= high;
      }
      product[byte * 8U + bit] = accum;
    }
    product[byte * 8U] = byte_bins[1];
  }
}

static void nla_matrix_mul(const nla_matrix_t *matrix,
                           const uint64_t *input,
                           uint64_t *output,
                           uint64_t *table) {
  unsigned long column;
  unsigned int dense_rows, sparse_row, sparse_rows;
  const size_t *row_offsets;
  const uint16_t *row_columns;

  if (!matrix->packed) {
    const la_col_t *cols = matrix->cols;
    unsigned long input_rows = matrix->input_rows;
    unsigned long input_dense_rows = matrix->input_dense_rows;
    unsigned long ncols = matrix->ncols;
    memset(output, 0, (size_t)input_rows * sizeof(*output));
    for (column = 0; column < ncols; column++) {
      const la_col_t *c = cols + column;
      const unsigned long *data = c->data;
      unsigned long weight = c->weight;
      const unsigned long *dense = input_dense_rows != 0 ?
                                    data + weight : NULL;
      uint64_t value = input[column];
      unsigned long i, row;
      for (i = 0; i < weight; i++)
        output[data[i]] ^= value;
      for (row = 0; row < input_dense_rows; row++)
        if (dense[row >> 5] & ((unsigned long)1 << (row & 31UL)))
          output[row] ^= value;
    }
    return;
  }

  if (matrix->packed_dense_rows != 0)
    nla_inner_product(matrix->dense_bits, input, output,
                      matrix->ncols, table);

  dense_rows = matrix->packed_dense_rows;
  sparse_rows = matrix->sparse_rows;
  row_offsets = matrix->row_offsets;
  row_columns = matrix->row_columns;
  for (sparse_row = 0; sparse_row < sparse_rows; sparse_row++) {
    size_t i = row_offsets[sparse_row];
    size_t end = row_offsets[sparse_row + 1U];
    uint64_t accum = 0;
    for (; i < end; i++)
      accum ^= input[row_columns[i]];
    output[dense_rows + sparse_row] = accum;
  }
}

static void nla_matrix_mul_transpose(const nla_matrix_t *matrix,
                                     const uint64_t *input,
                                     uint64_t *output,
                                     uint64_t *table) {
  unsigned long column;
  unsigned int dense_rows, sparse_row, sparse_rows;
  const size_t *row_offsets;
  const uint16_t *row_columns;

  if (!matrix->packed) {
    const la_col_t *cols = matrix->cols;
    unsigned long input_dense_rows = matrix->input_dense_rows;
    unsigned long ncols = matrix->ncols;
    for (column = 0; column < ncols; column++) {
      const la_col_t *c = cols + column;
      const unsigned long *data = c->data;
      unsigned long weight = c->weight;
      const unsigned long *dense = input_dense_rows != 0 ?
                                    data + weight : NULL;
      uint64_t accum = 0;
      unsigned long i, row;
      for (i = 0; i < weight; i++)
        accum ^= input[data[i]];
      for (row = 0; row < input_dense_rows; row++)
        if (dense[row >> 5] & ((unsigned long)1 << (row & 31UL)))
          accum ^= input[row];
      output[column] = accum;
    }
    return;
  }

  if (matrix->packed_dense_rows != 0) {
    nla_precompute_small(input, table);
    for (column = 0; column < matrix->ncols; column++)
      output[column] = nla_apply_small(matrix->dense_bits[column], table);
  } else {
    memset(output, 0, (size_t)matrix->ncols * sizeof(*output));
  }

  dense_rows = matrix->packed_dense_rows;
  sparse_rows = matrix->sparse_rows;
  row_offsets = matrix->row_offsets;
  row_columns = matrix->row_columns;
  for (sparse_row = 0; sparse_row < sparse_rows; sparse_row++) {
    size_t i = row_offsets[sparse_row];
    size_t end = row_offsets[sparse_row + 1U];
    uint64_t value = input[dense_rows + sparse_row];
    for (; i < end; i++)
      output[row_columns[i]] ^= value;
  }
}

static void nla_matrix_mul_symmetric(const nla_matrix_t *matrix,
                                     const uint64_t *input,
                                     uint64_t *output,
                                     uint64_t *row_scratch,
                                     uint64_t *table) {
  nla_matrix_mul(matrix, input, row_scratch, table);
  nla_matrix_mul_transpose(matrix, row_scratch, output, table);
}

/* Invert a maximal nonsingular 64x64 submatrix, preferring new columns. */
static unsigned int nla_find_nonsingular(const uint64_t *input,
                                         unsigned int *selected,
                                         const unsigned int *previous,
                                         unsigned int previous_count,
                                         uint64_t *inverse) {
  uint64_t augmented[64][2];
  unsigned int order[64];
  uint64_t used = 0;
  unsigned int i, j, count = 0;

  if (previous_count > 64U)
    return 0;
  for (i = 0; i < 64U; i++) {
    augmented[i][0] = input[i];
    augmented[i][1] = NLA_BIT(i);
  }
  for (i = 0; i < previous_count; i++) {
    order[63U - i] = previous[i];
    used |= NLA_BIT(previous[i]);
  }
  for (i = j = 0; i < 64U; i++)
    if ((used & NLA_BIT(i)) == 0)
      order[j++] = i;

  for (i = 0; i < 64U; i++) {
    unsigned int pivot_column = order[i];
    uint64_t *pivot_row = augmented[pivot_column];

    for (j = i; j < 64U; j++) {
      uint64_t *candidate = augmented[order[j]];
      if (candidate[0] & NLA_BIT(pivot_column)) {
        uint64_t a = candidate[0], b = candidate[1];
        candidate[0] = pivot_row[0];
        candidate[1] = pivot_row[1];
        pivot_row[0] = a;
        pivot_row[1] = b;
        break;
      }
    }

    if (j < 64U) {
      for (j = 0; j < 64U; j++) {
        uint64_t *candidate = augmented[order[j]];
        if (candidate != pivot_row &&
            (candidate[0] & NLA_BIT(pivot_column))) {
          candidate[0] ^= pivot_row[0];
          candidate[1] ^= pivot_row[1];
        }
      }
      selected[count++] = pivot_column;
      continue;
    }

    /* Complete the inverse even when this column is outside the submatrix. */
    for (j = i; j < 64U; j++) {
      uint64_t *candidate = augmented[order[j]];
      if (candidate[1] & NLA_BIT(pivot_column)) {
        uint64_t a = candidate[0], b = candidate[1];
        candidate[0] = pivot_row[0];
        candidate[1] = pivot_row[1];
        pivot_row[0] = a;
        pivot_row[1] = b;
        break;
      }
    }
    if (j == 64U)
      return 0;
    for (j = 0; j < 64U; j++) {
      uint64_t *candidate = augmented[order[j]];
      if (candidate != pivot_row &&
          (candidate[1] & NLA_BIT(pivot_column))) {
        candidate[0] ^= pivot_row[0];
        candidate[1] ^= pivot_row[1];
      }
    }
    pivot_row[0] = 0;
    pivot_row[1] = 0;
  }

  for (i = 0; i < 64U; i++)
    inverse[i] = augmented[i][1];
  return count;
}

static void nla_transpose_candidates(unsigned long rows,
                                     const uint64_t *input,
                                     uint64_t **output) {
  unsigned long row;
  for (row = 0; row < rows; row++) {
    uint64_t bits = input[row];
    unsigned long word = row >> 6;
    uint64_t mask = NLA_BIT(row & 63UL);
    while (bits != 0) {
      unsigned int candidate = nla_ctz64(bits);
      output[candidate][word] |= mask;
      bits &= bits - 1;
    }
  }
}

/*
 * Eliminate the images of 128 candidate vectors, applying the same row
 * operations to the candidates themselves.  Every zero-image, nonzero
 * candidate after elimination is an exact dependency of the represented
 * matrix.  Return up to 64 of them in low packed bits.
 */
static unsigned int nla_combine_candidates(unsigned long ncols,
                                            unsigned long image_rows,
                                            uint64_t *x,
                                            const uint64_t *v,
                                            const uint64_t *bx,
                                            const uint64_t *bv) {
  const unsigned int candidate_count = 128U;
  unsigned long vector_words = nla_ceil_div(ncols, 64UL);
  unsigned long image_words = nla_ceil_div(image_rows, 64UL);
  uint64_t *vector_store, *image_store;
  uint64_t *vectors[128], *images[128];
  uint64_t *dependencies[64];
  unsigned long dependency_pivots[64];
  unsigned int rank = 0, dependency_count = 0;
  unsigned int i, j;
  unsigned long bit_position, column, word;

  if (vector_words != 0 &&
      candidate_count > (size_t)-1 / sizeof(uint64_t) / vector_words)
    croak("lanczos: candidate vector size overflow");
  if (image_words != 0 &&
      candidate_count > (size_t)-1 / sizeof(uint64_t) / image_words)
    croak("lanczos: candidate image size overflow");
  vector_store = (uint64_t *)nla_calloc(
      (size_t)candidate_count * vector_words, sizeof(*vector_store));
  image_store = (uint64_t *)nla_calloc(
      (size_t)candidate_count * image_words, sizeof(*image_store));
  for (i = 0; i < candidate_count; i++) {
    vectors[i] = vector_store + (size_t)i * vector_words;
    images[i] = image_store + (size_t)i * image_words;
  }

  nla_transpose_candidates(ncols, x, vectors);
  nla_transpose_candidates(ncols, v, vectors + 64);
  nla_transpose_candidates(image_rows, bx, images);
  nla_transpose_candidates(image_rows, bv, images + 64);

  for (bit_position = 0;
       bit_position < image_rows && rank < candidate_count;
       bit_position++) {
    unsigned long pivot_word = bit_position >> 6;
    uint64_t pivot_bit = NLA_BIT(bit_position & 63UL);
    for (j = rank; j < candidate_count; j++)
      if (images[j][pivot_word] & pivot_bit)
        break;
    if (j == candidate_count)
      continue;
    if (j != rank) {
      uint64_t *temporary = images[rank];
      images[rank] = images[j];
      images[j] = temporary;
      temporary = vectors[rank];
      vectors[rank] = vectors[j];
      vectors[j] = temporary;
    }
    for (j = rank + 1U; j < candidate_count; j++) {
      if (images[j][pivot_word] & pivot_bit) {
        for (word = 0; word < image_words; word++)
          images[j][word] ^= images[rank][word];
        for (word = 0; word < vector_words; word++)
          vectors[j][word] ^= vectors[rank][word];
      }
    }
    rank++;
  }

  /* The zero-image candidates need not themselves be independent.  Reduce
   * them once more so all returned bits represent a useful basis vector. */
  for (i = rank; i < candidate_count && dependency_count < 64U; i++) {
    for (j = 0; j < dependency_count; j++) {
      unsigned long pivot = dependency_pivots[j];
      if (vectors[i][pivot >> 6] & NLA_BIT(pivot & 63UL))
        for (word = 0; word < vector_words; word++)
          vectors[i][word] ^= dependencies[j][word];
    }
    for (word = 0; word < vector_words; word++) {
      if (vectors[i][word] != 0) {
        dependency_pivots[dependency_count] =
            word * 64UL + nla_ctz64(vectors[i][word]);
        break;
      }
    }
    if (word == vector_words)
      continue;
    dependencies[dependency_count] = vectors[i];
    dependency_count++;
  }

  memset(x, 0, (size_t)ncols * sizeof(*x));
  for (i = 0; i < dependency_count; i++) {
    for (column = 0; column < ncols; column++) {
      if (dependencies[i][column >> 6] & NLA_BIT(column & 63UL))
        x[column] |= NLA_BIT(i);
    }
  }

  free(image_store);
  free(vector_store);
  return dependency_count;
}

static int nla_verify_input_matrix(const nla_matrix_t *matrix,
                                   const uint64_t *dependencies) {
  uint64_t *parity = (uint64_t *)nla_calloc((size_t)matrix->input_rows,
                                             sizeof(*parity));
  unsigned long column, row, i;
  int valid = 1;
  for (column = 0; column < matrix->ncols; column++) {
    const la_col_t *c = matrix->cols + column;
    uint64_t bits = dependencies[column];
    if (bits == 0)
      continue;
    for (i = 0; i < c->weight; i++)
      parity[c->data[i]] ^= bits;
    for (row = 0; row < matrix->input_dense_rows; row++)
      if (nla_input_dense_bit(c, row))
        parity[row] ^= bits;
  }
  for (row = 0; row < matrix->input_rows; row++) {
    if (parity[row] != 0) {
      valid = 0;
      break;
    }
  }
  free(parity);
  return valid;
}

static int nla_all_zero(const uint64_t *matrix) {
  unsigned int i;
  for (i = 0; i < 64U; i++)
    if (matrix[i] != 0)
      return 0;
  return 1;
}

static uint64_t *nla_block_lanczos_once(const nla_matrix_t *matrix,
                                        uint32_t *seed1,
                                        uint32_t *seed2,
                                        uint64_t *result_mask) {
  uint64_t *v[3], *vnext, *x, *initial;
  uint64_t *row_scratch, *table;
  uint64_t winv_store[3][64], vt_a_v_store[2][64];
  uint64_t vt_a2_v_store[2][64], vt_v0_store[4][64];
  uint64_t *winv[3], *vt_a_v[2], *vt_a2_v[2], *vt_v0[3], *vt_v0_next;
  uint64_t d[64], e[64], f[64], temporary[64];
  unsigned int selected[2][64];
  unsigned int dim0 = 0, dim1 = 64U;
  uint64_t mask0 = 0, mask1 = UINT64_MAX;
  unsigned long i, iteration = 0, dimensions_solved = 0;
  int failed = 0;

  *result_mask = 0;
  winv[0] = winv_store[0]; winv[1] = winv_store[1];
  winv[2] = winv_store[2];
  vt_a_v[0] = vt_a_v_store[0]; vt_a_v[1] = vt_a_v_store[1];
  vt_a2_v[0] = vt_a2_v_store[0]; vt_a2_v[1] = vt_a2_v_store[1];
  vt_v0[0] = vt_v0_store[0]; vt_v0[1] = vt_v0_store[1];
  vt_v0[2] = vt_v0_store[2]; vt_v0_next = vt_v0_store[3];
  memset(winv_store, 0, sizeof(winv_store));
  memset(vt_a_v_store, 0, sizeof(vt_a_v_store));
  memset(vt_a2_v_store, 0, sizeof(vt_a2_v_store));
  memset(vt_v0_store, 0, sizeof(vt_v0_store));
  for (i = 0; i < 64UL; i++)
    selected[1][i] = (unsigned int)i;

  v[0] = (uint64_t *)nla_malloc((size_t)matrix->ncols, sizeof(*v[0]));
  v[1] = (uint64_t *)nla_calloc((size_t)matrix->ncols, sizeof(*v[1]));
  v[2] = (uint64_t *)nla_calloc((size_t)matrix->ncols, sizeof(*v[2]));
  vnext = (uint64_t *)nla_malloc((size_t)matrix->ncols, sizeof(*vnext));
  x = (uint64_t *)nla_malloc((size_t)matrix->ncols, sizeof(*x));
  initial = (uint64_t *)nla_malloc((size_t)matrix->ncols, sizeof(*initial));
  row_scratch = (uint64_t *)nla_malloc((size_t)matrix->iteration_rows,
                                        sizeof(*row_scratch));
  table = (uint64_t *)nla_malloc(8U * 256U, sizeof(*table));

  for (i = 0; i < matrix->ncols; i++) {
    uint64_t high = nla_rand32(seed1, seed2);
    x[i] = (high << 32) | nla_rand32(seed1, seed2);
  }
  nla_matrix_mul_symmetric(matrix, x, v[0], row_scratch, table);
  memcpy(initial, v[0], (size_t)matrix->ncols * sizeof(*initial));

  for (;;) {
    uint64_t *swap;
    iteration++;
    nla_matrix_mul_symmetric(matrix, v[0], vnext, row_scratch, table);
    nla_inner_product(v[0], vnext, vt_a_v[0], matrix->ncols, table);
    if (nla_all_zero(vt_a_v[0]))
      break;
    nla_inner_product(vnext, vnext, vt_a2_v[0], matrix->ncols, table);

    dim0 = nla_find_nonsingular(vt_a_v[0], selected[0], selected[1],
                                dim1, winv[0]);
    if (dim0 == 0) {
      failed = 1;
      break;
    }
    mask0 = 0;
    for (i = 0; i < dim0; i++)
      mask0 |= NLA_BIT(selected[0][i]);

    /* The coverage condition is unreliable only in the terminal region. */
    if (matrix->active_rows > matrix->post_rows + 64UL &&
        dimensions_solved < matrix->active_rows - matrix->post_rows - 64UL &&
        (mask0 | mask1) != UINT64_MAX) {
      failed = 1;
      break;
    }
    dimensions_solved += dim0;

    if (iteration <= 3)
      nla_inner_product(v[0], initial, vt_v0[0], matrix->ncols, table);

    for (i = 0; i < 64UL; i++)
      d[i] = vt_a_v[0][i] ^ (vt_a2_v[0][i] & mask0);
    nla_small_multiply(winv[0], d, d, table);
    for (i = 0; i < 64UL; i++)
      d[i] ^= NLA_BIT(i);
    nla_vector_small_mask_acc(v[0], d, vnext, matrix->ncols,
                              mask0, table);
    nla_small_transpose(d, temporary);
    nla_small_multiply(temporary, vt_v0[0], vt_v0_next, table);

    nla_small_multiply(winv[1], vt_a_v[0], e, table);
    for (i = 0; i < 64UL; i++)
      e[i] &= mask0;
    nla_vector_small_mask_acc(v[1], e, vnext, matrix->ncols,
                              UINT64_MAX, table);
    nla_small_transpose(e, temporary);
    nla_small_multiply(temporary, vt_v0[1], e, table);
    for (i = 0; i < 64UL; i++)
      vt_v0_next[i] ^= e[i];

    if (mask1 != UINT64_MAX) {
      nla_small_multiply(vt_a_v[1], winv[1], f, table);
      for (i = 0; i < 64UL; i++)
        f[i] ^= NLA_BIT(i);
      nla_small_multiply(winv[2], f, f, table);
      for (i = 0; i < 64UL; i++)
        temporary[i] = mask0 &
            (vt_a_v[1][i] ^ (vt_a2_v[1][i] & mask1));
      nla_small_multiply(f, temporary, f, table);
      nla_vector_small_mask_acc(v[2], f, vnext, matrix->ncols,
                                UINT64_MAX, table);
      nla_small_transpose(f, temporary);
      nla_small_multiply(temporary, vt_v0[2], f, table);
      for (i = 0; i < 64UL; i++)
        vt_v0_next[i] ^= f[i];
    }

    nla_small_multiply(winv[0], vt_v0[0], d, table);
    nla_vector_small_mask_acc(v[0], d, x, matrix->ncols,
                              UINT64_MAX, table);

    swap = v[2]; v[2] = v[1]; v[1] = v[0]; v[0] = vnext; vnext = swap;
    swap = winv[2]; winv[2] = winv[1]; winv[1] = winv[0]; winv[0] = swap;
    swap = vt_v0[2]; vt_v0[2] = vt_v0[1]; vt_v0[1] = vt_v0[0];
    vt_v0[0] = vt_v0_next; vt_v0_next = swap;
    swap = vt_a_v[1]; vt_a_v[1] = vt_a_v[0]; vt_a_v[0] = swap;
    swap = vt_a2_v[1]; vt_a2_v[1] = vt_a2_v[0]; vt_a2_v[0] = swap;
    memcpy(selected[1], selected[0], sizeof(selected[0]));
    dim1 = dim0;
    mask1 = mask0;

    if (iteration == 3) {
      free(initial);
      initial = NULL;
    }
    if (iteration > matrix->ncols + 64UL) {
      failed = 1;
      break;
    }
  }

  if (get_verbose_level() > 3)
    printf("Lanczos halted after %lu iterations (dimension %lu)%s\n",
           iteration, dimensions_solved, failed ? ", retrying" : "");

  free(initial);
  free(row_scratch);
  free(table);
  free(vnext);
  if (failed) {
    free(x);
    free(v[0]); free(v[1]); free(v[2]);
    return NULL;
  }

  {
    uint64_t *bx = (uint64_t *)nla_calloc((size_t)matrix->image_rows,
                                           sizeof(*bx));
    uint64_t *bv = (uint64_t *)nla_calloc((size_t)matrix->image_rows,
                                           sizeof(*bv));
    uint64_t post_x[64], post_v[64];
    unsigned int dependencies;
    uint64_t actual_mask = 0;

    table = (uint64_t *)nla_malloc(8U * 256U, sizeof(*table));
    nla_matrix_mul(matrix, x, bx, table);
    nla_matrix_mul(matrix, v[0], bv, table);
    if (matrix->post_rows != 0) {
      nla_inner_product(matrix->post_bits, x, post_x,
                        matrix->ncols, table);
      nla_inner_product(matrix->post_bits, v[0], post_v,
                        matrix->ncols, table);
      memcpy(bx + matrix->iteration_rows, post_x,
             matrix->post_rows * sizeof(*post_x));
      memcpy(bv + matrix->iteration_rows, post_v,
             matrix->post_rows * sizeof(*post_v));
    }
    dependencies = nla_combine_candidates(matrix->ncols,
                                           matrix->image_rows,
                                           x, v[0], bx, bv);
    free(table);
    free(bv);
    free(bx);
    free(v[0]); free(v[1]); free(v[2]);

    if (dependencies == 0) {
      free(x);
      return NULL;
    }
    for (i = 0; i < matrix->ncols; i++)
      actual_mask |= x[i];
    if (actual_mask == 0) {
      free(x);
      return NULL;
    }
    if (!nla_verify_input_matrix(matrix, x))
      croak("lanczos: computed dependencies failed verification");
    *result_mask = actual_mask;
    return x;
  }
}

uint64_t *la_block_lanczos(unsigned long nrows,
                           unsigned long dense_rows,
                           unsigned long ncols,
                           la_col_t *cols,
                           uint32_t seed1,
                           uint32_t seed2,
                           uint64_t *mask) {
  nla_matrix_t matrix;
  uint64_t *result = NULL;
  unsigned int attempt;

  *mask = 0;
  if (ncols == 0)
    return NULL;
  if ((seed1 | seed2) == 0 ||
      (seed1 == UINT32_MAX && seed2 == NLA_RAND_MULT - 1U)) {
    seed1 = 11111111U;
    seed2 = 22222222U;
  }

  nla_matrix_init(&matrix, nrows, dense_rows, ncols, cols);
  for (attempt = 0; attempt < NLA_MAX_ATTEMPTS; attempt++) {
    result = nla_block_lanczos_once(&matrix, &seed1, &seed2, mask);
    if (result != NULL && *mask != 0)
      break;
    free(result);
    result = NULL;
    if (get_verbose_level() > 3)
      printf("linear algebra retry %u\n", attempt + 1U);
  }
  nla_matrix_clear(&matrix);
  return result;
}
