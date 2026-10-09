/* Matrix suite, included by siqs-check.c; not a separate compilation unit.
 * Including Lanczos exposes packing/kernels without production test hooks.
 * References operate on immutable original columns, not packed storage. */
#ifndef SIQS_CHECK_LANCZOS_SOURCE
# define SIQS_CHECK_LANCZOS_SOURCE "../lanczos.c"
#endif
#include "siqs-check-thread-seams.h"
#include SIQS_CHECK_LANCZOS_SOURCE
#undef malloc
#undef calloc
#undef pthread_create
#undef pthread_join
#undef pthread_mutex_init
#undef pthread_mutex_destroy
#undef pthread_cond_init
#undef pthread_cond_destroy
#undef pthread_mutex_unlock
#undef pthread_cond_signal
#undef fprintf
#undef exit

static uint32_t matrix_random_state = UINT32_C(0x183a97bd);
static uint32_t matrix_random(void) {
  uint32_t x = matrix_random_state;
  x ^= x << 13; x ^= x >> 17; x ^= x << 5;
  matrix_random_state = x;
  return x;
}
static uint64_t matrix_random64(void) {
  uint64_t high = matrix_random();
  return (high << 32) | matrix_random();
}

static void matrix_free(la_col_t *cols, unsigned long count) {
  unsigned long i;
  for (i = 0; i < count; i++) free(cols[i].data);
  free(cols);
}

static la_col_t *matrix_clone(const la_col_t *cols, unsigned long count) {
  unsigned long i;
  la_col_t *copy = (la_col_t *)check_allocate(count, sizeof(*copy));
  for (i = 0; i < count; i++) {
    copy[i] = cols[i];
    copy[i].data = (uint32_t *)check_allocate(cols[i].weight, sizeof(uint32_t));
    memcpy(copy[i].data, cols[i].data, (size_t)cols[i].weight * sizeof(uint32_t));
  }
  return copy;
}

static la_col_t *matrix_fixture(unsigned long rows, unsigned long count,
                                unsigned long dense) {
  unsigned long c, j;
  la_col_t *cols = (la_col_t *)check_allocate(count, sizeof(*cols));
  CHECK(rows > dense && rows <= UINT32_MAX && count <= UINT32_MAX);
  for (c = 0; c < count; c++) {
    uint32_t weight = 1U + matrix_random() % 9U;
    if (weight > rows - dense) weight = (uint32_t)(rows - dense);
    /* Uneven weights, plus zero-weight columns, exercise pool partitioning. */
    if (c >= rows && c % 17U == 0) weight = 0;
    if (c + 1U == count && rows > 100) weight = (uint32_t)(rows - dense);
    cols[c].orig = (uint32_t)c;
    cols[c].weight = weight;
    cols[c].data = (uint32_t *)check_allocate(weight + (dense + 31U) / 32U,
                                               sizeof(uint32_t));
    for (j = 0; j < weight; j++) {
      uint32_t row;
      unsigned long k;
      do {
        row = weight == rows - dense ? (uint32_t)(dense + j)
            : j == 0 && c >= dense && c < rows ? (uint32_t)c
            : (uint32_t)(dense + matrix_random() % (rows - dense));
        for (k = 0; k < j && cols[c].data[k] != row; k++) {}
      } while (k < j);
      cols[c].data[j] = row;
    }
    for (j = 0; j < dense; j++)
      if (c == j || matrix_random() % 5U == 0)
        cols[c].data[weight + j / 32U] |= UINT32_C(1) << (j % 32U);
  }
  return cols;
}

/* Simple scalar column traversal, also handles appended dense input words. */
static void matrix_reference_mul(unsigned long rows, unsigned long count,
    unsigned long dense, const la_col_t *cols, const uint64_t *input,
    uint64_t *output) {
  unsigned long c, i;
  memset(output, 0, (size_t)rows * sizeof(*output));
  for (c = 0; c < count; c++) {
    for (i = 0; i < cols[c].weight; i++) {
      CHECK(cols[c].data[i] >= dense && cols[c].data[i] < rows);
      output[cols[c].data[i]] ^= input[c];
    }
    for (i = 0; i < dense; i++)
      if (cols[c].data[cols[c].weight + i / 32U] & (UINT32_C(1) << (i % 32U)))
        output[i] ^= input[c];
  }
}

static void matrix_reference_transpose(unsigned long count, unsigned long dense,
    const la_col_t *cols, const uint64_t *input, uint64_t *output) {
  unsigned long c, i;
  for (c = 0; c < count; c++) {
    uint64_t value = 0;
    for (i = 0; i < cols[c].weight; i++) value ^= input[cols[c].data[i]];
    for (i = 0; i < dense; i++)
      if (cols[c].data[cols[c].weight + i / 32U] & (UINT32_C(1) << (i % 32U)))
        value ^= input[i];
    output[c] = value;
  }
}

/* Packing's documented order: decreasing incidence, ties by original row.
 * Use an independent selection sort, not nla_compare_rows or the packed map. */
static uint32_t *matrix_reference_order(unsigned long rows, unsigned long count,
    unsigned long dense, const la_col_t *cols, unsigned long *active) {
  unsigned long c, i, j;
  uint32_t *counts = (uint32_t *)check_allocate(rows, sizeof(uint32_t));
  uint32_t *order = (uint32_t *)check_allocate(rows, sizeof(uint32_t));
  for (c = 0; c < count; c++) {
    for (i = 0; i < cols[c].weight; i++) counts[cols[c].data[i]]++;
    for (i = 0; i < dense; i++)
      if (cols[c].data[cols[c].weight + i / 32U] & (UINT32_C(1) << (i % 32U)))
        counts[i]++;
  }
  for (i = *active = 0; i < rows; i++) if (counts[i]) order[(*active)++] = (uint32_t)i;
  for (i = 0; i < *active; i++) {
    unsigned long best = i;
    for (j = i + 1U; j < *active; j++)
      if (counts[order[j]] > counts[order[best]] ||
          (counts[order[j]] == counts[order[best]] && order[j] < order[best])) best = j;
    { uint32_t swap = order[i]; order[i] = order[best]; order[best] = swap; }
  }
  free(counts);
  return order;
}

static void matrix_verify_dependencies(unsigned long rows, unsigned long count,
    unsigned long dense, const la_col_t *original, const uint64_t *nullrows,
    uint64_t mask) {
  unsigned long c, row;
  uint64_t used = 0, basis[64] = {0};
  uint64_t *parity = (uint64_t *)check_allocate(rows, sizeof(uint64_t));
  CHECK(nullrows != NULL && mask != 0);
  for (c = 0; c < count; c++) {
    uint64_t bits = nullrows[c];
    CHECK((bits & ~mask) == 0);
    used |= bits;
    /* Check lane independence as well as nonempty kernel vectors. */
    for (row = 0; row < 64 && bits; row++)
      if (bits & (UINT64_C(1) << row)) {
        if (basis[row]) bits ^= basis[row];
        else { basis[row] = bits; break; }
      }
  }
  CHECK(used == mask);
  for (row = 0; row < 64; row++) CHECK((basis[row] != 0) == ((mask >> row) & 1U));
  matrix_reference_mul(rows, count, dense, original, nullrows, parity);
  for (row = 0; row < rows; row++) CHECK(parity[row] == 0);
  free(parity);
}

static void matrix_kernel_case(unsigned long rows, unsigned long count,
                                unsigned long dense, unsigned int post) {
  la_col_t *cols = matrix_fixture(rows, count, dense);
  nla_matrix_t matrix;
  unsigned long active, i, c, trial;
  uint32_t *order = matrix_reference_order(rows, count, dense, cols, &active);
  uint64_t *input = (uint64_t *)check_allocate(count, sizeof(uint64_t));
  uint64_t *want = (uint64_t *)check_allocate(count, sizeof(uint64_t));
  uint64_t *got = (uint64_t *)check_allocate(count, sizeof(uint64_t));
  uint64_t *row_values = (uint64_t *)check_allocate(rows, sizeof(uint64_t));
  uint64_t *row_actual = (uint64_t *)check_allocate(rows, sizeof(uint64_t));
  uint64_t *row_input = (uint64_t *)check_allocate(rows, sizeof(uint64_t));
  uint64_t *table = (uint64_t *)check_allocate(8U * 256U, sizeof(uint64_t));
  case_name = "original-matrix-kernels";
  nla_matrix_init(&matrix, rows, dense, count, cols, post, 0);
  CHECK(matrix.packed == (active >= 1024 && count <= 32768));
  CHECK(matrix.active_rows == active);
  if (detailed) printf("  matrix %lu x %lu, dense=%lu post=%u packed=%d\n",
                       rows, count, dense, matrix.post_rows, matrix.packed);
  /* Independently check post-row bitmaps (the iterative operator omits them). */
  if (matrix.packed) {
    uint32_t *rank = (uint32_t *)check_allocate(rows, sizeof(uint32_t));
    for (i = 0; i < active; i++) rank[order[i]] = (uint32_t)i;
    for (c = 0; c < count; c++) {
      uint64_t bits = 0;
      for (i = 0; i < cols[c].weight; i++)
        if (rank[cols[c].data[i]] < post) bits ^= UINT64_C(1) << rank[cols[c].data[i]];
      for (i = 0; i < dense; i++)
        if ((cols[c].data[cols[c].weight + i / 32U] & (UINT32_C(1) << (i % 32U))) && rank[i] < post)
          bits ^= UINT64_C(1) << rank[i];
      CHECK(matrix.post_bits[c] == bits);
    }
    free(rank);
  }
  for (trial = 0; trial < (extended ? 5U : 2U); trial++) {
    for (c = 0; c < count; c++) input[c] = matrix_random64();
    matrix_reference_mul(rows, count, dense, cols, input, row_values);
    nla_matrix_mul(&matrix, input, row_actual, table);
    for (i = 0; i < matrix.iteration_rows; i++)
      CHECK(row_actual[i] == row_values[matrix.packed ? order[i + post] : i]);
    for (i = 0; i < rows; i++) row_values[i] = matrix_random64();
    if (matrix.packed) {
      for (i = 0; i < post; i++) row_values[order[i]] = 0;
      for (i = 0; i < matrix.iteration_rows; i++) row_input[i] = row_values[order[i + post]];
    } else memcpy(row_input, row_values, (size_t)rows * sizeof(uint64_t));
    matrix_reference_transpose(count, dense, cols, row_values, want);
    nla_matrix_mul_transpose(&matrix, row_input, got, table);
    CHECK(memcmp(want, got, (size_t)count * sizeof(uint64_t)) == 0);
    matrix_reference_mul(rows, count, dense, cols, input, row_values);
    if (matrix.packed) for (i = 0; i < post; i++) row_values[order[i]] = 0;
    matrix_reference_transpose(count, dense, cols, row_values, want);
    nla_matrix_mul_symmetric(&matrix, input, got, row_actual, table);
    CHECK(memcmp(want, got, (size_t)count * sizeof(uint64_t)) == 0);
    memcpy(got, input, (size_t)count * sizeof(uint64_t));
    nla_matrix_mul_symmetric(&matrix, got, got, row_actual, table);
    CHECK(memcmp(want, got, (size_t)count * sizeof(uint64_t)) == 0);
#ifdef PSIQS
    if (!matrix.packed && count > 32768) {
      uint32_t threads;
      for (threads = 2; threads <= 4; threads++) {
        uint64_t small[64], inner[64] = {0}, actual_inner[64];
        uint64_t keep = threads == 2 ? 0 : threads == 3 ? UINT64_MAX
                                                        : UINT64_C(0x55aaff00aa55ff00);
        unsigned int bit;
        matrix.pool = nla_pool_create(&matrix, threads);
        CHECK(matrix.pool != NULL && matrix.pool->nthreads == threads);
        nla_matrix_mul_symmetric(&matrix, input, got, row_actual, table);
        CHECK(memcmp(want, got, (size_t)count * sizeof(uint64_t)) == 0);
        memcpy(got, input, (size_t)count * sizeof(uint64_t));
        nla_matrix_mul_symmetric(&matrix, got, got, row_actual, table);
        CHECK(memcmp(want, got, (size_t)count * sizeof(uint64_t)) == 0);
        for (bit = 0; bit < 64; bit++) small[bit] = matrix_random64();
        for (c = 0; c < count; c++) {
          uint64_t value = 0, bits = input[c];
          for (bit = 0; bit < 64; bit++) {
            if (bits & (UINT64_C(1) << bit)) value ^= small[bit];
            if (input[c] & (UINT64_C(1) << bit)) inner[bit] ^= want[c];
          }
          got[c] = (input[c] & keep) ^ value;
        }
        nla_solver_inner_product(&matrix, input, want, actual_inner, table);
        CHECK(memcmp(inner, actual_inner, sizeof(inner)) == 0);
        /* Aliased input/output, to expose an overwrite-before-read mistake. */
        memcpy(want, input, (size_t)count * sizeof(uint64_t));
        nla_solver_vector_acc(&matrix, want, small, want, keep, table);
        CHECK(memcmp(want, got, (size_t)count * sizeof(uint64_t)) == 0);
        nla_pool_destroy(matrix.pool); matrix.pool = NULL;
        /* Restore the symmetric reference for the next thread count. */
        matrix_reference_transpose(count, dense, cols, row_values, want);
      }
    }
#endif
  }
  nla_matrix_clear(&matrix);
  free(order); free(input); free(want); free(got); free(row_values);
  free(row_actual); free(row_input); free(table); matrix_free(cols, count);
}

/* Slow fixed-point pruning oracle. No incidence lists, queue, or production
 * sorting/pruning helpers; ties discard the largest original column first. */
static unsigned char *matrix_reference_reduce(unsigned long rows,
    unsigned long count, const la_col_t *cols) {
  unsigned long c, row, i, live, active;
  unsigned char *alive = (unsigned char *)check_allocate(count, 1);
  uint32_t *counts = (uint32_t *)check_allocate(rows, sizeof(uint32_t));
  memset(alive, 1, (size_t)count);
  for (;;) {
    int peeled = 0;
    memset(counts, 0, (size_t)rows * sizeof(uint32_t));
    for (c = live = 0; c < count; c++) if (alive[c]) {
      live++;
      for (i = 0; i < cols[c].weight; i++) counts[cols[c].data[i]]++;
    }
    for (row = active = 0; row < rows; row++) {
      if (counts[row]) active++;
      if (counts[row] == 1) {
        for (c = 0; c < count; c++) if (alive[c]) {
          for (i = 0; i < cols[c].weight; i++) if (cols[c].data[i] == row) break;
          if (i < cols[c].weight) { alive[c] = 0; peeled = 1; break; }
        }
        break; /* Recount before another peel. */
      }
    }
    if (peeled) continue;
    if (live <= active + 64U) break;
    for (i = live - active - 64U; i != 0; i--) {
      unsigned long worst = count;
      for (c = 0; c < count; c++) if (alive[c] &&
          (worst == count || cols[c].weight > cols[worst].weight ||
           (cols[c].weight == cols[worst].weight && cols[c].orig > cols[worst].orig))) worst = c;
      CHECK(worst < count); alive[worst] = 0;
    }
  }
  free(counts);
  return alive;
}

static void matrix_reduction(void) {
  unsigned long c, i, rows = 9, count = 100, reduced_rows = rows, reduced_cols = count;
  la_col_t *original = (la_col_t *)check_allocate(count, sizeof(*original)), *copy;
  unsigned char *alive;
  uint64_t mask, *nullrows, *mapped;
  case_name = "singleton-trim-and-original-mapping";
  for (c = 0; c < count; c++) {
    original[c].orig = (uint32_t)c;
    original[c].weight = c < 3 ? (c == 2 ? 1U : 2U) : c % 11U == 0 ? 0U : 2U;
    original[c].data = (uint32_t *)check_allocate(original[c].weight, sizeof(uint32_t));
    if (c < 3) {
      original[c].data[0] = (uint32_t)c;
      if (c < 2) original[c].data[1] = (uint32_t)c + 1U;
    } else if (original[c].weight) {
      original[c].data[0] = 3U + (uint32_t)c % 4U;
      original[c].data[1] = 3U + ((uint32_t)c + 1U) % 4U;
    }
  }
  alive = matrix_reference_reduce(rows, count, original);
  copy = matrix_clone(original, count);
  la_reduce_matrix(&reduced_rows, &reduced_cols, copy, 0);
  CHECK(reduced_rows == rows && reduced_cols < count);
  for (c = i = 0; c < count; c++) if (alive[c]) i++;
  CHECK(i == reduced_cols && !alive[0] && !alive[1] && !alive[2]);
  for (c = 0; c < reduced_cols; c++) {
    uint32_t orig = copy[c].orig;
    CHECK(orig < count && alive[orig]); alive[orig] = 0;
    CHECK(copy[c].weight == original[orig].weight);
    CHECK(memcmp(copy[c].data, original[orig].data,
                 (size_t)copy[c].weight * sizeof(uint32_t)) == 0);
  }
  nullrows = la_dense_nullspace(rows, reduced_cols, copy, &mask);
  mapped = (uint64_t *)check_allocate(count, sizeof(uint64_t));
  CHECK(nullrows != NULL);
  for (c = 0; c < reduced_cols; c++) mapped[copy[c].orig] = nullrows[c];
  matrix_verify_dependencies(rows, count, 0, original, mapped, mask);
  for (c = reduced_cols; c < count; c++) CHECK(copy[c].data == NULL);
  free(alive); free(nullrows); free(mapped);
  matrix_free(copy, count); matrix_free(original, count);
}

static void matrix_from_relations(void) {
  siqs_ctx_t ctx;
  siqs_factor_array_t result;
  unsigned long rows, count, c, i;
  la_col_t *cols;
  siqs_factor_t odd[3] = {{0, 1}, {1, 3}, {2, 5}}, even[1] = {{1, 2}};
  uint64_t mask, *nullrows;
  case_name = "relation-parity-to-column-mapping";
  relation_context(&ctx, &result, 1);
  siqs_materialize_smooth_relation(&ctx, relation_raw(&ctx, 1, 1, odd, 3));
  siqs_materialize_smooth_relation(&ctx, relation_raw(&ctx, 1, 1, even, 1));
  siqs_materialize_smooth_relation(&ctx, relation_raw(&ctx, 1, 1, odd, 3));
  cols = siqs_build_matrix(&ctx, &rows, &count);
  CHECK(rows == 3 && count == 3);
  for (c = 0; c < count; c++) {
    CHECK(cols[c].orig == c && cols[c].weight == (c == 1 ? 0U : 3U));
    for (i = 0; i < cols[c].weight; i++) CHECK(cols[c].data[i] == i);
  }
  nullrows = la_dense_nullspace(rows, count, cols, &mask);
  CHECK(mask == UINT64_C(3)); /* rank 1, nullity 2 */
  matrix_verify_dependencies(rows, count, 0, cols, nullrows, mask);
  free(nullrows); matrix_free(cols, count); relation_finish(&ctx, &result);
}

static void matrix_rng(void) {
  static const struct {
    uint64_t seed;
    uint64_t output[4];
  } vectors[] = {
    { UINT64_C(0), {
      UINT64_C(0xe220a8397b1dcdaf), UINT64_C(0x6e789e6aa1b965f4),
      UINT64_C(0x06c45d188009454f), UINT64_C(0xf88bb8a8724c81ec) } },
    { UINT64_C(0x83d2e5b79a4c610f), {
      UINT64_C(0xbfb8539d14b5da28), UINT64_C(0xa44562599cc2f0c2),
      UINT64_C(0x8842d535ffbed2c4), UINT64_C(0xd644e10e3fd85c53) } },
    { UINT64_MAX, {
      UINT64_C(0xe4d971771b652c20), UINT64_C(0xe99ff867dbf682c9),
      UINT64_C(0x382ff84cb27281e9), UINT64_C(0x6d1db36ccba982d2) } }
  };
  const uint64_t increment = UINT64_C(0x9e3779b97f4a7c15);
  unsigned int i, j;
  uint64_t state;
  nla_matrix_t matrix;
  la_col_t *cols;

  case_name = "splitmix64-known-answers-and-wraparound";
  for (i = 0; i < sizeof(vectors) / sizeof(vectors[0]); i++) {
    state = vectors[i].seed;
    for (j = 0; j < 4; j++) {
      CHECK(nla_rand64(&state) == vectors[i].output[j]);
      CHECK(state == vectors[i].seed + (uint64_t)(j + 1U) * increment);
    }
  }
  state = UINT64_C(0x61c8864680b583eb);
  CHECK(nla_rand64(&state) == 0 && state == 0);
  CHECK(nla_rand64(&state) == vectors[0].output[0]);

  case_name = "splitmix64-state-continues-across-solver-attempts";
  cols = matrix_fixture(101, 165, 0);
  nla_matrix_init(&matrix, 101, 0, 165, cols, NLA_POST_ROWS, 0);
  for (i = 0; i < sizeof(vectors) / sizeof(vectors[0]); i++) {
    state = vectors[i].seed;
    for (j = 0; j < 2; j++) {
      uint64_t mask, *deps = nla_block_lanczos_once(&matrix, &state, &mask, 0);
      CHECK(state == vectors[i].seed +
            (uint64_t)(j + 1U) * matrix.ncols * increment);
      if (deps != NULL)
        matrix_verify_dependencies(101, 165, 0, cols, deps, mask);
      free(deps);
    }
  }
  nla_matrix_clear(&matrix);
  matrix_free(cols, 165);
  puts("PASS matrix: SplitMix64 known answers, zero/wrapping seeds, and state across retries");
}

static void matrix_solver_case(unsigned long rows, unsigned long count,
                                unsigned long dense) {
  la_col_t *cols = matrix_fixture(rows, count, dense);
  uint64_t *nullrows, mask;
  unsigned int wide;
  case_name = "solver-dependencies-in-original-matrix";
  if (detailed) printf("  solve %lu x %lu, dense=%lu\n", rows, count, dense);
  if (dense == 0 && count <= 1537) {
    nullrows = la_dense_nullspace(rows, count, cols, &mask);
    matrix_verify_dependencies(rows, count, dense, cols, nullrows, mask);
    free(nullrows);
  }
  for (wide = 0; wide < 2; wide++) {
    nullrows = wide ? la_block_lanczos_wide(rows, dense, count, cols,
                                          UINT64_C(0x83d2e5b79a4c610f), &mask, 0)
                   : la_block_lanczos(rows, dense, count, cols,
                                     UINT64_C(0x83d2e5b79a4c610f), &mask, 0);
    matrix_verify_dependencies(rows, count, dense, cols, nullrows, mask);
#ifdef PSIQS
    {
      uint32_t threads;
      for (threads = 2; threads <= 4; threads++) {
        uint64_t threaded_mask, *threaded = la_block_lanczos_threaded(
            rows, dense, count, cols, UINT64_C(0x83d2e5b79a4c610f),
            &threaded_mask, wide, 0, threads);
        matrix_verify_dependencies(rows, count, dense, cols, threaded, threaded_mask);
        CHECK(threaded_mask == mask);
        CHECK(memcmp(threaded, nullrows, (size_t)count * sizeof(uint64_t)) == 0);
        free(threaded);
      }
    }
#endif
    free(nullrows);
  }
  matrix_free(cols, count);
}

static void matrix_degenerate(void) {
  la_col_t cols[3];
  uint32_t row = 0;
  uint64_t mask = UINT64_MAX, *nullrows;
  case_name = "zero-rank-and-no-kernel";
  memset(cols, 0, sizeof(cols));
  nullrows = la_dense_nullspace(0, 0, cols, &mask);
  CHECK(nullrows == NULL && mask == 0);
  nullrows = la_block_lanczos(0, 0, 0, cols, 0, &mask, 0);
  CHECK(nullrows == NULL && mask == 0);
  nullrows = la_dense_nullspace(0, 3, cols, &mask);
  CHECK(mask == UINT64_C(7));
  matrix_verify_dependencies(0, 3, 0, cols, nullrows, mask); free(nullrows);
  cols[0].weight = 1; cols[0].data = &row;
  nullrows = la_dense_nullspace(1, 1, cols, &mask);
  CHECK(nullrows == NULL && mask == 0);
}

static void suite_matrix(void) {
  matrix_rng();
  matrix_from_relations();
  matrix_reduction();
  matrix_degenerate();
  CHECK(siqs_use_dense_solver(1536) && !siqs_use_dense_solver(1537));
  matrix_kernel_case(1023, 1087, 0, 48);
  matrix_kernel_case(1024, 1088, 0, 48);
  matrix_kernel_case(1100, 1164, 37, 0);
  matrix_kernel_case(1100, 1164, 37, 48);
  matrix_kernel_case(1024, 32768, 0, 48);
  matrix_kernel_case(1024, 32769, 37, 48);
  matrix_solver_case(101, 165, 0);
  matrix_solver_case(1023, 1087, 0);
  matrix_solver_case(1024, 1088, 0);
  matrix_solver_case(1100, 1164, 37);
  matrix_solver_case(1472, 1536, 0);
  matrix_solver_case(1473, 1537, 0);
  if (extended) {
    matrix_kernel_case(7, 32769, 0, 48);
    matrix_solver_case(129, 32769, 37);
    matrix_solver_case(2050, 2114, 0);
  }
#ifdef PSIQS
  puts("PASS matrix: serial/pthread kernels and fixed-seed solver results agree (2/3/4 threads)");
#else
  puts("PASS matrix: serial kernels/solvers (use --threaded to check pthread kernels)");
#endif
  puts("PASS matrix: parity mapping, singleton/trim mapping, dense/packed boundaries, original-matrix dependencies");
  fflush(stdout);
}
