/*
 * Parallel block-Lanczos kernels for lanczos.c, enabled by PSIQS.
 * This is not a separate compilation unit.  The recurrence, random seeds,
 * packed matrices, and dependency verification stay in the serial solver.
 *
 * A pool lives for the whole solve, including retries.  The calling thread
 * is worker zero.  Forward multiplication scatters only into private row
 * buffers, then a row-partitioned XOR reduction precedes the transpose.
 * Every other kernel writes disjoint ranges or private partial products.
 * No atomics or floating-point reductions change the serial result.
 *
 * Copyright (c) 2026 Dana Jacobsen.  See LICENSE for redistribution terms.
 */

#include <pthread.h>

/* Sparse gathers/scatters soon become bandwidth limited.  Bound the initial
 * pool and its per-worker row storage independently of the sieve threads. */
#ifndef NLA_MAX_THREADS
/* Provisional cap; fewer workers may be faster. Revisit after scaling tests. */
#define NLA_MAX_THREADS 32U
#endif
#if NLA_MAX_THREADS < 1
#error "NLA_MAX_THREADS must be positive"
#endif

typedef enum {
  NLA_JOB_FORWARD, NLA_JOB_REDUCE, NLA_JOB_TRANSPOSE,
  NLA_JOB_INNER, NLA_JOB_VECTOR
} nla_job_t;

typedef struct {
  nla_pool_t *pool;
  pthread_t thread;
  unsigned long col_begin, col_end;
  unsigned long row_begin, row_end;
  unsigned long vector_begin, vector_end;
  uint64_t *rows;
  uint64_t bins[8U * 256U];
  uint64_t product[64];
} nla_worker_t;

struct nla_pool_t {
  const nla_matrix_t *matrix;
  nla_worker_t *workers;
  uint32_t nthreads, created, pending;
  pthread_mutex_t mutex;
  pthread_cond_t work, done;
  uint64_t generation;
  int stop;
  nla_job_t job;
  const uint64_t *left, *right, *table;
  uint64_t *output;
  uint64_t mask;
};

/* Execute one phase; all job inputs remain immutable until every worker ends. */
static void nla_pool_execute(nla_worker_t *worker) {
  nla_pool_t *pool = worker->pool;
  const nla_matrix_t *matrix = pool->matrix;
  unsigned long i, column;

  switch (pool->job) {
    case NLA_JOB_FORWARD:
      memset(worker->rows, 0, (size_t)matrix->input_rows * sizeof(uint64_t));
      for (column = worker->col_begin; column < worker->col_end; column++) {
        const la_col_t *c = matrix->cols + column;
        const uint32_t *dense = matrix->input_dense_rows != 0 ?
                                 c->data + c->weight : NULL;
        uint64_t value = pool->left[column];
        for (i = 0; i < c->weight; i++)
          worker->rows[c->data[i]] ^= value;
        for (i = 0; i < matrix->input_dense_rows; i++)
          if (dense[i >> 5] & ((uint32_t)1U << (i & 31UL)))
            worker->rows[i] ^= value;
      }
      break;
    case NLA_JOB_REDUCE:
      for (i = worker->row_begin; i < worker->row_end; i++) {
        uint32_t t;
        uint64_t value = pool->workers[0].rows[i];
        for (t = 1; t < pool->nthreads; t++)
          value ^= pool->workers[t].rows[i];
        pool->output[i] = value;
      }
      break;
    case NLA_JOB_TRANSPOSE:
      for (column = worker->col_begin; column < worker->col_end; column++) {
        const la_col_t *c = matrix->cols + column;
        const uint32_t *dense = matrix->input_dense_rows != 0 ?
                                 c->data + c->weight : NULL;
        uint64_t value = 0;
        for (i = 0; i < c->weight; i++)
          value ^= pool->left[c->data[i]];
        for (i = 0; i < matrix->input_dense_rows; i++)
          if (dense[i >> 5] & ((uint32_t)1U << (i & 31UL)))
            value ^= pool->left[i];
        pool->output[column] = value;
      }
      break;
    case NLA_JOB_INNER:
      nla_inner_product(pool->left + worker->vector_begin,
                         pool->right + worker->vector_begin,
                         worker->product,
                         worker->vector_end - worker->vector_begin,
                         worker->bins);
      break;
    case NLA_JOB_VECTOR:
      for (i = worker->vector_begin; i < worker->vector_end; i++)
        pool->output[i] = (pool->output[i] & pool->mask) ^
                           nla_apply_small(pool->left[i], pool->table);
      break;
  }
}

/* Generation predicates handle both spurious wakeups and fast successive jobs. */
static void *nla_pool_worker(void *argument) {
  nla_worker_t *worker = (nla_worker_t *)argument;
  nla_pool_t *pool = worker->pool;
  uint64_t generation = 0;
  pthread_mutex_lock(&pool->mutex);
  for (;;) {
    while (!pool->stop && generation == pool->generation)
      pthread_cond_wait(&pool->work, &pool->mutex);
    if (pool->stop)
      break;
    generation = pool->generation;
    pthread_mutex_unlock(&pool->mutex);
    nla_pool_execute(worker);
    pthread_mutex_lock(&pool->mutex);
    if (--pool->pending == 0)
      pthread_cond_signal(&pool->done);
  }
  pthread_mutex_unlock(&pool->mutex);
  return NULL;
}

/* Join every created worker before freeing buffers, including setup failures. */
static void nla_pool_destroy(nla_pool_t *pool) {
  uint32_t t;
  if (pool == NULL)
    return;
  pthread_mutex_lock(&pool->mutex);
  pool->stop = 1;
  pthread_cond_broadcast(&pool->work);
  pthread_mutex_unlock(&pool->mutex);
  for (t = 1; t <= pool->created; t++)
    pthread_join(pool->workers[t].thread, NULL);
  pthread_cond_destroy(&pool->done);
  pthread_cond_destroy(&pool->work);
  pthread_mutex_destroy(&pool->mutex);
  for (t = 0; t < pool->nthreads; t++)
    free(pool->workers[t].rows);
  free(pool->workers);
  free(pool);
}

/* Optional resources must fail back to serial, not abort an otherwise valid solve. */
static nla_pool_t *nla_pool_create(const nla_matrix_t *matrix, uint32_t nthreads) {
  nla_pool_t *pool;
  uint64_t total = 0, cumulative = 0;
  unsigned long column = 0;
  uint32_t t;
  if (nthreads < 2U || matrix->packed || matrix->ncols <= NLA_PACK_MAX_COLS ||
      matrix->input_rows == 0)
    return NULL;
  if (nthreads > NLA_MAX_THREADS)
    nthreads = NLA_MAX_THREADS;
  if (nthreads < 2U || matrix->input_rows > (size_t)-1 / sizeof(uint64_t) ||
      sizeof(nla_worker_t) > (size_t)-1 / nthreads)
    return NULL;

  for (column = 0; column < matrix->ncols; column++) {
    uint64_t cost = (uint64_t)matrix->cols[column].weight +
                     matrix->input_dense_rows + 1U;
    if (cost > UINT64_MAX - total)
      return NULL;
    total += cost;
  }
  pool = (nla_pool_t *)calloc(1, sizeof(*pool));
  if (pool == NULL)
    return NULL;
  pool->workers = (nla_worker_t *)calloc(nthreads, sizeof(*pool->workers));
  if (pool->workers == NULL) {
    free(pool);
    return NULL;
  }
  if (pthread_mutex_init(&pool->mutex, NULL) != 0)
    goto fail_mutex;
  if (pthread_cond_init(&pool->work, NULL) != 0)
    goto fail_work;
  if (pthread_cond_init(&pool->done, NULL) != 0)
    goto fail_done;
  pool->matrix = matrix;
  pool->nthreads = nthreads;

  /* Reduction sorts columns by weight.  Equal column counts would leave the
   * last worker most of the sparse entries; balance cumulative work instead. */
  column = 0;
  for (t = 0; t < nthreads; t++) {
    nla_worker_t *worker = pool->workers + t;
    uint64_t target = (total / nthreads) * (t + 1U) +
                       ((total % nthreads) * (t + 1U)) / nthreads;
    worker->pool = pool;
    worker->col_begin = column;
    while (column < matrix->ncols && cumulative < target) {
      cumulative += (uint64_t)matrix->cols[column].weight +
                      matrix->input_dense_rows + 1U;
      column++;
    }
    worker->col_end = column;
    worker->row_begin = (matrix->input_rows / nthreads) * t +
                         ((matrix->input_rows % nthreads) * t) / nthreads;
    worker->row_end = (matrix->input_rows / nthreads) * (t + 1U) +
                       ((matrix->input_rows % nthreads) * (t + 1U)) / nthreads;
    worker->vector_begin = (matrix->ncols / nthreads) * t +
                            ((matrix->ncols % nthreads) * t) / nthreads;
    worker->vector_end = (matrix->ncols / nthreads) * (t + 1U) +
                          ((matrix->ncols % nthreads) * (t + 1U)) / nthreads;
    worker->rows = (uint64_t *)malloc((size_t)matrix->input_rows * sizeof(uint64_t));
    if (worker->rows == NULL) {
      nla_pool_destroy(pool);
      return NULL;
    }
  }
  for (t = 1; t < nthreads; t++) {
    if (pthread_create(&pool->workers[t].thread, NULL, nla_pool_worker,
                        pool->workers + t) != 0) {
      nla_pool_destroy(pool);
      return NULL;
    }
    pool->created++;
  }
  return pool;

fail_done:
  pthread_cond_destroy(&pool->work);
fail_work:
  pthread_mutex_destroy(&pool->mutex);
fail_mutex:
  free(pool->workers);
  free(pool);
  return NULL;
}

/* The caller participates, then waits for all private results to be published. */
static void nla_pool_dispatch(nla_pool_t *pool, nla_job_t job) {
  pthread_mutex_lock(&pool->mutex);
  pool->job = job;
  pool->pending = pool->nthreads - 1U;
  pool->generation++;
  pthread_cond_broadcast(&pool->work);
  pthread_mutex_unlock(&pool->mutex);
  nla_pool_execute(pool->workers);
  pthread_mutex_lock(&pool->mutex);
  while (pool->pending != 0)
    pthread_cond_wait(&pool->done, &pool->mutex);
  pthread_mutex_unlock(&pool->mutex);
}

static void nla_pool_mul_symmetric(nla_pool_t *pool, const uint64_t *input,
                                    uint64_t *output, uint64_t *row_scratch) {
  pool->left = input;
  nla_pool_dispatch(pool, NLA_JOB_FORWARD);
  pool->output = row_scratch;
  nla_pool_dispatch(pool, NLA_JOB_REDUCE);
  pool->left = row_scratch;
  pool->output = output;
  nla_pool_dispatch(pool, NLA_JOB_TRANSPOSE);
}

static void nla_pool_inner_product(nla_pool_t *pool, const uint64_t *left,
                                    const uint64_t *right, uint64_t *product) {
  unsigned int bit;
  uint32_t t;
  pool->left = left;
  pool->right = right;
  nla_pool_dispatch(pool, NLA_JOB_INNER);
  for (bit = 0; bit < 64U; bit++) {
    uint64_t value = pool->workers[0].product[bit];
    for (t = 1; t < pool->nthreads; t++)
      value ^= pool->workers[t].product[bit];
    product[bit] = value;
  }
}

static void nla_pool_vector_acc(nla_pool_t *pool, const uint64_t *vector,
                                 uint64_t *output, uint64_t mask,
                                 const uint64_t *table) {
  pool->left = vector;
  pool->output = output;
  pool->mask = mask;
  pool->table = table;
  nla_pool_dispatch(pool, NLA_JOB_VECTOR);
}
