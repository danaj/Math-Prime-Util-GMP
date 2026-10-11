/*============================================================================

  Standalone driver and host dependencies for MPU-SIQS.

  Direct build:

    cc -O3 -DSTANDALONE -o msiqs \
      mpu-siqs.c siqs.c lanczos.c prime_iterator.c squfof126.c -lgmp -lm

  Add -march=native when building a binary for the local machine only.
  Add -DPSIQS -pthread for the parallel build.  After perl Makefile.PL,
  make siqs detects pthread support automatically; make siqs-serial forces
  a serial-only build.  Both targets produce msiqs.

  This driver owns one prime-cache startup/shutdown around all its inputs.
  Custom C hosts must own the same lifetime and supply the siqs_dep.h adapters;
  see tools/README-siqs-embedding.txt. Do not link this main into another driver.

  Copyright (c) 2026 Dana Jacobsen

============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "ptypes.h"
#include "prime_iterator.h"
#include "siqs.h"
#include "siqs_dep.h"
#define SIQS_PROGRAM_NAME "msiqs"

static int verbose_level = 0;  /* -v summaries, -vv progress, -vvv diagnostics. */

int siqs_is_prob_prime(const mpz_t n) {
  return mpz_probab_prime_p(n, 25);
}


static void print_usage(FILE *stream, const char *program) {
  fprintf(stream,
#ifdef PSIQS
          "usage: %s [-v|--verbose] [-threads N] [--] [INTEGER ...]\n"
          "Factor positive decimal integers with MPU-PSIQS.\n"
          "Use -threads N to select 1-%u workers (default 1); "
          "one worker uses serial SIQS.\n"
          "-t and --threads are aliases for -threads.\n"
#else
          "usage: %s [-v|--verbose] [--] [INTEGER ...]\n"
          "Factor positive decimal integers with MPU-SIQS.\n"
#endif
          "Repeat -v for increasingly detailed progress output.\n"
          "With no INTEGER arguments, read whitespace-separated integers "
          "from stdin.\n",
#ifdef PSIQS
          program, (unsigned)PSIQS_MAX_THREADS);
#else
          program);
#endif
}

#ifdef PSIQS
/* Reject signs, trailing text and out-of-range counts without narrowing. */
static int parse_thread_count(const char *text, uint32_t *nthreads) {
  uint32_t count = 0;
  const char *p;
  if (*text == '\0')
    return 0;
  for (p = text; *p != '\0'; p++) {
    if (*p < '0' || *p > '9')
      return 0;
    count = 10U * count + (uint32_t)(*p - '0');
    if (count > PSIQS_MAX_THREADS)
      return 0;
  }
  if (count == 0)
    return 0;
  *nthreads = count;
  return 1;
}
#endif

typedef struct {
  mpz_t *values;
  size_t count;
  size_t allocated;
} siqs_factor_list_t;

static void append_factor(siqs_factor_list_t *list, const mpz_t factor) {
  if (list->count == list->allocated) {
    size_t next = list->allocated ? list->allocated * 2 : 16;
    mpz_t *values;
    if (next < list->allocated ||
        next > ((size_t)-1) / sizeof(*values)) {
      fprintf(stderr, SIQS_PROGRAM_NAME ": too many factors\n");
      exit(3);
    }
    values = (mpz_t *)realloc(list->values, next * sizeof(*values));
    if (values == NULL) {
      fprintf(stderr, SIQS_PROGRAM_NAME ": unable to allocate factors\n");
      exit(3);
    }
    list->values = values;
    list->allocated = next;
  }
  mpz_init_set(list->values[list->count++], factor);
}

/* Revisit only strictly smaller composites, never an unsplit cofactor. */
static int collect_factors(const mpz_t n, siqs_factor_list_t *output,
                           uint32_t nthreads) {
  siqs_factor_list_t pending = {NULL, 0, 0};
  mpz_t current;
  int complete = 1;

#ifndef PSIQS
  (void)nthreads;
#endif
  mpz_init(current);
  append_factor(&pending, n);
  while (pending.count != 0) {
    mpz_t *partition;
    uint32_t count, i;

    mpz_swap(current, pending.values[--pending.count]);
    mpz_clear(pending.values[pending.count]);
#ifdef PSIQS
    partition = gmp_psiqs(current, &count, 2, verbose_level, nthreads);
#else
    partition = gmp_siqs(current, &count, 2, verbose_level);
#endif
    for (i = 0; i < count; i++) {
      if (siqs_is_prob_prime(partition[i])) {
        append_factor(output, partition[i]);
      } else if (mpz_cmp_ui(partition[i], 1) > 0 &&
                 mpz_cmp(partition[i], current) < 0) {
        append_factor(&pending, partition[i]);
      } else {
        append_factor(output, partition[i]);
        complete = 0;
      }
    }
    gmp_siqs_free(partition, count);
  }
  mpz_clear(current);
  free(pending.values);
  return complete;
}

static void sort_factors(siqs_factor_list_t *output) {
  size_t i, j;
  for (i = 1; i < output->count; i++) {
    for (j = i; j > 0 &&
         mpz_cmp(output->values[j], output->values[j - 1]) < 0; j--)
      mpz_swap(output->values[j], output->values[j - 1]);
  }
}

static int factor_number(const mpz_t n, uint32_t nthreads) {
  siqs_factor_list_t output = {NULL, 0, 0};
  size_t i;
  int complete;

  if (mpz_sgn(n) <= 0) {
    gmp_fprintf(stderr, SIQS_PROGRAM_NAME ": input must be positive: %Zd\n", n);
    return 1;
  }
  if (mpz_cmp_ui(n, 1) == 0) {
    puts("1:");
    return 0;
  }

  complete = collect_factors(n, &output, nthreads);
  sort_factors(&output);
  /* Print only after processing every cofactor; SIQS's verbose lines therefore
   * cannot interrupt the final result line. */
  gmp_printf("%Zd:", n);
  for (i = 0; i < output.count; i++)
    gmp_printf(" %Zd", output.values[i]);
  if (!complete)
    printf(" [incomplete]");
  putchar('\n');
  fflush(stdout);
  for (i = 0; i < output.count; i++)
    mpz_clear(output.values[i]);
  free(output.values);

  if (!complete) {
    gmp_fprintf(stderr, SIQS_PROGRAM_NAME ": incomplete factorization of %Zd\n", n);
    return 1;
  }
  return 0;
}

int main(int argc, char **argv) {
  mpz_t n;
  int argument = 1, status = 0;
  uint32_t nthreads = 1;

  while (argument < argc && argv[argument][0] == '-') {
    const char *option = argv[argument];
    if (strcmp(option, "--") == 0) {
      argument++;
      break;
    }
    if (strcmp(option, "-h") == 0 || strcmp(option, "--help") == 0) {
      print_usage(stdout, argv[0]);
      return 0;
    }
#ifdef PSIQS
    if (strcmp(option, "-t") == 0 ||
        strcmp(option, "-threads") == 0 ||
        strcmp(option, "--threads") == 0) {
      if (argument + 1 == argc ||
          !parse_thread_count(argv[argument + 1], &nthreads)) {
        fprintf(stderr, SIQS_PROGRAM_NAME ": %s requires a worker count "
                "between 1 and %u\n", option, (unsigned)PSIQS_MAX_THREADS);
        return 2;
      }
      argument += 2;
      continue;
    }
#endif
    if (strcmp(option, "--verbose") == 0) {
      verbose_level++;
    } else if (option[1] == 'v' && option[2] != '\0') {
      const char *p;
      for (p = option + 1; *p == 'v'; p++)
        verbose_level++;
      if (*p != '\0') {
        fprintf(stderr, SIQS_PROGRAM_NAME ": unknown option: %s\n", option);
        print_usage(stderr, argv[0]);
        return 2;
      }
    } else if (strcmp(option, "-v") == 0) {
      verbose_level++;
    } else {
      fprintf(stderr, SIQS_PROGRAM_NAME ": unknown option: %s\n", option);
      print_usage(stderr, argv[0]);
      return 2;
    }
    argument++;
  }

  /* One host-owned cache lifetime, shared by every input/worker in this run. */
  prime_iterator_global_startup();
  mpz_init(n);
  if (argument < argc) {
    for (; argument < argc; argument++) {
      if (mpz_set_str(n, argv[argument], 10) != 0) {
        fprintf(stderr, SIQS_PROGRAM_NAME ": invalid decimal integer: %s\n",
                argv[argument]);
        status = 2;
        continue;
      }
      if (factor_number(n, nthreads) != 0 && status == 0)
        status = 1;
    }
  } else {
    int scanned;
    while ((scanned = gmp_scanf("%Zd", n)) == 1)
      if (factor_number(n, nthreads) != 0 && status == 0)
        status = 1;
    if (scanned != EOF) {
      fprintf(stderr, SIQS_PROGRAM_NAME ": invalid input on stdin\n");
      status = 2;
    }
  }
  mpz_clear(n);
  prime_iterator_global_shutdown();
  return status;
}
