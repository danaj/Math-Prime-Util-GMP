/*============================================================================

  Standalone driver and host dependencies for MPU-SIQS.

  Direct build:

    cc -O3 -DSTANDALONE -o mpu-siqs \
      mpu-siqs.c siqs.c lanczos.c prime_iterator.c squfof126.c pbrent63.c \
      -lgmp -lm

  Add -march=native when building a binary for the local machine only.

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

static int verbose_level = 2;  /* SIQS starts verbose at level 3+ */

int siqs_verbose_level(void) {
  return verbose_level;
}

int siqs_is_prob_prime(const mpz_t n) {
  return mpz_probab_prime_p(n, 25);
}

#define SIQS_TEST_FOR_2357(n, f) \
  do { \
    if (mpz_divisible_ui_p((n), 2)) { mpz_set_ui((f), 2); return 1; } \
    if (mpz_divisible_ui_p((n), 3)) { mpz_set_ui((f), 3); return 1; } \
    if (mpz_divisible_ui_p((n), 5)) { mpz_set_ui((f), 5); return 1; } \
    if (mpz_divisible_ui_p((n), 7)) { mpz_set_ui((f), 7); return 1; } \
    if (mpz_cmp_ui((n), 121) < 0) return 0; \
  } while (0)

/* Portable GMP Pollard--Brent fallback for residual cofactor splitting.
 * Keep this algorithm aligned with _GMP_pbrent_factor until the cofactor
 * splitters have a shared lower-level interface. */
int siqs_pbrent_factor(const mpz_t n, mpz_t f, UV a, UV rounds) {
  mpz_t xi, xm, saved_xi, product, temporary;
  UV i, r;
  const UV inner = 256;

  SIQS_TEST_FOR_2357(n, f);
  mpz_init_set_ui(xi, 2);
  mpz_init_set_ui(xm, 2);
  mpz_init(product);
  mpz_init(temporary);
  mpz_init(saved_xi);

  r = 1;
  mpz_set_ui(f, 1);
  while (rounds > 0) {
    UV rleft = r > rounds ? rounds : r;
    while (rleft > 0) {
      UV dorounds = rleft > inner ? inner : rleft;
      mpz_set_ui(product, 1);
      mpz_set(saved_xi, xi);
      for (i = 0; i < dorounds; i++) {
        mpz_mul(temporary, xi, xi);
        mpz_add_ui(temporary, temporary, a);
        mpz_tdiv_r(xi, temporary, n);
        mpz_sub(f, xm, xi);
        mpz_mul(product, product, f);
        if ((i % 4) == ((dorounds - 1) % 4))
          mpz_tdiv_r(product, product, n);
      }
      rleft -= dorounds;
      rounds -= dorounds;
      mpz_gcd(f, product, n);
      if (mpz_cmp_ui(f, 1) != 0)
        break;
    }
    if (mpz_cmp_ui(f, 1) == 0) {
      r *= 2;
      mpz_set(xm, xi);
      continue;
    }
    if (mpz_cmp(f, n) == 0) {
      mpz_set(xi, saved_xi);
      do {
        mpz_mul(temporary, xi, xi);
        mpz_add_ui(temporary, temporary, a);
        mpz_tdiv_r(xi, temporary, n);
        mpz_sub(f, xm, xi);
        if (mpz_sgn(f) < 0)
          mpz_add(f, f, n);
        mpz_gcd(f, f, n);
      } while (mpz_cmp_ui(f, 1) == 0 && r-- != 0);
    }
    break;
  }

  mpz_clear(xi);
  mpz_clear(xm);
  mpz_clear(saved_xi);
  mpz_clear(product);
  mpz_clear(temporary);
  if (mpz_cmp_ui(f, 1) == 0 || mpz_cmp(f, n) == 0) {
    mpz_set(f, n);
    return 0;
  }
  return 1;
}

static void print_usage(FILE *stream, const char *program) {
  fprintf(stream,
          "usage: %s [-v|--verbose] [--] [INTEGER ...]\n"
          "Factor positive decimal integers with MPU-SIQS.\n"
          "Repeat -v for increasingly detailed progress output.\n"
          "With no INTEGER arguments, read whitespace-separated integers "
          "from stdin.\n",
          program);
}

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
      fprintf(stderr, "mpu-siqs: too many factors\n");
      exit(3);
    }
    values = (mpz_t *)realloc(list->values, next * sizeof(*values));
    if (values == NULL) {
      fprintf(stderr, "mpu-siqs: unable to allocate factors\n");
      exit(3);
    }
    list->values = values;
    list->allocated = next;
  }
  mpz_init_set(list->values[list->count++], factor);
}

/* Revisit only strictly smaller composites, never an unsplit cofactor. */
static int collect_factors(const mpz_t n, siqs_factor_list_t *output) {
  siqs_factor_list_t pending = {NULL, 0, 0};
  mpz_t current;
  int complete = 1;

  mpz_init(current);
  append_factor(&pending, n);
  while (pending.count != 0) {
    mpz_t *partition;
    uint32_t count, i;

    mpz_swap(current, pending.values[--pending.count]);
    mpz_clear(pending.values[pending.count]);
    partition = _GMP_siqs(current, &count, 2);
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
    _GMP_siqs_free(partition, count);
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

static int factor_number(const mpz_t n) {
  siqs_factor_list_t output = {NULL, 0, 0};
  size_t i;
  int complete;

  if (mpz_sgn(n) <= 0) {
    gmp_fprintf(stderr, "mpu-siqs: input must be positive: %Zd\n", n);
    return 1;
  }
  if (mpz_cmp_ui(n, 1) == 0) {
    puts("1:");
    return 0;
  }

  complete = collect_factors(n, &output);
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
    gmp_fprintf(stderr, "mpu-siqs: incomplete factorization of %Zd\n", n);
    return 1;
  }
  return 0;
}

int main(int argc, char **argv) {
  mpz_t n;
  int argument = 1, status = 0;

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
    if (strcmp(option, "--verbose") == 0) {
      verbose_level++;
    } else if (option[1] == 'v' && option[2] != '\0') {
      const char *p;
      for (p = option + 1; *p == 'v'; p++)
        verbose_level++;
      if (*p != '\0') {
        fprintf(stderr, "mpu-siqs: unknown option: %s\n", option);
        print_usage(stderr, argv[0]);
        return 2;
      }
    } else if (strcmp(option, "-v") == 0) {
      verbose_level++;
    } else {
      fprintf(stderr, "mpu-siqs: unknown option: %s\n", option);
      print_usage(stderr, argv[0]);
      return 2;
    }
    argument++;
  }

  prime_iterator_global_startup();
  mpz_init(n);
  if (argument < argc) {
    for (; argument < argc; argument++) {
      if (mpz_set_str(n, argv[argument], 10) != 0) {
        fprintf(stderr, "mpu-siqs: invalid decimal integer: %s\n",
                argv[argument]);
        status = 2;
        continue;
      }
      if (factor_number(n) != 0 && status == 0)
        status = 1;
    }
  } else {
    int scanned;
    while ((scanned = gmp_scanf("%Zd", n)) == 1)
      if (factor_number(n) != 0 && status == 0)
        status = 1;
    if (scanned != EOF) {
      fprintf(stderr, "mpu-siqs: invalid input on stdin\n");
      status = 2;
    }
  }
  mpz_clear(n);
  prime_iterator_global_shutdown();
  return status;
}
