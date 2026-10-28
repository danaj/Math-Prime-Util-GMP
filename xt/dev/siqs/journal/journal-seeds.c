/* A context-only journal must restore its own starting streams, even when
 * the caller now supplies different states. No relation checkpoint can hide
 * a failure to restore these states before the first A is selected. */
#define main journal_frontend_main
#include "../../../../siqs/mpu-siqs.c"
#undef main
#include "../../../../siqs/siqs.c"

int main(int argc, char **argv) {
  const uint64_t la = UINT64_C(0x1a6b50f924d38ce7);
  const uint64_t cof = UINT64_C(0x8e721c06a9b5d43f);
  const uint64_t choice = UINT64_C(0xc3f9821d4076ab5e);
  siqs_journal_t *journal;
  siqs_factor_array_t result;
  siqs_ctx_t ctx;
  siqs_poly_t poly;
  mpz_t n, q;
  uint64_t families = 0, polynomials = 0;
  uint32_t i;
  int reading;
  if (argc != 2) return 2;
  mpz_init_set_ui(n, 1); mpz_mul_2exp(n, n, 89); mpz_sub_ui(n, n, 1);
  mpz_init_set_ui(q, 1); mpz_mul_2exp(q, q, 107); mpz_sub_ui(q, q, 1);
  mpz_mul(n, n, q); mpz_clear(q);
  journal = gmp_siqs_journal_open(argv[1], n);
  reading = journal->reading;
  siqs_factor_array_init(&result, n); result.journal = journal;
  siqs_ctx_init(&ctx, n, n, &result, NULL, 0);
  if (!siqs_ctx_allocate(&ctx)) return 3;
  siqs_poly_init(&ctx, &poly);
  ctx.la_rng.state = reading ? 1U : la;
  ctx.cofactor_rng.state = reading ? 2U : cof;
  ctx.poly_rng.state = reading ? 3U : choice;
  if (sj_attach(&ctx, &poly, &families, &polynomials)) return 4;
  if (ctx.la_rng.state != la || ctx.cofactor_rng.state != cof ||
      ctx.poly_rng.state != choice || families != 0 || polynomials != 0)
    return 5;
  if (!siqs_choose_A(&ctx, &poly)) return 6;
  printf("Restored starting streams; first A indices:");
  for (i = 0; i < poly.q_count; i++) printf(" %u", poly.a_index[i]);
  putchar('\n');
  siqs_poly_clear(&ctx, &poly); siqs_ctx_clear(&ctx);
  gmp_siqs_journal_close(journal);
  gmp_siqs_free(result.values, result.count); free(result.primality);
  mpz_clear(n);
  return 0;
}
