/*
 * Build from the FLINT build directory:
 *   gcc -O3 -I. fqsieve.c -o fqsieve -L. -lflint -lgmp -lmpfr
 *
 * With no N argument, read decimal inputs from stdin and emit one flushed
 * "N: F1 F2 ..." line per input (exponents expanded). Blank lines are ignored.
 * The benchmark queue supplies exactly one input per fresh invocation. FLINT's
 * global thread pool and caches are cleaned up at process exit.
 * Use -q with an N argument for the same machine-readable output.
 *
 * qsieve_factor requires a composite input that is not a perfect power.
 * Returned factors are not guaranteed to be prime or distinct.
 */
#define _POSIX_C_SOURCE 200809L
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include <flint/qsieve.h>

static void usage(FILE *out, const char *program)
{
    fprintf(out, "Usage: %s [-t threads] [-q] [-v] [N]\n"
                 "  -t threads  Positive thread count (default: 1)\n"
                 "  -q          Print only N: F1 F2 ...\n"
                 "  -v          Print default QS tuning parameters\n"
                 "  -h          Show this help\n"
                 "Without N, read stdin; one quiet result per input line.\n"
                 "Verbose batch diagnostics go to stderr.\n"
                 "N must be composite and not a perfect power.\n",
            program);
}

static int parse_threads(const char *text, int *threads)
{
    const char *p;
    char *end;
    long value;

    if (*text == '\0')
        return 0;
    for (p = text; *p; ++p)
        if (*p < '0' || *p > '9')
            return 0;

    errno = 0;
    value = strtol(text, &end, 10);
    if (errno == ERANGE || *end != '\0' || value < 1 || value > INT_MAX)
        return 0;
    *threads = (int) value;
    return 1;
}

static void print_parameters(FILE *out, const fmpz_t n)
{
    flint_bitcnt_t bits = fmpz_bits(n);

    flint_fprintf(out, "Input size: %wu bits\n", (ulong) bits);
#ifdef QS_TUNE_SIZE
    {
        size_t i = 0;
        size_t rows = sizeof(qsieve_tune) / sizeof(qsieve_tune[0]);

        while (i + 1 < rows && qsieve_tune[i + 1][0] <= bits)
            ++i;

        fputs("Default tuning parameters from the installed FLINT header:\n", out);
        flint_fprintf(out, "  Tuning row:                 %wu bits\n",
                     (ulong) qsieve_tune[i][0]);
        flint_fprintf(out, "  Knuth-Schroeppel primes:     %wu\n",
                     (ulong) qsieve_tune[i][1]);
        flint_fprintf(out, "  Factor-base entries:        %wu\n",
                     (ulong) qsieve_tune[i][2]);
        flint_fprintf(out, "  Small-prime setting:        %wu\n",
                     (ulong) qsieve_tune[i][3]);
        flint_fprintf(out, "  Sieve size:                 %wu\n",
                     (ulong) qsieve_tune[i][4]);
        flint_fprintf(out, "  Sieve threshold setting:    %wu\n",
                     (ulong) qsieve_tune[i][5]);
        fputs("These are table defaults, not live engine measurements.\n", out);
    }
#else
    fputs("This FLINT header does not expose its QS tuning table.\n", out);
#endif
}

static int parse_number(fmpz_t n, const char *text)
{
    const char *p;

    if (*text == '\0')
        return 0;
    for (p = text; *p; ++p)
        if (*p < '0' || *p > '9')
            return 0;
    return fmpz_set_str(n, text, 10) == 0 && fmpz_cmp_ui(n, 2) >= 0;
}

static int factor_one(const fmpz_t n, int quiet, int verbose)
{
    fmpz_factor_t factors;
    fmpz_t product;
    slong i;
    ulong j;
    int success;

    if (!quiet || verbose)
    {
        FILE *out = quiet ? stderr : stdout;
        fprintf(out, "Running FLINT qsieve with %d thread(s)...\n",
                flint_get_num_threads());
        fputs("Factoring n = ", out);
        fmpz_fprint(out, n);
        fputc('\n', out);
        if (verbose)
            print_parameters(out, n);
        fflush(out);
    }

    /* A fresh factor container avoids retaining factors from earlier inputs. */
    fmpz_factor_init(factors);
    fmpz_init(product);
    qsieve_factor(factors, n);
    fmpz_factor_expand(product, factors);
    success = factors->num > 0 && fmpz_equal(product, n) &&
              (factors->num > 1 || factors->exp[0] > 1);
    for (i = 0; success && i < factors->num; ++i)
        if (factors->exp[i] == 0 || fmpz_cmp_ui(factors->p + i, 2) < 0 ||
            fmpz_cmp(factors->p + i, n) >= 0)
            success = 0;
    if (!success)
    {
        fputs("FLINT qsieve returned no complete split for ", stderr);
        fmpz_fprint(stderr, n);
        fputc('\n', stderr);
    }
    else if (quiet)
    {
        fmpz_print(n);
        putchar(':');
        for (i = 0; i < factors->num; ++i)
            for (j = 0; j < factors->exp[i]; ++j)
            {
                putchar(' ');
                fmpz_print(factors->p + i);
            }
        putchar('\n');
    }
    else
    {
        puts("\nFactors returned:");
        fmpz_factor_print(factors);
        putchar('\n');
    }

    fmpz_clear(product);
    fmpz_factor_clear(factors);
    /* Flush after cleanup so the caller's timing includes result release. */
    if (fflush(stdout) != 0)
    {
        perror("fqsieve stdout");
        success = 0;
    }
    return success;
}

int main(int argc, char **argv)
{
    int option, threads = 1, verbose = 0, quiet = 0;
    int result = EXIT_SUCCESS;
    fmpz_t n;

    opterr = 0;
    while ((option = getopt(argc, argv, "t:qvh")) != -1)
    {
        switch (option)
        {
            case 't':
                if (!parse_threads(optarg, &threads))
                {
                    fprintf(stderr, "Invalid thread count: %s\n", optarg);
                    return EXIT_FAILURE;
                }
                break;
            case 'q':
                quiet = 1;
                break;
            case 'v':
                if (++verbose > 1)
                {
                    fputs("-vv is unavailable: qsieve_factor has no runtime "
                          "progress-reporting interface.\n", stderr);
                    return EXIT_FAILURE;
                }
                break;
            case 'h':
                usage(stdout, argv[0]);
                return EXIT_SUCCESS;
            default:
                usage(stderr, argv[0]);
                return EXIT_FAILURE;
        }
    }

    if (argc - optind > 1)
    {
        usage(stderr, argv[0]);
        return EXIT_FAILURE;
    }

    fmpz_init(n);
    flint_set_num_threads(threads);
    if (optind < argc)
    {
        if (!parse_number(n, argv[optind]))
        {
            fputs("N must be a decimal integer greater than 1.\n", stderr);
            result = EXIT_FAILURE;
        }
        else if (!factor_one(n, quiet, verbose))
            result = EXIT_FAILURE;
    }
    else
    {
        char *line = NULL;
        size_t capacity = 0, line_number = 0;
        ssize_t length;

        while ((length = getline(&line, &capacity, stdin)) >= 0)
        {
            char *text = line, *end;

            ++line_number;
            if (memchr(line, '\0', (size_t) length) != NULL)
            {
                fprintf(stderr, "NUL byte on stdin line %zu.\n", line_number);
                result = EXIT_FAILURE;
                break;
            }
            while (isspace((unsigned char) *text))
                ++text;
            end = text + strlen(text);
            while (end > text && isspace((unsigned char) end[-1]))
                *--end = '\0';
            if (*text == '\0')
                continue;
            if (!parse_number(n, text))
            {
                fprintf(stderr, "Invalid decimal integer on stdin line %zu.\n",
                        line_number);
                result = EXIT_FAILURE;
                break;
            }
            if (!factor_one(n, 1, verbose))
            {
                result = EXIT_FAILURE;
                break;
            }
        }
        if (result == EXIT_SUCCESS && !feof(stdin))
        {
            perror("fqsieve stdin");
            result = EXIT_FAILURE;
        }
        free(line);
    }
    fmpz_clear(n);
    flint_cleanup();
    return result;
}
