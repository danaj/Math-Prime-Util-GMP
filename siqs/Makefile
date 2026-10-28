# Standalone SIQS; independent of Perl's configured compiler and linker.
# Run here, or use: make -C siqs
CC ?= cc
CPPFLAGS ?=
# 32-byte loop alignment avoids placement-sensitive GCC/x86 sieve slowdowns.
CFLAGS ?= -O3 -march=native -falign-loops=32
LDFLAGS ?=
LDLIBS ?= -lgmp -lm
PERL ?= perl
PYTHON ?= python3
EXTRA_FLAGS ?=
SIQS_CHECK_ARGS ?=
SIQS_TUNE_ARGS ?=

SOURCES = mpu-siqs.c siqs.c lanczos.c pr32_iter.c ../squfof126.c
HEADERS = siqs.h siqs_dep.h lanczos.h siqs_choose_inc.c siqs_sieve_kernel_inc.c siqs_journal_inc.c \
	psiqs_inc.c planczos_inc.c ../ptypes.h ../mont64.h pr32_iter.h \
	../squfof126.h ../misc_ui.h
BUILD = $(PERL) tools/build-siqs.pl --cc "$(CC)" \
	--cppflags "$(CPPFLAGS)" --cflags "$(CFLAGS) $(EXTRA_FLAGS)" \
	--ldflags "$(LDFLAGS)" --ldlibs "$(LDLIBS)"
BUILD_DEPS = tools/build-siqs.pl tools/pthread-probe.c
BENCH_DEPS = $(SOURCES) $(HEADERS)

.PHONY: all siqs serial siqs-serial tuned siqs-tuned siqs-tuned-8 \
	check check-siqs test-siqs check-pr32-iter check-siqs-autotune check-siqs-journal clean

all: msiqs
siqs: msiqs

# Probe by linking, never running. Standalone objects never enter the module build.
msiqs: $(SOURCES) $(HEADERS) $(BUILD_DEPS) Makefile
	@$(BUILD) --output $@ $(SOURCES)

# Explicit overrides rebuild even an otherwise up-to-date executable.
serial: siqs-serial
siqs-serial: $(SOURCES) $(HEADERS) $(BUILD_DEPS)
	@$(BUILD) --output msiqs --serial $(SOURCES)

tuned: siqs-tuned
siqs-tuned: siqs-sieve-bench tools/siqs-autotune.pl
	@siqs_tuned_flags=`$(PERL) tools/siqs-autotune.pl --flags $(SIQS_TUNE_ARGS)` && \
	$(BUILD) --cflags "$(CFLAGS) $(EXTRA_FLAGS) $$siqs_tuned_flags" --output msiqs $(SOURCES)

siqs-tuned-8:
	$(MAKE) siqs-tuned SIQS_TUNE_ARGS='--threads 8'

lanczos-bench: tools/lanczos-bench.c lanczos.c lanczos.h planczos_inc.c siqs_dep.h ../ptypes.h $(BUILD_DEPS) Makefile
	@$(BUILD) --ldlibs "-lm" --output $@ tools/lanczos-bench.c

siqs-la-bench: tools/siqs-la-bench.c $(BENCH_DEPS) $(BUILD_DEPS) Makefile
	@$(BUILD) --serial --output $@ tools/siqs-la-bench.c siqs.c pr32_iter.c ../squfof126.c

siqs-sieve-bench: tools/siqs-sieve-bench.c tools/siqs-autotune.inc.c $(BENCH_DEPS) $(BUILD_DEPS) Makefile
	@$(BUILD) --output $@ tools/siqs-sieve-bench.c lanczos.c pr32_iter.c ../squfof126.c

# Repository-only checks; not dependencies of the Perl module's test target.
# To force serial checks: make check-siqs SIQS_CHECK_ARGS=--serial
check: check-siqs
test-siqs: check-siqs
check-siqs:
	$(PERL) ../xt/dev/siqs/siqs-check.pl --cc "$(CC)" --cflags "$(CPPFLAGS) $(CFLAGS) $(EXTRA_FLAGS)" --ldflags "$(LDFLAGS) $(LDLIBS)" $(SIQS_CHECK_ARGS)

pr32-iter-check: ../xt/dev/siqs/pr32-iter-check.c pr32_iter.c pr32_iter.h Makefile
	$(CC) $(CPPFLAGS) $(CFLAGS) $(EXTRA_FLAGS) -o $@ ../xt/dev/siqs/pr32-iter-check.c $(LDFLAGS)

check-pr32-iter: pr32-iter-check
	@trap 'rm -f pr32-iter-check' 0; ./pr32-iter-check

check-siqs-journal:
	$(PYTHON) ../xt/dev/siqs/journal/run-all.py

siqs-autotune-check: ../xt/dev/siqs/siqs-autotune-check.c tools/siqs-sieve-bench.c tools/siqs-autotune.inc.c $(BENCH_DEPS) $(BUILD_DEPS) Makefile
	@$(BUILD) --output $@ ../xt/dev/siqs/siqs-autotune-check.c lanczos.c pr32_iter.c ../squfof126.c

check-siqs-autotune: siqs-autotune-check
	./siqs-autotune-check
	$(PERL) ../xt/dev/siqs/siqs-autotune.t

clean:
	rm -f msiqs lanczos-bench siqs-la-bench siqs-sieve-bench siqs-autotune-check pr32-iter-check
