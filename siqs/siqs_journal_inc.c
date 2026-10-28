/* Coordinator-owned relation journal. Included after the graph/chooser code.
 * See JOURNAL.md. Nothing here is called by sieve workers. */
#ifdef SIQS_NO_JOURNAL

/* Keep the front end and factor loop independent of build-mode conditionals.
 * In particular, close(NULL) is an ordinary no-op in a journal-free build. */
siqs_journal_t *gmp_siqs_journal_open(const char *path, mpz_srcptr n) {
  (void)path; (void)n;
  croak("SIQS journal support disabled in this build");
  return NULL;
}
void gmp_siqs_journal_input(const siqs_journal_t *j, mpz_t n) {
  (void)j; (void)n;
  croak("SIQS journal support disabled in this build");
}
void gmp_siqs_journal_close(siqs_journal_t *j) {
  if (j != NULL) croak("SIQS journal support disabled in this build");
}
void gmp_siqs_journal_finish(siqs_journal_t *j, mpz_t *values,
                            size_t count, int complete) {
  (void)values; (void)count; (void)complete;
  if (j != NULL) croak("SIQS journal support disabled in this build");
}
static void siqs_journal_relation(siqs_ctx_t *ctx, const siqs_raw_relation_t *r) {
  (void)ctx; (void)r;
}
static void siqs_journal_family(siqs_ctx_t *ctx, const siqs_poly_t *poly) {
  (void)ctx; (void)poly;
}
static void sj_flush(siqs_ctx_t *ctx) { (void)ctx; }
static void sj_maybe_flush(siqs_ctx_t *ctx) { (void)ctx; }
static void siqs_journal_full(siqs_ctx_t *ctx, uint32_t length) {
  (void)ctx; (void)length;
}
static int sj_attach(siqs_ctx_t *ctx, siqs_poly_t *poly, uint64_t *families,
                     uint64_t *polynomials) {
  (void)ctx; (void)poly; (void)families; (void)polynomials;
  return 0;
}
static void sj_begin_attempt(siqs_ctx_t *ctx) { (void)ctx; }
static void sj_end_attempt(siqs_ctx_t *ctx) { (void)ctx; }
static mpz_t *sj_cached_call(siqs_journal_t *j, const mpz_t n, uint32_t *count) {
  (void)j; (void)n; (void)count;
  return NULL;
}
static void sj_finish_call(siqs_journal_t *j, const mpz_t n, mpz_t *values,
                            uint32_t count) {
  (void)j; (void)n; (void)values; (void)count;
}

#else

#include <errno.h>
#include <fcntl.h>
#include <time.h>
#ifdef _WIN32
# include <windows.h>
# include <io.h>
# include <sys/stat.h>
# include <share.h>
#else
# include <sys/types.h>
# include <unistd.h>
#endif

#ifndef SIQS_JOURNAL_BATCH_BYTES
# define SIQS_JOURNAL_BATCH_BYTES 262144U
#endif
#define SJ_MAX_FRAME (64U * 1024U * 1024U)
#define SJ_HEADER 56U
#define SJ_MAGIC UINT32_C(0x4b4c424a)
#define SJ_VERSION 1U
enum { SJ_INPUT = 1, SJ_CONTEXT, SJ_BATCH, SJ_ATTEMPT_END, SJ_RETURN, SJ_JOB_END };
enum { SJ_RELATION = 1, SJ_FAMILY };
static const unsigned char sj_preamble[16] = {
  'M','P','U','S','I','Q','J',0,
  (unsigned char)SJ_VERSION, (unsigned char)(SJ_VERSION >> 8),
  (unsigned char)(SJ_VERSION >> 16), (unsigned char)(SJ_VERSION >> 24),
  4,3,2,1
};
typedef struct { unsigned char *data; size_t size, alloc; } sj_buffer;
typedef struct { const unsigned char *data; size_t left; } sj_reader;
struct siqs_journal_t {
  FILE *file;
  char *path;
  uint64_t sequence, append_offset, families;
  uint64_t run, mono_start;
  int have_run, writing_run, clock_started;
#ifdef _WIN32
  uint64_t clock_frequency;
#endif
  unsigned char *frame;
  uint32_t frame_type, frame_size, frame_count;
  int reading, have_frame, replaying, need_position;
  sj_buffer batch;
  uint32_t events;
  unsigned lp_bytes;
  mpz_t original;
};

static void sj_error(siqs_journal_t *j, const char *what) {
  croak("SIQS journal %s: %s", j->path, what);
}
static uint64_t sj_clock_us(siqs_journal_t *j, int monotonic) {
#ifdef _WIN32
  if (monotonic) {
    LARGE_INTEGER count, frequency;
    uint64_t seconds, remainder;
    if (j->clock_frequency == 0) {
      if (!QueryPerformanceFrequency(&frequency) || frequency.QuadPart <= 0)
        sj_error(j, "cannot read performance-counter frequency");
      j->clock_frequency = (uint64_t)frequency.QuadPart;
    }
    if (!QueryPerformanceCounter(&count) || count.QuadPart < 0)
      sj_error(j, "cannot read monotonic clock");
    seconds = (uint64_t)count.QuadPart / j->clock_frequency;
    remainder = (uint64_t)count.QuadPart % j->clock_frequency;
    if (seconds >= UINT64_MAX / UINT64_C(1000000) ||
        remainder > UINT64_MAX / UINT64_C(1000000))
      sj_error(j, "performance counter exceeds clock range");
    return seconds * UINT64_C(1000000) +
           remainder * UINT64_C(1000000) / j->clock_frequency;
  } else {
    FILETIME time;
    uint64_t ticks;
    GetSystemTimeAsFileTime(&time);
    ticks = ((uint64_t)time.dwHighDateTime << 32) | time.dwLowDateTime;
    if (ticks < UINT64_C(116444736000000000)) sj_error(j, "UTC clock predates Unix epoch");
    return (ticks - UINT64_C(116444736000000000)) / 10U;
  }
#else
  struct timespec time;
  if (clock_gettime(monotonic ? CLOCK_MONOTONIC : CLOCK_REALTIME, &time) != 0 ||
      time.tv_sec < 0 || (uint64_t)time.tv_sec >= UINT64_MAX / UINT64_C(1000000))
    sj_error(j, "cannot read journal clock");
  return (uint64_t)time.tv_sec * UINT64_C(1000000) + (uint64_t)time.tv_nsec / 1000U;
#endif
}
static void sj_lock(siqs_journal_t *j) {
#ifdef _WIN32
  /* The stream is opened with _SH_DENYWR: another writer is excluded while
   * read-only inspectors can still read every byte of the committed prefix. */
  (void)j;
#else
  struct flock lock;
  memset(&lock, 0, sizeof(lock));
  lock.l_type = F_WRLCK; lock.l_whence = SEEK_SET;
  if (fcntl(fileno(j->file), F_SETLK, &lock) != 0)
    sj_error(j, "cannot obtain exclusive journal lock");
#endif
}
static uint64_t sj_tell(siqs_journal_t *j) {
#ifdef _WIN32
  __int64 at = _ftelli64(j->file);
#else
  off_t at = ftello(j->file);
#endif
  if (at < 0) sj_error(j, "cannot obtain file position");
  return (uint64_t)at;
}
static void sj_position(siqs_journal_t *j) {
  if (!j->need_position) return;
  /* Discard only an incomplete final frame, never a complete bad checksum.
   * Seeking also supplies the required read/write stream transition. */
#ifdef _WIN32
  if (j->append_offset > INT64_MAX ||
      _chsize_s(_fileno(j->file), j->append_offset) != 0 ||
      _fseeki64(j->file, (__int64)j->append_offset, SEEK_SET) != 0)
#else
  if (j->append_offset > INT64_MAX ||
      ftruncate(fileno(j->file), (off_t)j->append_offset) != 0 ||
      fseeko(j->file, (off_t)j->append_offset, SEEK_SET) != 0)
#endif
    sj_error(j, "cannot position journal for append");
  j->need_position = 0;
}
static uint32_t sj_crc(uint32_t crc, const unsigned char *p, size_t n) {
  /* The immutable table permits independent concurrent journal sessions. */
  static const uint32_t table[16] = {
    0x00000000U,0x1db71064U,0x3b6e20c8U,0x26d930acU,
    0x76dc4190U,0x6b6b51f4U,0x4db26158U,0x5005713cU,
    0xedb88320U,0xf00f9344U,0xd6d6a3e8U,0xcb61b38cU,
    0x9b64c2b0U,0x86d3d2d4U,0xa00ae278U,0xbdbdf21cU
  };
  while (n-- != 0) {
    crc ^= *p++;
    crc = (crc >> 4) ^ table[crc & 15U];
    crc = (crc >> 4) ^ table[crc & 15U];
  }
  return crc;
}
static void sj_reserve(sj_buffer *b, size_t n) {
  size_t need, alloc;
  if (n > SJ_MAX_FRAME || b->size > SJ_MAX_FRAME - n)
    croak("SIQS journal: frame exceeds format size limit");
  need = b->size + n;
  if (need <= b->alloc) return;
  alloc = b->alloc ? b->alloc : 256U;
  while (alloc < need) alloc = alloc <= SJ_MAX_FRAME / 2 ? alloc * 2 : SJ_MAX_FRAME;
  b->data = (unsigned char *)siqs_realloc(b->data, alloc);
  b->alloc = alloc;
}
static void sj_bytes(sj_buffer *b, const void *p, size_t n) {
  sj_reserve(b, n);
  if (n != 0) memcpy(b->data + b->size, p, n);
  b->size += n;
}
static void sj_u32(sj_buffer *b, uint32_t x) {
  unsigned char p[4]; unsigned i;
  for (i = 0; i < 4; i++) p[i] = (unsigned char)(x >> (8U * i));
  sj_bytes(b, p, 4);
}
static void sj_u64(sj_buffer *b, uint64_t x) {
  unsigned char p[8]; unsigned i;
  for (i = 0; i < 8; i++) p[i] = (unsigned char)(x >> (8U * i));
  sj_bytes(b, p, 8);
}
static uint32_t sj_get32(sj_reader *r) {
  uint32_t x = 0; unsigned i;
  if (r->left < 4) croak("SIQS journal: short integer field");
  for (i = 0; i < 4; i++) x |= (uint32_t)r->data[i] << (8U * i);
  r->data += 4; r->left -= 4; return x;
}
static uint64_t sj_get64(sj_reader *r) {
  uint64_t x = 0; unsigned i;
  if (r->left < 8) croak("SIQS journal: short integer field");
  for (i = 0; i < 8; i++) x |= (uint64_t)r->data[i] << (8U * i);
  r->data += 8; r->left -= 8; return x;
}
/* Callers reserve the complete record once before using these encoders. */
static unsigned char *sj_put_uint(unsigned char *p, uint64_t x, unsigned bytes) {
  unsigned i;
  for (i = 0; i < bytes; i++) *p++ = (unsigned char)(x >> (8U * i));
  return p;
}
static uint64_t sj_get_uint(sj_reader *r, unsigned bytes) {
  uint64_t x = 0; unsigned i;
  if (r->left < bytes) croak("SIQS journal: short large-prime label");
  for (i = 0; i < bytes; i++) x |= (uint64_t)r->data[i] << (8U * i);
  r->data += bytes; r->left -= bytes; return x;
}
static unsigned char *sj_put_varint(unsigned char *p, uint32_t x) {
  while (x >= 128U) {
    *p++ = (unsigned char)((x & 127U) | 128U);
    x >>= 7;
  }
  *p++ = (unsigned char)x;
  return p;
}
static uint32_t sj_get_varint(sj_reader *r) {
  uint32_t x = 0; unsigned shift;
  for (shift = 0; shift <= 28U; shift += 7U) {
    unsigned byte;
    if (r->left == 0) croak("SIQS journal: short ULEB128 value");
    byte = *r->data++; r->left--;
    if (shift == 28U && byte > 15U)
      croak("SIQS journal: ULEB128 overflow");
    x |= (uint32_t)(byte & 127U) << shift;
    if ((byte & 128U) == 0) {
      if (shift != 0 && byte == 0) croak("SIQS journal: noncanonical ULEB128 value");
      return x;
    }
  }
  croak("SIQS journal: invalid ULEB128 value");
  return 0;
}
static void sj_mpz(sj_buffer *b, const mpz_t n) {
  size_t size = mpz_sgn(n) ? (mpz_sizeinbase(n, 2) + 7U) / 8U : 0, written = 0;
  if (mpz_sgn(n) < 0 || size > SJ_MAX_FRAME - 4U)
    croak("SIQS journal: invalid magnitude");
  sj_u32(b, (uint32_t)size); sj_reserve(b, size);
  if (size != 0) mpz_export(b->data + b->size, &written, 1, 1, 1, 0, n);
  if (size != written) croak("SIQS journal: magnitude size mismatch");
  b->size += size;
}
static void sj_get_mpz(sj_reader *r, mpz_t n, size_t maximum) {
  uint32_t size = sj_get32(r);
  if (size > maximum || size > r->left || (size != 0 && r->data[0] == 0))
    croak("SIQS journal: invalid magnitude length or encoding");
  mpz_import(n, size, 1, 1, 1, 0, r->data);
  r->data += size; r->left -= size;
}
static void sj_done(sj_reader *r) {
  if (r->left != 0) croak("SIQS journal: unexpected trailing payload");
}
static void sj_consume(siqs_journal_t *j) {
  free(j->frame); j->frame = NULL; j->have_frame = 0;
}
static int sj_peek(siqs_journal_t *j) {
  unsigned char header[SJ_HEADER];
  sj_reader h;
  uint64_t start, sequence, run;
  uint32_t magic, header_crc, crc;
  size_t got; int incomplete = 0;
  if (j->have_frame) return 1;
  if (!j->reading) return 0;
  start = sj_tell(j);
  got = fread(header, 1, SJ_HEADER, j->file);
  incomplete = got != 0;
  if (ferror(j->file)) sj_error(j, "read failed");
  if (got != SJ_HEADER) goto tail;
  h.data = header; h.left = SJ_HEADER;
  magic = sj_get32(&h); j->frame_type = sj_get32(&h);
  sequence = sj_get64(&h); j->frame_size = sj_get32(&h);
  j->frame_count = sj_get32(&h); run = sj_get64(&h);
  (void)sj_get64(&h); (void)sj_get64(&h); /* Clocks are informational. */
  header_crc = sj_get32(&h); crc = sj_get32(&h);
  if ((sj_crc(UINT32_MAX, header, 48) ^ UINT32_MAX) != header_crc)
    sj_error(j, "frame header checksum mismatch");
  if (magic != SJ_MAGIC || sequence != j->sequence ||
      j->frame_size > SJ_MAX_FRAME || j->frame_type < SJ_INPUT || j->frame_type > SJ_JOB_END)
    sj_error(j, "invalid frame header, version, or sequence");
  j->frame = (unsigned char *)siqs_malloc(j->frame_size ? j->frame_size : 1U);
  got = fread(j->frame, 1, j->frame_size, j->file);
  if (ferror(j->file)) sj_error(j, "read failed");
  if (got != j->frame_size) { free(j->frame); j->frame = NULL; goto tail; }
  if ((sj_crc(sj_crc(UINT32_MAX, header, 52), j->frame, j->frame_size) ^ UINT32_MAX) != crc)
    sj_error(j, "checksum mismatch in a complete frame");
  if ((!j->have_run && run != 0) ||
      (j->have_run && run != j->run && (j->run == UINT64_MAX || run != j->run + 1)))
    sj_error(j, "invalid run counter sequence");
  j->run = run; j->have_run = 1;
  if (j->sequence == UINT64_MAX) sj_error(j, "frame sequence exhausted");
  j->sequence++; j->have_frame = 1; return 1;
tail:
  if (incomplete) fprintf(stderr, "# siqs journal: ignoring incomplete final frame in %s\n", j->path);
  clearerr(j->file); j->reading = 0; j->append_offset = start; j->need_position = 1;
  return 0;
}
static void sj_write(siqs_journal_t *j, uint32_t type, uint32_t count, const sj_buffer *b) {
  sj_buffer header = {NULL, 0, 0}; uint32_t crc;
  uint64_t now, utc, elapsed;
  if (j->reading || j->have_frame) sj_error(j, "attempted append before replay completed");
  if (j->sequence == UINT64_MAX) sj_error(j, "frame sequence exhausted");
  sj_position(j);
  now = sj_clock_us(j, 1);
  if (!j->clock_started) {
    j->mono_start = now; j->clock_started = 1;
  }
  if (!j->writing_run) {
    if (j->have_run && j->run == UINT64_MAX) sj_error(j, "run counter exhausted");
    j->run = j->have_run ? j->run + 1 : 0;
    j->have_run = 1; j->writing_run = 1;
  }
  utc = sj_clock_us(j, 0);
  elapsed = now >= j->mono_start ? now - j->mono_start : 0;
  sj_u32(&header, SJ_MAGIC); sj_u32(&header, type); sj_u64(&header, j->sequence);
  sj_u32(&header, (uint32_t)b->size); sj_u32(&header, count);
  sj_u64(&header, j->run); sj_u64(&header, utc); sj_u64(&header, elapsed);
  sj_u32(&header, sj_crc(UINT32_MAX, header.data, header.size) ^ UINT32_MAX);
  crc = sj_crc(sj_crc(UINT32_MAX, header.data, header.size), b->data, b->size) ^ UINT32_MAX;
  sj_u32(&header, crc);
  if (fwrite(header.data, 1, header.size, j->file) != header.size ||
      fwrite(b->data, 1, b->size, j->file) != b->size || fflush(j->file) != 0)
    sj_error(j, "write or flush failed");
  free(header.data); j->sequence++;
}
static sj_reader sj_frame_reader(const siqs_journal_t *j) {
  sj_reader r; r.data = j->frame; r.left = j->frame_size; return r;
}

/* The readable prefix summarizes the same identity as the binary stream.
 * Use canonical LF line endings even on Windows (the file is open binary). */
static void sj_expect_text(siqs_journal_t *j, const char *text) {
  while (*text != '\0')
    if (fgetc(j->file) != (unsigned char)*text++)
      sj_error(j, "invalid or incomplete ASCII header");
}
static void sj_read_ascii(siqs_journal_t *j, mpz_srcptr n) {
  char version[32];
  sj_buffer digits = {NULL, 0, 0};
  mpz_t stored;
  int ch;
  sj_expect_text(j, "MPUSIQJ\n");
  snprintf(version, sizeof(version), "Version: %u\n", (unsigned)SJ_VERSION);
  sj_expect_text(j, version);
  sj_expect_text(j, "Byte-order: little-endian\nN: ");
  while ((ch = fgetc(j->file)) != '\n') {
    unsigned char digit;
    if (ch < '0' || ch > '9') sj_error(j, "invalid or incomplete ASCII input number");
    digit = (unsigned char)ch;
    sj_bytes(&digits, &digit, 1);
  }
  if (digits.size == 0) sj_error(j, "missing ASCII input number");
  sj_bytes(&digits, "", 1);
  mpz_init(stored);
  if (mpz_set_str(stored, (const char *)digits.data, 10) != 0)
    sj_error(j, "invalid ASCII input number");
  free(digits.data);
  if (mpz_sgn(stored) <= 0) sj_error(j, "ASCII input must be positive");
  if (n != NULL && mpz_cmp(stored, n) != 0)
    sj_error(j, "input does not match ASCII journal header");
  mpz_set(j->original, stored);
  mpz_clear(stored);
  /* Optional notes are uninterpreted; only an exact -- line ends the header.
   * Scan without allocating a buffer for comments or long metadata lines. */
  {
    size_t bytes = 0;
    unsigned dashes = 0;
    int matches = 1;
    for (;;) {
      ch = fgetc(j->file);
      if (ch == EOF) sj_error(j, "incomplete ASCII header: missing -- line");
      if (++bytes > SJ_MAX_FRAME) sj_error(j, "ASCII header exceeds size limit");
      if (ch == '\n') {
        if (matches && dashes == 2U) break;
        dashes = 0; matches = 1;
      } else if (matches) {
        if (ch != '-' || dashes == 2U) matches = 0;
        else dashes++;
      }
    }
  }
}

siqs_journal_t *gmp_siqs_journal_open(const char *path, mpz_srcptr n) {
  siqs_journal_t *j; sj_buffer b = {NULL, 0, 0}; unsigned char preamble[16];
  if (path == NULL || *path == 0 || (n != NULL && mpz_sgn(n) <= 0))
    croak("SIQS journal: expected a path and an optional positive input");
  j = (siqs_journal_t *)siqs_calloc(1, sizeof(*j));
  j->path = (char *)siqs_malloc(strlen(path) + 1); strcpy(j->path, path);
  mpz_init(j->original);
  if (n != NULL) mpz_set(j->original, n);
#ifdef _WIN32
  j->file = _fsopen(path, "r+b", _SH_DENYWR);
#else
  j->file = fopen(path, "r+b");
#endif
  if (j->file != NULL) {
    mpz_t stored; sj_reader r;
    sj_lock(j);
    j->reading = 1;
    sj_read_ascii(j, n);
    if (fread(preamble, 1, 16, j->file) != 16 || memcmp(preamble, sj_preamble, 16) != 0)
      sj_error(j, "invalid or incomplete binary header");
    if (!sj_peek(j) || j->frame_type != SJ_INPUT || j->frame_count != 1)
      sj_error(j, "missing input header");
    mpz_init(stored); r = sj_frame_reader(j);
    sj_get_mpz(&r, stored, SJ_MAX_FRAME); sj_done(&r);
    if (mpz_cmp(stored, j->original) != 0)
      sj_error(j, "ASCII and binary input headers do not match");
    mpz_clear(stored); sj_consume(j);
  } else {
    int fd;
    if (errno != ENOENT) sj_error(j, "cannot open journal");
    if (n == NULL) sj_error(j, "input omitted and journal does not exist");
#ifdef _WIN32
    fd = _sopen(path, _O_CREAT | _O_EXCL | _O_RDWR | _O_BINARY,
                 _SH_DENYWR, _S_IREAD | _S_IWRITE);
    if (fd < 0) sj_error(j, "cannot exclusively create journal");
    j->file = _fdopen(fd, "w+b");
#else
    fd = open(path, O_CREAT | O_EXCL | O_RDWR, 0644);
    if (fd < 0) sj_error(j, "cannot exclusively create journal");
    j->file = fdopen(fd, "w+b");
#endif
    if (j->file == NULL) sj_error(j, "cannot create journal stream");
    sj_lock(j);
    if (gmp_fprintf(j->file, "MPUSIQJ\nVersion: %u\nByte-order: little-endian\nN: %Zd\n--\n",
                    (unsigned)SJ_VERSION, n) < 0)
      sj_error(j, "cannot write ASCII header");
    if (fwrite(sj_preamble, 1, 16, j->file) != 16) sj_error(j, "cannot write file header");
    sj_mpz(&b, n); sj_write(j, SJ_INPUT, 1, &b); free(b.data);
  }
  return j;
}
void gmp_siqs_journal_input(const siqs_journal_t *j, mpz_t n) {
  mpz_set(n, j->original);
}
void gmp_siqs_journal_close(siqs_journal_t *j) {
  if (j == NULL) return;
  if (j->batch.size != 0) sj_error(j, "unflushed relation batch at close");
  if (sj_peek(j)) sj_error(j, "unconsumed calls remain in journal");
  if (fclose(j->file) != 0) sj_error(j, "close failed");
  mpz_clear(j->original); free(j->path); free(j->batch.data); free(j);
}

/* Fields are encoded individually: never write C structs, pointers, padding,
 * native mpz limbs, or native integer byte order. Floating fields use binary64. */
#define SJ_PARAMETER_U(X) \
 X(bits) X(policy_first_bits) X(policy_last_bits) X(max_large_primes) \
 X(fb_size) X(half_interval) X(q_count) X(poly_d) X(target_relations) \
 X(sieve_start) X(sieve_free_units) X(sieve_start_prime_floor) X(fb_floor) \
 X(relation_extra) X(multiplier_refine_divisor) X(fixed_half_interval) \
 X(sieve_initial) X(stage1_bias) X(sieve_raw_scores) X(sieve_saturating) \
 X(sieve_score_bound) X(large_prime_bound) X(smooth_bound)
#define SJ_PARAMETER_D(X) \
 X(lp_multiplier) X(residual_multiplier) X(sieve_start_index_exponent) \
 X(fb_coefficient) X(interval_scale) X(sieve_hit_exponent) \
 X(sieve_hit_bound_nominal) X(sieve_hit_residual_multiplier) X(sieve_hit_bound) X(log_scale)
#define SJ_COUNTERS(X) \
 X(total_candidates) X(split_attempts) X(split_square) X(split_siqs) \
 X(split_squfof) X(split_rho) X(split_fail) X(split_rejected)
static void sj_context(siqs_ctx_t *ctx, sj_buffer *b) {
  uint32_t i; size_t name = strlen(ctx->params.policy_name);
  if (sizeof(double) != sizeof(uint64_t) || DBL_MANT_DIG != 53 || DBL_MAX_EXP != 1024)
    croak("SIQS journal: binary64 doubles required by this journal format");
  sj_mpz(b, ctx->n); sj_u64(b, ctx->multiplier);
  sj_u32(b, SIQS_PACKED_FACTOR_ROW_BITS); sj_u32(b, (uint32_t)name);
  sj_bytes(b, ctx->params.policy_name, name);
#define SJ_PUT_U(field) sj_u64(b, (uint64_t)ctx->params.field);
  SJ_PARAMETER_U(SJ_PUT_U)
#undef SJ_PUT_U
#define SJ_PUT_D(field) { uint64_t word; memcpy(&word, &ctx->params.field, 8); sj_u64(b, word); }
  SJ_PARAMETER_D(SJ_PUT_D)
#undef SJ_PUT_D
  sj_u32(b, ctx->params.fb_size);
  for (i = 0; i < ctx->params.fb_size; i++) {
    sj_u32(b, ctx->fb[i].p); sj_u32(b, ctx->fb[i].sqrt_kn);
    sj_u32(b, ctx->fb[i].logp);
  }
  /* Collection-start checkpoints must remain explicit if seeding changes.
   * Replay restores them before regenerating the first A selection. */
  sj_u64(b, ctx->la_rng.state); sj_u64(b, ctx->cofactor_rng.state);
  sj_u64(b, ctx->poly_rng.state);
}
static void sj_flush(siqs_ctx_t *ctx) {
  siqs_journal_t *j = ctx->journal; sj_buffer b = {NULL, 0, 0};
  if (j == NULL || j->replaying || j->batch.size == 0) return;
  sj_u64(&b, j->families); sj_u64(&b, *ctx->journal_poly_count);
#define SJ_PUT_COUNTER(field) sj_u64(&b, ctx->field);
  SJ_COUNTERS(SJ_PUT_COUNTER)
#undef SJ_PUT_COUNTER
  sj_u64(&b, ctx->la_rng.state); sj_u64(&b, ctx->cofactor_rng.state);
  sj_u64(&b, ctx->poly_rng.state);
  sj_u64(&b, ctx->full_count);
  sj_u64(&b, (uint64_t)ctx->raw_count + ctx->one_lp.count);
  sj_u64(&b, ctx->matrix_next_target);
  {
    unsigned i;
    for (i = 0; i < 8; i++) sj_u64(&b, ctx->journal_cycles[i]);
  }
  sj_bytes(&b, j->batch.data, j->batch.size);
  sj_write(j, SJ_BATCH, j->events, &b); free(b.data);
  j->batch.size = 0; j->events = 0;
}
static void sj_maybe_flush(siqs_ctx_t *ctx) {
  siqs_journal_t *j = ctx->journal;
  if (!j->replaying && (j->batch.size >= SIQS_JOURNAL_BATCH_BYTES || j->events >= 256U))
    sj_flush(ctx);
}
static void siqs_journal_full(siqs_ctx_t *ctx, uint32_t length) {
  unsigned bin;
  if (ctx->journal == NULL) return;
  if (length == 0) sj_error(ctx->journal, "zero cycle length");
  bin = length >= 8 ? 7U : length - 1U;
  if (ctx->journal_cycles[bin] == UINT64_MAX)
    sj_error(ctx->journal, "invalid or exhausted cycle-length counter");
  ctx->journal_cycles[bin]++;
}
static void siqs_journal_family(siqs_ctx_t *ctx, const siqs_poly_t *poly) {
  siqs_journal_t *j = ctx->journal; uint32_t i;
  unsigned char event[2] = { SJ_FAMILY, (unsigned char)poly->q_count };
  if (j->replaying) return;
  if (j->families == UINT64_MAX) sj_error(j, "family counter exhausted");
  sj_bytes(&j->batch, event, sizeof(event));
  for (i = 0; i < poly->q_count; i++) sj_u32(&j->batch, poly->a_index[i]);
  j->families++; j->events++;
  if (j->batch.size >= SIQS_JOURNAL_BATCH_BYTES || j->events >= 256U) sj_flush(ctx);
}
static void siqs_journal_relation(siqs_ctx_t *ctx, const siqs_raw_relation_t *r) {
  siqs_journal_t *j = ctx->journal;
  sj_buffer *b = &j->batch;
  uint32_t i, previous = 0;
  unsigned labels = r->lp2 == 1 ? 0U : r->lp1 == 1 ? 1U : 2U;
  size_t size, written = 0, fixed;
  unsigned char *p;
  if (j->replaying) return;
  if (r->lp1 == 0 || r->lp1 > r->lp2 || r->lp2 > ctx->params.large_prime_bound)
    sj_error(j, "invalid large-prime labels");
  size = mpz_sgn(r->y) ? (mpz_sizeinbase(r->y, 2) + 7U) / 8U : 0;
  if (mpz_sgn(r->y) < 0 || size > SJ_MAX_FRAME - 26U)
    sj_error(j, "invalid relation magnitude");
  fixed = 10U + labels * j->lp_bytes + size;
  /* A 21-bit row delta plus its flag takes at most four ULEB128 bytes;
   * an explicit 11-bit exponent takes at most two. No per-factor reserves. */
  if (r->nfactors > (SJ_MAX_FRAME - fixed) / 6U)
    sj_error(j, "factor vector exceeds frame size limit");
  sj_reserve(b, fixed + 6U * r->nfactors);
  p = b->data + b->size;
  *p++ = SJ_RELATION; *p++ = (unsigned char)labels;
  if (labels == 2U) p = sj_put_uint(p, r->lp1, j->lp_bytes);
  if (labels != 0) p = sj_put_uint(p, r->lp2, j->lp_bytes);
  p = sj_put_uint(p, (uint32_t)size, 4);
  if (size != 0) mpz_export(p, &written, 1, 1, 1, 0, r->y);
  if (written != size) sj_error(j, "magnitude size mismatch");
  p += size;
  p = sj_put_uint(p, r->nfactors, 4);
  for (i = 0; i < r->nfactors; i++) {
    uint32_t row = siqs_raw_factor_row(r, i);
    uint32_t exponent = siqs_raw_factor_exponent(r, i);
    (void)siqs_pack_raw_factor(row, exponent);
    if (exponent == 0 || (i != 0 && row <= previous))
      sj_error(j, "invalid factor vector");
    p = sj_put_varint(p, ((row - previous) << 1) | (exponent != 1U));
    if (exponent != 1U) p = sj_put_varint(p, exponent);
    previous = row;
  }
  b->size = (size_t)(p - b->data);
  j->events++;
}

static void sj_partition(sj_buffer *b, const mpz_t n, mpz_t *values, uint32_t count) {
  uint32_t i; sj_mpz(b, n); sj_u32(b, count);
  for (i = 0; i < count; i++) sj_mpz(b, values[i]);
}
static mpz_t *sj_get_partition(siqs_journal_t *j, sj_reader *r, const mpz_t expected,
                              uint32_t *count) {
  mpz_t n, product; mpz_t *values; uint32_t i;
  mpz_init(n); mpz_init_set_ui(product, 1);
  sj_get_mpz(r, n, SJ_MAX_FRAME);
  if (mpz_cmp(n, expected) != 0) sj_error(j, "cofactor call or attempt does not match");
  *count = sj_get32(r);
  /* Every factor is >1 except the permitted partition [1] of input 1. */
  if (*count == 0 || *count > mpz_sizeinbase(expected, 2) || *count > r->left / 5U)
    sj_error(j, "invalid partition count");
  values = (mpz_t *)siqs_malloc(siqs_array_size(*count, sizeof(mpz_t)));
  for (i = 0; i < *count; i++) {
    mpz_init(values[i]); sj_get_mpz(r, values[i], (mpz_sizeinbase(expected, 2) + 7U) / 8U);
    if (mpz_cmp_ui(values[i], 1) < 0 || (mpz_cmp_ui(values[i], 1) == 0 && mpz_cmp_ui(expected, 1) != 0))
      sj_error(j, "invalid partition factor");
    mpz_mul(product, product, values[i]);
    if (mpz_cmp(product, expected) > 0) sj_error(j, "partition product exceeds input");
  }
  if (mpz_cmp(product, expected) != 0) sj_error(j, "partition product mismatch");
  mpz_clear(n); mpz_clear(product); return values;
}
static void sj_apply_end(siqs_ctx_t *ctx) {
  siqs_journal_t *j = ctx->journal; sj_reader r = sj_frame_reader(j);
  uint32_t found = sj_get32(&r), count;
  mpz_t *values = sj_get_partition(j, &r, ctx->original_n, &count);
  if (found > 1 || j->frame_count != 1) sj_error(j, "invalid attempt status");
  sj_done(&r);
  /* Preserve the recorded partition order as well as its product. The front
   * end's pending-cofactor stack must make the same subsequent calls. */
  gmp_siqs_free(ctx->result->values, ctx->result->count);
  free(ctx->result->primality);
  ctx->result->values = values;
  ctx->result->count = ctx->result->alloc = count;
  ctx->result->primality = (uint8_t *)siqs_calloc(count, sizeof(uint8_t));
  ctx->factor_found = (int)found; sj_consume(j);
}
static void sj_replay_batch(siqs_ctx_t *ctx, siqs_poly_t *poly, uint64_t *families,
                             uint64_t *polynomials) {
  siqs_journal_t *j = ctx->journal; sj_reader r = sj_frame_reader(j);
  uint64_t expected_families = sj_get64(&r), polys = sj_get64(&r);
  uint32_t event;
  mpz_t y;
  /* Counter snapshots precede events, but are applied after graph replay. */
  uint64_t counters[8], la, cof, poly_rng, full, partials, target, cycles[8];
  unsigned k = 0;
#define SJ_GET_COUNTER(field) counters[k++] = sj_get64(&r);
  SJ_COUNTERS(SJ_GET_COUNTER)
#undef SJ_GET_COUNTER
  la = sj_get64(&r); cof = sj_get64(&r); poly_rng = sj_get64(&r);
  if (la == 0 || cof == 0 || poly_rng == 0) sj_error(j, "invalid zero RNG state");
  full = sj_get64(&r); partials = sj_get64(&r); target = sj_get64(&r);
  for (k = 0; k < 8; k++) cycles[k] = sj_get64(&r);
  mpz_init(y);
  if (polys < *polynomials || polys > siqs_polynomial_limit(&ctx->params) ||
      j->frame_count > r.left / 2U)
    sj_error(j, "invalid progress or event count");
  for (event = 0; event < j->frame_count; event++) {
    unsigned kind = (unsigned)sj_get_uint(&r, 1);
    unsigned info = (unsigned)sj_get_uint(&r, 1);
    uint32_t i;
    if (kind == SJ_FAMILY) {
      uint32_t q = info;
      if (q != poly->q_count || q > r.left / 4U) sj_error(j, "invalid family factor count");
      /* Advance the chooser itself, without B/roots, sieving, or factoring.
       * This restores its cursor/RNG and A history, including q=1/2 cases. */
      for (i = 0; i < q; i++)
        if (poly->a_index[i] < ctx->params.fb_size) ctx->fb_in_a[poly->a_index[i]] = 0;
      if (!siqs_choose_A(ctx, poly)) sj_error(j, "recorded A sequence cannot be regenerated");
      for (i = 0; i < q; i++)
        if (sj_get32(&r) != poly->a_index[i]) sj_error(j, "A chooser differs from recorded sequence");
      if (j->families == UINT64_MAX) sj_error(j, "family counter exhausted");
      j->families++;
    } else if (kind == SJ_RELATION) {
      unsigned labels = info;
      uint64_t lp1, lp2;
      uint32_t count; siqs_factor_t *factors; siqs_raw_relation_t *raw;
      if (labels > 4U) sj_error(j, "invalid large-prime count");
      if (labels > 2U) sj_error(j, "3LP/4LP relation replay is not supported");
      if (labels > ctx->params.max_large_primes)
        sj_error(j, "relation exceeds the context's large-prime count");
      lp1 = labels == 2U ? sj_get_uint(&r, j->lp_bytes) : 1;
      lp2 = labels != 0 ? sj_get_uint(&r, j->lp_bytes) : 1;
      sj_get_mpz(&r, y, (mpz_sizeinbase(ctx->n, 2) + 7U) / 8U);
      if (mpz_cmp(y, ctx->n) >= 0 || lp1 > lp2 || lp2 > ctx->params.large_prime_bound ||
          (labels == 1U && lp2 <= 1) || (labels == 2U && lp1 <= 1))
        sj_error(j, "invalid relation magnitude or large-prime labels");
      count = sj_get32(&r);
      if (count > ctx->params.fb_size + 1U || count > r.left)
        sj_error(j, "invalid factor count");
      factors = (siqs_factor_t *)siqs_malloc(siqs_array_size(count ? count : 1U, sizeof(*factors)));
      for (i = 0; i < count; i++) {
        uint32_t code = sj_get_varint(&r), delta = code >> 1;
        uint32_t previous = i == 0 ? 0 : factors[i-1].row;
        uint32_t exponent = (code & 1U) ? sj_get_varint(&r) : 1;
        if (delta > ctx->params.fb_size - previous || (i != 0 && delta == 0) ||
            exponent == 0 || exponent > SIQS_PACKED_FACTOR_EXP_MAX ||
            ((code & 1U) != 0 && exponent == 1))
          sj_error(j, "invalid factor vector");
        factors[i].row = previous + delta;
        factors[i].exponent = exponent;
      }
      raw = siqs_raw_relation_new(ctx, y, factors, count, lp1, lp2);
      free(factors); siqs_accept_raw_relation(ctx, raw);
    } else sj_error(j, "unknown relation-batch event");
  }
  mpz_clear(y); sj_done(&r);
  if (full != ctx->full_count || partials != (uint64_t)ctx->raw_count + ctx->one_lp.count ||
      target > UINT32_MAX || memcmp(cycles, ctx->journal_cycles, sizeof(cycles)) != 0)
    sj_error(j, "graph statistics do not match replayed relations");
  ctx->matrix_next_target = (uint32_t)target;
  if (j->families != expected_families) sj_error(j, "family progress does not match events");
  *families = j->families; *polynomials = polys; k = 0;
#define SJ_SET_COUNTER(field) ctx->field = counters[k++];
  SJ_COUNTERS(SJ_SET_COUNTER)
#undef SJ_SET_COUNTER
  ctx->la_rng.state = la; ctx->cofactor_rng.state = cof;
  ctx->poly_rng.state = poly_rng; sj_consume(j);
}
/* Return 1 only for a previously completed attempt; otherwise EOF has been
 * reached and the caller can continue collection and append new frames. */
static int sj_attach(siqs_ctx_t *ctx, siqs_poly_t *poly, uint64_t *families,
                      uint64_t *polynomials) {
  siqs_journal_t *j = ctx->journal; sj_buffer b = {NULL, 0, 0};
  int ended = 0, existing; uint64_t relations_before;
  if (j == NULL) return 0;
  j->lp_bytes = 1;
  {
    uint64_t bound = ctx->params.large_prime_bound;
    while ((bound >>= 8) != 0) j->lp_bytes++;
  }
  existing = j->reading;
  ctx->journal_poly_count = polynomials; j->families = 0;
  sj_context(ctx, &b);
  if (sj_peek(j)) {
    sj_reader rng;
    if (j->frame_type != SJ_CONTEXT || j->frame_count != 1 ||
        b.size != j->frame_size || memcmp(b.data, j->frame, b.size - 24U) != 0)
      sj_error(j, "resolved parameters, factor base, or attempt differs");
    rng.data = j->frame + b.size - 24U; rng.left = 24U;
    ctx->la_rng.state = sj_get64(&rng);
    ctx->cofactor_rng.state = sj_get64(&rng);
    ctx->poly_rng.state = sj_get64(&rng);
    if (ctx->la_rng.state == 0 || ctx->cofactor_rng.state == 0 ||
        ctx->poly_rng.state == 0) sj_error(j, "invalid zero RNG state");
    sj_consume(j);
  } else sj_write(j, SJ_CONTEXT, 1, &b);
  free(b.data); j->replaying = 1;
  relations_before = ctx->accepted_smooth + ctx->accepted_one_lp + ctx->accepted_two_lp;
  while (sj_peek(j)) {
    if (j->frame_type == SJ_BATCH) sj_replay_batch(ctx, poly, families, polynomials);
    else if (j->frame_type == SJ_ATTEMPT_END) { sj_apply_end(ctx); ended = 1; break; }
    else sj_error(j, "unexpected frame inside collection attempt");
  }
  j->replaying = 0;
  if (ctx->verbose > 0 && existing) {
    printf("# siqs journal replayed %llu accepted relations, %llu A families%s\n",
      (unsigned long long)(ctx->accepted_smooth + ctx->accepted_one_lp + ctx->accepted_two_lp - relations_before),
      (unsigned long long)*families, ended ? " (attempt complete)" : "; appending new work");
    fflush(stdout);
  }
  return ended;
}
static void sj_begin_attempt(siqs_ctx_t *ctx) {
  siqs_journal_t *j = ctx->journal;
  if (j == NULL || j->clock_started) return;
  /* Time new work after replay, including the first batch. This only starts
   * the run clock: it neither appends nor consumes a run number, and later
   * recovery attempts/cofactors retain the same time origin. */
  j->mono_start = sj_clock_us(j, 1); j->clock_started = 1;
}
static void sj_end_attempt(siqs_ctx_t *ctx) {
  siqs_journal_t *j = ctx->journal; sj_buffer b = {NULL, 0, 0};
  if (j == NULL) return;
  sj_flush(ctx); sj_u32(&b, (uint32_t)ctx->factor_found);
  sj_partition(&b, ctx->original_n, ctx->result->values, ctx->result->count);
  sj_write(j, SJ_ATTEMPT_END, 1, &b); free(b.data);
}
static mpz_t *sj_cached_call(siqs_journal_t *j, const mpz_t n, uint32_t *count) {
  sj_reader r; mpz_t *values;
  if (j == NULL || !sj_peek(j) || j->frame_type != SJ_RETURN) return NULL;
  if (j->frame_count != 1) sj_error(j, "invalid return frame count");
  r = sj_frame_reader(j); values = sj_get_partition(j, &r, n, count);
  sj_done(&r); sj_consume(j); return values;
}
static void sj_finish_call(siqs_journal_t *j, const mpz_t n, mpz_t *values, uint32_t count) {
  sj_buffer b = {NULL, 0, 0}; uint32_t old_count, i; mpz_t *old;
  if (j == NULL) return;
  if (sj_peek(j)) {
    old = sj_cached_call(j, n, &old_count);
    if (old == NULL) sj_error(j, "expected completed cofactor call");
    if (old_count != count) sj_error(j, "return partition differs from completed attempt");
    for (i = 0; i < count; i++)
      if (mpz_cmp(old[i], values[i]) != 0)
        sj_error(j, "return partition order differs from completed attempt");
    gmp_siqs_free(old, old_count);
  } else {
    sj_partition(&b, n, values, count); sj_write(j, SJ_RETURN, 1, &b); free(b.data);
  }
}
void gmp_siqs_journal_finish(siqs_journal_t *j, mpz_t *values,
                            size_t count, int complete) {
  sj_buffer b = {NULL, 0, 0};
  if (j == NULL) return;
  if (count > UINT32_MAX) sj_error(j, "too many final factors");
  if (j->batch.size != 0) sj_error(j, "unflushed relation batch at job end");
  sj_u32(&b, complete != 0);
  sj_partition(&b, j->original, values, (uint32_t)count);
  if (sj_peek(j)) {
    if (j->frame_type != SJ_JOB_END || j->frame_count != 1 ||
        j->frame_size != b.size || memcmp(j->frame, b.data, b.size) != 0)
      sj_error(j, "final job status or partition differs");
    sj_consume(j);
  } else sj_write(j, SJ_JOB_END, 1, &b);
  free(b.data);
}
#undef SJ_PARAMETER_U
#undef SJ_PARAMETER_D
#undef SJ_COUNTERS
#endif /* !SIQS_NO_JOURNAL */
