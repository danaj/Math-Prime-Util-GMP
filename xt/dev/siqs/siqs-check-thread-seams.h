/* Test-only interposition, included before siqs.c and lanczos.c.  Definitions
 * are guarded, but the macros intentionally activate on each inclusion.
 * Callers undefine them immediately after including production source.
 * Native operations inside this header are never interposed. */
#ifdef PSIQS
#ifndef SIQS_CHECK_THREAD_SEAMS_H
#define SIQS_CHECK_THREAD_SEAMS_H
#include <errno.h>
#include <pthread.h>
#ifndef _WIN32
#include <unistd.h>
#include <signal.h>
#include <sys/wait.h>
#endif

typedef struct {
  void *(*start)(void *);
  void *argument;
} check_thread_call_t;

static pthread_mutex_t seam_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_key_t seam_key;
static int seam_key_ready;
static int seam_alloc_after = -1, seam_mutex_after = -1, seam_cond_after = -1;
static int seam_thread_after = -1, seam_thread_alternate, seam_quiet, seam_exit_status;
static unsigned long seam_thread_calls, seam_created, seam_joined, seam_warnings;
static long seam_threads_live, seam_mutexes_live, seam_conditions_live;
static void *(*seam_controlled_start)(void *);
static void (*seam_enter)(void *);
static int (*seam_unlock_before)(pthread_mutex_t *, void *);
static void (*seam_unlock_after)(void *);
static void (*seam_signal)(pthread_cond_t *, void *);

static int seam_fail(int *after) {
  int fail = *after == 0;
  if (*after >= 0) (*after)--;
  return fail;
}

static void *check_thread_malloc(size_t size) {
  int fail;
  /* Disabled during ordinary jobs: don't introduce allocator locks that
   * could accidentally hide a production race from the sanitizer. */
  if (seam_alloc_after < 0) return malloc(size);
  pthread_mutex_lock(&seam_mutex);
  fail = seam_fail(&seam_alloc_after);
  pthread_mutex_unlock(&seam_mutex);
  return fail ? NULL : malloc(size);
}
static void *check_thread_calloc(size_t count, size_t size) {
  int fail;
  if (seam_alloc_after < 0) return calloc(count, size);
  pthread_mutex_lock(&seam_mutex);
  fail = seam_fail(&seam_alloc_after);
  pthread_mutex_unlock(&seam_mutex);
  return fail ? NULL : calloc(count, size);
}
static int check_thread_mutex_init(pthread_mutex_t *mutex,
                                    const pthread_mutexattr_t *attr) {
  int fail, status;
  pthread_mutex_lock(&seam_mutex); fail = seam_fail(&seam_mutex_after);
  pthread_mutex_unlock(&seam_mutex);
  status = fail ? EAGAIN : pthread_mutex_init(mutex, attr);
  pthread_mutex_lock(&seam_mutex);
  if (status == 0) seam_mutexes_live++;
  pthread_mutex_unlock(&seam_mutex);
  return status;
}
static int check_thread_cond_init(pthread_cond_t *cond,
                                   const pthread_condattr_t *attr) {
  int fail, status;
  pthread_mutex_lock(&seam_mutex); fail = seam_fail(&seam_cond_after);
  pthread_mutex_unlock(&seam_mutex);
  status = fail ? EAGAIN : pthread_cond_init(cond, attr);
  pthread_mutex_lock(&seam_mutex);
  if (status == 0) seam_conditions_live++;
  pthread_mutex_unlock(&seam_mutex);
  return status;
}
static int check_thread_mutex_destroy(pthread_mutex_t *mutex) {
  int status = pthread_mutex_destroy(mutex);
  pthread_mutex_lock(&seam_mutex);
  if (status == 0) seam_mutexes_live--;
  pthread_mutex_unlock(&seam_mutex);
  return status;
}
static int check_thread_cond_destroy(pthread_cond_t *cond) {
  int status = pthread_cond_destroy(cond);
  pthread_mutex_lock(&seam_mutex);
  if (status == 0) seam_conditions_live--;
  pthread_mutex_unlock(&seam_mutex);
  return status;
}
static void *check_thread_entry(void *argument) {
  check_thread_call_t *call = (check_thread_call_t *)argument;
  void *result;
  if (pthread_setspecific(seam_key, call) != 0) abort();
  if (seam_enter != NULL) seam_enter(call->argument);
  result = call->start(call->argument);
  if (pthread_setspecific(seam_key, NULL) != 0) abort();
  free(call);
  return result;
}
static int check_thread_create(pthread_t *thread, const pthread_attr_t *attr,
                                void *(*start)(void *), void *argument) {
  int fail, status;
  check_thread_call_t *call = NULL;
  pthread_mutex_lock(&seam_mutex);
  fail = seam_thread_after == 0;
  if (seam_thread_after > 0) seam_thread_after--;
  if (seam_thread_alternate && (seam_thread_calls & 1U)) fail = 1;
  seam_thread_calls++;
  pthread_mutex_unlock(&seam_mutex);
  if (fail) return EAGAIN;
  if (start == seam_controlled_start && seam_key_ready) {
    call = (check_thread_call_t *)malloc(sizeof(*call));
    if (call == NULL) return ENOMEM;
    call->start = start; call->argument = argument;
    status = pthread_create(thread, attr, check_thread_entry, call);
    if (status != 0) free(call);
  } else status = pthread_create(thread, attr, start, argument);
  pthread_mutex_lock(&seam_mutex);
  if (status == 0) { seam_created++; seam_threads_live++; }
  pthread_mutex_unlock(&seam_mutex);
  return status;
}
static int check_thread_join(pthread_t thread, void **value) {
  int status = pthread_join(thread, value);
  pthread_mutex_lock(&seam_mutex);
  if (status == 0) { seam_joined++; seam_threads_live--; }
  pthread_mutex_unlock(&seam_mutex);
  return status;
}
static int check_thread_mutex_unlock(pthread_mutex_t *mutex) {
  check_thread_call_t *call = seam_key_ready
      ? (check_thread_call_t *)pthread_getspecific(seam_key) : NULL;
  int pause = call != NULL && seam_unlock_before != NULL
            ? seam_unlock_before(mutex, call->argument) : 0;
  int status = pthread_mutex_unlock(mutex);
  if (status == 0 && pause && seam_unlock_after != NULL)
    seam_unlock_after(call->argument);
  return status;
}
static int check_thread_cond_signal(pthread_cond_t *cond) {
  check_thread_call_t *call = seam_key_ready
      ? (check_thread_call_t *)pthread_getspecific(seam_key) : NULL;
  if (seam_signal != NULL) seam_signal(cond, call != NULL ? call->argument : NULL);
  return pthread_cond_signal(cond);
}
static int check_thread_fprintf(FILE *stream, const char *format, ...) {
  int quiet, result;
  va_list args;
  pthread_mutex_lock(&seam_mutex);
  quiet = seam_quiet && stream == stderr &&
    (strncmp(format, "PSIQS: pthread_create:", 22) == 0 ||
     strncmp(format, "PSIQS: worker pool unavailable;", 30) == 0);
  if (quiet) seam_warnings++;
  pthread_mutex_unlock(&seam_mutex);
  if (quiet) return 0;
  va_start(args, format); result = vfprintf(stream, format, args); va_end(args);
  return result;
}
static void check_thread_exit(int status) {
#ifndef _WIN32
  /* Expected fatal exits are exercised in a child. Avoid inherited/leaked
   * child allocations being diagnosed by LSan's normal-exit handler. */
  if (status == seam_exit_status && status != 0) _exit(status);
#endif
  exit(status);
}
#endif /* definitions */

#define malloc check_thread_malloc
#define calloc check_thread_calloc
#define pthread_create check_thread_create
#define pthread_join check_thread_join
#define pthread_mutex_init check_thread_mutex_init
#define pthread_mutex_destroy check_thread_mutex_destroy
#define pthread_cond_init check_thread_cond_init
#define pthread_cond_destroy check_thread_cond_destroy
#define pthread_mutex_unlock check_thread_mutex_unlock
#define pthread_cond_signal check_thread_cond_signal
#define fprintf check_thread_fprintf
#define exit check_thread_exit
#endif /* PSIQS */
