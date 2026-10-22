/* Lightweight Linux/macOS CLI adapter. No factor arithmetic or dependencies
 * beyond libc: the queue checks proper factors and their exact product.
 * Native stdout is copied to stderr; stdout contains only "N: factors".
 * Children stay in the queue's process group for interruption and timeouts.
 * The native style instead execs a normalized-stdin frontend directly.
 */
#define _POSIX_C_SOURCE 200809L
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <spawn.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

extern char **environ;
enum style { YAMAQUASI, RUSQSIEVE, MSIEVE, YAFU };
struct buffer { char *text; size_t used, allocated; };

static int error(const char *message)
{
  fprintf(stderr, "adapter: %s\n", message);
  return 1;
}

static int system_error(const char *message)
{
  fprintf(stderr, "adapter: %s: %s\n", message, strerror(errno));
  return 1;
}

static char *trim(char *s)
{
  char *end;
  while (isspace((unsigned char)*s)) s++;
  end = s + strlen(s);
  while (end > s && isspace((unsigned char)end[-1])) *--end = '\0';
  return s;
}

static int decimal(const char *s)
{
  if (*s < '1' || *s > '9') return 0;
  while (*++s) if (*s < '0' || *s > '9') return 0;
  return 1;
}

static char *read_number(void)
{
  char *line = NULL, *number = NULL;
  size_t capacity = 0;
  ssize_t length;
  while ((length = getline(&line, &capacity, stdin)) >= 0) {
    char *s;
    if (memchr(line, '\0', (size_t)length)) goto invalid;
    s = trim(line);
    if (!*s) continue;
    if (number || !decimal(s)) goto invalid;
    number = strdup(s);
    if (!number) { system_error("allocating input"); goto failed; }
  }
  if (ferror(stdin)) { system_error("reading input"); goto failed; }
  if (!number) goto invalid;
  free(line);
  return number;
invalid:
  error("expected exactly one positive decimal input");
failed:
  free(line);
  free(number);
  return NULL;
}

static int append(struct buffer *b, const char *s)
{
  size_t length = strlen(s), needed, size;
  char *p;
  if (length > SIZE_MAX - b->used - 1) {
    error("output size overflow");
    return 0;
  }
  needed = b->used + length + 1;
  if (needed > b->allocated) {
    size = b->allocated ? b->allocated : 256;
    while (size < needed) {
      if (size > SIZE_MAX / 2) { size = needed; break; }
      size *= 2;
    }
    p = realloc(b->text, size);
    if (!p) { system_error("allocating output"); return 0; }
    b->text = p;
    b->allocated = size;
  }
  memcpy(b->text + b->used, s, length + 1);
  b->used += length;
  return 1;
}

static char *factor_value(char *s, enum style style)
{
  char *p = s;
  if (style == YAMAQUASI || style == RUSQSIEVE)
    return decimal(s) ? s : NULL;
  if (style == MSIEVE) {
    if (!strncmp(p, "prp", 3)) p += 3;
    else if (*p == 'p') p++;
    else return NULL;
  } else {
    if (!strncmp(p, "PRP", 3)) p += 3;
    else if (*p == 'P') p++;
    else return NULL;
  }
  if (!isdigit((unsigned char)*p)) return NULL;
  while (isdigit((unsigned char)*p)) p++;
  if (style == YAFU) while (isspace((unsigned char)*p)) p++;
  if (*p != (style == MSIEVE ? ':' : '=')) return NULL;
  p++;
  while (isspace((unsigned char)*p)) p++;
  return decimal(p) ? p : NULL;
}

static int run(char *executable, char *threads, char *n, enum style style)
{
  FILE *input = NULL, *output = NULL;
  char *expression = NULL, *line = NULL;
  char *command[12] = { executable, NULL };
  struct buffer result = { NULL, 0, 0 };
  posix_spawn_file_actions_t actions;
  int rc = 1, err, status, parsed = 1;
  size_t capacity = 0;
  ssize_t length;
  pid_t pid, waited;

  /* Also enforced for direct adapter use, not just private queue workdirs. */
  if (style == MSIEVE && unlink("msieve.dat") && errno != ENOENT)
    return system_error("removing msieve.dat");
  if (style == YAFU) {
    if (unlink("siqs.dat") && errno != ENOENT)
      return system_error("removing siqs.dat");
    if (strlen(n) > SIZE_MAX - 7) return error("input size overflow");
    expression = malloc(strlen(n) + 7);
    if (!expression) return system_error("allocating expression");
    sprintf(expression, "siqs(%s)", n);
  }
  if (!(input = tmpfile()) || !(output = tmpfile())) {
    system_error("opening native I/O files");
    goto done;
  }
  switch (style) {
    case YAMAQUASI:
      command[1] = "--mode"; command[2] = "siqs";
      command[3] = "--verbose"; command[4] = "silent";
      command[5] = "--threads"; command[6] = threads; command[7] = n;
      command[8] = NULL;
      break;
    case RUSQSIEVE:
      command[1] = "--threads"; command[2] = threads;
      command[3] = "--progress"; command[4] = "never"; command[5] = NULL;
      fprintf(input, "%s\n", n);
      break;
    case MSIEVE:
      command[1] = "-q"; command[2] = "-t"; command[3] = threads;
      command[4] = "-s"; command[5] = "msieve.dat"; command[6] = n;
      command[7] = NULL;
      break;
    case YAFU:
      command[1] = expression; command[2] = "-threads"; command[3] = threads;
      command[4] = NULL;
      /* Piped stdin overrides YAFU's CLI expression; supply both. */
      fprintf(input, "%s\n", expression);
      break;
  }
  if (ferror(input) || fflush(input) || fseek(input, 0, SEEK_SET)) {
    system_error("preparing native input");
    goto done;
  }
  err = posix_spawn_file_actions_init(&actions);
  if (!err) {
    err = posix_spawn_file_actions_adddup2(&actions, fileno(input), STDIN_FILENO);
    if (!err) err = posix_spawn_file_actions_adddup2(&actions, fileno(output), STDOUT_FILENO);
    if (!err && fileno(input) > STDERR_FILENO)
      err = posix_spawn_file_actions_addclose(&actions, fileno(input));
    if (!err && fileno(output) > STDERR_FILENO)
      err = posix_spawn_file_actions_addclose(&actions, fileno(output));
    if (!err) err = posix_spawn(&pid, executable, &actions, NULL, command, environ);
    posix_spawn_file_actions_destroy(&actions);
  }
  if (err) {
    errno = err;
    system_error("launching executable");
    rc = 127;
    goto done;
  }
  do { waited = waitpid(pid, &status, 0); } while (waited < 0 && errno == EINTR);
  if (waited < 0) { system_error("waiting for executable"); goto done; }
  rc = WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
  if (fseek(output, 0, SEEK_SET)) { system_error("reading native output"); rc = 1; goto done; }
  if (rc == 0 && (!append(&result, n) || !append(&result, ":"))) parsed = 0;
  while ((length = getline(&line, &capacity, output)) >= 0) {
    char *s, *factor;
    if (fwrite(line, 1, (size_t)length, stderr) != (size_t)length) parsed = 0;
    if (rc || !parsed) continue;
    if (memchr(line, '\0', (size_t)length)) {
      error("NUL byte in native output"); parsed = 0; continue;
    }
    s = trim(line);
    if (!*s) continue;
    factor = factor_value(s, style);
    if (factor) {
      if (!append(&result, " ") || !append(&result, factor)) parsed = 0;
    } else if (style == YAMAQUASI || style == RUSQSIEVE) {
      error("unexpected output in factor-only stream"); parsed = 0;
    }
  }
  if (ferror(output) || fflush(stderr)) {
    system_error("copying native output"); parsed = 0;
  }
  if (!rc) {
    if (!parsed) rc = 1;
    else if (printf("%s\n", result.text) < 0 || fflush(stdout))
      rc = system_error("writing normalized output");
  }
done:
  if (input) fclose(input);
  if (output) fclose(output);
  free(expression);
  free(line);
  free(result.text);
  return rc;
}

int main(int argc, char **argv)
{
  char *executable = NULL, *threads = NULL, *style_name = NULL, *n, *end;
  enum style style;
  unsigned long thread_count;
  int i, rc;
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--help")) {
      puts("usage: cli --executable PATH --threads N --style native|yamaquasi|rusqsieve|msieve|yafu");
      return 0;
    }
    if (i + 1 == argc) return error("missing option value");
    if (!strcmp(argv[i], "--executable")) executable = argv[++i];
    else if (!strcmp(argv[i], "--threads")) threads = argv[++i];
    else if (!strcmp(argv[i], "--style")) style_name = argv[++i];
    else return error("unknown option");
  }
  if (!executable || !*executable || !threads || !style_name)
    return error("--executable, --threads and --style are required");
  errno = 0;
  thread_count = strtoul(threads, &end, 10);
  if (errno || *threads < '0' || *threads > '9' || *end ||
      !thread_count || thread_count > INT_MAX)
    return error("expected positive thread count within INT_MAX");
  if (!strcmp(style_name, "native")) {
    char *command[] = { executable, "-t", threads, NULL };
    /* MPU-SIQS and our FLINT frontend already implement the queue protocol.
     * Do not consume stdin, redirect output, or create another process. */
    execv(executable, command);
    system_error("launching native executable");
    return 127;
  }
  if (!strcmp(style_name, "yamaquasi")) style = YAMAQUASI;
  else if (!strcmp(style_name, "rusqsieve")) style = RUSQSIEVE;
  else if (!strcmp(style_name, "msieve")) style = MSIEVE;
  else if (!strcmp(style_name, "yafu")) style = YAFU;
  else return error("unknown style");
  n = read_number();
  if (!n) return 1;
  rc = run(executable, threads, n, style);
  free(n);
  return rc;
}
