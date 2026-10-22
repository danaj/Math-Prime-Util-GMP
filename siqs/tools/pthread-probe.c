/* Compile/link probe only; the standalone build does not run this program. */
#include <pthread.h>

static void *probe_worker(void *argument) {
  return argument;
}

int main(void) {
  pthread_t thread;
  if (pthread_create(&thread, 0, probe_worker, 0) != 0)
    return 1;
  return pthread_join(thread, 0) != 0;
}
