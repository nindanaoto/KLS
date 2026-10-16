/* Optional diagnostic link object for Linux x86-64. Set KLS_CPU_SAMPLE_FILE
 * to collect process-CPU instruction samples and a companion .maps file.
 * Link with -ldl. This is not part of the production solver library.
 *
 * ITIMER_PROF samples include system CPU time and can coalesce under load.
 * Counts locate hot code; they are not wall-time or speedup measurements.
 */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <signal.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <ucontext.h>

#if !defined(__linux__) || !defined(__x86_64__)
#error "This diagnostic requires Linux x86-64 ucontext instruction pointers"
#endif
_Static_assert(ATOMIC_INT_LOCK_FREE == 2,
               "Signal-handler counters must always be lock-free");

enum { SAMPLE_CAPACITY = 1 << 20, SAMPLE_PERIOD_US = 1000 };
static uintptr_t sample_pc[SAMPLE_CAPACITY];
static atomic_uint sample_count, handlers_active, closing;
static char *output_path;
static uintptr_t executable_base;
static int enabled;

static void capture_sample(int signal_number, siginfo_t *info, void *context) {
  (void)signal_number;
  (void)info;
  atomic_fetch_add(&handlers_active, 1);
  if (!atomic_load(&closing)) {
    /* Saturation bounds the counter as well as the storage. Each successful
     * reservation has one writer; serialization with shutdown precedes I/O. */
    unsigned at = atomic_load(&sample_count);
    while (at < SAMPLE_CAPACITY &&
           !atomic_compare_exchange_weak(&sample_count, &at, at + 1)) {}
    if (at < SAMPLE_CAPACITY)
      sample_pc[at] = (uintptr_t)
        ((ucontext_t *)context)->uc_mcontext.gregs[REG_RIP];
  }
  atomic_fetch_sub(&handlers_active, 1);
}

__attribute__((constructor)) static void start_sampling(void) {
  const char *path = getenv("KLS_CPU_SAMPLE_FILE");
  if (!path || !*path) return;
  struct sigaction previous;
  struct itimerval previous_timer;
  sigset_t mask;
  if (sigaction(SIGPROF, NULL, &previous) ||
      getitimer(ITIMER_PROF, &previous_timer) ||
      sigprocmask(SIG_SETMASK, NULL, &mask) ||
      previous.sa_handler != SIG_DFL || sigismember(&mask, SIGPROF) ||
      previous_timer.it_value.tv_sec || previous_timer.it_value.tv_usec ||
      previous_timer.it_interval.tv_sec || previous_timer.it_interval.tv_usec) {
    fputs("cpu_sampler: existing SIGPROF policy; sampling disabled\n", stderr);
    return;
  }
  output_path = strdup(path);
  if (!output_path) return;
  Dl_info module;
  if (!dladdr((void *)(uintptr_t)&start_sampling, &module)) return;
  executable_base = (uintptr_t)module.dli_fbase;
  struct sigaction action = {0};
  action.sa_sigaction = capture_sample;
  action.sa_flags = SA_SIGINFO | SA_RESTART;
  sigemptyset(&action.sa_mask);
  if (sigaction(SIGPROF, &action, NULL)) return;
  struct itimerval timer = {
    .it_interval = {.tv_usec = SAMPLE_PERIOD_US},
    .it_value = {.tv_usec = SAMPLE_PERIOD_US}
  };
  if (setitimer(ITIMER_PROF, &timer, NULL)) {
    sigaction(SIGPROF, &previous, NULL);
    return;
  }
  enabled = 1;
}

__attribute__((destructor)) static void finish_sampling(void) {
  if (!enabled) { free(output_path); return; }
  atomic_store(&closing, 1);
  const struct itimerval stopped = {0};
  setitimer(ITIMER_PROF, &stopped, NULL);
  while (atomic_load(&handlers_active)) {}
  /* Leave the now-inert handler installed until process exit: a pending
   * SIGPROF must not acquire its terminating default disposition. */
  const unsigned count = atomic_load(&sample_count);
  FILE *out = fopen(output_path, "wx");
  if (out) {
    fprintf(out, "{\"version\":1,\"period_us\":%d,\"capacity\":%d,"
            "\"samples\":%u,\"saturated\":%s,\"executable_base\":\"0x%lx\"}\n",
            SAMPLE_PERIOD_US, SAMPLE_CAPACITY, count,
            count == SAMPLE_CAPACITY ? "true" : "false",
            (unsigned long)executable_base);
    for (unsigned i = 0; i < count; ++i)
      fprintf(out, "{\"pc\":\"0x%lx\"}\n", (unsigned long)sample_pc[i]);
    if (fclose(out)) fputs("cpu_sampler: sample output failed\n", stderr);
  } else fputs("cpu_sampler: cannot create sample output\n", stderr);
  const size_t length = strlen(output_path);
  char *maps_path = malloc(length + sizeof(".maps"));
  if (maps_path) {
    memcpy(maps_path, output_path, length);
    memcpy(maps_path + length, ".maps", sizeof(".maps"));
    FILE *maps = fopen("/proc/self/maps", "r");
    FILE *copy = fopen(maps_path, "wx");
    if (maps && copy) {
      char line[4096];
      while (fgets(line, sizeof(line), maps)) fputs(line, copy);
    }
    if (maps) fclose(maps);
    if (copy) fclose(copy);
    free(maps_path);
  }
  free(output_path);
}
