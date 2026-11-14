#ifndef KLS_ENVIRONMENT_H
#define KLS_ENVIRONMENT_H
/* Shared by KLS and its vendored KLU sources. */
#include <stdlib.h>
#include <string.h>

/* Diagnostic switches are KLS_* environment variables, and several are
   consulted on every numeric call.  A lookup walks the whole environment, so
   on a small system the lookups cost more than the factorization.  Almost no
   process sets any of them.  That case is recognised in constant time: once a
   scan has found no KLS_ name, a new variable can only appear as a new last
   entry (possibly in a moved array), so an unchanged array, length and last
   entry mean there is still none.  If
   any KLS_ variable exists, every lookup goes to getenv as before, so a switch
   set or cleared between calls still takes effect on the next call. */
extern char **environ;
static inline char *kls_getenv(const char *name) {
  static _Thread_local char **clean_environment;
  static _Thread_local size_t clean_count;
  static _Thread_local const char *clean_last;
  char **const entries = environ;
  if (entries == NULL) return NULL;
  if (entries == clean_environment && entries[clean_count] == NULL &&
      (clean_count == 0u || entries[clean_count - 1u] == clean_last)) return NULL;
  size_t count = 0u;
  for (; entries[count] != NULL; ++count)
    if (!strncmp(entries[count], "KLS_", 4u)) {
      clean_environment = NULL;
      return getenv(name);
    }
  clean_environment = entries;
  clean_count = count;
  clean_last = count ? entries[count - 1u] : NULL;
  return NULL;
}
#endif
