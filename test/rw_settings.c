#include "problems/test_rw_settings.h"

int main(void) {
  const char *result = test_rw_settings();
  puts(result ? result : "rw_settings passed");
  return result != SCS_NULL;
}
