#include "glbopts.h"

const char *SCS(version)(void) {
  return SCS_VERSION;
}
size_t SCS(sizeof_int)(void) {
  return sizeof(scs_int);
}
size_t SCS(sizeof_float)(void) {
  return sizeof(scs_float);
}
