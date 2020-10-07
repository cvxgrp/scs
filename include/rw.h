#ifndef RW_H_GUARD
#define RW_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "scs.h"
#include "util.h"

void SCS(write_data)(const ScsData *d, const ScsCone *k);
scs_int SCS(read_data)(const char * filename, ScsData **d, ScsCone **k);
void SCS(log_data_to_csv)(const ScsData *d, const ScsCone *k, const ScsWork *w,
                          const ScsResiduals * r, scs_int iter,
                          SCS(timer) * solve_timer);

#ifdef __cplusplus
}
#endif
#endif
