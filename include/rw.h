/*
 * Read/write utilities for serializing SCS problem data to binary files
 * and logging solve progress to CSV. Disabled when compiled with
 * NO_READ_WRITE=1.
 */

#ifndef RW_H_GUARD
#define RW_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "util.h"

void SCS(write_data)(const ScsData *d, const ScsCone *k,
                     const ScsSettings *stgs);
scs_int SCS(read_data)(const char *filename, ScsData **d, ScsCone **k,
                       ScsSettings **stgs);
scs_int SCS(open_csv_log_file)(ScsWork *w);
void SCS(close_csv_log_file)(ScsWork *w);
void SCS(log_data_to_csv)(const ScsCone *k, const ScsWork *w, scs_int iter,
                          SCS(timer) * solve_timer);

#ifdef __cplusplus
}
#endif
#endif
