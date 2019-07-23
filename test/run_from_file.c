#include "rw.h"
#include "scs.h"
#include "util.h"

scs_int override_setting(ScsSettings *s, char *param, char *val) {
  scs_printf("Attempting to override %s with value %s.\n", param, val);
  if (strcmp(param, "normalize") == 0) {
    s->normalize = atoi(val);
  } else if (strcmp(param, "scale") == 0) {
    s->scale = atof(val);
  } else if (strcmp(param, "rho_x") == 0) {
    s->rho_x = atof(val);
  } else if (strcmp(param, "max_iters") == 0) {
    s->max_iters = atoi(val);
  } else if (strcmp(param, "eps") == 0) {
    s->eps = atof(val);
  } else if (strcmp(param, "alpha") == 0) {
    s->alpha = atof(val);
  } else if (strcmp(param, "cg_rate") == 0) {
    s->cg_rate = atof(val);
  } else if (strcmp(param, "verbose") == 0) {
    s->verbose = atoi(val);
  } else if (strcmp(param, "acceleration_lookback") == 0) {
    s->acceleration_lookback = atoi(val);
  } else {
    return -1;
  }
  scs_printf("Success.\n");
  return 0;
}

/* Simple helper function to run problems from data files */
/* Mostly useful for debugging */
int main(int argc, char **argv) {
  char *filename;
  scs_int read_status;
  ScsData *d;
  ScsCone *k;
  ScsSolution *sol;
  ScsInfo info = {0};
  scs_int i;
  if (argc < 2) {
    scs_printf("Need to specify a filename, exit.\n");
    return -1;
  }
  filename = argv[1];
  read_status = SCS(read_data)(filename, &d, &k);
  if (read_status < 0) {
    scs_printf("Data read failure, exit.\n");
    return -1;
  }
  for (i = 2; i < argc; i += 2) {
    if (argc < i + 2) {
      scs_printf("Incorrect number of arguments supplied\n.");
      SCS(free_data)(d, k);
      return -1;
    }
    if (override_setting(d->stgs, argv[i], argv[i + 1]) < 0) {
      scs_printf("Unrecognized setting %s\n", argv[i]);
      SCS(free_data)(d, k);
      return -1;
    }
  }
  sol = scs_calloc(1, sizeof(ScsSolution));
  scs(d, k, sol, &info);
  SCS(free_data)(d, k);
  SCS(free_sol)(sol);
  return 0;
}
