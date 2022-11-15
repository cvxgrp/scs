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
  } else if (strcmp(param, "eps_abs") == 0) {
    s->eps_abs = atof(val);
  } else if (strcmp(param, "eps_rel") == 0) {
    s->eps_rel = atof(val);
  } else if (strcmp(param, "eps_infeas") == 0) {
    s->eps_infeas = atof(val);
  } else if (strcmp(param, "alpha") == 0) {
    s->alpha = atof(val);
  } else if (strcmp(param, "verbose") == 0) {
    s->verbose = atoi(val);
  } else if (strcmp(param, "acceleration_lookback") == 0) {
    s->acceleration_lookback = atoi(val);
  } else if (strcmp(param, "acceleration_interval") == 0) {
    s->acceleration_interval = atoi(val);
  } else if (strcmp(param, "adaptive_scale") == 0) {
    s->adaptive_scale = atoi(val);
  } else if (strcmp(param, "log_csv_filename") == 0) {
    s->log_csv_filename = val;
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
  ScsSettings *stgs;
  ScsSolution *sol;
  ScsInfo info = {0};
  scs_int i;
  if (argc < 2) {
    scs_printf("Need to specify a filename, exit.\n");
    return -1;
  }
  filename = argv[1];
  read_status = SCS(read_data)(filename, &d, &k, &stgs);
  if (read_status < 0) {
    scs_printf("Data read failure, exit.\n");
    return -1;
  }
  for (i = 2; i < argc; i += 2) {
    if (argc < i + 2) {
      scs_printf("Incorrect number of arguments supplied.\n");

      SCS(free_data)(d);
      SCS(free_cone)(k);
      scs_free(stgs);

      return -1;
    }
    if (override_setting(stgs, argv[i], argv[i + 1]) < 0) {
      scs_printf("Unrecognized setting %s\n", argv[i]);

      SCS(free_data)(d);
      SCS(free_cone)(k);
      scs_free(stgs);

      return -1;
    }
  }
  if (!stgs->verbose) {
    scs_printf(
        "File data set `verbose` to 0, SCS will not output information. Add "
        "`verbose 1` to call to override.\n");
  }
  scs_printf("Solving problem.\n");
  sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  scs(d, k, stgs, sol, &info);

  SCS(free_data)(d);
  SCS(free_cone)(k);
  SCS(free_sol)(sol);
  scs_free(stgs);

  return 0;
}
