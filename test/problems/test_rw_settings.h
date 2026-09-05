#include <stdint.h>
#include <string.h>

#include "cones.h"
#include "glbopts.h"
#include "minunit.h"
#include "rw.h"
#include "scs.h"
#include "util.h"

/* Round-trip nondefault settings through write_data/read_data, then
 * verify the reader is robust to files stamped by OTHER 3.3.x builds:
 * the settings-block layout is keyed off a numeric version comparison
 * (>= 3.3.0), not exact string equality, and the extension section
 * carries its own schema version. */

static const char *test_rw_settings(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsCone *k2 = SCS_NULL;
  ScsData *d2 = SCS_NULL;
  ScsSettings *stgs2 = SCS_NULL;
  const char *fname = "test_rw_settings_data";
  scs_int pass;

  /* tiny LP: 1 var, 1 nonneg row */
  scs_float Ax[1] = {1.0};
  scs_int Ai[1] = {0};
  scs_int Ap[2] = {0, 1};
  scs_float b[1] = {1.0};
  scs_float c[1] = {1.0};
  ScsMatrix A = {Ax, Ai, Ap, 1, 1};
  d->m = 1;
  d->n = 1;
  d->A = &A;
  d->b = b;
  d->c = c;
  k->l = 1;

  scs_set_default_settings(stgs);
  /* nondefaults across both layout regions and the extension section */
  stgs->adaptive_diag_scale = 0; /* v2 extension field */
  stgs->time_limit_secs = 123.5; /* v1 extension field */
  stgs->acceleration_type_1 = 0; /* post-3.3.0 stgs-block field */
  stgs->acceleration_regularization = 3e-9;
  stgs->eps_abs = 2.5e-7;
  stgs->write_data_filename = fname;

  SCS(write_data)(d, k, stgs);

  mu_assert("rw_settings: read failed",
            SCS(read_data)(fname, &d2, &k2, &stgs2) == 0);
  /* Binary serialization preserves the stored values exactly in either
   * precision; comparing against double literals would fail with SFLOAT. */
  pass = stgs2->adaptive_diag_scale == stgs->adaptive_diag_scale &&
         stgs2->time_limit_secs == stgs->time_limit_secs &&
         stgs2->acceleration_type_1 == stgs->acceleration_type_1 &&
         stgs2->acceleration_regularization == stgs->acceleration_regularization &&
         stgs2->eps_abs == stgs->eps_abs;
  mu_assert("rw_settings: nondefault round-trip lost a setting", pass);
  SCS(free_data)(d2);
  SCS(free_cone)(k2);
  scs_free(stgs2);
  d2 = SCS_NULL;
  k2 = SCS_NULL;
  stgs2 = SCS_NULL;

  /* Patch the embedded version string to an adjacent 3.3.x build (3.3.9,
   * padded) and re-read: everything must still parse, since the settings
   * layout is keyed off the numeric major.minor, not the exact string. */
  {
    FILE *f = fopen(fname, "r+b");
    long pos = 2 * (long)sizeof(uint32_t); /* int_sz, float_sz */
    uint32_t ver_len;
    char patched[16];
    mu_assert("rw_settings: reopen failed", f != SCS_NULL);
    fseek(f, pos, SEEK_SET);
    mu_assert("rw_settings: read ver len",
              fread(&ver_len, sizeof(uint32_t), 1, f) == 1);
    mu_assert("rw_settings: unexpected version length",
              ver_len >= 5 && ver_len < 16);
    /* C requires a positioning call between a read and a write on an
     * update stream; Windows enforces it. */
    fseek(f, 0, SEEK_CUR);
    memset(patched, '9', ver_len); /* "3.3.99...": an adjacent 3.3.x */
    patched[0] = '3';
    patched[1] = '.';
    patched[2] = '3';
    patched[3] = '.';
    mu_assert("rw_settings: patch write",
              fwrite(patched, 1, ver_len, f) == ver_len);
    fclose(f);
  }
  mu_assert("rw_settings: adjacent-version read failed",
            SCS(read_data)(fname, &d2, &k2, &stgs2) == 0);
  pass = stgs2->adaptive_diag_scale == stgs->adaptive_diag_scale &&
         stgs2->time_limit_secs == stgs->time_limit_secs &&
         stgs2->acceleration_type_1 == stgs->acceleration_type_1 &&
         stgs2->acceleration_regularization == stgs->acceleration_regularization &&
         stgs2->eps_abs == stgs->eps_abs;
  mu_assert("rw_settings: adjacent-version read lost a setting", pass);
  SCS(free_data)(d2);
  SCS(free_cone)(k2);
  scs_free(stgs2);

  remove(fname);
  scs_free(stgs);
  scs_free(d);
  scs_free(k);
  return 0;
}
