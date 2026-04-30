/*
 * Read/write routines for serializing SCS problem data to binary files
 * and logging solve progress to CSV. Compiled as no-ops when NO_READ_WRITE=1.
 */

#include "rw.h"

#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cones.h"
#include "linalg.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

#if NO_READ_WRITE > 0 /* Disables all read / write functionality */

void SCS(write_data)(const ScsData *d, const ScsCone *k,
                     const ScsSettings *stgs) {
  /* Do nothing */
}
scs_int SCS(read_data)(const char *filename, ScsData **d, ScsCone **k,
                       ScsSettings **stgs) {
  /* Failure */
  return -1;
}
void SCS(log_data_to_csv)(const ScsCone *k, const ScsSettings *stgs,
                          const ScsWork *w, scs_int iter,
                          SCS(timer) * solve_timer) {
  /* Do nothing */
}

#else

/* writes/reads problem data to/from filename */
/* This is a VERY naive implementation, doesn't care about portability etc */

#define RW_EXT_MAGIC ((uint32_t)0x53435345u) /* "SCSE" */
#define RW_EXT_VERSION ((uint32_t)1u)

static scs_int checked_fread(void *ptr, size_t size, size_t nmemb, FILE *fin) {
  size_t ret;
  if (nmemb == 0) {
    return 0;
  }
  ret = fread(ptr, size, nmemb, fin);
  if (ret != nmemb) {
    scs_printf("Error: fread expected %lu items, got %lu\n",
               (unsigned long)nmemb, (unsigned long)ret);
    return -1;
  }
  return 0;
}

static scs_int checked_fwrite(const void *ptr, size_t size, size_t nmemb,
                              FILE *fout) {
  size_t ret;
  if (nmemb == 0) {
    return 0;
  }
  if (!ptr) {
    scs_printf("Error: fwrite passed NULL data for %lu items\n",
               (unsigned long)nmemb);
    return -1;
  }
  ret = fwrite(ptr, size, nmemb, fout);
  if (ret != nmemb) {
    scs_printf("Error: fwrite expected %lu items, wrote %lu\n",
               (unsigned long)nmemb, (unsigned long)ret);
    return -1;
  }
  return 0;
}

static scs_int write_int_array(const scs_int *x, scs_int n, FILE *fout) {
  if (n < 0) {
    scs_printf("Error: cannot write negative array length %li\n", (long)n);
    return -1;
  }
  return checked_fwrite(x, sizeof(scs_int), (size_t)n, fout);
}

static scs_int write_float_array(const scs_float *x, scs_int n, FILE *fout) {
  if (n < 0) {
    scs_printf("Error: cannot write negative array length %li\n", (long)n);
    return -1;
  }
  return checked_fwrite(x, sizeof(scs_float), (size_t)n, fout);
}

static scs_int write_scs_cone(const ScsCone *k, FILE *fout) {
  scs_int box_len = MAX(k->bsize - 1, 0);
  if (checked_fwrite(&(k->z), sizeof(scs_int), 1, fout) < 0 ||
      checked_fwrite(&(k->l), sizeof(scs_int), 1, fout) < 0 ||
      checked_fwrite(&(k->bsize), sizeof(scs_int), 1, fout) < 0 ||
      write_float_array(k->bl, box_len, fout) < 0 ||
      write_float_array(k->bu, box_len, fout) < 0 ||
      checked_fwrite(&(k->qsize), sizeof(scs_int), 1, fout) < 0 ||
      write_int_array(k->q, k->qsize, fout) < 0 ||
      checked_fwrite(&(k->ssize), sizeof(scs_int), 1, fout) < 0 ||
      write_int_array(k->s, k->ssize, fout) < 0 ||
      checked_fwrite(&(k->ep), sizeof(scs_int), 1, fout) < 0 ||
      checked_fwrite(&(k->ed), sizeof(scs_int), 1, fout) < 0 ||
      checked_fwrite(&(k->psize), sizeof(scs_int), 1, fout) < 0 ||
      write_float_array(k->p, k->psize, fout) < 0) {
    return -1;
  }
  return 0;
}

/*
 * Read integer data from file. If the integer width on file is
 * different to scs_int then it will cast the ints after reading
 * to be compatible with the SCS data types.
 */
static scs_int read_int(scs_int *dest, size_t file_int_sz, size_t nitems,
                        FILE *fin) {
  size_t val;
  size_t i;
  if (nitems == 0) {
    return 0;
  }
  if (file_int_sz == sizeof(scs_int)) {
    return checked_fread(dest, sizeof(scs_int), nitems, fin);
  }
  void *ptr = scs_calloc(nitems, file_int_sz);
  if (!ptr)
    return -1;
  val = fread(ptr, file_int_sz, nitems, fin);
  if (val != nitems) {
    scs_printf("Error: fread expected %lu integer items, got %lu\n",
               (unsigned long)nitems, (unsigned long)val);
    scs_free(ptr);
    return -1;
  }
  switch (file_int_sz) {
  case 4:
    for (i = 0; i < nitems; ++i) {
      dest[i] = (scs_int)(((int *)ptr)[i]);
    }
    break;
  case 8:
    for (i = 0; i < nitems; ++i) {
      dest[i] = (scs_int)(((long long *)ptr)[i]);
    }
    break;
  default:
    scs_printf("Error: unsupported integer size %lu in data file\n",
               (unsigned long)file_int_sz);
    scs_free(ptr);
    return -1;
  }
  scs_free(ptr);
  return 0;
}

static scs_int read_int_array(scs_int **dest, scs_int n, size_t file_int_sz,
                              FILE *fin) {
  *dest = SCS_NULL;
  if (n < 0) {
    scs_printf("Error: negative array length %li in data file\n", (long)n);
    return -1;
  }
  if (n == 0) {
    return 0;
  }
  *dest = (scs_int *)scs_calloc(n, sizeof(scs_int));
  if (!*dest) {
    return -1;
  }
  if (read_int(*dest, file_int_sz, (size_t)n, fin) < 0) {
    scs_free(*dest);
    *dest = SCS_NULL;
    return -1;
  }
  return 0;
}

static scs_int read_float_array(scs_float **dest, scs_int n, FILE *fin) {
  *dest = SCS_NULL;
  if (n < 0) {
    scs_printf("Error: negative array length %li in data file\n", (long)n);
    return -1;
  }
  if (n == 0) {
    return 0;
  }
  *dest = (scs_float *)scs_calloc(n, sizeof(scs_float));
  if (!*dest) {
    return -1;
  }
  if (checked_fread(*dest, sizeof(scs_float), (size_t)n, fin) < 0) {
    scs_free(*dest);
    *dest = SCS_NULL;
    return -1;
  }
  return 0;
}

static scs_int skip_bytes(FILE *fin, uint64_t bytes) {
  unsigned char buf[4096];
  while (bytes > 0) {
    size_t chunk = bytes < sizeof(buf) ? (size_t)bytes : sizeof(buf);
    if (checked_fread(buf, 1, chunk, fin) < 0) {
      return -1;
    }
    bytes -= chunk;
  }
  return 0;
}

static scs_int skip_int_array(FILE *fin, size_t file_int_sz, scs_int n) {
  if (n < 0) {
    scs_printf("Error: negative array length %li in data file\n", (long)n);
    return -1;
  }
  return skip_bytes(fin, (uint64_t)n * (uint64_t)file_int_sz);
}

static ScsCone *free_cone_null(ScsCone *k) {
  SCS(free_cone)(k);
  return SCS_NULL;
}

static ScsSettings *free_stgs_null(ScsSettings *s) {
  scs_free(s);
  return SCS_NULL;
}

static ScsMatrix *free_matrix_null(ScsMatrix *A) {
  SCS(free_scs_matrix)(A);
  return SCS_NULL;
}

static ScsData *free_data_null(ScsData *d) {
  SCS(free_data)(d);
  return SCS_NULL;
}

static scs_int read_data_cleanup(FILE *fin, ScsData **d, ScsCone **k,
                                 ScsSettings **stgs) {
  fclose(fin);
  SCS(free_data)(*d);
  SCS(free_cone)(*k);
  scs_free(*stgs);
  *d = SCS_NULL;
  *k = SCS_NULL;
  *stgs = SCS_NULL;
  return -1;
}

static ScsCone *read_scs_cone(FILE *fin, size_t file_int_sz) {
  scs_int box_len;
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  if (!k)
    return SCS_NULL;
  if (read_int(&(k->z), file_int_sz, 1, fin) < 0 ||
      read_int(&(k->l), file_int_sz, 1, fin) < 0 ||
      read_int(&(k->bsize), file_int_sz, 1, fin) < 0) {
    return free_cone_null(k);
  }
  if (k->bsize < 0) {
    scs_printf("Error: negative box cone size in data file\n");
    return free_cone_null(k);
  }
  box_len = MAX(k->bsize - 1, 0);
  if (read_float_array(&(k->bl), box_len, fin) < 0 ||
      read_float_array(&(k->bu), box_len, fin) < 0 ||
      read_int(&(k->qsize), file_int_sz, 1, fin) < 0 ||
      read_int_array(&(k->q), k->qsize, file_int_sz, fin) < 0 ||
      read_int(&(k->ssize), file_int_sz, 1, fin) < 0 ||
      read_int_array(&(k->s), k->ssize, file_int_sz, fin) < 0 ||
      read_int(&(k->ep), file_int_sz, 1, fin) < 0 ||
      read_int(&(k->ed), file_int_sz, 1, fin) < 0 ||
      read_int(&(k->psize), file_int_sz, 1, fin) < 0 ||
      read_float_array(&(k->p), k->psize, fin) < 0) {
    return free_cone_null(k);
  }
  return k;
}

static scs_int write_scs_stgs(const ScsSettings *s, FILE *fout) {
  /* Warm start to false for now */
  scs_int warm_start = 0;
  if (checked_fwrite(&(s->normalize), sizeof(scs_int), 1, fout) < 0 ||
      checked_fwrite(&(s->scale), sizeof(scs_float), 1, fout) < 0 ||
      checked_fwrite(&(s->rho_x), sizeof(scs_float), 1, fout) < 0 ||
      checked_fwrite(&(s->max_iters), sizeof(scs_int), 1, fout) < 0 ||
      checked_fwrite(&(s->eps_abs), sizeof(scs_float), 1, fout) < 0 ||
      checked_fwrite(&(s->eps_rel), sizeof(scs_float), 1, fout) < 0 ||
      checked_fwrite(&(s->eps_infeas), sizeof(scs_float), 1, fout) < 0 ||
      checked_fwrite(&(s->alpha), sizeof(scs_float), 1, fout) < 0 ||
      checked_fwrite(&(s->verbose), sizeof(scs_int), 1, fout) < 0 ||
      checked_fwrite(&warm_start, sizeof(scs_int), 1, fout) < 0 ||
      checked_fwrite(&(s->acceleration_lookback), sizeof(scs_int), 1, fout) <
          0 ||
      checked_fwrite(&(s->acceleration_interval), sizeof(scs_int), 1, fout) <
          0 ||
      checked_fwrite(&(s->acceleration_type_1), sizeof(scs_int), 1, fout) <
          0 ||
      checked_fwrite(&(s->acceleration_regularization), sizeof(scs_float), 1,
                     fout) < 0 ||
      checked_fwrite(&(s->acceleration_relaxation), sizeof(scs_float), 1,
                     fout) < 0 ||
      checked_fwrite(&(s->adaptive_scale), sizeof(scs_int), 1, fout) < 0) {
    return -1;
  }
  /* Do not write the write_data_filename */
  /* Do not write the log_csv_filename */
  return 0;
}

static ScsSettings *read_scs_stgs(FILE *fin, size_t file_int_sz,
                                  scs_int legacy_settings) {
  ScsSettings *s = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  if (!s)
    return SCS_NULL;
  scs_set_default_settings(s);
  if (read_int(&(s->normalize), file_int_sz, 1, fin) < 0 ||
      checked_fread(&(s->scale), sizeof(scs_float), 1, fin) < 0 ||
      checked_fread(&(s->rho_x), sizeof(scs_float), 1, fin) < 0 ||
      read_int(&(s->max_iters), file_int_sz, 1, fin) < 0 ||
      checked_fread(&(s->eps_abs), sizeof(scs_float), 1, fin) < 0 ||
      checked_fread(&(s->eps_rel), sizeof(scs_float), 1, fin) < 0 ||
      checked_fread(&(s->eps_infeas), sizeof(scs_float), 1, fin) < 0 ||
      checked_fread(&(s->alpha), sizeof(scs_float), 1, fin) < 0 ||
      read_int(&(s->verbose), file_int_sz, 1, fin) < 0 ||
      read_int(&(s->warm_start), file_int_sz, 1, fin) < 0 ||
      read_int(&(s->acceleration_lookback), file_int_sz, 1, fin) < 0 ||
      read_int(&(s->acceleration_interval), file_int_sz, 1, fin) < 0) {
    return free_stgs_null(s);
  }
  if (legacy_settings) {
    if (read_int(&(s->adaptive_scale), file_int_sz, 1, fin) < 0) {
      return free_stgs_null(s);
    }
  } else if (read_int(&(s->acceleration_type_1), file_int_sz, 1, fin) < 0 ||
             checked_fread(&(s->acceleration_regularization),
                           sizeof(scs_float), 1, fin) < 0 ||
             checked_fread(&(s->acceleration_relaxation), sizeof(scs_float), 1,
                           fin) < 0 ||
             read_int(&(s->adaptive_scale), file_int_sz, 1, fin) < 0) {
    return free_stgs_null(s);
  }
  return s;
}

static scs_int write_amatrix(const ScsMatrix *A, FILE *fout) {
  scs_int Anz = A->p[A->n];
  if (Anz < 0) {
    scs_printf("Error: matrix has negative nonzero count %li\n", (long)Anz);
    return -1;
  }
  if (checked_fwrite(&(A->m), sizeof(scs_int), 1, fout) < 0 ||
      checked_fwrite(&(A->n), sizeof(scs_int), 1, fout) < 0 ||
      write_int_array(A->p, A->n + 1, fout) < 0 ||
      write_float_array(A->x, Anz, fout) < 0 ||
      write_int_array(A->i, Anz, fout) < 0) {
    return -1;
  }
  return 0;
}

static ScsMatrix *read_amatrix(FILE *fin, size_t file_int_sz) {
  scs_int Anz;
  ScsMatrix *A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  if (!A)
    return SCS_NULL;
  if (read_int(&(A->m), file_int_sz, 1, fin) < 0 ||
      read_int(&(A->n), file_int_sz, 1, fin) < 0) {
    return free_matrix_null(A);
  }
  if (A->m < 0 || A->n < 0) {
    scs_printf("Error: negative matrix dimensions in data file\n");
    return free_matrix_null(A);
  }
  A->p = (scs_int *)scs_calloc(A->n + 1, sizeof(scs_int));
  if (!A->p || read_int(A->p, file_int_sz, A->n + 1, fin) < 0) {
    return free_matrix_null(A);
  }
  Anz = A->p[A->n];
  if (Anz < 0) {
    scs_printf("Error: negative matrix nonzero count in data file\n");
    return free_matrix_null(A);
  }
  A->x = (scs_float *)scs_calloc(Anz, sizeof(scs_float));
  A->i = (scs_int *)scs_calloc(Anz, sizeof(scs_int));
  if (!A->x || !A->i ||
      checked_fread(A->x, sizeof(scs_float), Anz, fin) < 0 ||
      read_int(A->i, file_int_sz, Anz, fin) < 0) {
    return free_matrix_null(A);
  }
  return A;
}

static scs_int write_scs_data(const ScsData *d, FILE *fout) {
  scs_int has_p = d->P ? 1 : 0;
  if (checked_fwrite(&(d->m), sizeof(scs_int), 1, fout) < 0 ||
      checked_fwrite(&(d->n), sizeof(scs_int), 1, fout) < 0 ||
      write_float_array(d->b, d->m, fout) < 0 ||
      write_float_array(d->c, d->n, fout) < 0 ||
      write_amatrix(d->A, fout) < 0) {
    return -1;
  }
  /* write has P bit */
  if (checked_fwrite(&has_p, sizeof(scs_int), 1, fout) < 0) {
    return -1;
  }
  if (d->P) {
    return write_amatrix(d->P, fout);
  }
  return 0;
}

static ScsData *read_scs_data(FILE *fin, size_t file_int_sz) {
  scs_int has_p = 0;
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  if (!d)
    return SCS_NULL;
  if (read_int(&(d->m), file_int_sz, 1, fin) < 0 ||
      read_int(&(d->n), file_int_sz, 1, fin) < 0) {
    return free_data_null(d);
  }
  if (d->m < 0 || d->n < 0) {
    scs_printf("Error: negative problem dimensions in data file\n");
    return free_data_null(d);
  }
  d->b = (scs_float *)scs_calloc(d->m, sizeof(scs_float));
  d->c = (scs_float *)scs_calloc(d->n, sizeof(scs_float));
  if (!d->b || !d->c ||
      checked_fread(d->b, sizeof(scs_float), d->m, fin) < 0 ||
      checked_fread(d->c, sizeof(scs_float), d->n, fin) < 0) {
    return free_data_null(d);
  }
  d->A = read_amatrix(fin, file_int_sz);
  if (!d->A) {
    return free_data_null(d);
  }
  /* If has_p bit is not set or this hits end of file then has_p = 0 */
  if (read_int(&has_p, file_int_sz, 1, fin) < 0) {
    return free_data_null(d);
  }
  d->P = has_p ? read_amatrix(fin, file_int_sz) : SCS_NULL;
  if (has_p && !d->P) {
    return free_data_null(d);
  }
  return d;
}

static scs_int write_ext_int_array(const scs_int *x, scs_int n, FILE *fout) {
  if (checked_fwrite(&n, sizeof(scs_int), 1, fout) < 0) {
    return -1;
  }
  return write_int_array(x, n, fout);
}

static scs_int write_scs_extensions(const ScsCone *k, const ScsSettings *s,
                                    FILE *fout) {
  scs_int zero = 0;
  uint32_t magic = RW_EXT_MAGIC;
  uint32_t version = RW_EXT_VERSION;
  if (checked_fwrite(&magic, sizeof(uint32_t), 1, fout) < 0 ||
      checked_fwrite(&version, sizeof(uint32_t), 1, fout) < 0 ||
      write_ext_int_array(k->cs, k->cssize, fout) < 0) {
    return -1;
  }
#ifdef USE_SPECTRAL_CONES
  if (write_ext_int_array(k->d, k->dsize, fout) < 0 ||
      checked_fwrite(&(k->nucsize), sizeof(scs_int), 1, fout) < 0 ||
      write_int_array(k->nuc_m, k->nucsize, fout) < 0 ||
      write_int_array(k->nuc_n, k->nucsize, fout) < 0 ||
      write_ext_int_array(k->ell1, k->ell1_size, fout) < 0 ||
      checked_fwrite(&(k->sl_size), sizeof(scs_int), 1, fout) < 0 ||
      write_int_array(k->sl_n, k->sl_size, fout) < 0 ||
      write_int_array(k->sl_k, k->sl_size, fout) < 0) {
    return -1;
  }
#else
  if (write_ext_int_array(SCS_NULL, zero, fout) < 0 ||
      checked_fwrite(&zero, sizeof(scs_int), 1, fout) < 0 ||
      write_int_array(SCS_NULL, zero, fout) < 0 ||
      write_int_array(SCS_NULL, zero, fout) < 0 ||
      write_ext_int_array(SCS_NULL, zero, fout) < 0 ||
      checked_fwrite(&zero, sizeof(scs_int), 1, fout) < 0 ||
      write_int_array(SCS_NULL, zero, fout) < 0 ||
      write_int_array(SCS_NULL, zero, fout) < 0) {
    return -1;
  }
#endif
  return checked_fwrite(&(s->time_limit_secs), sizeof(scs_float), 1, fout);
}

static scs_int read_ext_int_array(scs_int **dest, scs_int *n,
                                  size_t file_int_sz, FILE *fin) {
  if (read_int(n, file_int_sz, 1, fin) < 0) {
    return -1;
  }
  return read_int_array(dest, *n, file_int_sz, fin);
}

static scs_int read_scs_extensions(FILE *fin, size_t file_int_sz, ScsCone *k,
                                   ScsSettings *s) {
  unsigned char magic_buf[sizeof(uint32_t)];
  uint32_t magic, version;
  scs_int dsize, nucsize, ell1_size, sl_size;
  size_t ret = fread(magic_buf, 1, sizeof(magic_buf), fin);
  if (ret == 0 && feof(fin)) {
    return 0;
  }
  if (ret != sizeof(magic_buf)) {
    scs_printf("Error: incomplete read/write extension header\n");
    return -1;
  }
  memcpy(&magic, magic_buf, sizeof(magic));
  if (magic != RW_EXT_MAGIC) {
    scs_printf("Warning: ignoring unrecognized trailing data in SCS file\n");
    return 0;
  }
  if (checked_fread(&version, sizeof(uint32_t), 1, fin) < 0) {
    return -1;
  }
  if (version != RW_EXT_VERSION) {
    scs_printf("Error: unsupported SCS read/write extension version %lu\n",
               (unsigned long)version);
    return -1;
  }
  if (read_ext_int_array(&(k->cs), &(k->cssize), file_int_sz, fin) < 0 ||
      read_int(&dsize, file_int_sz, 1, fin) < 0) {
    return -1;
  }
#ifdef USE_SPECTRAL_CONES
  k->dsize = dsize;
  if (read_int_array(&(k->d), k->dsize, file_int_sz, fin) < 0 ||
      read_int(&nucsize, file_int_sz, 1, fin) < 0) {
    return -1;
  }
  k->nucsize = nucsize;
  if (read_int_array(&(k->nuc_m), k->nucsize, file_int_sz, fin) < 0 ||
      read_int_array(&(k->nuc_n), k->nucsize, file_int_sz, fin) < 0 ||
      read_ext_int_array(&(k->ell1), &(k->ell1_size), file_int_sz, fin) < 0 ||
      read_int(&sl_size, file_int_sz, 1, fin) < 0) {
    return -1;
  }
  k->sl_size = sl_size;
  if (read_int_array(&(k->sl_n), k->sl_size, file_int_sz, fin) < 0 ||
      read_int_array(&(k->sl_k), k->sl_size, file_int_sz, fin) < 0) {
    return -1;
  }
#else
  if (skip_int_array(fin, file_int_sz, dsize) < 0 ||
      read_int(&nucsize, file_int_sz, 1, fin) < 0 ||
      skip_int_array(fin, file_int_sz, nucsize) < 0 ||
      skip_int_array(fin, file_int_sz, nucsize) < 0 ||
      read_int(&ell1_size, file_int_sz, 1, fin) < 0 ||
      skip_int_array(fin, file_int_sz, ell1_size) < 0 ||
      read_int(&sl_size, file_int_sz, 1, fin) < 0 ||
      skip_int_array(fin, file_int_sz, sl_size) < 0 ||
      skip_int_array(fin, file_int_sz, sl_size) < 0) {
    return -1;
  }
#endif
  return checked_fread(&(s->time_limit_secs), sizeof(scs_float), 1, fin);
}

void SCS(write_data)(const ScsData *d, const ScsCone *k,
                     const ScsSettings *stgs) {
  FILE *fout = fopen(stgs->write_data_filename, "wb");
  scs_int status = 0;
  if (!fout) {
    scs_printf("Error: could not open %s for writing\n",
               stgs->write_data_filename);
    return;
  }
  uint32_t scs_int_sz = (uint32_t)sizeof(scs_int);
  uint32_t scs_float_sz = (uint32_t)sizeof(scs_float);
  const char *scs_version = SCS_VERSION;
  uint32_t scs_version_sz = (uint32_t)strlen(scs_version);
  if (checked_fwrite(&(scs_int_sz), sizeof(uint32_t), 1, fout) < 0 ||
      checked_fwrite(&(scs_float_sz), sizeof(uint32_t), 1, fout) < 0 ||
      checked_fwrite(&(scs_version_sz), sizeof(uint32_t), 1, fout) < 0 ||
      checked_fwrite(scs_version, 1, scs_version_sz, fout) < 0 ||
      write_scs_cone(k, fout) < 0 || write_scs_data(d, fout) < 0 ||
      write_scs_stgs(stgs, fout) < 0 ||
      write_scs_extensions(k, stgs, fout) < 0) {
    status = -1;
  }
  if (fclose(fout) != 0) {
    status = -1;
  }
  if (status < 0) {
    scs_printf("Error: failed writing SCS data to %s\n",
               stgs->write_data_filename);
  }
}

scs_int SCS(read_data)(const char *filename, ScsData **d, ScsCone **k,
                       ScsSettings **stgs) {
  uint32_t file_int_sz;
  uint32_t file_float_sz;
  uint32_t file_version_sz;
  char file_version[16];
  scs_int legacy_settings;
  errno = 0;
  FILE *fin = fopen(filename, "rb");
  *d = SCS_NULL;
  *k = SCS_NULL;
  *stgs = SCS_NULL;
  if (!fin) {
    scs_printf("Error reading file %s\n", filename);
    scs_printf("errno:%i:%s\n", errno, strerror(errno));
    return -1;
  }
  scs_printf("Reading data from %s\n", filename);
  if (checked_fread(&(file_int_sz), sizeof(uint32_t), 1, fin) < 0 ||
      checked_fread(&(file_float_sz), sizeof(uint32_t), 1, fin) < 0) {
    return read_data_cleanup(fin, d, k, stgs);
  }
  if (file_int_sz != 4 && file_int_sz != 8) {
    scs_printf("Error: unsupported file integer size %lu\n",
               (unsigned long)file_int_sz);
    return read_data_cleanup(fin, d, k, stgs);
  }
  if (file_int_sz != (uint32_t)sizeof(scs_int)) {
    scs_printf(
        "Warning, sizeof(file int) is %lu, but scs expects sizeof(int) %lu. "
        "SCS will attempt to cast the data, which may be slow. "
        "This message can be avoided by recompiling with the correct flags.\n",
        (unsigned long)file_int_sz, (unsigned long)sizeof(scs_int));
  }
  if (file_float_sz != (uint32_t)sizeof(scs_float)) {
    scs_printf(
        "Error, sizeof(file float) is %lu, but scs expects sizeof(float) %lu, "
        "scs should be recompiled with the correct flags.\n",
        (unsigned long)file_float_sz, (unsigned long)sizeof(scs_float));
    return read_data_cleanup(fin, d, k, stgs);
  }
  if (checked_fread(&(file_version_sz), sizeof(uint32_t), 1, fin) < 0) {
    return read_data_cleanup(fin, d, k, stgs);
  }
  if (file_version_sz >= sizeof(file_version)) {
    scs_printf("Error: file version string length %lu exceeds buffer size\n",
               (unsigned long)file_version_sz);
    return read_data_cleanup(fin, d, k, stgs);
  }
  if (checked_fread(file_version, 1, file_version_sz, fin) < 0) {
    return read_data_cleanup(fin, d, k, stgs);
  }
  file_version[file_version_sz] = '\0';
  legacy_settings = strcmp(file_version, SCS_VERSION) != 0;
  if (strcmp(file_version, SCS_VERSION) != 0) {
    scs_printf("************************************************************\n"
               "Warning: SCS file version %s, this is SCS version %s.\n"
               "The file reading / writing logic might have changed.\n"
               "************************************************************\n",
               file_version, SCS_VERSION);
  }
  *k = read_scs_cone(fin, file_int_sz);
  if (!*k) {
    return read_data_cleanup(fin, d, k, stgs);
  }
  *d = read_scs_data(fin, file_int_sz);
  if (!*d) {
    return read_data_cleanup(fin, d, k, stgs);
  }
  *stgs = read_scs_stgs(fin, file_int_sz, legacy_settings);
  if (!*stgs) {
    return read_data_cleanup(fin, d, k, stgs);
  }
  if (read_scs_extensions(fin, file_int_sz, *k, *stgs) < 0) {
    return read_data_cleanup(fin, d, k, stgs);
  }
  scs_printf("Finished reading data.\n");
  fclose(fin);
  return 0;
}

void SCS(log_data_to_csv)(const ScsCone *k, const ScsSettings *stgs,
                          const ScsWork *w, scs_int iter,
                          SCS(timer) * solve_timer) {
  ScsResiduals *r = w->r_orig;
  ScsResiduals *r_n = w->r_normalized;
  ScsSolution *sol = w->xys_orig;
  ScsSolution *sol_n = w->xys_normalized;
  /* if iter 0 open to write, else open to append */
  FILE *fout = fopen(stgs->log_csv_filename, iter == 0 ? "w" : "a");
  if (!fout) {
    scs_printf("Error: Could not open %s for writing\n",
               stgs->log_csv_filename);
    return;
  }
  scs_int l = w->d->m + w->d->n + 1;
  if (iter == 0) {
    /* need to end in comma so that csv parsing is correct */
    fprintf(fout, "iter,"
                  "res_pri,"
                  "res_dual,"
                  "gap,"
                  "x_nrm_inf,"
                  "y_nrm_inf,"
                  "s_nrm_inf,"
                  "x_nrm_2,"
                  "y_nrm_2,"
                  "s_nrm_2,"
                  "x_nrm_inf_normalized,"
                  "y_nrm_inf_normalized,"
                  "s_nrm_inf_normalized,"
                  "x_nrm_2_normalized,"
                  "y_nrm_2_normalized,"
                  "s_nrm_2_normalized,"
                  "ax_s_btau_nrm_inf,"
                  "px_aty_ctau_nrm_inf,"
                  "ax_s_btau_nrm_2,"
                  "px_aty_ctau_nrm_2,"
                  "res_infeas,"
                  "res_unbdd_a,"
                  "res_unbdd_p,"
                  "pobj,"
                  "dobj,"
                  "tau,"
                  "kap,"
                  "res_pri_normalized,"
                  "res_dual_normalized,"
                  "gap_normalized,"
                  "ax_s_btau_nrm_inf_normalized,"
                  "px_aty_ctau_nrm_inf_normalized,"
                  "ax_s_btau_nrm_2_normalized,"
                  "px_aty_ctau_nrm_2_normalized,"
                  "res_infeas_normalized,"
                  "res_unbdd_a_normalized,"
                  "res_unbdd_p_normalized,"
                  "pobj_normalized,"
                  "dobj_normalized,"
                  "tau_normalized,"
                  "kap_normalized,"
                  "ax_nrm_inf,"
                  "ax_s_nrm_inf,"
                  "px_nrm_inf,"
                  "aty_nrm_inf,"
                  "xt_p_x,"
                  "xt_p_x_tau,"
                  "ctx,"
                  "ctx_tau,"
                  "bty,"
                  "bty_tau,"
                  "b_nrm_inf,"
                  "c_nrm_inf,"
                  "scale,"
                  "diff_u_ut_nrm_2,"
                  "diff_v_v_prev_nrm_2,"
                  "diff_u_ut_nrm_inf,"
                  "diff_v_v_prev_nrm_inf,"
                  "aa_norm,"
                  "accepted_accel_steps,"
                  "rejected_accel_steps,"
                  "time,"
#ifdef USE_LAPACK
                  "spectral_Newton_iter,"
                  "plain_Newton_success,"
                  "res_dual_spectral,"
                  "res_pri_spectral,"
                  "comp_spectral,"
#endif
                  "\n");
  }
  fprintf(fout, "%li,", (long)iter);
  fprintf(fout, "%.16e,", r->res_pri);
  fprintf(fout, "%.16e,", r->res_dual);
  fprintf(fout, "%.16e,", r->gap);
  fprintf(fout, "%.16e,", SCS(norm_inf)(sol->x, w->d->n));
  fprintf(fout, "%.16e,", SCS(norm_inf)(sol->y, w->d->m));
  fprintf(fout, "%.16e,", SCS(norm_inf)(sol->s, w->d->m));
  fprintf(fout, "%.16e,", SCS(norm_2)(sol->x, w->d->n));
  fprintf(fout, "%.16e,", SCS(norm_2)(sol->y, w->d->m));
  fprintf(fout, "%.16e,", SCS(norm_2)(sol->s, w->d->m));
  fprintf(fout, "%.16e,", SCS(norm_inf)(sol_n->x, w->d->n));
  fprintf(fout, "%.16e,", SCS(norm_inf)(sol_n->y, w->d->m));
  fprintf(fout, "%.16e,", SCS(norm_inf)(sol_n->s, w->d->m));
  fprintf(fout, "%.16e,", SCS(norm_2)(sol_n->x, w->d->n));
  fprintf(fout, "%.16e,", SCS(norm_2)(sol_n->y, w->d->m));
  fprintf(fout, "%.16e,", SCS(norm_2)(sol_n->s, w->d->m));
  fprintf(fout, "%.16e,", SCS(norm_inf)(r->ax_s_btau, w->d->m));
  fprintf(fout, "%.16e,", SCS(norm_inf)(r->px_aty_ctau, w->d->n));
  fprintf(fout, "%.16e,", SCS(norm_2)(r->ax_s_btau, w->d->m));
  fprintf(fout, "%.16e,", SCS(norm_2)(r->px_aty_ctau, w->d->n));
  fprintf(fout, "%.16e,", r->res_infeas);
  fprintf(fout, "%.16e,", r->res_unbdd_a);
  fprintf(fout, "%.16e,", r->res_unbdd_p);
  fprintf(fout, "%.16e,", r->pobj);
  fprintf(fout, "%.16e,", r->dobj);
  fprintf(fout, "%.16e,", r->tau);
  fprintf(fout, "%.16e,", r->kap);
  fprintf(fout, "%.16e,", r_n->res_pri);
  fprintf(fout, "%.16e,", r_n->res_dual);
  fprintf(fout, "%.16e,", r_n->gap);
  fprintf(fout, "%.16e,", SCS(norm_inf)(r_n->ax_s_btau, w->d->m));
  fprintf(fout, "%.16e,", SCS(norm_inf)(r_n->px_aty_ctau, w->d->n));
  fprintf(fout, "%.16e,", SCS(norm_2)(r_n->ax_s_btau, w->d->m));
  fprintf(fout, "%.16e,", SCS(norm_2)(r_n->px_aty_ctau, w->d->n));
  fprintf(fout, "%.16e,", r_n->res_infeas);
  fprintf(fout, "%.16e,", r_n->res_unbdd_a);
  fprintf(fout, "%.16e,", r_n->res_unbdd_p);
  fprintf(fout, "%.16e,", r_n->pobj);
  fprintf(fout, "%.16e,", r_n->dobj);
  fprintf(fout, "%.16e,", r_n->tau);
  fprintf(fout, "%.16e,", r_n->kap);
  fprintf(fout, "%.16e,", SCS(norm_inf)(r->ax, w->d->m));
  fprintf(fout, "%.16e,", SCS(norm_inf)(r->ax_s, w->d->m));
  fprintf(fout, "%.16e,", SCS(norm_inf)(r->px, w->d->n));
  fprintf(fout, "%.16e,", SCS(norm_inf)(r->aty, w->d->n));
  fprintf(fout, "%.16e,", r->xt_p_x);
  fprintf(fout, "%.16e,", r->xt_p_x_tau);
  fprintf(fout, "%.16e,", r->ctx);
  fprintf(fout, "%.16e,", r->ctx_tau);
  fprintf(fout, "%.16e,", r->bty);
  fprintf(fout, "%.16e,", r->bty_tau);
  fprintf(fout, "%.16e,", SCS(norm_inf)(w->b_orig, w->d->m));
  fprintf(fout, "%.16e,", SCS(norm_inf)(w->c_orig, w->d->n));
  fprintf(fout, "%.16e,", w->stgs->scale);
  fprintf(fout, "%.16e,", SCS(norm_diff)(w->u, w->u_t, l));
  fprintf(fout, "%.16e,", SCS(norm_diff)(w->v, w->v_prev, l));
  fprintf(fout, "%.16e,", SCS(norm_inf_diff)(w->u, w->u_t, l));
  fprintf(fout, "%.16e,", SCS(norm_inf_diff)(w->v, w->v_prev, l));
  fprintf(fout, "%.16e,", w->aa_norm);
  fprintf(fout, "%li,", (long)w->accepted_accel_steps);
  fprintf(fout, "%li,", (long)w->rejected_accel_steps);
  fprintf(fout, "%.16e,", SCS(tocq)(solve_timer) / 1e3);
#ifdef USE_SPECTRAL_CONES
  fprintf(fout, "%li,", (long)w->cone_work->newton_stats.iter);
  fprintf(fout, "%li,", (long)w->cone_work->newton_stats.newton_success);
  fprintf(fout, "%.16e,", w->cone_work->newton_stats.residuals[0]);
  fprintf(fout, "%.16e,", w->cone_work->newton_stats.residuals[1]);
  fprintf(fout, "%.16e,", w->cone_work->newton_stats.residuals[2]);
#endif
  fprintf(fout, "\n");
  fclose(fout);
}

#endif
