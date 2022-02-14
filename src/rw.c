#include "rw.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "linalg.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

/* writes/reads problem data to/from filename */
/* This is a VERY naive implementation, doesn't care about portability etc */

static void write_scs_cone(const ScsCone *k, FILE *fout) {
  fwrite(&(k->z), sizeof(scs_int), 1, fout);
  fwrite(&(k->l), sizeof(scs_int), 1, fout);
  fwrite(&(k->bsize), sizeof(scs_int), 1, fout);
  fwrite(k->bl, sizeof(scs_float), MAX(k->bsize - 1, 0), fout);
  fwrite(k->bu, sizeof(scs_float), MAX(k->bsize - 1, 0), fout);
  fwrite(&(k->qsize), sizeof(scs_int), 1, fout);
  fwrite(k->q, sizeof(scs_int), k->qsize, fout);
  fwrite(&(k->ssize), sizeof(scs_int), 1, fout);
  fwrite(k->s, sizeof(scs_int), k->ssize, fout);
  fwrite(&(k->ep), sizeof(scs_int), 1, fout);
  fwrite(&(k->ed), sizeof(scs_int), 1, fout);
  fwrite(&(k->psize), sizeof(scs_int), 1, fout);
  fwrite(k->p, sizeof(scs_float), k->psize, fout);
}

static ScsCone *read_scs_cone(FILE *fin) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  fread(&(k->z), sizeof(scs_int), 1, fin);
  fread(&(k->l), sizeof(scs_int), 1, fin);
  fread(&(k->bsize), sizeof(scs_int), 1, fin);
  if (k->bsize > 1) {
    k->bl = (scs_float *)scs_calloc(MAX(k->bsize - 1, 0), sizeof(scs_float));
    k->bu = (scs_float *)scs_calloc(MAX(k->bsize - 1, 0), sizeof(scs_float));
    fread(k->bl, sizeof(scs_float), MAX(k->bsize - 1, 0), fin);
    fread(k->bu, sizeof(scs_float), MAX(k->bsize - 1, 0), fin);
  }
  fread(&(k->qsize), sizeof(scs_int), 1, fin);
  if (k->qsize) {
    k->q = (scs_int *)scs_calloc(k->qsize, sizeof(scs_int));
    fread(k->q, sizeof(scs_int), k->qsize, fin);
  }
  fread(&(k->ssize), sizeof(scs_int), 1, fin);
  if (k->ssize) {
    k->s = (scs_int *)scs_calloc(k->ssize, sizeof(scs_int));
    fread(k->s, sizeof(scs_int), k->ssize, fin);
  }
  fread(&(k->ep), sizeof(scs_int), 1, fin);
  fread(&(k->ed), sizeof(scs_int), 1, fin);
  fread(&(k->psize), sizeof(scs_int), 1, fin);
  if (k->psize) {
    k->p = (scs_float *)scs_calloc(k->psize, sizeof(scs_float));
    fread(k->p, sizeof(scs_float), k->psize, fin);
  }
  return k;
}

static void write_scs_stgs(const ScsSettings *s, FILE *fout) {
  /* Warm start to false for now */
  scs_int warm_start = 0;
  fwrite(&(s->normalize), sizeof(scs_int), 1, fout);
  fwrite(&(s->scale), sizeof(scs_float), 1, fout);
  fwrite(&(s->rho_x), sizeof(scs_float), 1, fout);
  fwrite(&(s->max_iters), sizeof(scs_int), 1, fout);
  fwrite(&(s->eps_abs), sizeof(scs_float), 1, fout);
  fwrite(&(s->eps_rel), sizeof(scs_float), 1, fout);
  fwrite(&(s->eps_infeas), sizeof(scs_float), 1, fout);
  fwrite(&(s->alpha), sizeof(scs_float), 1, fout);
  fwrite(&(s->verbose), sizeof(scs_int), 1, fout);
  fwrite(&warm_start, sizeof(scs_int), 1, fout);
  fwrite(&(s->acceleration_lookback), sizeof(scs_int), 1, fout);
  fwrite(&(s->acceleration_interval), sizeof(scs_int), 1, fout);
  fwrite(&(s->adaptive_scale), sizeof(scs_int), 1, fout);
  /* Do not write the write_data_filename */
  /* Do not write the log_csv_filename */
}

static ScsSettings *read_scs_stgs(FILE *fin) {
  ScsSettings *s = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  fread(&(s->normalize), sizeof(scs_int), 1, fin);
  fread(&(s->scale), sizeof(scs_float), 1, fin);
  fread(&(s->rho_x), sizeof(scs_float), 1, fin);
  fread(&(s->max_iters), sizeof(scs_int), 1, fin);
  fread(&(s->eps_abs), sizeof(scs_float), 1, fin);
  fread(&(s->eps_rel), sizeof(scs_float), 1, fin);
  fread(&(s->eps_infeas), sizeof(scs_float), 1, fin);
  fread(&(s->alpha), sizeof(scs_float), 1, fin);
  fread(&(s->verbose), sizeof(scs_int), 1, fin);
  fread(&(s->warm_start), sizeof(scs_int), 1, fin);
  fread(&(s->acceleration_lookback), sizeof(scs_int), 1, fin);
  fread(&(s->acceleration_interval), sizeof(scs_int), 1, fin);
  fread(&(s->adaptive_scale), sizeof(scs_int), 1, fin);
  return s;
}

static void write_amatrix(const ScsMatrix *A, FILE *fout) {
  scs_int Anz = A->p[A->n];
  fwrite(&(A->m), sizeof(scs_int), 1, fout);
  fwrite(&(A->n), sizeof(scs_int), 1, fout);
  fwrite(A->p, sizeof(scs_int), A->n + 1, fout);
  fwrite(A->x, sizeof(scs_float), Anz, fout);
  fwrite(A->i, sizeof(scs_int), Anz, fout);
}

static ScsMatrix *read_amatrix(FILE *fin) {
  scs_int Anz;
  ScsMatrix *A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  fread(&(A->m), sizeof(scs_int), 1, fin);
  fread(&(A->n), sizeof(scs_int), 1, fin);
  A->p = (scs_int *)scs_calloc(A->n + 1, sizeof(scs_int));
  fread(A->p, sizeof(scs_int), A->n + 1, fin);
  Anz = A->p[A->n];
  A->x = (scs_float *)scs_calloc(Anz, sizeof(scs_float));
  A->i = (scs_int *)scs_calloc(Anz, sizeof(scs_int));
  fread(A->x, sizeof(scs_float), Anz, fin);
  fread(A->i, sizeof(scs_int), Anz, fin);
  return A;
}

static void write_scs_data(const ScsData *d, FILE *fout) {
  scs_int has_p = d->P ? 1 : 0;
  fwrite(&(d->m), sizeof(scs_int), 1, fout);
  fwrite(&(d->n), sizeof(scs_int), 1, fout);
  fwrite(d->b, sizeof(scs_float), d->m, fout);
  fwrite(d->c, sizeof(scs_float), d->n, fout);
  write_amatrix(d->A, fout);
  /* write has P bit */
  fwrite(&has_p, sizeof(scs_int), 1, fout);
  if (d->P) {
    write_amatrix(d->P, fout);
  }
}

static ScsData *read_scs_data(FILE *fin) {
  scs_int has_p = 0;
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  fread(&(d->m), sizeof(scs_int), 1, fin);
  fread(&(d->n), sizeof(scs_int), 1, fin);
  d->b = (scs_float *)scs_calloc(d->m, sizeof(scs_float));
  d->c = (scs_float *)scs_calloc(d->n, sizeof(scs_float));
  fread(d->b, sizeof(scs_float), d->m, fin);
  fread(d->c, sizeof(scs_float), d->n, fin);
  d->A = read_amatrix(fin);
  /* If has_p bit is not set or this hits end of file then has_p = 0 */
  has_p &= fread(&has_p, sizeof(scs_int), 1, fin);
  d->P = has_p ? read_amatrix(fin) : SCS_NULL;
  return d;
}

void SCS(write_data)(const ScsData *d, const ScsCone *k,
                     const ScsSettings *stgs) {
  FILE *fout = fopen(stgs->write_data_filename, "wb");
  uint32_t scs_int_sz = (uint32_t)sizeof(scs_int);
  uint32_t scs_float_sz = (uint32_t)sizeof(scs_float);
  const char *scs_version = SCS_VERSION;
  uint32_t scs_version_sz = (uint32_t)strlen(scs_version);
  scs_printf("writing data to %s\n", stgs->write_data_filename);
  fwrite(&(scs_int_sz), sizeof(uint32_t), 1, fout);
  fwrite(&(scs_float_sz), sizeof(uint32_t), 1, fout);
  fwrite(&(scs_version_sz), sizeof(uint32_t), 1, fout);
  fwrite(scs_version, 1, scs_version_sz, fout);
  write_scs_cone(k, fout);
  write_scs_data(d, fout);
  write_scs_stgs(stgs, fout);
  fclose(fout);
}

scs_int SCS(read_data)(const char *filename, ScsData **d, ScsCone **k,
                       ScsSettings **stgs) {
  uint32_t file_int_sz;
  uint32_t file_float_sz;
  uint32_t file_version_sz;
  char file_version[16];
  FILE *fin = fopen(filename, "rb");
  if (!fin) {
    scs_printf("Error reading file %s\n", filename);
    return -1;
  }
  scs_printf("Reading data from %s\n", filename);
  fread(&(file_int_sz), sizeof(uint32_t), 1, fin);
  fread(&(file_float_sz), sizeof(uint32_t), 1, fin);
  if (file_int_sz != (uint32_t)sizeof(scs_int)) {
    scs_printf(
        "Error, sizeof(file int) is %lu, but scs expects sizeof(int) %lu, "
        "scs should be recompiled with correct flags.\n",
        (unsigned long)file_int_sz, (unsigned long)sizeof(scs_int));
    fclose(fin);
    return -1;
  }
  if (file_float_sz != (uint32_t)sizeof(scs_float)) {
    scs_printf(
        "Error, sizeof(file float) is %lu, but scs expects sizeof(float) %lu, "
        "scs should be recompiled with the correct flags.\n",
        (unsigned long)file_float_sz, (unsigned long)sizeof(scs_float));
    fclose(fin);
    return -1;
  }
  fread(&(file_version_sz), sizeof(uint32_t), 1, fin);
  fread(file_version, 1, file_version_sz, fin);
  file_version[file_version_sz] = '\0';
  if (strcmp(file_version, SCS_VERSION) != 0) {
    scs_printf("************************************************************\n"
               "Warning: SCS file version %s, this is SCS version %s.\n"
               "The file reading / writing logic might have changed.\n"
               "************************************************************\n",
               file_version, SCS_VERSION);
  }
  *k = read_scs_cone(fin);
  *d = read_scs_data(fin);
  *stgs = read_scs_stgs(fin);
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
                  "px_nrm_inf,"
                  "aty_nrm_inf,"
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
  fprintf(fout, "%.16e,", SCS(norm_inf)(r->px, w->d->n));
  fprintf(fout, "%.16e,", SCS(norm_inf)(r->aty, w->d->n));
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
  fprintf(fout, "\n");
  fclose(fout);
}
