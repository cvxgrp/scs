#ifndef NORMALIZE_H_GUARD
#define NORMALIZE_H_GUARD

#define MIN_SCALE (1e-3)
#define MAX_SCALE (1e3)

void normalizeBC(Work *w) {
    scs_int i;
    scs_float nm, *D = w->scal->D, *E = w->scal->E, *b = w->b, *c = w->c;
    /* scale b */
    for (i = 0; i < w->m; ++i) {
        b[i] /= D[i];
    }
    nm = calcNorm(b, w->m);
    w->sc_b = w->scal->meanNormColA / MAX(nm, MIN_SCALE);
    /* scale c */
    for (i = 0; i < w->n; ++i) {
        c[i] /= E[i];
    }
    nm = calcNorm(c, w->n);
    w->sc_c = w->scal->meanNormRowA / MAX(nm, MIN_SCALE);
    scaleArray(b, w->sc_b * w->stgs->scale, w->m);
    scaleArray(c, w->sc_c * w->stgs->scale, w->n);
}

void calcScaledResids(Work *w, struct residuals *r) {
    scs_float *D = w->scal->D;
    scs_float *E = w->scal->E;
    scs_float *u = w->u;
    scs_float *u_t = w->u_t;
    scs_float *u_prev = w->u_prev;
    scs_float tmp;
    scs_int i, n = w->n, m = w->m;

    r->resPri = 0;
    for (i = 0; i < n; ++i) {
        tmp = (u[i] - u_t[i]) / (E[i] * w->sc_b);
        r->resPri += tmp * tmp;
    }
    for (i = 0; i < m; ++i) {
        tmp = (u[i + n] - u_t[i + n]) / (D[i] * w->sc_c);
        r->resPri += tmp * tmp;
    }
    tmp = u[n + m] - u_t[n + m];
    r->resPri += tmp * tmp;
    r->resPri = sqrt(r->resPri);

    r->resDual = 0;
    for (i = 0; i < n; ++i) {
        tmp = (u[i] - u_prev[i]) * E[i] / w->sc_b;
        r->resDual += tmp * tmp;
    }
    for (i = 0; i < m; ++i) {
        tmp = (u[i + n] - u_prev[i + n]) * D[i] / w->sc_c;
        r->resDual += tmp * tmp;
    }
    tmp = u[n + m] - u_t[n + m];
    r->resDual += tmp * tmp;
    r->resDual = sqrt(r->resDual);
}

void normalizeWarmStart(Work *w) {
    scs_int i;
    scs_float *D = w->scal->D;
    scs_float *E = w->scal->E;
    scs_float *x = w->u;
    scs_float *y = &(w->u[w->n]);
    scs_float *s = &(w->v[w->n]);
    for (i = 0; i < w->n; ++i) {
        x[i] *= (E[i] * w->sc_b);
    }
    for (i = 0; i < w->m; ++i) {
        y[i] *= (D[i] * w->sc_c);
    }
    for (i = 0; i < w->m; ++i) {
        s[i] /= (D[i] / (w->sc_b * w->stgs->scale));
    }
}

void unNormalizeSol(Work *w, Sol *sol) {
    scs_int i;
    scs_float *D = w->scal->D;
    scs_float *E = w->scal->E;
    for (i = 0; i < w->n; ++i) {
        sol->x[i] /= (E[i] * w->sc_b);
    }
    for (i = 0; i < w->m; ++i) {
        sol->y[i] /= (D[i] * w->sc_c);
    }
    for (i = 0; i < w->m; ++i) {
        sol->s[i] *= D[i] / (w->sc_b * w->stgs->scale);
    }
}

/* unused:
void unNormalizeBC(Work * w, Sol * sol) {
        scs_int i;
        scs_float * D = w->D;
        scs_float * E = w->E;
        for (i = 0; i < w->n; ++i) {
                w->c[i] *= E[i] / (w->sc_c * w->scale);
        }
        for (i = 0; i < w->m; ++i) {
                w->b[i] *= D[i] / (w->sc_b * w->scale);
        }
}
*/
#endif
