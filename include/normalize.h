#ifndef NORMALIZE_H_GUARD
#define NORMALIZE_H_GUARD

#define MIN_SCALE (1e-3)
#define MAX_SCALE (1e3)

void normalizeBC(Data * d, Work * w) {
	scs_int i;
	scs_float *D = w->D, *E = w->E;
	/*
    scs_printf("norm b = %4f\n", calcNorm(d->b, d->m));
    scs_printf("norm c = %4f\n", calcNorm(d->b, d->n));
    */
    /* scale b */
	for (i = 0; i < d->m; ++i) {
		d->b[i] /= D[i];
	}
	w->sc_b = w->meanNormColA / MAX(calcNorm(d->b, d->m), MIN_SCALE);
	/* scale c */
	for (i = 0; i < d->n; ++i) {
		d->c[i] /= E[i];
	}
	w->sc_c = w->meanNormRowA / MAX(calcNorm(d->c, d->n), MIN_SCALE);
	scaleArray(d->b, w->sc_b * d->scale, d->m);
	scaleArray(d->c, w->sc_c * d->scale, d->n);
}

void calcScaledResids(Data * d, Work * w, struct residuals * r) {
	scs_float * D = w->D;
	scs_float * E = w->E;
	scs_float * u = w->u;
	scs_float * u_t = w->u_t;
	scs_float * u_prev = w->u_prev;
	scs_float tmp;
	scs_int i, n = d->n, m = d->m;

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

void normalizeWarmStart(Data *d, Work * w) {
	scs_int i;
	scs_float * D = w->D;
	scs_float * E = w->E;
	scs_float * x = w->u;
	scs_float * y = &(w->u[d->n]);
	scs_float * s = &(w->v[d->n]);
	for (i = 0; i < d->n; ++i) {
		x[i] *= (E[i] * w->sc_b);
	}
	for (i = 0; i < d->m; ++i) {
		y[i] *= (D[i] * w->sc_c);
	}
	for (i = 0; i < d->m; ++i) {
		s[i] /= D[i] / (w->sc_b * d->scale);
	}
}

void unNormalizeSolBC(Data *d, Work * w, Sol * sol) {
	scs_int i;
	scs_float * D = w->D;
	scs_float * E = w->E;
	for (i = 0; i < d->n; ++i) {
		sol->x[i] /= (E[i] * w->sc_b);
	}
	for (i = 0; i < d->m; ++i) {
		sol->y[i] /= (D[i] * w->sc_c);
	}
	for (i = 0; i < d->m; ++i) {
		sol->s[i] *= D[i] / (w->sc_b * d->scale);
	}
	for (i = 0; i < d->n; ++i) {
		d->c[i] *= E[i] / (w->sc_c * d->scale);
	}
	for (i = 0; i < d->m; ++i) {
		d->b[i] *= D[i] / (w->sc_b * d->scale);
	}
}

#endif
