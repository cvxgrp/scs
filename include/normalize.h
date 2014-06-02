#ifndef NORMALIZE_H_GUARD
#define NORMALIZE_H_GUARD

#define MIN_SCALE 1e-3
#define MAX_SCALE 1e3

void normalizeBC(Data * d, Work * w) {
	idxint i;
	pfloat *D = w->D, *E = w->E;
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
	scaleArray(d->b, w->sc_b * d->SCALE, d->m);
	scaleArray(d->c, w->sc_c * d->SCALE, d->n);
}

void calcScaledResids(Data * d, Work * w, struct residuals * r) {
	pfloat * D = w->D;
	pfloat * E = w->E;
	pfloat * u = w->u;
	pfloat * u_t = w->u_t;
	pfloat * u_prev = w->u_prev;
	pfloat tmp;
	idxint i, n = d->n, m = d->m;

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
	idxint i;
	pfloat * D = w->D;
	pfloat * E = w->E;
	pfloat * x = w->u;
	pfloat * y = &(w->u[d->n]);
	pfloat * s = &(w->v[d->n]);
	for (i = 0; i < d->n; ++i) {
		x[i] *= (E[i] * w->sc_b);
	}
	for (i = 0; i < d->m; ++i) {
		y[i] *= (D[i] * w->sc_c);
	}
	for (i = 0; i < d->m; ++i) {
		s[i] /= D[i] / (w->sc_b * d->SCALE);
	}
}

void unNormalizeSolBC(Data *d, Work * w, Sol * sol) {
	idxint i;
	pfloat * D = w->D;
	pfloat * E = w->E;
	for (i = 0; i < d->n; ++i) {
		sol->x[i] /= (E[i] * w->sc_b);
	}
	for (i = 0; i < d->m; ++i) {
		sol->y[i] /= (D[i] * w->sc_c);
	}
	for (i = 0; i < d->m; ++i) {
		sol->s[i] *= D[i] / (w->sc_b * d->SCALE);
	}
	for (i = 0; i < d->n; ++i) {
		d->c[i] *= E[i] / (w->sc_c * d->SCALE);
	}
	for (i = 0; i < d->m; ++i) {
		d->b[i] *= D[i] / (w->sc_b * d->SCALE);
	}
}

#endif
