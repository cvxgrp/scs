#ifndef NORMAL_H_GUARD
#define NORMAL_H_GUARD

#define MIN_SCALE 1e-2
#define MAX_SCALE 1e3

void normalizeA(Data * d, Work * w, Cone * k){
    
	pfloat * D = scs_calloc(d->m, sizeof(pfloat));
	pfloat * E = scs_calloc(d->n, sizeof(pfloat));
	
	idxint i, j, count;
	pfloat wrk, *nms;

	/* heuristic rescaling, seems to do well with a scaling of about 1 */
	w->scale = 1; /*MAX( MIN( sqrt( d->n * ((pfloat) d->m / d->Ap[d->n] ) , MAX_SCALE) , 1); */

	/* calculate row norms */
	for(i = 0; i < d->n; ++i){
		for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
			wrk = d->Ax[j];
			D[d->Ai[j]] += wrk*wrk;
		}
	}
	for (i=0; i < d->m; ++i){
		D[i] = sqrt(D[i]); /* just the norms */
	}
    /* mean of norms of rows across each cone  */
    count = k->l+k->f;
	for(i = 0; i < k->qsize; ++i)
    {
		wrk = 0;
	    /*
        for (j = count; j < count + k->q[i]; ++j){
        	wrk = MAX(wrk,D[j]);
		}
        */
        for (j = count; j < count + k->q[i]; ++j){
        	wrk += D[j];
		}
        wrk /= k->q[i];
        for (j = count; j < count + k->q[i]; ++j){
        	D[j] = wrk;
		}
		count += k->q[i];
    }
    for (i=0; i < k->ssize; ++i)
	{
 		wrk = 0;
        /*
        for (j = count; j < count + (k->s[i])*(k->s[i]); ++j){
        	wrk = MAX(wrk,D[j]);
		}
        */
		for (j = count; j < count + (k->s[i])*(k->s[i]); ++j){
        	wrk += D[j]*D[j];
		}
        wrk = sqrt(wrk);
        wrk /= k->s[i];
        for (j = count; j < count + (k->s[i])*(k->s[i]); ++j){
        	D[j] = wrk;
		}
		count += (k->s[i])*(k->s[i]);
    }
    
    for(i = 0; i < k->ep + k->ed; ++i)
    {
        wrk = D[count]/3 + D[count + 1]/3 + D[count + 2]/3;
        D[count] = wrk;
        D[count + 1] = wrk;
        D[count + 2] = wrk;
        count += 3;
    }

    for (i=0; i<d->m; ++i){
        if (D[i] < MIN_SCALE) D[i] = MIN_SCALE;
        else if (D[i] > MAX_SCALE) D[i] = MAX_SCALE;

    }
	/* scale the rows with D */
	for(i = 0; i < d->n; ++i){
		for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
			d->Ax[j] /= D[d->Ai[j]];
		}
	}
	/* calculate and scale by col norms, E */
	for (i = 0; i < d->n; ++i){
		E[i] = calcNorm(&(d->Ax[d->Ap[i]]),d->Ap[i+1] - d->Ap[i]);
        if (E[i] < MIN_SCALE) E[i] = MIN_SCALE;
		else if (E[i] > MAX_SCALE) E[i] = MAX_SCALE;
        scaleArray(&(d->Ax[d->Ap[i]]), 1.0/E[i], d->Ap[i+1] - d->Ap[i]);
	}

    nms = scs_calloc(d->m,sizeof(pfloat));
    for(i = 0; i < d->n; ++i){
        for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
            wrk = d->Ax[j];
            nms[d->Ai[j]] += wrk*wrk;
        }
    }
    w->meanNormRowA = 0.0;
    for (i=0; i < d->m; ++i){
        w->meanNormRowA += sqrt(nms[i])/d->m;
    }
    scs_free(nms);
    
    w->D = D;
    w->E = E;
    scaleArray(d->Ax, w->scale, d->Ap[d->n]);

    /*
       scs_printf("norm D is %4f\n", calcNorm(D,d->m));
       scs_printf("norm E is %4f\n", calcNorm(E,d->n));
       scs_printf("norm A is %4f\n", calcNorm(d->Ax, d->Ap[d->n]));
       scs_printf("norm b is %4f\n", calcNorm(d->b,d->m));
       scs_printf("norm c is %4f\n", calcNorm(d->c,d->n));
     */
}

void normalizeBC(Data * d, Work * w) {
	idxint i;
	pfloat *D = w->D, *E = w->E;
    /* scale b */
    for (i = 0; i < d->m; ++i){
        d->b[i] /= D[i];
    }
    w->sc_b = 1/MAX(calcNorm(d->b,d->m),MIN_SCALE);
    /* scale c */
    for (i = 0; i < d->n; ++i){
        d->c[i] /= E[i];
    }
    w->sc_c = w->meanNormRowA/MAX(calcNorm(d->c,d->n),MIN_SCALE);
    
    scaleArray(d->b, w->sc_b * w->scale, d->m);
    scaleArray(d->c, w->sc_c * w->scale, d->n);
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
    for (i = 0; i < n; ++i){
        tmp = (u[i] - u_t[i])/(E[i] * w->sc_b);
        r->resPri += tmp * tmp;        
    }
    for (i = 0; i < m; ++i){
        tmp = (u[i + n] - u_t[i + n])/(D[i] * w->sc_c);
        r->resPri += tmp * tmp;
    }
    tmp = u[n + m] - u_t[n + m];
    r->resPri += tmp * tmp;
    r->resPri = sqrt(r->resPri);

    r->resDual = 0;
    for (i = 0; i < n; ++i){
        tmp = (u[i] - u_prev[i]) * E[i] / w->sc_b;
        r->resDual += tmp * tmp;        
    }
    for (i = 0; i < m; ++i){
        tmp = (u[i + n] - u_prev[i + n]) * D[i] / w->sc_c;
        r->resDual += tmp * tmp;
    }
    tmp = u[n + m] - u_t[n + m];
    r->resDual += tmp * tmp;
    r->resDual = sqrt(r->resDual);
}

pfloat calcScaledNormInf(const pfloat *a, const pfloat * s, idxint l){ 
    pfloat tmp, max = 0.0;
    idxint i;
    for ( i=0; i<l; ++i){
        tmp = ABS(s[i] * a[i]);
        if(tmp > max) max = tmp;
    }   
    return max;
}

void normalizeWarmStart(Data *d, Work * w){
	idxint i;
	pfloat * D = w->D;
	pfloat * E = w->E;
    pfloat * x = w->u;
    pfloat * y = &(w->u[d->n]);
    pfloat * s = &(w->v[d->n]);
	for (i = 0; i < d->n; ++i){
		x[i] *= (E[i] * w->sc_b);
	}
	for (i = 0; i < d->m; ++i){
		y[i] *= (D[i] * w->sc_c);
	}
	for (i = 0; i < d->m; ++i){
		s[i] /= D[i]/(w->sc_b * w->scale);
	}
}

void unNormalizeA(Data *d, Work * w) {
    idxint i, j;
    pfloat * D = w->D;
    pfloat * E = w->E;
    scaleArray(d->Ax, 1.0/w->scale, d->Ap[d->n]);
    for (i = 0; i < d->n; ++i){
        scaleArray(&(d->Ax[d->Ap[i]]), E[i], d->Ap[i+1] - d->Ap[i]);
    }   
    for(i = 0; i < d->n; ++i){
        for(j = d->Ap[i]; j < d->Ap[i+1]; ++j) {
            d->Ax[j] *= D[d->Ai[j]];
        } 
    }
}


void unNormalizeSolBC(Data *d, Work * w, Sol * sol){
    idxint i;
    pfloat * D = w->D;
    pfloat * E = w->E;
    for (i = 0; i < d->n; ++i){
        sol->x[i] /= (E[i] * w->sc_b);
    }
    for (i = 0; i < d->m; ++i){
        sol->y[i] /= (D[i] * w->sc_c);
    }
    for (i = 0; i < d->m; ++i){
        sol->s[i] *= D[i]/(w->sc_b * w->scale);
    }
    for (i = 0; i < d->n; ++i){
        d->c[i] *= E[i]/(w->sc_c * w->scale);
    }
    for (i = 0; i < d->m; ++i){
        d->b[i] *= D[i]/(w->sc_b * w->scale);
    }
}

#endif
