#include "util.h"
#include <sys/time.h>

static struct timeval tic_timestart;

void tic(void) {
	gettimeofday(&tic_timestart, NULL);
}

pfloat tocq(void) {
	struct timeval tic_timestop;
	gettimeofday(&tic_timestop, NULL);
	/* scs_printf("time: %8.4f seconds.\n", (float)(tic_timestop - tic_timestart)); */
	return tic_timestop.tv_sec * 1e3 + tic_timestop.tv_usec / 1e3 - tic_timestart.tv_sec * 1e3
			- tic_timestart.tv_usec / 1e3;
}

pfloat toc(void) {
	pfloat time = tocq();
	scs_printf("time: %8.4f milli-seconds.\n", time);
	return time;
}

void printConeData(Cone * k) {
	idxint i;
	scs_printf("num zeros = %i\n", (int) k->f);
	scs_printf("num LP = %i\n", (int) k->l);
	scs_printf("num SOCs = %i\n", (int) k->qsize);
	scs_printf("soc array:\n");
	for (i = 0; i < k->qsize; i++) {
		scs_printf("%i\n", (int) k->q[i]);
	}
	scs_printf("num SDCs = %i\n", (int) k->ssize);
	scs_printf("sdc array:\n");
	for (i = 0; i < k->ssize; i++) {
		scs_printf("%i\n", (int) k->s[i]);
	}
	scs_printf("num ep = %i\n", (int) k->ep);
	scs_printf("num ed = %i\n", (int) k->ed);
}

void printData(Data * d) {
	scs_printf("d->n is %i\n", (int) d->n);
	scs_printf("d->m is %i\n", (int) d->m);
	scs_printf("d->b[0] is %4f\n", d->b[0]);
	scs_printf("d->c[0] is %4f\n", d->c[0]);
	scs_printf("d->Ax[0] is %4f\n", d->Ax[0]);
	scs_printf("Anz is %li\n", (long) d->Ap[d->n]);
	scs_printf("d->MAX_ITERS is %i\n", (int) d->MAX_ITERS);
	scs_printf("d->VERBOSE is %i\n", (int) d->VERBOSE);
	scs_printf("d->NORMALIZE is %i\n", (int) d->VERBOSE);
	scs_printf("d->ALPHA is %6f\n", d->ALPHA);
	scs_printf("d->EPS is %6f\n", d->EPS);
	scs_printf("d->EPS is %6f\n", d->EPS);
}

void printWork(Data * d, Work * w) {
	idxint i, l = d->n + d->m;
	scs_printf("\n u_t is \n");
	for (i = 0; i < l; i++) {
		scs_printf("%f\n", w->u_t[i]);
	}
	scs_printf("\n u is \n");
	for (i = 0; i < l; i++) {
		scs_printf("%f\n", w->u[i]);
	}
	scs_printf("\n v is \n");
	for (i = 0; i < l; i++) {
		scs_printf("%f\n", w->v[i]);
	}
}

