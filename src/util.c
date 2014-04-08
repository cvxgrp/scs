#include "util.h"
/* return milli-seconds */
#if (defined WIN32 || _WIN64 || defined _WINDLL)
void tic(timer* t)
{
	QueryPerformanceFrequency(&t->freq);
	QueryPerformanceCounter(&t->tic);
}

pfloat tocq(timer* t)
{
	QueryPerformanceCounter(&t->toc);
	return (1e3 * (t->toc.QuadPart - t->tic.QuadPart) / (pfloat) t->freq.QuadPart);
}
#elif (defined __APPLE__)
void tic(timer* t) {
	/* read current clock cycles */
	t->tic = mach_absolute_time();
}

pfloat tocq(timer* t) {

	uint64_t duration; /* elapsed time in clock cycles*/

	t->toc = mach_absolute_time();
	duration = t->toc - t->tic;

	/*conversion from clock cycles to nanoseconds*/
	mach_timebase_info(&(t->tinfo));
	duration *= t->tinfo.numer;
	duration /= t->tinfo.denom;

	return (pfloat) duration / 1e6;
}
#else
void tic(timer* t)
{
	clock_gettime(CLOCK_MONOTONIC, &t->tic);
}

pfloat tocq(timer* t)
{
	struct timespec temp;

	clock_gettime(CLOCK_MONOTONIC, &t->toc);

	if ((t->toc.tv_nsec - t->tic.tv_nsec)<0) {
		temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec-1;
		temp.tv_nsec = 1e9+t->toc.tv_nsec - t->tic.tv_nsec;
	} else {
		temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec;
		temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
	}
	return (pfloat)temp.tv_sec * 1e3 + (pfloat)temp.tv_nsec / 1e6;
}
#endif

pfloat toc(timer * t) {
	pfloat time = tocq(t);
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

void printData(Data * d) {
	scs_printf("m = %i\n", (int) d->m);
	scs_printf("n = %i\n", (int) d->n);

	scs_printf("b[0] = %4f\n", d->b[0]);
	scs_printf("c[0] = %4f\n", d->c[0]);

	scs_printf("MAX_ITERS = %i\n", (int) d->MAX_ITERS);
	scs_printf("VERBOSE = %i\n", (int) d->VERBOSE);
	scs_printf("NORMALIZE = %i\n", (int) d->NORMALIZE);
	scs_printf("WARM_START = %i\n", (int) d->WARM_START);
	scs_printf("EPS = %4f\n", d->EPS);
	scs_printf("ALPHA = %4f\n",  d->ALPHA);
	scs_printf("RHO_X = %4f\n", d->RHO_X);
	scs_printf("CG_RATE = %4f\n", d->CG_RATE);
	scs_printf("SCALE = %4f\n", d->SCALE);
}

