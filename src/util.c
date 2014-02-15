#include "util.h"
/* return milli-seconds */
#if (defined WIN32 || _WIN64 || defined _WINDLL)
void tic(timer* t)
{
	QueryPerformanceFrequency(&t->freq);
	QueryPerformanceCounter(&t->tic);
}

double tocq(timer* t)
{
	QueryPerformanceCounter(&t->toc);
	return (1e3 * (t->toc.QuadPart - t->tic.QuadPart) / (double) t->freq.QuadPart);
}
#elif (defined __APPLE__)
void tic(timer* t)
{
    /* read current clock cycles */
    t->tic = mach_absolute_time();
}

double tocq(timer* t)
{

    uint64_t duration; /* elapsed time in clock cycles*/

    t->toc = mach_absolute_time();
    duration = t->toc - t->tic;

    /*conversion from clock cycles to nanoseconds*/
    mach_timebase_info(&(t->tinfo));
    duration *= t->tinfo.numer;
    duration /= t->tinfo.denom;

    return (double)duration / 1e6;
}
#else
void tic(timer* t)
{
	clock_gettime(CLOCK_MONOTONIC, &t->tic);
}

double tocq(timer* t)
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
	return (double)temp.tv_sec * 1e3 + (double)temp.tv_nsec / 1e6;
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

