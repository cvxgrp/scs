/*
 * Interface for SCS signal handling.
 */

#ifndef CTRLC_H_GUARD
#define CTRLC_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#if CTRLC > 0

void scs_start_interrupt_listener(void);
void scs_end_interrupt_listener(void);
int scs_is_interrupted(void);

#else /* CTRLC = 0 */

/* Simply to suppress empty translation unit warnings. */
typedef int scs_make_iso_compilers_happy;

/* No signal handling. */
#define scs_start_interrupt_listener()
#define scs_end_interrupt_listener()
#define scs_is_interrupted() 0

#endif /* END IF CTRLC > 0 */

#ifdef __cplusplus
}
#endif
#endif
