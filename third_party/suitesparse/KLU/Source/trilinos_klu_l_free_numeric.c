/* ========================================================================== */
/* === TRILINOS_KLU_free_numeric ===================================================== */
/* ========================================================================== */

/* Free the KLU Numeric object. */

/* This file should make the long int version of KLU */
#ifndef KLS_KLU_INT_VARIANT
#define DLONG 1
#endif

#include "trilinos_klu_internal.h"
#include <stdlib.h>

Int TRILINOS_KLU_free_numeric
(
    TRILINOS_KLU_numeric **NumericHandle,
    TRILINOS_KLU_common	*Common
)
{
    TRILINOS_KLU_numeric *Numeric ;
    Unit **LUbx ;
    size_t *LUsize ;
    Int block, n, nzoff, nblocks ;

    if (Common == NULL)
    {
	return (FALSE) ;
    }
    if (NumericHandle == NULL || *NumericHandle == NULL)
    {
	return (TRUE) ;
    }

    Numeric = *NumericHandle ;

    /* task #21 free-tracking shim: record the freed numeric's array
       pointers with a backtrace so a later solve on a stale numeric
       can name the premature-free site (KLS_TRACK_NUMERIC_FREE) */
    {
	static int kls_track = -1 ;
	if (kls_track == -1)
	{
	    kls_track = getenv ("KLS_TRACK_NUMERIC_FREE") != NULL ;
	}
	if (kls_track)
	{
	    extern void kls_numeric_free_log (const void *numeric,
					      const void *lip,
					      const void *llen,
					      const void *uip,
					      const void *ulen) ;
	    kls_numeric_free_log (Numeric, Numeric->Lip, Numeric->Llen,
				  Numeric->Uip, Numeric->Ulen) ;
	}
    }

    n = Numeric->n ;
    nzoff = Numeric->nzoff ;
    nblocks = Numeric->nblocks ;
    LUsize = Numeric->LUsize ;

    LUbx = (Unit **) Numeric->LUbx ;
    if (LUbx != NULL)
    {
	for (block = 0 ; block < nblocks ; block++)
	{
	    TRILINOS_KLU_free (LUbx [block], LUsize ? LUsize [block] : 0,
		sizeof (Unit), Common) ;
	}
    }

    TRILINOS_KLU_free (Numeric->Pnum, n, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Numeric->Offp, n+1, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Numeric->Offi, nzoff+1, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Numeric->Offx, nzoff+1, sizeof (Entry), Common) ;

    TRILINOS_KLU_free (Numeric->Lip,  n, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Numeric->Llen, n, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Numeric->Uip,  n, sizeof (Int), Common) ;
    TRILINOS_KLU_free (Numeric->Ulen, n, sizeof (Int), Common) ;

    TRILINOS_KLU_free (Numeric->LUsize, nblocks, sizeof (size_t), Common) ;

    TRILINOS_KLU_free (Numeric->LUbx, nblocks, sizeof (Unit *), Common) ;

    TRILINOS_KLU_free (Numeric->Udiag, n, sizeof (Entry), Common) ;

    TRILINOS_KLU_free (Numeric->Rs,   n, sizeof (double), Common) ;
    TRILINOS_KLU_free (Numeric->Pinv, n, sizeof (Int), Common) ;

    TRILINOS_KLU_free (Numeric->Work, Numeric->worksize, 1, Common) ;

    TRILINOS_KLU_free (Numeric, 1, sizeof (TRILINOS_KLU_numeric), Common) ;

    *NumericHandle = NULL ;
    return (TRUE) ;
}
