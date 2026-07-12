#define _GNU_SOURCE
/* ========================================================================== */
/* === TRILINOS_KLU_kernel =========================================================== */
/* ========================================================================== */

/* Sparse left-looking LU factorization, with partial pivoting.  Based on
 * Gilbert & Peierl's method, with a non-recursive DFS and with Eisenstat &
 * Liu's symmetric pruning.  No user-callable routines are in this file.
 */

/* This file should make the long int version of KLU */
#define DLONG 1

#include "trilinos_klu_internal.h"

/* ========================================================================== */
/* === dfs ================================================================== */
/* ========================================================================== */

/* Does a depth-first-search, starting at node j. */

#include <pthread.h>
#include <stdatomic.h>
#include <time.h>
#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#define KLS_KLU_CPU_RELAX() _mm_pause ()
#else
#define KLS_KLU_CPU_RELAX() do { } while (0)
#endif

#include <stdint.h>
#include <string.h>
/* units occupied by an int32 index array in the chunked column layout */
#define KLS_UNITS32(len) \
    (((size_t) (len) * sizeof (int32_t) + sizeof (Unit) - 1) / sizeof (Unit))

_Thread_local double kls_pipe_t_sym ;
_Thread_local double kls_pipe_t_num ;
_Thread_local long kls_pipe_copy_bytes ;
_Thread_local long kls_pipe_madds ;
static int kls_pipe_phase_prof ;
static int kls_pipe_dense_panels ;
Int KLS_SN_PANEL_FACTOR (double *A, Int R, Int W, Int rowids [ ],
			 const Int diagrows [ ], double tol) ;

_Thread_local double kls_construct_secs ;
_Thread_local long kls_construct_calls ;
_Thread_local long kls_construct_entries ;
static int kls_construct_prof ;

static double kls_klu_now (void)
{
    struct timespec ts ;
    clock_gettime (CLOCK_MONOTONIC, &ts) ;
    return ((double) ts.tv_sec + 1e-9 * (double) ts.tv_nsec) ;
}
#include <stdlib.h>

/* chunked mode: stable per-column storage; classic mode: base+offset */
#define KLS_COL_BASE(Colptr,LU,Lip,j) \
    ((Colptr) != NULL ? (Colptr)[j] : (LU) + (Lip)[j])

static Int dfs
(
    /* input, not modified on output: */
    Int j,		/* node at which to start the DFS */
    Int k,		/* mark value, for the Flag array */
    Int Pinv [ ],	/* Pinv [i] = k if row i is kth pivot row, or TRILINOS_KLU_EMPTY if
			 * row i is not yet pivotal.  */
    Int Llen [ ],	/* size n, Llen [k] = # nonzeros in column k of L */
    Int Lip [ ],	/* size n, Lip [k] is position in LU of column k of L */

    /* workspace, not defined on input or output */
    Int Stack [ ],	/* size n */

    /* input/output: */
    Int Flag [ ],	/* Flag [i] == k means i is marked */
    Int Lpend [ ],	/* for symmetric pruning */
    Int top,		/* top of stack on input*/
    Unit LU [],
    Unit *const *Colptr,	/* chunked per-column storage, or NULL */
    Int *Lik,		/* Li row index array of the kth column */
    Int *plength,

    /* other, not defined on input or output */
    Int Ap_pos [ ]	/* keeps track of position in adj list during DFS */
)
{
    Int i, pos, jnew, head, l_length ;
    Int *Li ;

    l_length = *plength ;

    head = 0 ;
    Stack [0] = j ;
    ASSERT (Flag [j] != k) ;

    while (head >= 0)
    {
	j = Stack [head] ;
	jnew = Pinv [j] ;
	ASSERT (jnew >= 0 && jnew < k) ;	/* j is pivotal */

	if (Flag [j] != k)	    /* a node is not yet visited */
	{
	    /* first time that j has been visited */
	    Flag [j] = k ;
	    PRINTF (("[ start dfs at %d : new %d\n", j, jnew)) ;
	    /* set Ap_pos [head] to one past the last entry in col j to scan */
	    Ap_pos [head] =
		(Lpend [jnew] == TRILINOS_KLU_EMPTY) ?  Llen [jnew] : Lpend [jnew] ;
	}

	/* add the adjacent nodes to the recursive stack by iterating through
	 * until finding another non-visited pivotal node */
	Li = (Int *) KLS_COL_BASE (Colptr, LU, Lip, jnew) ;
	for (pos = --Ap_pos [head] ; pos >= 0 ; --pos)
	{
	    i = Li [pos] ;
	    if (Flag [i] != k)
	    {
		/* node i is not yet visited */
		if (Pinv [i] >= 0)
		{
		    /* keep track of where we left off in the scan of the
		     * adjacency list of node j so we can restart j where we
		     * left off. */
		    Ap_pos [head] = pos ;

		    /* node i is pivotal; push it onto the recursive stack
		     * and immediately break so we can recurse on node i. */
		    Stack [++head] = i ;
		    break ;
		}
		else
		{
		    /* node i is not pivotal (no outgoing edges). */
		    /* Flag as visited and store directly into L,
		     * and continue with current node j. */
    	 	    Flag [i] = k ;
		    Lik [l_length] = i ;
	   	    l_length++ ;
		}
	    }
	}

	if (pos == -1)
	{
	    /* if all adjacent nodes of j are already visited, pop j from
	     * recursive stack and push j onto output stack */
	    head-- ;
	    Stack[--top] = j ;
	    PRINTF (("  end   dfs at %d ] head : %d\n", j, head)) ;
	}
    }

    *plength = l_length ;
    return (top) ;
}


/* ========================================================================== */
/* === lsolve_symbolic ====================================================== */
/* ========================================================================== */

/* Finds the pattern of x, for the solution of Lx=b */

static Int lsolve_symbolic
(
    /* input, not modified on output: */
    Int n,              /* L is n-by-n, where n >= 0 */
    Int k,		/* also used as the mark value, for the Flag array */
    Int Ap [ ],
    Int Ai [ ],
    Int Q [ ],
    Unit *const *Colptr,	/* chunked per-column storage, or NULL */
    Int Pinv [ ],	/* Pinv [i] = k if i is kth pivot row, or TRILINOS_KLU_EMPTY if row i
			 * is not yet pivotal.  */

    /* workspace, not defined on input or output */
    Int Stack [ ],	/* size n */

    /* workspace, defined on input and output */
    Int Flag [ ],	/* size n.  Initially, all of Flag [0..n-1] < k.  After
			 * lsolve_symbolic is done, Flag [i] == k if i is in
			 * the pattern of the output, and Flag [0..n-1] <= k. */

    /* other */
    Int Lpend [ ],	/* for symmetric pruning */
    Int Ap_pos [ ],	/* workspace used in dfs */

    Unit LU [ ],	/* LU factors (pattern and values) */
    Int lup,		/* pointer to free space in LU */
    Int Llen [ ],	/* size n, Llen [k] = # nonzeros in column k of L */
    Int Lip [ ],	/* size n, Lip [k] is position in LU of column k of L */

    /* ---- the following are only used in the TRILINOS_BTF case --- */

    Int k1,		/* the block of A is from k1 to k2-1 */
    Int PSinv [ ]	/* inverse of P from symbolic factorization */
)
{
    Int *Lik ;
    Int i, p, pend, oldcol, kglobal, top, l_length ;

    top = n ;
    l_length = 0 ;
    Lik = (Int *) (LU + lup);

    /* ---------------------------------------------------------------------- */
    /* TRILINOS_BTF factorization of A (k1:k2-1, k1:k2-1) */
    /* ---------------------------------------------------------------------- */

    kglobal = k + k1 ;	/* column k of the block is col kglobal of A */
    oldcol = Q [kglobal] ;	/* Q must be present for TRILINOS_BTF case */
    pend = Ap [oldcol+1] ;
    for (p = Ap [oldcol] ; p < pend ; p++)
    {
	i = PSinv [Ai [p]] - k1 ;
	if (i < 0) continue ;	/* skip entry outside the block */

	/* (i,k) is an entry in the block.  start a DFS at node i */
	PRINTF (("\n ===== DFS at node %d in b, inew: %d\n", i, Pinv [i])) ;
	if (Flag [i] != k)
	{
	    if (Pinv [i] >= 0)
	    {
		top = dfs (i, k, Pinv, Llen, Lip, Stack, Flag,
			   Lpend, top, LU, Colptr, Lik, &l_length, Ap_pos) ;
	    }
	    else
	    {
		/* i is not pivotal, and not flagged. Flag and put in L */
		Flag [i] = k ;
		Lik [l_length] = i ;
		l_length++;
	    }
	}
    }

    /* If Llen [k] is zero, the matrix is structurally singular */
    Llen [k] = l_length ;
    return (top) ;
}


/* ========================================================================== */
/* === construct_column ===================================================== */
/* ========================================================================== */

/* Construct the kth column of A, and the off-diagonal part, if requested.
 * Scatter the numerical values into the workspace X, and construct the
 * corresponding column of the off-diagonal matrix. */

static void construct_column
(
    /* inputs, not modified on output */
    Int k,	    /* the column of A (or the column of the block) to get */
    Int Ap [ ],
    Int Ai [ ],
    Entry Ax [ ],
    Int Q [ ],	    /* column pre-ordering */

    /* zero on input, modified on output */
    Entry X [ ],

    /* ---- the following are only used in the TRILINOS_BTF case --- */

    /* inputs, not modified on output */
    Int k1,	    /* the block of A is from k1 to k2-1 */
    Int PSinv [ ],  /* inverse of P from symbolic factorization */
    double Rs [ ],  /* scale factors for A */
    Int scale,	    /* 0: no scaling, nonzero: scale the rows with Rs */

    /* inputs, modified on output */
    Int Offp [ ],   /* off-diagonal matrix (modified by this routine) */
    Int Offi [ ],
    Entry Offx [ ]
)
{
    Entry aik ;
    Int i, p, pend, oldcol, kglobal, poff, oldrow ;
    const double kls_cc_t0 = kls_construct_prof ? kls_klu_now () : 0.0 ;

    /* ---------------------------------------------------------------------- */
    /* Scale and scatter the column into X. */
    /* ---------------------------------------------------------------------- */

    kglobal = k + k1 ;		/* column k of the block is col kglobal of A */
    poff = Offp [kglobal] ;	/* start of off-diagonal column */
    oldcol = Q [kglobal] ;
    pend = Ap [oldcol+1] ;

    if (scale <= 0)
    {
	/* no scaling */
	for (p = Ap [oldcol] ; p < pend ; p++)
	{
	    oldrow = Ai [p] ;
	    i = PSinv [oldrow] - k1 ;
	    aik = Ax [p] ;
	    if (i < 0)
	    {
		/* this is an entry in the off-diagonal part */
		Offi [poff] = oldrow ;
		Offx [poff] = aik ;
		poff++ ;
	    }
	    else
	    {
		/* (i,k) is an entry in the block.  scatter into X */
		X [i] = aik ;
	    }
	}
    }
    else
    {
	/* row scaling */
	for (p = Ap [oldcol] ; p < pend ; p++)
	{
	    oldrow = Ai [p] ;
	    i = PSinv [oldrow] - k1 ;
	    aik = Ax [p] ;
	    SCALE_DIV (aik, Rs [oldrow]) ;
	    if (i < 0)
	    {
		/* this is an entry in the off-diagonal part */
		Offi [poff] = oldrow ;
		Offx [poff] = aik ;
		poff++ ;
	    }
	    else
	    {
		/* (i,k) is an entry in the block.  scatter into X */
		X [i] = aik ;
	    }
	}
    }

    Offp [kglobal+1] = poff ;   /* start of the next col of off-diag part */
    if (kls_construct_prof)
    {
	kls_construct_secs += kls_klu_now () - kls_cc_t0 ;
	kls_construct_calls++ ;
	kls_construct_entries += pend - Ap [oldcol] ;
    }
}


/* ========================================================================== */
/* === lsolve_numeric ======================================================= */
/* ========================================================================== */

/* Computes the numerical values of x, for the solution of Lx=b.  Note that x
 * may include explicit zeros if numerical cancelation occurs.  L is assumed
 * to be unit-diagonal, with possibly unsorted columns (but the first entry in
 * the column must always be the diagonal entry). */

static void lsolve_numeric
(
    /* input, not modified on output: */
    Int Pinv [ ],	/* Pinv [i] = k if i is kth pivot row, or TRILINOS_KLU_EMPTY if row i
			 * is not yet pivotal.  */
    Unit *LU,		/* LU factors (pattern and values) */
    Unit *const *Colptr,	/* chunked per-column storage, or NULL */
    Int Stack [ ],	/* stack for dfs */
    Int Lip [ ],	/* size n, Lip [k] is position in LU of column k of L */
    Int top,		/* top of stack on input */
    Int n,		/* A is n-by-n */
    Int Llen [ ],	/* size n, Llen [k] = # nonzeros in column k of L */

    /* output, must be zero on input: */
    Entry X [ ]	/* size n, initially zero.  On output,
		 * X [Ui [up1..up-1]] and X [Li [lp1..lp-1]]
		 * contains the solution. */

)
{
    Entry xj ;
    Entry *Lx ;
    Int *Li ;
    Int p, s, j, jnew, len ;

    /* solve Lx=b */
    for (s = top ; s < n ; s++)
    {
	/* forward solve with column j of L */
 	j = Stack [s] ;
	jnew = Pinv [j] ;
	ASSERT (jnew >= 0) ;
	xj = X [j] ;
	{
	    Unit *xp = KLS_COL_BASE (Colptr, LU, Lip, jnew) ;
	    len = Llen [jnew] ;
	    Li = (Int *) xp ;
	    Lx = (Entry *) (xp + UNITS (Int, len)) ;
	}
	ASSERT (Lip [jnew] <= Lip [jnew+1]) ;
	for (p = 0 ; p < len ; p++)
	{
	    /*X [Li [p]] -= Lx [p] * xj ; */
	    MULT_SUB (X [Li [p]], Lx [p], xj) ;
	}
    }
}


/* ========================================================================== */
/* === lpivot =============================================================== */
/* ========================================================================== */

/* Find a pivot via partial pivoting, and scale the column of L. */

static Int lpivot
(
    Int diagrow,
    Int *p_pivrow,
    Entry *p_pivot,
    double *p_abs_pivot,
    double tol,
    Entry X [ ],
    Unit *LU,		/* LU factors (pattern and values) */
    Int Lip [ ],
    Int Llen [ ],
    Int k,
    Int n,

    Int Pinv [ ],	/* Pinv [i] = k if row i is kth pivot row, or TRILINOS_KLU_EMPTY if
			 * row i is not yet pivotal.  */

    Int *p_firstrow,
    TRILINOS_KLU_common *Common
)
{
    Entry x, pivot, *Lx ;
    double abs_pivot, xabs ;
    Int p, i, ppivrow, pdiag, pivrow, *Li, last_row_index, firstrow, len ;

    pivrow = TRILINOS_KLU_EMPTY ;
    if (Llen [k] == 0)
    {
	/* matrix is structurally singular */
	if (Common->halt_if_singular)
	{
	    return (FALSE) ;
	}
	for (firstrow = *p_firstrow ; firstrow < n ; firstrow++)
	{
	    PRINTF (("check %d\n", firstrow)) ;
	    if (Pinv [firstrow] < 0)
	    {
		/* found the lowest-numbered non-pivotal row.  Pick it. */
		pivrow = firstrow ;
		PRINTF (("Got pivotal row: %d\n", pivrow)) ;
		break ;
	    }
	}
	ASSERT (pivrow >= 0 && pivrow < n) ;
	CLEAR (pivot) ;
	*p_pivrow = pivrow ;
	*p_pivot = pivot ;
	*p_abs_pivot = 0 ;
	*p_firstrow = firstrow ;
	return (FALSE) ;
    }

    pdiag = TRILINOS_KLU_EMPTY ;
    ppivrow = TRILINOS_KLU_EMPTY ;
    abs_pivot = TRILINOS_KLU_EMPTY ;
    i = Llen [k] - 1 ;
    GET_POINTER (LU, Lip, Llen, Li, Lx, k, len) ;
    last_row_index = Li [i] ;

    /* decrement the length by 1 */
    Llen [k] = i ;
    GET_POINTER (LU, Lip, Llen, Li, Lx, k, len) ;

    /* look in Li [0 ..Llen [k] - 1 ] for a pivot row */
    for (p = 0 ; p < len ; p++)
    {
	/* gather the entry from X and store in L */
	i = Li [p] ;
	x = X [i] ;
	CLEAR (X [i]) ;

	Lx [p] = x ;
	/* xabs = ABS (x) ; */
	ABS (xabs, x) ;

	/* find the diagonal */
	if (i == diagrow)
	{
	    pdiag = p ;
	}

	/* find the partial-pivoting choice */
	if (xabs > abs_pivot)
	{
	    abs_pivot = xabs ;
	    ppivrow = p ;
	}
    }

    /* xabs = ABS (X [last_row_index]) ;*/
    ABS (xabs, X [last_row_index]) ;
    if (xabs > abs_pivot)
    {
        abs_pivot = xabs ;
        ppivrow = TRILINOS_KLU_EMPTY ;
    }

    /* compare the diagonal with the largest entry */
    if (last_row_index == diagrow)
    {
	if (xabs >= tol * abs_pivot)
	{
    	    abs_pivot = xabs ;
            ppivrow = TRILINOS_KLU_EMPTY ;
        }
    }
    else if (pdiag != TRILINOS_KLU_EMPTY)
    {
	/* xabs = ABS (Lx [pdiag]) ;*/
	ABS (xabs, Lx [pdiag]) ;
	if (xabs >= tol * abs_pivot)
	{
	    /* the diagonal is large enough */
	    abs_pivot = xabs ;
	    ppivrow = pdiag ;
	}
    }

    if (ppivrow != TRILINOS_KLU_EMPTY)
    {
        pivrow = Li [ppivrow] ;
        pivot  = Lx [ppivrow] ;
	/* overwrite the ppivrow values with last index values */
        Li [ppivrow] = last_row_index ;
        Lx [ppivrow] = X [last_row_index] ;
    }
    else
    {
        pivrow = last_row_index ;
        pivot = X [last_row_index] ;
    }
    CLEAR (X [last_row_index]) ;

    if (Common->kls_static_perturb > 0 && !IS_ZERO (pivot) &&
	abs_pivot < Common->kls_static_perturb)
    {
	/* hard-static small-pivot perturbation (SuperLU_DIST style):
	   bounded growth in exchange for a factorization of a nearby
	   matrix; the enrolled solve refinement recovers the contract */
	pivot = pivot < 0 ? -Common->kls_static_perturb
			  : Common->kls_static_perturb ;
	abs_pivot = Common->kls_static_perturb ;
	Common->kls_perturb_count++ ;
    }
    *p_pivrow = pivrow ;
    *p_pivot = pivot ;
    *p_abs_pivot = abs_pivot ;
    ASSERT (pivrow >= 0 && pivrow < n) ;

    if (IS_ZERO (pivot) && Common->halt_if_singular)
    {
	/* numerically singular case */
	return (FALSE) ;
    }

    /* divide L by the pivot value */
    for (p = 0 ; p < Llen [k] ; p++)
    {
	/* Lx [p] /= pivot ; */
	DIV (Lx [p], Lx [p], pivot) ;
    }

    return (TRUE) ;
}


/* ========================================================================== */
/* === prune ================================================================ */
/* ========================================================================== */

/* Prune the columns of L to reduce work in subsequent depth-first searches */
static void prune
(
    /* input/output: */
    Int Lpend [ ],	/* Lpend [j] marks symmetric pruning point for L(:,j) */

    /* input: */
    Int Pinv [ ],	/* Pinv [i] = k if row i is kth pivot row, or TRILINOS_KLU_EMPTY if
			 * row i is not yet pivotal.  */
    Int k,		/* prune using column k of U */
    Int pivrow,		/* current pivot row */

    /* input/output: */
    Unit *LU,		/* LU factors (pattern and values) */

    /* input */
    Int Uip [ ],	/* size n, column pointers for U */
    Int Lip [ ],	/* size n, column pointers for L */
    Int Ulen [ ],	/* size n, column length of U */
    Int Llen [ ]	/* size n, column length of L */
)
{
    Entry x ;
    Entry *Lx, *Ux ;
    Int *Li, *Ui ;
    Int p, i, j, p2, phead, ptail, llen, ulen ;

    /* check to see if any column of L can be pruned */
    GET_POINTER (LU, Uip, Ulen, Ui, Ux, k, ulen) ;
    for (p = 0 ; p < ulen ; p++)
    {
	j = Ui [p] ;
	ASSERT (j < k) ;
	PRINTF (("%d is pruned: %d. Lpend[j] %d Lip[j+1] %d\n",
	    j, Lpend [j] != TRILINOS_KLU_EMPTY, Lpend [j], Lip [j+1])) ;
	if (Lpend [j] == TRILINOS_KLU_EMPTY)
	{
	    /* scan column j of L for the pivot row */
            GET_POINTER (LU, Lip, Llen, Li, Lx, j, llen) ;
	    for (p2 = 0 ; p2 < llen ; p2++)
	    {
		if (pivrow == Li [p2])
		{
		    /* found it!  This column can be pruned */
#ifndef NDEBUG
		    PRINTF (("==== PRUNE: col j %d of L\n", j)) ;
		    {
			Int p3 ;
			for (p3 = 0 ; p3 < Llen [j] ; p3++)
			{
			    PRINTF (("before: %i  pivotal: %d\n", Li [p3],
					Pinv [Li [p3]] >= 0)) ;
			}
		    }
#endif

		    /* partition column j of L.  The unit diagonal of L
		     * is not stored in the column of L. */
		    phead = 0 ;
		    ptail = Llen [j] ;
		    while (phead < ptail)
		    {
			i = Li [phead] ;
			if (Pinv [i] >= 0)
			{
			    /* leave at the head */
			    phead++ ;
			}
			else
			{
			    /* swap with the tail */
			    ptail-- ;
			    Li [phead] = Li [ptail] ;
			    Li [ptail] = i ;
			    x = Lx [phead] ;
			    Lx [phead] = Lx [ptail] ;
			    Lx [ptail] = x ;
			}
		    }

		    /* set Lpend to one past the last entry in the
		     * first part of the column of L.  Entries in
		     * Li [0 ... Lpend [j]-1] are the only part of
		     * column j of L that needs to be scanned in the DFS.
		     * Lpend [j] was TRILINOS_KLU_EMPTY; setting it >= 0 also flags
		     * column j as pruned. */
		    Lpend [j] = ptail ;

#ifndef NDEBUG
		    {
			Int p3 ;
			for (p3 = 0 ; p3 < Llen [j] ; p3++)
			{
			    if (p3 == Lpend [j]) PRINTF (("----\n")) ;
			    PRINTF (("after: %i  pivotal: %d\n", Li [p3],
					Pinv [Li [p3]] >= 0)) ;
			}
		    }
#endif

		    break ;
		}
	    }
	}
    }
}


/* ========================================================================== */
/* === TRILINOS_KLU_kernel =========================================================== */
/* ========================================================================== */


/* ========================================================================== */
/* Diagonal-claim pivot for the parallel phase: identical selection to
 * lpivot, but if the threshold test would choose an off-diagonal pivot
 * (or the column is structurally singular), all effects are rolled back
 * and RETRY is returned so the column defers to the serial cleanup.
 * Distinct columns own distinct diagonal rows, so accepted claims are
 * collision-free without atomics. */
/* colptr-aware symmetric pruning: identical to prune() but addresses
 * columns through the chunked per-column pointers.  Callers must be
 * single-threaded phases (serial super-levels with the other workers
 * parked at the barrier, or the cleanup): pruning reorders finalized
 * columns in place. */
static void kls_prune_chunked
(
    Unit *const *Colptr,
    Int Lpend [ ],
    Int Pinv [ ],
    Int k,
    Int pivrow,
    Int Uip [ ],
    Int Lip [ ],
    Int Ulen [ ],
    Int Llen [ ]
)
{
    Entry x ;
    Entry *Lx, *Ux ;
    Int *Li, *Ui ;
    Int p, i, j, p2, phead, ptail, llen, ulen ;
    {
	Unit *xp = Colptr [k] + Uip [k] ;
	ulen = Ulen [k] ;
	Ui = (Int *) xp ;
	Ux = (Entry *) (xp + UNITS (Int, ulen)) ;
    }
    (void) Ux ;
    for (p = 0 ; p < ulen ; p++)
    {
	j = Ui [p] ;
	if (Lpend [j] == TRILINOS_KLU_EMPTY)
	{
	    Unit *xp = Colptr [j] + Lip [j] ;
	    llen = Llen [j] ;
	    Li = (Int *) xp ;
	    Lx = (Entry *) (xp + UNITS (Int, llen)) ;
	    for (p2 = 0 ; p2 < llen ; p2++)
	    {
		if (pivrow == Li [p2])
		{
		    phead = 0 ;
		    ptail = llen ;
		    while (phead < ptail)
		    {
			i = Li [phead] ;
			if (Pinv [i] >= 0)
			{
			    phead++ ;
			}
			else
			{
			    ptail-- ;
			    Li [phead] = Li [ptail] ;
			    Li [ptail] = i ;
			    x = Lx [phead] ;
			    Lx [phead] = Lx [ptail] ;
			    Lx [ptail] = x ;
			}
		    }
		    Lpend [j] = ptail ;
		    break ;
		}
	    }
	}
    }
}

/* Per-thread routing request from the solver layer: heavy-column
 * foreground first factors opt in; background/race threads see 0.
 * The KLS_KLU_PIPE env overrides for experiments. */
_Thread_local int kls_klu_pipe_threads = 0 ;

#define KLS_PIVOT_RETRY (-2)
static Int kls_lpivot_diag_claim
(
    Int diagrow,
    Int *p_pivrow,
    Entry *p_pivot,
    double *p_abs_pivot,
    double tol,
    Entry X [ ],
    Unit *LU,
    Int Lip [ ],
    Int Llen [ ],
    Int k,
    Int n
)
{
    Entry x, pivot, *Lx ;
    double abs_pivot, xabs ;
    Int p, i, pdiag, *Li, last_row_index, len ;

    if (Llen [k] == 0)
    {
	return (KLS_PIVOT_RETRY) ;   /* structurally singular: serial path */
    }
    i = Llen [k] - 1 ;
    GET_POINTER (LU, Lip, Llen, Li, Lx, k, len) ;
    last_row_index = Li [i] ;
    Llen [k] = i ;
    GET_POINTER (LU, Lip, Llen, Li, Lx, k, len) ;

    pdiag = TRILINOS_KLU_EMPTY ;
    abs_pivot = 0 ;
    for (p = 0 ; p < len ; p++)
    {
	i = Li [p] ;
	x = X [i] ;
	Lx [p] = x ;    /* gather (X preserved for possible rollback) */
	ABS (xabs, x) ;
	if (i == diagrow)
	{
	    pdiag = p ;
	}
	if (xabs > abs_pivot)
	{
	    abs_pivot = xabs ;
	}
    }
    ABS (xabs, X [last_row_index]) ;
    if (xabs > abs_pivot)
    {
	abs_pivot = xabs ;
    }

    if (last_row_index == diagrow)
    {
	ABS (xabs, X [last_row_index]) ;
	if (!(xabs >= tol * abs_pivot) || xabs == 0)
	{
	    /* the gather's Lx[0] clobbered Li[len] (same unit): restore */
	    Li [len] = last_row_index ;
	    Llen [k] = len + 1 ;
	    return (KLS_PIVOT_RETRY) ;
	}
	pivot = X [last_row_index] ;
	abs_pivot = xabs ;
	/* diagonal is the dropped last entry: keep L as gathered */
	for (p = 0 ; p < len ; p++)
	{
	    CLEAR (X [Li [p]]) ;
	}
	CLEAR (X [last_row_index]) ;
    }
    else if (pdiag != TRILINOS_KLU_EMPTY)
    {
	ABS (xabs, Lx [pdiag]) ;
	if (!(xabs >= tol * abs_pivot) || xabs == 0)
	{
	    /* the gather's Lx[0] clobbered Li[len] (same unit): restore */
	    Li [len] = last_row_index ;
	    Llen [k] = len + 1 ;
	    return (KLS_PIVOT_RETRY) ;
	}
	pivot = Lx [pdiag] ;
	abs_pivot = xabs ;
	/* swap the dropped last entry into the diagonal slot */
	Li [pdiag] = last_row_index ;
	Lx [pdiag] = X [last_row_index] ;
	for (p = 0 ; p < len ; p++)
	{
	    CLEAR (X [Li [p]]) ;
	}
	CLEAR (X [diagrow]) ;
	CLEAR (X [last_row_index]) ;
    }
    else
    {
	/* diagonal not present in the column */
	Li [len] = last_row_index ;
	Llen [k] = len + 1 ;
	return (KLS_PIVOT_RETRY) ;
    }

    for (p = 0 ; p < len ; p++)
    {
	DIV (Lx [p], Lx [p], pivot) ;
    }
    *p_pivrow = diagrow ;
    *p_pivot = pivot ;
    *p_abs_pivot = abs_pivot ;
    return (TRUE) ;
}

/* ========================================================================== */
/* === KLS per-column step API ============================================== */
/* ========================================================================== */

/* The factorization loop body, exposed one column at a time so an external
 * (level-scheduled, eventually parallel) driver can run it.  The classic
 * TRILINOS_KLU_kernel below is rewritten as init + loop{step} + finish over
 * this state, so every existing caller validates the extraction. */

typedef struct KLS_KLU_KERNEL_STATE_STRUCT
{
    Int n ;
    Int *Ap ;
    Int *Ai ;
    Entry *Ax ;
    Int *Q ;
    size_t lusize ;
    Int *Pinv ;
    Int *P ;
    Unit *LU ;
    Entry *Udiag ;
    Int *Llen ;
    Int *Ulen ;
    Int *Lip ;
    Int *Uip ;
    Int lnz ;
    Int unz ;
    Entry *X ;
    Int *Stack ;
    Int *Flag ;
    Int *Ap_pos ;
    Int *Lpend ;
    Int k1 ;
    Int *PSinv ;
    double *Rs ;
    Int *Offp ;
    Int *Offi ;
    Entry *Offx ;
    TRILINOS_KLU_common *Common ;
    Int firstrow ;
    Int lup ;
    Int scale ;
    double tol ;
    double memgrow ;
    Int no_prune ;   /* parallel drivers: finalized columns stay immutable */
    Int diag_claim ; /* parallel phase: only diagonal pivots; else defer */
    Int chunked_prune ; /* single-threaded phases may prune chunked cols */
    Int idx32 ;         /* chunked columns store int32 row/pivot indices
                           (pipe-published layout; KERNEL_FINISH widens
                           to Int while packing) */
    /* chunked storage (parallel drivers): stable per-column pointers,
       arena chunks never realloc once a column is published */
    Unit **colptr ;         /* size n, or NULL for classic base+offset */
    Unit *scratch ;         /* chunked: per-worker column build buffer
			       (size 2n+4 units); columns copy to exact-
			       size chunks at publish, so the arena holds
			       only real fill instead of the O(n^2) dense
			       bound */
    Int cols_done ;         /* columns completed (chunked pack bound) */
    Int pack_keep_row_indices ; /* early-return pack: match the classic
				   mid-factor layout (no Pinv convert) */
    Unit *chunk_head ;      /* current chunk (first Unit links to prev) */
    size_t chunk_used ;
    size_t chunk_size ;
} KLS_KLU_KERNEL_STATE ;

/* allocate nunits from the state's chunk arena (never moves memory) */
static Unit *kls_klu_chunk_alloc (KLS_KLU_KERNEL_STATE *S, size_t nunits)
{
    Unit *col ;
    if (S->chunk_head == NULL || S->chunk_used + nunits > S->chunk_size)
    {
	size_t want = S->chunk_size > 0 ? S->chunk_size : (size_t) 1 << 18 ;
	Unit *chunk ;
	while (want < nunits + 1) want *= 2 ;
	/* raw malloc: concurrent workers must not touch Common's
	   memusage counters */
	chunk = (Unit *) malloc (want * sizeof (Unit)) ;
	if (chunk == NULL)
	{
	    return (NULL) ;
	}
	*((Unit **) chunk) = S->chunk_head ;   /* link previous */
	S->chunk_head = chunk ;
	S->chunk_used = 1 ;
	S->chunk_size = want ;
    }
    col = S->chunk_head + S->chunk_used ;
    S->chunk_used += nunits ;
    return (col) ;
}

void KLS_KLU_KERNEL_CHUNKS_FREE (KLS_KLU_KERNEL_STATE *S)
{
    Unit *chunk = S->chunk_head ;
    while (chunk != NULL)
    {
	Unit *prev = *((Unit **) chunk) ;
	free (chunk) ;
	chunk = prev ;
    }
    S->chunk_head = NULL ;
    S->chunk_used = 0 ;
}

void KLS_KLU_KERNEL_INIT
(
    KLS_KLU_KERNEL_STATE *S
)
{
    Int k ;
    S->cols_done = 0 ;
    S->firstrow = 0 ;
    S->lup = 0 ;
    S->lnz = 0 ;
    S->unz = 0 ;
    for (k = 0 ; k < S->n ; k++)
    {
	CLEAR (S->X [k]) ;
	S->Flag [k] = TRILINOS_KLU_EMPTY ;
	S->Lpend [k] = TRILINOS_KLU_EMPTY ;
    }
    for (k = 0 ; k < S->n ; k++)
    {
	S->P [k] = k ;
	S->Pinv [k] = FLIP (k) ;
    }
    S->Offp [0] = 0 ;
}

/* returns 0 on success, -1 on fatal (status set), 1 on singular-halt */
Int KLS_KLU_KERNEL_STEP
(
    KLS_KLU_KERNEL_STATE *S,
    Int k
)
{
    Entry pivot ;
    double abs_pivot, xsize, nunits ;
    Entry *Ux ;
    Int *Li, *Ui ;
    Unit *LU = S->LU ;
    Int p, i, j, pivrow = TRILINOS_KLU_EMPTY, kbar, diagrow, top, len ;
    size_t newlusize ;
    const Int n = S->n ;
    Int *Pinv = S->Pinv ;
    Int *Llen = S->Llen ;
    Int *Lip = S->Lip ;
    TRILINOS_KLU_common *Common = S->Common ;

    nunits = DUNITS (Int, n - k) + DUNITS (Int, k) +
	     DUNITS (Entry, n - k) + DUNITS (Entry, k) ;
    if (S->colptr != NULL)
    {
	/* chunked: build the column in the worker's scratch buffer and
	   copy it to an exact-size chunk at publish time (the dense
	   bound would make the arena O(n^2) on large blocks) */
	(void) nunits ;
	Lip [k] = 0 ;
	LU = S->scratch ;
    }
    else
    {
    xsize = ((double) S->lup) + nunits ;
    if (xsize > (double) S->lusize)
    {
	xsize = (S->memgrow * ((double) S->lusize) + 4*n + 1) ;
	if (INT_OVERFLOW (xsize))
	{
	    Common->status = TRILINOS_KLU_TOO_LARGE ;
	    return (-1) ;
	}
	newlusize = S->memgrow * S->lusize + 2*n + 1 ;
	LU = (Unit*) TRILINOS_KLU_realloc (newlusize, S->lusize, sizeof (Unit),
					   LU, Common) ;
	Common->nrealloc++ ;
	S->LU = LU ;
	if (Common->status == TRILINOS_KLU_OUT_OF_MEMORY)
	{
	    return (-1) ;
	}
	S->lusize = newlusize ;
    }

    Lip [k] = S->lup ;
    }

    top = lsolve_symbolic (n, k, S->Ap, S->Ai, S->Q, S->colptr, Pinv,
		S->Stack, S->Flag, S->Lpend, S->Ap_pos, LU,
		S->colptr != NULL ? 0 : S->lup, Llen, Lip, S->k1, S->PSinv) ;

    construct_column (k, S->Ap, S->Ai, S->Ax, S->Q, S->X,
	S->k1, S->PSinv, S->Rs, S->scale, S->Offp, S->Offi, S->Offx) ;

    lsolve_numeric (Pinv, LU, S->colptr, S->Stack, Lip, top, n, Llen,
		    S->X) ;

    diagrow = S->P [k] ;

    if (S->diag_claim)
    {
	Int claim = kls_lpivot_diag_claim (diagrow, &pivrow, &pivot,
					   &abs_pivot, S->tol, S->X, LU, Lip,
					   Llen, k, n) ;
	if (claim == KLS_PIVOT_RETRY)
	{
	    /* roll back: clear the U-part scatter (L part is preserved by
	       the claim's own rollback) and defer to the serial cleanup */
	    for (p = top ; p < n ; p++)
	    {
		CLEAR (S->X [S->Stack [p]]) ;
	    }
	    {
		/* the column is not published yet: it lives in the local
		   build buffer (scratch in chunked mode) */
		Int *Lik = (Int *) (LU + Lip [k]) ;
		for (p = 0 ; p < Llen [k] ; p++)
		{
		    CLEAR (S->X [Lik [p]]) ;
		}
	    }
	    return (2) ;
	}
	/* accepted diagonal claim: distinct rows, no race */
	S->Uip [k] = Lip [k] + UNITS (Int, Llen [k]) +
		     UNITS (Entry, Llen [k]) ;
	S->Ulen [k] = n - top ;
	GET_POINTER (LU, S->Uip, S->Ulen, Ui, Ux, k, len) ;
	for (p = top, i = 0 ; p < n ; p++, i++)
	{
	    j = S->Stack [p] ;
	    Ui [i] = Pinv [j] ;
	    Ux [i] = S->X [j] ;
	    CLEAR (S->X [j]) ;
	}
	S->Udiag [k] = pivot ;
	{
	    /* publish: exact-size copy out of the scratch buffer */
	    size_t used = (size_t) (UNITS (Int, Llen [k]) +
				    UNITS (Entry, Llen [k]) +
				    UNITS (Int, S->Ulen [k]) +
				    UNITS (Entry, S->Ulen [k])) ;
	    Unit *col = kls_klu_chunk_alloc (S, used) ;
	    if (col == NULL)
	    {
		Common->status = TRILINOS_KLU_OUT_OF_MEMORY ;
		return (-1) ;
	    }
	    memcpy (col, S->scratch, used * sizeof (Unit)) ;
	    S->colptr [k] = col ;
	}
	S->P [k] = diagrow ;
	Pinv [diagrow] = k ;
	if (S->chunked_prune)
	{
	    kls_prune_chunked (S->colptr, S->Lpend, Pinv, k, diagrow,
			       S->Uip, Lip, S->Ulen, Llen) ;
	}
	S->lnz += Llen [k] + 1 ;
	S->unz += S->Ulen [k] + 1 ;
	S->cols_done = k + 1 ;
	return (0) ;
    }

    if (!lpivot (diagrow, &pivrow, &pivot, &abs_pivot, S->tol, S->X, LU, Lip,
		Llen, k, n, Pinv, &S->firstrow, Common))
    {
	Common->status = TRILINOS_KLU_SINGULAR ;
	if (Common->numerical_rank == TRILINOS_KLU_EMPTY)
	{
	    Common->numerical_rank = k + S->k1 ;
	    Common->singular_col = S->Q [k + S->k1] ;
	}
	if (Common->halt_if_singular)
	{
	    return (1) ;
	}
    }

    S->Uip [k] = Lip [k] + UNITS (Int, Llen [k]) + UNITS (Entry, Llen [k]) ;
    if (S->colptr == NULL)
    {
	S->lup += UNITS (Int, Llen [k]) + UNITS (Entry, Llen [k]) ;
    }
    S->Ulen [k] = n - top ;

    GET_POINTER (LU, S->Uip, S->Ulen, Ui, Ux, k, len) ;
    {
	/* accumulate the TRILINOS_KLU_flops sum while the U pattern is
	   hot: Llen of every source column is final once it pivots */
	double kls_fl = 0.0 ;
	for (p = top, i = 0 ; p < n ; p++, i++)
	{
	    j = S->Stack [p] ;
	    Ui [i] = Pinv [j] ;
	    Ux [i] = S->X [j] ;
	    CLEAR (S->X [j]) ;
	    kls_fl += (double) Llen [Ui [i]] ;
	}
	Common->kls_kernel_flops += 2.0 * kls_fl + (double) Llen [k] ;
    }
    if (S->colptr == NULL)
    {
	S->lup += UNITS (Int, S->Ulen [k]) + UNITS (Entry, S->Ulen [k]) ;
    }

    S->Udiag [k] = pivot ;

    if (pivrow != diagrow)
    {
	Common->noffdiag++ ;
	if (Pinv [diagrow] < 0)
	{
	    kbar = FLIP (Pinv [pivrow]) ;
	    S->P [kbar] = diagrow ;
	    Pinv [diagrow] = FLIP (kbar) ;
	}
    }
    S->P [k] = pivrow ;
    Pinv [pivrow] = k ;

    if (!S->no_prune && S->colptr == NULL)
    {
	prune (S->Lpend, Pinv, k, pivrow, LU, S->Uip, Lip, S->Ulen, Llen) ;
    }

    if (S->colptr != NULL)
    {
	Int do_prune_after = S->chunked_prune ;
	(void) do_prune_after ;
	/* publish: exact-size copy out of the scratch buffer */
	size_t used = (size_t) (UNITS (Int, Llen [k]) +
				UNITS (Entry, Llen [k]) +
				UNITS (Int, S->Ulen [k]) +
				UNITS (Entry, S->Ulen [k])) ;
	Unit *col = kls_klu_chunk_alloc (S, used) ;
	if (col == NULL)
	{
	    Common->status = TRILINOS_KLU_OUT_OF_MEMORY ;
	    return (-1) ;
	}
	memcpy (col, S->scratch, used * sizeof (Unit)) ;
	S->colptr [k] = col ;
	if (S->chunked_prune)
	{
	    kls_prune_chunked (S->colptr, S->Lpend, Pinv, k, pivrow,
			       S->Uip, Lip, S->Ulen, Llen) ;
	}
    }

    S->lnz += Llen [k] + 1 ;
    S->unz += S->Ulen [k] + 1 ;
    S->cols_done = k + 1 ;
    return (0) ;
}

size_t KLS_KLU_KERNEL_FINISH   /* returns final LU size */
(
    KLS_KLU_KERNEL_STATE *S
)
{
    Int p, i ;
    Int *Li ;
    size_t newlusize ;
    if (S->colptr != NULL)
    {
	/* pack the chunked columns into one contiguous LU in column
	   order, rebasing Lip/Uip and converting L row indices to
	   pivotal indices */
	size_t total = 0, off ;
	Unit *packed ;
	const Int ncols = S->cols_done ;
	for (p = 0 ; p < ncols ; p++)
	{
	    if (S->colptr [p] == NULL)
	    {
		continue ;   /* never processed (parallel failure path) */
	    }
	    total += UNITS (Int, S->Llen [p]) + UNITS (Entry, S->Llen [p]) +
		     UNITS (Int, S->Ulen [p]) + UNITS (Entry, S->Ulen [p]) ;
	}
	packed = (Unit *) TRILINOS_KLU_realloc (total, S->lusize,
						sizeof (Unit), S->LU,
						S->Common) ;
	if (packed == NULL || S->Common->status == TRILINOS_KLU_OUT_OF_MEMORY)
	{
	    return (S->lusize) ;
	}
	off = 0 ;
	for (p = 0 ; p < ncols ; p++)
	{
	    size_t lunits ;
	    if (S->colptr [p] == NULL)
	    {
		S->Lip [p] = 0 ;
		S->Uip [p] = 0 ;
		S->Llen [p] = 0 ;
		S->Ulen [p] = 0 ;
		continue ;
	    }
	    lunits = UNITS (Int, S->Llen [p]) +
			    UNITS (Entry, S->Llen [p]) ;
	    size_t uunits = UNITS (Int, S->Ulen [p]) +
			    UNITS (Entry, S->Ulen [p]) ;
	    if (S->idx32)
	    {
		/* widen the pipe's int32 column to the standard Int
		   layout while packing (converting L rows to pivotal
		   indices in the same pass) */
		const Unit *src = S->colptr [p] ;
		const int32_t *sLi = (const int32_t *) src ;
		const Entry *sLx =
		    (const Entry *) (src + KLS_UNITS32 (S->Llen [p])) ;
		const size_t slpart =
		    KLS_UNITS32 (S->Llen [p]) + UNITS (Entry, S->Llen [p]) ;
		const int32_t *sUi = (const int32_t *) (src + slpart) ;
		const Entry *sUx =
		    (const Entry *) (src + slpart +
				     KLS_UNITS32 (S->Ulen [p])) ;
		Int *dLi = (Int *) (packed + off) ;
		Entry *dLx =
		    (Entry *) (packed + off + UNITS (Int, S->Llen [p])) ;
		Int *dUi = (Int *) (packed + off + lunits) ;
		Entry *dUx = (Entry *) (packed + off + lunits +
					UNITS (Int, S->Ulen [p])) ;
		for (i = 0 ; i < S->Llen [p] ; i++)
		{
		    dLi [i] = S->pack_keep_row_indices
			? (Int) sLi [i] : S->Pinv [sLi [i]] ;
		    dLx [i] = sLx [i] ;
		}
		for (i = 0 ; i < S->Ulen [p] ; i++)
		{
		    dUi [i] = (Int) sUi [i] ;
		    dUx [i] = sUx [i] ;
		}
		S->Lip [p] = (Int) off ;
		S->Uip [p] = (Int) (off + lunits) ;
		off += lunits + uunits ;
		continue ;
	    }
	    memcpy (packed + off, S->colptr [p],
		    (lunits + uunits) * sizeof (Unit)) ;
	    S->Lip [p] = (Int) off ;
	    S->Uip [p] = (Int) (off + lunits) ;
	    if (!S->pack_keep_row_indices)
	    {
		Li = (Int *) (packed + off) ;
		for (i = 0 ; i < S->Llen [p] ; i++)
		{
		    Li [i] = S->Pinv [Li [i]] ;
		}
	    }
	    off += lunits + uunits ;
	}
	KLS_KLU_KERNEL_CHUNKS_FREE (S) ;
	S->LU = packed ;
	S->lusize = total ;
	return (total) ;
    }
    for (p = 0 ; p < S->n ; p++)
    {
	Li = (Int *) (S->LU + S->Lip [p]) ;
	for (i = 0 ; i < S->Llen [p] ; i++)
	{
	    Li [i] = S->Pinv [Li [i]] ;
	}
    }
    newlusize = S->lup ;
    S->LU = (Unit*) TRILINOS_KLU_realloc (newlusize, S->lusize, sizeof (Unit),
					  S->LU, S->Common) ;
    S->lusize = newlusize ;
    return (newlusize) ;
}


/* ========================================================================== */
/* === KLS level-scheduled kernel =========================================== */
/* ========================================================================== */

/* Column elimination tree of the block view (Liu's algorithm on A'A,
 * CSparse cs_etree ata-style): bounds every pivoting-time dependency
 * (George/Ng), so a level schedule over it stays valid regardless of the
 * runtime pivot choices. */
static int kls_klu_block_coletree
(
    Int n,
    Int Ap [ ],
    Int Ai [ ],
    Int Q [ ],
    Int k1,
    Int PSinv [ ],
    Int parent [ ],     /* size n out */
    Int ancestor [ ],   /* size n workspace */
    Int prev [ ]        /* size n workspace (per block row) */
)
{
    Int k, p, i, jroot, next ;
    for (k = 0 ; k < n ; k++)
    {
	parent [k] = TRILINOS_KLU_EMPTY ;
	ancestor [k] = TRILINOS_KLU_EMPTY ;
	prev [k] = TRILINOS_KLU_EMPTY ;
    }
    for (k = 0 ; k < n ; k++)
    {
	Int oldcol = Q [k + k1] ;
	for (p = Ap [oldcol] ; p < Ap [oldcol+1] ; p++)
	{
	    i = PSinv [Ai [p]] - k1 ;
	    if (i < 0 || i >= n)
	    {
		continue ;      /* off-block entry */
	    }
	    jroot = prev [i] ;
	    prev [i] = k ;
	    while (jroot != TRILINOS_KLU_EMPTY && jroot < k)
	    {
		next = ancestor [jroot] ;
		ancestor [jroot] = k ;
		if (next == TRILINOS_KLU_EMPTY)
		{
		    parent [jroot] = k ;
		    break ;
		}
		jroot = next ;
	    }
	}
    }
    return (1) ;
}

typedef struct kls_klu_par_shared_s
{
    Int *level_ptr ;
    Int *level_cols ;
    Int nlevels ;
    Int *super_ptr ;     /* super-level boundaries (level indices) */
    char *super_serial ; /* 1: narrow batch, tid 0 runs it alone */
    Int nsuper ;
    Int n ;
    Int *parent ;
    _Atomic char *defer ;
    _Atomic int abort_flag ;
    _Atomic long defer_count ;
    long defer_limit ;   /* heavy deferral: not a diagonal-claimable
			    matrix; abort and let the caller retry
			    serially instead of grinding the no-prune
			    cleanup (mac_econ: 37s wasted attempt) */
    pthread_barrier_t barrier ;
    int nthreads ;
    int prof ;           /* KLS_KLU_LEVELS_PROF: worker 0 buckets wall */
    double wide_secs, serial_secs ;
} kls_klu_par_shared ;

typedef struct kls_klu_par_worker_s
{
    KLS_KLU_KERNEL_STATE S ;
    kls_klu_par_shared *sh ;
    int tid ;
} kls_klu_par_worker ;

static void *kls_klu_par_worker_main (void *arg)
{
    kls_klu_par_worker *W = (kls_klu_par_worker *) arg ;
    kls_klu_par_shared *sh = W->sh ;
    Int sl, pos, k, r, j ;
    for (sl = 0 ; sl < sh->nsuper ; sl++)
    {
	const Int pos0 = sh->level_ptr [sh->super_ptr [sl]] ;
	const Int pos1 = sh->level_ptr [sh->super_ptr [sl + 1]] ;
	const int serial = sh->super_serial [sl] ;
	const double t0 = (sh->prof && W->tid == 0) ? kls_klu_now () : 0 ;
	if (serial && W->tid != 0)
	{
	    pthread_barrier_wait (&sh->barrier) ;
	    continue ;
	}
	W->S.chunked_prune = serial ;
	if (!atomic_load_explicit (&sh->abort_flag, memory_order_acquire))
	{
	    for (pos = serial ? pos0 : pos0 + W->tid ;
		 pos < pos1 ; pos += serial ? 1 : sh->nthreads)
	    {
		k = sh->level_cols [pos] ;
		if (atomic_load_explicit (&sh->defer [k],
					  memory_order_acquire))
		{
		    continue ;
		}
		r = KLS_KLU_KERNEL_STEP (&W->S, k) ;
		if (r == 2)
		{
		    long marked = 0 ;
		    /* defer this column and its etree ancestor path:
		       every dependent of k is one of k's ancestors */
		    atomic_store_explicit (&sh->defer [k], 1,
					   memory_order_release) ;
		    marked++ ;
		    for (j = sh->parent [k] ; j != TRILINOS_KLU_EMPTY ;
			 j = sh->parent [j])
		    {
			if (atomic_load_explicit (&sh->defer [j],
						  memory_order_relaxed))
			{
			    break ;
			}
			atomic_store_explicit (&sh->defer [j], 1,
					       memory_order_release) ;
			marked++ ;
		    }
		    if (atomic_fetch_add_explicit (&sh->defer_count, marked,
						   memory_order_relaxed) +
			marked > sh->defer_limit)
		    {
			atomic_store_explicit (&sh->abort_flag, 1,
					       memory_order_release) ;
			break ;
		    }
		}
		else if (r != 0)
		{
		    atomic_store_explicit (&sh->abort_flag, 1,
					   memory_order_release) ;
		    break ;
		}
	    }
	}
	pthread_barrier_wait (&sh->barrier) ;
	if (sh->prof && W->tid == 0)
	{
	    double dt = kls_klu_now () - t0 ;
	    if (serial) { sh->serial_secs += dt ; }
	    else        { sh->wide_secs += dt ; }
	}
    }
    return (NULL) ;
}

/* ========================================================================== */
/* === pipelined pivoted factorization ===================================== */
/* ========================================================================== */

/* Software-pipelined left-looking factorization with real partial
 * pivoting.  Columns are claimed in order; PUBLICATION IS IN ORDER (a
 * column pivots only once every earlier column has published), which
 * keeps pivot choices collision-free and makes the final phase of each
 * column exactly the serial algorithm.  The overlap is the update
 * pre-accumulation: while waiting for its predecessors, a column
 * eliminates through already-published pivot columns in rounds.  Round
 * r eliminates pivots in [P_{r-1}, P_r); contributions to a pivotal
 * row only ever come from strictly smaller pivot columns, so per-round
 * immediate elimination and U extraction preserve serial semantics.
 * No pruning (concurrent readers), as in the level-scheduled kernel. */

typedef struct kls_klu_pipe_shared_s
{
    _Atomic Int next_col ;
    _Atomic Int prefix ;      /* all columns < prefix are published */
    _Atomic int abort_flag ;
    Int n ;
    int nthreads ;
    /* concurrent-safe symmetric pruning (KLS_KLU_PIPE_PRUNE): finals
       are serialized by in-order publication, so pruning is single-
       writer; readers use per-column seqlocks.  Lpend[j] != EMPTY
       bounds the DFS scan of column j to its pivotal prefix. */
    _Atomic Int *lpend ;
    _Atomic unsigned *colver ;
    /* panel claiming (supernodal-pipe step a): workers claim whole
       fundamental supernodes; the per-column processing and in-order
       publication are unchanged, so results are bit-identical to
       column claiming - panels only group the claim order. */
    _Atomic Int next_panel ;
    const Int *panel_start ;  /* npanels+1 offsets, or NULL */
    Int npanels ;
} kls_klu_pipe_shared ;

typedef struct kls_klu_pipe_worker_s
{
    KLS_KLU_KERNEL_STATE S ;
    kls_klu_pipe_shared *sh ;
    Int *ubuf_i ;             /* per-round U pattern accumulation */
    Entry *ubuf_x ;
    double kflops ;           /* assembly-time flop sum of the columns
                                 this worker published */
    unsigned *ap_ver ;        /* seqlock snapshots for resumable scans */
    Unit *copybuf ;           /* seqlock column copies for the numeric */
    int tid ;
    double t_work, t_spin, t_final ;   /* KLS_KLU_PIPE_PROF buckets */
    long n_cols, n_rounds ;
    /* supernodal-pipe step (b): per-panel multi-column workspace.
       pX[w] is the wth panel column's dense accumulator (pX[0] aliases
       S.X); pLik/pUbuf_i/pUbuf_x are the per-column pattern and
       U-segment builds; pW is the panel width cap this worker can
       afford (memory-bounded at allocation). */
    Entry *pX [32] ;
    Int *pLik [32] ;
    Int *pUbuf_i [32] ;
    Entry *pUbuf_x [32] ;
    Int *pPFlag ;             /* union dedup marks (panel generation) */
    Int *pUlist ;             /* union source list (the DFS owns Stack) */
    Int pColCap ;             /* per-column pattern/ubuf capacity */
    Int *pFlagW [32] ;        /* per-column DFS marks: sibling columns'
                                 interleaved collections must not
                                 overwrite each other's marks */
    /* panel-major dense accumulator: B[rowpos*W + w] keeps the panel's
       numeric updates contiguous in w (one-two cachelines per source
       entry instead of W scattered full-size vectors) */
    Entry *pB ;
    size_t pBcap ;            /* rows currently allocated */
    Int *pRowPos ;            /* row -> buffer position (gen-tagged) */
    Int *pRowGen ;
    Int *pRowList ;           /* buffer position -> row (dense finalize) */
    Int pGen ;
    Int pW ;
} kls_klu_pipe_worker ;

/* DFS from node i (row of the block), eliminating only through pivot
 * columns < plimit; rows pivotal at >= plimit (or not pivotal) are L
 * candidates.  Identical in shape to dfs() but with the prefix filter.
 * With sh non-NULL, scans are bounded by the shared Lpend and guarded
 * by per-column seqlocks: a scan (or a resume after descending) whose
 * column version changed restarts that column's scan - marked rows
 * dedup, so restarts only recover misses.  Returns the new top;
 * Stack[top..oldtop-1] is the round's topological segment. */
static Int kls_pipe_dfs
(
    Int j,                    /* node at which to start the DFS */
    Int k,                    /* mark value for Flag */
    Int plimit,               /* eliminate only through pivots < plimit */
    Int Pinv [ ],
    Int Llen [ ],
    Int Lip [ ],
    Int Stack [ ],
    Int Flag [ ],
    Int top,
    Unit LU [ ],
    Unit *const *Colptr,
    Int *Lik,
    Int *plength,
    Int Ap_pos [ ],
    const kls_klu_pipe_shared *sh,
    unsigned Ap_ver [ ]
)
{
    Int i, pos, jnew, head, l_length ;
    Int *Li ;

    l_length = *plength ;

    head = 0 ;
    Stack [0] = j ;
    ASSERT (Flag [j] != k) ;

    while (head >= 0)
    {
	j = Stack [head] ;
	jnew = Pinv [j] ;
	if (Flag [j] != k)
	{
	    /* first time j has been visited */
	    Flag [j] = k ;
	    /* set Ap_pos [head] to one past the last entry in col j to scan */
	    if (jnew < 0 || jnew >= plimit)
	    {
		Ap_pos [head] = 0 ;
	    }
	    else if (sh != NULL && sh->lpend != NULL)
	    {
		unsigned v ;
		Int bound ;
		do
		{
		    v = atomic_load_explicit (&sh->colver [jnew],
					      memory_order_acquire) ;
		    KLS_KLU_CPU_RELAX () ;
		} while (v & 1u) ;
		bound = atomic_load_explicit (&sh->lpend [jnew],
					      memory_order_acquire) ;
		Ap_pos [head] = (bound == TRILINOS_KLU_EMPTY)
		    ? Llen [jnew] : bound ;
		Ap_ver [head] = v ;
	    }
	    else
	    {
		Ap_pos [head] = Llen [jnew] ;
	    }
	}
	else if (sh != NULL && sh->lpend != NULL &&
		 jnew >= 0 && jnew < plimit &&
		 atomic_load_explicit (&sh->colver [jnew],
				       memory_order_acquire) !=
		   Ap_ver [head])
	{
	    /* resumed a column whose rows were pruned mid-descent:
	       restart its scan under the new version (dedup recovers) */
	    unsigned v ;
	    Int bound ;
	    do
	    {
		v = atomic_load_explicit (&sh->colver [jnew],
					  memory_order_acquire) ;
		KLS_KLU_CPU_RELAX () ;
	    } while (v & 1u) ;
	    bound = atomic_load_explicit (&sh->lpend [jnew],
					  memory_order_acquire) ;
	    Ap_pos [head] = (bound == TRILINOS_KLU_EMPTY)
		? Llen [jnew] : bound ;
	    Ap_ver [head] = v ;
	}

	/* add the adjacent nodes to the recursive stack by iterating through
	 * until finding another non-visited pivotal node */
	if (jnew >= 0 && jnew < plimit)
	{
	    const int kls_pipe_cols = sh != NULL && sh->lpend != NULL ;
	    Unit *xp = kls_pipe_cols
		? (Unit *) __atomic_load_n ((Unit **) &Colptr [jnew],
					    __ATOMIC_ACQUIRE)
		: KLS_COL_BASE (Colptr, LU, Lip, jnew) ;
	    const int32_t *Li32 = (const int32_t *) xp ;
	    Li = (Int *) xp ;
	    for (pos = --Ap_pos [head] ; pos >= 0 ; --pos)
	    {
		i = kls_pipe_cols ? (Int) Li32 [pos] : Li [pos] ;
		if (Flag [i] != k)
		{
		    /* node i is not yet visited */
		    if (Pinv [i] >= 0 && Pinv [i] < plimit)
		    {
			/* keep track of where we left off in the scan of the
			 * adjacency list of node j so we can restart j where we
			 * left off. */
			Ap_pos [head] = pos ;

			/* node i is pivotal; push onto recursive stack */
			Stack [++head] = i ;
			break ;
		    }
		    else
		    {
			/* i is an L candidate for column k of the block */
			Flag [i] = k ;
			Lik [l_length] = i ;
			l_length++ ;
		    }
		}
	    }
	    if (pos >= 0)
	    {
		continue ;
	    }
	}

	/* pop off the recursive stack and push j in the output stack */
	if (sh != NULL && sh->lpend != NULL && jnew >= 0 && jnew < plimit &&
	    atomic_load_explicit (&sh->colver [jnew],
				  memory_order_acquire) != Ap_ver [head])
	{
	    /* the column was pruned during its final partial scan:
	       rescan under the new version before retiring it */
	    continue ;
	}
	head-- ;
	if (jnew >= 0 && jnew < plimit)
	{
	    Stack [--top] = j ;
	}
	else
	{
	    /* j itself is an L candidate (promoted-root case is handled by
	       the caller; roots reaching here were classified above) */
	}
    }

    *plength = l_length ;
    return (top) ;
}

static void kls_pipe_lsolve_numeric
(
    Int Pinv [ ], Unit *LU, Unit *const *Colptr, Int Stack [ ],
    Int Lip [ ], Int top, Int n, Int Llen [ ], Entry X [ ],
    const kls_klu_pipe_shared *sh, Unit *copybuf
) ;

/* One elimination round for column k of the block: symbolic (from the
 * given roots) + numeric + U-segment extraction, eliminating only
 * through pivots < plimit.  root_rows==NULL means round 1 (roots are
 * A(:,k)'s block rows); otherwise roots are promoted former candidates
 * (already flagged, already pivotal).  Returns 0 on success. */
static void kls_pipe_round
(
    KLS_KLU_KERNEL_STATE *S,
    Int k,
    Int plimit,
    const Int *root_rows,
    Int nroots,
    Unit *LU,                 /* local scratch base */
    Int *Lik,
    Int *plength,
    Int *ubuf_i,
    Entry *ubuf_x,
    Int *pucount,
    const kls_klu_pipe_shared *sh,
    unsigned *Ap_ver,
    Unit *copybuf
)
{
    const Int n = S->n ;
    Int *Pinv = S->Pinv ;
    Int top = n ;
    Int r, i, p, s ;
    const double kls_rd_t0 = kls_pipe_phase_prof ? kls_klu_now () : 0.0 ;

    if (root_rows == NULL)
    {
	/* round 1: roots are the block rows of A(:,k) */
	Int kglobal = k + S->k1 ;
	Int oldcol = S->Q [kglobal] ;
	Int pend = S->Ap [oldcol+1] ;
	for (p = S->Ap [oldcol] ; p < pend ; p++)
	{
	    i = S->PSinv [S->Ai [p]] - S->k1 ;
	    if (i < 0) continue ;
	    if (S->Flag [i] != k)
	    {
		if (Pinv [i] >= 0 && Pinv [i] < plimit)
		{
		    top = kls_pipe_dfs (i, k, plimit, Pinv, S->Llen, S->Lip,
					S->Stack, S->Flag, top, LU, S->colptr,
					Lik, plength, S->Ap_pos, sh, Ap_ver) ;
		}
		else
		{
		    S->Flag [i] = k ;
		    Lik [*plength] = i ;
		    (*plength)++ ;
		}
	    }
	}
    }
    else
    {
	/* later round: roots are promoted candidates (flagged, pivotal).
	   Push each through the DFS expansion of its pivot column. */
	for (r = 0 ; r < nroots ; r++)
	{
	    i = root_rows [r] ;
	    /* the root is already flagged; expand it like dfs would after
	       marking: temporarily unflag so kls_pipe_dfs takes it */
	    S->Flag [i] = TRILINOS_KLU_EMPTY ;
	    top = kls_pipe_dfs (i, k, plimit, Pinv, S->Llen, S->Lip,
				S->Stack, S->Flag, top, LU, S->colptr,
				Lik, plength, S->Ap_pos, sh, Ap_ver) ;
	}
    }

    if (kls_pipe_phase_prof)
    {
	kls_pipe_t_sym += kls_klu_now () - kls_rd_t0 ;
    }
    /* numeric for this round's topological segment */
    {
	const double kls_ph_t0 = kls_pipe_phase_prof ? kls_klu_now () : 0.0 ;
	kls_pipe_lsolve_numeric (Pinv, LU, S->colptr, S->Stack, S->Lip, top,
				 n, S->Llen, S->X, sh, copybuf) ;
	if (kls_pipe_phase_prof)
	{
	    const double kls_ph_t1 = kls_klu_now () ;
	    kls_pipe_t_num += kls_ph_t1 - kls_ph_t0 ;
	}
    }

    /* extract this round's U segment: values are final (contributions
       to a pivotal row come only from smaller pivot columns, all of
       which are eliminated by this or an earlier round) */
    for (s = top ; s < n ; s++)
    {
	Int j = S->Stack [s] ;
	ubuf_i [*pucount] = Pinv [j] ;
	ubuf_x [*pucount] = S->X [j] ;
	CLEAR (S->X [j]) ;
	(*pucount)++ ;
    }
}


/* Publish one dense-panel column: identical to kls_pipe_finalize_column
 * after lpivot, but the pivot and the L rows/values come from the dense
 * within-panel factorization (KLS_SN_PANEL_FACTOR) instead of a scalar
 * candidate scan.  U entries are already in W->ubuf_i/x (externals plus
 * intra-panel siblings, unsorted).  Returns 1 on abort. */
static int kls_pipe_emit_dense_column
(
    kls_klu_pipe_worker *W,
    Int k,
    const Int *Lrows,       /* llen global row ids (below the pivot) */
    const double *Lvals,    /* llen multipliers, stride PW apart */
    Int lstride,
    Int llen,
    Int pivrow,
    Entry pivot,
    Int ucount
)
{
    kls_klu_pipe_shared *sh = W->sh ;
    KLS_KLU_KERNEL_STATE *S = &W->S ;
    const Int n = S->n ;
    Int diagrow, p, done_limit ;

    S->Llen [k] = llen ;
    diagrow = S->P [k] ;
    {
	size_t lpart = KLS_UNITS32 (llen) + (size_t) UNITS (Entry, llen) ;
	size_t used = lpart + KLS_UNITS32 (ucount) +
		      (size_t) UNITS (Entry, ucount) ;
	Unit *col = kls_klu_chunk_alloc (S, used) ;
	int32_t *Ui ;
	Entry *Ux ;
	if (col == NULL)
	{
	    S->Common->status = TRILINOS_KLU_OUT_OF_MEMORY ;
	    atomic_store_explicit (&sh->abort_flag, 1,
				   memory_order_release) ;
	    return (1) ;
	}
	{
	    int32_t *dLi = (int32_t *) col ;
	    Entry *dLx = (Entry *) (col + KLS_UNITS32 (llen)) ;
	    for (p = 0 ; p < llen ; p++)
	    {
		dLi [p] = (int32_t) Lrows [p] ;
		dLx [p] = Lvals [(size_t) p * lstride] ;
	    }
	}
	S->Uip [k] = (Int) lpart ;
	S->Ulen [k] = ucount ;
	Ui = (int32_t *) (col + lpart) ;
	Ux = (Entry *) (col + lpart + KLS_UNITS32 (ucount)) ;
	{
	    Int a, b ;
	    for (a = 1 ; a < ucount ; a++)
	    {
		Int vi = W->ubuf_i [a] ;
		Entry vx = W->ubuf_x [a] ;
		for (b = a ; b > 0 && W->ubuf_i [b-1] > vi ; b--)
		{
		    W->ubuf_i [b] = W->ubuf_i [b-1] ;
		    W->ubuf_x [b] = W->ubuf_x [b-1] ;
		}
		W->ubuf_i [b] = vi ;
		W->ubuf_x [b] = vx ;
	    }
	}
	for (p = 0 ; p < ucount ; p++)
	{
	    Ui [p] = (int32_t) W->ubuf_i [p] ;
	    Ux [p] = W->ubuf_x [p] ;
	}
	S->colptr [k] = col ;
    }
    S->Udiag [k] = pivot ;
    if (pivrow != diagrow)
    {
	S->Common->noffdiag++ ;
	if (S->Pinv [diagrow] < 0)
	{
	    Int kbar = FLIP (S->Pinv [pivrow]) ;
	    if (kbar < 0 || kbar >= n)
	    {
		fprintf (stderr, "KLS pipe: BAD dense kbar=%ld k=%ld"
			 " pivrow=%ld tid=%d\n", (long) kbar, (long) k,
			 (long) pivrow, W->tid) ;
		atomic_store_explicit (&sh->abort_flag, 1,
				       memory_order_release) ;
		return (1) ;
	    }
	    S->P [kbar] = diagrow ;
	    S->Pinv [diagrow] = FLIP (kbar) ;
	}
    }
    S->P [k] = pivrow ;
    S->Pinv [pivrow] = k ;
    if (sh->lpend != NULL)
    {
	Int up, jcol ;
	for (up = 0 ; up < ucount ; up++)
	{
	    jcol = W->ubuf_i [up] ;
	    if (jcol >= k - 0)
	    {
		continue ;   /* intra-panel U columns: prune later cols only */
	    }
	    if (atomic_load_explicit (&sh->lpend [jcol],
				      memory_order_relaxed) !=
		TRILINOS_KLU_EMPTY)
	    {
		continue ;
	    }
	    {
		Unit *xp = (Unit *) __atomic_load_n (
		    (Unit **) &S->colptr [jcol], __ATOMIC_ACQUIRE) ;
		Int jlen = S->Llen [jcol] ;
		const int32_t *Lij = (const int32_t *) xp ;
		Int p2, found = 0 ;
		for (p2 = 0 ; p2 < jlen ; p2++)
		{
		    if ((Int) Lij [p2] == pivrow)
		    {
			found = 1 ;
			break ;
		    }
		}
		if (found)
		{
		    Int julen = S->Ulen [jcol] ;
		    size_t lpart = KLS_UNITS32 (jlen) +
			(size_t) UNITS (Entry, jlen) ;
		    size_t used = lpart + KLS_UNITS32 (julen) +
			(size_t) UNITS (Entry, julen) ;
		    Unit *newcol = kls_klu_chunk_alloc (S, used) ;
		    if (newcol != NULL)
		    {
			int32_t *nLi = (int32_t *) newcol ;
			Entry *nLx =
			    (Entry *) (newcol + KLS_UNITS32 (jlen)) ;
			Int phead = 0, ptail = jlen ;
			memcpy (newcol, xp, used * sizeof (Unit)) ;
			while (phead < ptail)
			{
			    int32_t row = nLi [phead] ;
			    if (S->Pinv [row] >= 0)
			    {
				phead++ ;
			    }
			    else
			    {
				Entry xtmp ;
				ptail-- ;
				nLi [phead] = nLi [ptail] ;
				nLi [ptail] = row ;
				xtmp = nLx [phead] ;
				nLx [phead] = nLx [ptail] ;
				nLx [ptail] = xtmp ;
			    }
			}
			atomic_fetch_add_explicit (&sh->colver [jcol], 1,
						   memory_order_acq_rel) ;
			__atomic_store_n ((Unit **) &S->colptr [jcol],
					  newcol, __ATOMIC_RELEASE) ;
			atomic_store_explicit (&sh->lpend [jcol], ptail,
					       memory_order_release) ;
			atomic_fetch_add_explicit (&sh->colver [jcol], 1,
						   memory_order_acq_rel) ;
		    }
		}
	    }
	}
    }
    {
	double kfl = 0.0 ;
	Int up2 ;
	for (up2 = 0 ; up2 < ucount ; up2++)
	{
	    if (W->ubuf_i [up2] < k)
	    {
		kfl += (double) S->Llen [W->ubuf_i [up2]] ;
	    }
	}
	W->kflops += 2.0 * kfl + (double) S->Llen [k] ;
    }
    S->lnz += S->Llen [k] + 1 ;
    S->unz += S->Ulen [k] + 1 ;
    S->cols_done = k + 1 ;
    done_limit = k + 1 ;
    atomic_store_explicit (&sh->prefix, done_limit, memory_order_release) ;
    return (0) ;
}

/* lsolve_numeric with seqlock column copies: x -= Lx*xj is not
 * retryable, so each applied column is copied to scratch under its
 * version check and applied from the copy. */
static void kls_pipe_lsolve_numeric
(
    Int Pinv [ ],
    Unit *LU,
    Unit *const *Colptr,
    Int Stack [ ],
    Int Lip [ ],
    Int top,
    Int n,
    Int Llen [ ],
    Entry X [ ],
    const kls_klu_pipe_shared *sh,
    Unit *copybuf
)
{
    Entry xj ;
    Entry *Lx ;
    const int32_t *Li32 ;
    Int p, s, j, jnew, len ;

    if (sh == NULL || sh->lpend == NULL)
    {
	lsolve_numeric (Pinv, LU, Colptr, Stack, Lip, top, n, Llen, X) ;
	return ;
    }
    for (s = top ; s < n ; s++)
    {
	j = Stack [s] ;
	jnew = Pinv [j] ;
	ASSERT (jnew >= 0) ;
	xj = X [j] ;
	{
	    /* copy-on-prune makes every published column snapshot
	       immutable: one atomic pointer load pins a consistent
	       (index,value) pairing - the pre-prune and post-prune
	       arrangements are permutations of the same multiset, and
	       the scatter below is order-independent (distinct rows).
	       No copy, no seqlock retry. */
	    Unit *xp = (Unit *) __atomic_load_n (
		(Unit **) &Colptr [jnew], __ATOMIC_ACQUIRE) ;
	    len = Llen [jnew] ;
	    Li32 = (const int32_t *) xp ;
	    Lx = (Entry *) (xp + KLS_UNITS32 (len)) ;
	    if (kls_pipe_phase_prof)
	    {
		kls_pipe_madds += len ;
		kls_construct_calls++ ;      /* source applies (avg len) */
	    }
	}
	for (p = 0 ; p < len ; p++)
	{
	    MULT_SUB (X [Li32 [p]], Lx [p], xj) ;
	}
    }
}


/* ========================================================================== */
/* === supernodal-pipe lockstep panel body (step b) ========================= */
/* ========================================================================== */


/* Final phase of one pipeline column: pivot, prune, assemble, publish.
 * Candidates are in the scratch head (l_length of them), the U segment
 * in W->ubuf_i/x (ucount), values in S->X.  Shared by the per-column
 * body and the panel cascade.  Returns 1 on abort (flag set). */
static int kls_pipe_finalize_column
(
    kls_klu_pipe_worker *W,
    Int k,
    Int l_length,
    Int ucount
)
{
    kls_klu_pipe_shared *sh = W->sh ;
    KLS_KLU_KERNEL_STATE *S = &W->S ;
    const Int n = S->n ;
    Unit *LU = S->scratch ;
    Int *Lik = (Int *) LU ;
    Entry pivot ;
    double abs_pivot ;
    Int pivrow = TRILINOS_KLU_EMPTY, diagrow, i, p, done_limit ;
    (void) Lik ;
    (void) i ;

    S->Llen [k] = l_length ;   /* the body sets this above its call;
                                  the cascade relies on it here */
	/* prefix == k: the final state is exactly the serial algorithm's */
	diagrow = S->P [k] ;
	if (!lpivot (diagrow, &pivrow, &pivot, &abs_pivot, S->tol, S->X, LU,
		     S->Lip, S->Llen, k, n, S->Pinv, &S->firstrow,
		     S->Common))
	{
	    /* singular: match the serial kernel's bookkeeping, then abort
	       to the serial fallback (halt semantics are the caller's) */
	    S->Common->status = TRILINOS_KLU_SINGULAR ;
	    if (S->Common->numerical_rank == TRILINOS_KLU_EMPTY)
	    {
		S->Common->numerical_rank = k + S->k1 ;
		S->Common->singular_col = S->Q [k + S->k1] ;
	    }
	    if (getenv ("KLS_KLU_PIPE_TRACE") != NULL)
	    {
		fprintf (stderr, "KLS pipe: singular k=%ld k1=%ld llen=%ld "
			 "tid=%d\n", (long) k, (long) S->k1,
			 (long) S->Llen [k], W->tid) ;
	    }
	    atomic_store_explicit (&sh->abort_flag, 1, memory_order_release) ;
	    return (1) ;
	}

	/* assemble and publish the exact-size column: [Li|Lx] is already
	   in the scratch (lpivot gathered it); append [Ui|Ux] */
	{
	    Int llen = S->Llen [k] ;
	    size_t lpart = KLS_UNITS32 (llen) + (size_t) UNITS (Entry, llen) ;
	    size_t used = lpart + KLS_UNITS32 (ucount) +
			  (size_t) UNITS (Entry, ucount) ;
	    Unit *col = kls_klu_chunk_alloc (S, used) ;
	    int32_t *Ui ;
	    Entry *Ux ;
	    if (col == NULL)
	    {
		S->Common->status = TRILINOS_KLU_OUT_OF_MEMORY ;
		atomic_store_explicit (&sh->abort_flag, 1,
				       memory_order_release) ;
		return (1) ;
	    }
	    /* narrow the L indices while copying out of the scratch
	       (block-local rows always fit int32: the pipe requires
	       n < 2^31) */
	    {
		Int *sLi = (Int *) S->scratch ;
		Entry *sLx = (Entry *) (S->scratch + UNITS (Int, llen)) ;
		int32_t *dLi = (int32_t *) col ;
		Entry *dLx = (Entry *) (col + KLS_UNITS32 (llen)) ;
		for (p = 0 ; p < llen ; p++)
		{
		    dLi [p] = (int32_t) sLi [p] ;
		    dLx [p] = sLx [p] ;
		}
	    }
	    S->Uip [k] = (Int) lpart ;
	    S->Ulen [k] = ucount ;
	    Ui = (int32_t *) (col + lpart) ;
	    Ux = (Entry *) (col + lpart + KLS_UNITS32 (ucount)) ;
	    /* sort ascending by pivot index: the multi-round concatenation
	       is a valid topological order but downstream consumers of the
	       stored pattern (refactor map/schedule builders) mis-handle
	       it on multi-round columns; ascending is the consumer-blessed
	       canonical order (the predicted builder's convention) */
	    {
		Int a, b ;
		for (a = 1 ; a < ucount ; a++)
		{
		    Int vi = W->ubuf_i [a] ;
		    Entry vx = W->ubuf_x [a] ;
		    for (b = a ; b > 0 && W->ubuf_i [b-1] > vi ; b--)
		    {
			W->ubuf_i [b] = W->ubuf_i [b-1] ;
			W->ubuf_x [b] = W->ubuf_x [b-1] ;
		    }
		    W->ubuf_i [b] = vi ;
		    W->ubuf_x [b] = vx ;
		}
	    }
	    for (p = 0 ; p < ucount ; p++)
	    {
		Ui [p] = (int32_t) W->ubuf_i [p] ;
		Ux [p] = W->ubuf_x [p] ;
	    }
	    S->colptr [k] = col ;
	}
	S->Udiag [k] = pivot ;
	if (pivrow != diagrow)
	{
	    S->Common->noffdiag++ ;
	    if (S->Pinv [diagrow] < 0)
	    {
		Int kbar = FLIP (S->Pinv [pivrow]) ;
		if (kbar < 0 || kbar >= n)
		{
		    fprintf (stderr, "KLS pipe: BAD kbar=%ld k=%ld pivrow=%ld"
			     " Pinv[pivrow]=%ld diagrow=%ld tid=%d\n",
			     (long) kbar, (long) k, (long) pivrow,
			     (long) S->Pinv [pivrow], (long) diagrow,
			     W->tid) ;
		    atomic_store_explicit (&sh->abort_flag, 1,
					   memory_order_release) ;
		    return (1) ;
		}
		S->P [kbar] = diagrow ;
		S->Pinv [diagrow] = FLIP (kbar) ;
	    }
	}
	S->P [k] = pivrow ;
	S->Pinv [pivrow] = k ;
	if (sh->lpend != NULL)
	{
	    /* symmetric pruning, serialized by in-order finals: for each
	       U column j of k not yet pruned, if pivrow appears in j's
	       published L rows, partition them pivotal-first under j's
	       seqlock and set the shared scan bound. */
	    Int up, jcol ;
	    for (up = 0 ; up < ucount ; up++)
	    {
		jcol = W->ubuf_i [up] ;
		if (atomic_load_explicit (&sh->lpend [jcol],
					  memory_order_relaxed) !=
		    TRILINOS_KLU_EMPTY)
		{
		    continue ;
		}
		{
		    Unit *xp = (Unit *) __atomic_load_n (
			(Unit **) &S->colptr [jcol], __ATOMIC_ACQUIRE) ;
		    Int jlen = S->Llen [jcol] ;
		    const int32_t *Lij = (const int32_t *) xp ;
		    Int p2, found = 0 ;
		    for (p2 = 0 ; p2 < jlen ; p2++)
		    {
			if ((Int) Lij [p2] == pivrow)
			{
			    found = 1 ;
			    break ;
			}
		    }
		    if (found)
		    {
			/* copy-on-prune: partition a fresh copy of the
			   whole packed column and swap the pointer.
			   Readers holding the old pointer keep a
			   consistent (index,value) pairing forever (the
			   chunk arenas outlive the factorization), so
			   the numeric consumers need no copy and no
			   seqlock retry.  The swap loop below replays
			   the exact in-place partition order so the
			   published arrangement stays bit-identical to
			   the previous algorithm. */
			Int julen = S->Ulen [jcol] ;
			size_t lpart = KLS_UNITS32 (jlen) +
			    (size_t) UNITS (Entry, jlen) ;
			size_t used = lpart + KLS_UNITS32 (julen) +
			    (size_t) UNITS (Entry, julen) ;
			Unit *newcol = kls_klu_chunk_alloc (S, used) ;
			if (newcol != NULL)
			{
			    int32_t *nLi = (int32_t *) newcol ;
			    Entry *nLx =
				(Entry *) (newcol + KLS_UNITS32 (jlen)) ;
			    Int phead = 0, ptail = jlen ;
			    memcpy (newcol, xp, used * sizeof (Unit)) ;
			    while (phead < ptail)
			    {
				int32_t row = nLi [phead] ;
				if (S->Pinv [row] >= 0)
				{
				    phead++ ;
				}
				else
				{
				    Entry xtmp ;
				    ptail-- ;
				    nLi [phead] = nLi [ptail] ;
				    nLi [ptail] = row ;
				    xtmp = nLx [phead] ;
				    nLx [phead] = nLx [ptail] ;
				    nLx [ptail] = xtmp ;
				}
			    }
			    atomic_fetch_add_explicit (&sh->colver [jcol],
						       1,
						       memory_order_acq_rel) ;
			    __atomic_store_n ((Unit **) &S->colptr [jcol],
					      newcol, __ATOMIC_RELEASE) ;
			    atomic_store_explicit (&sh->lpend [jcol], ptail,
						   memory_order_release) ;
			    atomic_fetch_add_explicit (&sh->colver [jcol],
						       1,
						       memory_order_acq_rel) ;
			}
			/* allocation failure: skip the prune (it is an
			   optimization; unpruned scans stay correct) */
		    }
		}
	    }
	}
	{
	    /* assembly-time flop sum: in-order publication makes
	       S->Llen of every U column final here (same invariant as
	       the serial kernel's accumulation) */
	    double kfl = 0.0 ;
	    Int up2 ;
	    for (up2 = 0 ; up2 < ucount ; up2++)
	    {
		kfl += (double) S->Llen [W->ubuf_i [up2]] ;
	    }
	    W->kflops += 2.0 * kfl + (double) S->Llen [k] ;
	}
	S->lnz += S->Llen [k] + 1 ;
	S->unz += S->Ulen [k] + 1 ;
	S->cols_done = k + 1 ;

	/* publish in order */
	done_limit = k + 1 ;
	atomic_store_explicit (&sh->prefix, done_limit,
			       memory_order_release) ;
    return (0) ;
}

/* Processes the panel's W columns in lockstep: external sources are
 * applied through ONE masked union pass per round (each source column
 * streamed once for all panel columns that discovered it - the
 * membership bitmask in W->pPFlag guards against partially-updated
 * future dependencies), then the columns run their normal final
 * phases in order as the prefix advances through the panel.  Returns
 * 0 on success, 1 on abort (abort_flag already set). */
static int kls_pipe_panel_lockstep
(
    kls_klu_pipe_worker *W,
    Int k0,
    Int kend,
    Int *promoted,
    Int l_out [ ],
    Int u_out [ ]
)
{
    kls_klu_pipe_shared *sh = W->sh ;
    KLS_KLU_KERNEL_STATE *S = &W->S ;
    const Int n = S->n ;
    const Int PW = kend - k0 ;
    const int use_buf = W->pB != NULL ;
    Int nbrows = 0 ;
    Int w, k, j, p ;
#define KLS_PANEL_ROWPOS(r, out_pos)                                   \
    do                                                                 \
    {                                                                  \
	if (W->pRowGen [r] != W->pGen)                                 \
	{                                                              \
	    Int zw ;                                                   \
	    W->pRowGen [r] = W->pGen ;                                 \
	    if ((size_t) nbrows >= W->pBcap)                           \
	    {                                                          \
		size_t ncap = W->pBcap * 2 ;                           \
		Entry *nb = (Entry *) realloc (W->pB,                  \
		    ncap * (size_t) W->pW * sizeof (Entry)) ;          \
		if (nb == NULL)                                        \
		{                                                      \
		    atomic_store_explicit (&sh->abort_flag, 1,         \
					   memory_order_release) ;     \
		    return (1) ;                                       \
		}                                                      \
		W->pB = nb ;                                           \
		W->pBcap = ncap ;                                      \
	    }                                                          \
	    W->pRowList [nbrows] = (r) ;                               \
	    W->pRowPos [r] = nbrows++ ;                                \
	    for (zw = 0 ; zw < PW ; zw++)                              \
	    {                                                          \
		W->pB [(size_t) W->pRowPos [r] * PW + zw] = 0.0 ;      \
	    }                                                          \
	}                                                              \
	(out_pos) = W->pRowPos [r] ;                                   \
    } while (0)
    Int l_len [32] ;
    Int u_cnt [32] ;
    Int seg_done [32] ;
    Int plimit, newlimit ;
    Int *ulist = W->pUlist ;
    Int ulen = 0 ;

    /* slot 0 aliases the worker's own arrays */
    W->pX [0] = S->X ;
    W->pLik [0] = NULL ;      /* column 0 pattern lives in scratch head */

    /* phase A: construct every panel column into its own accumulator */
    for (w = 0 ; w < PW ; w++)
    {
	l_len [w] = 0 ;
	u_cnt [w] = 0 ;
	seg_done [w] = 0 ;
	S->Lip [k0 + w] = 0 ;
	construct_column (k0 + w, S->Ap, S->Ai, S->Ax, S->Q, W->pX [w],
			  S->k1, S->PSinv, S->Rs, S->scale,
			  S->Offp, S->Offi, S->Offx) ;
	if (use_buf)
	{
	    /* migrate the in-block seeds into the panel-major buffer */
	    Int kglobal = k0 + w + S->k1 ;
	    Int oldcol = S->Q [kglobal] ;
	    Int pend = S->Ap [oldcol+1] ;
	    for (p = S->Ap [oldcol] ; p < pend ; p++)
	    {
		Int i = S->PSinv [S->Ai [p]] - S->k1 ;
		Int pos ;
		if (i < 0) continue ;
		KLS_PANEL_ROWPOS (i, pos) ;
		W->pB [(size_t) pos * PW + w] = W->pX [w][i] ;
		W->pX [w][i] = 0.0 ;
	    }
	}
    }

    /* pattern lists: column 0 uses the scratch head (the final block
       below expects Lik there); columns 1.. use the panel arrays */
    {
	Int *Lik0 = (Int *) S->scratch ;
	W->pLik [0] = Lik0 ;
    }

    /* unified rounds + cascade: windows advance with the prefix up to
       the next unfinalized column; when the window reaches it, that
       column finalizes (publish advances the boundary) and one more
       cascade round lets every remaining sibling absorb the fresh
       column and the external subtree it unlocks - all through the
       buffered union apply. */
    {
	Int next_final = 0 ;
	plimit = atomic_load_explicit (&sh->prefix,
				       memory_order_acquire) ;
	if (plimit > k0)
	{
	    plimit = k0 ;
	}
	for ( ; ; )
	{
		/* per-column DFS collections for this window */
		ulen = 0 ;
		for (w = next_final ; w < PW ; w++)
		{
		    Int top = n ;
		    Int kw = k0 + w ;
		    if (seg_done [w] == 0)
		    {
			/* round-1 roots: A(:,kw) block rows */
			Int kglobal = kw + S->k1 ;
			Int oldcol = S->Q [kglobal] ;
			Int pend = S->Ap [oldcol+1] ;
			for (p = S->Ap [oldcol] ; p < pend ; p++)
			{
			    Int i = S->PSinv [S->Ai [p]] - S->k1 ;
			    if (i < 0) continue ;
			    if (W->pFlagW [w][i] != kw)
			    {
				if (S->Pinv [i] >= 0 && S->Pinv [i] < plimit)
				{
				    top = kls_pipe_dfs (i, kw, plimit, S->Pinv,
							S->Llen, S->Lip, S->Stack,
							W->pFlagW [w], top, NULL,
							S->colptr, W->pLik [w],
							&l_len [w], S->Ap_pos, sh,
							W->ap_ver) ;
				}
				else
				{
				    W->pFlagW [w][i] = kw ;
				    W->pLik [w][l_len [w]] = i ;
				    l_len [w]++ ;
				}
			    }
			}
			seg_done [w] = 1 ;
		    }
		    else
		    {
			/* later round: promoted candidates of column w */
			Int npromoted = 0, cw = 0 ;
			for (p = 0 ; p < l_len [w] ; p++)
			{
			    Int i = W->pLik [w][p] ;
			    if (S->Pinv [i] >= 0 && S->Pinv [i] < plimit)
			    {
				promoted [npromoted++] = i ;
			    }
			    else
			    {
				W->pLik [w][cw++] = i ;
			    }
			}
			l_len [w] = cw ;
			for (p = 0 ; p < npromoted ; p++)
			{
			    Int i = promoted [p] ;
			    W->pFlagW [w][i] = TRILINOS_KLU_EMPTY ;
			    top = kls_pipe_dfs (i, kw, plimit, S->Pinv, S->Llen,
						S->Lip, S->Stack, W->pFlagW [w],
						top, NULL, S->colptr, W->pLik [w],
						&l_len [w], S->Ap_pos, sh,
						W->ap_ver) ;
			}
		    }
		    /* merge this column's topological segment into the union
		       (first occurrence keeps a globally topological order) */
		    for (p = top ; p < n ; p++)
		    {
			j = S->Stack [p] ;
			if (W->pPFlag [j] == TRILINOS_KLU_EMPTY)
			{
			    W->pPFlag [j] = (Int) 1 << w ;
			    ulist [ulen++] = j ;
			}
			else
			{
			    W->pPFlag [j] |= (Int) 1 << w ;
			}
		    }
		}
		/* pivot order is a topological order of the L-DAG (an edge
		   l->j implies Pinv[l] < Pinv[j]); the concatenation of
		   per-column segments is not, so sort the round's union */
		for (p = 1 ; p < ulen ; p++)
		{
		    Int jj = ulist [p] ;
		    Int pv = S->Pinv [jj] ;
		    Int q2 = p ;
		    while (q2 > 0 && S->Pinv [ulist [q2-1]] > pv)
		    {
			ulist [q2] = ulist [q2-1] ;
			q2-- ;
		    }
		    ulist [q2] = jj ;
		}
		/* union apply: each source streamed once, executing only the
		   member lanes (compressed) - useful madds only */
		for (p = 0 ; p < ulen ; p++)
		{
		    Int mask ;
		    Int jnew ;
		    Entry xj [32] ;
		    j = ulist [p] ;
		    mask = W->pPFlag [j] ;
		    jnew = S->Pinv [j] ;
		    {
			Unit *xp = (Unit *) __atomic_load_n (
			    (Unit **) &S->colptr [jnew], __ATOMIC_ACQUIRE) ;
			Int len = S->Llen [jnew] ;
			const int32_t *Li32 = (const int32_t *) xp ;
			const Entry *Lx = (const Entry *) (xp +
							   KLS_UNITS32 (len)) ;
			Int q ;
			if (use_buf)
			{
			    Int jpos ;
			    Int mw [32] ;
			    Int nm = 0 ;
			    KLS_PANEL_ROWPOS (j, jpos) ;
			    for (w = 0 ; w < PW ; w++)
			    {
				if (mask & ((Int) 1 << w))
				{
				    xj [nm] = W->pB [(size_t) jpos * PW + w] ;
				    mw [nm++] = w ;
				}
			    }
			    for (q = 0 ; q < len ; q++)
			    {
				const Int r = (Int) Li32 [q] ;
				const Entry v = Lx [q] ;
				Entry *brow ;
				Int rpos, m ;
				KLS_PANEL_ROWPOS (r, rpos) ;
				brow = W->pB + (size_t) rpos * PW ;
				for (m = 0 ; m < nm ; m++)
				{
				    brow [mw [m]] -= v * xj [m] ;
				}
			    }
			}
			else
			{
			    Int mw [32] ;
			    Int nm = 0 ;
			    for (w = 0 ; w < PW ; w++)
			    {
				if (mask & ((Int) 1 << w))
				{
				    xj [nm] = W->pX [w][j] ;
				    mw [nm++] = w ;
				}
			    }
			    for (q = 0 ; q < len ; q++)
			    {
				const Int r = (Int) Li32 [q] ;
				const Entry v = Lx [q] ;
				Int m ;
				for (m = 0 ; m < nm ; m++)
				{
				    W->pX [mw [m]][r] -= v * xj [m] ;
				}
			    }
			}
		    }
		    if (kls_pipe_phase_prof)
		    {
			Int nb = 0 ;
			for (w = 0 ; w < PW ; w++) nb += (mask >> w) & 1 ;
			kls_pipe_madds += (long) S->Llen [jnew] * nb ;
			kls_pipe_copy_bytes += (long) S->Llen [jnew] * 12 ;
			/* streamed source bytes (12B/entry, once per source) */
			kls_construct_calls++ ;          /* union sources */
			kls_construct_entries += nb ;    /* popcount sum */
		    }
		}
		/* U extraction per member column, then clear */
		for (p = 0 ; p < ulen ; p++)
		{
		    Int mask ;
		    Int jpos = -1 ;
		    j = ulist [p] ;
		    mask = W->pPFlag [j] ;
		    if (use_buf)
		    {
			KLS_PANEL_ROWPOS (j, jpos) ;
		    }
		    for (w = 0 ; w < PW ; w++)
		    {
			if (mask & ((Int) 1 << w))
			{
			    Int *ubi = w == 0 ? W->ubuf_i : W->pUbuf_i [w] ;
			    Entry *ubx = w == 0 ? W->ubuf_x : W->pUbuf_x [w] ;
			    ubi [u_cnt [w]] = S->Pinv [j] ;
			    if (use_buf)
			    {
				ubx [u_cnt [w]] = W->pB [(size_t) jpos * PW + w] ;
				W->pB [(size_t) jpos * PW + w] = 0.0 ;
			    }
			    else
			    {
				ubx [u_cnt [w]] = W->pX [w][j] ;
				W->pX [w][j] = 0.0 ;
			    }
			    u_cnt [w]++ ;
			}
		    }
		    W->pPFlag [j] = TRILINOS_KLU_EMPTY ;
		}

	    if (kls_pipe_dense_panels && use_buf && next_final == 0 &&
		plimit == k0)
	    {
		/* dense within-panel finalize: all external sources are
		   applied and extracted; the remaining work is the
		   intra-panel elimination, done as one dense LU over the
		   panel's non-pivotal row set (threshold pivoting with
		   the same tol and per-column diagonal preference).  On
		   any failure fall through to the scalar cascade - the
		   buffer is copied, not consumed. */
		Int R2 = 0, pos ;
		for (pos = 0 ; pos < nbrows ; pos++)
		{
		    if (S->Pinv [W->pRowList [pos]] < 0)
		    {
			R2++ ;
		    }
		}
		if (R2 >= PW)
		{
		    double *A = (double *) malloc ((size_t) R2 * PW *
						   sizeof (double)) ;
		    Int *rowids = (Int *) malloc ((size_t) R2 *
						  sizeof (Int)) ;
		    Int diagrows [32] ;
		    if (A != NULL && rowids != NULL)
		    {
			Int c = 0, fac ;
			for (pos = 0 ; pos < nbrows ; pos++)
			{
			    Int r = W->pRowList [pos] ;
			    if (S->Pinv [r] < 0)
			    {
				memcpy (A + (size_t) c * PW,
					W->pB + (size_t) pos * PW,
					(size_t) PW * sizeof (double)) ;
				rowids [c++] = r ;
			    }
			}
			for (w = 0 ; w < PW ; w++)
			{
			    diagrows [w] = S->P [k0 + w] ;
			}
			/* threshold pivoting with the kernel tol; panel
			   tol=1.0 measured 4x WORSE on mac_econ (off-
			   diagonal pivots cascade into later panels'
			   candidate sets - the "no fill cost" intuition
			   only holds within one panel) */
			fac = KLS_SN_PANEL_FACTOR (A, R2, PW, rowids,
						   diagrows,
						   getenv ("KLS_KLU_PIPE_DENSE_PP")
						   != NULL ? 1.0 : S->tol) ;
			if (fac == PW)
			{
			    /* structural per-column emission: pattern_w =
			       own candidates + patterns of earlier siblings
			       whose pivot row lies in pattern_w (left-
			       looking symbolic).  Every numerically nonzero
			       dense multiplier is inside the structural
			       pattern (numeric values only arise along
			       structural paths), so the dropped
			       outside-pattern lanes are exact zeros and the
			       emitted factor equals the union emission with
			       the padding removed. */
			    int aborted = 0 ;
			    char *inpat = (char *) malloc ((size_t) R2) ;
			    Int *pat = (Int *) malloc ((size_t) PW * R2 *
						       sizeof (Int)) ;
			    Int *plen = (Int *) malloc ((size_t) PW *
							sizeof (Int)) ;
			    Int *tmpr = (Int *) malloc ((size_t) R2 *
							sizeof (Int)) ;
			    Entry *tmpx = (Entry *) malloc ((size_t) R2 *
							    sizeof (Entry)) ;
			    /* pRowPos is generation-tagged and this panel
			       is done with the buffer map: safe scratch
			       (pFlagW would collide with future DFS marks) */
			    Int *rpos = W->pRowPos ;
			    if (inpat == NULL || pat == NULL ||
				plen == NULL || tmpr == NULL || tmpx == NULL)
			    {
				free (inpat) ; free (pat) ; free (plen) ;
				free (tmpr) ; free (tmpx) ;
				free (A) ; free (rowids) ;
				goto kls_dense_bail ;
			    }
			    {
				Int c ;
				for (c = 0 ; c < R2 ; c++)
				{
				    rpos [rowids [c]] = c ;
				}
			    }
			    for (w = 0 ; w < PW && !aborted ; w++)
			    {
				Int uc, jj, len2 = 0, q2 ;
				Int *pw = pat + (size_t) w * R2 ;
				Int *lik = w == 0 ? (Int *) S->scratch
						  : W->pLik [w] ;
				memset (inpat, 0, (size_t) R2) ;
				for (q2 = 0 ; q2 < l_len [w] ; q2++)
				{
				    Int c = rpos [lik [q2]] ;
				    if (!inpat [c])
				    {
					inpat [c] = 1 ;
					pw [len2++] = c ;
				    }
				}
				for (jj = 0 ; jj < w ; jj++)
				{
				    if (inpat [jj])
				    {
					const Int *pj = pat +
					    (size_t) jj * R2 ;
					for (q2 = 0 ; q2 < plen [jj] ; q2++)
					{
					    Int c = pj [q2] ;
					    if (!inpat [c])
					    {
						inpat [c] = 1 ;
						pw [len2++] = c ;
					    }
					}
				    }
				}
				plen [w] = len2 ;
				if (w > 0)
				{
				    memcpy (W->ubuf_i, W->pUbuf_i [w],
					    (size_t) u_cnt [w] *
					    sizeof (Int)) ;
				    memcpy (W->ubuf_x, W->pUbuf_x [w],
					    (size_t) u_cnt [w] *
					    sizeof (Entry)) ;
				}
				uc = u_cnt [w] ;
				for (jj = 0 ; jj < w ; jj++)
				{
				    if (inpat [jj])
				    {
					W->ubuf_i [uc] = k0 + jj ;
					W->ubuf_x [uc] =
					    A [(size_t) jj * PW + w] ;
					uc++ ;
				    }
				}
				{
				    Int llen2 = 0 ;
				    for (q2 = 0 ; q2 < len2 ; q2++)
				    {
					Int c = pw [q2] ;
					if (c > w)
					{
					    tmpr [llen2] = rowids [c] ;
					    tmpx [llen2] =
						A [(size_t) c * PW + w] ;
					    llen2++ ;
					}
				    }
				    if (kls_pipe_emit_dense_column (W,
					    k0 + w, tmpr, tmpx, 1, llen2,
					    rowids [w],
					    A [(size_t) w * PW + w], uc))
				    {
					aborted = 1 ;
				    }
				    l_len [w] = llen2 ;
				    u_cnt [w] = uc ;
				}
			    }
			    free (inpat) ; free (pat) ; free (plen) ;
			    free (tmpr) ; free (tmpx) ;
			    free (A) ;
			    free (rowids) ;
			    if (aborted)
			    {
				return (1) ;
			    }
			    S->Common->kls_dense_panels = 1 ;
			    next_final = PW ;
			    break ;
			}
		    }
		    free (A) ;
		    free (rowids) ;
		}
	    }
	    kls_dense_bail: ;
	    if (plimit == k0 + next_final)
	    {
		/* this column's window is complete: finalize it */
		const Int kf = k0 + next_final ;
		const Int wf = next_final ;
		Int *lik_src = wf == 0 ? (Int *) S->scratch
				       : W->pLik [wf] ;
		if (use_buf)
		{
		    for (p = 0 ; p < l_len [wf] ; p++)
		    {
			Int i = lik_src [p] ;
			if (W->pRowGen [i] == W->pGen)
			{
			    S->X [i] =
				W->pB [(size_t) W->pRowPos [i] * PW + wf] ;
			}
		    }
		}
		if (wf > 0)
		{
		    memcpy ((Int *) S->scratch, W->pLik [wf],
			    (size_t) l_len [wf] * sizeof (Int)) ;
		    memcpy (W->ubuf_i, W->pUbuf_i [wf],
			    (size_t) u_cnt [wf] * sizeof (Int)) ;
		    memcpy (W->ubuf_x, W->pUbuf_x [wf],
			    (size_t) u_cnt [wf] * sizeof (Entry)) ;
		}
		if (kls_pipe_finalize_column (W, kf, l_len [wf],
					      u_cnt [wf]))
		{
		    return (1) ;
		}
		next_final++ ;
		if (next_final == PW)
		{
		    break ;
		}
		plimit = k0 + next_final ;
		continue ;
	    }
	    /* wait for the prefix to advance the external window */
	    {
		Int newlimit2 ;
		do
		{
		    if (atomic_load_explicit (&sh->abort_flag,
					      memory_order_acquire))
		    {
			return (1) ;
		    }
		    KLS_KLU_CPU_RELAX () ;
		    newlimit2 = atomic_load_explicit (&sh->prefix,
						      memory_order_acquire) ;
		} while (newlimit2 <= plimit) ;
		if (newlimit2 > k0 + next_final)
		{
		    newlimit2 = k0 + next_final ;
		}
		plimit = newlimit2 ;
		if (kls_pipe_phase_prof) { W->n_rounds++ ; }
	    }
	}
    }
    if (use_buf)
    {
	W->pGen++ ;
    }
    for (w = 0 ; w < PW ; w++)
    {
	l_out [w] = l_len [w] ;
	u_out [w] = u_cnt [w] ;
    }
    return (0) ;
}
#undef KLS_PANEL_ROWPOS


static void *kls_klu_pipe_worker_main (void *arg)
{
    kls_klu_pipe_worker *W = (kls_klu_pipe_worker *) arg ;
    kls_klu_pipe_shared *sh = W->sh ;
    KLS_KLU_KERNEL_STATE *S = &W->S ;
    const Int n = S->n ;
    const int prof = getenv ("KLS_KLU_PIPE_PROF") != NULL ;
    double t0 = 0 ;
    Int panel_l [32], panel_u [32] ;
    Int panel_k0 = 0 ;
    int panel_inject = 0 ;
    Int *promoted = (Int *) malloc ((size_t) n * sizeof (Int)) ;
    if (promoted == NULL)
    {
	atomic_store_explicit (&sh->abort_flag, 1, memory_order_release) ;
	return (NULL) ;
    }

    for ( ; ; )
    {
	Int k, kp_end ;
	if (sh->panel_start != NULL)
	{
	    /* claim a whole fundamental supernode; the per-column body
	       below is unchanged (bit-identical results - panels only
	       group the claim order) */
	    Int pnl = atomic_fetch_add_explicit (&sh->next_panel, 1,
						 memory_order_relaxed) ;
	    if (pnl >= sh->npanels ||
		atomic_load_explicit (&sh->abort_flag,
				      memory_order_acquire))
	    {
		break ;
	    }
	    k = sh->panel_start [pnl] ;
	    kp_end = sh->panel_start [pnl + 1] ;
	    panel_inject = 0 ;
	    panel_k0 = k ;
	    if (kp_end - k > 1 && kp_end - k <= W->pW &&
		sh->lpend != NULL && W->pPFlag != NULL)
	    {
		if (kls_pipe_panel_lockstep (W, k, kp_end, promoted,
					     panel_l, panel_u))
		{
		    free (promoted) ;
		    return (NULL) ;
		}
		/* the cascade finalized and published every panel
		   column: claim the next panel */
		if (prof) { W->n_cols += kp_end - k ; }
		continue ;
	    }
	}
	else
	{
	    k = atomic_fetch_add_explicit (&sh->next_col, 1,
					   memory_order_relaxed) ;
	    if (k >= n ||
		atomic_load_explicit (&sh->abort_flag,
				      memory_order_acquire))
	    {
		break ;
	    }
	    kp_end = k + 1 ;
	}
	for ( ; k < kp_end ; k++)
	{
	if (atomic_load_explicit (&sh->abort_flag, memory_order_acquire))
	{
	    free (promoted) ;
	    return (NULL) ;
	}

	Unit *LU = S->scratch ;
	Int *Lik ;
	Int l_length = 0, ucount = 0, plimit, done_limit ;
	if (prof) { t0 = kls_klu_now () ; W->n_cols++ ; }
	Entry pivot ;
	double abs_pivot ;
	Int pivrow = TRILINOS_KLU_EMPTY, diagrow, i, p ;

	S->Lip [k] = 0 ;
	Lik = (Int *) LU ;

	if (panel_inject)
	{
	    /* lockstep already constructed, ran the external rounds and
	       extracted their U segments for this column: inject that
	       state and let the normal promotion rounds + final phase
	       finish the column (intra-panel sources arrive through the
	       prefix like any published column) */
	    const Int w = k - panel_k0 ;
	    l_length = panel_l [w] ;
	    ucount = panel_u [w] ;
	    /* transfer this column's lockstep marks into the shared
	       Flag: candidates from the pattern list, applied sources
	       via their pivot rows (P[u]) - the continuation's rounds
	       must not re-visit (candidates would duplicate, sources
	       would double-apply) */
	    {
		Int q4 ;
		for (q4 = 0 ; q4 < l_length ; q4++)
		{
		    S->Flag [w > 0 ? W->pLik [w][q4]
				   : ((Int *) S->scratch) [q4]] = k ;
		}
		for (q4 = 0 ; q4 < ucount ; q4++)
		{
		    Int uv = w > 0 ? W->pUbuf_i [w][q4] : W->ubuf_i [q4] ;
		    S->Flag [S->P [uv]] = k ;
		}
	    }
	    if (w > 0)
	    {
		memcpy (Lik, W->pLik [w], (size_t) l_length * sizeof (Int)) ;
		memcpy (W->ubuf_i, W->pUbuf_i [w],
			(size_t) ucount * sizeof (Int)) ;
		memcpy (W->ubuf_x, W->pUbuf_x [w],
			(size_t) ucount * sizeof (Entry)) ;
		/* the column's accumulator: swap the panel X in */
		{
		    Entry *tmp = S->X ;
		    S->X = W->pX [w] ;
		    W->pX [w] = tmp ;
		}
	    }
	    /* the lockstep covered exactly the external window [0, k0):
	       resume from there so the promotion rounds process the
	       intra-panel window (k0, k) - the earlier panel columns'
	       pivots and updates arrive like any published column */
	    plimit = panel_k0 ;
	}
	else
	{
	construct_column (k, S->Ap, S->Ai, S->Ax, S->Q, S->X,
			  S->k1, S->PSinv, S->Rs, S->scale,
			  S->Offp, S->Offi, S->Offx) ;

	/* round 1 against the prefix at claim time */
	plimit = atomic_load_explicit (&sh->prefix, memory_order_acquire) ;
	if (plimit > k)
	{
	    plimit = k ;
	}
	kls_pipe_round (S, k, plimit, NULL, 0, LU, Lik, &l_length,
			W->ubuf_i, W->ubuf_x, &ucount, sh, W->ap_ver,
			W->copybuf) ;
	}

	/* later rounds as the prefix advances toward k */
	while (plimit < k)
	{
	    Int newlimit ;
	    if (prof)
	    {
		double ts = kls_klu_now () ;
		W->t_work += ts - t0 ;
		t0 = ts ;
	    }
	    do
	    {
		if (atomic_load_explicit (&sh->abort_flag,
					  memory_order_acquire))
		{
		    free (promoted) ;
		    return (NULL) ;
		}
		KLS_KLU_CPU_RELAX () ;
		newlimit = atomic_load_explicit (&sh->prefix,
						 memory_order_acquire) ;
	    } while (newlimit <= plimit) ;
	    if (prof)
	    {
		double ts = kls_klu_now () ;
		W->t_spin += ts - t0 ;
		t0 = ts ;
		W->n_rounds++ ;
	    }
	    if (newlimit > k)
	    {
		newlimit = k ;
	    }
	    /* collect candidates promoted into [plimit, newlimit) and
	       compact the candidate list */
	    {
		Int npromoted = 0, w = 0 ;
		for (p = 0 ; p < l_length ; p++)
		{
		    i = Lik [p] ;
		    if (S->Pinv [i] >= 0 && S->Pinv [i] < newlimit)
		    {
			promoted [npromoted++] = i ;
		    }
		    else
		    {
			Lik [w++] = i ;
		    }
		}
		l_length = w ;
		if (npromoted > 0)
		{
		    kls_pipe_round (S, k, newlimit, promoted, npromoted, LU,
				    Lik, &l_length, W->ubuf_i, W->ubuf_x,
				    &ucount, sh, W->ap_ver, W->copybuf) ;
		}
	    }
	    plimit = newlimit ;
	}

	if (prof)
	{
	    double ts = kls_klu_now () ;
	    W->t_work += ts - t0 ;
	    t0 = ts ;
	}
	S->Llen [k] = l_length ;

	if (kls_pipe_finalize_column (W, k, l_length, ucount))
	{
	    break ;
	}
	if (prof)
	{
	    W->t_final += kls_klu_now () - t0 ;
	}
	if (panel_inject && k - panel_k0 > 0)
	{
	    /* restore the ownership swap: the worker's own X returns to
	       S->X so the cleanup frees exactly the buffers it
	       allocated (the caller owns the template X) */
	    Entry *tmp = S->X ;
	    S->X = W->pX [k - panel_k0] ;
	    W->pX [k - panel_k0] = tmp ;
	}
	}   /* panel column loop */
    }
    if (prof)
    {
	fprintf (stderr, "KLS pipe prof tid=%d cols=%ld rounds=%ld "
		 "work=%.2fs spin=%.2fs final=%.2fs construct=%.2fs/%ld/%ld"
		 " sym=%.2fs num=%.2fs copyMB=%ld madds=%ld\n",
		 W->tid, W->n_cols, W->n_rounds,
		 W->t_work, W->t_spin, W->t_final,
		 kls_construct_secs, kls_construct_calls,
		 kls_construct_entries,
		 kls_pipe_t_sym, kls_pipe_t_num,
		 kls_pipe_copy_bytes >> 20, kls_pipe_madds) ;
    }
    free (promoted) ;
    return (NULL) ;
}

/* ========================================================================== */
/* === supernodal panel factor (phase 1b) ================================== */
/* ========================================================================== */

/* Dense panel LU with threshold row pivoting over a gathered row set:
 * A is row-major R x W (R >= W), rows carry global identities in
 * rowids[].  On return the first W entries of rowids are the chosen
 * pivot rows (panel-internal partial pivoting with the klu-style
 * diagonal preference: prefer diagrow when |x_diag| >= tol * |x_max|),
 * A holds L\U packed in place (unit lower within the top W x W, the
 * multipliers below), and the routine reports the number of columns
 * factored (== W unless a column of the remaining row set is exactly
 * zero).  diagrows[j] names the preferred pivot row of column j
 * (TRILINOS_KLU_EMPTY for none). */
Int KLS_SN_PANEL_FACTOR
(
    double *A,          /* row-major R x W, modified in place */
    Int R,
    Int W,
    Int rowids [ ],     /* size R, permuted alongside */
    const Int diagrows [ ],  /* size W or NULL */
    double tol
)
{
    Int j, r, p, best ;
    for (j = 0 ; j < W ; j++)
    {
	double amax = 0.0, adiag = -1.0 ;
	Int rdiag = -1 ;
	best = -1 ;
	for (r = j ; r < R ; r++)
	{
	    double v = A [r * W + j] ;
	    v = v < 0 ? -v : v ;
	    if (v > amax)
	    {
		amax = v ;
		best = r ;
	    }
	    if (diagrows != NULL && diagrows [j] >= 0 &&
		rowids [r] == diagrows [j])
	    {
		rdiag = r ;
		adiag = v ;
	    }
	}
	if (best < 0 || amax == 0.0)
	{
	    return (j) ;    /* structurally/numerically dead column */
	}
	if (rdiag >= 0 && adiag >= tol * amax)
	{
	    best = rdiag ;
	}
	if (best != j)
	{
	    Int tid = rowids [best] ;
	    rowids [best] = rowids [j] ;
	    rowids [j] = tid ;
	    for (p = 0 ; p < W ; p++)
	    {
		double tv = A [best * W + p] ;
		A [best * W + p] = A [j * W + p] ;
		A [j * W + p] = tv ;
	    }
	}
	{
	    double piv = A [j * W + j] ;
	    for (r = j + 1 ; r < R ; r++)
	    {
		double m = A [r * W + j] / piv ;
		A [r * W + j] = m ;
		if (m != 0.0)
		{
		    for (p = j + 1 ; p < W ; p++)
		    {
			A [r * W + p] -= m * A [j * W + p] ;
		    }
		}
	    }
	}
    }
    return (W) ;
}

/* Pipelined factorization entry: same contract as TRILINOS_KLU_kernel.
 * On any worker abort (singular, OOM) the caller retries serially. */
size_t KLS_KLU_KERNEL_PIPE
(
    Int n, Int Ap [ ], Int Ai [ ], Entry Ax [ ], Int Q [ ], size_t lusize,
    Int Pinv [ ], Int P [ ], Unit **p_LU, Entry Udiag [ ],
    Int Llen [ ], Int Ulen [ ], Int Lip [ ], Int Uip [ ],
    Int *lnz, Int *unz, Entry X [ ],
    Int Stack [ ], Int Flag [ ], Int Ap_pos [ ], Int Lpend [ ],
    Int k1, Int PSinv [ ], double Rs [ ],
    Int Offp [ ], Int Offi [ ], Entry Offx [ ],
    TRILINOS_KLU_common *Common,
    int nthreads
)
{
    KLS_KLU_KERNEL_STATE S ;
    kls_klu_pipe_shared sh ;
    kls_klu_pipe_worker workers [16] ;
    pthread_t tids [16] ;
    int t, spawn_failed = 0 ;
    Int k, poff ;
    size_t final_size ;
    Int *panel_start = NULL ;
    Int npanels = 0 ;
    const int kls_wall_prof = getenv ("KLS_KLU_PIPE_PROF") != NULL ;
    double kls_wall_t0 = kls_wall_prof ? kls_klu_now () : 0.0 ;
    double kls_wall_t1 = 0.0, kls_wall_t2 = 0.0 ;

    if (nthreads < 1) nthreads = 1 ;
    if (nthreads > 16) nthreads = 16 ;

    /* Offp prefix so construct_column's off-diagonal writes are disjoint
       and idempotent under any completion order */
    poff = Offp [k1] ;
    for (k = 0 ; k < n ; k++)
    {
	Int oldcol = Q [k + k1] ;
	Int p2 ;
	Offp [k + k1] = poff ;
	for (p2 = Ap [oldcol] ; p2 < Ap [oldcol+1] ; p2++)
	{
	    if (PSinv [Ai [p2]] - k1 < 0)
	    {
		poff++ ;
	    }
	}
    }
    Offp [n + k1] = poff ;

    S.n = n ;
    S.Ap = Ap ; S.Ai = Ai ; S.Ax = Ax ; S.Q = Q ;
    S.lusize = lusize ;
    S.Pinv = Pinv ; S.P = P ;
    S.LU = *p_LU ;
    S.Udiag = Udiag ;
    S.Llen = Llen ; S.Ulen = Ulen ; S.Lip = Lip ; S.Uip = Uip ;
    S.X = X ; S.Stack = Stack ; S.Flag = Flag ;
    S.Ap_pos = Ap_pos ; S.Lpend = Lpend ;
    S.k1 = k1 ; S.PSinv = PSinv ; S.Rs = Rs ;
    S.Offp = Offp ; S.Offi = Offi ; S.Offx = Offx ;
    S.Common = Common ;
    S.scale = Common->scale ;
    S.tol = Common->tol ;
    S.memgrow = Common->memgrow ;
    S.no_prune = 1 ;
    S.diag_claim = 0 ;
    S.chunked_prune = 0 ;
    S.idx32 = 1 ;  /* pipe columns publish int32 indices */
    S.pack_keep_row_indices = 0 ;
    S.chunk_head = NULL ; S.chunk_used = 0 ; S.chunk_size = 0 ;
    S.scratch = (Unit *) malloc ((2 * (size_t) n + 4) * sizeof (Unit)) ;
    S.colptr = (Unit **) calloc ((size_t) n, sizeof (Unit *)) ;
    if (S.scratch == NULL || S.colptr == NULL)
    {
	free (S.scratch) ;
	free (S.colptr) ;
	Common->status = TRILINOS_KLU_OUT_OF_MEMORY ;
	return (lusize) ;
    }

    KLS_KLU_KERNEL_INIT (&S) ;

    /* fundamental-supernode panels over the column etree, processed
       lockstep with the panel-major buffer.  The partition width is
       capped at the per-worker workspace width so NO panel falls back
       to serialized per-column claims: fallback panels were the
       measured catastrophe (pre2 42s with 62% fallback vs 13s at full
       coverage; column mode 17.7s).  The win is the buffer: union
       applies land in a compact L2-resident panel-major block instead
       of the 5MB X vector (pre2 workers 12.0 -> 7.2s). */
    if (getenv ("KLS_KLU_PIPE_NOPANELS") == NULL &&
	(size_t) n * sizeof (Entry) > ((size_t) 1 << 20))
    {
	/* the buffer's win is X-footprint compression into L2; when X
	   already fits L2 the lockstep only adds synchronization
	   (onetone1 measured +5%) */
	Int *sp = (Int *) malloc ((size_t) (4 * n + 2) * sizeof (Int)) ;
	if (sp != NULL)
	{
	    Int *par = sp, *anc = sp + n, *prv = sp + 2*n ;
	    Int *pst = sp + 3*n ;   /* worst case n+1 panel offsets */
	    Int kk, width = 1, np = 0 ;
	    Int wcap_env ;
	    {
		/* same budget formula as the worker workspace below */
		size_t col_cap = (size_t) n / 4 + 16 ;
		size_t per_col = (size_t) n * sizeof (Entry) +
		    col_cap * (sizeof (Entry) + 2 * sizeof (Int)) ;
		size_t budget = (size_t) 128 << 20 ;
		{
		    const char *be = getenv ("KLS_KLU_PIPE_BUDGET_MB") ;
		    if (be != NULL && be [0] != '\0')
		    {
			long bm = atol (be) ;
			if (bm >= 64 && bm <= 4096)
			{
			    budget = (size_t) bm << 20 ;
			}
		    }
		}
		wcap_env = (Int) (budget / (per_col > 0 ? per_col : 1)) ;
		if (wcap_env > 32) wcap_env = 32 ;
		if (wcap_env < 2) wcap_env = 2 ;
	    }
	    {
		const char *we = getenv ("KLS_KLU_PIPE_PANELW") ;
		if (we != NULL)
		{
		    wcap_env = atol (we) ;
		    if (wcap_env < 2) wcap_env = 2 ;
		    if (wcap_env > 32) wcap_env = 32 ;
		}
	    }
	    kls_klu_block_coletree (n, Ap, Ai, Q, k1, PSinv, par, anc,
				    prv) ;
	    pst [0] = 0 ;
	    {
	    /* relaxed merge: continue a panel across an etree-chain
	       break when the adjacent columns' A patterns differ by at
	       most KLS_KLU_PIPE_RELAX rows - grid-like matrices carry
	       near-identical adjacent columns whose fundamental chains
	       are shorter than the profitable dense-panel width */
	    Int relax = 0 ;
	    Int *rflag = NULL ;
	    {
		const char *re = getenv ("KLS_KLU_PIPE_RELAX") ;
		if (re != NULL && re [0] != '\0')
		{
		    relax = atol (re) ;
		    if (relax > 0)
		    {
			rflag = (Int *) calloc ((size_t) n, sizeof (Int)) ;
			if (rflag == NULL)
			{
			    relax = 0 ;
			}
		    }
		}
	    }
	    for (kk = 1 ; kk <= n ; kk++)
	    {
		if (kk < n && par [kk-1] == kk && width < wcap_env)
		{
		    width++ ;
		    continue ;
		}
		if (kk < n && relax > 0 && width < wcap_env)
		{
		    /* symmetric difference of block-local A patterns */
		    Int oc1 = Q [kk - 1 + k1], oc2 = Q [kk + k1] ;
		    Int p1, d = 0, n1 = 0 ;
		    for (p1 = Ap [oc1] ; p1 < Ap [oc1+1] ; p1++)
		    {
			Int i2 = PSinv [Ai [p1]] - k1 ;
			if (i2 >= 0) { rflag [i2] = kk ; n1++ ; }
		    }
		    for (p1 = Ap [oc2] ; p1 < Ap [oc2+1] ; p1++)
		    {
			Int i2 = PSinv [Ai [p1]] - k1 ;
			if (i2 >= 0)
			{
			    if (rflag [i2] == kk) { n1-- ; }
			    else { d++ ; }
			}
		    }
		    if (d + n1 <= relax)
		    {
			width++ ;
			continue ;
		    }
		}
		np++ ;
		pst [np] = kk ;
		width = 1 ;
	    }
	    free (rflag) ;
	    }
	    panel_start = (Int *) malloc ((size_t) (np + 1) * sizeof (Int)) ;
	    if (panel_start != NULL)
	    {
		memcpy (panel_start, pst, (size_t) (np + 1) * sizeof (Int)) ;
		npanels = np ;
	    }
	    if (getenv ("KLS_KLU_PIPE_PROF") != NULL)
	    {
		Int s1 = 0, wmax = 0, pp ;
		for (pp = 0 ; pp < np ; pp++)
		{
		    Int wd = pst [pp+1] - pst [pp] ;
		    if (wd == 1) s1++ ;
		    if (wd > wmax) wmax = wd ;
		}
		fprintf (stderr, "KLS pipe panels: n=%ld npanels=%ld"
			 " avgw=%.1f singles=%ld (%.0f%%) maxw=%ld\n",
			 (long) n, (long) np,
			 np > 0 ? (double) n / np : 0.0, (long) s1,
			 np > 0 ? 100.0 * s1 / np : 0.0, (long) wmax) ;
	    }
	    free (sp) ;
	}
    }

    /* shared struct is initialized before worker setup */
    atomic_init (&sh.next_col, 0) ;
    atomic_init (&sh.prefix, 0) ;
    atomic_init (&sh.abort_flag, 0) ;
    atomic_init (&sh.next_panel, 0) ;
    sh.panel_start = panel_start ;
    sh.npanels = npanels ;
    sh.n = n ;
    sh.nthreads = nthreads ;
    sh.lpend = NULL ;
    sh.colver = NULL ;
    if (getenv ("KLS_KLU_PIPE_NOPRUNE") == NULL)
    {
	Int q ;
	sh.lpend = (_Atomic Int *) malloc ((size_t) n * sizeof (Int)) ;
	sh.colver = (_Atomic unsigned *) calloc ((size_t) n,
						 sizeof (unsigned)) ;
	if (sh.lpend == NULL || sh.colver == NULL)
	{
	    free ((void *) sh.lpend) ;
    free (panel_start) ;
	    free ((void *) sh.colver) ;
	    sh.lpend = NULL ;
	    sh.colver = NULL ;
	}
	else
	{
	    for (q = 0 ; q < n ; q++)
	    {
		atomic_init (&sh.lpend [q], TRILINOS_KLU_EMPTY) ;
	    }
	}
    }

    for (t = 0 ; t < nthreads ; t++)
    {
	kls_klu_pipe_worker *W = &workers [t] ;
	W->S = S ;
	W->sh = &sh ;
	W->tid = t ;
	W->S.chunk_head = NULL ;
	W->S.chunk_used = 0 ;
	W->S.chunk_size = 0 ;
	W->S.lnz = 0 ;
	W->S.unz = 0 ;
	W->S.firstrow = 0 ;
	W->ubuf_i = (Int *) malloc ((size_t) n * sizeof (Int)) ;
	W->ubuf_x = (Entry *) malloc ((size_t) n * sizeof (Entry)) ;
	W->ap_ver = (unsigned *) malloc ((size_t) n * sizeof (unsigned)) ;
	W->copybuf = sh.lpend != NULL
	    ? (Unit *) malloc ((2 * (size_t) n + 4) * sizeof (Unit))
	    : NULL ;
	W->t_work = 0 ; W->t_spin = 0 ; W->t_final = 0 ;
	W->n_cols = 0 ; W->n_rounds = 0 ;
	W->kflops = 0.0 ;
	/* panel-mode workspace: width capped by a ~128MB/worker budget;
	   slot 0 aliases the worker's own X/ubuf/scratch pattern space */
	{
	    Int w, wcap = 1 ;
	    memset (W->pX, 0, sizeof (W->pX)) ;
	    memset (W->pLik, 0, sizeof (W->pLik)) ;
	    memset (W->pUbuf_i, 0, sizeof (W->pUbuf_i)) ;
	    memset (W->pUbuf_x, 0, sizeof (W->pUbuf_x)) ;
	    memset (W->pFlagW, 0, sizeof (W->pFlagW)) ;
	    W->pPFlag = NULL ;
	    W->pUlist = NULL ;
	    W->pB = NULL ;
	    W->pBcap = 0 ;
	    W->pRowPos = NULL ;
	    W->pRowGen = NULL ;
	    W->pRowList = NULL ;
	    W->pGen = 0 ;
	    if (panel_start != NULL)
	    {
size_t col_cap = (size_t) n / 4 + 16 ;
		size_t per_col = (size_t) n * sizeof (Entry) +
		    col_cap * (sizeof (Entry) + 2 * sizeof (Int)) ;
		size_t budget = (size_t) 128 << 20 ;
		{
		    const char *be = getenv ("KLS_KLU_PIPE_BUDGET_MB") ;
		    if (be != NULL && be [0] != '\0')
		    {
			long bm = atol (be) ;
			if (bm >= 64 && bm <= 4096)
			{
			    budget = (size_t) bm << 20 ;
			}
		    }
		}
		wcap = (Int) (budget / (per_col > 0 ? per_col : 1)) ;
		if (wcap > 32) wcap = 32 ;
		if (wcap < 1) wcap = 1 ;
		for (w = 0 ; w < wcap ; w++)
		{
		    if (w > 0)
		    {
			W->pX [w] = (Entry *) calloc ((size_t) n,
						      sizeof (Entry)) ;
			W->pLik [w] = (Int *) malloc (col_cap *
						      sizeof (Int)) ;
			W->pUbuf_i [w] = (Int *) malloc (col_cap *
							 sizeof (Int)) ;
			W->pUbuf_x [w] = (Entry *) malloc (col_cap *
							   sizeof (Entry)) ;
			if (W->pX [w] == NULL || W->pLik [w] == NULL ||
			    W->pUbuf_i [w] == NULL ||
			    W->pUbuf_x [w] == NULL)
			{
			    wcap = w ;
			    break ;
			}
		    }
		    W->pFlagW [w] = (Int *) malloc ((size_t) n *
						    sizeof (Int)) ;
		    if (W->pFlagW [w] == NULL)
		    {
			wcap = w > 0 ? w : 1 ;
			break ;
		    }
		    {
			Int q3 ;
			for (q3 = 0 ; q3 < n ; q3++)
			{
			    W->pFlagW [w][q3] = TRILINOS_KLU_EMPTY ;
			}
		    }
		}
		W->pPFlag = (Int *) malloc ((size_t) n * sizeof (Int)) ;
		W->pUlist = (Int *) malloc ((size_t) n * sizeof (Int)) ;
		W->pColCap = (Int) col_cap ;
		if (getenv ("KLS_KLU_PIPE_NOBUF") == NULL)
		{
		    W->pBcap = 4096 ;
		    W->pB = (Entry *) malloc (W->pBcap * (size_t) wcap *
					      sizeof (Entry)) ;
		    W->pRowPos = (Int *) malloc ((size_t) n * sizeof (Int)) ;
		    W->pRowGen = (Int *) calloc ((size_t) n, sizeof (Int)) ;
		    W->pRowList = (Int *) malloc ((size_t) n * sizeof (Int)) ;
		    W->pGen = 1 ;   /* calloc'd generations are 0 */
		    if (W->pB == NULL || W->pRowPos == NULL ||
			W->pRowGen == NULL || W->pRowList == NULL)
		    {
			free (W->pB) ; W->pB = NULL ;
			free (W->pRowPos) ; W->pRowPos = NULL ;
			free (W->pRowGen) ; W->pRowGen = NULL ;
			free (W->pRowList) ; W->pRowList = NULL ;
		    }
		}
		if (W->pUlist == NULL)
		{
		    free (W->pPFlag) ;
		    W->pPFlag = NULL ;
		}
		if (W->pPFlag != NULL)
		{
		    Int q2 ;
		    for (q2 = 0 ; q2 < n ; q2++)
		    {
			W->pPFlag [q2] = TRILINOS_KLU_EMPTY ;
		    }
		}
		else
		{
		    wcap = 1 ;
		}
	    }
	    W->pW = wcap ;
	    if (t == 0 && panel_start != NULL &&
		getenv ("KLS_KLU_PIPE_PROF") != NULL)
	    {
		Int fb = 0, pp ;
		for (pp = 0 ; pp < npanels ; pp++)
		{
		    if (panel_start [pp+1] - panel_start [pp] > wcap) fb++ ;
		}
		fprintf (stderr, "KLS pipe pW=%ld fallback_panels=%ld/%ld"
			 " buf=%d\n", (long) wcap, (long) fb,
			 (long) npanels, W->pB != NULL) ;
	    }
	}
	if (W->ubuf_i == NULL || W->ubuf_x == NULL || W->ap_ver == NULL ||
	    (sh.lpend != NULL && W->copybuf == NULL))
	{
	    spawn_failed = 1 ;
	}
	if (t > 0 && !spawn_failed)
	{
	    Int q ;
	    W->S.scratch = (Unit *) malloc ((2 * (size_t) n + 4) *
					    sizeof (Unit)) ;
	    W->S.X = (Entry *) malloc ((size_t) n * sizeof (Entry)) ;
	    W->S.Stack = (Int *) malloc ((size_t) n * sizeof (Int)) ;
	    W->S.Flag = (Int *) malloc ((size_t) n * sizeof (Int)) ;
	    W->S.Ap_pos = (Int *) malloc ((size_t) n * sizeof (Int)) ;
	    if (W->S.scratch == NULL || W->S.X == NULL ||
		W->S.Stack == NULL || W->S.Flag == NULL ||
		W->S.Ap_pos == NULL)
	    {
		spawn_failed = 1 ;
	    }
	    else
	    {
		for (q = 0 ; q < n ; q++)
		{
		    CLEAR (W->S.X [q]) ;
		    W->S.Flag [q] = TRILINOS_KLU_EMPTY ;
		}
	    }
	}
    }
    if (!spawn_failed)
    {
	if (kls_wall_prof) { kls_wall_t1 = kls_klu_now () ; }
	for (t = 1 ; t < nthreads ; t++)
	{
	    if (pthread_create (&tids [t], NULL, kls_klu_pipe_worker_main,
				&workers [t]) != 0)
	    {
		atomic_store (&sh.abort_flag, 1) ;
		nthreads = t ;
		spawn_failed = 1 ;
		break ;
	    }
	}
	kls_klu_pipe_worker_main (&workers [0]) ;
	for (t = 1 ; t < nthreads ; t++)
	{
	    pthread_join (tids [t], NULL) ;
	}
    }
    else
    {
	atomic_store (&sh.abort_flag, 1) ;
    }

    if (!spawn_failed &&
	!atomic_load_explicit (&sh.abort_flag, memory_order_acquire) &&
	atomic_load_explicit (&sh.prefix, memory_order_acquire) == n)
    {
	/* merge tallies; worker 0's state owns the template arena */
	S = workers [0].S ;
	for (t = 1 ; t < nthreads ; t++)
	{
	    S.lnz += workers [t].S.lnz ;
	    S.unz += workers [t].S.unz ;
	}
	for (t = 0 ; t < nthreads ; t++)
	{
	    Common->kls_kernel_flops += workers [t].kflops ;
	}
	S.cols_done = n ;
	{
	    const char *dump = getenv ("KLS_KLU_PIPE_DUMP") ;
	    if (dump != NULL)
	    {
		FILE *f = fopen (dump, "a") ;
		if (f != NULL)
		{
		    Int dk, dp ;
		    for (dk = 0 ; dk < n ; dk++)
		    {
			long lsum = 0, usum = 0 ;
			Unit *xp = S.colptr [dk] ;
			const int32_t *dLi = (const int32_t *) xp ;
			const int32_t *dUi = (const int32_t *)
			    (xp + KLS_UNITS32 (S.Llen [dk]) +
			     UNITS (Entry, S.Llen [dk])) ;
			for (dp = 0 ; dp < S.Llen [dk] ; dp++)
			{
			    lsum += dLi [dp] ;
			}
			for (dp = 0 ; dp < S.Ulen [dk] ; dp++)
			{
			    usum += dUi [dp] ;
			}
			fprintf (f, "%ld p=%ld l=%ld u=%ld ls=%ld us=%ld\n",
				 (long) (dk + k1), (long) S.P [dk],
				 (long) S.Llen [dk], (long) S.Ulen [dk],
				 lsum, usum) ;
		    }
		    fclose (f) ;
		}
	    }
	}
	if (kls_wall_prof) { kls_wall_t2 = kls_klu_now () ; }
	final_size = KLS_KLU_KERNEL_FINISH (&S) ;
	if (kls_wall_prof)
	{
	    fprintf (stderr, "KLS pipe wall: setup=%.2fs workers=%.2fs"
		     " finish=%.2fs\n", kls_wall_t1 - kls_wall_t0,
		     kls_wall_t2 - kls_wall_t1, kls_klu_now () - kls_wall_t2) ;
	}
	for (t = 0 ; t < nthreads ; t++)
	{
	    if (t > 0)
	    {
		KLS_KLU_KERNEL_CHUNKS_FREE (&workers [t].S) ;
		free (workers [t].S.scratch) ;
		free (workers [t].S.X) ;
		free (workers [t].S.Stack) ;
		free (workers [t].S.Flag) ;
		free (workers [t].S.Ap_pos) ;
	    }
	    free (workers [t].ubuf_i) ;
	    free (workers [t].ubuf_x) ;
	    free (workers [t].ap_ver) ;
	    free (workers [t].copybuf) ;
	{
	    Int w ;
	    for (w = 0 ; w < workers [t].pW ; w++)
	    {
		if (w > 0)
		{
		    free (workers [t].pX [w]) ;
		    free (workers [t].pLik [w]) ;
		    free (workers [t].pUbuf_i [w]) ;
		    free (workers [t].pUbuf_x [w]) ;
		}
		free (workers [t].pFlagW [w]) ;
	    }
	    free (workers [t].pPFlag) ;
	    free (workers [t].pUlist) ;
	    free (workers [t].pB) ;
	    free (workers [t].pRowPos) ;
	    free (workers [t].pRowGen) ;
	    free (workers [t].pRowList) ;
	}
	}
	free ((void *) sh.lpend) ;
    free (panel_start) ;
	free ((void *) sh.colver) ;
	free (S.colptr) ;
	*p_LU = S.LU ;
	*lnz = S.lnz ;
	*unz = S.unz ;
	return (final_size) ;
    }

    /* abort: free everything and signal the caller to retry serially */
    for (t = 0 ; t < nthreads ; t++)
    {
	if (t > 0)
	{
	    KLS_KLU_KERNEL_CHUNKS_FREE (&workers [t].S) ;
	    free (workers [t].S.scratch) ;
	    free (workers [t].S.X) ;
	    free (workers [t].S.Stack) ;
	    free (workers [t].S.Flag) ;
	    free (workers [t].S.Ap_pos) ;
	}
	free (workers [t].ubuf_i) ;
	free (workers [t].ubuf_x) ;
	free (workers [t].ap_ver) ;
	free (workers [t].copybuf) ;
	{
	    Int w ;
	    for (w = 1 ; w < workers [t].pW ; w++)
	    {
		free (workers [t].pX [w]) ;
		free (workers [t].pLik [w]) ;
		free (workers [t].pUbuf_i [w]) ;
		free (workers [t].pUbuf_x [w]) ;
	    }
	    free (workers [t].pPFlag) ;
	}
    }
    KLS_KLU_KERNEL_CHUNKS_FREE (&workers [0].S) ;
    free (S.scratch) ;
    free (S.colptr) ;
    if (Common->status == TRILINOS_KLU_OK)
    {
	Common->status = TRILINOS_KLU_SINGULAR ;
    }
    return (lusize) ;
}

/* Level-scheduled factorization: same contract as TRILINOS_KLU_kernel.
 * v0 runs the schedule serially through the chunked step (validates the
 * schedule and storage); the worker pool lands on top of this. */
size_t KLS_KLU_KERNEL_LEVELS
(
    Int n, Int Ap [ ], Int Ai [ ], Entry Ax [ ], Int Q [ ], size_t lusize,
    Int Pinv [ ], Int P [ ], Unit **p_LU, Entry Udiag [ ],
    Int Llen [ ], Int Ulen [ ], Int Lip [ ], Int Uip [ ],
    Int *lnz, Int *unz, Entry X [ ],
    Int Stack [ ], Int Flag [ ], Int Ap_pos [ ], Int Lpend [ ],
    Int k1, Int PSinv [ ], double Rs [ ],
    Int Offp [ ], Int Offi [ ], Entry Offx [ ],
    TRILINOS_KLU_common *Common
)
{
    KLS_KLU_KERNEL_STATE S ;
    Int k, lev, result ;
    size_t final_size ;
    Int *parent, *anc, *prevrow, *level, *level_ptr, *level_cols ;
    Int nlevels, poff ;
    const int prof = (getenv ("KLS_KLU_LEVELS_PROF") != NULL) ;
    double t_sched0 = prof ? kls_klu_now () : 0, t_sched = 0 ;

    parent = (Int *) TRILINOS_KLU_malloc ((size_t) (6*n + 2), sizeof (Int),
					  Common) ;
    if (parent == NULL)
    {
	Common->status = TRILINOS_KLU_OUT_OF_MEMORY ;
	return (lusize) ;
    }
    anc = parent + n ;
    prevrow = anc + n ;
    level = prevrow + n ;
    level_ptr = level + n ;         /* n+1 slots */
    level_cols = level_ptr + n + 1 ;

    kls_klu_block_coletree (n, Ap, Ai, Q, k1, PSinv, parent, anc, prevrow) ;

    /* levels: parent[k] > k, one ascending pass */
    nlevels = 0 ;
    for (k = 0 ; k < n ; k++)
    {
	level [k] = 0 ;
    }
    for (k = 0 ; k < n ; k++)
    {
	if (parent [k] != TRILINOS_KLU_EMPTY &&
	    level [k] + 1 > level [parent [k]])
	{
	    level [parent [k]] = level [k] + 1 ;
	}
	if (level [k] + 1 > nlevels)
	{
	    nlevels = level [k] + 1 ;
	}
    }
    for (lev = 0 ; lev <= n ; lev++)
    {
	level_ptr [lev] = 0 ;
    }
    for (k = 0 ; k < n ; k++)
    {
	level_ptr [level [k] + 1]++ ;
    }
    for (lev = 0 ; lev < n ; lev++)
    {
	level_ptr [lev + 1] += level_ptr [lev] ;
    }
    {
	Int *cursor = anc ;   /* reuse workspace */
	for (lev = 0 ; lev < nlevels ; lev++)
	{
	    cursor [lev] = level_ptr [lev] ;
	}
	for (k = 0 ; k < n ; k++)
	{
	    level_cols [cursor [level [k]]++] = k ;
	}
    }

    /* Offp prefix so column writes are disjoint under any order */
    poff = Offp [k1] ;
    for (k = 0 ; k < n ; k++)
    {
	Int oldcol = Q [k + k1] ;
	Int p2 ;
	Offp [k + k1] = poff ;
	for (p2 = Ap [oldcol] ; p2 < Ap [oldcol+1] ; p2++)
	{
	    if (PSinv [Ai [p2]] - k1 < 0)
	    {
		poff++ ;
	    }
	}
    }
    Offp [n + k1] = poff ;

    S.n = n ;
    S.Ap = Ap ; S.Ai = Ai ; S.Ax = Ax ; S.Q = Q ;
    S.lusize = lusize ;
    S.Pinv = Pinv ; S.P = P ;
    S.LU = *p_LU ;
    S.Udiag = Udiag ;
    S.Llen = Llen ; S.Ulen = Ulen ; S.Lip = Lip ; S.Uip = Uip ;
    S.X = X ; S.Stack = Stack ; S.Flag = Flag ;
    S.Ap_pos = Ap_pos ; S.Lpend = Lpend ;
    S.k1 = k1 ; S.PSinv = PSinv ; S.Rs = Rs ;
    S.Offp = Offp ; S.Offi = Offi ; S.Offx = Offx ;
    S.Common = Common ;
    S.scale = Common->scale ;
    S.tol = Common->tol ;
    S.memgrow = Common->memgrow ;
    S.no_prune = 1 ;
    S.diag_claim = 0 ;
    S.chunked_prune = 0 ;
    S.idx32 = 0 ;
    S.pack_keep_row_indices = 0 ;
    S.chunk_head = NULL ; S.chunk_used = 0 ; S.chunk_size = 0 ;
    if (prof)
    {
	t_sched = kls_klu_now () - t_sched0 ;
    }
    S.scratch = (Unit *) malloc ((2 * (size_t) n + 4) * sizeof (Unit)) ;
    S.colptr = (Unit **) calloc ((size_t) n, sizeof (Unit *)) ;
    if (S.colptr == NULL)
    {
	TRILINOS_KLU_free (parent, (size_t) (6*n + 2), sizeof (Int), Common) ;
	Common->status = TRILINOS_KLU_OUT_OF_MEMORY ;
	return (lusize) ;
    }

    KLS_KLU_KERNEL_INIT (&S) ;

    {
	int nthreads = 1 ;
	{
	    const char *tenv = getenv ("KLS_KLU_LEVELS") ;
	    if (tenv != NULL && tenv [0] != '\0')
	    {
		nthreads = atoi (tenv) ;
	    }
	    if (nthreads < 1)
	    {
		nthreads = 1 ;
	    }
	    if (nthreads > 16)
	    {
		nthreads = 16 ;
	    }
	}
	if (nthreads > 1)
	{
	    /* parallel phase: diagonal claims only; deferred columns and
	       their ancestor paths fall to the serial cleanup below */
	    kls_klu_par_shared sh ;
	    kls_klu_par_worker workers [16] ;
	    pthread_t tids [16] ;
	    int t, spawn_failed = 0 ;
	    Int cleanup_failed = 0 ;
	    Int *super_ptr = NULL ;
	    char *super_serial = NULL ;
	    sh.level_ptr = level_ptr ;
	    sh.level_cols = level_cols ;
	    sh.nlevels = nlevels ;
	    sh.n = n ;
	    sh.parent = parent ;
	    /* batch consecutive narrow levels into single-thread
	       super-levels: chain-heavy etrees otherwise pay a barrier
	       per near-empty level (rajat25: parallel 1.55s vs serial
	       0.87 was almost pure barrier churn) */
	    super_ptr = (Int *) malloc (((size_t) nlevels + 1) *
					sizeof (Int)) ;
	    super_serial = (char *) malloc ((size_t) nlevels + 1) ;
	    if (super_ptr != NULL && super_serial != NULL)
	    {
		const Int wide = 4 * nthreads ;
		Int lev2 = 0, ns = 0 ;
		while (lev2 < nlevels)
		{
		    Int width = level_ptr [lev2 + 1] - level_ptr [lev2] ;
		    super_ptr [ns] = lev2 ;
		    if (width >= wide)
		    {
			super_serial [ns] = 0 ;
			lev2++ ;
		    }
		    else
		    {
			super_serial [ns] = 1 ;
			while (lev2 < nlevels &&
			       level_ptr [lev2 + 1] - level_ptr [lev2] <
				 wide)
			{
			    lev2++ ;
			}
		    }
		    ns++ ;
		}
		super_ptr [ns] = nlevels ;
		sh.super_ptr = super_ptr ;
		sh.super_serial = super_serial ;
		sh.nsuper = ns ;
	    }
	    else
	    {
		free (super_ptr) ;
		free (super_serial) ;
		super_ptr = NULL ;
		super_serial = NULL ;
	    }
	    sh.defer = (_Atomic char *) calloc ((size_t) n, 1) ;
	    atomic_init (&sh.abort_flag, 0) ;
	    atomic_init (&sh.defer_count, 0) ;
	    sh.defer_limit = (long) n / 10 + 64 ;
	    sh.nthreads = nthreads ;
	    sh.prof = prof ;
	    sh.wide_secs = 0 ;
	    sh.serial_secs = 0 ;
	    if (prof && super_ptr != NULL)
	    {
		Int nw = 0, cw = 0, maxw = 0, sl2 ;
		for (sl2 = 0 ; sl2 < sh.nsuper ; sl2++)
		{
		    Int c = level_ptr [super_ptr [sl2 + 1]] -
			    level_ptr [super_ptr [sl2]] ;
		    if (!super_serial [sl2]) { nw++ ; cw += c ; }
		}
		for (sl2 = 0 ; sl2 < nlevels ; sl2++)
		{
		    Int w = level_ptr [sl2 + 1] - level_ptr [sl2] ;
		    if (w > maxw) { maxw = w ; }
		}
		fprintf (stderr, "LEVELSPROF n=%ld nlev=%ld nsuper=%ld "
			 "wide_supers=%ld wide_cols=%ld serial_cols=%ld "
			 "maxwidth=%ld sched=%.3fs\n",
			 (long) n, (long) nlevels, (long) sh.nsuper,
			 (long) nw, (long) cw, (long) (n - cw), (long) maxw,
			 t_sched) ;
	    }
	    if (sh.defer == NULL || super_ptr == NULL ||
		pthread_barrier_init (&sh.barrier, NULL,
				      (unsigned) nthreads) != 0)
	    {
		free ((void *) sh.defer) ;
		free (super_ptr) ;
		free (super_serial) ;
		nthreads = 1 ;
	    }
	    else
	    {
	    for (t = 0 ; t < nthreads ; t++)
	    {
		kls_klu_par_worker *W = &workers [t] ;
		W->S = S ;
		W->sh = &sh ;
		W->tid = t ;
		W->S.diag_claim = 1 ;
		W->S.chunk_head = NULL ;
		W->S.chunk_used = 0 ;
		W->S.chunk_size = 0 ;
		W->S.lnz = 0 ;
		W->S.unz = 0 ;
		W->S.firstrow = 0 ;
		if (t > 0)
		{
		    W->S.scratch = (Unit *) malloc ((2 * (size_t) n + 4) *
						    sizeof (Unit)) ;
		    if (W->S.scratch == NULL)
		    {
			spawn_failed = 1 ;
			break ;
		    }
		    Int q ;
		    W->S.X = (Entry *) malloc ((size_t) n * sizeof (Entry)) ;
		    W->S.Stack = (Int *) malloc ((size_t) n * sizeof (Int)) ;
		    W->S.Flag = (Int *) malloc ((size_t) n * sizeof (Int)) ;
		    W->S.Ap_pos = (Int *) malloc ((size_t) n * sizeof (Int)) ;
		    W->S.Lpend = (Int *) malloc ((size_t) n * sizeof (Int)) ;
		    if (W->S.X == NULL || W->S.Stack == NULL ||
			W->S.Flag == NULL || W->S.Ap_pos == NULL ||
			W->S.Lpend == NULL)
		    {
			spawn_failed = 1 ;
			break ;
		    }
		    for (q = 0 ; q < n ; q++)
		    {
			CLEAR (W->S.X [q]) ;
			W->S.Flag [q] = TRILINOS_KLU_EMPTY ;
			W->S.Lpend [q] = TRILINOS_KLU_EMPTY ;
		    }
		}
	    }
	    if (!spawn_failed)
	    {
		for (t = 1 ; t < nthreads ; t++)
		{
		    if (pthread_create (&tids [t], NULL,
					kls_klu_par_worker_main,
					&workers [t]) != 0)
		    {
			atomic_store (&sh.abort_flag, 1) ;
			spawn_failed = 1 ;
			nthreads = t ;   /* join only the created ones */
			break ;
		    }
		}
		kls_klu_par_worker_main (&workers [0]) ;
		for (t = 1 ; t < nthreads ; t++)
		{
		    pthread_join (tids [t], NULL) ;
		}
	    }
	    pthread_barrier_destroy (&sh.barrier) ;
	    if (!spawn_failed &&
		!atomic_load_explicit (&sh.abort_flag, memory_order_acquire))
	    {
		/* merge worker tallies; the template arena is worker 0's */
		S = workers [0].S ;
		S.diag_claim = 0 ;
		for (t = 1 ; t < nthreads ; t++)
		{
		    S.lnz += workers [t].S.lnz ;
		    S.unz += workers [t].S.unz ;
		}
		/* serial cleanup: deferred columns in ascending order.
		   Diagonal claims can assign rows differently than the
		   classic pivot order would, so a deferred column may
		   find its candidates taken; that is pivot STARVATION,
		   not singularity - take the zero-pivot rescue (the
		   caller's probe/nudge/refinement machinery arbitrates,
		   as on the predicted path) instead of halting. */
		{
		    Int saved_halt = Common->halt_if_singular ;
		    double t_cl0 = prof ? kls_klu_now () : 0 ;
		    Common->halt_if_singular = 0 ;
		    S.chunked_prune = 1 ;
		    for (k = 0 ; k < n ; k++)
		    {
			if (!atomic_load_explicit (&sh.defer [k],
						   memory_order_acquire))
			{
			    continue ;
			}
			result = KLS_KLU_KERNEL_STEP (&S, k) ;
			if (result != 0)
			{
			    cleanup_failed = 1 ;
			    break ;
			}
		    }
		    Common->halt_if_singular = saved_halt ;
		    if (prof)
		    {
			fprintf (stderr, "LEVELSPROF wide=%.3fs serial=%.3fs "
				 "cleanup=%.3fs deferred=%ld failed=%ld\n",
				 sh.wide_secs, sh.serial_secs,
				 kls_klu_now () - t_cl0,
				 (long) atomic_load (&sh.defer_count),
				 (long) cleanup_failed) ;
		    }
		}
		S.cols_done = cleanup_failed ? S.cols_done : n ;
		for (t = 1 ; t < nthreads ; t++)
		{
		    if (workers [t].S.X != NULL)
		    {
			/* pack first needs the chunks alive; only the
			   scratch arrays can go now */
		    }
		}
		if (!cleanup_failed)
		{
		    final_size = KLS_KLU_KERNEL_FINISH (&S) ;
		}
		else
		{
		    S.pack_keep_row_indices = 1 ;
		    final_size = KLS_KLU_KERNEL_FINISH (&S) ;
		}
		for (t = 1 ; t < nthreads ; t++)
		{
		    KLS_KLU_KERNEL_CHUNKS_FREE (&workers [t].S) ;
		    free (workers [t].S.X) ;
		    free (workers [t].S.Stack) ;
		    free (workers [t].S.Flag) ;
		    free (workers [t].S.Ap_pos) ;
		    free (workers [t].S.Lpend) ;
		    free (workers [t].S.scratch) ;
		}
		free ((void *) sh.defer) ;
		free (super_ptr) ;
		free (super_serial) ;
		free (S.colptr) ;
		TRILINOS_KLU_free (parent, (size_t) (6*n + 2), sizeof (Int),
				   Common) ;
		*p_LU = S.LU ;
		*lnz = S.lnz ;
		*unz = S.unz ;
		if (cleanup_failed)
		{
		    return (final_size) ;
		}
		return (final_size) ;
	    }
	    if (prof)
	    {
		fprintf (stderr, "LEVELSPROF ABORT spawn_failed=%d "
			 "deferred=%ld limit=%ld wide=%.3fs serial=%.3fs\n",
			 spawn_failed, (long) atomic_load (&sh.defer_count),
			 sh.defer_limit, sh.wide_secs, sh.serial_secs) ;
	    }
	    /* spawn/abort: free worker scratch and fall through serial */
	    for (t = 1 ; t < nthreads ; t++)
	    {
		KLS_KLU_KERNEL_CHUNKS_FREE (&workers [t].S) ;
		free (workers [t].S.X) ;
		free (workers [t].S.Stack) ;
		free (workers [t].S.Flag) ;
		free (workers [t].S.Ap_pos) ;
		free (workers [t].S.Lpend) ;
		free (workers [t].S.scratch) ;
	    }
	    KLS_KLU_KERNEL_CHUNKS_FREE (&workers [0].S) ;
	    free ((void *) sh.defer) ;
	    free (super_ptr) ;
	    free (super_serial) ;
	    KLS_KLU_KERNEL_INIT (&S) ;
	    S.chunk_head = NULL ;
	    S.chunk_used = 0 ;
	    S.chunk_size = 0 ;
	    }
	}
    }

    for (lev = 0 ; lev < nlevels ; lev++)
    {
	Int pos ;
	for (pos = level_ptr [lev] ; pos < level_ptr [lev + 1] ; pos++)
	{
	    k = level_cols [pos] ;
	    result = KLS_KLU_KERNEL_STEP (&S, k) ;
	    if (result != 0)
	    {
		S.pack_keep_row_indices = 1 ;
		final_size = KLS_KLU_KERNEL_FINISH (&S) ;
		free (S.scratch) ;
		free (S.colptr) ;
		TRILINOS_KLU_free (parent, (size_t) (6*n + 2), sizeof (Int),
				   Common) ;
		*p_LU = S.LU ;
		*lnz = S.lnz ;
		*unz = S.unz ;
		return (final_size) ;
	    }
	}
    }

    final_size = KLS_KLU_KERNEL_FINISH (&S) ;
    free (S.scratch) ;
    free (S.colptr) ;
    TRILINOS_KLU_free (parent, (size_t) (6*n + 2), sizeof (Int), Common) ;
    *p_LU = S.LU ;
    *lnz = S.lnz ;
    *unz = S.unz ;
    return (final_size) ;
}

size_t TRILINOS_KLU_kernel   /* final size of LU on output */
(
    /* input, not modified */
    Int n,	    /* A is n-by-n */
    Int Ap [ ],	    /* size n+1, column pointers for A */
    Int Ai [ ],	    /* size nz = Ap [n], row indices for A */
    Entry Ax [ ],   /* size nz, values of A */
    Int Q [ ],	    /* size n, optional input permutation */
    size_t lusize,  /* initial size of LU on input */

    /* output, not defined on input */
    Int Pinv [ ],   /* size n, inverse row permutation, where Pinv [i] = k if
		     * row i is the kth pivot row */
    Int P [ ],	    /* size n, row permutation, where P [k] = i if row i is the
		     * kth pivot row. */
    Unit **p_LU,	/* LU array, size lusize on input */
    Entry Udiag [ ],	/* size n, diagonal of U */
    Int Llen [ ],       /* size n, column length of L */
    Int Ulen [ ],	/* size n, column length of U */
    Int Lip [ ],	/* size n, column pointers for L */
    Int Uip [ ],	/* size n, column pointers for U */
    Int *lnz,		/* size of L*/
    Int *unz,		/* size of U*/
    /* workspace, not defined on input */
    Entry X [ ],    /* size n, undefined on input, zero on output */

    /* workspace, not defined on input or output */
    Int Stack [ ],  /* size n */
    Int Flag [ ],   /* size n */
    Int Ap_pos [ ],	/* size n */

    /* other workspace: */
    Int Lpend [ ],		    /* size n workspace, for pruning only */

    /* inputs, not modified on output */
    Int k1,	    	/* the block of A is from k1 to k2-1 */
    Int PSinv [ ],  	/* inverse of P from symbolic factorization */
    double Rs [ ],  	/* scale factors for A */

    /* inputs, modified on output */
    Int Offp [ ],   /* off-diagonal matrix (modified by this routine) */
    Int Offi [ ],
    Entry Offx [ ],
    /* --------------- */
    TRILINOS_KLU_common *Common
)
{
    KLS_KLU_KERNEL_STATE S ;
    Int k, result ;
    size_t final_size ;

    ASSERT (Common != NULL) ;
    kls_construct_prof = getenv ("KLS_CONSTRUCT_PROF") != NULL ;
    kls_pipe_phase_prof = getenv ("KLS_KLU_PIPE_PHASES") != NULL ;
    kls_pipe_dense_panels = getenv ("KLS_KLU_PIPE_DENSE") != NULL ;
    Common->kls_dense_panels = 0 ;
    if (n >= 4096 && getenv ("KLS_SN_STATS") != NULL)
    {
	/* supernodal-first-factor feasibility probe: fundamental
	   supernode partition over the column etree (parent[j]==j+1
	   chains, width cap 128).  Statistics only; no behavior. */
	Int *sp = (Int *) malloc ((size_t) (3 * n) * sizeof (Int)) ;
	if (sp != NULL)
	{
	    Int *par = sp, *anc = sp + n, *prv = sp + 2*n ;
	    Int k, nsn = 0, width = 1, singles = 0, maxw = 1 ;
	    long wsum = 0, w2sum = 0 ;
	    kls_klu_block_coletree (n, Ap, Ai, Q, k1, PSinv, par, anc, prv) ;
	    for (k = 1 ; k <= n ; k++)
	    {
		if (k < n && par [k-1] == k && width < 128)
		{
		    width++ ;
		    continue ;
		}
		nsn++ ;
		wsum += width ;
		w2sum += (long) width * width ;
		if (width > maxw) maxw = width ;
		if (width == 1) singles++ ;
		width = 1 ;
	    }
	    fprintf (stderr, "KLS snstats: n=%ld nsn=%ld avgw=%.1f "
		     "rmsw=%.1f maxw=%ld singles=%ld (%.0f%%)\n",
		     (long) n, (long) nsn, (double) wsum / (nsn > 0 ? nsn : 1),
		     nsn > 0 ? __builtin_sqrt ((double) w2sum / nsn) : 0.0,
		     (long) maxw, (long) singles,
		     100.0 * singles / (nsn > 0 ? nsn : 1)) ;
	    /* A'A-bound fill count (Gilbert-Ng-Peyton over the column
	       etree with row subtrees): the panel-pivoting bound's
	       affordability vs the true LU fill.  Skeleton-free upper
	       variant: count col j's bound rows as the union of row
	       subtrees - via the classic prevleaf/first-descendant LCA
	       skip counts. */
	    {
		Int *first = (Int *) malloc ((size_t) (4 * n) * sizeof (Int)) ;
		if (first != NULL)
		{
		    Int *prevleaf = first + n ;
		    Int *setparent = first + 2*n ;
		    Int *count = first + 3*n ;
		    Int j, kk, p2 ;
		    double bound_fill = 0.0 ;
		    /* first descendant via postorder-free approximation:
		       process columns ascending (etree children < parent
		       for the column etree of an ordered matrix) */
		    for (j = 0 ; j < n ; j++)
		    {
			first [j] = j ;
			prevleaf [j] = TRILINOS_KLU_EMPTY ;
			setparent [j] = j ;
			count [j] = 1 ;  /* diagonal */
		    }
		    for (j = 0 ; j < n ; j++)
		    {
			if (par [j] != TRILINOS_KLU_EMPTY &&
			    first [par [j]] > first [j])
			{
			    first [par [j]] = first [j] ;
			}
		    }
		    /* rows: each structural entry (i in column j of the
		       permuted block) contributes count increments at
		       LCA skips along row i's subtree leaves */
		    for (j = 0 ; j < n ; j++)
		    {
			Int kglobal = j + k1 ;
			Int oldcol = Q [kglobal] ;
			for (p2 = Ap [oldcol] ; p2 < Ap [oldcol+1] ; p2++)
			{
			    Int i = PSinv [Ai [p2]] - k1 ;
			    if (i < 0 || j <= i)
			    {
				continue ;   /* upper/diag or off-block */
			    }
			    /* column j sees row i (i < j): leaf j in row
			       i's subtree.  Skip-count: add path from j
			       up to the previous leaf's LCA */
			    {
				Int q = prevleaf [i] ;
				Int lca ;
				if (q == TRILINOS_KLU_EMPTY)
				{
				    lca = i ;
				}
				else
				{
				    /* find root of q's set */
				    Int r0 = q ;
				    while (setparent [r0] != r0)
				    {
					r0 = setparent [r0] ;
				    }
				    lca = r0 ;
				    /* path compress */
				    while (setparent [q] != r0)
				    {
					Int nx = setparent [q] ;
					setparent [q] = r0 ;
					q = nx ;
				    }
				}
				{
				    Int t ;
				    for (t = j ; t != lca &&
					 t != TRILINOS_KLU_EMPTY ;
					 t = par [t])
				    {
					count [t]++ ;
				    }
				}
				prevleaf [i] = j ;
			    }
			}
			/* union j into parent's set after processing */
			if (par [j] != TRILINOS_KLU_EMPTY)
			{
			    setparent [j] = par [j] ;
			}
		    }
		    for (j = 0 ; j < n ; j++)
		    {
			bound_fill += (double) count [j] ;
		    }
		    fprintf (stderr,
			     "KLS snstats: ata-bound fill %.3e\n",
			     bound_fill) ;
		    free (first) ;
		}
	    }
	    free (sp) ;
	}
    }
    if (n >= 512 && n < 2147483647 &&
	(getenv ("KLS_KLU_PIPE") != NULL || kls_klu_pipe_threads > 0))
    {
	int pipe_threads = getenv ("KLS_KLU_PIPE") != NULL
	    ? atoi (getenv ("KLS_KLU_PIPE"))
	    : kls_klu_pipe_threads ;
	size_t pipe_size = KLS_KLU_KERNEL_PIPE (n, Ap, Ai, Ax, Q, lusize,
	    Pinv, P, p_LU, Udiag, Llen, Ulen, Lip, Uip, lnz, unz, X, Stack,
	    Flag, Ap_pos, Lpend, k1, PSinv, Rs, Offp, Offi, Offx, Common,
	    pipe_threads) ;
	if (Common->status == TRILINOS_KLU_OK)
	{
	    return (pipe_size) ;
	}
	/* failed pipe run: drop any partial worker sums; the serial
	   retry below re-accumulates this block completely */
	/* singular or resource failure in the pipeline: retry the block
	   with the classic serial kernel (X/Flag state is clean; the
	   pipeline frees its own arenas on abort) */
	Common->status = TRILINOS_KLU_OK ;
	Common->numerical_rank = TRILINOS_KLU_EMPTY ;
	Common->singular_col = TRILINOS_KLU_EMPTY ;
	/* the serial path below re-runs KLS_KLU_KERNEL_INIT, which
	   restores X/Flag/Lpend and the P/Pinv flipped identity */
    }
    if (n >= 512 && getenv ("KLS_KLU_LEVELS") != NULL)
    {
	size_t par_size = KLS_KLU_KERNEL_LEVELS (n, Ap, Ai, Ax, Q, lusize,
	    Pinv, P, p_LU, Udiag, Llen, Ulen, Lip, Uip, lnz, unz, X, Stack,
	    Flag, Ap_pos, Lpend, k1, PSinv, Rs, Offp, Offi, Offx, Common) ;
	if (Common->status == TRILINOS_KLU_OK)
	{
	    return (par_size) ;
	}
	/* pivot starvation (diagonal claims can strand a deferred
	   column) or another parallel-path failure: retry the block
	   with the classic serial kernel on the same (possibly
	   repacked) LU allocation */
	Common->status = TRILINOS_KLU_OK ;
	Common->numerical_rank = TRILINOS_KLU_EMPTY ;
	Common->singular_col = TRILINOS_KLU_EMPTY ;
	lusize = par_size ;
    }
    S.n = n ;
    S.Ap = Ap ;
    S.Ai = Ai ;
    S.Ax = Ax ;
    S.Q = Q ;
    S.lusize = lusize ;
    S.Pinv = Pinv ;
    S.P = P ;
    S.LU = *p_LU ;
    S.Udiag = Udiag ;
    S.Llen = Llen ;
    S.Ulen = Ulen ;
    S.Lip = Lip ;
    S.Uip = Uip ;
    S.X = X ;
    S.Stack = Stack ;
    S.Flag = Flag ;
    S.Ap_pos = Ap_pos ;
    S.Lpend = Lpend ;
    S.k1 = k1 ;
    S.PSinv = PSinv ;
    S.Rs = Rs ;
    S.Offp = Offp ;
    S.Offi = Offi ;
    S.Offx = Offx ;
    S.Common = Common ;
    S.scale = Common->scale ;
    S.tol = Common->tol ;
    S.memgrow = Common->memgrow ;
    S.no_prune = 0 ;
    S.diag_claim = 0 ;
    S.chunked_prune = 0 ;
    S.idx32 = 0 ;
    S.pack_keep_row_indices = 0 ;
    S.colptr = NULL ;
    S.scratch = NULL ;
    S.chunk_head = NULL ;
    S.chunk_used = 0 ;
    S.chunk_size = 0 ;
    if (getenv ("KLS_KLU_CHUNKED") != NULL)
    {
	/* validation mode: run every factorization through the chunked
	   storage + pack path the parallel driver will use */
	S.colptr = (Unit **) TRILINOS_KLU_malloc ((size_t) n,
						  sizeof (Unit *), Common) ;
	S.scratch = (Unit *) malloc ((2 * (size_t) n + 4) * sizeof (Unit)) ;
	if (S.colptr != NULL && S.scratch != NULL)
	{
	    S.no_prune = 1 ;
	}
	else
	{
	    if (S.colptr != NULL)
	    {
		TRILINOS_KLU_free (S.colptr, (size_t) n, sizeof (Unit *),
				   Common) ;
	    }
	    free (S.scratch) ;
	    S.colptr = NULL ;
	    S.scratch = NULL ;
	}
    }

    KLS_KLU_KERNEL_INIT (&S) ;

    for (k = 0 ; k < n ; k++)
    {
	result = KLS_KLU_KERNEL_STEP (&S, k) ;
	if (result != 0)
	{
	    size_t partial = S.lusize ;
	    if (S.colptr != NULL)
	    {
		/* leave a classic-layout partial LU for the restart and
		   diagnostic paths (row indices, as mid-factor); the
		   failed column's L pattern is valid and inspected by
		   the restart heuristics - include it with an empty U */
		if (S.cols_done == k && k < n)
		{
		    S.Ulen [k] = 0 ;
		    S.Uip [k] = S.Lip [k] + UNITS (Int, S.Llen [k]) +
				UNITS (Entry, S.Llen [k]) ;
		    S.cols_done = k + 1 ;
		}
		S.pack_keep_row_indices = 1 ;
		partial = KLS_KLU_KERNEL_FINISH (&S) ;
		TRILINOS_KLU_free (S.colptr, (size_t) n, sizeof (Unit *),
				   Common) ;
	    }
	    *p_LU = S.LU ;
	    *lnz = S.lnz ;
	    *unz = S.unz ;
	    return (partial) ;
	}
    }

    final_size = KLS_KLU_KERNEL_FINISH (&S) ;
    if (S.colptr != NULL)
    {
	TRILINOS_KLU_free (S.colptr, (size_t) n, sizeof (Unit *), Common) ;
	free (S.scratch) ;
    }
    *p_LU = S.LU ;
    *lnz = S.lnz ;
    *unz = S.unz ;
    return (final_size) ;
}
