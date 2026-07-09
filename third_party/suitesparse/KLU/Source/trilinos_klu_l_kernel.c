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
    for (p = top, i = 0 ; p < n ; p++, i++)
    {
	j = S->Stack [p] ;
	Ui [i] = Pinv [j] ;
	Ux [i] = S->X [j] ;
	CLEAR (S->X [j]) ;
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
    }
    return (NULL) ;
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
    S.pack_keep_row_indices = 0 ;
    S.chunk_head = NULL ; S.chunk_used = 0 ; S.chunk_size = 0 ;
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
