/* ========================================================================== */
/* === BTF_ORDER ============================================================ */
/* ========================================================================== */

/* Find a permutation P and Q to permute a square sparse matrix into upper block
 * triangular form.  A(P,Q) will contain a zero-free diagonal if A has
 * structural full-rank.  Otherwise, the number of nonzeros on the diagonal of
 * A(P,Q) will be maximized, and will equal the structural rank of A.
 *
 * Q[k] will be "flipped" if a zero-free diagonal was not found.  Q[k] will be
 * negative, and j = TRILINOS_BTF_UNFLIP (Q [k]) gives the corresponding permutation.
 *
 * R defines the block boundaries of A(P,Q).  The kth block consists of rows
 * and columns R[k] to R[k+1]-1.
 *
 * If maxwork > 0 on input, then the work performed in btf_maxtrans is limited
 * to maxwork*nnz(A) (excluding the "cheap match" phase, which can take another
 * nnz(A) work).  On output, the work parameter gives the actual work performed,
 * or -1 if the limit was reached.  In the latter case, the diagonal of A(P,Q)
 * might not be zero-free, and the number of nonzeros on the diagonal of A(P,Q)
 * might not be equal to the structural rank.
 *
 * See btf.h for more details.
 *
 * Copyright (c) 2004-2007.  Tim Davis, University of Florida,
 * with support from Sandia National Laboratories.  All Rights Reserved.
 */

/* This file should make the long int version of TRILINOS_BTF */
#define DLONG 1

#include <stdlib.h>
#include "trilinos_btf_decl.h"
#include "trilinos_btf_internal.h"

/* This function only operates on square matrices (either structurally full-
 * rank, or structurally rank deficient). */

/* Hopcroft-Karp maximum bipartite matching: BFS phases + layered DFS,
 * O(sqrt(n) * nnz) worst case vs the depth-first augmenting maxtrans,
 * whose pathological chains cost ~9s on mac_econ-class patterns
 * (89% structurally missing diagonal).  Produces an identical-size
 * maximum matching with the same output contract as maxtrans:
 * Match [i] = column matched to row i, EMPTY if unmatched; returns
 * the match count.  Workspace: 4n from Work plus one allocation. */
static Int kls_hk_maxtrans
(
    Int n,
    Int Ap [ ],
    Int Ai [ ],
    Int Match [ ],      /* size n: row -> col */
    Int Work [ ]        /* size >= 4n */
)
{
    Int *match_col = Work ;          /* col -> row */
    Int *dist = Work + n ;           /* BFS layer, read-only in the DFS */
    Int *queue = Work + 2*n ;
    Int *dfs_pos = Work + 3*n ;
    Int *used = Work + 4*n ;         /* per-phase visited marker */
    Int i, j, p, nmatch, head, tail, found, top ;
    Int *stack = (Int *) malloc ((size_t) n * sizeof (Int)) ;
    if (stack == NULL)
    {
	return (TRILINOS_BTF_EMPTY) ;
    }
    for (i = 0 ; i < n ; i++)
    {
	Match [i] = TRILINOS_BTF_EMPTY ;
	match_col [i] = TRILINOS_BTF_EMPTY ;
    }
    /* greedy warm start */
    nmatch = 0 ;
    for (j = 0 ; j < n ; j++)
    {
	for (p = Ap [j] ; p < Ap [j+1] ; p++)
	{
	    i = Ai [p] ;
	    if (Match [i] == TRILINOS_BTF_EMPTY)
	    {
		Match [i] = j ;
		match_col [j] = i ;
		nmatch++ ;
		break ;
	    }
	}
    }
    for ( ; ; )
    {
	/* BFS phase: layer the unmatched columns */
	head = 0 ;
	tail = 0 ;
	found = 0 ;
	for (j = 0 ; j < n ; j++)
	{
	    used [j] = 0 ;
	    if (match_col [j] == TRILINOS_BTF_EMPTY)
	    {
		dist [j] = 0 ;
		queue [tail++] = j ;
	    }
	    else
	    {
		dist [j] = -1 ;
	    }
	}
	while (head < tail)
	{
	    j = queue [head++] ;
	    for (p = Ap [j] ; p < Ap [j+1] ; p++)
	    {
		Int jm = Match [Ai [p]] ;
		if (jm == TRILINOS_BTF_EMPTY)
		{
		    found = 1 ;
		}
		else if (dist [jm] < 0)
		{
		    dist [jm] = dist [j] + 1 ;
		    queue [tail++] = jm ;
		}
	    }
	}
	if (!found)
	{
	    break ;
	}
	/* layered DFS from each unmatched column */
	for (j = 0 ; j < n ; j++)
	{
	    if (match_col [j] != TRILINOS_BTF_EMPTY)
	    {
		continue ;
	    }
	    top = 0 ;
	    stack [top] = j ;
	    dfs_pos [j] = Ap [j] ;
	    while (top >= 0)
	    {
		Int jc = stack [top] ;
		Int advanced = 0 ;
		for (p = dfs_pos [jc] ; p < Ap [jc+1] ; p++)
		{
		    Int r = Ai [p] ;
		    Int jm = Match [r] ;
		    if (jm == TRILINOS_BTF_EMPTY)
		    {
			/* augmenting path found: flip along the stack */
			Int t ;
			dfs_pos [jc] = p + 1 ;
			Match [r] = jc ;
			match_col [jc] = r ;
			for (t = top - 1 ; t >= 0 ; t--)
			{
			    Int jprev = stack [t] ;
			    Int rprev = match_col [stack [t + 1]] ;
			    (void) jprev ; (void) rprev ;
			}
			/* rebuild the flips pairwise down the stack */
			for (t = top ; t > 0 ; t--)
			{
			    Int jchild = stack [t] ;
			    Int jparent = stack [t - 1] ;
			    Int pparent = dfs_pos [jparent] - 1 ;
			    Int rlink = Ai [pparent] ;
			    Match [rlink] = jparent ;
			    match_col [jparent] = rlink ;
			    (void) jchild ;
			}
			nmatch++ ;
			top = -1 ;   /* done with this root */
			advanced = 1 ;
			break ;
		    }
		    else if (!used [jm] && dist [jm] == dist [jc] + 1)
		    {
			dfs_pos [jc] = p + 1 ;
			used [jm] = 1 ;    /* visited this phase */
			stack [++top] = jm ;
			dfs_pos [jm] = Ap [jm] ;
			advanced = 1 ;
			break ;
		    }
		}
		if (!advanced && top >= 0)
		{
		    dfs_pos [jc] = Ap [jc+1] ;
		    used [jc] = 1 ;
		    top-- ;
		}
	    }
	}
    }
    free (stack) ;
    return (nmatch) ;
}

Int TRILINOS_BTF(order)	    /* returns number of blocks found */
(
    /* input, not modified: */
    Int n,	    /* A is n-by-n in compressed column form */
    Int Ap [ ],	    /* size n+1 */
    Int Ai [ ],	    /* size nz = Ap [n] */
    double maxwork, /* do at most maxwork*nnz(A) work in the maximum
		     * transversal; no limit if <= 0 */

    /* output, not defined on input */
    double *work,   /* work performed in maxtrans, or -1 if limit reached */
    Int P [ ],	    /* size n, row permutation */
    Int Q [ ],	    /* size n, column permutation */
    Int R [ ],	    /* size n+1.  block b is in rows/cols R[b] ... R[b+1]-1 */
    Int *nmatch,    /* # nonzeros on diagonal of P*A*Q */

    /* workspace, not defined on input or output */
    Int Work [ ]    /* size 5n */
)
{
    Int *Flag ;
    Int nblocks, i, j, nbadcol ;

    /* ---------------------------------------------------------------------- */
    /* compute the maximum matching */
    /* ---------------------------------------------------------------------- */

    /* Zero-free diagonal quick check: a full structural diagonal is already
     * a maximum matching, and every maximum matching yields the same block
     * triangular form.  The cheap-match phase in maxtrans assigns the first
     * unmatched row of each column, so full-diagonal circuit matrices would
     * otherwise pay an augmenting-path repair for nearly every column. */
    {
	Int diag_full = 1 ;
	for (j = 0 ; diag_full && j < n ; j++)
	{
	    Int p ;
	    Int found = 0 ;
	    for (p = Ap [j] ; p < Ap [j+1] ; p++)
	    {
		if (Ai [p] == j)
		{
		    found = 1 ;
		    break ;
		}
	    }
	    diag_full = found ;
	}
	if (diag_full)
	{
	    for (i = 0 ; i < n ; i++)
	    {
		Q [i] = i ;
	    }
	    *nmatch = n ;
	    if (work != NULL)
	    {
		*work = 0 ;
	    }
	    nblocks = TRILINOS_BTF(strongcomp) (n, Ap, Ai, Q, P, R, Work) ;
	    return (nblocks) ;
	}
    }

    /* if maxwork > 0, then a maximum matching might not be found */

    if (n >= 30000 && getenv ("KLS_NO_HK_MAXTRANS") == NULL)
    {
	Int hk = kls_hk_maxtrans (n, Ap, Ai, Q, Work) ;
	if (hk != TRILINOS_BTF_EMPTY)
	{
	    *nmatch = hk ;
	    if (work != NULL)
	    {
		*work = 0 ;
	    }
	}
	else
	{
	    *nmatch = TRILINOS_BTF(maxtrans) (n, n, Ap, Ai, maxwork, work,
					      Q, Work) ;
	}
    }
    else
    {
	*nmatch = TRILINOS_BTF(maxtrans) (n, n, Ap, Ai, maxwork, work, Q,
					  Work) ;
    }

    /* ---------------------------------------------------------------------- */
    /* complete permutation if the matrix is structurally singular */
    /* ---------------------------------------------------------------------- */

    /* Since the matrix is square, ensure TRILINOS_BTF_UNFLIP(Q[0..n-1]) is a
     * permutation of the columns of A so that A has as many nonzeros on the
     * diagonal as possible.
     */

    if (*nmatch < n)
    {
	/* get a size-n work array */
	Flag = Work + n ;
	for (j = 0 ; j < n ; j++)
	{
	    Flag [j] = 0 ;
	}

	/* flag all matched columns */
	for (i = 0 ; i < n ; i++)
	{
	    j = Q [i] ;
	    if (j != TRILINOS_BTF_EMPTY)
	    {
		/* row i and column j are matched to each other */
		Flag [j] = 1 ;
	    }
	}

	/* make a list of all unmatched columns, in Work [0..nbadcol-1]  */
	nbadcol = 0 ;
	for (j = n-1 ; j >= 0 ; j--)
	{
	    if (!Flag [j])
	    {
		/* j is matched to nobody */
		Work [nbadcol++] = j ;
	    }
	}
	ASSERT (*nmatch + nbadcol == n) ;

	/* make an assignment for each unmatched row */
	for (i = 0 ; i < n ; i++)
	{
	    if (Q [i] == TRILINOS_BTF_EMPTY && nbadcol > 0)
	    {
		/* get an unmatched column j */
		j = Work [--nbadcol] ;
		/* assign j to row i and flag the entry by "flipping" it */
		Q [i] = TRILINOS_BTF_FLIP (j) ;
	    }
	}
    }

    /* The permutation of a square matrix can be recovered as follows: Row i is
     * matched with column j, where j = TRILINOS_BTF_UNFLIP (Q [i]) and where j
     * will always be in the valid range 0 to n-1.  The entry A(i,j) is zero
     * if TRILINOS_BTF_ISFLIPPED (Q [i]) is true, and nonzero otherwise.  nmatch
     * is the number of entries in the Q array that are non-negative. */

    /* ---------------------------------------------------------------------- */
    /* find the strongly connected components */
    /* ---------------------------------------------------------------------- */

    nblocks = TRILINOS_BTF(strongcomp) (n, Ap, Ai, Q, P, R, Work) ;
    return (nblocks) ;
}
