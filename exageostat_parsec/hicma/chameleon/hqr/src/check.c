/**
 *
 * @file check.c
 *
 * @copyright 2010-2017 The University of Tennessee and The University
 *                      of Tennessee Research Foundation.  All rights
 *                      reserved.
 *
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Raphael Boucherie
 * @author Mathieu Faverge
 * @date 2017-03-21
 *
 */
#include "libhqr_internal.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* static int libhqr_qrtree_getinon0( const qr_piv_t *arg,  */
/*                                 const int k, int i, int mt ); */

#define ENDCHECK( test, ret )                   \
    if ( !test ) {                              \
        assert( ret == 0 );                     \
        return ret;                             \
    }

int
libhqr_check( const libhqr_matrix_t *A, const libhqr_tree_t *qrtree)
{
    int minMN = libhqr_imin(A->mt, A->nt );
    int i, m, k, nb;
    int check;

    int a = qrtree->a;
    int p = qrtree->p;

    /*
     * Check Formula for NB geqrt
     */
    {
	/* libhqr_print_type( A, qrtree ); */
	/* libhqr_print_nbgeqrt( A, qrtree ); */
	check = 1;
	for (k=0; k<minMN; k++) {
	    nb = 0;
	    for (m=k; m < A->mt; m++) {
		if ( qrtree->gettype( qrtree, k, m ) > 0 )
		    nb++;
	    }

	    if ( nb != qrtree->getnbgeqrf( qrtree, k ) ) {
		check = 0;
		printf(" ----------------------------------------------------\n"
		       "  - a = %d, p = %d, M = %d, N = %d\n"
		       "     Check number of geqrt:\n"
		       "       For k=%d => return %d instead of %d",
		       a, p, A->mt, A->nt, k, qrtree->getnbgeqrf( qrtree, k ), nb );
	    }
	}

	ENDCHECK( check, 1 );
    }

    /*
     * Check indices of geqrt
     */
    {
	int prevm = -1;
	check = 1;
	for (k=0; k<minMN; k++) {
	    /* libhqr_print_geqrt_k( A, qrtree, k ); */
	    nb = qrtree->getnbgeqrf( qrtree, k );
	    prevm = -1;
	    for (i=0; i < nb; i++) {

		m = qrtree->getm( qrtree, k, i );

		/*
		 * getm has to be the inverse of geti
		 */
		if ( i != qrtree->geti( qrtree, k, m) ) {
		    check = 0;
		    printf(" ----------------------------------------------------\n"
			   "  - a = %d, p = %d, M = %d, N = %d\n"
			   "     Check indices of geqrt:\n"
			   "        getm( k=%d, i=%d ) => m = %d && geti( k=%d, m=%d ) => i = %d\n",
			   a, p, A->mt, A->nt,
			   k, i, m, k, m, qrtree->geti( qrtree, k, m));
		}
		/* tile before the diagonal are factorized and
		 * the m is a growing list (not true with round-robin inside TS)
		 */
		else if ( (a == 1) && (( m < k ) || ( m < prevm )) ) {
		    check = 0;
		    printf(" ----------------------------------------------------\n"
			   "  - a = %d, p = %d, M = %d, N = %d\n"
			   "     Check indices of geqrt:\n"
			   "        getm( k=%d, i=%d ) => m = %d",
			   a, p, A->mt, A->nt, k, i, m);
		}
#if 0
		else if ( m != qrtree->getinon0( qrtree, k, i, A->mt ) ) {
		    check = 0;
		    printf(" ----------------------------------------------------\n"
			   "  - a = %d, p = %d, M = %d, N = %d\n"
			   "     Check indices of geqrt:\n"
			   "        getm( k=%d, i=%d ) => m = %d but should be %d",
			   a, p, A->mt, A->nt, k, i, m, qrtree->getinon0( qrtree, k, i, A->mt));
		}
#endif
		prevm = m;
	    }
	}
	ENDCHECK( check, 2 );
    }

    /*
     * Check number of exit in next
     */
    {
	int s;
	check = 1;

	for (k=0; k<minMN; k++) {
	    for(m=k; m<A->mt; m++) {
		nb = 0;
		for(s=A->mt; s>k; s--) {
		    if ( qrtree->nextpiv(qrtree, k, m, s) == A->mt )
			nb++;
		}
		if ( nb > 1 ) {
		    libhqr_print_next_k( A, qrtree, k);
		    libhqr_print_prev_k( A, qrtree, k);

		    printf(" ----------------------------------------------------\n"
			   "  - a = %d, p = %d, M = %d, N = %d\n"
			   "     Next of line %d for step %d contains more than one exit:\n",
			   a, p, A->mt, A->nt,
			   m, k);
		    check = 0;
		    return 3;
		}
		else if ( nb == 0 ) {
		    libhqr_print_next_k( A, qrtree, k);
		    libhqr_print_prev_k( A, qrtree, k);

		    printf(" ----------------------------------------------------\n"
			   "  - a = %d, p = %d, M = %d, N = %d\n"
			   "     Next of line %d for step %d needs one exit:\n",
			   a, p, A->mt, A->nt,
			   m, k);
		    check = 0;
                    ENDCHECK( check, 3 );
		}
	    }
	}
	ENDCHECK( check, 3 );
    }

    /*
     * Check number of exit in prev
     */
    {
	int s;
	check = 1;

	for (k=0; k<minMN; k++) {
	    for(m=k; m<A->mt; m++) {
		nb = 0;
		for(s=k; s<A->mt; s++) {
		    if ( qrtree->prevpiv(qrtree, k, m, s) == A->mt )
			nb++;
		}
		if ( nb > 1 ) {
		    libhqr_print_next_k( A, qrtree, k);
		    libhqr_print_prev_k( A, qrtree, k);

		    printf(" ----------------------------------------------------\n"
			   "  - a = %d, p = %d, M = %d, N = %d\n"
			   "     Prev of line %d for step %d contains more than one exit:\n",
			   a, p, A->mt, A->nt,
			   m, k);
		    check = 0;
		    return 3;
		}
		else if ( nb == 0 ) {
		    libhqr_print_next_k( A, qrtree, k);
		    libhqr_print_prev_k( A, qrtree, k);

		    printf(" ----------------------------------------------------\n"
			   "  - a = %d, p = %d, M = %d, N = %d\n"
			   "     Prev of line %d for step %d needs one exit:\n",
			   a, p, A->mt, A->nt,
			   m, k);
		    check = 0;
                    ENDCHECK( check, 3 );
		}
	    }
	}
	ENDCHECK( check, 3 );
    }

    /*
     * Check next/prev
     */
    {
	int start, next, prev;
	check = 1;

	for (k=0; k<minMN; k++) {
	    start = A->mt;
	    for(m=k; m<A->mt; m++) {

		do {
		    next = qrtree->nextpiv(qrtree, k, m, start);
		    if ( next == A->mt )
			prev = qrtree->prevpiv(qrtree, k, m, m);
		    else
			prev = qrtree->prevpiv(qrtree, k, m, next);

		    if ( start != prev ) {
			libhqr_print_next_k( A, qrtree, k);
			libhqr_print_prev_k( A, qrtree, k);

			printf(" ----------------------------------------------------\n"
			       "  - a = %d, p = %d, M = %d, N = %d\n"
			       "     Check next/prev:\n"
			       "       next( k=%d, m=%d, start=%d ) => %d && prev( k=%d, m=%d, start=%d ) => %d (instead of %d)\n",
			       a, p, A->mt, A->nt,
			       k, m, start, next, k, m, next, prev, start);
                        check = 0;
			ENDCHECK( check, 3 );
		    }
		    start = next;
		} while ( start != A->mt );
	    }
	}
	ENDCHECK( check, 3 );
    }

    return 0;
}

void
libhqr_print_type( const libhqr_matrix_t *A, const libhqr_tree_t *qrtree )
{
    int minMN = libhqr_imin(A->mt, A->nt );
    int m, k;
    int lm = 0;
    int lmg = 0;
    int rank = 0;

    printf("\n------------ Localization = Type of pivot --------------\n");
    for(m=0; m<A->mt; m++) {
	printf("%3d | ", m);
	for (k=0; k<libhqr_imin(minMN, m+1); k++) {
	    printf( "%3d ", qrtree->gettype( qrtree, k, m ) );
	}
	for (k=libhqr_imin(minMN, m+1); k<minMN; k++) {
	    printf( "    " );
	}

	printf("    ");
	printf("%2d,%3d | ", rank, lmg);
	for (k=0; k<libhqr_imin(minMN, lmg+1); k++) {
	    printf( "%3d ", qrtree->gettype( qrtree, k, lmg) );
	}
	for (k=libhqr_imin(minMN, lmg+1); k<minMN; k++) {
	    printf( "    " );
	}
	lm++; lmg+=qrtree->p;
	if ( lmg >= A->mt ) {
	    rank++;
	    lmg = rank;
	    lm = 0;
	}
	printf("\n");
    }
}

void
libhqr_print_pivot( const libhqr_matrix_t *A, const libhqr_tree_t *qrtree )
{
    int minMN = libhqr_imin(A->mt, A->nt );
    int m, k;
    int lm = 0;
    int lmg = 0;
    int rank = 0;
    printf("\n------------ Current Pivot--------------\n");
    for(m=0; m<A->mt; m++) {
	printf("%3d | ", m);
	for (k=0; k<libhqr_imin(minMN, m+1); k++) {
	    printf( "%3d ", qrtree->currpiv(qrtree, k, m) );
	}
	for (k=libhqr_imin(minMN, m+1); k<minMN; k++) {
	    printf( "    " );
	}

	printf("    ");
	printf("%2d,%3d | ", rank, lmg);
	for (k=0; k<libhqr_imin(minMN, lmg+1); k++) {
	    printf( "%3d ", qrtree->currpiv(qrtree, k, lmg) );
	}
	for (k=libhqr_imin(minMN, lmg+1); k<minMN; k++) {
	    printf( "    " );
	}
	lm++; lmg+=qrtree->p;
	if ( lmg >= A->mt ) {
	    rank++;
	    lmg = rank;
	    lm = 0;
	}
	printf("\n");
    }
}

void
libhqr_print_next_k( const libhqr_matrix_t *A, const libhqr_tree_t *qrtree, int k )
{
    int m, s;
    printf("\n------------ Next (k = %d)--------------\n", k);

    printf( "      " );
    for(s=A->mt; s>0; s--)
	printf( "%3d ", s );
    printf( "\n" );

    for(m=0; m<A->mt; m++) {
	printf("%3d | ", m);
	for(s=A->mt; s>0; s--) {
	    printf( "%3d ", qrtree->nextpiv(qrtree, k, m, s) );
	}
	printf("\n");
    }
}

void
libhqr_print_prev_k( const libhqr_matrix_t *A, const libhqr_tree_t *qrtree, int k )
{
    int m, s;
    printf("\n------------ Prev (k = %d)--------------\n", k);

    printf( "      " );
    for(s=A->mt; s>-1; s--)
	printf( "%3d ", s );
    printf( "\n" );

    for(m=0; m<A->mt; m++) {
	printf("%3d | ", m);
	for(s=A->mt; s>-1; s--) {
	    printf( "%3d ", qrtree->prevpiv(qrtree, k, m, s) );
	}
	printf("\n");
    }
}

void
libhqr_print_perm( const libhqr_matrix_t *A, const libhqr_tree_t *qrtree, int *perm )
{
    int minMN = libhqr_imin(A->mt, A->nt );
    int m, k;
    (void)qrtree;

    printf("\n------------ Permutation --------------\n");
    for (k=0; k<minMN; k++) {
	printf( "%3d ", k );
    }
    printf( "\n" );
    for (k=0; k<minMN; k++) {
	printf( "----" );
    }
    printf( "\n" );

    for (m=0; m < A->mt+1; m++) {
	for (k=0; k<minMN; k++) {
	    printf( "%3d ", perm[ k*(A->mt+1) + m ] );
	}
	printf( "\n" );
    }
    printf( "\n" );
}

void
libhqr_print_nbgeqrt( const libhqr_matrix_t *A, const libhqr_tree_t *qrtree )
{
    int minMN = libhqr_imin(A->mt, A->nt );
    int m, k, nb;

    printf("\n------------ Nb GEQRT per k --------------\n");
    printf(" k      : ");
    for (k=0; k<minMN; k++) {
	printf( "%3d ", k );
    }
    printf( "\n" );
    printf(" Compute: ");
    for (k=0; k<minMN; k++) {
	nb = 0;
	for (m=k; m < A->mt; m++) {
	    if ( qrtree->gettype(qrtree, k, m) > 0 )
		nb++;
	}
	printf( "%3d ", nb );
    }
    printf( "\n" );
    printf(" Formula: ");
    for (k=0; k<minMN; k++) {
	printf( "%3d ", qrtree->getnbgeqrf( qrtree, k ) );
    }
    printf( "\n" );
}

void
libhqr_print_geqrt_k( const libhqr_matrix_t *A, const libhqr_tree_t *qrtree, int k )
{
    int i, m, nb;
    (void)A;

    printf("\n------------ Liste of geqrt for k = %d --------------\n", k);

    printf( "  m:");
    nb = qrtree->getnbgeqrf( qrtree, k );
    for (i=0; i < nb; i++) {
	m = qrtree->getm( qrtree, k, i );
	if ( i == qrtree->geti( qrtree, k, m) )
	    printf( "%3d ", m );
	else
	    printf( "x%2d ", qrtree->geti( qrtree, k, m) );
    }
    printf( "\n" );
}
