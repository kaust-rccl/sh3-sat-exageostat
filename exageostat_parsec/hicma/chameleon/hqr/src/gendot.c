/**
 *
 * @file gendot.c
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

#define DAG_HEADER        "digraph G { orientation=portrait; \n"
#define DAG_FOOTER        "} // close graph\n"
#define DAG_LABELNODE     "%d [label=\"%d\",color=white,pos=\"-1.,-%d.!\"]\n"
#define DAG_LENGTHNODE    "l%d [label=\"%d\",color=white,pos=\"%d.,0.5!\"]\n"
#define DAG_INVISEDGE     "%d->%d [style=\"invis\"];\n"
#define DAG_STARTNODE     "p%d_m%d_k%d [shape=point,width=0.1, pos=\"%d.,-%d.!\",color=\"%s\"];\n"
#define DAG_NODE          "p%d_m%d_k%d [shape=point,width=0.1, pos=\"%d.,-%d.!\",color=\"%s\"];\n"
#define DAG_FIRSTEDGE_PIV "%d->p%d_m%d_k0\n"
#define DAG_FIRSTEDGE_TS  "%d->p%d_m%d_k0 [style=dotted,width=0.1]\n"
#define DAG_FIRSTEDGE_TT  "%d->p%d_m%d_k0 [style=dotted,width=0.1]\n"
#define DAG_EDGE_PIV      "p%d_m%d_k%d->p%d_m%d_k%d [width=0.1,color=\"%s\"]\n"
#define DAG_EDGE_TS       "p%d_m%d_k%d->p%d_m%d_k%d [style=dotted, width=0.1,color=\"%s\"]\n"
#define DAG_EDGE_TT       "p%d_m%d_k%d->p%d_m%d_k%d [style=dashed, width=0.1,color=\"%s\"]\n"

char *color[] = {
    "red",
    "blue",
    "green",
    "orange",
    "cyan",
    "purple",
    "yellow",
};
#define DAG_NBCOLORS 7

void
libhqr_print_dot( const libhqr_tree_t *qrtree,
                  const char          *filename )
{
    int *pos, *next, *done;
    int k, m, n, lpos, prev, length;
    int minMN = libhqr_imin( qrtree->mt, qrtree->nt );
    FILE *f = fopen( filename, "w" );

    done = (int*)malloc( qrtree->mt * sizeof(int) );
    pos  = (int*)malloc( qrtree->mt * sizeof(int) );
    next = (int*)malloc( qrtree->mt * sizeof(int) );
    memset(pos,  0, qrtree->mt * sizeof(int) );
    memset(next, 0, qrtree->mt * sizeof(int) );

    /* Print header */
    fprintf(f, DAG_HEADER ); /*, A->mt+2, minMN+2 );*/
    for(m=0; m < qrtree->mt; m++) {
	fprintf(f, DAG_LABELNODE, m, m, m);
    }

    for(k=0; k<minMN; k++ ) {
	int nb2reduce = qrtree->mt - k - 1;

	for(m=k; m < qrtree->mt; m++) {
	    fprintf(f, DAG_STARTNODE, m, qrtree->mt, k, pos[m], m, color[ (m%qrtree->p) % DAG_NBCOLORS ]);
	    next[m] = qrtree->nextpiv( qrtree, k, m, qrtree->mt);
	}

	while( nb2reduce > 0 ) {
	    memset(done, 0, qrtree->mt * sizeof(int) );
	    for(m=qrtree->mt-1; m > (k-1); m--) {
		n = next[m];
		if ( next[n] != qrtree->mt )
		    continue;
		if ( n != qrtree->mt ) {
		    lpos = libhqr_imax( pos[m], pos[n] );
		    lpos++;
		    pos[m] = lpos;
		    pos[n] = lpos;

		    fprintf(f, DAG_NODE, m, n, k, pos[m], m, color[ (m%qrtree->p) % DAG_NBCOLORS ]);

		    prev = qrtree->prevpiv( qrtree, k, m, n );
		    fprintf(f, DAG_EDGE_PIV,
			    m, prev, k,
			    m, n,    k,
			    color[ (m%qrtree->p) % DAG_NBCOLORS ]);

		    prev = qrtree->prevpiv( qrtree, k, n, n );
		    if ( qrtree->gettype(qrtree, k, n) == 0 )
			fprintf(f, DAG_EDGE_TS,
				n, prev, k,
				m, n, k,
				color[ (m%qrtree->p) % DAG_NBCOLORS ]);
		    else
			fprintf(f, DAG_EDGE_TT,
				n, prev, k,
				m, n, k,
				color[ (m%qrtree->p) % DAG_NBCOLORS ]);

		    next[m] = qrtree->nextpiv( qrtree, k, m, n);
		    done[m] = done[n] = 1;
		    nb2reduce--;
		}
	    }
	}
    }

    length = 0;
    for(m=0; m < qrtree->mt; m++) {
	length = libhqr_imax(length, pos[m]);
    }
    length++;
    for(k=0; k<length; k++)
	fprintf(f, DAG_LENGTHNODE, k, k, k);
    fprintf(f, DAG_FOOTER);

    printf("Tic Max = %d\n", length-1);

    fclose( f );
    free(pos);
    free(next);
}
