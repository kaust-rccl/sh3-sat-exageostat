/**
 *
 * @file gensvg.c
 *
 * File for algorithm of treewalking.
 *
 * @copyright 2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                 Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Raphael Boucherie
 * @author Mathieu Faverge
 * @date 2017-03-21
 *
 */
#include "libhqr_internal.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#define WIDTH  50
#define HEIGHT 50
#define SIZE   100

/*
 * Global array for color
 */
char *colortree[] = {"red", "blue", "green", "orange", "cyan", "purple", "yellow" };
#define NBCOLORS (sizeof( colortree ) / sizeof( char* ))

/*
 * functions writing in the svg file
 */
static void
drawsvg_header( FILE *file )
{
    int rc;

    rc = fprintf(file,
                 "<?xml version=\"1.0\" standalone=\"no\"?>\n"
                 "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
                 "<svg width=\"2000\" height=\"2000\" version=\"1.1\" \n xmlns=\"http://www.w3.org/2000/svg\">\n");

    if (rc < 0) {
        fprintf(stderr, "Something wrong happened while writing the file (drawsvg_header)\n");
    }
    return;
}

static void
drawsvg_top_TS( FILE *file, int k, int x, int y, int w, int h )
{
    int rc;
    rc = fprintf(file,"<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"%s\" /> \n", x, y, w, h, colortree[k%NBCOLORS]);

    if (rc < 0) {
        fprintf(stderr, "Something wrong happened while writing the file (drawsvg_top_TS)\n");
    }
    return;
}

static void
drawsvg_bot_TS( FILE *file, int k, int x, int y, int w, int h )
{
    int rc, x2, y2, w2, h2;

    rc = fprintf(file,"<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"%s\" /> \n", x, y, w, h, colortree[k%NBCOLORS]);
    if (rc < 0) {
        fprintf(stderr, "Something wrong happened while writing the file (drawsvg_bot_TS)\n");
        return;
    }

    x2 = x + (w / 4);
    y2 = y + (h / 4);
    w2 = (w / 2);
    h2 = (h / 2);

    rc = fprintf(file,"<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill =\"white\"/> \n", x2, y2, w2, h2);
    if (rc < 0) {
        fprintf(stderr, "Something wrong happened while writing the file (drawsvg_bot_TS)\n");
    }
    return;
}

static void
drawsvg_top_TT( FILE *file, int k, int x, int y, int w, int h )
{
    int rc;
    rc = fprintf( file,"<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" /> \n",
                  x + w / 2, y + h / 2, w / 2, colortree[k%NBCOLORS] );
    if (rc < 0) {
        fprintf(stderr, "Something wrong happened while writing the file (drawsvg_top_TT)\n");
    }
    return;
}

static void
drawsvg_bot_TT( FILE *file, int k, int x, int y, int w, int h )
{
    int rc;
    rc = fprintf( file,"<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"%s\" /> \n",
                  x + w / 2, y + h / 2, w / 2, colortree[k%NBCOLORS] );
    if (rc < 0) {
        fprintf(stderr, "Something wrong happened while writing the file (drawsvg_bot_TT)\n");
        return;
    }

    rc = fprintf( file,"<circle cx=\"%d\" cy=\"%d\" r=\"%d\" fill=\"white\" /> \n",
                  x + w / 2, y + h / 2, w / 4 );
    if (rc < 0) {
        fprintf(stderr, "Something wrong happened while writing the file (drawsvg_bot_TT)\n");
    }
    return;
}

static void
drawsvg_line( FILE *file, int k, int x1, int y1, int x2, int y2 )
{
    int rc;
    rc = fprintf(file,"<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" style=\"fill:none;stroke:%s;stroke-width:2px;\"/> \n", x1, y1, x2, y2, colortree[k%NBCOLORS]);

    if (rc < 0) {
        fprintf(stderr, "Something wrong happened while writing the file (drawsvg_line)\n");
    }
    return;
}

static void
drawsvg_footer( FILE *file )
{
    int rc;
    rc = fprintf(file, "</svg>\n");
    if (rc < 0) {
        fprintf(stderr, "Something wrong happened while writing the file (drawsvg_footer)\n");
    }
    return;
}


static void
drawsvg_lines_rowm( FILE *file, int k,
                    int p, int m, int beg_p, int beg_m, int end )
{
    int yp, ym;
    int x, xp, xm;

    /* Row of the tiles */
    ym = SIZE + SIZE * m;
    yp = SIZE + SIZE * p;

    /* Starting position of the tiles */
    xm = (SIZE + (SIZE / 4)) + SIZE * beg_m;
    xp = (SIZE + (SIZE / 4)) + SIZE * beg_p;

    /* Final position of the tiles */
    x = SIZE + SIZE * end;

    /* Horizontal lines */
    drawsvg_line( file, k, xm, ym, x + (SIZE / 4), ym );
    drawsvg_line( file, k, xp, yp, x + (SIZE / 4), yp );

    /* Vertical line */
    drawsvg_line( file, k, x, ym, x, yp );
}

static void
drawsvg_lines_stepk( const libhqr_tree_t *qrtree, FILE *file,
                     int k, int *tiles, int *step )
{
    int i, m, p, end;

    /* Get order for step k */
    libhqr_walk_stepk( qrtree, k, tiles+(k+1) );

    for(i = k+1; i < qrtree->mt; i++){
        m = tiles[i];
        p = qrtree->currpiv(qrtree, k, m);

        end = libhqr_imax( step[p], step[m] ) + 1;

        /* Draw horizontal lines for rows p and m */
        drawsvg_lines_rowm( file, k, p, m, step[p], step[m], end );

        /* Update last time the rows p and m have been modified for future lines */
        step[m] = end;
        step[p] = end;
    }
}

static void
drawsvg_nodes_rowm( FILE *file, int k,
                    int type, int p, int m, int step_m )
{
    int x, yp, ym;
    x  = ((SIZE * 3) / 4) + SIZE * step_m;
    ym = ((SIZE * 3) / 4) + SIZE * m;
    yp = ((SIZE * 3) / 4) + SIZE * p;

    if ( type == LIBHQR_KILLED_BY_TS ) {
        drawsvg_top_TS(file, k, x, yp, WIDTH, HEIGHT );
        drawsvg_bot_TS(file, k, x, ym, WIDTH, HEIGHT );
    }
    else {
        drawsvg_top_TT(file, k, x, yp, WIDTH, HEIGHT );
        drawsvg_bot_TT(file, k, x, ym, WIDTH, HEIGHT );
    }
}

void
libhqr_print_svg( const libhqr_tree_t *qrtree,
                  const char          *filename )
{
    FILE *file;
    int  *tiles, *steps;
    int   minMN, k, i;

    minMN = libhqr_imin( qrtree->mt, qrtree->nt );

    tiles = (int*)calloc( qrtree->mt, sizeof(int));
    steps = (int*)calloc( qrtree->mt, sizeof(int));

    file = fopen(filename,"w+");

    drawsvg_header(file);
    for (k = 0; k < minMN; k++) {
        /* Drawing the lines */
        drawsvg_lines_stepk( qrtree, file, k, tiles, steps );

        /* Drawing the rectangles */
        for (i = k+1; i < qrtree->mt; i++) {
            drawsvg_nodes_rowm( file, k,
                                qrtree->gettype(qrtree, k, i),
                                qrtree->currpiv(qrtree, k, i), i, steps[i] );
        }
    }
    drawsvg_footer(file);
    fclose(file);

    free(tiles);
    free(steps);
}
