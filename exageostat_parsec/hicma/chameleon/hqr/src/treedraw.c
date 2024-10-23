/**
 *
 * @file treedraw.c
 *
 * All the functions required for drawing tree are here.
 *
 * @copyright 2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                 Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Raphael Boucherie
 * @author Mathieu Faverge
 * @date 2017-04-04
 *
 */
#include <string.h>
#include <stdio.h>
#include "libhqr_draw.h"

/*
 * Global array for color
 */

char *colortree[4] = {"red", "blue", "green", "purple"};

/*
 * functions writing in the svg file
 */

void libhqr_writeheader(FILE *file){
    if(fprintf(file, "<?xml version=\"1.0\" standalone=\"no\"?>\n"
                     "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
                     "<svg width=\"2000\" height=\"2000\" version=\"1.1\" \n xmlns=\"http://www.w3.org/2000/svg\">\n") <0)
        return;
}

/*
 * Common prameters to the 2 following functions:
 *    x - Parameter x for the x-axis
 *    y - Parameter y for the y-axis
 *    w - Parameter w for the width
 *    h - Parameter h for the height
 *    k - Factorization step for the color
 */

void libhqr_drawTT(int x, int y, int w, int h, int k, FILE *file){
    if(fprintf(file,"<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"%s\" /> \n", x, y, w, h, colortree[k%4]) < 0 )
        return;
}

void libhqr_drawTS(int x, int y, int w, int h, int k, FILE *file){
    if(fprintf(file,"<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"%s\" /> \n", x, y, w, h, colortree[k%4]) < 0 )
        return;
    int x2 = x + (w / 4);
    int y2 = y + (h / 4);
    int w2 = (w / 2);
    int h2 = (h / 2);
    if(fprintf(file,"<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill =\"white\"/> \n", x2, y2, w2, h2) < 0 ) return;
}

void libhqr_drawline(int x1, int y1, int x2, int y2, int k, FILE *file){
    if(fprintf(file,"<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" style=\"fill:none;stroke:%s;stroke-width:2px;\"/> \n", x1, y1, x2, y2, colortree[k%4]) < 0 ) return;
}

void libhqr_writeend(FILE *file){
    if(fprintf(file, "</svg>") < 0) return;
}
