/*****************************************************************************/
/* weight.h								     */
/*****************************************************************************/

#ifndef	__WEIGHT_H_INCLUDED
#define	__WEIGHT_H_INCLUDED	1

#include <fits/fits.h>

#include "stars.h"
#include "psf.h"

/*****************************************************************************/

typedef struct
 {	double	x,y;
	int	ix,iy;
	double	flux;
	double	**iarr;
 } weight;

typedef struct
 {	double	***zdata;
	weight	*weights;
	int	nweight;
	int	hsize,grid;
 } weightlist;

/* weight-star.c */ /*********************************************************/

int	weight_draw(weightlist *wl,star *stars,int nstar,
		int hsize,int grid,psf *tpd);

/* weight-io.c */ /***********************************************************/

fits *	weight_fits_create(weightlist *wl);
int	weight_parse_fits(fits *img,weightlist *wl);

/* weight-gen.c */ /**********************************************************/

int	weight_sort(weightlist *wl);
weight *weight_get_closest(weightlist *wl,double x0,double y0);
int	weight_get_intersec_list(weightlist *wl,int x0,int y0,int sx,int sy,
		int *list,int maxcnt);

/*****************************************************************************/

#endif

/*****************************************************************************/
