/*****************************************************************************/
/* transform.h								     */
/*****************************************************************************/

#ifndef	__TRANSFORM_H_INCLUDED
#define	__TRANSFORM_H_INCLUDED	1

#include <stdio.h>

#define		TRANS_POLYNOMIAL	1
#define		TRANS_ZOOM		2
#define		TRANS_SHRINK		3

#define		TRANS_WR_COMMENT	0x01
#define		TRANS_WR_DXDY		0x02
#define		TRANS_WR_IEEE_32	0x04
#define		TRANS_WR_IEEE_64	0x08

/*****************************************************************************/

typedef struct 
 {	int	type;
	int	order;
	int	nval;
	double	**vfits;
	double	ox,oy,scale;
	double	bshx,bshy;
 } transformation;

/*****************************************************************************/

int	transformation_free(transformation *tf);
int	transformation_read_data(FILE *fr,transformation *tf);
int	transformation_parse_params(char *params,transformation *tf);
int	transformation_write_data(FILE *fr,transformation *tf,int flags);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int	transformation_check_if_null(transformation *tf);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int	transformation_get_jacobi(transformation *tf,double **rjxx,double **rjxy,double **rjyx,double **rjyy);
int	transformation_eval_normal_2d(double x,double y,transformation *tf,double *rx,double *ry);
int	transformation_eval_invert_2d(double x,double y,transformation *tf,double *rx,double *ry,double *jxx,double *jxy,double *jyx,double *jyy);

/*****************************************************************************/

#endif
