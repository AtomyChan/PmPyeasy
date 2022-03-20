/*****************************************************************************/
/* wcs.h								     */
/*****************************************************************************/

#ifndef	__WCS_H_INCLUDED
#define	__WCS_H_INCLUDED	1

/*****************************************************************************/

#define		M_R2D		(180.0/M_PI)
#define		M_D2R		(M_PI/180.0)

/*****************************************************************************/

#define		WCS_SIN		0
#define		WCS_ARC		1
#define		WCS_TAN		2

typedef	double	vector[3];
typedef double	matrix[3][3];

typedef struct
 {	double	ra0,de0;
	int	type,order;
	matrix	mproj;
	double	zfactor;
 } wcsinit;

typedef struct
 {	wcsinit	init;
 	double	*pixpoly[2];
	double	*prjpoly[2];
	double	crota;
	double	crpix1,crpix2;
	double	cdelt1,cdelt2;
 } wcsdata;

int	wcs_get_projection_matrix(double ra0,double de0,matrix mproj);
int	wcs_get_projected_coords_matrix(matrix mproj,
		double ra,double de,double *rx,double *ry);
int	wcs_invert_projected_coords_matrix(matrix mproj,
		double x,double y,double *rra,double *rde);

int	wcs_project_distort(int type,double *rx,double *ry);
int	wcs_invert_project_distort(int type,double *rx,double *ry);

int	wcs_get_projected_coords(double ra0,double de0,
		double ra,double de,double *rx,double *ry);

/*****************************************************************************/

#endif

/*****************************************************************************/
                         
