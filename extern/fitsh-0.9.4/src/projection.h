/*****************************************************************************/
/* projection.h								     */
/*****************************************************************************/

#ifndef	__PROJECTION_H_INCLUDED
#define	__PROJECTION_H_INCLUDED	1

/*****************************************************************************/

#define		M_R2D		(180.0/M_PI)
#define		M_D2R		(M_PI/180.0)

/*****************************************************************************/

#define		PROJECTION_SIN		0
#define		PROJECTION_ARC		1
#define		PROJECTION_TAN		2

typedef	double	vector[3];
typedef double	matrix[3][3];

#define		PROJECTION_BCR_NCOEFF	4
#define		PROJECTION_BCT_NCOEFF	4

typedef struct
 {	int	brown_conrady_radial_ncoeff;
	double	brown_conrady_radial_coeffs[PROJECTION_BCR_NCOEFF];
	int	brown_conrady_tangential_ncoeff;
	double	brown_conrady_tangential_coeffs[PROJECTION_BCT_NCOEFF];
 } projection_distort;

int	projection_get_matrix_rdr(double ra0,double de0,double roll0,matrix mproj);
int	projection_get_matrix_quat(double qr,double qi,double qj,double qk,matrix mproj);
int	projection_do_matrix_coord(matrix mproj,
		double ra,double de,double *rx,double *ry,double *rz);
int	projection_do_inverse_matrix_coord(matrix mproj,
		double x,double y,double *rra,double *rde);

int	projection_do_distortion(int type,projection_distort *dist,double *rx,double *ry,double *rz);
int	projection_do_inverse_distortion(int type,projection_distort *dist,double *rx,double *ry,double *rz);

int	projection_do_coord(double ra0,double de0,double roll0,
		double ra,double de,double *rx,double *ry,double *rz);

/*****************************************************************************/

#endif

/*****************************************************************************/
                    
