/*****************************************************************************/
/* optcalc.h								     */
/*****************************************************************************/

#ifndef	__OPTCALC_H_INCLUDED	
#define	__OPTCALC_H_INCLUDED	1

/*****************************************************************************/

#include <stdio.h>
#include "math/matrixvector.h"

struct	glass
 {	char	g_name[32];
	double	g_n;
	int	g_nncoeff;
	double	g_ncoeffs[8];
 };

struct	surface
 {	double	s_curvature;
	double	s_conic;
	int	s_nalpha;
	double	s_alphas[8];
 };

struct lens
 {	double	l_offset;
	double	l_thickness;
	double	l_radius1;
	double	l_radius2;
	struct	surface	l_s1;
	struct	surface	l_s2;
	int	l_index_glass;
	double	l_n_lambda;
 };

struct optics
 {	struct	glass	*opt_glasses;
	int		opt_nglass;
	struct	lens	*opt_lenses;
	int		opt_nlens;
	double		opt_z_focal;
	struct	surface	opt_s_focal;
 };

struct raytrace
 {	int	rt_npoint;
	vector	*rt_points;
 };

typedef	double	transfer_matrix[2][2];

/*****************************************************************************/

double	optcalc_get_refraction_index(struct glass *g,double lambda);
int	optcalc_snell_descartes(vector norm,vector s1,vector s2,double n1,double n2);
int	optcalc_reflection(vector norm,vector s1,vector s2);
int	optcalc_surface_aspheric_eval(double c,double k,int nalpha,double *alphas,double r2,double *rz);
int	optcalc_surface_aspheric_diff(double c,double k,int nalpha,double *alphas,double r2,double *rdz);

/* optcalc_surface_aspheric_ray_trace():
   computes the intersection of the aspheric surface described by the 
   parameter set of (c,k,alphas[nalpha]) and the ray going from x0[] to the
   direction n0[]. 0 is returned if there is an intersection and x_int[]
   and n_int[] is set appropriaterly (intersection point and surface normal
   vector, respectively). Here the aspheric surface is directed to tangent the 
   (x,y) plane at (0,0,0) while positive curvature referes to positive z 
   values.
*/
int	optcalc_surface_aspheric_ray_trace(double c,double k,int nalpha,double *alphas,vector x0,vector n0,vector x_int,vector n_int);
int	optcalc_refract_quadratic_ray_trace(matrix A,vector B,double C,vector x0,vector n0,double rn1,double rn2);
int	optcalc_refract_aspheric_ray_trace(double c,double k,int nalpha,double *alphas,double z0,vector x0,vector n0,double rn1,double rn2);

int	optcalc_surface_reset(struct surface *s);
int	optcalc_surface_is_quadratic(struct surface *s);

int	optcalc_lens_ray_trace(struct lens *l,double n_refr,vector x0,vector n0,struct raytrace *rt);

int	optcalc_raytrace_reset(struct raytrace *rt);
int	optcalc_raytrace_free(struct raytrace *rt);

int	optcalc_glass_refraction_precompute(struct optics *opt,double lambda);
int	optcalc_ray_trace(struct optics *opt,double lambda,vector v0,vector n0,struct raytrace *rt);

int	optcalc_compute_transfer_matrix(struct optics *opt,double lambda,transfer_matrix m,double z0,double *rzend);

int	optcalc_reset(struct optics *opt);

int	optcalc_read_optics(FILE *fr,struct optics *opt);

/*****************************************************************************/

#endif

/*****************************************************************************/
