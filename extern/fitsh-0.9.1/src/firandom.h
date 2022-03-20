/*****************************************************************************/
/* firandom.h								     */
/*****************************************************************************/

#ifndef	__FIRANDOM_H_INCLUDED
#define	__FIRANDOM_H_INCLUDED	1

#include "magnitude.h"

/*****************************************************************************/

#define		VAR_X		0
#define		VAR_Y		1
#define		VAR_AX		2
#define		VAR_AY		3
#define		VAR_I		4
#define		VAR_M		5
#define		VAR_SH_F	6
#define		VAR_SH_E	7
#define		VAR_SH_P	8
#define		VAR_SN_S	9
#define		VAR_SN_D	10
#define		VAR_SN_K	11
#define		VAR_SI_S	12
#define		VAR_SI_D	13
#define		VAR_SI_K	14

#define		MSK_X		(1<<VAR_X)
#define		MSK_Y		(1<<VAR_Y)
#define		MSK_AX		(1<<VAR_AX)
#define		MSK_AY		(1<<VAR_AY)
#define		MSK_COORD	(MSK_X|MSK_Y|MSK_AX|MSK_AY)
#define		MSK_I		(1<<VAR_I)
#define		MSK_M		(1<<VAR_M)
#define		MSK_INTENSITY	(MSK_I|MSK_M)
#define		MSK_SH_F	(1<<VAR_SH_F)
#define		MSK_SH_E	(1<<VAR_SH_E)
#define		MSK_SH_P	(1<<VAR_SH_P)
#define		MSK_SH		(MSK_SH_F|MSK_SH_E|MSK_SH_P)
#define		MSK_SN_S	(1<<VAR_SN_S)
#define		MSK_SN_D	(1<<VAR_SN_D)
#define		MSK_SN_K	(1<<VAR_SN_K)
#define		MSK_SN		(MSK_SN_S|MSK_SN_D|MSK_SN_K)
#define		MSK_SI_S	(1<<VAR_SI_S)
#define		MSK_SI_D	(1<<VAR_SI_D)
#define		MSK_SI_K	(1<<VAR_SI_K)
#define		MSK_SI		(MSK_SI_S|MSK_SI_D|MSK_SI_K)
#define		MSK_SHAPE	(MSK_SH|MSK_SN|MSK_SI)

#define		NUM_SET		15

#define		PARAM_CNTR	15
#define		RND_1		16
#define		RND_2		17
#define		RND_3		18
#define		RND_4		19
#define		PARAM_SX	20
#define		PARAM_SY	21

#define		NUM_TOTAL	22

/*****************************************************************************/

typedef struct
 {	int	subg;
	double	**subpixeldata;
	double	gain;
	int	is_photnoise,method;
	int	is_intinelect,dontquantize;
	double	nsuppress;
	psf	*tpd;
 } stargenparam;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define		LISTTYPE_FEP	0   /* --fep : x,y,flux,FWHM,ELL,PA	   */
#define		LISTTYPE_SIG	1   /* --sdk : x,y,flux,sigma,delta,kappa  */
#define		LISTTYPE_SDK	2   /* --SDK : x,y,flux,S/D/K momenta      */
#define		LISTTYPE_SMM	3   /* --Smom: x,y,flux,S, M20/M11/M02/... */
#define		LISTTYPE_PSF	4   /* --psf : x,y,flux + external PSF!    */
#define		LISTTYPE_PDS	5   /* --psfdst : x,y,flux + PSF + distort */

typedef struct
 {	magflux	mf0;
	double	ox,oy,scale;
	int	sx,sy,basetype;
 } starlistparam;

typedef struct
 {	int	colx,coly;
	int	colflux,colmag;
	int	colsh1,colsh2,colsh3,colsh4;
	int	listtype;
	magflux	mf0;
 } inlistparam;

/*****************************************************************************/

double	get_gaussian(double mean,double stddev);
int	get_gaussian_2d(double x0,double y0,double is,double id,double ik,double *rx,double *ry);

int	fep_to_sdk(double f,double e,double p,double *s,double *d,double *k);
int	sdk_to_fep(double s,double d,double k,double *f,double *e,double *p);
int	sdk_to_isdk(double s,double d,double k,double *is,double *id,double *ik);
int	isdk_to_sdk(double is,double id,double ik,double *s,double *d,double *k);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int	create_background(fitsimage *img,char *bgarg,double stddev,double ox,double oy,double scale,int zoom);

int	replace_limiters(char *buff);
int	create_input_list(char *buff,starlistparam *lp,star **rstars,int *rnstar,int basistype,int sseed);

/*****************************************************************************/

#endif
