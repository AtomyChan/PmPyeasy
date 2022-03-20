/*****************************************************************************/
/* background.h								     */
/*****************************************************************************/

#ifndef	__BACKGROUND_H_DEFINED
#define	__BACKGROUND_H_DEFINED	1

/*****************************************************************************/

#define		BACKGROUND_POLYNOMIAL		1
#define		BACKGROUND_SPLINE		2

/*****************************************************************************/

typedef struct
 {	int	type,blocks,order;
	double	sx,sy,scale;
	double	*xfit,*yfit;
 } background;

typedef struct
 {      int     order;
        double  ox,oy,scale;
        double  *coeff;
 } spatial;

int	determine_background(fitsimage *img,spatial *bg,int gx,int gy,int order);

/*****************************************************************************/

#endif

/*****************************************************************************/
                                                                
