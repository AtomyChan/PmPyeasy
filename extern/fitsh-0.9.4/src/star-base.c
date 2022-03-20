/*****************************************************************************/
/* star-base.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Base module (common functions) for star searching and fitting...	     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "fitsh.h"
#include "common.h"
#include "stars.h"

/*****************************************************************************/

/* fit_small_parabola_point(), fit_small_parabola_block(),
   Fits the A+B*(dx)+C*(dy)+1/2*D*(dx)^2+E*(dx)*(dy)+1/2*F*(dy)^2 parabola 
   to the nine points / blocks of the image section [x-1:x+1,y-1:y+1] 
   (dx=-1, 0, 1 and dy=-1, 0, 1 for this nine point/block). The fit parameters 
   A, B, C, D, E and F are stored in the array fit: in fit[0], fit[1], fit[2], 
   fit[3], fit[4], and fit[5] respectively. If the point (x,y) is out of 
   the image boundary, a nonzero value is returned indicating the error, 
   otherwise the functions return 0.					     */

int fit_small_parabola_point(fitsimage *img,int x,int y,double *fit)
{
 double	bv[6],xls,xhs,yls,yhs;

 if ( x<1 || y<1 || x>=img->sx-1 || y>=img->sy-1 )	return(1);

 xls=img->data[y-1][x-1]+img->data[y  ][x-1]+img->data[y+1][x-1];
 xhs=img->data[y-1][x+1]+img->data[y  ][x+1]+img->data[y+1][x+1];

 yls=img->data[y-1][x-1]+img->data[y-1][x  ]+img->data[y-1][x+1];
 yhs=img->data[y+1][x-1]+img->data[y+1][x  ]+img->data[y+1][x+1];

 bv[0]=xls+(img->data[y-1][x]+img->data[y][x]+img->data[y+1][x])+xhs;
 bv[1]=xhs-xls;
 bv[2]=yhs-yls;
 bv[3]=0.5*(xhs+xls);
 bv[4]=	(img->data[y-1][x-1]+img->data[y+1][x+1])-
	(img->data[y+1][x-1]+img->data[y-1][x+1]);
 bv[5]=0.5*(yhs+yls);

 fit[0]=(5.0*bv[0]-6.0*(bv[3]+bv[5]))/9.0;
 fit[1]=bv[1]/6.0;
 fit[2]=bv[2]/6.0;
 fit[3]=(-2.0*bv[0]+6.0*bv[3])/3.0;
 fit[4]=0.25*bv[4];
 fit[5]=(-2.0*bv[0]+6.0*bv[5])/3.0;
 return(0);
}

int fit_small_parabola_block(fitsimage *img,int x,int y,double *fit)
{
 double	bv[6],xls,xms,xhs,yls,yms,yhs;

 if ( x<1 || y<1 || x>=img->sx-1 || y>=img->sy-1 )	return(1);

 xls=img->data[y-1][x-1]+img->data[y  ][x-1]+img->data[y+1][x-1];
 xms=img->data[y-1][x  ]+img->data[y  ][x  ]+img->data[y+1][x  ];
 xhs=img->data[y-1][x+1]+img->data[y  ][x+1]+img->data[y+1][x+1];

 yls=img->data[y-1][x-1]+img->data[y-1][x  ]+img->data[y-1][x+1];
 yms=img->data[y  ][x-1]+img->data[y  ][x  ]+img->data[y  ][x+1];
 yhs=img->data[y+1][x-1]+img->data[y+1][x  ]+img->data[y+1][x+1];

 bv[0]=xls+xms+xhs;
 bv[1]=xhs-xls;
 bv[2]=yhs-yls;
 bv[3]=(13.0*xhs+xms+13.0*xls)/24.0;
 bv[4]=	(img->data[y-1][x-1]+img->data[y+1][x+1])-
	(img->data[y+1][x-1]+img->data[y-1][x+1]);
 bv[5]=(13.0*yhs+yms+13.0*yls)/24.0;

 fit[0]=(97.0/144.0)*bv[0]-0.75*(bv[3]+bv[5]);
 fit[1]=(1.0/6.0)*bv[1];
 fit[2]=(1.0/6.0)*bv[2];
 fit[3]=-0.75*bv[0]+2*bv[3];
 fit[4]= 0.25*bv[4];
 fit[5]=-0.75*bv[0]+2*bv[5];
 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int fit_small_parabola_block_param(fitsimage *img,int j,int i,
double *rcx,double *rcy,double *raxx,double *raxy,double *rayy,double *rpeak)
{
 double	pfit[6],a,ax,ay,axx,axy,ayy,cx,cy,det,tr,peak;

 if ( fit_small_parabola_block(img,j,i,pfit) )	return(1);

 a=pfit[0],ax=pfit[1],ay=pfit[2],
 axx=pfit[3],axy=pfit[4],ayy=pfit[5];
 det=axx*ayy-axy*axy; tr=axx+ayy;
 if ( det<=0.0 || tr>=0.0 )		return(1);
 cx=-(+ayy*ax-axy*ay)/det,
 cy=-(-axy*ax+axx*ay)/det;
 if ( fabs(cx)>1 || fabs(cy)>1 )	return(1);

 peak=a+ax*cx+ay*cy+0.5*(axx*cx*cx+2.0*axy*cx*cy+ayy*cy*cy);

 *rcx=cx;
 *rcy=cy;
 *raxx=axx,*raxy=axy,*rayy=ayy;
 *rpeak=peak;

 return(0);
}

/*****************************************************************************/

static int order_candidates_by_peak_compare(const void *vc1,const void *vc2)
{
 candidate *c1,*c2;
 c1=(candidate *)vc1,c2=(candidate *)vc2;
      if ( c1->peak > c2->peak )	return(1);
 else if ( c1->peak < c2->peak )	return(-1);
 else					return(0);
}
int order_candidates_by_peak(candidate *cands,int ncand)
{
 qsort(cands,ncand,sizeof(candidate),order_candidates_by_peak_compare);
 return(0);
}

/*****************************************************************************/

int cleanup_starlist(star **rstars,int *rnstar)
{
 star	*stars,*tstars;
 int	nstar,i,nn;

 stars=*rstars,nstar=*rnstar;
 
 tstars=(star *)malloc(sizeof(star)*nstar);
 memcpy(tstars,stars,sizeof(star)*nstar);

 nn=0;
 for ( i=0 ; i<nstar ; i++ )
  {	if ( ! tstars[i].marked )
	 {	memcpy(&stars[nn],&tstars[i],sizeof(star));
		nn++;
	 }
  }
 free(tstars);

 nstar=nn;
 stars=realloc(stars,sizeof(star)*nstar);

 *rstars=stars,*rnstar=nstar;
 return(0);
}

int cleanup_candlist(candidate **rcands,int *rncand)
{
 candidate	*cands;
 int		ncand,i;

 cands=*rcands,ncand=*rncand;

 for ( i=0 ; i<ncand ; )
  {	if ( cands[i].marked )
	 {	if ( cands[i].nipoint > 0 && cands[i].ipoints != NULL )
			free(cands[i].ipoints);
		if ( i<ncand-1 )
			memmove(cands+i,cands+i+1,sizeof(candidate)*(ncand-i-1));
		ncand--;
	 }
	else	i++;
  }
 cands=realloc(cands,sizeof(candidate)*ncand);

 *rcands=cands,*rncand=ncand;
 return(0);
}

/*****************************************************************************/

int free_stars(star *stars,int nstar)
{
 free(stars);
 return(0);
}

int free_candidates(candidate *cands,int ncand)
{
 int	i;
 for ( i=0 ; i<ncand ; i++ )
  {	if ( cands[i].ipoints != NULL && cands[i].nipoint>0 )
		free(cands[i].ipoints);
  }
 free(cands);
 return(0);
}

/*****************************************************************************/

#define		MAX_FIT		(MAX_DEVIATION_COEFF+3)

/* F=B+A1*I1(x1,y1)+A2*I2(x2,y2)+...					     */
/* a[] = { B, x1, y1, A1, <shape1>, x2, y2, A2, <shape2>, ... }		     */
/* xpnt = { ix, iy }							     */
void model_merge(void *xpnt,double *a,double *yy,double *dyda,void *p)
{
 modelparam	*mp=p;
 double		wa[MAX_FIT],wdyda[MAX_FIT],wy;
 int		i;

 *yy=a[0];
 dyda[0]=1.0;
 a++,dyda++;

 while ( mp->funct != NULL && mp->nshape>0 )
  {	wa[0]=a[2];
	wa[1]=0.0;
	wa[2]=a[0];
	wa[3]=a[1];
	a+=3;
	for ( i=0 ; i<mp->nshape ; i++,a++ )
	 {	wa[i+4]=*a;			}
	mp->funct(xpnt,wa,&wy,wdyda,mp->param);
	*yy+=wy;
	dyda[0]=wdyda[2];
	dyda[1]=wdyda[3];
	dyda[2]=wdyda[0];
	dyda+=3;
	for ( i=0 ; i<mp->nshape ; i++,dyda++ )
	 {	*dyda=wdyda[i+4];		}
	mp++;
  }
}

/* F=B+A1*I1(x1,y1)+A2*I2(x2,y2)+...					     */
/* a[] = { B, x0, y0, A1, <shape1>, A2, <shape2>, A3, <shape3>, ... }	     */
/* xpnt = { ix, iy }							     */
void model_combine(void *xpnt,double *a,double *yy,double *dyda,void *p)
{
 modelparam	*mp=p;
 double		wa[MAX_FIT],wdyda[MAX_FIT],wy,*dyda0,*a0,amp;
 int		i;

 dyda0=dyda,a0=a;

 *yy=a[0];
 a+=3;

 if ( dyda != NULL )
  {	dyda0[0]=1.0;
	dyda0[1]=0.0;
	dyda0[2]=0.0;
	dyda+=3;
  }

 while ( mp->funct != NULL && mp->nshape>0 )
  {	amp  =a[0],a++;
	wa[0]=amp;
	wa[1]=0.0;
	wa[2]=a0[1];
	wa[3]=a0[2];

	for ( i=0 ; i<mp->nshape ; i++,a++ )
	 {	wa[i+4]=*a;			}

	if ( dyda != NULL )	mp->funct(xpnt,wa,&wy,wdyda,mp->param);
	else			mp->funct(xpnt,wa,&wy,NULL ,mp->param);

	*yy+=wy;

	if ( dyda != NULL )
	 {	dyda0[1]+=wdyda[2];
		dyda0[2]+=wdyda[3];
		dyda [0] =wdyda[0];
		dyda++;
		for ( i=0 ; i<mp->nshape ; i++,dyda++ )
		 {	*dyda=wdyda[i+4];		}
	 }

	mp++;
  }
}

/*****************************************************************************/

int refine_candidate_data(fitsimage *img,candidate *wc)
{
 ipoint	*wi;
 double	cx,cy,cxx,cxy,cyy,m,mx,my,mxx,mxy,myy,w,x,y,det,x0,y0;
 int	i;

 if ( wc->ipoints==NULL || wc->nipoint<=0 )	return(1);

 x0=wc->cx,
 y0=wc->cy;

 m=mx=my=0.0;
 mxx=mxy=myy=0.0;

 for ( i=0 ; i<wc->nipoint ; i++ )
  {	wi=&wc->ipoints[i];
	w=(img->data[wi->y][wi->x])-(wc->bg);
	x=wi->x-x0,
	y=wi->y-y0;
	if ( w<=0.0 )	continue;
	m  +=w;
	mx +=w*(x+0.5);
	my +=w*(y+0.5);
	mxx+=w*(x*x+x+1.0/3.0);
	mxy+=w*(x+0.5)*(y+0.5);
	myy+=w*(y*y+y+1.0/3.0);
  }

 if ( m<=0.0 )	return(0);
 
 cx=mx/m;
 cy=my/m;

 wc->cx=x0+cx;
 wc->cy=y0+cy;

 cxx=mxx/m-cx*cx;
 cxy=mxy/m-cx*cy;
 cyy=myy/m-cy*cy;

 det=cxx*cyy-cxy*cxy;
 wc->sxx=+cyy/det;
 wc->sxy=-cxy/det;
 wc->syy=+cxx/det;

 return(0);
}

int refine_candidate_params(fitsimage *img,candidate *cands,int ncand)
{
 int	i;

 if ( img==NULL || img->data==NULL )	return(1);

 for ( i=0 ; i<ncand ; i++ )
  {	refine_candidate_data(img,&cands[i]);		}

 return(0);
}

/*****************************************************************************/

int star_set_common_shape_params(double gs,double gd,double gk,star *ws)
{
 double	gg,gm;
 double	gsig,gdev,gdel,gkap,gpa,gellip;

 gg=gs*gs-gd*gd-gk*gk;
 if ( gg<=0.0 )		return(1);
 gm=sqrt(gg);

 gsig=sqrt((gs+gm)/(2*gg));

 if ( gs-gm>0.0 )	gdev=sqrt((gs-gm)/(2*gg));
 else			gdev=0.0;

 gdel=-gd/(2*gsig*gg);
 gkap=-gk/(2*gsig*gg);
 gellip=1.0-(gsig-gdev)/(gsig+gdev);
 if ( gd*gd+gk*gk == 0.0 )	gpa=0.0;
 else				gpa=0.5*atan2(gkap,gdel)*180.0/M_PI;

 ws->gfwhm =gsig*SIG_FWHM;
 ws->gellip=gellip;
 ws->gpa   =gpa;

 ws->gsig=gsig;
 ws->gdel=gdel;
 ws->gkap=gkap;

 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

double star_get_unity_flux(starshape *sw)
{
 double	flux,gg,c,*mom,is;

 switch ( sw->model )
  {  case SHAPE_GAUSS:
	flux=2.0*M_PI/sw->gs;
	break;
     case SHAPE_ELLIPTIC:
	gg=sw->gs*sw->gs-sw->gd*sw->gd-sw->gk*sw->gk;
	if ( gg > 0.0 )	flux=2.0*M_PI/sqrt(gg);
	else		flux=0.0;
	break;
     case SHAPE_DEVIATED:
	mom=sw->mom;
	is=1.0/sw->gs;
	switch ( sw->order )
	 {   case 2: case 3:
		c=1.0+0.5*(mom[0]+mom[2])*is;
		break;
	     case 4: case 5:
		c=1.0+(0.5*(mom[0]+mom[2])+
		  0.25*(0.5*mom[7]+mom[9]+0.5*mom[11])*is)*is;
		break;
	     default:
		c=1.0;
		break;
	 }
	flux=c*2*M_PI*is;
	break;
     default:
	flux=0.0;
  }

 return(flux);
}

/*****************************************************************************/

int convert_candidates(candidate *cands,int ncand,star **rstars,int *rnstar)
{
 star		*stars,*ws;
 candidate	*wc;
 int		nstar,i,j;
 double		gs,gd,gk;

 stars=(star *)malloc(sizeof(star)*ncand);
 nstar=ncand;

 for ( i=0 ; i<ncand ; i++ )
  {	ws=&stars[i];
	wc=&cands[i];

	ws->location.gamp=wc->amp;
	ws->location.gbg =wc->bg;
	ws->location.gcx =wc->cx;
	ws->location.gcy =wc->cy;

	ws->shape.model=SHAPE_ELLIPTIC;
	ws->shape.order=0;
	gs=ws->shape.gs=0.5*(wc->sxx+wc->syy);
	gd=ws->shape.gd=0.5*(wc->sxx-wc->syy);
	gk=ws->shape.gk=wc->sxy;
	ws->shape.gl=0.0;
	for ( j=0 ; j<MAX_DEVIATION_COEFF ; j++ )
	 {	ws->shape.mom[j]=0.0;			}

	#ifdef	STAR_MULTIMODEL
	ws->mshapes=NULL;
	ws->nmshape=0;
	#endif

	star_set_common_shape_params(gs,gd,gk,ws);
	ws->flux=ws->location.gamp*star_get_unity_flux(&ws->shape);
	
	ws->marked=0;

	ws->cand=wc;
  }

 if ( rstars != NULL )	*rstars=stars;
 if ( rnstar != NULL )	*rnstar=nstar;

 return(0);
}

/*****************************************************************************/
