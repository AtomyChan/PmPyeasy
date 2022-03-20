/*****************************************************************************/
/* grmatch.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line utility to match files based on:			     */
/* 	- two-dimensional point sets (point matching),			     */
/*	- any dimensional coordinate sets (coordinate matching), and	     */
/*	- identifiers (identifier matching).				     */
/*****************************************************************************/
#define	FI_GRMATCH_VERSION	"0.9d6"
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>
#include <time.h>

#include "longhelp.h"
#include "fi.h"

#include "io/iof.h"
#include "io/scanarg.h"
#include "io/tokenize.h"
#include "math/fit/lmfit.h"
#include "math/poly.h"
#include "math/polyfit.h"
#include "math/tpoint.h"
#include "math/trimatch.h"
#include "math/cpmatch.h"

#include "transform.h"
#include "common.h"

#ifdef  HAVE_NO_CC_EXTENSION 
#define __extension__ 
#endif 

/*****************************************************************************/

#define		MAX_COORDMATCH_DIM	5	/* To make sizeof(iline)=64, */
						/* however, can be changed   */
						/* to arbitrary dimension... */

#define		MATCH_POINTS		1
#define		MATCH_IDS		2
#define		MATCH_COORDS		3

#define		MATCH_READ_COORDS	0x01
#define		MATCH_READ_ID		0x02

/*****************************************************************************/

int	is_verbose,is_comment;
char	*progbasename;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int fprint_error(char *expr,...)
{
 va_list	ap;
 fprintf(stderr,"%s: error: ",progbasename);
 va_start(ap,expr);
 vfprintf(stderr,expr,ap);
 va_end(ap);
 fprintf(stderr,"\n");
 return(0);
}

int fprint_warning(char *expr,...)
{
 va_list	ap;
 fprintf(stderr,"%s: warning: ",progbasename);
 va_start(ap,expr);
 vfprintf(stderr,expr,ap);
 va_end(ap);
 fprintf(stderr,"\n");
 return(0);
}


/*****************************************************************************/

typedef struct
 {	int	*colcoords;
	int	ncolcoord;

	int	colord;
	int	colwgh;

	int	*colids;
	int	ncolid;
	
	int	neg_ordering;
 } colinfo;

typedef struct
 {	double	x,y;
 } ilinepoint2d;

typedef struct
 {	double	x[MAX_COORDMATCH_DIM];
 } ilinepointnd;

typedef union 
 {	ilinepoint2d	d;
	ilinepointnd	n;
 } ilinepoint;

typedef struct
 {	ilinepoint	p;
	double		ordseq,weight;
	char		*line,*id;
 } iline;

/*****************************************************************************/

int get_number_list(char *strcolid,int **rids,int *rnid)
{
 int	*ids,nid,c;

 ids=NULL;
 nid=0;
 while ( *strcolid )
  {	if ( sscanf(strcolid,"%d",&c)<1 || c<=0 )
	 {	if ( ids != NULL )	free(ids);
		return(1);
	 }
	ids=(int *)realloc(ids,sizeof(int)*(nid+1));
	ids[nid]=c-1;
	nid++;
	while ( *strcolid && *strcolid != ',' )	strcolid++;
	if ( *strcolid==',' )			strcolid++;
  };
 if ( rids != NULL )	*rids=ids;
 if ( rnid != NULL )	*rnid=nid;
 
 return(0);
}

int normalize_columns(colinfo *col,char *strcolcoord,char *strcolid)
{
 if ( col->colord>0 )		col->colord--,col->neg_ordering=0;
 else if ( col->colord<0 )	col->colord=-col->colord-1,col->neg_ordering=1;
 else				col->colord=-1;

 if ( strcolcoord==NULL )
  {	col->ncolcoord=2;
	col->colcoords=(int *)malloc(sizeof(int)*2);
	col->colcoords[0]=0;
	col->colcoords[1]=1;
  }
 else
  {	if ( get_number_list(strcolcoord,&col->colcoords,&col->ncolcoord) )
		return(1);
  }

 if ( strcolid==NULL )
  {	col->ncolid=1;
	col->colids=(int *)malloc(sizeof(int));
	col->colids[0]=0;
  }
 else
  {	if ( get_number_list(strcolid,&col->colids,&col->ncolid) )
		return(1); 
  }

 if ( col->colids==NULL )	return(1);
 else				return(0);
}

static char *wrninappr="Warning: inappropriate content in line %d, skipped.\n";

int cut_newline(char *buff)
{
 for ( ; *buff ; buff++ )
  {	if ( *buff==10 )	*buff=0;	}
 return(0);
}

int read_match_data_points(FILE *fr,colinfo *col,iline **rils,int *rnil,int mtype)
{
 char	*rbuff,*sbuff,**cmd,*id;
 iline	*ils;
 int	nil,n,i,k,l,ln;
 double	xrs[MAX_COORDMATCH_DIM],wo,ww;

 ils=NULL;
 nil=0;

 ln=0;
 rbuff=sbuff=NULL;cmd=NULL;
 while ( ! feof(fr) )
  {	if ( sbuff != NULL )	{ free(sbuff);sbuff=NULL; }
	if ( rbuff != NULL )	{ free(rbuff);rbuff=NULL; }
	if ( cmd   != NULL )	{ free(cmd);  cmd  =NULL; }
	
	rbuff=freadline(fr);
	if ( rbuff==NULL )	break;
	ln++;
	sbuff=strdup(rbuff);
	remove_newlines_and_comments(rbuff);
	cmd=tokenize_spaces_dyn(rbuff);
	if ( cmd==NULL )	continue;
	for ( n=0 ; cmd[n] != NULL ; )	n++;
	if ( n==0 )		continue;

	id=NULL;

	if ( mtype & MATCH_READ_COORDS )
	 {	k=0;
		for ( i=0 ; i<col->ncolcoord ; i++ )
		 {	xrs[i]=0.0;
			if ( col->colcoords[i]<n )
				k+=sscanf(cmd[col->colcoords[i]],"%lg",&xrs[i]);
		 }
		if ( k<col->ncolcoord )
		 {	if ( is_verbose )	fprintf(stderr,wrninappr,ln);
			continue;
		 }
	 }
	if ( mtype & MATCH_READ_ID )
	 {	id=NULL;l=0;
		for ( i=0 ; i<col->ncolid ; i++ )
		 {	k=col->colids[i];
			if ( 0<=k && k<n )
			 {	l+=strlen(cmd[k]);
				if ( id==NULL )
					id=strdup(cmd[k]);
				else
				 {	id=realloc(id,l+2);
					strcat(id,"\n");
					strcat(id,cmd[k]);
				 }
			 }
			else	break;
		 }
		if ( i<col->ncolid )
		 {	if ( id != NULL )	free(id);
			id=NULL;
		 }

		if ( id==NULL )
		 {	if ( is_verbose )	fprintf(stderr,wrninappr,ln);
			continue;
		 }
	 }

	if ( 0<=col->colord && col->colord<n )
	 {	k=sscanf(cmd[col->colord],"%lg",&wo);
		if ( k<1 || ! isfinite(wo) )
		 {	if ( is_verbose )	fprintf(stderr,wrninappr,ln);
			continue;
		 }
		if ( col->neg_ordering )	wo=-wo;
	 }
	else	wo=1.0;

	if ( 0<=col->colwgh && col->colwgh<n )
	 {	k=sscanf(cmd[col->colwgh],"%lg",&ww);
		if ( k<1 || ! isfinite(ww) )
		 {	if ( is_verbose )	fprintf(stderr,wrninappr,ln);
			continue;
		 }
	 }
	else	ww=0.0;

	ils=(iline *)realloc(ils,sizeof(iline)*(nil+1));

	for ( i=0 ; i<col->ncolcoord ; i++ )
	 {	ils[nil].p.n.x[i]=xrs[i];		}
	ils[nil].ordseq=wo;
	ils[nil].weight=ww;

	cut_newline(sbuff);
	ils[nil].line=sbuff,sbuff=NULL;
	ils[nil].id=id;

	nil++;
  };
 if ( sbuff != NULL )	free(sbuff);
 if ( rbuff != NULL )	free(rbuff);
 if ( cmd   != NULL )	free(cmd);
 *rils=ils;
 *rnil=nil;

 return(0);
}

/*****************************************************************************/

/* get_unitarity(): deprecated, use calc_2d_unitarity() from poly.c instead */
/*
double get_unitarity(double ma,double mb,double mc,double md)
{
 double	n1,n2,nn,dd;

 n1=(ma-md)*(ma-md)+(mb+mc)*(mb+mc);
 n2=(ma+md)*(ma+md)+(mb-mc)*(mb-mc);
 nn=(n1<n2?n1:n2);
 dd=ma*ma+mb*mb+mc*mc+md*md;
 if ( dd<=0.0 )		return(-1.0);
 else if ( nn<=0.0 )	return(0.0);
 else			return(sqrt(nn/dd));
}
*/

/*****************************************************************************/

/* debug after trimatch() */ 
/***
 for ( i=0 ; i<nvar ; i++ )
  {	fprintf(stderr,"%12g ",xfit[i]);		}
 fprintf(stderr,"\n");
 for ( i=0 ; i<nvar ; i++ )
  {	fprintf(stderr,"%12g ",yfit[i]);		}
 fprintf(stderr,"\n");
***/

typedef struct 	/* fine-tune parameters of point matching */
 {	int		nmiter,friter;
	double		rejlevel;

	double		maxdist,unitarity;
	int		parity;
	int		ttype,use_ordering,
			maxnum_ref,maxnum_inp;

	int		wcat,w_magnitude;
	double		wpower;

	int		is_centering;
	double		refcx,refcy;
	double		inpcx,inpcy;
	double		maxcenterdist;

	transformation	*htf;
	int		hintorder;

 } matchpointtune;

typedef struct
 {	double	wsigma;		/* weighted RMS of the matched points	     */
	double	nsigma;		/* normal RMS of the matched points	     */
	double	unitarity;	/* unitarity of the fitted polynomial	     */

	double	time_total;	/* total time (in seconds) of point matching */
	double	time_trimatch;	/* time of triangle matching 		     */
	double	time_symmatch;	/* time of symmetric point matching	     */
	int	nmiter;		/* total iterations			     */

	int	tri_level;	/* maximal triangulation level used	     */

 } matchpointstat;

int match_compare_ordering(const void *vi1,const void *vi2)
{
 if ( ((iline*)vi1)->ordseq < ((iline*)vi2)->ordseq ) 	return(1);
 else							return(-1);
}

int do_pointmatch(iline *refls,int nref,iline *inpls,int ninp,
	matchpointtune *mptp,cphit **rhits,int *rnhit,int order,double **vfits,
	matchpointstat *mps)
{
 tpoint		*refps,*rftps,*inpps,*wp;
 point		*pfits;
 tpointarr	arrref,arrinp;
 cphit		*hits;
 int		nhit,i,j,k,nvar,nmax,nmin,ntriref,ntriinp;
 int		i_iter,nmiter,tri_level;
 double		rejlevel,x,y,w;
 double		*xfit,*yfit;
 trimatchpar	tmp;
 trimatchlog	tml;
 iline		*wi;
 double		time_trimatch,time_symmatch,time_total;
 clock_t	t0,t1;

 t0=clock();

 nmiter=mptp->nmiter;
 if ( nmiter<=0 )	nmiter=0;
 rejlevel=mptp->rejlevel;
 if ( rejlevel<=0.0 )	nmiter=0;

 if ( mptp->use_ordering )
  {	qsort(inpls,ninp,sizeof(iline),match_compare_ordering);
	qsort(refls,nref,sizeof(iline),match_compare_ordering);
  } 

 nvar=(order+1)*(order+2)/2;
 xfit=vfits[0];
 yfit=vfits[1];

 refps=(tpoint *)malloc(sizeof(tpoint)*nref);
 inpps=(tpoint *)malloc(sizeof(tpoint)*ninp);
 rftps=(tpoint *)malloc(sizeof(tpoint)*nref);

 hits=NULL;

 time_trimatch=0.0;
 time_symmatch=0.0;

 tri_level=0;

 for ( i_iter=0 ; i_iter<=nmiter ; i_iter++ )
  {	clock_t	t0,t1;

	if ( i_iter<=0 )
	 {	for ( i=0 ; i<nref ; i++ )
		 {	wp=&refps[i];
 			wp->id=i;		/* trivial id */
			wp->xcoord=refls[i].p.d.x,
			wp->ycoord=refls[i].p.d.y;
		 }	
	 }
	else
	 {	for ( i=0 ; i<nref ; i++ )
		 {	wp=&refps[i];
			wp->id=i;
			x=refls[i].p.d.x,
			y=refls[i].p.d.y;
			wp->xcoord=eval_2d_poly(x,y,order,xfit,0,0,1);
			wp->ycoord=eval_2d_poly(x,y,order,yfit,0,0,1);
		 }
	 }

	for ( i=0 ; i<ninp ; i++ )
	 {	wp=&inpps[i];
	 	wp->id=i;		/* trivial id */
		wp->xcoord=inpls[i].p.d.x,
		wp->ycoord=inpls[i].p.d.y;
	 }

	hits=NULL;

	tmp.level    =mptp->ttype;
	tmp.maxdist  =-1.0;
	tmp.unitarity=mptp->unitarity;
	tmp.parity   =mptp->parity;

	if ( ninp>nref )	nmax=ninp,nmin=nref;
	else			nmax=nref,nmin=ninp;

	if ( mptp->use_ordering && mptp->maxnum_inp>0 && mptp->maxnum_ref>0 )
	 {	if ( nref>=mptp->maxnum_ref && ninp>=mptp->maxnum_inp )
		 {	ntriref=mptp->maxnum_ref;
			ntriinp=mptp->maxnum_inp;
		 }
		else
		 {	double	nimr,nrmi;
			ntriref=nref;
			ntriinp=ninp;
			nimr=(double)ntriinp*(double)mptp->maxnum_ref;
			nrmi=(double)ntriref*(double)mptp->maxnum_inp;
			if ( nimr<=nrmi )
			 {	ntriref=(int)(nimr/(double)mptp->maxnum_inp);	}
			else
			 {	ntriinp=(int)(nrmi/(double)mptp->maxnum_ref);	}
			if ( ntriref>mptp->maxnum_ref )	ntriref=mptp->maxnum_ref;
			if ( ntriinp>mptp->maxnum_inp )	ntriinp=mptp->maxnum_inp;
			if ( ntriref>nref )		ntriref=nref;
			if ( ntriinp>ninp )		ntriinp=ninp;
		 }
	 }
	else if ( mptp->use_ordering )
	 {	ntriref=nmin;
		ntriinp=nmin;
	 }
	else
	 {	ntriref=nref,
		ntriinp=ninp;
	 }

	/*fprintf(stderr,"nref=%d,ninp=%d;mxref=%d,mxinp=%d;ntriref=%d,ntriinp=%d\n",nref,ninp,mptp->maxnum_ref,mptp->maxnum_inp,ntriref,ntriinp);*/

	if ( hits != NULL )
	 {	free(hits);
		hits=NULL;
	 }

	t0=clock();
	trimatch(refps,ntriref,inpps,ntriinp,order,&tmp,&hits,&nhit,xfit,yfit,&tml);
	t1=clock();

	time_trimatch+=(double)(t1-t0)/(double)CLOCKS_PER_SEC;
	if ( tml.level_used > tri_level )	tri_level=tml.level_used;

	/*
	fprintf(stderr,"x: -> %12g %12g %12g\n",xfit[0],xfit[1],xfit[2]);
	fprintf(stderr,"y: -> %12g %12g %12g\n",yfit[0],yfit[1],yfit[2]);
	*/

	for ( i=0 ; i<nref ; i++ )
	 {	x=refps[i].xcoord,
		y=refps[i].ycoord;
		rftps[i].xcoord=eval_2d_poly(x,y,order,xfit,0,0,1);
 		rftps[i].ycoord=eval_2d_poly(x,y,order,yfit,0,0,1);
		rftps[i].id=refps[i].id;
	 }
	if ( hits != NULL )	free(hits);

	arrref.points=rftps,arrref.length=nref;
	arrinp.points=inpps,arrinp.length=ninp;

	t0=clock();
	hits=cpmatch_symmetric(&arrref,&arrinp,&nhit,0.0,mptp->maxdist);
	/*hits=cpmatch_symmetric(&arrref,&arrinp,&nhit,3.0,0.0);*/ /* old */
	t1=clock();
	time_symmatch+=(double)(t1-t0)/(double)CLOCKS_PER_SEC;

	pfits=(point *)malloc(sizeof(point)*nhit);
	for ( i=0 ; i<nhit ; i++ )
	 {	j=hits[i].idx[0];
	 	pfits[i].x=refls[j].p.d.x,
		pfits[i].y=refls[j].p.d.y,
		pfits[i].weight=1.0;
	 }
	if ( mptp->wcat>=0 )
	 {	for ( i=0 ; i<nhit ; i++ )
		 {	if ( ! mptp->wcat )	wi=&refls[hits[i].idx[0]];
			else			wi=&inpls[hits[i].idx[1]];
			w=wi->weight;
			if ( mptp->w_magnitude )
			 {	w=exp(-0.4*M_LN10*(w-20.0));	}
			if ( mptp->wpower != 1.0 )
			 {	w=pow(w,mptp->wpower);		}
			pfits[i].weight=w;
		 }
	 }
	for ( i=0 ; i<nhit ; i++ )
	 {	pfits[i].value=inpls[hits[i].idx[1]].p.d.x;		}
	fit_2d_poly(pfits,nhit,order,xfit,0,0,1);
	for ( i=0 ; i<nhit ; i++ )
	 {	pfits[i].value=inpls[hits[i].idx[1]].p.d.y;		}
	fit_2d_poly(pfits,nhit,order,yfit,0,0,1);

	if ( mptp->is_centering && mptp->maxcenterdist>0.0 )
	 {	x=eval_2d_poly(mptp->refcx,mptp->refcy,order,xfit,0,0,1)-mptp->inpcx;
		y=eval_2d_poly(mptp->refcx,mptp->refcy,order,yfit,0,0,1)-mptp->inpcy;
		if ( x*x+y*y > mptp->maxcenterdist*mptp->maxcenterdist )
			break;
	 }
	
  }

 t1=clock();
 time_total=(double)(t1-t0)/(double)CLOCKS_PER_SEC;

 if ( mps != NULL )
  {	mps->time_total=time_total;
	mps->time_trimatch=time_trimatch;
	mps->time_symmatch=time_symmatch;
	mps->nmiter=nmiter+1;
	mps->tri_level=tri_level;
  }

 if ( mps != NULL && nhit>0 )
  {	double	ws,ns,wdd,ndd,nx,ny,dx,dy,dd;

	ws=wdd=0.0;
	ns=ndd=0.0;
	for ( i=0 ; i<nhit ; i++ )
	 {	j=hits[i].idx[0];
		k=hits[i].idx[1];

		if ( mptp->wcat<0 )		wi=NULL;
		else if ( ! mptp->wcat )	wi=&refls[j];
		else				wi=&inpls[k];
		if ( wi != NULL )
		 {	w=wi->weight;
			if ( mptp->w_magnitude )
			 {	w=exp(-0.4*M_LN10*(w-20.0));	}
			if ( mptp->wpower != 1.0 )
			 {	w=pow(w,mptp->wpower);		}		
		 }
		else
			w=1.0;

		x=refls[j].p.d.x,
		y=refls[j].p.d.y;
		nx=eval_2d_poly(x,y,order,xfit,0,0,1);
		ny=eval_2d_poly(x,y,order,yfit,0,0,1);
		dx=nx-inpls[k].p.d.x;
		dy=ny-inpls[k].p.d.y;
		dd=dx*dx+dy*dy;
	
		ws+=w,  wdd+=dd*w;
		ns+=1.0,ndd+=dd;
	 }
	if ( ws>0.0 )	mps->wsigma=sqrt(wdd/ws);
	else		mps->wsigma=0.0;
	if ( ns>0.0 )	mps->nsigma=sqrt(ndd/ns);
	else		mps->nsigma=0.0;

	mps->unitarity=calc_2d_unitarity(xfit,yfit,order);
  }
 else if ( mps != NULL )
  {	mps->wsigma=0.0;
	mps->nsigma=0.0;
	mps->unitarity=-1.0;
  }

 free(rftps);
 free(inpps);
 free(refps);

 if ( rhits != NULL )	*rhits=hits;
 if ( rnhit != NULL )	*rnhit=nhit;

 return(0); 
}

/*****************************************************************************/

int match_compare_firstcoord(const void *vi1,const void *vi2)
{
 if ( ((iline*)vi1)->p.n.x[0] < ((iline*)vi2)->p.n.x[0] ) 	return(-1);
 else								return(1);
}

int coordmatch_search_nearest(iline *ils,int nil,double *x0,int dim)
{
 int	best,min,max,mid,i,j,brk;
 double	*xc,xc0[MAX_COORDMATCH_DIM],
	xx,cc,dd,dist[MAX_COORDMATCH_DIM],cdist,edist,mindist,maxdist;

 int	outl,outr,pl,pr;

 if ( ils==NULL || nil<=0 || x0==NULL || dim<=0 || dim>MAX_COORDMATCH_DIM )
	return(0);
 
 min=mid=0;
 max=nil;
 while ( max>min )
  {	mid=(min+max)/2;
	cdist=x0[0]-ils[mid].p.n.x[0];
	if ( fabs(cdist)<1e-10 )	break;
	else if ( cdist>0.0 )		min=mid+1;
	else				max=mid;
  }

 /* one dimension: it is a special (and simple) case... */
 if ( dim==1 )
  {	best=mid;
	mindist=fabs(x0[0]-ils[mid].p.n.x[0]);
	if ( mid>0 )
	 {	cdist=fabs(x0[0]-ils[mid-1].p.n.x[0]);
		if ( cdist<mindist )	mindist=cdist,best=mid-1;
	 }
	if ( mid<nil-1 )
	 {	cdist=fabs(x0[0]-ils[mid+1].p.n.x[0]);
		if ( cdist<mindist )	mindist=cdist,best=mid+1;
	 }
	return(best);
  }

 /* two and more dimension: do the 'real' coordinate matching */
 mindist=maxdist=0.0;
 for ( i=0 ; i<dim ; i++ )
  {	xc0[i]=ils[mid].p.n.x[i];
	dist[i]=fabs(x0[i]-xc0[i]);
	mindist+=dist[i]*dist[i];
	maxdist+=dist[i];
  }

 best=mid;
 for ( outl=outr=1,pl=best-1,pr=best+1 ; outl || outr ; pl--,pr++ )
  {	if ( outl )
	 {	if ( pl>=0 )
		 {	xc=&ils[pl].p.n.x[0];
			xx=fabs(x0[0]-xc[0]);
			if ( xx<maxdist )
			 {	brk=0;
				cdist=xx;
				for ( i=1 ; i<dim ; i++ )
			 	 {	cc=fabs(xc[i]-x0[i]);
					cdist+=cc;
					if ( cc<=dist[i] )
					 {	dist[i]=cc,brk=1;	}
				 }
				if ( brk )
				 {	edist=xx*xx;
					for ( j=1 ; j<dim ; j++ )
					 {	dd=xc[j]-x0[j];
						edist+=dd*dd;
				 	 }
					if ( edist<=mindist )
					 {	best=pl,mindist=edist;
						if ( cdist<maxdist )
							maxdist=cdist;
					 }
				 }
			 }
			else	outl=0;
		 }
		else	outl=0;
	 }
	if ( outr )
	 {	if ( pr<nil )
		 {	xc=&ils[pr].p.n.x[0];
			xx=fabs(x0[0]-xc[0]);
			if ( xx<maxdist )
			 {	brk=0;
				cdist=xx;
				for ( i=1 ; i<dim ; i++ )
			 	 {	cc=fabs(xc[i]-x0[i]);
					cdist+=cc;
					if ( cc<=dist[i] )
					 {	dist[i]=cc,brk=1;	}
				 }
				if ( brk )
				 {	edist=xx*xx;
					for ( j=1 ; j<dim ; j++ )
					 {	dd=xc[j]-x0[j];
						edist+=dd*dd;
				 	 }
					if ( edist<=mindist )
					 {	best=pr,mindist=edist;
						if ( cdist<maxdist )
							maxdist=cdist;
					 }
				 }
			 }
			else	outr=0;
		 }
		else	outr=0;
	 }
  }

 return(best);
}

static double get_distance(double *x1,double *x2,int dim)
{
 double	r;
 for ( r=0.0 ; dim>0 ; dim--,x1++,x2++ )
  {	r+=((*x1)-(*x2))*((*x1)-(*x2));		}
 return(r);
}

int do_coordmatch(iline *refls,int nref,iline *inpls,int ninp,cphit **rhits,int *rnhit,int dim,double maxdist)
{
 cphit	*hits;
 int	nmax,nmin,*refhs,*inphs,i,j,nhit;
 double	mxd2;

 if ( refls==NULL || nref<=0 )	return(-1);
 if ( inpls==NULL || ninp<=0 )	return(-1);

 qsort(refls,nref,sizeof(iline),match_compare_firstcoord);
 qsort(inpls,ninp,sizeof(iline),match_compare_firstcoord);

 nmax=(nref>ninp?nref:ninp);
 nmin=(nref>ninp?ninp:nref);
 
 hits=(cphit *)malloc(sizeof(cphit)*nmin);
 
 refhs=(int *)malloc(sizeof(int)*nref);
 inphs=(int *)malloc(sizeof(int)*ninp);

 for ( i=0 ; i<nref ; i++ )
  {	refhs[i]=coordmatch_search_nearest(inpls,ninp,refls[i].p.n.x,dim);	}
 for ( i=0 ; i<ninp ; i++ )
  {	inphs[i]=coordmatch_search_nearest(refls,nref,inpls[i].p.n.x,dim);	}

 nhit=0;
 mxd2=maxdist*maxdist;
 for ( i=0 ; i<nref ; i++ )
  {	j=refhs[i];
	if ( inphs[j]==i && ( maxdist<0 || get_distance(refls[i].p.n.x,inpls[j].p.n.x,dim)<=mxd2 ) )
	 {	hits[nhit].idx[0]=i;
		hits[nhit].idx[1]=j;
		nhit++;
	 }
  }

 free(inphs);
 free(refhs);

 hits=(cphit *)realloc(hits,sizeof(cphit)*nhit);

 if ( rhits != NULL )	*rhits=hits;
 if ( rnhit != NULL )	*rnhit=nhit;
 
 return(0);
}

/*****************************************************************************/

int stridcmp(char *id1,char *id2)
{
 if ( id1==NULL && id2==NULL )	return(0);
 else if ( id1==NULL )	return(1);
 else if ( id2==NULL )	return(-1);
 else return(strcmp(id1,id2));
}

int id_compare(const void *vi1,const void *vi2)
{ 
 return ( stridcmp(((iline*)vi1)->id,((iline*)vi2)->id) );
}

/*
int search_id(iline *inpls,int n,char *id)
{
 int	f,k,w;
 f=0;
 while ( n )
  {	k=f+n/2;
	w=stridcmp(inpls[k].id,id);
	if ( w==0 )	return(k);
	else if ( w<0 )	f+=n/2+1,n-=n/2+1;
	else		n=n/2;
  };
 return(-1);
}
*/

int search_id_limiters(iline *inpls,int n,char *id,int *rleft,int *rright)
{
 int	min,mid,max;

 if ( n<=0 )
	return(1);
 if ( rleft==NULL || rright==NULL )
	return(-1);

 min=0,max=n;
 if ( ! ( 0<=stridcmp(inpls[n-1].id,id) ) )
  {	*rleft=1,*rright=0;
	return(1);
  }
 while ( max>min+1 )
  {	mid=(min+max)/2;
	if ( 0<=stridcmp(inpls[mid-1].id,id) )
		max=mid;
	else
		min=mid;
  };
 *rleft=min;

 min=0,max=n;
 if ( ! ( stridcmp(inpls[0].id,id)<=0 ) )
  {	*rleft=1,*rright=0;
	return(1);
  }
 while ( max>min+1 )
  {	mid=(min+max)/2;
	if ( stridcmp(inpls[mid].id,id)<=0 )
		min=mid;
	else
		max=mid;
  };	
 *rright=min;

 if ( *rleft>*rright )	return(1);
 else			return(0);
}

#define		AMBIG_NONE	0
#define		AMBIG_FIRST	1
#define		AMBIG_ANY	2
#define		AMBIG_FULL	3

int do_idmatch(iline *refls,int nref,iline *inpls,int ninp,cphit **rhits,int *rnhit,int ambig)
{
 int	i,j,k,nhit,ahit,hr,hi,rl,rr,il,ir;
 cphit	*hits;
 char	*id;

 qsort(inpls,ninp,sizeof(iline),id_compare);
 qsort(refls,nref,sizeof(iline),id_compare);

 ahit=nref;
 hits=(cphit *)malloc(sizeof(cphit)*ahit);
 nhit=0;

 for ( i=0 ; i<nref ; )
  {	id=refls[i].id;
	if ( id==NULL )
	 {	i++;
		continue;
	 }
	hr=search_id_limiters(refls,nref,id,&rl,&rr);
	if ( hr )
	 {	i++;
		continue;
	 }
	hi=search_id_limiters(inpls,ninp,id,&il,&ir);
	if ( hi )
	 {	i+=rr-rl+1;
		continue;
	 }
	
	hr=rr-rl+1;
	hi=ir-il+1;

	switch ( ambig )
	 {   case AMBIG_FIRST:
		k=1;
		break;
	     case AMBIG_ANY:
		k=(hr<hi?hr:hi);
		break;
	     case AMBIG_FULL:
		k=hr*hi;
		break;
	     case AMBIG_NONE:
	     default:
		if ( hr==1 && hi==1 )	k=1;
		else			k=0;
	 }

	if ( k<=0 )
	 {	i+=rr-rl+1;
		continue;
	 }

	if ( nhit+k>ahit )
	 {	ahit=nhit+k;
		hits=(cphit *)realloc(hits,sizeof(cphit)*ahit);
	 }

	switch ( ambig )
	 {   case AMBIG_FIRST:
		hits[nhit].idx[0]=rl;
		hits[nhit].idx[1]=il;
		nhit++;
		break;
	     case AMBIG_ANY:
		for ( j=0 ; j<hr && j<hi ; j++ )
		 {	hits[nhit].idx[0]=rl+j;
			hits[nhit].idx[1]=il+j;
			nhit++;
		 }
		break;
	     case AMBIG_FULL:
		for ( j=0 ; j<hr ; j++ )
		 {	for ( k=0 ; k<hi ; k++ )
			 {	hits[nhit].idx[0]=rl+j;
				hits[nhit].idx[1]=il+k;
				nhit++;
			 }
		 }
		break;

	     case AMBIG_NONE:
	     default:
		if ( hr==1 && hi==1 )
	 	 {	hits[nhit].idx[0]=rl;
			hits[nhit].idx[1]=il;
			nhit++;
		 }
		break;
	 }

	i+=rr-rl+1;
	
  }
			
 if ( rhits != NULL )	*rhits=hits;
 if ( rnhit != NULL )	*rnhit=nhit;
 
 return(0);
}

/*****************************************************************************/

int fprint_grmatch_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\tgrmatch [-h|--help] [-C|--comment] [-V|--verbose]\n"
"\t-r|--input-reference <reffile> [[-i|--input] <in>] [-o|--output <out>]\n"
"Point matching (default, --match-points can be omitted):\n"
"\t--match-points [--col-ref <>,<>] [--col-inp <>,<>] \n"
"\t[--col-ref-ordering [-]<> --col-inp-ordering [-]<>]\n"
"\t[--output-transformation <output-transformation>]\n"
"\t[-a|--order <order>] [-f|--offset <ox>,<oy>] [--scale <scale>]\n");
 fprintf(fw,
"Identifier matching:\n"
"\t--match-id [--col-ref-id <>[,<>[,...]]] [--col-inp-id <>[,<>[,...]]]\n"
"\t[--{no|first|any|full}-ambiguity]\n");
 fprintf(fw,
"Coordinate matching:\n"
"\t--match-coords [--col-ref <>[,<>[,...]]] [--col-inp <>[,<>[,...]]]\n"
"\t[-m|--max-distance <maxdist>]\n");
 fprintf(fw,
"Outputs:\n"
"\t[--output|--output-matched|-o <line-pairs>] [--output-id <id-pairs>]\n"
"\t[--output-excluded-{reference|input} <excluded-lines>]\n");
 fprintf(fw,
"Fine-tuning of point matching:\n"
"\t[-H|--hint-transformation <transformation> [-b|--hint-order <order>]]\n"
"\t[--triangulation delaunay|level=<level>|full|auto[,unitarity=<factor>],\n"
"\t                 mixed|conformable|reverse \n"
"\t                 max{number|ref|inp}=<max.number>]\n"
"\t[-m|--max-distance <maxdist>]\n"
"\t[--fit iterations=<num.iter.>,firstrejection=<1st.rej.>,sigma=<sigma>]\n"
"\t[--weight [reference|input],column=<>,[magnitude],[power=<power>]\n"
"\t[--center-{reference|input} <x>,<y> [--center-maxdist <mxdist>]]\n");
 return(0);
}

longhelp_entry grmatch_long_help[] = 
{
 LONGHELP_OPTIONS,

 { "General options:", NULL },
 { "-h, --help",
	"Give general summary about the command line options." },
 { "--long-help",
	"Give a detailed list of command line options." },
 { "--version",
	"Give some version information about the program." },
 { "-C, --comment",
	"Comment the output (both the transformation file and the match file)." },

 { "Options for input/output specifications:", NULL },
 { "-r <referencefile>, --input-reference <referencefile>",
	"Mandatory, name of the reference file." },
 { "<inputfile>, -i <inputile>, --input <inputfile>",
	"Name of the input file. If this switch is omitted, the input is"
	"read from stdin (specifying some input is mandatory)." },
 { "-o <output>, --output <output>, --output-matched <output>",
	"Name of the output file, containing the matched lines. The "
	"matched lines are pasted lines, the first part is from the "
	"reference file and the second part is from the input file, these "
	"two parts are concatenated by a TAB character. This switch is "
	" optional, if it is not specified, no such output will "
	"be generated. " },
 { "--output-excluded-reference <out>, --output-excluded-input <out>",
	"Names of the files which contain the valid but excluded lines "
	"from the reference and from the input. These outputs are "
	"disjoint from the previous output and altogether contaions all "
	"valid lines." },
 { "--output-id <out>",
	"Name of the file which contaions only the identifiers of the "
	"matched lines. If the primary matching method was not "
	"identifier matching, one should specify the column indices of the "
	"identifiers by --col-ref-id and --col-inp-id also." },
 { "--output-transformation <output-transformation-file>",
	"Name of the output file containing the  geometrical  transformation, "
	"in human-readable format, if the matching method was point "
	"matching (in other case, this option has no  effect).  The "
	" commented  version  of this file includes some statistics about "
	"the matching (the total  number  of  lines  used  and  matched,  "
	"the required CPU time, the final triangulation level, the fit "
	"residuals and other things like these)." },
 { "In all of the above input/output file specifications,  the  replacement "
   "of  the  file name by ``-'' (a single minus sign) forces the reading from "
   "stdin or writing to stdout. Note that all parts of the any  line  after "
   "``#'' (hashmark) are treated as a comment, therefore "
   "ignored. ", NULL }, { "", NULL },

 { "General options for point matching:", NULL },
 { "--match-points",
	"This  switch  forces  the usage of the point matching method. By "
	"default, this method is  assumed  to  be  used,  therefore  this "
	"switch can be omitted." },
 { "--col-ref <x>,<y>, --col-inp <x>,<y>",
	"The  column  indices containing the X and Y coordinates, for the "
	"reference and for the input file, respectively. The index of the "
	"first  column  is  always 1, the index of the second is 2 and so "
	"on. Lines in which these columns do not contain valid real numbers "
	"bers are omitted." },
 { "-a <order>, --order <order>",
	"This  switch specifies the polynomial order of the resulted "
	"geometrical transformation. It can be arbitrary  positive  integer. "
	"Note that if the order is A, at least (A+1)*(A+2)/2 valid points "
	"are needed both from the reference and both from the input  file "
	"to fit the transformation. " },
 { "--max-distance <maxdist>",
	__extension__
	"The  maximal accepted distance between the matched points in the "
	"coordinate frame of the input coordinate list (and not in the "
	"coordinate frame of the reference coordinate list). Possible pairs "
	"(which are valid pairs due to "
	"the  symmetric  coordinate  matching  algorihms) are excluded if "
	"their Eucledian distance is larger than maxdist. Note that  this "
	"option has no initial value, therefore, if omitted, all possible "
	"pairs due to the symmetric matching are resulted, which, in certain "
	"cases  in  practice,  can result unexpected behaviour. One "
	"should always specify a reasonable maximal distance which can be "
	"estimated  only  by  the  knowledge  of the physics of the input "
	"files. " },
 { "See more options concerning to point  matching  in  the  section "
   "``Fine-Tuning   of  Point  Matching'' below. That  section  also "
   "describes the tuning of the  triangulation  used  by  the  point "
   "matching  algorithm.  For  a more detailed description about the "
   "point matching algorithms based on pattern and triangle matching "
   "see [1], [2] or [3].", NULL }, { "", NULL },

 { "General options for coordinate matching:", NULL },
 { "--match-coord, --match-coords",
	"This  switch forces the usage of the coordinate matching method. "
	"Note that because of the common options with the point  matching "
	"method, one should specify this switch to force the usage of the "
	"coordinate matching method (the default method is  point  matching, "
	"see above)." },
 { "--col-ref <x>[,<y>,[<z>...]] --col-inp <x>[,<y>,[<z>...]]",
	__extension__
	"The  column  indices containing the spatial coordinates, for the "
	"reference and for the input file, respectively. The index of the "
	"first  column  is  always 1, the index of the second is 2 and so "
	"on. Lines in which these columns do not contain valid real  numbers "
	"are  omitted.  Note  that  the dimension of the coordinate "
	"matching space is specified indirectly, by the number of  column "
	"indices  listed  here.  Because  of  this,  the number of column "
	"indices should be the same for the reference and input, in other "
	"case,  when  the  dimensions  are  mismatched, the program exits "
	"unsuccessfully." },
 { "--max-distance <maxdist>",
	"The maximal accepted distance between the matched points. "
	"Possible  pairs (which are valid pairs due to the symmetric "
	"coordinate matching algorihms) are excluded if  their  Eucledian "
	"distance  is larger than maxdist. Note that this option has no "
	"initial value, therefore, if omitted, all possible pairs due to "
	"the  symmetric  matching  are  resulted (see also point matching, "
	"above)." }, { "", NULL },
 { "General options for identifier matching:", NULL },
 { "--match-id, --match-identifiers",
	"This switch forces the usage of the identifier matching  method." },
 { "--col-ref-id <i>[,<j>,[<k>...]] --col-inp-id <i>[,<j>,[<k>...]]",
	"Column  index  or  indices  containing the identifiers, from the "
	"reference and from the input file, respectively." },
 { "--no-ambiguity, --first-ambiguity, --any-ambiguity, --full-ambiguity",
	__extension__
	"These options tune the behaviour of the matching when  there is "
	"more  than one occurrence of a given identifier in the reference "
	"and/or input file.  If --no-ambiguity is specified, these  "
	"identifiers  are discarded, this is the default method.  "
	"If --first-ambiguity is specified, only the first occurence is "
	"treated as a matched  line, independently from the number of "
	"occurrences.  If the switch --any-ambiguity is specified, the "
	"lines  are  paired sequentally, until there is any left from the "
	"reference and from the input.  For example, if there is 4 "
	"occurrences in the reference  and  6  in the input file of a given "
	"identifier, 4 matched pairs are returned.  Otherwise, if  "
	"--full-ambiguity  is  specified,  all  possible  combinations  "
	"of  the lines are treated as matched lines. For example, if "
	"there is  4  occurrences  in  the reference  and  6  in  the "
	"input file of a given identifier, all 4*6=24 combinations are "
	"returned as matched pairs." },
 
 { "Fine-tuning of point matching:", NULL },
 { "--triangulation <parameters>",
	"This switch  is  followed  by  comma-separated  directives, which "
	"specify the parameters of the triangulation-based point matching "
	"algorithm: " },
 { "delaunay, level=<level>, full, auto, unitarity=<U>",
	__extension__
	"These  directives specify the triangulation level used for point "
	"matching. ``delaunay'' forces the usage only of the "
	"Delaunay-triangles.  This is the fastest method, however, it is "
	"only working if the points in the reference and input lists are "
	"almost  competely  overlapping  and  describe  almost  the  same "
	"point sets (within a ratio of common  points  above  60-70%).  "
	"The  ``level'' specifies  the level of the expansion of the "
	"Delaunay-triangulation (see [1] for more details).  In  practice,  "
	"the  lower  the ratio  of common points and/or the ratio of the "
	"overlapping, the higher level should be used.  Specifying "
	"``level=1'' or  ``level=2'' gives  a  robust  but  still  fast "
	"method for general usage. The directive ``full'' forces full "
	"triangulation.  This can  be  overwhelmingly  slow  and  annoying  "
	"and  requires tons of memory if there are more than 40-50 points "
	"(the amounts of these resources are  proportional  to  the 6th(!) "
	"and 3rd power of the number of the points, respectively). The "
	"directive  ``auto''  increases  the level  of  the  triangulation  "
	"expansion  automatically  until a proper match is found. A match "
	"is considered as a good match if the unitarity of the transformation "
	"is less than the unitarity U specified by the ``unitarity=U'' "
	"directive (see also  the  section Notes/Unitarity below). " },
 { "mixed, conformable, reverse",
	__extension__
	"These directives define the chirality of the triangle spaces to "
	"be used.  Practically, it means the following. If we don't  know "
	"whether the input and reference lists are inverted respecting to "
	"each other, one should use ``mixed'' triangle  space.  If  we  "
	"are sure  about that the input and reference lists are not inverted, "
	"we can use ``conformable'' triangle space. If  we  know  that  the "
	"input  and  reference  lists  are inverted, we can use ``reverse'' "
	"space. Note that although  ``mixed''  triangle  space  can  always "
	"result  a  good match, it is a wise idea to fix the chirality by "
	"specifying ``conformable'' or ``reverse'' if we really know that "
	"the point  sets  are  not  inverted  or  inverted respecting "
	"to each other. If the  chirality  is  fixed,  the  program  yields "
	"more matched  pairs,  the  appropriate  triangulation  level  can  "
	"be smaller and in ``auto'' mode, the program returns the match  "
	"definitely faster." },
 { "maxnumber=<max>, maxref=<mr>, maxinp=<mi>",
	__extension__
	"These directives specify the maximal number of points which are "
	"used for triangulation (for  any  type  of  triangulation). "
	"If ``maxnumber''  is  specified,  it is equivalent to define "
	"``maxref'' and ``maxinp'' with the same values. Then,  the  "
	"first  <mr>  points from  the  reference and the first <mi> points "
	"from the input list are used to generate the triangle sets. "
	"The ``first''  points  are selected  using  the  optional  "
	"information  found in one of the columns, see the "
	"following switches." },

 { "(Note that there should be only one --triangulation switch, all desired "
   "directives  should  be  written in the same argument, "
   "separated by commas.)", NULL },

 { "--col-ref-ordering [-]<w>, --col-inp-ordering [-]<w>.",
	__extension__
	"These switches specify one-one column index from  the  reference "
	"and from the input files which are used to order these lists and "
	"select the first ``maxref'' and ``maxinp'' points  (see  above)  "
	"for the  generation  of the two triangle meshes. Both columns "
	"should contain valid real  numbers,  otherwise  the  whole(!)  "
	"line  is excluded (not only from sorting but from the whole matching "
	"procedure). If there is no negative sign before the  column  index, "
	"the  data are sorted in descending(!) order, therefore the lines "
	"with the lines with the highest(!) values are selected for  "
	"triangulation.  If  there  is a negative sign before the index, "
	"the data are sorted in ascending order by  these  values,  "
	"therefore the  lines with the smallest(!) values are selected for "
	"triangulation. For example, if we want to match star  lists,  we  "
	"might want  to  use  only  the brightest ones to generate the "
	"triangle sets. If the brightnesses of the stars are  specified  "
	"by  their fluxes,  we should not use the negative sign (the list "
	"should be sorted in descending order to select the first few "
	"lines as  the brightest  stars),  and if the brightness is known "
	"by the magnitude, we have to use the negative sign." },
 
 { "--fit iterations=<N>,firstrejection=<F>,sigma=<S>",
	__extension__
	"Like --triangulation, this switch is  followed  by  some  "
	"directives.  These  directives  specify  the  number <N> of "
	"iterations (``iterations=<N>'')  for  point  matching.  The  "
	"``firstrejection'' directive  speciy  the  serial  number <F> "
	"of the first iteration where points farer than <S> ``sigma'' level "
	"are excluded in the next iteration.  Note  that  in  practice  "
	"these type of iteration is really not important (due to, for instance, "
	"the limitations of the outliers by the --max-distance switch), "
	"however, some suspicious users can be convinced by such arguments." },
 { "--weight reference|input,column=<wi>,[magnitude],[power=<p>]",
	__extension__
	"These  directives  specify the weights which are used during the "
	"fit of the geometrical transformation. For example, in  practice "
	"it  is  useful  in the following situation. We try to match star "
	"lists, then the fainter stars are believed to have higher "
	"astrometrical errors, therefore they should have smaller influence "
	"in the fit. We can take the weights  from  the  reference  "
	"(specify ``reference'') and from the input (specify ``input''), "
	"from the column specified by the weight-index. The weights  can  "
	"be  derived from  stellar  magnitudes, if so, specify ``magnitude'' "
	"to convert the read values in magnitude to flux. The real weights  "
	"then  is the  ``power''th  power  of  the  flux.  The  default "
	"value of the ``power'' is 1, however, for the maximum-likelihood "
	"estimation  of an assumed Gaussian distribution, the weights "
	"should be the second power of the fluxes." },

 { __extension__
   "Some notes on unitarity.  The unitarity of a geometrical transformation "
   "measures  how it differs from the closest transformation which is affine "
   "and a combination of dilation, rotation and shift. For such a "
   "transformation  the unitarity  is  0 and if the second-order terms in "
   "a transformation distort a such unitary transformation, the unitarity "
   "will  have  the  same magnitude  like the magnitude of this second-order "
   "effect. For example, to map a part of a sphere with the size of d degrees "
   "will have an  unitarity of 1-cos(d). Therefore, for astrometrical "
   "purposes, a reasonable value of the critical unitarity in ``auto'' "
   "triangulation  mode  can  be estimated  as  2 or 3 times 1-cos(d/2) "
   "where d is the size of the field in which astrometry "
   "should be performed.", NULL },
 
 { NULL, NULL }
 
};

int fprint_grmatch_long_help(FILE *fw)
{
 fprintf(fw,
"Usage:\tgrmatch [options] -r <reference> -i <input> [-o <output>]\n");
 fprintf(fw,
"The program `grmatch` matches lines read from two input files, namely\n"
"from a reference and from an input file. All implemented algorithms are\n"
"symmetric, in the manner that the result should be the same if these\n"
"two files are\n"
"swapped. The only case when the order of these files is important is when\n"
"a geometrical transformation is also returned (see point matching below),\n"
"in this case the swapping of the files results the inverse form of the\n"
"original transformation.\n");
 fprintf(fw,
"The lines (rows) can be matched using various criteria.\n\n");
 fprintf(fw,
"1. Lines can be matched by identifier, where the identifier can\n"
"be any concatenation of arbitrary, space-separated columns found in the\n"
"files. Generally, the identifier is represented by a single column\n"
"(e.g. it is an astronomical catalog identifier). The behaviour of the\n"
"program can be tuned for the cases when there are more than one rows\n"
"with the same identifier.\n\n");
 fprintf(fw, __extension__
"2. Lines can be matched using a 2-dimensional point matchig algorithm.\n"
"In this method, the program expects two-two columns both from the \n"
"reference and input files which can be treated as X and Y coordinates.\n"
"If both point lists are known, the program tries to find the appropriate\n"
"geometrical transformation which transforms the points from the \n"
"frame of the reference list to the frame of the input list and, \n"
"simultaneously, tries to find as many pairs as possible. The \n"
"parameters of the geometrical transformation and the whole algorithm\n"
"can be fine-tuned.\n\n");
 fprintf(fw,
"3. Lines can be matched using arbitrary- (N-) dimensional coordinate matching\n"
"algorithm. This method expects N-N columns both from the\n"
"reference and input files which can be treated as X_1, ..., X_N Cartesian\n"
"coordinates and the method assumes both of the point sets in the\n"
"same reference frame. The point 'A' from the reference list and \n"
"the point 'P' from the input list forms a pair if the closest point\n"
"to 'A' from the input list is 'P' and vice versa.\n\n");

 longhelp_fprint(fw,grmatch_long_help,0,-1);

 fprintf(fw,"\n");
 fprintf(fw,"Report bugs to <%s>\n",FI_MAINT_EMAIL);

 return(0);
}

/*****************************************************************************/

int fprint_pointmatch_stat(FILE *fw,int nhit,int nref,int ninp,matchpointstat *mps)
{
 double	rat;
 fprintf(fw,"# Residual: %g (native)\n",mps->nsigma);
 fprintf(fw,"# Residual: %g (weighted)\n",mps->wsigma);
 fprintf(fw,"# Unitarity: %g\n",mps->unitarity);
 fprintf(fw,"# Points: %d %d %d (number of: matched, reference, input)\n",nhit,nref,ninp);
 if ( ninp<nref )	rat=(double)nhit/(double)ninp;
 else			rat=(double)nhit/(double)nref;
 fprintf(fw,"# Ratio: %.2f (percent)\n",rat*100.0);
 fprintf(fw,"# Timing: %.3f %.3f %.3f %d (total, triangle match, coordinate match: in seconds; iterations)\n",
		mps->time_total,mps->time_trimatch,mps->time_symmatch,mps->nmiter);
 fprintf(fw,"# Triangulation: %d (maximal triangulation level)\n",mps->tri_level);
 fprintf(fw,"# All: %g %g %g %d %d %d %.2f %.3f %.3f %.3f %d %d\n",
	mps->nsigma,mps->wsigma,mps->unitarity,
	nhit,nref,ninp,rat*100.0,
	mps->time_total,mps->time_trimatch,mps->time_symmatch,
	mps->nmiter,mps->tri_level);
 return(0);
}

int fprint_general_stat(FILE *fw,int nhit,int nref,int ninp)
{
 double	rat;

 if ( ninp<nref )	rat=(double)nhit/(double)ninp;
 else			rat=(double)nhit/(double)nref;
 fprintf(fw,"# Ratio: %.2f (percent)\n",rat*100.0);
 fprintf(fw,"# Points: %d %d %d (number of: matched, reference, input)\n",nhit,nref,ninp);

 return(0);
}


int fwrite_excluded_lines(FILE *fw,iline *refls,int nref,cphit *hits,int nhit,int t)
{
 char	*excflag;
 int	i,l;

 excflag=(char *)malloc(nref);
 memset(excflag,0,nref);

 if ( t )	t=1;

 for ( i=0 ; i<nhit && hits != NULL ; i++ )
  {	l=hits[i].idx[t];
	excflag[l]=1;
  }
 for ( i=0 ; i<nref ; i++ )
  {	if ( excflag[i] )	continue;
	fprintf(fw,"%s\n",refls[i].line);
  }

 free(excflag);

 return(0);
}

int main(int argc,char *argv[])
{
 FILE		*fr,*fw;
 int		i,order,hintorder,nvar,is_help;
 char		*infile,*outfile,*reffile,*outlinefile,*outidfile,
		*outxreffile,*outxinpfile,
		*inpcolids,*inpcolcoords,
		*refcolids,*refcolcoords,
		*outtransfile,*hntransfile,*intransfile,*intransparam,
		*triparam,*fitparam,*wghparam;
 colinfo 	colref,colinp;
 int		matchtype,matchread,defined_id,defined_coord,ambig;
 matchpointtune	mptp;
 matchpointstat	mps;

 iline		*refls,*inpls;
 cphit		*hits;
 int		nref,ninp,nhit;
 int		is_failed;		/* true if the matching is failed */

 double		ox,oy,scale;
 double		*xfit,*yfit,**vfits;
 transformation	itf_data,*itf=&itf_data,
		otf_data,*otf=&otf_data,
		htf_data,*htf=&htf_data;

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 infile=outfile=reffile=NULL;
 outxreffile=outxinpfile=NULL;
 outidfile=outlinefile=NULL;
 intransfile=outtransfile=intransparam=triparam=fitparam=wghparam=NULL;
 is_comment=is_verbose=is_help=0;
 order=1; ox=oy=0.0,scale=1.0; 

 hintorder=-1;
 hntransfile=NULL;
 
 colref.colord=0,colref.neg_ordering=0;
 colref.ncolid=colref.ncolcoord=0;refcolids=refcolcoords=NULL;colref.colwgh=-1;
 colref.colids=colref.colcoords=NULL;
 colinp.colord=0,colinp.neg_ordering=0;
 colinp.ncolid=colinp.ncolcoord=0;inpcolids=inpcolcoords=NULL;colinp.colwgh=-1;
 colref.colids=colinp.colcoords=NULL;

 matchtype=MATCH_POINTS; matchread=0;
 ambig=AMBIG_NONE;

 defined_id=0,defined_coord=0;

 mptp.ttype=2;
 mptp.maxdist=-1.0;mptp.use_ordering=0;
 mptp.maxnum_ref=mptp.maxnum_inp=0;
 mptp.unitarity=0.0;
 mptp.parity=0;

 mptp.nmiter=0;
 mptp.friter=0;
 mptp.rejlevel=3;

 mptp.maxcenterdist=0.0;
 mptp.is_centering=0;

 i=scanarg(argc,argv,SCANARG_ALLOW_FLAGS,
	"--version:%SN-1f%q",&is_help,
	"--version-short|--short-version:%SN-2f%q",&is_help,
	"-h|--help|--short-help|--help-short:%f%q",&is_help,
	"--long-help|--help-long:%SN2f%q",&is_help,

	"-i|--input:%s",&infile,
	"-r|--reference|--input-reference:%s",&reffile,
	"-o|--output|-om|--output-match|--output-matched:%s",&outfile,
	"--output-id|--output-ids:%s",&outidfile,
	"-L|--output-line|--output-lines:%s",&outlinefile,
	"--output-excluded-reference:%s",&outxreffile,
	"--output-excluded-input:%s",&outxinpfile,

	"-T|--input-transformation:%s",&intransfile,
	"-t|--transformation:%s",&intransparam,

	"-H|--hint-transformation:%s",&hntransfile,
	"-b|--hint-order:%d",&hintorder,

	"--match-point|--match-points:" SNf(MATCH_POINTS),&matchtype,
	"--colr|--col-ref:%f%s",&defined_coord,&refcolcoords,
	"--coli|--col-inp:%f%s",&defined_coord,&inpcolcoords,
	"--col-ref-ordering:%d",&colref.colord,
	"--col-inp-ordering:%d",&colinp.colord,

	"--center-reference:%0f%g,%g",&mptp.is_centering,&mptp.refcx,&mptp.refcy,
	"--center-input:%1f%g,%g",&mptp.is_centering,&mptp.inpcx,&mptp.inpcy,
	"--center-maxdist:%g",&mptp.maxcenterdist,

	"--triangulate|--triangulation:%s",&triparam,
	"-m|--maxdist|--max-distance:%g",&mptp.maxdist,
	"--fit:%s",&fitparam,
	"--weight:%s",&wghparam,
	"-a|--order:%d",&order,
	"-f|--offset:%g,%g",&ox,&oy,
	"-s|--scale:%g",&scale,
	"--output-transformation:%s",&outtransfile,

	"--match-coord|--match-coords:"SNf(MATCH_COORDS),&matchtype,

	"--match-id|--match-identifiers:" SNf(MATCH_IDS),&matchtype,
	"--match-ids|--match-identifier:" SNf(MATCH_IDS),&matchtype,
	"--col-ref-id:%f%s",&defined_id,&refcolids,
	"--col-inp-id:%f%s",&defined_id,&inpcolids,
	"--no-ambiguity:"	SNf(AMBIG_NONE ), &ambig,
	"--first-ambiguity:"	SNf(AMBIG_FIRST), &ambig,
	"--any-ambiguity:"	SNf(AMBIG_ANY  ), &ambig,
	"--full-ambiguity:"	SNf(AMBIG_FULL ), &ambig,

	"--comment:%i",&is_comment,"(C):%i",&is_comment,
	"--verbose:%i",&is_verbose,"(V):%i",&is_verbose,
	
	"-:%w",&infile,
	"-*|+*:%e",
	"*:%w",&infile,
	NULL);

 if ( i )		
  {	fprint_error("invalid command line argument near '%s'",argv[i]);
	return(1);
  }
 if ( is_help<0 )
  {	fprint_generic_version(stdout,argv[0],"grmatch",FI_GRMATCH_VERSION,is_help);
	return(0);
  }
 else if ( is_help>1 )
  {	fprint_grmatch_long_help(stdout);
	return(0);
  } 
 else if ( is_help )
  {	fprint_grmatch_usage(stdout);
	return(0);
  }
 
 if ( defined_id && ! defined_coord )
	matchtype=MATCH_IDS;

 if ( matchtype==MATCH_POINTS || matchtype==MATCH_COORDS )
	matchread |= MATCH_READ_COORDS;
 if ( defined_id )
	matchread |= MATCH_READ_ID;

 if ( normalize_columns(&colref,refcolcoords,refcolids) )
  {	fprint_error("invalid column specification");
	return(1);
  }
 if ( normalize_columns(&colinp,inpcolcoords,inpcolids) )
  {	fprint_error("invalid column specification");
	return(1);
  }

 if ( (matchread & MATCH_READ_COORDS) && colref.ncolcoord != colinp.ncolcoord )
  {	fprint_error("reference and input dimension mismatched");
	return(1);
  }
 if ( matchtype==MATCH_POINTS && colref.ncolcoord != 2 )
  {	fprint_error("only two-dimensional point matching is implemented");
	return(1);
  }

 if ( colref.colord>=0 && colinp.colord>=0 )	mptp.use_ordering=1; 
 else						mptp.use_ordering=0;

 if ( triparam != NULL )
  {	int	is_auto;
	int	maxnum;

	is_auto=0;
	maxnum=0;
	i=scanpar(triparam,SCANPAR_DEFAULT,
		"full:%SN-1f",&mptp.ttype,
		"delaunay:%SN0f",&mptp.ttype,
		"level:%d",&mptp.ttype,
		"auto:%f",&is_auto,
		"maxnumber:%d",&maxnum,
		"maxref|maxnumref|maxrefnum:%d",&mptp.maxnum_ref,
		"maxinp|maxnuminp|maxinpnum:%d",&mptp.maxnum_inp,
		"unitarity:%g",&mptp.unitarity,
		"mixed:%SN0f",&mptp.parity,
		"conformable:%SN1f",&mptp.parity,
		"reverse:%SN-1f",&mptp.parity,		
		NULL);
	
	if ( i )
	 {	fprint_error("invalid triangulation parameter in '%s'",triparam);
		return(1);
	 }

	if ( maxnum>0 )
	 {	mptp.maxnum_ref=maxnum;
		mptp.maxnum_inp=maxnum;
	 }

	if ( ! is_auto )		mptp.unitarity=0.0;
	else if ( mptp.unitarity<=0.0 )	mptp.unitarity=0.005;
  }

 if ( fitparam != NULL )
  {	i=scanpar(fitparam,SCANPAR_DEFAULT,
		"iterations:%d",&mptp.nmiter,
		"firstrejection:%d",&mptp.friter,
		"sigma:%g",&mptp.rejlevel,
		NULL);
	
	if ( i )
	 {	fprint_error("invalid fit parameter in '%s'",fitparam);
		return(1);
	 }
  }

 if ( wghparam != NULL )
  {	int	wcat,wcol,is_mag;
	double	power;

	wcat=0;wcol=1;
	is_mag=0;
	power=1.0;

	i=scanpar(wghparam,SCANPAR_DEFAULT,
		"reference:" SNf(0), &wcat,
		"input:"     SNf(1), &wcat,
		"column:%d",&wcol,
		"magnitude:%f",&is_mag,
		"power:%g",&power,
		NULL);
	if ( i || wcol<=0 )
	 {	fprint_error("invalid weight parameter in '%s'",wghparam);
		return(1);
	 }

	mptp.wcat=wcat;
	mptp.wpower=power;
	mptp.w_magnitude=is_mag;

	if ( ! wcat )	colref.colwgh=wcol-1;
	else		colinp.colwgh=wcol-1;
  }
 else
  {	mptp.wcat=-1;
	mptp.wpower=0.0;
	mptp.w_magnitude=0;
  }

 if ( mptp.is_centering != 3 )	mptp.is_centering=0;
 else				mptp.is_centering=1;
 if ( mptp.maxcenterdist<0.0 )	mptp.maxcenterdist=0.0;

 if ( reffile==NULL )	
  {	fprint_error("no reference file has been specified");
	return(1);
  }
 fr=fopenread(reffile);
 if ( fr==NULL )
  {	fprint_error("unable to open reference file '%s'",reffile);
	return(1);
  }
 read_match_data_points(fr,&colref,&refls,&nref,matchread);
 fcloseread(fr);
 if ( infile==NULL )	fr=stdin;
 else			fr=fopenread(infile);
 if ( fr==NULL )	
  {	fprint_error("unable to open input file '%s'",infile);
	return(1);
  }
 read_match_data_points(fr,&colinp,&inpls,&ninp,matchread);
 fcloseread(fr);

 if ( intransfile != NULL )
  {	FILE    *ft;
	ft=fopenread(intransfile);
	if ( ft==NULL ) 
	 {	fprint_error("unable to open input transformation file '%s'",intransfile);
		return(1);
	 }
	i=transformation_read_data(ft,itf);
	if ( i )
	 {	fprint_error("unable to parse transformation data in '%s'",intransfile);
		return(1);
	 }
	if ( itf->nval != 2 )
	 {	fprint_error("transformation is not a 2D -> 2D one");
		return(1);
	 }
	fcloseread(ft);
  }
 else if ( intransparam != NULL )
  {	i=transformation_parse_params(intransparam,itf);
	if ( i )		
	 {	fprint_error("unable to parse transformation string '%s'",intransparam);
		return(1);
	 }
	if ( itf->nval != 2 )
	 {	fprint_error("transformation is not a 2D -> 2D one");
		return(1);
	 }
  }
 else   itf=NULL;

 if ( hntransfile != NULL )
  {	FILE    *ft;
	ft=fopenread(hntransfile);
	if ( ft==NULL ) 
	 {	fprint_error("unable to open hint transformation file '%s'",hntransfile);
		return(1);
	 }
	i=transformation_read_data(ft,htf);
	if ( i )
	 {	fprint_error("unable to parse transformation data in '%s'",intransfile);
		return(1);
	 }
	if ( htf->nval != 2 )
	 {	fprint_error("hint transformation is not a 2D -> 2D one");
		return(1);
	 }
	fcloseread(ft);
	if ( hintorder <= 0 )
		hintorder=order;
  }
 else
	htf=NULL;

 nvar=(order+1)*(order+2)/2;
 vfits=(double **)malloc(sizeof(double *)*2);
 vfits[0]=xfit=(double *)malloc(sizeof(double)*nvar);
 vfits[1]=yfit=(double *)malloc(sizeof(double)*nvar);

 hits=NULL;
 nhit=0;

 switch ( matchtype )
  {  
     case MATCH_POINTS:

	mptp.htf=htf;
	mptp.hintorder=hintorder;

 	do_pointmatch(refls,nref,inpls,ninp,&mptp,&hits,&nhit,order,vfits,&mps);

	if ( nhit<=0 )				is_failed=1;
	else if ( ! isfinite(mps.unitarity) )	is_failed=1,mps.unitarity=-1;
	else if ( mps.unitarity<0.0 )		is_failed=1;
	else					is_failed=0;

	for ( i=0 ; i<(order+1)*(order+2)/2 ; i++ )
	 {	if ( ! isfinite(vfits[0][i]) )	is_failed=1;
		if ( ! isfinite(vfits[1][i]) )	is_failed=1;
	 }
	if ( is_failed )
	 {	for ( i=0 ; i<(order+1)*(order+2)/2 ; i++ )
		 {	vfits[0][i]=0.0;
			vfits[1][i]=0.0;
		 }
		mps.wsigma=0.0;
		mps.nsigma=0.0;
	 }
	break;

     case MATCH_COORDS:
	do_coordmatch(refls,nref,inpls,ninp,&hits,&nhit,colref.ncolcoord,mptp.maxdist);
	if ( nhit<=0 )	is_failed=1;
	else		is_failed=0;	/* almost always successful */
	break;

     case MATCH_IDS:
 	do_idmatch(refls,nref,inpls,ninp,&hits,&nhit,ambig);	
	if ( nhit<=0 )	is_failed=1;
	else		is_failed=0;	/* almost always successful */
	break;

     default:
	hits=NULL;
	nhit=0;
	is_failed=1;
	break;

  }

 if ( hits==NULL )	nhit=0;

 if ( outfile != NULL && hits != NULL && nhit>0 )
  {	fw=fopenwrite(outfile);
	if ( fw==NULL )	
	 {	fprint_error("unable to create output file '%s'",outfile);
		return(1);
	 } 

	for ( i=0 ; i<nhit ; i++ )
	 {	int	l1,l2;
		l1=hits[i].idx[0],
		l2=hits[i].idx[1];
		fprintf(fw,"%s\t%s\n",refls[l1].line,inpls[l2].line);
	 }

	if ( is_comment )
	 {	switch ( matchtype )
	 	 {   case MATCH_POINTS:
			fprint_pointmatch_stat(fw,nhit,nref,ninp,&mps);
			if ( is_failed )	fprintf(fw,"# Match failed.\n");
			break;
		     case MATCH_COORDS:
			fprint_general_stat(fw,nhit,nref,ninp);
			break;
		     case MATCH_IDS:
			fprint_general_stat(fw,nhit,nref,ninp);
			break;
		 }
	 }
	fclosewrite(fw);
  }

 if ( outxreffile != NULL )
  {	fw=fopenwrite(outxreffile);
	if ( fw==NULL )	
	 {	fprint_error("unable to create output file '%s' for excluded reference lines",outxreffile);
		return(1);
	 }
	fwrite_excluded_lines(fw,refls,nref,hits,nhit,0);
	fclosewrite(fw);
  }
 if ( outxinpfile != NULL )
  {	fw=fopenwrite(outxinpfile);
	if ( fw==NULL )	
	 {	fprint_error("unable to create output file '%s' for excluded input lines",outxinpfile);
		return(1);
	 }
	fwrite_excluded_lines(fw,inpls,ninp,hits,nhit,1);
	fclosewrite(fw);
  }
	

 if ( outidfile != NULL && hits != NULL && nhit>0 )
  {	int	mxr,mxi,l;
	char	sbuff[32],*ir,*ii;

	fw=fopenwrite(outidfile);
	if ( fw==NULL )	
	 {	fprint_error("unable to create output identifier list file '%s'",outidfile);
		return(1);
	 }

	mxr=mxi=1;
	for ( i=0 ; i<nhit ; i++ )
	 {	int	l1,l2;
		l1=hits[i].idx[0],
		l2=hits[i].idx[1];
		if ( refls[l1].id != NULL )
		 {	l=strlen(refls[l1].id);if ( l>mxr ) mxr=l;	}
		if ( inpls[l2].id != NULL )
		 {	l=strlen(inpls[l2].id);if ( l>mxi ) mxi=l;	}
	 }
	sprintf(sbuff,"%%%ds\t%%%ds\n",mxr,mxi);
	for ( i=0 ; i<nhit ; i++ )
	 {	int	l1,l2,j;
		l1=hits[i].idx[0],
		l2=hits[i].idx[1];
		if ( refls[l1].id != NULL )	ir=refls[l1].id;
		else				ir="-";
		if ( inpls[l2].id != NULL )	ii=inpls[l2].id;
		else				ii="-";
		for ( j=0 ; ir[j] ; j++ )  { if ( ir[j]=='\n' ) ir[j]=' '; }
		for ( j=0 ; ii[j] ; j++ )  { if ( ii[j]=='\n' ) ii[j]=' '; }
		fprintf(fw,sbuff,ir,ii);
	 }
	fclosewrite(fw);
  }

 if ( outtransfile != NULL && matchtype==MATCH_POINTS )
  {	
	fw=fopenwrite(outtransfile);
	if ( fw==NULL )	
	 {	fprint_error("unable to create output transformation file '%s'",outidfile);
		return(1);
	 }

	if ( is_comment )
	 {	fprintf(fw,"# Created by grmatch %s (fi: %s)\n",FI_GRMATCH_VERSION,FI_VERSION);
		fprintf(fw,"# Invoked command:");
		for ( i=0 ; i<argc ; i++ )
		 {	if ( is_any_nasty_char(argv[i]) )
				fprintf(fw," \"%s\"",argv[i]);
			else
				fprintf(fw," %s",argv[i]);
		 }
		fprintf(fw,"\n");
	 }

	otf->type=TRANS_POLYNOMIAL;
	otf->order=order;
	otf->ox=otf->oy=0.0,otf->scale=1.0;
	otf->nval=2;
	otf->vfits=vfits;
	otf->bshx=otf->bshy=0.0;

	transformation_write_data(fw,otf,(is_comment?TRANS_WR_COMMENT:0)|TRANS_WR_DXDY);

	if ( is_comment )
	 {	fprint_pointmatch_stat(fw,nhit,nref,ninp,&mps);
		if ( is_failed )	fprintf(fw,"# Match failed.\n");
	 }

	fclosewrite(fw);
  } 

 /* if matching is failed, errorlevel will be 2: we follow DOS conventions,  */
 /* where higher return code represent less-critical error (invalid command  */
 /* line argument, invalid input/output file or such critical errors will    */
 /* return 1 as exit status).						     */
 if ( is_failed )	return(2);
 else			return(0);

}
