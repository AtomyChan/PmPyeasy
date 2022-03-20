/***********************b*****************************************************/
/* trimatch.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Triangle match of 2D pointsets. Originally written by Istvan Domsa (see   */
/* copyright below), minor modifications (to be 'fi'-compatible) and some    */
/* extensions (see expand_triangulation(), new triangle space coordinates    */
/* and related functions) were made by A. Pal (apal@szofi.elte.hu).	     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2003; Domsa, I. (domsa@konkoly.hu); 2005-06, Pal, A. (apal@szofi)     */
/*****************************************************************************/

/* trimatch.c : triangle match of 2D pointsets (domsa@konkoly.hu (2002)) */
/* code based on Michael Richmond's "match" program (PASP, 107, 1119) */
/* written for the HAT Project */
/* trimatch.c,v 5.5 2003/05/13 16:05:26 domsa Exp */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "fit/lmfit.h"
#include "poly.h"
#include "polyfit.h"

#include "spmatrix.h"
#include "cpmatch.h"
#include "tpoint.h"
#include "delaunay.h"

#include "trimatch.h"

#define MINVOTE			1 	/* pairs get accounted with larger no. of votes */
#define MATCH_PERCENTILE	0.4 	/* best 70% used for transformation calc */
#define MATCH_MINSIGMA		1e-06	/* stop iterating if this sigma reached */
#define MATCH_BIGNSIGMA		36.0	/* no pairs above 6sigma distance */
#define MATCH_MAXITER		2	/* how many times to iterate to find trans */

#define MAX(x,y)      		((x) > (y) ? (x) : (y))

typedef struct
 {	tpoint *p1;		/* ptr to "obj" */
	tpoint *p2;		/* ptr to "cat" */
	union
	 {	double  votes;	/* confidence of matched pairs */
		double  dist;
	 } confidence;
 } pair;

/* turns triangles into points in "triangle space" */
tpoint *tri2point(int len, triangle *tris)
{
 int	ii;
 tpoint	*points,*curr;

 /* propagate internal error */
 if ( len < 1 || tris == NULL )
	return(NULL);
 if ( (points=malloc(len*sizeof(*points))) == NULL )
	return(NULL);

 for ( ii=0, curr=points ; ii<len ; ii++,curr++ )
  {	curr->id    =ii;
	curr->xcoord=tris[ii].trix;	/* "x" is the b/a ratio */
	curr->ycoord=tris[ii].triy;	/* "y" is the c/a ratio */
  }

 return(points);
}

static int triangle_sidesquares(triangle *t,double *ret)
{
 ret[0]=(t->vertices[1]->xcoord-t->vertices[2]->xcoord)*
	(t->vertices[1]->xcoord-t->vertices[2]->xcoord)+
 	(t->vertices[1]->ycoord-t->vertices[2]->ycoord)*
	(t->vertices[1]->ycoord-t->vertices[2]->ycoord);
 ret[1]=(t->vertices[2]->xcoord-t->vertices[0]->xcoord)*
	(t->vertices[2]->xcoord-t->vertices[0]->xcoord)+
 	(t->vertices[2]->ycoord-t->vertices[0]->ycoord)*
	(t->vertices[2]->ycoord-t->vertices[0]->ycoord);
 ret[2]=(t->vertices[0]->xcoord-t->vertices[1]->xcoord)*
	(t->vertices[0]->xcoord-t->vertices[1]->xcoord)+
 	(t->vertices[0]->ycoord-t->vertices[1]->ycoord)*
	(t->vertices[0]->ycoord-t->vertices[1]->ycoord);
 return(0);
}

static int triangle_side_sort(double *s,int *i)
{
 if ( s[0] <= s[1] && s[0] <= s[2] )
  {	if ( s[1] <= s[2] )
		i[0]=0,i[1]=1,i[2]=2;
	else
		i[0]=0,i[1]=2,i[2]=1;
  }
 else if ( s[1] <= s[0] && s[1] <= s[2] )
  {	if ( s[0] <= s[2] )
		i[0]=1,i[1]=0,i[2]=2;
	else
		i[0]=1,i[1]=2,i[2]=0;
  }
 else
  {	if ( s[0] <= s[1] )
		i[0]=2,i[1]=0,i[2]=1;
	else
		i[0]=2,i[1]=1,i[2]=0;
  }
 return(0);
}

/* create a match matrix */
spmatrix *getmatchmatrix(int nobjtri,triangle *objtri,int ncattri,triangle *cattri)
{
 int		ii,vote;
 int		hlen;
 spmatrix	*votematrix;
 tpoint		*catpoints, *objpoints;
 tpointarr	*catarr, *objarr;
 cphit		*hits;

 catpoints=tri2point(ncattri,cattri);
 objpoints=tri2point(nobjtri,objtri);
 if ( catpoints==NULL || objpoints==NULL )
  {	if ( catpoints != NULL )	free(catpoints);
	if ( objpoints != NULL )	free(objpoints);
	return(NULL);
  }

 catarr=tpoint_buildarr(ncattri,catpoints);
 objarr=tpoint_buildarr(nobjtri,objpoints); 
 if ( catarr==NULL || objarr==NULL )
  {	if ( catarr != NULL )	tpoint_destroyarr(catarr);
	if ( objarr != NULL )	tpoint_destroyarr(objarr);
	return(0);
  }

 if ( (votematrix= spm_create()) == NULL )
  {	tpoint_destroyarr(catarr);
	tpoint_destroyarr(objarr);
	return(NULL);
  }

 
/* 
 maxdist=1.0;
 if ( (hits=cpmatch(catarr,objarr,&hlen))==NULL ||
      (hlen=cphit_revise(hlen,&hits,MAX(ncattri,nobjtri),CPMTYPE_ACC,&maxdist)) == 0 )
*/

 /* try to match as many pairs as possible */
 if ( (hits=cpmatch_symmetric(catarr,objarr,&hlen,4.0,0.0)) == NULL )
  {	tpoint_destroyarr(catarr);
	tpoint_destroyarr(objarr);
	spm_destroy(votematrix);
	return(NULL);
  }

 qsort(hits,hlen,sizeof(cphit),cphit_sort_dist);

 /* fill up vote matrix */
 for ( ii=0 ; ii<hlen ; ii++ )
  {	
	double	sobj[3],scat[3];
	int	iobj[3],icat[3];
	int	xc,xo;

	vote=hlen-ii;
	xc=hits[ii].idx[0];
	xo=hits[ii].idx[1];

	triangle_sidesquares(&objtri[xo],sobj);
	triangle_side_sort(sobj,iobj);
	triangle_sidesquares(&cattri[xc],scat);
	triangle_side_sort(scat,icat);

	/* fprintf(stderr,"[%d%d%d][%d%d%d]\n",iobj[0],iobj[1],iobj[2],icat[0],icat[1],icat[2]); */

	spm_addval(votematrix,vote,objtri[xo].vertices[iobj[0]]->id,cattri[xc].vertices[icat[0]]->id);
	spm_addval(votematrix,vote,objtri[xo].vertices[iobj[1]]->id,cattri[xc].vertices[icat[1]]->id);
	spm_addval(votematrix,vote,objtri[xo].vertices[iobj[2]]->id,cattri[xc].vertices[icat[2]]->id);
  }

 /* release memory */
 tpoint_destroyarr(objarr);
 tpoint_destroyarr(catarr);

 return(votematrix);
}

/* get all non-zero elements from a vote matrix */
pair *get_pairs(spmatrix *votematrix,tpoint *obj,tpoint *cat,int *plen)
{
 int	npairs;
 pair	*pairs, *cp;
 spnode	node,witness;

 /* just to be sure */
 if ( votematrix == NULL || obj == NULL || cat == NULL )
	return(NULL);
 if ( (pairs=malloc(spm_numofelems(votematrix)*sizeof(*pairs))) == NULL )
	return(NULL);

 npairs=0;witness=NULL;
 node=spm_getnextnode(votematrix,NULL,&witness);
 while ( node != NULL )
  {	cp=pairs+npairs;
	cp->p1=obj+spm_getrow(node);
	cp->p2=cat+spm_getcol(node);
	cp->confidence.votes=spm_getnodeval(node);
	npairs++;
	node=spm_getnextnode(votematrix,node,&witness);
  }

 *plen=npairs;
 return(pairs);
}

/* sort a list of pairs according to votes (for qsort()) */
int sort_pairs(const void *p1,const void *p2)
{
 int vote1=((pair *)p1)->confidence.votes;
 int vote2=((pair *)p2)->confidence.votes;
 return(vote2-vote1);
}

/* sort a list of pairs according to distances (for qsort()) */
int sort_pairsdist(const void *p1, const void *p2)
{
 double dist1=((pair *)p1)->confidence.dist;
 double dist2=((pair *)p2)->confidence.dist;

 return(dist2>dist1 ? -1 : 1);
}

/* return the "perc"th percentile distance */
double find_percdist(int len, pair *pairs, double perc)
{
 int index=(int)(len*perc+0.5);

 if ( index >= len )
	index=len-1;
 return(pairs[index].confidence.dist);
}

/* that's the main stuff: get the linear transform from found pairs */
int calc_trans(int len,pair *pairs,int order,double *xfit,double *yfit,double *sigma)
{
 int    	ii,k,nvar;
 point		*fpoints;

 nvar=(order+1)*(order+2)/2;
 for ( ii=0 ; ii<nvar ; ii++ )
  {	xfit[ii]=0.0;
	yfit[ii]=0.0;
  }

 /* quit if we don't have enough points */
 if ( len<3 || pairs == NULL )	return(1);

 /* get transformation */
 fpoints=(point *)malloc(len*sizeof(point));
 for ( ii=0 ; ii<len ; ii++ )
  {	fpoints[ii].x=pairs[ii].p2->xcoord;
	fpoints[ii].y=pairs[ii].p2->ycoord;
	fpoints[ii].weight=1.0;
  }

 for ( ii=0 ; ii<len ; ii++ )
  {	fpoints[ii].value=pairs[ii].p1->xcoord;		}
 k=fit_2d_poly(fpoints,len,order,xfit,0,0,1);
 if ( ! k )
  {	for ( ii=0 ; ii<len ; ii++ )
	 {	fpoints[ii].value=pairs[ii].p1->ycoord;		}
	k=fit_2d_poly(fpoints,len,order,yfit,0,0,1);
  }
 free(fpoints);
	
 return(k);
}

/* iterate until the best matching transformation found */
int get_trans(int *len, pair *pairs, int order,double *xfit,double *yfit,
	      double *sumsig, double maxdist)
{
 int		length,ii,iter,badnumber;
 double		sigma,bigsigma,nx,ny,dx,dy;

 /* propagate internal error */
 if ( *len < 3 || pairs == NULL )	return(1);

 if ( maxdist > 0 )		maxdist=maxdist*maxdist;
 length=*len;

 /* start iteration */
 iter=0;
 do
  {	/* calculate transformation */
	if ( calc_trans(length,pairs,order,xfit,yfit,sumsig) )
		return(1);

	/* calculate euclidian distances of transformed points to their pairs */
	for ( ii=0; ii<length; ii++ )
	 {	nx=eval_2d_poly(pairs[ii].p2->xcoord,pairs[ii].p2->ycoord,order,xfit,0,0,1);
		ny=eval_2d_poly(pairs[ii].p2->xcoord,pairs[ii].p2->ycoord,order,yfit,0,0,1);
		dx=pairs[ii].p1->xcoord-nx,
		dy=pairs[ii].p1->ycoord-ny;
 		pairs[ii].confidence.dist=dx*dx+dy*dy;
	 }

	/* sort pairs along distance */
	qsort((void *)pairs,length,sizeof(*pairs),sort_pairsdist);

	/* no bad point so far */
	badnumber=0;

	/* get rid of pairs which are too separated */
	if ( maxdist > 0.0 )
	 {	for ( ii=0 ; ii<length ; ii++ )
		 {	if ( pairs[ii].confidence.dist > maxdist )
			 {	badnumber = 1;
				length    = ii;
				break;
			 }
		 }
	 }

	/* check if there are any points left */
	if ( length==0 )
		return(1);

	/* get sigma for MATCH_PERCENTILE best matches */
	sigma=find_percdist(length,pairs,MATCH_PERCENTILE);

	/* if sigma is too small then stop now */
	if ( sigma<MATCH_MINSIGMA )
		break;

	/* throw out points with too big sigmas */
	bigsigma=MATCH_BIGNSIGMA*sigma;

	if ( maxdist < 0.0 || (maxdist > 0.0 && bigsigma < maxdist) )
	 {	for ( ii=0 ; ii<length ; ii++ )
		 {	if ( pairs[ii].confidence.dist > bigsigma )
			 {	badnumber=1;
				length   =ii;
				break;
			 }
		 }
		if ( length==0 )	return(1);
	 }

  } while ( badnumber && ++iter < MATCH_MAXITER );

 *len=length;

 return(calc_trans(length,pairs,order,xfit,yfit,sumsig));
}

/* get those pairs with highest votes (ret 0 on error) */
int get_bigvotes(pair *matches, int len)
{
 int	ii,jj,min_vote,already;
 tpoint	*cur_obj,*cur_cat;
 tpoint	**already_obj,**already_cat;

 if ( matches==NULL || len<1 )
	return(0);

 already_obj=malloc(len*sizeof(tpoint *));
 already_cat=malloc(len*sizeof(tpoint *));
 if ( already_obj==NULL || already_cat==NULL )
  {	if ( already_obj != NULL )	free(already_obj);
	if ( already_cat != NULL )	free(already_cat);
	return(0);
  }
 memset(already_obj,0,len*sizeof(tpoint *));
 memset(already_cat,0,len*sizeof(tpoint *));

 min_vote=MINVOTE;
 already =0;

 for ( ii=0 ; ii<len ; ii++ )
  {	if ( matches[ii].confidence.votes<min_vote )
		break;

	cur_obj=matches[ii].p1;
	cur_cat=matches[ii].p2;

	for ( jj=0 ; jj<already ; jj++ )
	 {	if ( cur_obj==already_obj[jj] )
		        break;
	 }

	if ( jj != already )
		break;

	for ( jj=0 ; jj<already ; jj++ )
	 {	if ( cur_cat==already_cat[jj] )
			break;
	 }

	if ( jj != already )
		break;

	already_obj[already]=cur_obj;
	already_cat[already]=cur_cat;
	already++;
  }

 free(already_cat);
 free(already_obj);

 return(already);
}

/*****************************************************************************/

int triangle_sort(triangle *tt,int parity)
{
 double	x1,y1,x2,y2;
 tpoint	*p0,*p1,*p2,*pt;
 
 p0=tt->vertices[0];
 p1=tt->vertices[1];
 p2=tt->vertices[2];

 x1=p1->xcoord-p0->xcoord,y1=p1->ycoord-p0->ycoord;
 x2=p2->xcoord-p0->xcoord,y2=p2->ycoord-p0->ycoord;
 if ( (double)parity*(x1*y2-y1*x2) < 0.0 )
 	pt=tt->vertices[2],tt->vertices[2]=tt->vertices[1],tt->vertices[1]=pt;

 return(0);
}

int triangle_sort_ccw(triangle *tt)
{
 return ( triangle_sort(tt,+1) );
}
int triangle_sort_cw(triangle *tt)
{
 return ( triangle_sort(tt,-1) );
}

int triangle_space_coords_mixed(triangle *tt)
{
 tpoint	*p0,*p1,*p2,*sorted[3];
 double	sides[3],ab,ac;

 p0=tt->vertices[0];
 p1=tt->vertices[1];
 p2=tt->vertices[2];
 sides[0]=tpoint_eucdist(p1,p2),sorted[0]=p0;
 sides[1]=tpoint_eucdist(p2,p0),sorted[1]=p1;
 sides[2]=tpoint_eucdist(p0,p1),sorted[2]=p2;
 vertice_sort(sides,sorted);
 ab=sides[1]/sides[0];
 ac=sides[2]/sides[0];
 tt->trix=ab;
 tt->triy=ac;

 tt->vertices[0]=sorted[0],
 tt->vertices[1]=sorted[1],
 tt->vertices[2]=sorted[2];

 triangle_sort_ccw(tt);

 return(0);
}


int triangle_space_coords_continous(triangle *tt,int parity)
{
 tpoint	*p0,*p1,*p2;
 double	s0,s1,s2,alpha,beta,n2,x,y,r,px,py,ir3;

 if ( parity>0 )	triangle_sort_ccw(tt);
 else			triangle_sort_cw(tt);

 p0=tt->vertices[0];
 p1=tt->vertices[1];
 p2=tt->vertices[2];

 s0=tpoint_eucdist(p1,p2);
 s1=tpoint_eucdist(p2,p0);
 s2=tpoint_eucdist(p0,p1);

 if ( s0>=s1 && s0>=s2 )
  {	alpha=1.0-s1/s0,
	beta =1.0-s2/s0;
  }
 else if ( s1>=s0 && s1>=s2 )
  {	alpha=1.0-s2/s1,
	beta =1.0-s0/s1;
  }
 else
  {	alpha=1.0-s0/s2,
	beta =1.0-s1/s2;
  }

 n2=sqrt(alpha*alpha+beta*beta);
 if ( n2 <= 0.0 )
  {	tt->trix=0.0;
	tt->triy=0.0;
  }
 else
  {	x=alpha*(alpha+beta)/n2;
	y=beta *(alpha+beta)/n2;
	r=alpha+beta;
	ir3=1.0/(r*r*r);
	px=x*x-y*y;
	py=2*x*y;
	tt->trix=ir3*(px*px-py*py);
	tt->triy=ir3*(2*px*py);
  }

 return(0);
}

 
int set_triangle_space_coords(triangle *tris,int ntri,int parity)
{
 if ( ! parity )
  {	while ( ntri>0 )
	 {	triangle_space_coords_mixed(tris);
		tris++;ntri--;
	 };
  }
 else if ( parity>0 )
  {	while ( ntri>0 )
	 {	triangle_space_coords_continous(tris,+1);
		tris++;ntri--;
	 };
  }
 else
  {	while ( ntri>0 )
	 {	triangle_space_coords_continous(tris,-1);
		tris++;ntri--;
	 };
  }
 return(0);
}

/*****************************************************************************/

triangle *full_triangulation(tpoint *points,int len,int *ntri)
{
 int		numtri;
 tpoint		*p1,*p2,*p3,*e1,*e2,*e3;
 triangle	*triangles,*tt;

 *ntri=0;
 if ( len<3 || points==NULL )	return(NULL);

 /* set boundaries */
 e1=points+len-2;
 e2=points+len-1;
 e3=points+len;

 numtri=(len)*(len-1)*(len-2)/6;
 triangles=malloc(numtri*sizeof(triangle));
 if ( triangles==NULL )		return(NULL);

 for ( p1=points,tt=triangles ; p1<e1 ; p1++ )
  {	for ( p2=p1+1 ; p2<e2 ; p2++ )
	 {	for ( p3=p2+1 ; p3<e3 ; p3++,tt++ )
		 {	tt->vertices[0]=p1,
			tt->vertices[1]=p2,
			tt->vertices[2]=p3;
			tt->trix=tt->triy=0.0;
		 }
	 }
  }

 if ( ntri != NULL )	*ntri=numtri;
 return(triangles);
}

int copy_triangulation(triangle *tris,int ntri,triangle **rtris,int *rntri)
{
 if ( rtris==NULL || rntri==NULL )	return(1);
 *rtris=(triangle *)malloc(sizeof(triangle)*ntri);
 memcpy(*rtris,tris,sizeof(triangle)*ntri);
 *rntri=ntri;
 return(0);
}

int trimatch_int(tpoint *cat,int catlen,tpoint *obj,int objlen,
	triangle *cattri,int ncattri,triangle *objtri,int nobjtri,
	int order,double *xfit,double *yfit,double *rsigma,double maxdist,
	int parity,
	pair **rmatches,int *rnbigs)
{
 spmatrix	*matchmatrix;
 pair		*matches;
 int		nbigs,plen,ret;

 matchmatrix=NULL;
 matches=NULL;

 if ( ! parity )
  {	set_triangle_space_coords(cattri,ncattri,0);
	set_triangle_space_coords(objtri,nobjtri,0);
  }
 else
  {	set_triangle_space_coords(cattri,ncattri,1);
	set_triangle_space_coords(objtri,nobjtri,parity);
  }

 /* get match matrix */
 matchmatrix=getmatchmatrix(nobjtri,objtri,ncattri,cattri);

 if ( matchmatrix == NULL )	return(1);

 /* get matches */
 matches=get_pairs(matchmatrix,obj,cat,&plen);
 spm_destroy(matchmatrix);

 if ( matches == NULL )		return(1);

 /* sort pairs according to votes */
 qsort((void *)matches,plen,sizeof(*matches),sort_pairs);

 /* get pairs with biggest number of votes */
 nbigs=get_bigvotes(matches,plen);

 /* calculate transformation between catalog and image */
 ret=get_trans(&nbigs,matches,order,xfit,yfit,rsigma,maxdist);

 if ( ret )
  {	free(matches);
	if ( rmatches != NULL )	*rmatches=NULL;
	if ( rnbigs != NULL )	*rnbigs=0;
	return(1);
  }
 else
  {	if ( rmatches != NULL )	*rmatches=matches;
	if ( rnbigs != NULL )	*rnbigs=nbigs;
	return(0);
  }

}

/* trimatch_unitarity(): deprecated, use calc_2d_unitarity() instead. */

/*
double trimatch_unitarity(double *xfit,double *yfit)
{
 double	ma,mb,mc,md;
 double	n1,n2,nn,dd;
 double	unitarity;

 ma=xfit[1],mb=xfit[2];
 mc=yfit[1],md=yfit[2];
 n1=(ma-md)*(ma-md)+(mb+mc)*(mb+mc);
 n2=(ma+md)*(ma+md)+(mb-mc)*(mb-mc);
 nn=(n1<n2?n1:n2);
 dd=ma*ma+mb*mb+mc*mc+md*md;
 if ( dd<=0.0 )
  {	return(-1.0);		}
 else
  {	if ( nn<=0.0 )	unitarity=0.0;
	else		unitarity=sqrt(nn/dd);
	return(unitarity);
  }
}
*/

/* make a triangulate match of pointsets: return found transformation */
int trimatch(tpoint *cat,int catlen,tpoint *obj,int objlen,
	int order,trimatchpar *itmp,
        cphit **rhits, int *rnhit,double *xfit,double *yfit,
	trimatchlog *tml)
{
 int		nobjtri,ncattri,ncatdeltri,nobjdeltri;
 int		nbigs,ret,nhit,level,level_init,level_used;
 triangle	*cattris,*objtris,*catdeltris,*objdeltris;
 triinfo	cati,obji;
 cphit		*hits;
 trimatchpar	stmp,*tmp;
 pair		*matches=NULL;
 double		unitarity,sigma;

 /* just to be sure */
 if ( catlen <= 3 || cat == NULL || objlen <= 3 || obj == NULL )
	return(1);

 /* check order: transformation should be at least linear... */
 if ( order<1 )
	return(1);

 /* set return values */
 hits=NULL;
 nhit=0;

 if ( itmp != NULL )		/* inherit external fine-tune parameters*/
  {	tmp=itmp;	}
 else				/* Default fine-tune parameters are:	*/
  {	stmp.level=0;		/*  - Delaunay-triangulation		*/
	stmp.maxdist=-1.0;	/*  - no maximal distance constraint	*/
				/*    ^^^^^^ => this is f*** unclear    */
	stmp.unitarity=0.0;	/*  - no automatic unitarity usage	*/
	stmp.parity=0;		/*  - mixed left/right-handed triangles */
	tmp=&stmp;
  }

 unitarity=-1.0;

 if ( tmp->unitarity <= 0.0 )
  {	if ( tmp->level < 0 )	/* full triangulation */
	 {	
		objtris=full_triangulation(obj,objlen,&nobjtri);
		if ( objtris==NULL )	return(1);
		cattris=full_triangulation(cat,catlen,&ncattri);
		if ( cattris==NULL )	{ free(objtris);return(1); }
	
		/* no delaunay triangulation   */
		catdeltris=NULL,ncatdeltri=0;	      
		objdeltris=NULL,nobjdeltri=0;

		/* no extra triangulation info */
		cati.neigs=NULL,cati.neigpoints=NULL; 
		obji.neigs=NULL,obji.neigpoints=NULL;
		level_init=-1;
		level_used=-1;
	 }
	else
	 {	int	ro,rc;
		ro=delaunay_triangulation(obj,objlen,&objdeltris,&nobjdeltri,&obji);
		if ( ro )	return(1);
		rc=delaunay_triangulation(cat,catlen,&catdeltris,&ncatdeltri,&cati);
		if ( rc )
		 {	free(objdeltris);
			free(obji.neigs);free(obji.neigpoints);
			return(1);
	 	 }
		if ( tmp->level>0 )
		 {	expand_triangulation(&cati,&cattris,&ncattri,cat,catlen,tmp->level);
			expand_triangulation(&obji,&objtris,&nobjtri,obj,objlen,tmp->level);
		 }
		else
		 {	copy_triangulation(catdeltris,ncatdeltri,&cattris,&ncattri);
			copy_triangulation(objdeltris,nobjdeltri,&objtris,&nobjtri);
		 }
		level_init=tmp->level;
		level_used=tmp->level;

	 }
	ret=trimatch_int(cat,catlen,obj,objlen,cattris,ncattri,objtris,nobjtri,
		order,xfit,yfit,&sigma,tmp->maxdist,tmp->parity,&matches,&nbigs);

	free(objtris);
	free(cattris);
  }
 else
  {	int	ro,rc;
	ro=delaunay_triangulation(obj,objlen,&objdeltris,&nobjdeltri,&obji);
	if ( ro )	return(1);
	rc=delaunay_triangulation(cat,catlen,&catdeltris,&ncatdeltri,&cati);
	if ( rc )
	 {	free(objdeltris);
		free(obji.neigs);free(obji.neigpoints);
		return(1);
	 }
	matches=NULL;
	ret=trimatch_int(cat,catlen,obj,objlen,catdeltris,ncatdeltri,objdeltris,nobjdeltri,
		order,xfit,yfit,&sigma,tmp->maxdist,tmp->parity,&matches,&nbigs);

	/*fprintf(stderr,"ret_0=%d\n",ret);*/
	for ( level=1 ; level<=4 ; level++ )
	 {	unitarity=calc_2d_unitarity(xfit,yfit,order);
		/*fprintf(stderr,"unitarity: %g\n",unitarity);*/
		if ( unitarity < 0.0 )			ret=1;
		else if ( unitarity > tmp->unitarity )	ret=1;

		if ( ! ret )	break;
		if ( matches != NULL )
		 {	free(matches);
			matches=NULL;
		 }

		expand_triangulation(&cati,&cattris,&ncattri,cat,catlen,level);
		expand_triangulation(&obji,&objtris,&nobjtri,obj,objlen,level);
		ret=trimatch_int(cat,catlen,obj,objlen,cattris,ncattri,objtris,nobjtri,
			order,xfit,yfit,&sigma,tmp->maxdist,tmp->parity,&matches,&nbigs);
		free(objtris);
		free(cattris);
	 }
	level_init=0;
	level_used=level-1;
  }

 if ( ! ret && matches != NULL )
  {	if ( (hits=malloc(sizeof(cphit)*nbigs)) != NULL )
	 {	int  ii;
		for ( ii=0 ; ii<nbigs ; ii++ )
		 {	hits[ii].idx[0]=matches[ii].p2->id;
			hits[ii].idx[1]=matches[ii].p1->id;
		 }
	 }
	nhit=nbigs;
  }
 else	hits=NULL,nhit=0;

 if ( matches != NULL )		free(matches);

 if ( catdeltris != NULL )	free(catdeltris);
 if ( cati.neigs != NULL )	free(cati.neigs);
 if ( cati.neigpoints != NULL )	free(cati.neigpoints);
 if ( objdeltris != NULL )	free(objdeltris);
 if ( obji.neigs != NULL )	free(obji.neigs);
 if ( obji.neigpoints != NULL )	free(obji.neigpoints);

 if ( rhits != NULL )	*rhits=hits;
 if ( rnhit != NULL )	*rnhit=nhit;

 if ( tml != NULL )
  {	tml->residual  =sigma;
	tml->level_init=level_init;
	tml->level_used=level_used;
	if ( unitarity < 0.0 )
	 {	unitarity=calc_2d_unitarity(xfit,yfit,order);		}
	tml->unitarity=unitarity;

  }
 
 return(ret);
}

/*****************************************************************************/
                
