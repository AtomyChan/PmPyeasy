/*****************************************************************************/
/* cpmatch.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions related to point matching: 				     */
/*   -	matching between two 2d point sets;				     */
/*   -	searching for nearest point in a given point list;		     */
/*   -	searching for neighbouring points in a given point list.	     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Written by Istvan Domsa (see original copyright below) and 		     */
/* A. Pal (apal@szofi.elte.hu). (c) 2003, 2005-2006			     */
/*****************************************************************************/

/* cpmatch.c : match two 2D point array (domsa@konkoly.hu (2003)) */
/* cpmatch.c,v 5.5 2003/05/13 16:09:43 domsa Exp */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "tpoint.h"
#include "cpmatch.h"

#define	FLT_EPSILON	1e-10

#define SWAP(a,b,type) do {type t;t=a;a=b;b=t;} while (0)

/*****************************************************************************/

/* function for qsort() interface: sort hit array along distance */
int  cphit_sort_dist(const void *p1, const void *p2)
{
 double	d1=((cphit *)p1)->distance;
 double	d2=((cphit *)p2)->distance;
 return(d1>d2?1:-1);
}

/*****************************************************************************/

/* return the index of the point in the base array, whose x coord is
   just below or equals the xcoord of the probe point. base array should be
   sorted along its x coordinates. return a negative number on error */
int cpmatch_find_xidx(tpoint *probe,tpointarr *base)
{
 int	min,max,mid;
 double	xcoord,dist;
 tpoint	*arr;

 mid   = -1;    /* set return value to something invalid */

 min  =0;
 max  =base->length;
 arr  =base->points;
 xcoord=probe->xcoord;

 while ( max > min )
  {	mid = (max + min)/2;
	dist = xcoord - arr[mid].xcoord;

	if ( fabs(dist) < FLT_EPSILON )
		break;
	else if ( dist > 0.0 )
		min = mid + 1;
	else
		max = mid;
  }
 return(mid);
}

/* find the closest point in a direction. return negative number on error */
int cpmatch_find_nearest(tpoint *probe,tpointarr *base,int neighbour,
	double *rdistance,int is_exclude_self)
{
 int	posl,posr,best;
 int	length;
 double	xc,yc,xx,yy;
 double	xdist,ydist;
 double	xcoord,ycoord;
 double	maxdist,cdist;
 double	mindist,dist;
 tpoint	*basearr,*nbpoint;
 int	outl,outr;

 /* check if neighbour point is inside the array */
 length  = base->length;
 if ( neighbour < 0 || neighbour >= length )
	return(-1);

 /* x/y coordantes of the probe point */
 xcoord   = probe->xcoord;
 ycoord   = probe->ycoord;

 /* array of the base points */
 basearr = base->points;

 /* the point to start the search with: */
 /* if the index with 'neighbour' is not excluded... */
 if ( ! is_exclude_self )
  {	/* x/y distances from the probe and its neighbouring point */
	nbpoint = basearr + neighbour;
	xdist = fabs(nbpoint->xcoord-xcoord);
	ydist = fabs(nbpoint->ycoord-ycoord);

	/* calculate the initial squared minimal distance, and set best poss */
	best   =neighbour;
	mindist=xdist*xdist+ydist*ydist;
	maxdist=xdist+ydist;
  }
 /* ... or if the index with 'neighbour' is excluded: */
 else
  {	xdist=ydist=maxdist=mindist=0.0;
	if ( neighbour>0 )
	 {	nbpoint=basearr+neighbour-1;
		xdist=fabs(nbpoint->xcoord-xcoord);
		ydist=fabs(nbpoint->ycoord-ycoord);
		best   =neighbour-1;
		mindist=xdist*xdist+ydist*ydist;
		maxdist=xdist+ydist;
	 }
	else
	 {	best=-1;			}

	if ( neighbour<length-1 )
	 {	nbpoint=basearr+neighbour+1;
		xdist=fabs(nbpoint->xcoord-xcoord);
		ydist=fabs(nbpoint->ycoord-ycoord);
		if ( best<0 || xdist*xdist+ydist*ydist<mindist )
		 {	best   =neighbour+1;
			mindist=xdist*xdist+ydist*ydist;
			maxdist=xdist+ydist;
		 }
	 }

	if ( best<0 )	return(-1);
  }

 /* set maximal distance in the x axis to search for */

 outl = 1;
 outr = 1;

 for ( posl=neighbour-1,posr=neighbour+1 ; outl || outr ; --posl,++posr )
  {	if ( outl )
	 {	if ( posl >= 0 )
		 {	/* the nearest neighbour candidate on the left */
			nbpoint = basearr + posl;
			/* its coordinates */
			xc = nbpoint->xcoord;
			yc = nbpoint->ycoord;

			/* break out from loop if it's too far away */
			xx = xcoord - xc;
			if ( xx < maxdist )
			 {	yy = yc - ycoord;
				yy = fabs(yy);
				/* calculate distance only, if ycoord is less than before */
				if ( yy <= ydist )
				 {	ydist = yy;
					dist = xx*xx + yy*yy;
					if ( dist <= mindist )
					 {	/* a better match found */
						best    = posl;
						mindist = dist;
						/* xx should be positive */
						xx=fabs(xx);
						/* reset limit */
						cdist=xx+yy;
						if ( cdist < maxdist )
							maxdist = cdist;
					 }
				 }
			 }
			else
				outl=0;
		 }
		else
			outl=0;
	 }
	if ( outr )
	 {	if ( posr < length )
		 {	/* the nearest neighbour candidate on the right */
			nbpoint = basearr + posr;
			/* and its coordinates */
			xc = nbpoint->xcoord;
			yc = nbpoint->ycoord;
			/* break out from loop if it's too far away */
			xx = xc - xcoord;
			if ( xx < maxdist )
			 {	yy = yc - ycoord;
				yy = fabs(yy);
				/* calculate distance only, if ycoord is less than before */
				if ( yy <= ydist )
				 {	ydist = yy;
					dist = xx*xx + yy*yy;
					if ( dist <= mindist )
					 {	/* a better match found */
						best   =posr;
						mindist=dist;
						/* xx should be positive */
						xx=fabs(xx);
						/* reset limit */
						cdist = xx + yy;
						if ( cdist < maxdist )
							maxdist = cdist;
					 }
				 }
			 }
			else
				outr=0;
		 }
		else
			outr=0;
	 }
  }

 /* set minimal distance and return the closest point index */
 if ( rdistance != NULL )	*rdistance = mindist;

 return(best);
}

int cpmatch_find_neighbour(tpointarr *base,int indx)
{
 int	neighbour;
 tpoint	*probe;

 if ( base->length<=1 )			return(-1);
 if ( indx<0 || indx>=base->length )	return(-1);

 probe=&base->points[indx];

 neighbour=cpmatch_find_nearest(probe,base,indx,NULL,1);

 return(neighbour);
}

/* find closest point matches between two point arrays (see point.h)
 * returns matches in the hit structure (see cpmatch.h) */
cphit  *cpmatch(tpointarr *base,tpointarr *probe, int *hlen)
{
 int	ii,pos;
 int	baselen,probelen;
 int	bidx,pidx;
 tpoint	*basearr,*probearr,*pp;
 cphit	*hits;

 /* set return values */
 pos   = 0;
 *hlen = 0;

 /* check arguments */
 if ( base ==NULL || base->points ==NULL || base->length <1 ||
      probe==NULL || probe->points==NULL || probe->length<1 )
	return(NULL);

 bidx=0,pidx=1;

 /* shortcuts to variables in the structures */
 baselen  = base->length;
 basearr  = base->points;
 probelen = probe->length;
 probearr = probe->points;

 /* allocate return structure */
 if ( (hits = malloc(probelen*sizeof(*hits))) == NULL )
	return(NULL);

 /* sort the base point array along the x-axis */
 qsort((void *)basearr,baselen,sizeof(*basearr),tpoint_sortx);

 /* step through the probe points, and find its closest neighbour
    among the base points */

 for ( pp = probearr, ii = probelen ; ii-- ; pp++ )
  {	double	mindist=0.0;		/* cp distance squared */
	int	neighbour,best;

	/* neighbour is the index to start the search with */
	neighbour = cpmatch_find_xidx(pp, base);
	/* closest neighbour index in the base array */
	best  = cpmatch_find_nearest(pp, base, neighbour, &mindist,0);

	/* save hit to the return structure */
	hits[pos].idx[pidx] = pp - probearr;
	hits[pos].idx[bidx] = basearr[best].id;

	hits[pos].distance  = mindist;

	/* increment counter of the return length */
	pos++;
  }

 /* resize array if needed */
 if ( pos < probelen )
	hits = realloc((void *) hits, pos*sizeof(*hits));

 /* reset hits array to NULL, if it was freed by realloc() */
 if ( !pos )	hits=NULL;

 /* return length too */
 *hlen = pos;

 return(hits);
}

/* find closest point matches between two point arrays	*/
cphit  *cpmatch_symmetric(tpointarr *cat1,tpointarr *cat2,
		int *hlen,double rev_sig,double rev_maxdist)
{
 int	i,j,k,pos;
 int	cat1len,cat2len,maxlen;
 tpoint	*cat1arr,*cat2arr;
 cphit	*hits,*whits;
 int	neighbour,best;
 double	mindist=0.0,meddist,limdist;

 /* set return values */
 *hlen = 0;

 /* check arguments */
 if ( cat1==NULL || cat1->points==NULL || cat1->length<1 ||
      cat2==NULL || cat2->points==NULL || cat2->length<1 )
	return(NULL);

 /* shortcuts to variables in the structures */
 cat1len = cat1->length,cat1arr = cat1->points;
 cat2len = cat2->length,cat2arr = cat2->points;

/* fix ids */
/*
 for ( i=0 ; i<cat1len ; i++ )	cat1arr[i].id=i;
 for ( i=0 ; i<cat2len ; i++ )	cat2arr[i].id=i;
*/

 /* sort both arrays along the x-axis */
 qsort((void *)cat1arr,cat1len,sizeof(tpoint),tpoint_sortx);
 qsort((void *)cat2arr,cat2len,sizeof(tpoint),tpoint_sortx);

 maxlen=(cat1len>cat2len?cat1len:cat2len);
 hits =(cphit *)malloc(maxlen*sizeof(cphit));
 whits=(cphit *)malloc(maxlen*sizeof(cphit));
 for ( i=0 ; i<cat1len ; i++ )
  {	neighbour=cpmatch_find_xidx(&cat1arr[i],cat2);
	best=cpmatch_find_nearest(&cat1arr[i],cat2,neighbour,&mindist,0);
	whits[i].idx[0]=best;
	whits[i].distance=mindist;
  }
 for ( ; i<maxlen ; i++ ) 	whits[i].idx[0]=-1;
 for ( i=0 ; i<cat2len ; i++ )
  {	neighbour=cpmatch_find_xidx(&cat2arr[i],cat1);
	best=cpmatch_find_nearest(&cat2arr[i],cat1,neighbour,&mindist,0);
	whits[i].idx[1]=best;
  }
 for ( ; i<maxlen ; i++ ) 	whits[i].idx[1]=-1;
 pos=0;
 for ( i=0 ; i<maxlen ; i++ )
  {	j=whits[i].idx[0];
	if ( j<0 )	continue;
	k=whits[j].idx[1];
	if ( k<0 )	continue;
	if ( k==i )
	 {	hits[pos].idx[0]=cat1arr[i].id;
		hits[pos].idx[1]=cat2arr[j].id;
		hits[pos].distance=whits[i].distance;
		pos++;
	 }
  }
 free(whits);

 /* revise matched pairs a bit... */
 if ( rev_sig>0.0 )
  {	qsort(hits,pos,sizeof(cphit),cphit_sort_dist);
	meddist=0.5*(hits[(pos-1)/2].distance+hits[pos/2].distance);
	limdist=meddist*rev_sig*rev_sig;
	while ( pos>0 && hits[pos-1].distance>limdist )	pos--;
  }
 if ( rev_maxdist>0.0 )
  {	qsort(hits,pos,sizeof(cphit),cphit_sort_dist);
	limdist=rev_maxdist*rev_maxdist;
	while ( pos>0 && hits[pos-1].distance>limdist )	pos--;
  }

 hits=(cphit *)realloc(hits,pos*sizeof(cphit));

 qsort((void *)cat1arr,cat1len,sizeof(tpoint),tpoint_sortid);
 qsort((void *)cat2arr,cat2len,sizeof(tpoint),tpoint_sortid);
 
 /* reset hits array to NULL, if it was freed by realloc() */
 if ( ! pos )	hits=NULL;

 /* return length too */
 *hlen = pos;

 return(hits);
}

/* sort out "good" matches from a hits list. return the new length of
   the size of the input hit array */
int cphit_revise(int hlen,cphit **hits,int maxhlen,cpmtype mtype,double *maxdist)
{
 int	ii, newlen;
 int	*idx0,*idx1;
 int	probe;

 /* shortcut to the input array */
 cphit  *hh=*hits;

 /* check parameters: propagate error */
 if ( hlen < 1 || hh == NULL )	return(hlen);

 /* init return value */
 newlen = 0;

 /* allocate memory for temp index arrays */
 if ( (idx0 = malloc(maxhlen*sizeof(int))) == NULL ||
      (idx1 = malloc(maxhlen*sizeof(int))) == NULL )
  {	free(idx0);
	return(hlen);
  }
 memset(idx0,0,maxhlen*sizeof(int));
 memset(idx1,0,maxhlen*sizeof(int));

 /* sort input array according to distance */
 qsort((void *) hh, hlen, sizeof(*hh), cphit_sort_dist);

 if ( mtype == CPMTYPE_GUESS )
  {	for ( ; newlen < hlen ; newlen++ )
	 {	/* check the first index */
		probe = hh[newlen].idx[0];
		if ( idx0[probe] )	break;
		else			idx0[probe]=1;
		/* check the second index */
		probe = hh[newlen].idx[1];
		if ( idx1[probe] )	break;
		else			idx1[probe]=1;
	 }
  }
 else if ( mtype == CPMTYPE_ACC )
  {	int	keepon=1;
	double	misses=0.0;

	/* accumulate error: do not break out at the first unambigous match */
	for ( ii = 0 ; keepon && ii < hlen ; ii++ )
	 {	/* check the first index */
		probe = hh[ii].idx[0];
		if ( idx0[probe] )
		 {	misses += 1.0;
			if ( misses/newlen >= *maxdist )	keepon=0;
			continue;
		 }
		else	idx0[probe]=1;
		/* check the second index */
		probe = hh[ii].idx[1];
		if ( idx1[probe] )
		 {	misses += 1.0;
			if ( misses/newlen >= *maxdist )	keepon=0;
			continue;
		 }
		else	idx1[probe]=1;
		/* increase new length counter */
		newlen++;
	 }
  }
 else if ( mtype == CPMTYPE_DIST && maxdist != NULL )
  {	/* we have a user set maximal distance level */
	double	md;

	md=(*maxdist)*(*maxdist);
	for ( ii = 0 ; ii < hlen && hh[ii].distance < md ; ii++ )
	 {	/* check the first index */
		probe = hh[ii].idx[0];
		if ( idx0[probe] )	continue;
		else			idx0[probe]=1;
		/* check the second index */
		probe = hh[ii].idx[1];
		if(idx1[probe] )	continue;
		else			idx1[probe]=1;
		/* check the second index */
		newlen++;
	 }
  }

 /* write maximal distance of unambigous pairs if asked for */
 if ( maxdist != NULL )	*maxdist=hh[newlen-1].distance;

 /* release temp arrays */
 free(idx0);
 free(idx1);

 /* reallocate the hits array */
 if ( newlen != hlen )	*hits=realloc((void *)hh,newlen*sizeof(*hh));

 /* set original position if there was an error when reallocating */
 if ( *hits==NULL )	*hits=hh;

 return(newlen);
}

/*****************************************************************************/
                         
