/*****************************************************************************/
/* multiindex.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library for (currently only) 2-dimensional orthogonal range    */
/* search.								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2006; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "multiindex.h"

/*****************************************************************************/

int get_nbit(int x)
{
 int    nbit;
 for ( nbit=0 ; x>0 ; )
  {	x=x/2;
	nbit++;
  }
 return(nbit);
}

/*****************************************************************************/

int index_quicksort_rec(int *index,int l,int r,
	int (*compare)(int,int,void *),void *param)
{
 int	i,j,tmp,ch;

 ch=index[(l+r)/2];
 i=l-1,j=r+1;
 while ( 1 )
  {	do { i++; } while ( index[i] != ch && compare(index[i],ch,param)<0 ) ;
	do { j--; } while ( index[j] != ch && compare(index[j],ch,param)>0 ) ;
	if ( i<j )
		tmp=index[j],index[j]=index[i],index[i]=tmp;
	else if ( i==j )
	 {	i++;
		break;
	 }
	else	break;
  };
 if ( l<j )	index_quicksort_rec(index,l,j,compare,param);
 if ( i<r )	index_quicksort_rec(index,i,r,compare,param);
 return(0);
}

int index_quicksort(int *index,int size,
	int (*compare)(int,int,void *),void *param)
{
 index_quicksort_rec(index,0,size-1,compare,param);
 return(0);
}

int index_subsort(int *index,int (*compare)(int,int,void *),
	int level,int min,int max,void *param)
{
 int	mid;
 if ( level>0 )
  {	mid=(min+max)/2;
	index_subsort(index,compare,level-1,min,mid,param);
	index_subsort(index,compare,level-1,mid,max,param);
	return(0);
  }
 else if ( max>min+1 ) 
  {	index_quicksort(index+min,max-min,compare,param);
	return(0);
  }
 else
	return(0);
}

/*****************************************************************************/

/* functions related to 1d range searching */

int check_if_in_interval(void *param,int p,
	int (*compare)(int,double,void *),double *xl,double *xr)
{
 if ( xl != NULL && ! ( 0 <= compare(p,*xl,param) ) )		/* II */
	return(0);
 else if ( xr != NULL && ! ( compare(p,*xr,param)<=0 ) )	/* II */
	return(0);
 else
	return(1);
}

int search_index_boundaries(void *param,int *index,int size,
	int (*compare)(int,double,void *),
	double *rx0,int *rleft ,nodeindex *nleft,
	double *rx1,int *rright,nodeindex *nright)
{
 nodeindex	n;
 int		min,max,mid,mask;
 double		x0,x1;

 if ( rleft == NULL || rright == NULL )	return(-1);

 if ( rx0==NULL )
  {	*rleft=0;
	if ( nleft != NULL )
	 {	nleft->depth=get_nbit(size)-1;
		nleft->value=0;
	 }
  }
 else
  {	x0=*rx0;
	if ( ! ( 0 <= compare(index[size-1],x0,param) ) )	/* II */
	 {	*rleft=1;
		*rright=0;
		return(1);
	 }
	min=0;
	max=size;
	n.depth=0;
	n.value=0;
	mask=1;
	while ( max>min+1 )
	 {	mid=(min+max)/2;
		if ( 0 <= compare(index[mid-1],x0,param) )	/* II */
			max=mid;
		else
			min=mid,
			n.value|=mask;
		n.depth++;
		mask<<=1;
	 }
	*rleft=min;
	if ( nleft != NULL )
	 {	nleft->depth=n.depth;
		nleft->value=n.value;
	 }
  }

 if ( rx1==NULL )
  {	*rright=size-1;
	if ( nright != NULL )
	 {	nright->depth=get_nbit(size-1);
		nright->value=(1<<nright->depth)-1;
	 };
  }	
 else
  {	x1=*rx1;
	if ( ! ( compare(index[0],x1,param)<=0 ) )		/* II */
	 {	*rleft=1;
		*rright=0;
		return(1);
	 }
	min=0;
	max=size;
	n.depth=0;
	n.value=0;
	mask=1;
	while ( max>min+1 )
	 {	mid=(min+max)/2;
		if ( compare(index[mid],x1,param)<=0 )		/* II */
			min=mid,
			n.value|=mask;
		else	
			max=mid;
		n.depth++;
		mask<<=1;
	 };
	*rright=min;
	if ( nleft != NULL )
	 {	nright->depth=n.depth;
		nright->value=n.value;
	 }
  }

 if ( *rleft > *rright )	return(1);	/* epmty interval */
 else				return(0);	
}

/*****************************************************************************/

/* functions related to walking on a static binary tree */

int get_nodeinterval_rec(int ileft,int iright,int min,int max,
	int depth,int value,nodeinterval **nis)
{
 int	nnileft,nniright,mid;
 if ( iright < min || max <= ileft )
	return(0);
 else if ( ileft<=min && max-1<=iright )
  {	(*nis)->i.depth=depth;
	(*nis)->i.value=value;
	(*nis)->min=min;
	(*nis)->max=max;
	(*nis)++;
	return(1);
  }
 else
  {	mid=(min+max)/2;
	nnileft =get_nodeinterval_rec(ileft,iright,
			min,mid,depth+1,value,nis);
	nniright=get_nodeinterval_rec(ileft,iright,
			mid,max,depth+1,value|(1<<depth),nis);
	return(nnileft+nniright);
  }
}

int get_nodeintervals(int ileft,int iright,int size,nodeinterval *nis)
{
 int	nni;
 nni=get_nodeinterval_rec(ileft,iright,0,size,0,0,&nis);
 return(nni);
}

/*****************************************************************************/

int multiindex_create(void *param,int size,
	int (*compare_x)(int,int,void *),
	int (*compare_y)(int,int,void *),
	multiindex *mi)
{
 int	nbit,i;

 nbit=get_nbit(size);

 mi->param=param;

 mi->size=size;
 mi->nbit=nbit;

 mi->i_pri=(int *)malloc(sizeof(int)*size);
 mi->i_sec=(int *)malloc(sizeof(int)*size*nbit);

 mi->ind_x=mi->i_pri;
 for ( i=0 ; i<size ; i++ )
  {	mi->ind_x[i]=i;				}
 mi->ind_xy=(int **)malloc(sizeof(int *)*nbit);
 for ( i=0 ; i<nbit ; i++ )
  {	mi->ind_xy[i]=&mi->i_sec[size*i];	}

 index_quicksort(mi->ind_x,size,compare_x,param);
 for ( i=0 ; i<nbit ; i++ )
  {	memcpy(mi->ind_xy[i],mi->ind_x,sizeof(int)*size);
	index_subsort(mi->ind_xy[i],compare_y,i,0,size,param);
  }
 
 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int multiindex_range_query(multiindex *mi,
	int (*compare_value_x)(int,double,void *),
	int (*compare_value_y)(int,double,void *),
	double xl,double xr,double yl,double yr,
	int *ret,int maxret)
{
 int		il,ir,iyl,iyr,nni,maxdepth;
 int		*lindy,depth,min,max,i,j,k,nret;
 nodeinterval	*nis;

 if ( mi==NULL )
	return(-1);
 if ( mi->param==NULL || mi->ind_x==NULL || mi->ind_xy==NULL )
	return(-1);

 search_index_boundaries(mi->param,mi->ind_x,mi->size,compare_value_x,
	&xl,&il,NULL,&xr,&ir,NULL);

 nis=(nodeinterval *)malloc(sizeof(nodeinterval)*2*mi->nbit);
 nni=get_nodeintervals(il,ir,mi->size,nis);

 maxdepth=mi->nbit;

 nret=0;
 for ( i=0 ; i<nni ; i++ )
  {	depth=nis[i].i.depth;
	min=nis[i].min;
	max=nis[i].max;
	if ( depth >= maxdepth || depth >= mi->nbit )
	 {	for ( j=min ; j<max ; j++ )
		 {	k=mi->ind_x[j];
			if ( !  check_if_in_interval(mi->param,k,
				compare_value_y,&yl,&yr) )
					continue;
			else if ( nret<maxret )
			 {	ret[nret]=k;
				nret++;
			 }
		 }
	 }
	else
	 {	lindy=mi->ind_xy[depth];
		search_index_boundaries(mi->param,lindy+min,max-min,
			compare_value_y,&yl,&iyl,NULL,&yr,&iyr,NULL);
		iyl+=min;
		iyr+=min;
		for ( j=iyl ; j<=iyr ; j++ )
		 {	k=lindy[j];
			if ( nret<maxret )
			 {	ret[nret]=k;
				nret++;
			 }
		 }
	 }
  }

 free(nis);

 return(nret); 
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int multiindex_reset(multiindex *mi)
{
 if ( mi==NULL )	return(-1);

 mi->nbit=mi->size=0;
 mi->param=NULL;
 mi->i_pri=mi->i_sec=NULL;
 mi->ind_x=NULL;
 mi->ind_xy=NULL;
 
 return(0);
}

int multiindex_free(multiindex *mi)
{
 if ( mi==NULL )	return(-1);

 if ( mi->i_pri != NULL )	free(mi->i_pri);
 if ( mi->i_sec != NULL )	free(mi->i_sec);
 if ( mi->ind_xy != NULL )	free(mi->ind_xy);

 multiindex_reset(mi);
 
 return(0);
}

/*****************************************************************************/
                                                     
