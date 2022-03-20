/*****************************************************************************/
/* linkpoint.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions related to linking points (see linkpoint.h for a more detailed  */
/* description). This simple standalone library is the core of some 	     */
/* sophisticated object detection algorithms.				     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2006, 2007; Pal, A. (apal@szofi.elte.hu).			     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "linkpoint.h"

/*****************************************************************************/

static void *link_alloc(int typesize,int sx,int sy)
{
 size_t	tsize,psize,bsize,cd,cb,snext,pnext;
 size_t	i,j;
 int	arr[2],rank;
 void	*ret,**pret;

 rank=2;
 arr[0]=sx,arr[1]=sy;

 psize=0;bsize=1;
 for ( i=rank-1 ; i>=1 ; i-- )
  {	bsize=bsize*arr[i];
	psize+=bsize;
  }
 psize=psize*sizeof(void *);
 tsize=psize+typesize*bsize*arr[0];
 pret=(void **)malloc(tsize);
 if ( pret==NULL )	return(NULL);

 psize=0;bsize=1;
 for ( i=rank-1 ; i>=1 ; i-- )
  {	cd=arr[i];
	if ( i>1 )	snext=sizeof(void *);
	else		snext=typesize;
	cb=cd*bsize;
	pnext=psize+cb;
	for ( j=0 ; j<cb ; j++ )
		pret[psize+j]=(void *)( ((char *)(&pret[pnext]))+j*snext*arr[i-1] );
	psize=pnext;
	bsize=cb;
  }
 ret=(void *)pret;
 return(ret);
}

static int link_free(void *arr)
{
 free(arr);
 return(0);
}

/*****************************************************************************/

int linkpoint_local_extreme(double **arr,int ox,int oy,int sx,int sy,
	char **mask,int j,int i,
	int flag,int exclude_mask,int *rnx,int *rny)
{
 double	wmin,wmax,w;
 int	maxx,minx,maxy,miny,n;
 int	k0,k,k1,l0,l,l1,msi,msd,cmask;

 if ( arr==NULL || sx<=0 || sy<=0 )
	return(-1);
 if ( mask != NULL && mask[i][j] )
	return(0);

 k0=i-1;msi=0x01;if ( k0<oy  )	k0=oy,msi=0x08;
 l0=j-1;msd=0;   if ( l0<ox )	l0=ox,msd=1;
 k1=i+1; if ( k1>=sy )	k1=sy-1;
 l1=j+1; if ( l1>=sx )	l1=sx-1;

 if ( ! exclude_mask )
  {	wmin=wmax=arr[i][j];
	n=0;
	maxx=minx=j;
	maxy=miny=i;

	if ( mask != NULL )
	 {	for ( k=k0 ; k<=k1 ; k++ )
		 {	for ( l=l0 ; l<=l1 ; l++ )
			 {	if ( mask[k][l] ) continue;
				w=arr[k][l];
				if ( w>wmax )	   wmax=w,maxx=l,maxy=k;
				else if ( w<wmin ) wmin=w,minx=l,miny=k;
				n++;
			 }
		 }
	 }
	else
	 {	for ( k=k0 ; k<=k1 ; k++ )
		 {	for ( l=l0 ; l<=l1 ; l++ )
			 {	w=arr[k][l];
				if ( w>wmax )	   wmax=w,maxx=l,maxy=k;
				else if ( w<wmin ) wmin=w,minx=l,miny=k;
				n++;
			 }
		 }
	 }
  }
 else
  {	wmin=wmax=arr[i][j];
	n=0;
	maxx=minx=j;
	maxy=miny=i;

	if ( mask != NULL )
	 {	for ( k=k0 ; k<=k1 ; k++,msi<<=3 )
		 {	for ( l=l0,cmask=msi<<msd ; l<=l1 ; l++,cmask<<=1 )
			 {	if ( mask[k][l] || (exclude_mask & cmask) ) continue;
				w=arr[k][l];
				if ( ! n )	   	wmax=wmin=w,maxx=minx=l,maxy=miny=k;
				else if ( w>wmax )	wmax=w,maxx=l,maxy=k;
				else if ( w<wmin )	wmin=w,minx=l,miny=k;
				n++;
			 }
		 }
	 }
	else
	 {	for ( k=k0 ; k<=k1 ; k++,msi<<=3 )
		 {	for ( l=l0,cmask=msi<<msd ; l<=l1 ; l++,cmask<<=1 )
			 {	if ( exclude_mask & cmask )	continue;
				w=arr[k][l];
				if ( ! n )	   	wmax=wmin=w,maxx=minx=l,maxy=miny=k;
				else if ( w>wmax )	wmax=w,maxx=l,maxy=k;
				else if ( w<wmin )	wmin=w,minx=l,miny=k;
				n++;
			 }
		 }
	 }
  }

 if ( n>0 )
  {	if ( flag>0 )
	 {	*rnx=maxx;
		*rny=maxy;
	 }
	else
	 {	*rny=minx;
		*rny=miny;
	 }
  }

 return(n);
}

int linkpoint_mask_same_group(linkpoint **lparr,int sx,int sy,int j,int i)
{
 int	k0,k,k1,l0,l,l1,mx,my;
 int	msi,msd,cmask,rmask;

 k0=i-1;msi=0x01;if ( k0<0  )	k0=0,msi=0x08;
 l0=j-1;msd=0;   if ( l0<0 )	l0=0,msd=1;
 k1=i+1; if ( k1>=sy )	k1=sy-1;
 l1=j+1; if ( l1>=sx )	l1=sx-1;

 rmask=0;

 mx=lparr[i][j].mx;
 my=lparr[i][j].my;

 for ( k=k0 ; k<=k1 ; k++,msi<<=3 )
  {	for ( l=l0,cmask=msi<<msd ; l<=l1 ; l++,cmask<<=1 )
	 {	if ( lparr[k][l].nx<0 || lparr[k][l].ny<0 )
			rmask |= cmask;
		else if ( lparr[k][l].mx==mx && lparr[k][l].my==my )
			rmask |= cmask;
	 }
  }
 
 return(rmask);
}


/*****************************************************************************/

linkpoint **linkpoint_create(double **arr,int sx,int sy,linkrange *lr,char **mask,int flag)
{
 linkpoint	**lparr;
 int		i,j,imin,imax,jmin,jmax,nx,ny;
 int		n;

 if ( ! flag )	return(NULL);		/* should be negative or positive */

 if ( lr==NULL )
  {	imin=0,imax=sy-1,
	jmin=0,jmax=sx-1;
  }
 else
  {	imin=lr->ymin,imax=lr->ymax;
	jmin=lr->xmin,jmax=lr->xmax;
	if ( imin>imax )	i=imin,imin=imax,imax=i;
	if ( jmin>jmax )	j=jmin,jmin=jmax,jmax=j;
  }
 if ( imin<0 )		imin=0;
 if ( imax>=sy )	imax=sy-1;
 if ( jmin<0 )		jmin=0;
 if ( jmax>=sx )	jmax=sx-1;

 lparr=(linkpoint **)link_alloc(sizeof(linkpoint),sx,sy);

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	lparr[i][j].nx=lparr[i][j].ny=-1;
		lparr[i][j].mx=lparr[i][j].my=-1;
	 } 
  }

 for ( i=imin ; i<=imax ; i++ )
  {	for ( j=jmin ; j<=jmax ; j++ )
	 {	
		if ( mask != NULL && mask[i][j] )
			continue;

		nx=ny=0;
		n=linkpoint_local_extreme(arr,jmin,imin,jmax-1,imax-1,mask,j,i,flag,0,&nx,&ny);
		if ( n>0 )
		 {	lparr[i][j].nx=nx;
			lparr[i][j].ny=ny;
		 }
		/* lparr[i][j].mx=lparr[i][j].my=-1; */
	 }
  }

 return(lparr);
}

int linkpoint_get_link_end(linkpoint **lparr,int nx,int ny,int *rmx,int *rmy)
{
 int	wx,wy;
 
 if ( lparr[ny][nx].nx < 0 || lparr[ny][nx].ny < 0 )
	return(-1);

 while ( lparr[ny][nx].nx != nx || lparr[ny][nx].ny != ny )
  {	wx=lparr[ny][nx].nx,
	wy=lparr[ny][nx].ny;
	nx=wx,
	ny=wy;
  }

 *rmx=nx,
 *rmy=ny;

 return(0);
}

int linkpoint_is_same_endpoint(linkpoint **lparr,int x1,int y1,int x2,int y2)
{
 int	mx1,my1,mx2,my2;

 if ( linkpoint_get_link_end(lparr,x1,y1,&mx1,&my1) < 0 )
	return(-1);
 else if ( linkpoint_get_link_end(lparr,x2,y2,&mx2,&my2) < 0 )
	return(-1);
 else if ( mx1==mx2 && my1==my2 )
	return(1);
 else 
	return(0);
}


int linkpoint_reconnect(linkpoint **lparr,int sx,int sy)
{
 int	i,j,mx,my;

 if ( lparr==NULL || sx<=0 || sy<=0 )
	return(1);

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( lparr[i][j].nx<0 || lparr[i][j].ny<0 )	continue;
		if ( ! linkpoint_get_link_end(lparr,j,i,&mx,&my) )
		 {	lparr[i][j].mx=mx,
			lparr[i][j].my=my;
		 }
	 }
  }

 return(0);
}

int linkpoint_free(linkpoint **arr)
{
 link_free((void *)arr);
 return(0);
}

/*****************************************************************************/

linkreference **linkreference_create(linkpoint **lparr,int sx,int sy)
{
 linkreference	**lrarr;
 int		i,j,wx,wy;
 int		k0,k,k1,l0,l,l1,nngh,nsgn;

 if ( lparr==NULL || sx<=0 || sy<=0 )
	return(NULL);

 lrarr=(linkreference **)link_alloc(sizeof(linkreference),sx,sy);
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	lrarr[i][j].ncnt=0;
		lrarr[i][j].mcnt=0;
		lrarr[i][j].nngh=0;
		lrarr[i][j].nsgn=0;
		lrarr[i][j].identifier=-1;
		lrarr[i][j].aux=0;
	 }
  }

 for ( i=0 ; i<sy ; i++ )
  {	k0=i-1; if ( k0<0  )	k0=0;
	k1=i+1; if ( k1>=sy )	k1=sy-1;
	for ( j=0 ; j<sx ; j++ )
	 {	if ( lparr[i][j].nx<0 || lparr[i][j].ny<0 )
		 {	lrarr[i][j].ncnt=-1;
			lrarr[i][j].mcnt=-1;
			lrarr[i][j].nngh=-1;
			lrarr[i][j].nsgn=-1;
			continue;
		 }

		wx=lparr[i][j].nx;
		wy=lparr[i][j].ny;
		lrarr[wy][wx].ncnt++;

		wx=lparr[i][j].mx;
		wy=lparr[i][j].my;
		lrarr[wy][wx].mcnt++;

		l0=j-1; if ( l0<0 )	l0=0;
		l1=j+1; if ( l1>=sx )	l1=sx-1;
	
		nngh=0;
		nsgn=0;
		for ( k=k0 ; k<=k1 ; k++ )
		 {	for ( l=l0 ; l<=l1 ; l++ )
			 {	if ( lparr[k][l].nx<0 || lparr[k][l].ny<0 )
					continue;
				nngh++;
				if ( lparr[k][l].mx==wx && lparr[k][l].my==wy )
					nsgn++;
			 }
		 }
		lrarr[i][j].nngh=nngh-1;
		lrarr[i][j].nsgn=nsgn-1;		
	 } 
  }

 return(lrarr);
}

int linkreference_free(linkreference **arr)
{
 link_free((void *)arr);
 return(0);
}

/*****************************************************************************/
                                            
