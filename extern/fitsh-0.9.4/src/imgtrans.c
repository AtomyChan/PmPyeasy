/*****************************************************************************/
/* imgtrans.c [-> fits.c, sort.c]					     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Library for applying the discrete spatial operators on FITS images.	     */
/* The library is NOT standalone, it requires some functions from fits.a     */
/* and sort.c								     */
/* (c) 2004, Pal, A. (apal@szofi.elte.hu). 				     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in imgtrans.h      */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>

#include <fits/fits.h>

#include "statistics.h"

#include "index/sort.h"
#include "imgtrans.h"

/*****************************************************************************/

static int laplace_of_image_ign_flag(fitsimage *img,int flag)
{
 int 	i,j,sx,sy;
 double	*d0,*d1,*dd;

 if ( img==NULL )	return(1);
 if ( img->data==NULL )	return(1);

 sx=img->sx,sy=img->sy;

 if ( sy<3 || sx<3 )
  {	for ( i=0 ; i<sy ; i++ )
	 {  for ( j=0 ; j<sx ; j++ )
	     {	img->data[i][j]=0.0;		}
	 }
	return(0);
  }

 dd=(double *)malloc(sizeof(double)*sx*2);
 if ( dd==NULL )	return(-1);
 d0=dd,d1=d0+sx;
 for ( j=0 ; j<sx ; j++ )
  {	d0[j]=img->data[0][j],
	d1[j]=img->data[1][j];
  }
 for ( i=0 ; i<sy ; i++ )
  {	if ( i==0 || i==sy-1 )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	img->data[i][j]=0.0;		}
		continue;
	 }
	img->data[i][0]=img->data[i][sx-1]=0.0;

	if ( ! flag )
	 {	for ( j=1 ; j<sx-1 ; j++ )
		 {	img->data[i][j]=4*d1[j]-(d0[j]+img->data[i+1][j]+d1[j-1]+d1[j+1]);	}
	 }
	else
	 {	for ( j=1 ; j<sx-1 ; j++ )
		 {	if ( d0[j]<=0.0 || d1[j]<=0.0 || d1[j-1]<=0.0 || d1[j+1]<=0 || img->data[i+1][j]<=0.0 )
				img->data[i][j]=0.0;
			else
				img->data[i][j]=4*d1[j]-(d0[j]+img->data[i+1][j]+d1[j-1]+d1[j+1]);
		 }
	 }
	
	if ( i<sy-1 )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	d0[j]=d1[j];
			d1[j]=img->data[i+1][j];
		 }
	 }
  }
 free(dd);
 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int laplace_of_image(fitsimage *img)
{
 return ( laplace_of_image_ign_flag(img,0) );
}
int laplace_of_image_ign(fitsimage *img)
{
 return ( laplace_of_image_ign_flag(img,1) );
}

/*****************************************************************************/

static int cyclic_laplace_of_image_ign_flag(fitsimage *img,int flag)
{
 int 	i,j,sx,sy;
 double	*d0,*d1,*di,*dd;

 if ( img==NULL )	return(1);
 if ( img->data==NULL )	return(1);

 sx=img->sx,sy=img->sy;

 dd=(double *)malloc(sizeof(double)*sx*3);
 if ( dd==NULL )	return(-1);
 di=dd,d0=dd+sx,d1=d0+sx;
 for ( j=0 ; j<sx ; j++ )
  {	di[j]=d1[j]=img->data[0][j],
	d0[j]=img->data[sy-1][j];
  }
 for ( i=0 ; i<sy ; i++ )
  {	double	*dprev,*dcurr,*dnext;
	dprev=d0;
	dcurr=d1;
	if ( i<sy-1 )	dnext=img->data[i+1];
	else		dnext=di;

	if ( ! flag )
	 {	img->data[i][0   ]=4*dcurr[0   ]-(dprev[0   ]+dnext[0   ]+dcurr[sx-1]+dcurr[1]);
		for ( j=1 ; j<sx-1 ; j++ )
		 {	img->data[i][j]=4*dcurr[j]-(dprev[j]+dnext[j]+dcurr[j-1]+dcurr[j+1]);	}
		img->data[i][sx-1]=4*dcurr[sx-1]-(dprev[sx-1]+dnext[sx-1]+dcurr[sx-2]+dcurr[0]);
	 }
	else
	 {	if ( dcurr[0]<=0.0 || dprev[0]<=0.0 || dnext[0]<=0.0 || dcurr[sx-1]<=0.0 || dcurr[1]<=0.0 )
			img->data[i][0   ]=0.0;
		else 
			img->data[i][0   ]=4*dcurr[0   ]-(dprev[0   ]+dnext[0   ]+dcurr[sx-1]+dcurr[1]);
		for ( j=1 ; j<sx-1 ; j++ )
		 {	if ( dprev[j]<=0.0 || dcurr[j]<=0.0 || dcurr[j-1]<=0.0 || dcurr[j+1]<=0 || dnext[j]<=0.0 )
				img->data[i][j]=0.0;
			else
				img->data[i][j]=4*dcurr[j]-(dprev[j]+dnext[j]+dcurr[j-1]+dcurr[j+1]);
		 }
		if ( dcurr[sx-1]<=0.0 || dprev[sx-1]<=0.0 || dnext[sx-1]<=0.0 || dcurr[sx-2]<=0.0 || dcurr[0]<=0.0 )
			img->data[i][sx-1]=0.0;
		else 
			img->data[i][sx-1]=4*dcurr[sx-1]-(dprev[sx-1]+dnext[sx-1]+dcurr[sx-2]+dcurr[0]);
	 }
	
	if ( i<sy-1 )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	d0[j]=d1[j];
			d1[j]=img->data[i+1][j];
		 }
	 }
  }
 free(dd);
 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int cyclic_laplace_of_image(fitsimage *img)
{
 return ( cyclic_laplace_of_image_ign_flag(img,0) );
}
int cyclic_laplace_of_image_ign(fitsimage *img)
{
 return ( cyclic_laplace_of_image_ign_flag(img,1) );
}
        
/*****************************************************************************/

/*
static int median_block_filter_compare(int i1,int i2,void *param)
{
 double	*rawdata;
 rawdata=(double *)param;
 if ( rawdata[i1] < rawdata[i2] )	return(-1);
 else					return(1);
}

int median_block_filter(fitsimage *img,fitsimage *med,char **mask,int hb)
{
 int	sx,sy,i,j,hs;
 int	imin,imax,in,ii,
	jmin,jmax,jn,jj;
 int	*idx,n;
 double	*rawdata;

 if ( img==NULL || med==NULL )			return(1);
 if ( img->data==NULL || med->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;
 if ( sx != med->sx || sy != med->sy )		return(1);

 hs=2*hb+1;

 rawdata=(double *)malloc(sizeof(double)*hs*sx);
 if ( rawdata==NULL )	return(-1);
 idx=(int *)malloc(sizeof(int)*hs*hs);
 if ( idx==NULL )	return(-1);

 for ( i=0 ; i<sy ; i++ )
  {	imin=i-hb;if ( imin<0   ) imin=0;
	imax=i+hb;if ( imax>=sy ) imax=sy-1;
	in=imax-imin+1;
	for ( ii=0 ; ii<in ; ii++ )
	 {	memcpy(rawdata+ii*sx,img->data[imin+ii],sizeof(double)*sx);	}
	for ( j=0 ; j<sx ; j++ )
	 {	jmin=j-hb;if ( jmin<0   ) jmin=0;
		jmax=j+hb;if ( jmax>=sx ) jmax=sx-1;
		jn=jmax-jmin+1;
		if ( mask==NULL )
		 {	for ( ii=0,n=0 ; ii<in ; ii++ )
			 {	for ( jj=0 ; jj<jn ; jj++,n++ )
				 {	idx[n]=ii*sx+jmin+jj;		}
			 }
		 }
		else
		 {	for ( ii=0,n=0 ; ii<in ; ii++ )
			 {	for ( jj=0 ; jj<jn ; jj++ )
				 {	if ( mask[imin+ii][jmin+jj] )	continue;
					idx[n]=ii*sx+jmin+jj;
					n++;
				 }
			 }
		 }
		if ( n==0 )
		 {	med->data[i][j]=0.0;			}
		else if ( n==1 )
		 {	med->data[i][j]=rawdata[idx[0]];	}
		else
		 {	int	p1,p2;
			index_qsort(idx,n,median_block_filter_compare,(void *)rawdata);
			p1=(n-1)/2,p2=(n)/2;
			med->data[i][j]=0.5*(rawdata[idx[p1]]+rawdata[idx[p2]]);
		 }
	 }
  }
	

 free(idx);
 free(rawdata);

 return(0);
}
*/
