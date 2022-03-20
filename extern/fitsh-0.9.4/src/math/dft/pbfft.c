/*****************************************************************************/
/* pbfft.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library for Fast Fourier Transformation, based on the 	     */
/* prime-based FFT (PBFFT) algorithm.					     */
/* (c) 2001-2003, 2004, Pal A. (apal@szofi.elte.hu)			     */
/* Last modified: 2006.10.23						     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in pbfft.h.	     */
/*****************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "pbfft.h"

#ifndef	M_PI
#define	M_PI 3.1415926535897932384626433832795028841968
#endif

typedef struct
 {	int	prime[32];
	complex *exp,*work;
	int	size,type;
 } pbfftcache;

static pbfftcache cache={ {1}, NULL,NULL, 0, 0 };

/* TBD: should be obsoleted by the upcoming implementation, see the commented
   code of pbfft_getprimdivs() just right after this function. */
static void pbfft_getprimdivs(int n,int *r)
{
 int k;

 while ( n>1 )
  {	if ( n%2==0 )	*r=2,r++,n/=2;
	else		break;
  }

 for ( k=3 ; n>1 ; )
  {	if ( n%k==0 )	*r=k,r++,n=n/k;
	else		k+=2;
  }
 *r=0;
}

/*
static void pbfft_getprimdivs(int n,int *r)
{
 int k,w;

 while ( n>1 )
  {	if ( n%2==0 )	*r=2,r++,n/=2;
	else		break;
  }

 for ( k=3,w=9 ; n>1 && w<=n ; )
  {	if ( n%k==0 )	*r=k,r++,n=n/k;
	else		w+=4*k+4,k+=2;
  }
 if ( n>1 )		*r=n,r++;
 *r=0;
}
*/

static int pbfft_get_reversal(int n,int *pbb)
{
 int	i,j,wd[32];

 i=0;
 while ( *pbb )
  {	wd[i]=n%(*pbb),
	n=n/(*pbb);
	i++,pbb++;
  };
 j=1,n=0;
 while ( i )
  {	i--,pbb--;
	n+=j*wd[i],j=j*(*pbb);
  };
 return(n);
}

static void pbfft_setpwr(complex *x,double a)
{
 x->re=cos(a);
 x->im=sin(a);
}

static void pbfft_set_expbuffer(complex *eb,int n,int typ)
{
 int	i;

 if ( typ )	typ=1;
 else		typ=-1;

 for ( i=0 ; i<n ; i++ )
  {	pbfft_setpwr(eb,2.0*M_PI*(double)(i*typ)/(double)n);
	eb++;
  }
}

static void pbfft_permute(complex *dt,complex *dt2,int n,int *pbb)
{
 int	i,j;

 for ( i=0 ; i<n ; i++ )
  {	j=pbfft_get_reversal(i,pbb);
	dt2[i].re=dt[j].re;
	dt2[i].im=dt[j].im;
  }

 memcpy(dt,dt2,sizeof(complex)*n);
}

static void pbfft_conv_base2(complex *dt,complex *eb,int k0,int e0)
{		
 int	k,a;
 double	re,im,re1,im1,re2,im2;

 for ( k=0 ; k<k0 ; k++ )
  {	a=k*e0;
	re=eb[a].re*dt[k+k0].re-eb[a].im*dt[k+k0].im;
	im=eb[a].im*dt[k+k0].re+eb[a].re*dt[k+k0].im;
	re1=dt[k].re+re;
	im1=dt[k].im+im;
	re2=dt[k].re-re;
	im2=dt[k].im-im;
	dt[k].re=re1;
	dt[k].im=im1;
	dt[k+k0].re=re2;
	dt[k+k0].im=im2;
  }
}
static int pbfft_conv_base(complex *dt,complex *dt2,
			   complex *eb,int k0,int e0,int s0)
{
 int	k,a,b,n,s;
 double	re0,im0;

 n=k0*s0;
 for ( k=0 ; k<n ; k++ )
  {	re0=im0=0.0;
	for ( s=0 ; s<s0 ; s++ )
	 {	a=((s*k)%n)*e0;
		b=k%k0+s*k0;
		re0+=eb[a].re*dt[b].re-eb[a].im*dt[b].im;
		im0+=eb[a].im*dt[b].re+eb[a].re*dt[b].im;
	 }
	dt2[k].re=re0;
	dt2[k].im=im0;
  }

 memcpy((char *)dt,(char *)dt2,sizeof(complex)*n);
 return(0);
}
int pbfft_conv(complex *dt,int n,int type)
{
 int		i,k,s,e0,*pbb;
 complex	*dt2;

 if ( type )	type=1;

 pbb=cache.prime;

 if ( cache.size != n )
  {	if ( cache.exp==NULL )
		cache.exp =(complex *)malloc (sizeof(complex)*n);
	else	
		cache.exp =(complex *)realloc(cache.exp ,sizeof(complex)*n);
	if ( cache.work==NULL )
		cache.work=(complex *)malloc (sizeof(complex)*n);
	else	
		cache.work=(complex *)realloc(cache.work,sizeof(complex)*n);
	cache.size=n;
	pbfft_getprimdivs(n,pbb);
	pbfft_set_expbuffer(cache.exp,n,type);
	cache.type=type;
  }
 if ( cache.type != type )
  {	pbfft_set_expbuffer(cache.exp,n,type);
	cache.type=type;
  }

 dt2=cache.work;

 pbfft_permute(dt,dt2,n,pbb);

 k=1;e0=n;

 while ( *pbb==2 )
  {	e0/=2;
	for ( i=0 ; i<n ; i+=2*k )
	 {	pbfft_conv_base2(dt+i,cache.exp,k,e0);			}
	k=k*2;pbb++;
  };

 while ( *pbb )
  {	s=*pbb;e0/=s;
	for ( i=0 ; i<n ; i+=k*s )
	 {	pbfft_conv_base(dt+i,dt2+i,cache.exp,k,e0,s);		}
	k=k*s;
	pbb++;	
  };

 return(0);
}

/*****************************************************************************/
                                    
