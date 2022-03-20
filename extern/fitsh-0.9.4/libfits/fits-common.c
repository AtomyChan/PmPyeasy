/*****************************************************************************/
/* fits-common.c							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Another simple standalone library for manipulating FITS files:	     */
/* Declaration of some global constants and common allocation functions.     */
/* (c) 2004-06, Pal, A. (apal@szofi.elte.hu). 				     */
/* See reference(s) at the end of this source code.			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in fits.h	     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#define	_FITS_SOURCE

#include <fits/fits.h>
#include "fits-common.h"

/*****************************************************************************/

char *headers[]=
 {	"NAXIS"   ,		/*  0 */
	"NAXIS1"  ,		/*  1 */
	"NAXIS2"  ,		/*  2 */
	"NAXIS3"  ,		/*  3 */
	"NAXIS4"  ,		/*  4 */
	"BITPIX"  ,		/*  5 */
	"BSCALE"  ,		/*  6 */	
	"BZERO"   ,		/*  7 */
	"EXTEND"  ,		/*  8 */
	"XTENSION",		/*  9 */
	"INHERIT" ,		/* 10 */
	"SIMPLE"  ,		/* 11 */
	"ORIGIN"  ,		/* 12 */
	"GAIN"    ,		/* 13 */
	"END"	  ,		/* 14 */
	"TFIELDS" ,		/* 15 */
	"TBCOL"	  ,		/* 16 */
	"TFORM"	  ,		/* 17 */
	"TSCAL"	  ,		/* 18 */
	"TZERO"	  ,		/* 19 */
	"TNULL"	  ,		/* 20 */
	"TTYPE"	  ,		/* 21 */
	"TUNIT"	  ,		/* 22 */
	"PCOUNT"  ,		/* 23 */
	"GCOUNT"  ,		/* 24 */
	NULL
 };

char *comment_fits_standard="Fits standard";

int substantial_headers[]=
 {	HDR_SIMPLE,
	HDR_NAXIS,
	HDR_BITPIX,
	HDR_EXTEND,
	HDR_XTENSION,
	-1
 };

/*****************************************************************************/

void *fits_tensor_alloc_arr(int typesize,int rank,int *arr)
{
 size_t	tsize,psize,bsize,cd,cb,snext,pnext;
 size_t	i,j;
 void	*ret,**pret;

/*
 if ( 1 )
  {	int	i;
	fprintf(stderr,"fits_tensor_alloc_arr(): rank=%d",rank);
	for ( i=0; i<rank; i++ )
	 {	fprintf(stderr," arr[%d]=%d",i,arr[i]);		}
	fprintf(stderr,"\n");
  }
*/

 if ( rank<=0 )	
	return(NULL);

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

int fits_tensor_free(void *tensor)
{
 free(tensor);
 return(0);
}

/*****************************************************************************/

int fits_arch_is_swapped(void)
{
 short	a;
 a=1;
 if ( *(char *)(&a) == 1 )	return(1);
 else				return(0);
}

/*****************************************************************************/

int fits_swap_line_bytes(unsigned char *wbuff,int bs,int sx)
{
 int		j;
 unsigned char	w;

 if ( fits_arch_is_swapped() )
  {	switch ( bs )
	 {   case 2:
		for ( j=0 ; j<2*sx ; j+=2 )
		 {	w=wbuff[j+0],wbuff[j+0]=wbuff[j+1],wbuff[j+1]=w;	}
		break;
	     case 4:
		for ( j=0 ; j<4*sx ; j+=4 )
		 {	w=wbuff[j+0],wbuff[j+0]=wbuff[j+3],wbuff[j+3]=w;
			w=wbuff[j+1],wbuff[j+1]=wbuff[j+2],wbuff[j+2]=w;
		 }
		break;
	     case 8:
		for ( j=0 ; j<8*sx ; j+=8 )
		 {	w=wbuff[j+0],wbuff[j+0]=wbuff[j+7],wbuff[j+7]=w;
			w=wbuff[j+1],wbuff[j+1]=wbuff[j+6],wbuff[j+6]=w;
			w=wbuff[j+2],wbuff[j+2]=wbuff[j+5],wbuff[j+5]=w;
			w=wbuff[j+3],wbuff[j+3]=wbuff[j+4],wbuff[j+4]=w;
		 }
		break;
	 }
  }
 return(0);
}

/*****************************************************************************/

int fits_cb_skip(void *param,int length)
{
 if ( param==NULL )
	return(-1);
 else
  {	fseek((FILE *)param,(long)length,SEEK_CUR);
	return(length); 
  }
}
int fits_cb_is_end(void *param)
{
 if ( param==NULL )
	return(-1);
 else if ( feof((FILE *)param) )
	return(1);
 else
	return(0);
}
int fits_cb_read(void *param,void *dst,int length)
{
 if ( param==NULL )
	return(-1); 
 else if ( fits_cb_is_end(param) )
	return(0);
 else if ( dst != NULL )
	return(fread(dst,1,length,(FILE *)param));
 else
	return(fits_cb_skip(param,length));
}
int fits_cb_write(void *param,void *src,int length)
{
 if ( param==NULL )
	return(-1);
 else if ( src != NULL )
	return(fwrite(src,1,length,(FILE *)param));
 else
	return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int fits_cb_mem_write(void *param,void *src,int length)
{
 fitsmemwrite	*fmw;

 fmw=(fitsmemwrite *)param;

 if ( fmw==NULL )
	return(-1);
 else if ( src != NULL )
  {	fmw->buffer=(char *)realloc(fmw->buffer,fmw->length+length);
	memcpy(fmw->buffer+fmw->length,src,length);
	fmw->length+=length;
	return(length);
  }
 else
	return(0);
}

/******************************************************************************
 Reference:
  [1]:	Definition of the Flexible Image Transport System (FITS)
	NASA / Science Office of Standards and Technology (NOST)
	NOST 100-2.0 (March 29, 1999)
	NASA Goddard Space Flight Center
	http://fits.gsfc.nasa.gov/fits_documentation.html
******************************************************************************/
                                                      
