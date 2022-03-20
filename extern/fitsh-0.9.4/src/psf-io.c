/*****************************************************************************/
/* psf-io.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions related to PSF I/O (reading/writing/dumping PSF information).   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2004-2006; Pal, A. (apal@szofi.elte.hu)			             */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>

#include <fits/fits.h>

#include "tensor.h"

#include "psf.h"
#include "psf-io.h"

/*****************************************************************************/

int psf_write(FILE *fw,psf *p)
{
 int	i,j,k,wwd,nvar;
 double	dx,dy,igrid;

 wwd=p->grid*(2*p->hsize+1);
 nvar=(p->order+1)*(p->order+2)/2;
 igrid=1.0/(double)(p->grid);

 for ( i=0 ; i<wwd ; i++ )
  {	for ( j=0 ; j<wwd ; j++ )
	 {	dx=igrid*((double)j-0.5*(double)(wwd-1));
		dy=igrid*((double)i-0.5*(double)(wwd-1));
		fprintf(fw,"%11g %11g ",dx,dy);
		for ( k=0 ; k<nvar ; k++ )
			fprintf(fw,"%11g ",p->coeff[k][i][j]);
		fprintf(fw,"\n");
	 }
  }

 return(0);
}

/*****************************************************************************/

#define		PSFHDR_PSFHSIZE		0
#define		PSFHDR_PSFSGRID		1
#define		PSFHDR_PSFORDER		2
#define		PSFHDR_PSFOFFSX		3
#define		PSFHDR_PSFOFFSY		4
#define		PSFHDR_PSFSCALE		5

static char *psfheaders[] = 
 {	"PSFHSIZE",	/* 0,	PSFHDR_PSFHSIZE		*/
	"PSFSGRID",	/* 1,	PSFHDR_PSFSGRID		*/
	"PSFORDER",	/* 2,	PSFHDR_PSFORDER		*/
	"PSFOFFSX",	/* 3,	PSFHDR_PSFOFFSX		*/
	"PSFOFFSY",	/* 4,	PSFHDR_PSFOFFSY		*/
	"PSFSCALE",	/* 5,	PSFHDR_PSFSCALE		*/
	NULL
 };

int psf_write_fits(FILE *fw,psf *p)
{
 int		wwd,nvar;
 fits		*outfits;
 fitsimage	*out;

 wwd=p->grid*(2*p->hsize+1);
 nvar=(p->order+1)*(p->order+2)/2;

 outfits=fits_create();
 out=&outfits->i;

 out->dim=3;
 out->naxis[0]=wwd,
 out->naxis[1]=wwd,
 out->naxis[2]=nvar;
 out->vdata=(void *)p->coeff;
 out->bit=-32;
 fits_set_standard(outfits,NULL);/* SIMPLE */
 fits_set_image_params(outfits); /* NAXIS, NAXIS[123], BITPIX, BSCALE, BZERO */
 fits_set_header_integer(outfits,psfheaders[PSFHDR_PSFHSIZE],FITS_SH_FIRST,p->hsize,NULL);
 fits_set_header_integer(outfits,psfheaders[PSFHDR_PSFSGRID],FITS_SH_FIRST,p->grid ,NULL);
 fits_set_header_integer(outfits,psfheaders[PSFHDR_PSFORDER],FITS_SH_FIRST,p->order,NULL);
 fits_set_header_double (outfits,psfheaders[PSFHDR_PSFOFFSX],FITS_SH_FIRST,p->ox   ,NULL);
 fits_set_header_double (outfits,psfheaders[PSFHDR_PSFOFFSY],FITS_SH_FIRST,p->oy   ,NULL);
 fits_set_header_double (outfits,psfheaders[PSFHDR_PSFSCALE],FITS_SH_FIRST,p->scale,NULL);
 fits_set_origin(outfits,"fi/src/psf.c","PSF determination");

 fits_write(fw,outfits);

 fits_free(outfits);

 return(0);
}

int psf_parse_fits(fits *img,psf *p)
{
 int	k,wwd,nvar,order,hsize,grid,i,j;
 double	ox,oy,scale,w;

 if ( img==NULL )	return(1);

 if ( img->i.vdata==NULL || img->i.dim != 3 )	return(1);

 k=0;
 k|=fits_get_header_as_double(img,psfheaders[0],&w,0);hsize=(int)w;
 k|=fits_get_header_as_double(img,psfheaders[1],&w,0);grid =(int)w;
 k|=fits_get_header_as_double(img,psfheaders[2],&w,0);order=(int)w;
 k|=fits_get_header_as_double(img,psfheaders[3],&ox,0);
 k|=fits_get_header_as_double(img,psfheaders[4],&oy,0);
 k|=fits_get_header_as_double(img,psfheaders[5],&scale,0);
 if ( k || hsize<=0 || grid<=0 || order<0 )
	return(1);
 wwd=grid*(2*hsize+1);
 nvar=(order+1)*(order+2)/2;
 if ( img->i.naxis[0] != wwd || img->i.naxis[1] != wwd || img->i.naxis[2] != nvar )
	return(1);

 p->hsize=hsize;
 p->grid=grid;
 p->order=order;
 
 p->ox=ox,
 p->oy=oy;
 p->scale=scale;

 p->coeff=(double ***)tensor_alloc_3d(double,wwd,wwd,nvar);
 for ( i=0 ; i<wwd ; i++ )
  {	for ( j=0 ; j<wwd ; j++ )
	 {	for ( k=0 ; k<nvar ; k++ )
			p->coeff[k][i][j]=((double ***)img->i.vdata)[k][i][j];
	 }
  }

 return(0);
}

int psf_read_fits(FILE *fr,psf *p)
{
 fits	*img;

 img=fits_read(fr);
 if ( img==NULL )
	return(1);
 if ( psf_parse_fits(img,p) )
  {	fits_free(img);return(1);	}

 img->i.vdata=NULL;
 img->i.data=NULL;

 fits_free(img);
 
 return(0);
}

/*****************************************************************************/

