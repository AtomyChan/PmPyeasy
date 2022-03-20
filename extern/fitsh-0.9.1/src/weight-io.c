/*****************************************************************************/
/* weight-io.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions related to weight stamp handling.				     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2006; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <stdarg.h>

#include <fits/fits.h>

#include "fi.h"

#include "tensor.h"
#include "weight.h"

/*****************************************************************************/

#define		HDR_WGTDATA		0
#define		HDR_WGTTABLE		1
#define		HDR_WGTHSIZE		2
#define		HDR_WGTSGRID		3

static	char *wgthdr[]=
 {	"WGTDATA",
	"WGTTABLE",
	"WGTHSIZE",
	"WGTSGRID",
	NULL
 };

/*****************************************************************************/

fits *weight_fits_create(weightlist *wl)
{
 int		fsize,hsize,grid,n;
 fits		*img;
 fitsextension	*fx;
 unsigned char	*line;
 weight		*ww;

 img=fits_create();
 fits_set_standard(img,NULL);

 hsize=wl->hsize;
 grid =wl->grid;
 fsize=grid*(2*hsize+1);

 img->i.dim=3;
 img->i.sx=img->i.naxis[0]=fsize;
 img->i.sy=img->i.naxis[1]=fsize;
 img->i.naxis[2]=wl->nweight;
 img->i.vdata=wl->zdata;
 img->i.data =wl->zdata[0];
 img->i.allocdata=NULL;

 img->i.bit=-32;
 img->i.curr.bscale=1.0,img->i.curr.bzero=0.0;
 img->i.read.bscale=1.0,img->i.read.bzero=0.0;
 fits_set_image_params(img);
 fits_set_extend(img,1,NULL);
 fits_set_header_boolean(img,wgthdr[HDR_WGTDATA],FITS_SH_FIRST,1,"Weight metadata is present as a BINTABLE, ...");
 fits_set_header_integer(img,wgthdr[HDR_WGTTABLE],FITS_SH_FIRST,2,"... in sec [2] (1st extension after pri data)");
 fits_set_header_integer(img,wgthdr[HDR_WGTHSIZE],FITS_SH_FIRST,hsize,"Half-size of the stamps");
 fits_set_header_integer(img,wgthdr[HDR_WGTSGRID],FITS_SH_FIRST,grid,"Subpixel grid resolution of the stamps");

 fx=fits_extension_new(img,FITS_EXT_BINTABLE);
 fits_bintable_create_fields(&fx->x.b,wl->nweight,5,FTF_SHORT,1,FTF_SHORT,1,FTF_FLOAT,1,FTF_FLOAT,1,FTF_FLOAT,1);
 fits_bintable_set_xtr_params(&fx->x.b,0,"IX"  ,"pixel",NULL);
 fits_bintable_set_xtr_params(&fx->x.b,1,"IY"  ,"pixel",NULL);
 fits_bintable_set_xtr_params(&fx->x.b,2,"X"   ,"pixel",NULL);
 fits_bintable_set_xtr_params(&fx->x.b,3,"Y"   ,"pixel",NULL);
 fits_bintable_set_xtr_params(&fx->x.b,4,"flux","adu"  ,NULL);

 fits_bintable_alloc(&fx->x.b);
 fits_bintable_set_params(&fx->header,&fx->x.b);
 fits_headerset_set_string(&fx->header,"EXTNAME",FITS_SH_FIRST,"WEIGHTSTAMPS",NULL);

 for ( n=0 ; n<wl->nweight ; n++ )
  {	line=fx->x.b.data[n];
	ww=&wl->weights[n];
	*(short *)(line+fx->x.b.bfields[0].offset)=ww->ix;
	*(short *)(line+fx->x.b.bfields[1].offset)=ww->iy;
	*(float *)(line+fx->x.b.bfields[2].offset)=ww->x;
	*(float *)(line+fx->x.b.bfields[3].offset)=ww->y;
	*(float *)(line+fx->x.b.bfields[4].offset)=ww->flux;
  }

 return(img);
}	

/*****************************************************************************/

int weight_parse_fits(fits *img,weightlist *wl)
{
 fitsheader	*hd;
 fitsextension	*fx;
 fitsbtable	*fb;
 fitsimage	*fi;
 weight		*ww;
 int		i,k,l,hsize,grid,nweight,fsize;
 unsigned char	*line;
 double		***zdata;

 if ( wl==NULL || img==NULL )	return(0);

 wl->zdata=NULL;
 wl->weights=NULL;
 wl->nweight=0;
 
 hd=fits_headerset_get_uniq_header(&img->header,wgthdr[HDR_WGTDATA]);
 if ( hd==NULL || hd->vtype != FITS_VBOOLEAN || ! hd->vint )	return(1);
 hd=fits_headerset_get_uniq_header(&img->header,wgthdr[HDR_WGTTABLE]);
 if ( hd==NULL || hd->vtype != FITS_VINT || hd->vint<2 )	return(2);
 i=hd->vint-2;
 if ( i>=img->nxtn )			return(3);
 fx=&img->xtns[i];
 if ( fx->type != FITS_EXT_BINTABLE )	return(4);

 hd=fits_headerset_get_uniq_header(&img->header,wgthdr[HDR_WGTHSIZE]);
 if ( hd==NULL || hd->vtype != FITS_VINT || hd->vint<0 )	return(5);
 hsize=hd->vint;
 hd=fits_headerset_get_uniq_header(&img->header,wgthdr[HDR_WGTSGRID]);
 if ( hd==NULL || hd->vtype != FITS_VINT || hd->vint<=0 )	return(6);
 grid =hd->vint;

 fb=&fx->x.b;
 fi=&img->i;
 
 if ( fb->rowsize != 16 || fb->nrow<=0 )	return(7);
 i=fits_bintable_check_fields(fb,5,FTF_SHORT,1,FTF_SHORT,1,FTF_FLOAT,1,FTF_FLOAT,1,FTF_FLOAT,1);
 if ( i )	return(8);
 
 nweight=fb->nrow;
 fsize=grid*(2*hsize+1);
 if ( fi->dim != 3 || fi->vdata==NULL )	return(9);
 if ( fi->naxis[0] != fsize || fi->naxis[1] != fsize || fi->naxis[2] != nweight )	return(10);
 zdata=(double ***)fi->vdata;

 wl->zdata=(double ***)tensor_alloc_3d(double,fsize,fsize,nweight);
 wl->hsize=hsize; 
 wl->grid =grid;
 wl->nweight=nweight;
 wl->weights=(weight *)malloc(sizeof(weight)*nweight);

 for ( i=0 ; i<nweight ; i++ )
  {	ww=&wl->weights[i];
	line=fb->data[i];

	ww->ix	=*(short *)(line+fb->bfields[0].offset);
	ww->iy	=*(short *)(line+fb->bfields[1].offset);
	ww->x	=*(float *)(line+fb->bfields[2].offset);
	ww->y   =*(float *)(line+fb->bfields[3].offset);
	ww->flux=*(float *)(line+fb->bfields[4].offset);

	ww->iarr=wl->zdata[i];
	for ( k=0 ; k<fsize ; k++ )
	 {	for ( l=0 ; l<fsize ; l++ )
		 {	ww->iarr[k][l]=zdata[i][k][l];		}
	 }
  }

 return(0);
}

/*****************************************************************************/

