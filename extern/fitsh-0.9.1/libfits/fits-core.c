/*****************************************************************************/
/* fits-core.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Another simple standalone library for manipulating FITS files:	     */
/* Core procedures, functions and some wrappers for backward compatibility   */
/* (and of course for simpler usage...)					     */ 
/* (c) 2004-06, Pal, A. (apal@szofi.elte.hu). 				     */
/* See reference(s) at the end of this source code.			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in fits.h          */
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

int fits_tapeblock_size(int size)
{
 size=(size+FITS_TAPE_BLOCKSIZE-1)/FITS_TAPE_BLOCKSIZE;
 return(size*FITS_TAPE_BLOCKSIZE);
}

/*****************************************************************************/

int fits_alloc_image(fits *img,int sx,int sy)
{
 int		ret;
 ret=fits_image_alloc(&img->i,sx,sy);
 return(ret);
}
int fits_alloc_image_gen(fits *img,int dim,int *pnaxis)
{
 int		ret;
 ret=fits_image_alloc_gen(&img->i,dim,pnaxis);
 return(ret);
}

int fits_alloc_image_2d(fits *img,int sx,int sy)
{
 return ( fits_alloc_image(img,sx,sy) );
}
int fits_alloc_image_3d(fits *img,int sx,int sy,int sz)
{
 int	pnaxis[3];
 pnaxis[0]=sx,
 pnaxis[1]=sy,
 pnaxis[2]=sz;
 return ( fits_alloc_image_gen(img,3,pnaxis) );
}

/*****************************************************************************/

fits *fits_create(void)
{
 fits *ret;

 ret=(fits *)malloc(sizeof(fits));
 if ( ret==NULL )	return(NULL);

 fits_headerset_reset(&ret->header);	/* Empty header	*/

 ret->i.bit=0;			
 ret->i.sx=0;
 ret->i.sy=0;		/* No primary image	*/
 ret->i.data=NULL;

 ret->i.dim=0;		/* No primary image	*/
 ret->i.vdata=NULL;
 ret->i.allocdata=NULL;

 ret->i.read.bscale=ret->i.curr.bscale=1.0; /* default scaling parameters */
 ret->i.read.bzero =ret->i.curr.bzero =0.0; /* default scaling parameters */

 ret->length=0;
 ret->rawdata=NULL;	/* No raw data.	*/

 ret->xtns=NULL;	/* No extensions */
 ret->nxtn=0;	

 return(ret);
}

static fits *fits_duplicate_native(fits *img,int flag)
{
 fits		*ret;
 int		i;
 fitsimage	*si,*di;
 fitsttable	*st,*dt;
 fitsbtable	*sb,*db;

 if ( img==NULL )
  {	ret=fits_create();
	return(ret);
  }

 ret=(fits *)malloc(sizeof(fits));
 if ( ret==NULL )	return(NULL);

 fits_headerset_duplicate(&ret->header,&img->header);

 fits_image_duplicate(&ret->i,&img->i,flag);

 if ( img->rawdata != NULL )
  {	ret->length=img->length;
	ret->rawdata=(char *)malloc(img->length);
	if ( flag )	memcpy(ret->rawdata,img->rawdata,img->length);
	else		memset(ret->rawdata,0,img->length);
  }
 else
  {	ret->length=0;
	ret->rawdata=NULL;
  }

 if ( img->nxtn>0 && img->xtns != NULL )
  {	ret->nxtn=img->nxtn;
	ret->xtns=(fitsextension *)malloc(sizeof(fitsextension)*ret->nxtn);
	for ( i=0 ; i<img->nxtn ; i++ )
	 {	ret->xtns[i].type=img->xtns[i].type;
		fits_headerset_duplicate(&ret->xtns[i].header,&img->xtns[i].header);
		switch ( ret->xtns[i].type )
		 {   case FITS_EXT_IMAGE:
			si=&img->xtns[i].x.i;
			di=&ret->xtns[i].x.i;
			fits_image_duplicate(di,si,flag);
			break;
		     case FITS_EXT_TABLE:
			st=&img->xtns[i].x.t;
			dt=&ret->xtns[i].x.t;
			fits_table_duplicate(dt,st,flag);
			break;
		     case FITS_EXT_BINTABLE:
			sb=&img->xtns[i].x.b;
			db=&ret->xtns[i].x.b;
			fits_bintable_duplicate(db,sb,flag);
			break;
		 }
	 }
  }
 else
  {	ret->nxtn=0;
	ret->xtns=NULL;
  }

 return(ret);
}

fits *fits_duplicate(fits *img)
{
 fits	*ret;
 ret=fits_duplicate_native(img,1);
 return(ret);
}

fits *fits_duplicate_empty(fits *img)
{
 fits	*ret;
 ret=fits_duplicate_native(img,0);
 return(ret);
}

void fits_free_image(fits *img)
{
 fits_image_free(&img->i);
}

void fits_free(fits *img)
{
 int	i;

 fits_headerset_free(&img->header);

 fits_free_image(img);

 if ( img->rawdata != NULL )
  {	free(img->rawdata);
	img->length=0;
  }

 for ( i=0 ; i<img->nxtn && img->xtns != NULL ; i++ )
  {	fits_headerset_free(&img->xtns[i].header);
	switch ( img->xtns[i].type )
	 {   case FITS_EXT_IMAGE:
		fits_image_free(&img->xtns[i].x.i);
		break;
	     case FITS_EXT_TABLE:
		fits_table_free(&img->xtns[i].x.t);
		break;
	     case FITS_EXT_BINTABLE:
		fits_bintable_free(&img->xtns[i].x.b);
		break;
	 }
  }
 if ( img->xtns != NULL )	free(img->xtns);

 free(img);
}

/*****************************************************************************/

int fits_read_header(FILE *fr,fits *img)
{
 int	ret;

 ret=fits_headerset_read_cb(fits_cb_read,(void *)fr,&img->header);

 return(ret);
}

/*****************************************************************************/
 
int fits_get_header_count(fits *img,char *hdr)
{
 int	n;
 n=fits_headerset_get_count(&img->header,hdr);
 return(n);
}

int fits_get_header_id(fits *img,char *hdr,int cnt)
{
 int	i;
 i=fits_headerset_get_id(&img->header,hdr,cnt);
 return(i);
}

fitsheader * fits_get_header(fits *img,char *hdr,int cnt)
{
 fitsheader	*h;
 h=fits_headerset_get_header(&img->header,hdr,cnt);
 return(h);
}

int fits_get_header_as_double(fits *img,char *hdr,double *ret,int is_ambigous_allowed)
{
 int	r;
 r=fits_headerset_get_as_double(&img->header,hdr,ret,is_ambigous_allowed);
 return(r);
}

int fits_get_gain(fits *img,double *ret)
{
 int	r;
 r=fits_get_header_as_double(img,headers[HDR_GAIN],ret,0);
 return(r);
}

/*****************************************************************************/

fitsheader *fits_set_header_any(fits *img,char *hdr,int rule,char *comment)
{
 fitsheader	*hd;
 hd=fits_headerset_set_any(&img->header,hdr,rule,comment);
 return(hd);
}

int fits_set_header_integer(fits *img,char *hdr,int rule,int val,char *comment)
{
 int	ret;
 ret=fits_headerset_set_integer(&img->header,hdr,rule,val,comment);
 return(ret);
}
int fits_set_header_double(fits *img,char *hdr,int rule,double val,char *comment)
{
 int	ret;
 ret=fits_headerset_set_double(&img->header,hdr,rule,val,comment);
 return(ret);
}
int fits_set_header_string(fits *img,char *hdr,int rule,char *str,char *comment)
{
 int	ret;
 ret=fits_headerset_set_string(&img->header,hdr,rule,str,comment);
 return(ret);
}
int fits_set_header_boolean(fits *img,char *hdr,int rule,int vbool,char *comment)
{
 int	ret;
 ret=fits_headerset_set_boolean(&img->header,hdr,rule,vbool,comment);
 return(ret);
}

/*****************************************************************************/

int fits_delete_header(fits *img,char *hdr,int k)
{
 int	ret;
 ret=fits_headerset_delete(&img->header,hdr,k);
 return(ret);
}
int fits_delete_all_header(fits *img,char *hdr)
{
 int	ret;
 ret=fits_headerset_delete_all(&img->header,hdr);
 return(ret);
}

/*****************************************************************************/

void fits_copy_full_header(fits *im1,fits *im2)
{
 fits_headerset_copy(&im1->header,&im2->header);
}

/*****************************************************************************/

int fits_set_image(fits *img,double value)
{
 fits_image_set_value(&img->i,value);
 return(0);
}
int fits_reset_image(fits *img)
{
 fits_image_reset(&img->i);
 return(0);
}

/*****************************************************************************/

int fits_read_image_line(FILE *f,int sx,int bit,double *line)
{
 int	ret;
 ret=fits_image_read_line(f,sx,bit,line);
 return(ret);
}

/*****************************************************************************/

int fits_get_scale(fits *img,double *rbscale,double *rbzero)
{
 int	ret;
 ret=fits_image_get_scale(&img->header,&img->i,rbscale,rbzero);
 return(ret);
}

int fits_get_image_params(fits *img)
{
 int	ret;
 ret=fits_image_get_params(&img->header,&img->i);
 return(ret);
}

/*****************************************************************************/

int fits_rescale(fits *img)
{
 int	ret;
 ret=fits_image_rescale(&img->i);
 return(ret);
}

int fits_backscale(fits *img,double ibscale,double ibzero)
{
 int	ret;
 ret=fits_image_backscale(&img->i,ibscale,ibzero);
 return(ret);
}

int fits_set_image_params(fits *img)
{
 int	ret;
 ret=fits_image_set_params(&img->header,&img->i);
 return(ret);
}

/*****************************************************************************/

int fits_set_standard(fits *img,char *comment)
{
 if ( comment==NULL )	comment="FITS standard";
 fits_headerset_set_boolean(&img->header,headers[HDR_SIMPLE],
	FITS_SH_FORCEFIRST,1,comment);
 return(0);
}
int fits_set_origin(fits *img,char *origin,char *comment)
{
 fits_headerset_set_string(&img->header,headers[HDR_ORIGIN],
	FITS_SH_FIRST,origin,comment);
 return(0);
}
int fits_set_extend(fits *img,int flag,char *comment)
{
 fits_headerset_set_boolean(&img->header,headers[HDR_EXTEND],
	FITS_SH_ADD,flag,comment);
 return(0);
}

/*****************************************************************************/

int fits_read_image(FILE *fr,fits *img)
{
 int	ret;
 ret=fits_image_read(fr,&img->i);
 return(ret);
}

int fits_read_rawdata(FILE *fr,fits *img)
{
 int	rd,length;
 char	*data;

 data=(char *)malloc(BUFSIZE);
 length=0;
 while ( ! feof(fr) )
  {	rd=fread(data+length,1,BUFSIZE,fr);
	if ( rd==0 )	break;
	length+=rd;
	if ( ! feof(fr) )
		data=(char *)realloc(data,length+BUFSIZE);
  }
 data=(char *)realloc(data,length);
 
 img->i.data=NULL;
 img->i.sx=img->i.sy=0;
 img->i.dim=0;
 img->i.vdata=NULL;
 
 img->rawdata=data,
 img->length=length;

 return(0);
}

/*****************************************************************************/

fitsextension *fits_extension_add(fits *img,int cnt)
{
 fitsextension	*fx;

 img->xtns=(fitsextension *)realloc(img->xtns,sizeof(fitsextension)*(img->nxtn+cnt));

 fx=&img->xtns[img->nxtn]; 
 memset(fx,0,sizeof(fitsextension)*cnt);
 img->nxtn+=cnt;
 return(fx);

}
fitsextension *fits_extension_new(fits *img,int type)
{
 fitsextension	*fx;
 fx=fits_extension_add(img,1);
 fx->type=type;
 return(fx);
}

/*****************************************************************************/

int fits_read_extension_image_cb(int (*cb_read)(void *,void *,int),void *param,fitsextension *xtn)
{
 int	ret;

 ret=fits_image_get_params(&xtn->header,&xtn->x.i);
 if ( ret )	return(ret);
 ret=fits_image_alloc_gen(&xtn->x.i,xtn->x.i.dim,xtn->x.i.naxis);
 if ( ret )	return(ret);
 ret=fits_image_read_cb(cb_read,param,&xtn->x.i);
 return(ret);
}
int fits_read_extension_image(FILE *fr,fitsextension *xtn)
{
 return(fits_read_extension_image_cb(fits_cb_read,(void *)fr,xtn));
}

int fits_read_extension_table_cb(int (*cb_read)(void *,void *,int),void *param,fitsextension *xtn)
{
 int	ret;

 ret=fits_table_get_params(&xtn->header,&xtn->x.t);
 if ( ret )	return(ret);
 ret=fits_table_alloc(&xtn->x.t);
 if ( ret )	return(ret);
 ret=fits_table_read_cb(cb_read,param,&xtn->x.t);
 return(ret);
}
int fits_read_extension_table(FILE *fr,fitsextension *xtn)
{
 return(fits_read_extension_table_cb(fits_cb_read,(void *)fr,xtn));
}

int fits_read_extension_bintable_cb(int (*cb_read)(void *,void *,int),void *param,fitsextension *xtn)
{
 int	ret;

 ret=fits_bintable_get_params(&xtn->header,&xtn->x.b);
 if ( ret )	return(ret);
 ret=fits_bintable_alloc(&xtn->x.b);
 if ( ret )	return(ret);
 ret=fits_bintable_read_cb(cb_read,param,&xtn->x.b);
 return(ret);
}
int fits_read_extension_bintable(FILE *fr,fitsextension *xtn)
{
 return(fits_read_extension_bintable_cb(fits_cb_read,(void *)fr,xtn));
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int fits_skip_extension_image_cb(int (*cb_read)(void *,void *,int),void *param,fitsextension *xtn)
{
 int	ret;

 ret=fits_image_get_params(&xtn->header,&xtn->x.i);
 if ( ret )	return(ret);
 ret=fits_image_skip_cb(cb_read,param,&xtn->x.i);
 return(ret);
}

int fits_skip_extension_table_cb(int (*cb_read)(void *,void *,int),void *param,fitsextension *xtn)
{
 int	ret;

 ret=fits_table_get_params(&xtn->header,&xtn->x.t);
 if ( ret )	return(ret);
 ret=fits_table_skip_cb(cb_read,param,&xtn->x.t);
 return(ret);
}

int fits_skip_extension_bintable_cb(int (*cb_read)(void *,void *,int),void *param,fitsextension *xtn)
{
 int	ret;

 ret=fits_bintable_get_params(&xtn->header,&xtn->x.b);
 if ( ret )	return(ret);
 ret=fits_bintable_skip_cb(cb_read,param,&xtn->x.b);
 return(ret);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int fits_read_extensions_cb(int (*cb_read)(void *,void *,int),void *param,fits *img)
{
 fitsextension	xtn;
 int		type;

 while ( 1 )
  {	xtn.type=0;
	fits_headerset_reset(&xtn.header);
	fits_headerset_read_cb(cb_read,param,&xtn.header);

	if ( xtn.header.nhdr<=0 )	break;

	type=fits_headerset_is_extension(&xtn.header);
	xtn.type=type;

	switch ( type )
	 {   case FITS_EXT_IMAGE:
		fits_read_extension_image_cb(cb_read,param,&xtn);
		break;
	     case FITS_EXT_TABLE:
		fits_read_extension_table_cb(cb_read,param,&xtn);
		break;
	     case FITS_EXT_BINTABLE:
		fits_read_extension_bintable_cb(cb_read,param,&xtn);
		break;
	     default:
		fits_headerset_free(&xtn.header);
		return(1);
		break;
	 }

	img->xtns=(fitsextension *)realloc(img->xtns,sizeof(fitsextension)*(img->nxtn+1));
	memcpy(&img->xtns[img->nxtn],&xtn,sizeof(fitsextension));
	img->nxtn++;
  };

 return(0);
}

fits *fits_read_cb(int (*cb_read)(void *,void *,int),void *param)
{
 fits		*img;
 int		i;
 fitsimage	*fi;
 fitsheader	*hx;

 img=fits_create();
 if ( img==NULL )	return(NULL);

 fits_headerset_reset(&img->header);
 fits_headerset_read_cb(cb_read,param,&img->header);

 fi=&img->i;

 i=fits_image_get_params(&img->header,fi);
 if ( ! i )
  {	i=fits_image_alloc_gen(fi,fi->dim,fi->naxis);
	if ( i )
	 {	fits_free(img);
		return(NULL);
	 }
	fits_image_read_cb(cb_read,param,fi);
  }

/***
 else
  {	i=fits_read_rawdata(fr,img);
	if ( i )
	 {	fits_free(img);
		return(NULL);
	 }
  }
***/

 hx=fits_headerset_get_uniq_header(&img->header,headers[HDR_EXTEND]);
 if ( hx != NULL && hx->vtype==FITS_VBOOLEAN && hx->vint )
  {	fits_read_extensions_cb(cb_read,param,img);		}

 return(img);
}
fits *fits_read(FILE *fr)
{
 return(fits_read_cb(fits_cb_read,(void *)fr));
}

/*****************************************************************************/

int fits_skip_more_extensions_cb(int (*cb_read)(void *,void *,int),void *param,int nframe)
{
 fitsextension	xtn;
 int		type;
 
 while ( nframe>0 )
  {	
	fits_headerset_reset(&xtn.header);
	fits_headerset_read_cb(cb_read,param,&xtn.header);
	if ( xtn.header.nhdr<=0 )	break;

	type=fits_headerset_is_extension(&xtn.header);
	switch ( type )
	 {   case FITS_EXT_IMAGE:
		fits_skip_extension_image_cb(cb_read,param,&xtn);
		break;
	     case FITS_EXT_TABLE:
		fits_skip_extension_table_cb(cb_read,param,&xtn);
		break;
	     case FITS_EXT_BINTABLE:
		fits_skip_extension_bintable_cb(cb_read,param,&xtn);
		break;
	     default:
		fits_headerset_free(&xtn.header);
		break;
	 }

	nframe--;
  };

 return(0);
}
int fits_skip_more_extensions(FILE *fr,int nframe)
{
 return(fits_skip_more_extensions_cb(fits_cb_read,(void *)fr,nframe));
}

fits * fits_read_frame_as_extension_cb(int (*cb_read)(void *,void *,int),void *param,int frameno)
{
 fits		*img;
 int		i,type;
 fitsheaderset	header;
 fitsimage	priimg;
 fitsextension	*fx,xtn;

 if ( frameno<0 )	return(NULL);	/* invalid frame number */
 if ( cb_read==NULL )	return(NULL);	/* invalid callback     */

 img=fits_create();

 fits_headerset_reset(&header);
 fits_headerset_read_cb(cb_read,param,&header);
 i=fits_image_get_params(&header,&priimg);

 if ( ! i )	/* valid image */
  {	if ( frameno==0 )
	 {	fits_image_alloc_gen(&priimg,priimg.dim,priimg.naxis);
		fits_image_read_cb(cb_read,param,&priimg);
		fx=fits_extension_add(img,1);
		memcpy(&fx->header,&header,sizeof(fitsheaderset));
		fx->type=FITS_EXT_IMAGE;
		memcpy(&fx->x.i,&priimg,sizeof(fitsimage));
		return(img);
	 }
	else
		fits_image_skip_cb(cb_read,param,&priimg);

	frameno--;
  }

 fits_headerset_free(&header);

 fits_skip_more_extensions_cb(cb_read,param,frameno);

 fits_headerset_reset(&xtn.header);
 fits_headerset_read_cb(cb_read,param,&xtn.header);

 type=fits_headerset_is_extension(&xtn.header);
 switch ( type )
  {  case FITS_EXT_IMAGE:
	fits_read_extension_image_cb(cb_read,param,&xtn);
	break;
     case FITS_EXT_TABLE:
	fits_read_extension_table_cb(cb_read,param,&xtn);
	break;
     case FITS_EXT_BINTABLE:
	fits_read_extension_bintable_cb(cb_read,param,&xtn);
	break;
     default:
	fits_headerset_free(&xtn.header);
	return(NULL);
	break;
  }
 xtn.type=type;

 fx=fits_extension_add(img,1);
 memcpy(fx,&xtn,sizeof(fitsextension));

 return(img);
}
fits * fits_read_frame_as_extension(FILE *fr,int frameno)
{
 return(fits_read_frame_as_extension_cb(fits_cb_read,(void *)fr,frameno));
}

fits * fits_seek_frame_to_image_cb(int (*cb_read)(void *,void *,int),void *param,int frameno)
{
 fits		*img;
 int		i,type;
 fitsheaderset	header;

 if ( frameno<0 )	return(NULL);	/* invalid frame number */

 img=fits_create();
 fits_headerset_reset(&img->header);
 fits_headerset_read_cb(cb_read,param,&img->header);
 if ( img->header.nhdr<=0 || img->header.hdrs==NULL )
  {	fits_free(img);
	return(0);
  }
 i=fits_image_get_params(&img->header,&img->i);

 if ( ! i )	/* valid image */
  {	if ( frameno==0 )
		return(img);
	else
	 {	fits_image_skip_cb(cb_read,param,&img->i);
		frameno--;
	 }
  }

 fits_skip_more_extensions_cb(cb_read,param,frameno);

 fits_headerset_reset(&header);
 fits_headerset_read_cb(cb_read,param,&header);

 type=fits_headerset_is_extension(&header);

 if ( type==FITS_EXT_IMAGE )
  {	fits_image_get_params(&header,&img->i);
	fits_headerset_merge(&img->header,&header,0);
	return(img);
  }
 else
  {	fits_headerset_free(&header);
	fits_headerset_free(&img->header);
	fits_free(img);
	return(NULL);
  }

 return(img);
}
fits * fits_seek_frame_to_image(FILE *fr,int frameno)
{
 return(fits_seek_frame_to_image_cb(fits_cb_read,(void *)fr,frameno));
}

fits * fits_read_frame_to_image_cb(int (*cb_read)(void *,void *,int),void *param,int frameno)
{
 fits	*img;

 if ( frameno<0 )	frameno=0;

 img=fits_seek_frame_to_image_cb(cb_read,param,frameno);
 if ( img==NULL )	return(NULL);

 fits_image_alloc_gen(&img->i,img->i.dim,img->i.naxis);
 fits_image_read_cb(cb_read,param,&img->i);

 return(img);
}
fits * fits_read_frame_to_image(FILE *fr,int frameno)
{
 return(fits_read_frame_to_image_cb(fits_cb_read,(void *)fr,frameno));
}

/*****************************************************************************/

fits *fits_read_raw(FILE *fr)
{
 fits	*img;
 int	i;

 img=fits_create();
 if ( img==NULL )	return(NULL);

 fits_headerset_reset(&img->header);
 fits_headerset_read(fr,&img->header);

 i=fits_read_rawdata(fr,img);
 if ( i )
  {	fits_free(img);
	return(NULL);
  }

 return(img);
}

/*****************************************************************************/

int fits_write_header(FILE *fw,fits *img)
{
 int	ret;
 ret=fits_headerset_write(fw,&img->header);
 return(ret);
}

/*****************************************************************************/

int fits_write_image_line(FILE *fw,int sx,int bit,double *line)
{
 int	ret;
 ret=fits_image_write_line(fw,sx,bit,line);
 return(ret);
}
int fits_write_image(FILE *fw,fits *img)
{
 int	wr;
 wr=fits_image_write(fw,&img->i,0);
 return(wr);		/* total bytes written */
}

int fits_write_cb(int (*cb_write)(void *,void *,int),void *param,fits *img)
{
 int		wr,wn;
 char		*wbuff;
 int		i,j,validimg;

 wr=fits_headerset_write_cb(cb_write,param,&img->header);

 validimg=1; 
 for ( i=0 ; i<img->i.dim ; i++ )
  {	if ( img->i.naxis[i] <= 0 )	validimg=0;	}

 if ( img->i.vdata != NULL && img->i.dim>0 && validimg )
  {	wr+=fits_image_write_cb(cb_write,param,&img->i,1);	/* pad image */		}

 for ( i=0 ; i<img->nxtn && img->xtns != NULL ; i++ )
  {	wr+=fits_headerset_write_cb(cb_write,param,&img->xtns[i].header);
	switch ( img->xtns[i].type )
	 {   case FITS_EXT_IMAGE:
		wr+=fits_image_write_cb(cb_write,param,&img->xtns[i].x.i,1);
		break;
	     case FITS_EXT_TABLE:
		wr+=fits_table_write_cb(cb_write,param,&img->xtns[i].x.t,1);
		break;
	     case FITS_EXT_BINTABLE:
		wr+=fits_bintable_write_cb(cb_write,param,&img->xtns[i].x.b,1);
		break;
	 }
  }

 if ( img->rawdata != NULL && img->length>0 )
  {	for ( i=0 ; i<img->length ;  )
	 {	j=BUFSIZE;
		if ( img->length-i<BUFSIZE ) j=img->length-i;
		cb_write(param,img->rawdata+i,j);
		wr+=j;i+=j;
	 }
	wn=FITS_TAPE_BLOCKSIZE-(wr%FITS_TAPE_BLOCKSIZE);
	if ( wn==FITS_TAPE_BLOCKSIZE )	return(wr);
	wbuff=(char *)malloc(wn);
	memset(wbuff,0,wn);
	cb_write(param,wbuff,wn);
	wr+=wn;
	free(wbuff);
  }

 return(wr);
}
int fits_write(FILE *fw,fits *img)
{
 return(fits_write_cb(fits_cb_write,(void *)fw,img));
}

int fits_mem_write(void **rbuff,fits *img)
{
 fitsmemwrite	fmw;
 int		length;

 if ( rbuff==NULL )	return(0);
 
 fmw.buffer=NULL;
 fmw.length=0;
 length=fits_write_cb(fits_cb_mem_write,(void *)(&fmw),img);
 
 *rbuff=(void *)fmw.buffer;
 return(length);
}

/*****************************************************************************/

double * fits_dump_image_raw(fits *img)
{
 double	*ret;
 ret=fits_image_dump(&img->i);
 return(ret);
}

/******************************************************************************
 Reference:
  [1]:	Definition of the Flexible Image Transport System (FITS)
	NASA / Science Office of Standards and Technology (NOST)
	NOST 100-2.0 (March 29, 1999)
	NASA Goddard Space Flight Center
	http://fits.gsfc.nasa.gov/fits_documentation.html
******************************************************************************/
       
