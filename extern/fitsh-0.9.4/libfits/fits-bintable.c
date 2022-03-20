/*****************************************************************************/
/* fits-bintable.c							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Another simple standalone library for manipulating FITS files:	     */
/* Binary table extensions (XTENSION = 'BINTABLE')			     */
/* (c) 2006, Pal, A. (apal@szofi.elte.hu). 				     */
/* See reference(s) at the end of this source code.			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in fits.h	     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>

#define	_FITS_SOURCE

#include <fits/fits.h>
#include "fits-common.h"

/*****************************************************************************/

int fits_bintable_form_basesize(int form)
{
 int	basesize;

 switch ( form )	/* basesize = 0 means: basesize is a bit, 1/8 byte */
  {	case FTF_LOGICAL	: basesize= 1;break;
	case FTF_BIT		: basesize= 0;break; 	
	case FTF_BYTE		: basesize= 1;break; 
	case FTF_SHORT		: basesize= 2;break; 
	case FTF_LONG		: basesize= 4;break; 
	case FTF_LONGLONG	: basesize= 8;break;
	case FTF_CHAR		: basesize= 1;break; 
	case FTF_FLOAT		: basesize= 4;break; 
	case FTF_DOUBLE		: basesize= 8;break; 
	case FTF_FLOATCOMPLEX	: basesize= 8;break; 
	case FTF_DOUBLECOMPLEX	: basesize=16;break; 
	case FTF_ARRAY		: basesize= 8;break; 
	default : basesize=-1;break;
  }
 return(basesize);
}

char *fits_bintable_form_cname(int form)
{
 char	*ret;

 switch ( form )
  {	case FTF_LOGICAL	: ret="logical";	break;
	case FTF_BIT		: ret="bit";		break; 	
	case FTF_BYTE		: ret="byte";		break; 
	case FTF_SHORT		: ret="short";		break; 
	case FTF_LONG		: ret="long";		break; 
	case FTF_LONGLONG	: ret="long long";	break; 
	case FTF_CHAR		: ret="char*";		break; 
	case FTF_FLOAT		: ret="float";		break; 
	case FTF_DOUBLE		: ret="double";		break; 
	case FTF_FLOATCOMPLEX	: ret="float complex";	break; 
	case FTF_DOUBLECOMPLEX	: ret="double complex";	break; 
	case FTF_ARRAY		: ret="array";		break;
	default : ret=NULL; break;
  }
 
 return(ret);
}



/*****************************************************************************/

int fits_bintable_get_params(fitsheaderset *header,fitsbtable *fb)
{
 fitsheader	*hdr;
 fitsbfield	*tb;
 int		i,repeat,form,basesize,offset;
 char		chdr[16],c;

 memset(fb,0,sizeof(fitsbtable));

 hdr=fits_headerset_get_uniq_header(header,headers[HDR_NAXIS]);
 if ( hdr==NULL || hdr->vtype != FITS_VINT || hdr->vint != 2 )	return(1);
 hdr=fits_headerset_get_uniq_header(header,headers[HDR_BITPIX]);
 if ( hdr==NULL || hdr->vtype != FITS_VINT || hdr->vint != 8 )	return(1);
 hdr=fits_headerset_get_uniq_header(header,headers[HDR_NAXIS1]);
 if ( hdr==NULL || hdr->vtype != FITS_VINT || hdr->vint<=0 )	return(1);
 fb->rowsize=hdr->vint;
 hdr=fits_headerset_get_uniq_header(header,headers[HDR_NAXIS2]);
 if ( hdr==NULL || hdr->vtype != FITS_VINT || hdr->vint<=0 )	return(1);
 fb->nrow=hdr->vint;
 hdr=fits_headerset_get_uniq_header(header,headers[HDR_TFIELDS]);
 if ( hdr==NULL || hdr->vtype != FITS_VINT || hdr->vint<=0 )	return(1);
 fb->nbfield=hdr->vint;

 fb->bfields=(fitsbfield *)malloc(sizeof(fitsbfield)*fb->nbfield);
 
 offset=0;
 for ( i=0 ; i<fb->nbfield ; i++ )
  {	tb=&fb->bfields[i];

	sprintf(chdr,"%s%d",headers[HDR_TFORM],i+1);
	hdr=fits_headerset_get_uniq_header(header,chdr);
	if ( hdr==NULL || hdr->vtype != FITS_VSTR )
	 {	free(fb->bfields);return(1);				}
	if ( isdigit((int)hdr->vstr[0]) )
	 {	if ( sscanf(hdr->vstr,"%d%c",&repeat,&c)<2 )
		 {	free(fb->bfields);return(1);			}
	 }
	else
	 {	c=toupper((int)hdr->vstr[0]);
		repeat=1;
	 }
	form=(int)c;
	basesize=fits_bintable_form_basesize(form);
	if ( basesize<0 )
	 {	free(fb->bfields);return(1);				}
	tb->form=form;
	tb->basesize=basesize;
	tb->repeat=repeat;
	tb->offset=offset;

	if ( basesize>0 )	offset+=basesize*repeat;
	else			offset+=(repeat+7)/8;

	tb->scale.bscale=1.0;
	tb->scale.bzero =0.0;
	tb->null[0]=0;
	tb->type[0]=0;
	tb->unit[0]=0;

	sprintf(chdr,"%s%d",headers[HDR_TSCAL],i+1);
	hdr=fits_headerset_get_uniq_header(header,chdr);
	if ( hdr != NULL )
	 {	if ( hdr->vtype==FITS_VINT )	tb->scale.bscale=hdr->vint;
		if ( hdr->vtype==FITS_VDOUBLE )	tb->scale.bscale=hdr->vdouble;
	 }
	sprintf(chdr,"%s%d",headers[HDR_TZERO],i+1);
	hdr=fits_headerset_get_uniq_header(header,chdr);
	if ( hdr != NULL )
	 {	if ( hdr->vtype==FITS_VINT )	tb->scale.bzero=hdr->vint;
		if ( hdr->vtype==FITS_VDOUBLE )	tb->scale.bzero=hdr->vdouble;
	 }
	sprintf(chdr,"%s%d",headers[HDR_TNULL],i+1);
	hdr=fits_headerset_get_uniq_header(header,chdr);
	if ( hdr != NULL && hdr->vtype==FITS_VSTR )
	 {	strncpy(tb->null,hdr->vstr,31);tb->null[31]=0;		}
	sprintf(chdr,"%s%d",headers[HDR_TTYPE],i+1);
	hdr=fits_headerset_get_uniq_header(header,chdr);
	if ( hdr != NULL && hdr->vtype==FITS_VSTR )
	 {	strncpy(tb->type,hdr->vstr,31);tb->type[31]=0;		}
	sprintf(chdr,"%s%d",headers[HDR_TUNIT],i+1);
	hdr=fits_headerset_get_uniq_header(header,chdr);
	if ( hdr != NULL && hdr->vtype==FITS_VSTR )
	 {	strncpy(tb->unit,hdr->vstr,31);tb->unit[31]=0;		}
	
  }

 return(0);
}

int fits_bintable_create_fields(fitsbtable *fb,int nrow,int nbfield,...)
{
 va_list	ap;
 int		i,repeat,form,basesize,offset;

 /* check whether the field descriptor array is initialized or not: */
 if ( fb->bfields != NULL || fb->nbfield != 0 )	return(1);

 if ( nbfield<=0 )	return(0);

 va_start(ap,nbfield);

 fb->bfields=(fitsbfield *)malloc(sizeof(fitsbfield)*nbfield);
 fb->nbfield=nbfield;

 offset=0;
 for ( i=0 ; i<nbfield && offset>=0 ; i++ )
  {	form=va_arg(ap,int);
	repeat=va_arg(ap,int);
	basesize=fits_bintable_form_basesize(form);
	if ( basesize<0 || repeat<=0 )
 	 {	offset=-1;break;	}
	fb->bfields[i].form=form;
	fb->bfields[i].repeat=repeat;
	fb->bfields[i].basesize=basesize;
	fb->bfields[i].offset=offset;
	fb->bfields[i].scale.bscale=1.0;
	fb->bfields[i].scale.bzero =0.0;
	fb->bfields[i].null[0]=0;
	fb->bfields[i].type[0]=0;
	fb->bfields[i].unit[0]=0;
	if ( basesize>0 )	offset+=basesize*repeat;
	else			offset+=(repeat+7)/8;
  }

 if ( offset<0 )	free(fb->bfields);
 va_end(ap);

 fb->rowsize=offset;
 fb->nrow=nrow;

 if ( offset<0 )	return(1);
 else			return(0);
}

int fits_bintable_check_fields(fitsbtable *fb,int nbfield,...)
{
 va_list	ap;
 int		i,repeat,form,basesize;

 if ( fb->nbfield != nbfield )	return(1);
 if ( fb->bfields == NULL )	return(1);

 va_start(ap,nbfield);

 for ( i=0 ; i<nbfield ; i++ )
  {	form=va_arg(ap,int);
	repeat=va_arg(ap,int);
	basesize=fits_bintable_form_basesize(form);
	if ( basesize<0 || repeat<=0 )		return(-1);
	if ( fb->bfields[i].form != form )	return(1);
	if ( fb->bfields[i].repeat != repeat )	return(1);
  }
 return(0);
}

int fits_bintable_set_xtr_params(fitsbtable *fb,int n,char *type,char *unit,char *null)
{
 fitsbfield	*bb;
 if ( n<0 || n>=fb->nbfield )	return(1);
 if ( fb->bfields==NULL )	return(1);
 bb=&fb->bfields[n];
 if ( type != NULL )	{ strncpy(bb->type,type,31);bb->type[31]=0; }
 if ( unit != NULL )	{ strncpy(bb->unit,unit,31);bb->unit[31]=0; }
 if ( null != NULL )	{ strncpy(bb->null,null,31);bb->null[31]=0; }
 return(0);
}

int fits_bintable_set_xtr_scale(fitsbtable *fb,int n,double bscale,double bzero)
{
 fitsbfield	*bb;
 if ( n<0 || n>=fb->nbfield )	return(1);
 if ( fb->bfields==NULL )	return(1);
 bb=&fb->bfields[n];
 bb->scale.bscale=bscale;
 bb->scale.bzero =bzero ;
 return(0);
}

int fits_bintable_alloc(fitsbtable *fb)
{
 int		naxis[2];
 naxis[0]=fb->rowsize;
 naxis[1]=fb->nrow;
 fb->allocdata=fits_tensor_alloc_arr(sizeof(unsigned char),2,naxis);
 fb->data=fb->allocdata;
 return(0);
}

int fits_bintable_swap_line(unsigned char *data,int rowsize,fitsbfield *bfs,int nbf)
{
 int		i,n,basesize;
 unsigned char	w;

 for ( i=0 ; i<nbf ; i++ )
  {	n=bfs[i].repeat;

	if ( bfs[i].form == FTF_COMPLEX32 || bfs[i].form == FTF_COMPLEX64 )
	 {	basesize=bfs[i].basesize/2;
		n=n*2;
	 }
	else
	 {	basesize=bfs[i].basesize;		}

	if ( basesize>1 )
	 {	while ( n>0 )
		 {	switch ( basesize )
			 {   case 2:
				w=data[0],data[0]=data[1],data[1]=w;
				break;
			     case 4:
				w=data[0],data[0]=data[3],data[3]=w;
				w=data[1],data[1]=data[2],data[2]=w;
				break;
			     case 8:
				w=data[0],data[0]=data[7],data[7]=w;
				w=data[1],data[1]=data[6],data[6]=w;
				w=data[2],data[2]=data[5],data[5]=w;
				w=data[3],data[3]=data[4],data[4]=w;
				break;
			 }
			data+=basesize;
			n--;
		 }
	 }
	else
	 {	data+=basesize*n;		}
  }

 return(0);
}

int fits_bintable_read_line_cb(int (*cb_read)(void *,void *,int),void *param,unsigned char *data,int rowsize,fitsbfield *bfs,int nbf)
{
 cb_read(param,data,rowsize);

 if ( bfs==NULL || nbf<0 )		return(0);

 if ( ! fits_arch_is_swapped() )
  {	return(0);		}
 else
  {	fits_bintable_swap_line(data,rowsize,bfs,nbf);
	return(0);
  }
}
int fits_bintable_read_cb(int (*cb_read)(void *,void *,int),void *param,fitsbtable *fb)
{
 int	cread,i;

 cread=0;
 for ( i=0 ; i<fb->nrow ; i++ )
  {	fits_bintable_read_line_cb(cb_read,param,fb->data[i],fb->rowsize,fb->bfields,fb->nbfield);
	cread=(cread+fb->rowsize)%FITS_TAPE_BLOCKSIZE;
  }
 if ( cread>0 )	cb_read(param,NULL,FITS_TAPE_BLOCKSIZE-cread);
 
 return(0);
}
int fits_bintable_read(FILE *fr,fitsbtable *fb)
{
 return(fits_bintable_read_cb(fits_cb_read,(void *)fr,fb));
}

int fits_bintable_skip_cb(int (*cb_read)(void *,void *,int),void *param,fitsbtable *fb)
{
 int	bseek,btot;

 btot=fb->rowsize*fb->nrow;
 if ( btot<=0 ) return(1);

 bseek=(btot+FITS_TAPE_BLOCKSIZE-1)/FITS_TAPE_BLOCKSIZE;
 cb_read(param,NULL,bseek*FITS_TAPE_BLOCKSIZE);

 return(0);
}
int fits_bintable_skip(FILE *fr,fitsbtable *fb)
{
 return(fits_bintable_skip_cb(fits_cb_read,(void *)fr,fb));
}

/*****************************************************************************/

int fits_bintable_set_offsets(fitsbtable *fb)
{
 int		i,offset,basesize;
 fitsbfield	*tb;

 offset=0;
 for ( i=0 ; i<fb->nbfield && fb->bfields != NULL ; i++ )
  {	tb=&fb->bfields[i];
	basesize=fits_bintable_form_basesize(tb->form);
	if ( basesize<0 )	return(-1);
	tb->basesize=basesize;
	tb->offset=offset;
	if ( basesize>0 )	offset+=basesize*tb->repeat;
	else			offset+=(tb->repeat+7)/8;
  }

 if ( offset != fb->rowsize )
  {	fb->rowsize=offset;
	return(1);
  }
 else
	return(0);
}

int fits_bintable_set_params(fitsheaderset *header,fitsbtable *fb)
{
 int		ret,i;
 char		chdr[16],cfrm[16];
 fitsbfield	*tb;

 if ( fb==NULL )	return(1);

 ret=0;

 fits_bintable_set_offsets(fb);

 ret|=fits_headerset_set_string(header,headers[HDR_XTENSION],FITS_SH_FIRST,"BINTABLE",NULL);
 ret|=fits_headerset_set_integer(header,headers[HDR_BITPIX],FITS_SH_FIRST,8,NULL);
 ret|=fits_headerset_set_integer(header,headers[HDR_NAXIS],FITS_SH_FIRST,2,NULL);
 ret|=fits_headerset_set_integer(header,headers[HDR_NAXIS1],FITS_SH_FIRST,fb->rowsize,NULL);
 ret|=fits_headerset_set_integer(header,headers[HDR_NAXIS2],FITS_SH_FIRST,fb->nrow,NULL);
 ret|=fits_headerset_set_integer(header,headers[HDR_PCOUNT],FITS_SH_FIRST,0,NULL);
 ret|=fits_headerset_set_integer(header,headers[HDR_GCOUNT],FITS_SH_FIRST,1,NULL);
 ret|=fits_headerset_set_integer(header,headers[HDR_TFIELDS],FITS_SH_FIRST,fb->nbfield,NULL);
 for ( i=0 ; i<fb->nbfield && fb->bfields != NULL ; i++ )
  {	tb=&fb->bfields[i];

	sprintf(chdr,"%s%d",headers[HDR_TFORM],i+1);
	sprintf(cfrm,"%d%c",tb->repeat,tb->form);
	ret|=fits_headerset_set_string(header,chdr,FITS_SH_FIRST,cfrm,NULL);

	if ( tb->null[0] )
	 {	sprintf(chdr,"%s%d",headers[HDR_TNULL],i+1);
		ret|=fits_headerset_set_string(header,chdr,FITS_SH_FIRST,tb->null,NULL);
	 }
	if ( tb->type[0] )
	 {	sprintf(chdr,"%s%d",headers[HDR_TTYPE],i+1);
		ret|=fits_headerset_set_string(header,chdr,FITS_SH_FIRST,tb->type,NULL);
	 }
	if ( tb->unit[0] )
	 {	sprintf(chdr,"%s%d",headers[HDR_TUNIT],i+1);
		ret|=fits_headerset_set_string(header,chdr,FITS_SH_FIRST,tb->unit,NULL);
	 }
	if ( tb->scale.bscale != 1.0 || tb->scale.bzero != 0.0 )
	 {	sprintf(chdr,"%s%d",headers[HDR_TSCAL],i+1);
		ret|=fits_headerset_set_double(header,chdr,FITS_SH_FIRST,tb->scale.bscale,NULL);
		sprintf(chdr,"%s%d",headers[HDR_TZERO],i+1);
		ret|=fits_headerset_set_double(header,chdr,FITS_SH_FIRST,tb->scale.bzero ,NULL);
	 }
  }
 
 return(ret);
}
 
/*****************************************************************************/

int fits_bintable_write_cb(int (*cb_write)(void *,void *,int),void *param,fitsbtable *fb,int is_pad)
{
 int		wr;
 int		i,psize;
 unsigned char *data;
 
 if ( fb==NULL || fb->data==NULL )	return(0);

 wr=0;
 
 data=(unsigned char *)malloc(fb->rowsize);
 for ( i=0 ; i<fb->nrow ; i++ )
  {	memcpy(data,fb->data[i],fb->rowsize);
	fits_bintable_swap_line(data,fb->rowsize,fb->bfields,fb->nbfield);
 	cb_write(param,data,fb->rowsize);
	wr+=fb->rowsize;
  }
 free(data);

 if ( is_pad && wr%FITS_TAPE_BLOCKSIZE>0 )
  {	psize=FITS_TAPE_BLOCKSIZE-(wr%FITS_TAPE_BLOCKSIZE);
	data=(unsigned char *)malloc(psize);
	memset(data,0,psize);
	cb_write(param,data,psize);
	wr+=psize;
	free(data);
  }

 return(wr);
}
int fits_bintable_write(FILE *fw,fitsbtable *fb,int is_pad)
{
 return(fits_bintable_write_cb(fits_cb_write,(void *)fw,fb,is_pad));
}

/*****************************************************************************/

int fits_bintable_free(fitsbtable *fb)
{
 if ( fb->allocdata != NULL )   fits_tensor_free(fb->allocdata);
 fb->allocdata=NULL;
 fb->data=NULL;
 if ( fb->bfields != NULL && fb->nbfield>0 )	free(fb->bfields);
 fb->nbfield=0;
 fb->bfields=NULL;
 fb->nrow=0;
 fb->rowsize=0;
 fb->pcount=0;
 return(0);
}

/*****************************************************************************/

int fits_bintable_duplicate(fitsbtable *db,fitsbtable *sb,int flag)
{
 int	naxis[2]; 

 db->nrow=sb->nrow;
 db->rowsize=sb->rowsize;

 db->nbfield=sb->nbfield;
 db->bfields=(fitsbfield *)malloc(sizeof(fitsbfield)*db->nbfield);
 memcpy(db->bfields,sb->bfields,sizeof(fitsbfield)*db->nbfield);

 naxis[0]=db->rowsize;
 naxis[1]=db->nrow;
 if ( db->data != NULL )
  {	db->data=fits_tensor_alloc_arr(1,2,naxis);
	memcpy(&db->data[0][0],&sb->data[0][0],db->rowsize*db->nrow);
	db->allocdata=db->data;
  }
 else
  {	db->data=NULL;
	db->allocdata=NULL;
  }

 return(0);
}

/******************************************************************************
 Reference:
  [1]:	Definition of the Flexible Image Transport System (FITS)
	NASA / Science Office of Standards and Technology (NOST)
	NOST 100-2.0 (March 29, 1999)
	NASA Goddard Space Flight Center
	http://fits.gsfc.nasa.gov/fits_documentation.html
  [2]:  Author's comment: this library contains the implementation of the 
	unreferenced int64_t field 'K' as well.
******************************************************************************/
