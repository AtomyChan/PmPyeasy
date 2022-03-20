/*****************************************************************************/
/* fits-table.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Another simple standalone library for manipulating FITS files:	     */
/* ASCII table extensions (XTENSION = 'TABLE')				     */
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

#define	_FITS_SOURCE

#include <fits/fits.h>
#include "fits-common.h"

/*****************************************************************************/

int fits_table_get_params(fitsheaderset *header,fitsttable *ft)
{
 fitsheader	*hdr;
 fitstfield	*tt;
 int		i;
 char		chdr[16];

 memset(ft,0,sizeof(fitsttable));

 hdr=fits_headerset_get_uniq_header(header,headers[HDR_NAXIS]);
 if ( hdr==NULL || hdr->vtype != FITS_VINT || hdr->vint != 2 )	return(1);
 hdr=fits_headerset_get_uniq_header(header,headers[HDR_BITPIX]);
 if ( hdr==NULL || hdr->vtype != FITS_VINT || hdr->vint != 8 )	return(1);
 hdr=fits_headerset_get_uniq_header(header,headers[HDR_NAXIS1]);
 if ( hdr==NULL || hdr->vtype != FITS_VINT || hdr->vint<=0 )	return(1);
 ft->rowsize=hdr->vint;
 hdr=fits_headerset_get_uniq_header(header,headers[HDR_NAXIS2]);
 if ( hdr==NULL || hdr->vtype != FITS_VINT || hdr->vint<=0 )	return(1);
 ft->nrow=hdr->vint;
 hdr=fits_headerset_get_uniq_header(header,headers[HDR_TFIELDS]);
 if ( hdr==NULL || hdr->vtype != FITS_VINT || hdr->vint<=0 )	return(1);
 ft->ntfield=hdr->vint;

 ft->tfields=(fitstfield *)malloc(sizeof(fitstfield)*ft->ntfield);
 
 for ( i=0 ; i<ft->ntfield ; i++ )
  {	tt=&ft->tfields[i];

	sprintf(chdr,"%s%d",headers[HDR_TBCOL],i+1);
	hdr=fits_headerset_get_uniq_header(header,chdr);
	if ( hdr==NULL || hdr->vtype != FITS_VINT || hdr->vint<=0 )
	 {	free(ft->tfields);return(1);				}
	tt->colindex=hdr->vint-1;
	sprintf(chdr,"%s%d",headers[HDR_TFORM],i+1);
	hdr=fits_headerset_get_uniq_header(header,chdr);
	if ( hdr==NULL || hdr->vtype != FITS_VSTR )
	 {	free(ft->tfields);return(1);				}
	strncpy(tt->format,hdr->vstr,11);tt->format[11]=0;

	tt->scale.bscale=1.0;
	tt->scale.bzero =0.0;
	tt->null[0]=0;
	tt->type[0]=0;
	tt->unit[0]=0;

	sprintf(chdr,"%s%d",headers[HDR_TSCAL],i+1);
	hdr=fits_headerset_get_uniq_header(header,chdr);
	if ( hdr != NULL )
	 {	if ( hdr->vtype==FITS_VINT )	tt->scale.bscale=hdr->vint;
		if ( hdr->vtype==FITS_VDOUBLE )	tt->scale.bscale=hdr->vdouble;
	 }
	sprintf(chdr,"%s%d",headers[HDR_TZERO],i+1);
	hdr=fits_headerset_get_uniq_header(header,chdr);
	if ( hdr != NULL )
	 {	if ( hdr->vtype==FITS_VINT )	tt->scale.bzero=hdr->vint;
		if ( hdr->vtype==FITS_VDOUBLE )	tt->scale.bzero=hdr->vdouble;
	 }
	sprintf(chdr,"%s%d",headers[HDR_TNULL],i+1);
	hdr=fits_headerset_get_uniq_header(header,chdr);
	if ( hdr != NULL && hdr->vtype==FITS_VSTR )
	 {	strncpy(tt->null,hdr->vstr,31);tt->null[31]=0;		}
	sprintf(chdr,"%s%d",headers[HDR_TTYPE],i+1);
	hdr=fits_headerset_get_uniq_header(header,chdr);
	if ( hdr != NULL && hdr->vtype==FITS_VSTR )
	 {	strncpy(tt->type,hdr->vstr,31);tt->type[31]=0;		}
	sprintf(chdr,"%s%d",headers[HDR_TUNIT],i+1);
	hdr=fits_headerset_get_uniq_header(header,chdr);
	if ( hdr != NULL && hdr->vtype==FITS_VSTR )
	 {	strncpy(tt->unit,hdr->vstr,31);tt->unit[31]=0;		}
	
  }

 return(0);
}

int fits_table_alloc(fitsttable *ft)
{
 int	naxis[2];
 naxis[0]=ft->rowsize;
 naxis[1]=ft->nrow;
 ft->allocdata=fits_tensor_alloc_arr(sizeof(unsigned char),2,naxis);
 ft->data=ft->allocdata;
 return(0);
}

int fits_table_read_cb(int (*cb_read)(void *,void *,int),void *param,fitsttable *ft)
{
 int	cread,i;

 cread=0;
 for ( i=0 ; i<ft->nrow ; i++ )
  {	cb_read(param,ft->data[i],ft->rowsize);
	cread=(cread+ft->rowsize)%FITS_TAPE_BLOCKSIZE;
  }
 if ( cread>0 )	cb_read(param,NULL,FITS_TAPE_BLOCKSIZE-cread);
 
 return(0);
}
int fits_table_read(FILE *fr,fitsttable *ft)
{
 return(fits_table_read_cb(fits_cb_read,(void *)fr,ft));
}

int fits_table_skip_cb(int (*cb_read)(void *,void *,int),void *param,fitsttable *ft)
{
 int	bseek,btot;

 btot=ft->rowsize*ft->nrow;
 if ( btot<=0 )	return(1);

 bseek=(btot+FITS_TAPE_BLOCKSIZE-1)/FITS_TAPE_BLOCKSIZE;
 cb_read(param,NULL,bseek*FITS_TAPE_BLOCKSIZE);

 return(0);
}
int fits_table_skip(FILE *fr,fitsttable *ft)
{
 return(fits_table_skip_cb(fits_cb_read,(void *)fr,ft));
}

/*****************************************************************************/

int fits_table_set_params(fitsheaderset *header,fitsttable *ft)
{
 int		ret,i;
 char		chdr[16];
 fitstfield	*tt;

 if ( ft==NULL )	return(1);

 ret=0;

 ret|=fits_headerset_set_string(header,headers[HDR_XTENSION],FITS_SH_FIRST,"TABLE",NULL);
 ret|=fits_headerset_set_integer(header,headers[HDR_BITPIX],FITS_SH_FIRST,8,NULL);
 ret|=fits_headerset_set_integer(header,headers[HDR_NAXIS],FITS_SH_FIRST,2,NULL);
 ret|=fits_headerset_set_integer(header,headers[HDR_NAXIS1],FITS_SH_FIRST,ft->rowsize,NULL);
 ret|=fits_headerset_set_integer(header,headers[HDR_NAXIS2],FITS_SH_FIRST,ft->nrow,NULL);
 ret|=fits_headerset_set_integer(header,headers[HDR_PCOUNT],FITS_SH_FIRST,0,NULL);
 ret|=fits_headerset_set_integer(header,headers[HDR_GCOUNT],FITS_SH_FIRST,1,NULL);
 ret|=fits_headerset_set_integer(header,headers[HDR_TFIELDS],FITS_SH_FIRST,ft->ntfield,NULL);
 for ( i=0 ; i<ft->ntfield && ft->tfields != NULL ; i++ )
  {	tt=&ft->tfields[i];
	sprintf(chdr,"%s%d",headers[HDR_TBCOL],i+1);
	ret|=fits_headerset_set_integer(header,chdr,FITS_SH_FIRST,tt->colindex+1,NULL);
	sprintf(chdr,"%s%d",headers[HDR_TFORM],i+1);
	ret|=fits_headerset_set_string(header,chdr,FITS_SH_FIRST,tt->format,NULL);
	if ( tt->null[0] )
	 {	sprintf(chdr,"%s%d",headers[HDR_TNULL],i+1);
		ret|=fits_headerset_set_string(header,chdr,FITS_SH_FIRST,tt->null,NULL);
	 }
	if ( tt->type[0] )
	 {	sprintf(chdr,"%s%d",headers[HDR_TTYPE],i+1);
		ret|=fits_headerset_set_string(header,chdr,FITS_SH_FIRST,tt->type,NULL);
	 }
	if ( tt->unit[0] )
	 {	sprintf(chdr,"%s%d",headers[HDR_TUNIT],i+1);
		ret|=fits_headerset_set_string(header,chdr,FITS_SH_FIRST,tt->unit,NULL);
	 }
	if ( tt->scale.bscale != 1.0 || tt->scale.bzero != 0.0 )
	 {	sprintf(chdr,"%s%d",headers[HDR_TSCAL],i+1);
		ret|=fits_headerset_set_double(header,chdr,FITS_SH_FIRST,tt->scale.bscale,NULL);
		sprintf(chdr,"%s%d",headers[HDR_TZERO],i+1);
		ret|=fits_headerset_set_double(header,chdr,FITS_SH_FIRST,tt->scale.bzero ,NULL);
	 }
  }
 
 return(ret);
}

/*****************************************************************************/

int fits_table_write_cb(int (*cb_write)(void *,void *,int),void *param,fitsttable *ft,int is_pad)
{
 int	wr;
 int	i,psize;
 char	*padding;

 if ( ft==NULL || ft->data==NULL )	return(0);

 wr=0;
 for ( i=0 ; i<ft->nrow ; i++ )
  {	cb_write(param,ft->data[i],ft->rowsize);
	wr+=ft->rowsize;
  }
	
 if ( is_pad && wr%FITS_TAPE_BLOCKSIZE>0 )
  {	psize=FITS_TAPE_BLOCKSIZE-(wr%FITS_TAPE_BLOCKSIZE);
	padding=(char *)malloc(psize);
	memset(padding,32,psize);
	cb_write(param,padding,psize);
	wr+=psize;
	free(padding);
  }

 return(wr);
}
int fits_table_write(FILE *fw,fitsttable *ft,int is_pad)
{
 return(fits_table_write_cb(fits_cb_write,(void *)fw,ft,is_pad));
}

/*****************************************************************************/

int fits_table_free(fitsttable *ft)
{
 if ( ft->allocdata != NULL )	fits_tensor_free(ft->allocdata);
 ft->allocdata=NULL;
 ft->data=NULL;
 if ( ft->tfields != NULL && ft->ntfield>0 )	free(ft->tfields);
 ft->ntfield=0;
 ft->tfields=NULL;
 ft->nrow=0;
 ft->rowsize=0;
 return(0);
}

int fits_table_duplicate(fitsttable *dt,fitsttable *st,int flag)
{
 int	naxis[2];

 dt->nrow=st->nrow;
 dt->rowsize=st->rowsize;

 dt->ntfield=st->ntfield;
 dt->tfields=(fitstfield *)malloc(sizeof(fitstfield)*dt->ntfield);
 memcpy(dt->tfields,st->tfields,sizeof(fitstfield)*dt->ntfield);

 naxis[0]=dt->rowsize;
 naxis[1]=dt->nrow;
 if ( dt->data != NULL )
  {	dt->data=fits_tensor_alloc_arr(1,2,naxis);
	memcpy(&dt->data[0][0],&st->data[0][0],dt->rowsize*dt->nrow);
	dt->allocdata=dt->data;
  }
 else
  {	dt->data=NULL;
	dt->allocdata=NULL;
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
******************************************************************************/
                           
