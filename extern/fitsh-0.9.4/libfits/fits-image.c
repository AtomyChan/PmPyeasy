/*****************************************************************************/
/* fits-image.c 							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Another simple standalone library for manipulating FITS files:	     */
/* Image extensions (XTENSION = 'IMAGE') and primary images.		     */
/* (c) 2004-06, Pal, A. (apal@szofi.elte.hu).				     */
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

int fits_image_alloc(fitsimage *fi,int sx,int sy)
{
 double **r;
 int 	arr[2];

 arr[0]=sx,arr[1]=sy;
 r=(double **)fits_tensor_alloc_arr(sizeof(double),2,arr);

 fi->dim=2;
 fi->naxis[0]=sx;
 fi->naxis[1]=sy;

 fi->sx=sx,
 fi->sy=sy;

 fi->data=r;
 fi->allocdata=fi->vdata=(void *)r;

 return(0);
}

double ** fits_image_first_layer(void *data,int dim)
{
 while ( dim>2 )
  {	data=*(void **)data;
	dim--;
  };
 return((double **)data);
}
double * fits_image_first_pixel(void *data,int dim)
{
 while ( dim>1 )
  {	data=*(void **)data;
	dim--;
  };
 return((double *)data);
}
int fits_image_total_pixels(int dim,int *naxis)
{
 int	ntot;
 ntot=1;
 while ( dim>0 )
  {	if ( *naxis <= 0 )	return(-1);
	ntot=ntot*(*naxis);
	dim--,naxis++;
  };
 return(ntot);
}

int fits_image_alloc_gen(fitsimage *fi,int dim,int *pnaxis)
{
 void	*vdata;
 int	i,*naxis,cnaxis[2];

 if ( dim<1 || dim>FITS_MAX_NAXIS )	return(1);
 else if ( dim==1 )
  {	cnaxis[0]=pnaxis[0];
	cnaxis[1]=1;
	dim=2;
	naxis=cnaxis;
  }
 else	naxis=pnaxis;

 vdata=fits_tensor_alloc_arr(sizeof(double),dim,naxis);
 if ( vdata==NULL )			return(-1);

 fi->dim=dim;
 for ( i=0 ; i<dim ; i++ )		fi->naxis[i]=naxis[i];
 fi->data=fits_image_first_layer(vdata,dim);
 fi->allocdata=vdata;
 fi->vdata=vdata;
 fi->sx=naxis[0];
 fi->sy=naxis[1];
 
 return(0);
}

int fits_image_alloc_2d(fitsimage *fi,int sx,int sy)
{
 return ( fits_image_alloc(fi,sx,sy) );
}
int fits_image_alloc_3d(fitsimage *fi,int sx,int sy,int sz)
{
 int	pnaxis[3];
 pnaxis[0]=sx,
 pnaxis[1]=sy,
 pnaxis[2]=sz;
 return ( fits_image_alloc_gen(fi,3,pnaxis) );
}

/*****************************************************************************/

void fits_image_free(fitsimage *fi)
{
 if ( fi->allocdata != NULL )	fits_tensor_free(fi->allocdata);
 fi->allocdata=fi->vdata=NULL;
 fi->data=NULL;
 fi->sx=0,fi->sy=0;
 fi->dim=0;
}

/*****************************************************************************/

int fits_image_set_value(fitsimage *fi,double value)
{
 int	i,ntot;
 double	*ddata;
 if ( fi==NULL || fi->vdata==NULL )	return(1);
 ddata=fits_image_first_pixel (fi->vdata,fi->dim);
 ntot =fits_image_total_pixels(fi->dim,fi->naxis);
 for ( i=0 ; i<ntot ; i++ )
  {	ddata[i]=value;		}
 return(0);
}
int fits_image_reset(fitsimage *fi)
{
 int	r;
 r=fits_image_set_value(fi,0.0);
 return(r);
}

/*****************************************************************************/

int fits_image_read_line_cb(int (*cb_read)(void *,void *,int),void *param,int sx,int bit,double *line)
{
 int		j,bs;
 unsigned char	*wbuff,*w;

 bs=bit;
 if ( bs<0 )	bs=-bs;
 bs=bs/8;

 wbuff=(unsigned char *)malloc((unsigned)(sx*bs));
 if ( wbuff==NULL )	return(0);

 cb_read(param,wbuff,bs*sx);

 fits_swap_line_bytes(wbuff,bs,sx);

 switch ( bit )
  {  case 8:
	for ( j=0,w=wbuff ; j<sx ; j++,w++ )
	 {	line[j]=(double)(*(unsigned char *)(w));	}
	break;
     case 16:
	for ( j=0,w=wbuff ; j<sx ; j++,w+=2 )
	 {	line[j]=(double)(*(short *)(w));		}
	break;
     case 32:
	if ( sizeof(long)==4 )		/* i386 */
	 {	for ( j=0,w=wbuff ; j<sx ; j++,w+=4 )
		 {	line[j]=(double)(*(long *)(w));		}
	 }
	else if ( sizeof(int)==4 )	/* ia64 */
	 {	for ( j=0,w=wbuff ; j<sx ; j++,w+=4 )
		 {	line[j]=(double)(*(int  *)(w));		}
	 }
	else				/* ???? */
	 {	for ( j=0,w=wbuff ; j<sx ; j++,w+=4 )
		 {	line[j]=0.0;				}
	 }
	break;
     case -32:
	for ( j=0,w=wbuff ; j<sx ; j++,w+=4 )
	 {	line[j]=(double)(*(float *)(w));	}
	break;
     case -64:
	for ( j=0,w=wbuff ; j<sx ; j++,w+=8 )
	 {	line[j]=(double)(*(double *)(w));	}
	break;
     default:
	for ( j=0,w=wbuff ; j<sx ; j++,w+=4 )
	 {	line[j]=0.0;				}
	break;
  }
 free(wbuff);

 return(bs*sx);
}
int fits_image_read_line(FILE *fr,int sx,int bit,double *line)
{
 return(fits_image_read_line_cb(fits_cb_read,(void *)fr,sx,bit,line));
}

/*****************************************************************************/

int fits_image_get_scale(fitsheaderset *header,fitsimage *fi,double *rbscale,double *rbzero)
{
 fitsheader	*hdr;
 double		bscale,bzero;

 *rbscale=1.0;
 *rbzero =0.0;

 if ( fi==NULL )	return(1);

 hdr=fits_headerset_get_header(header,headers[HDR_BSCALE],0);
 if ( hdr==NULL )			bscale=1.0;			/* There's no BSCALE header, default value is 1.0 */
 else if ( hdr->vtype==FITS_VINT )	bscale=(double)hdr->vint;	/* Integer */
 else if ( hdr->vtype==FITS_VDOUBLE )	bscale=(double)hdr->vdouble;	/* Double */
 else					return(1);			/* Everything else, ignored */
 hdr=fits_headerset_get_header(header,headers[HDR_BZERO],0);
 if ( hdr==NULL )			bzero=0.0;			/* Default BZERO is 0 */
 else if ( hdr->vtype==FITS_VINT )	bzero=(double)hdr->vint;
 else if ( hdr->vtype==FITS_VDOUBLE )	bzero=(double)hdr->vdouble;
 else					return(1);

 if ( rbscale != NULL )	*rbscale=bscale;
 if ( rbzero  != NULL )	*rbzero =bzero;

 return(0);
}

int fits_image_get_params(fitsheaderset *header,fitsimage *fi)
{
 int		sx,sy,bit,naxis[FITS_MAX_NAXIS],dim,i;
 fitsheader	*hdr;
 char		naxishdr[16];

 hdr=fits_headerset_get_uniq_header(header,headers[HDR_NAXIS]);
 if ( hdr==NULL || hdr->vtype != FITS_VINT )	return(1);
 dim=hdr->vint;
 if ( dim<1 || dim>FITS_MAX_NAXIS )	return(1);

 for ( i=0 ; i<dim ; i++ )
  {	sprintf(naxishdr,"%s%d",headers[HDR_NAXIS],i+1);
	hdr=fits_headerset_get_uniq_header(header,naxishdr);
	if ( hdr==NULL || hdr->vtype != FITS_VINT )	return(1);
	naxis[i]=hdr->vint;
  }
 sx=naxis[0];
 if ( dim==1 )	sy=1;
 else		sy=naxis[1];

 hdr=fits_headerset_get_uniq_header(header,headers[HDR_BITPIX]);
 if ( hdr==NULL || hdr->vtype != FITS_VINT )	return(1);
 bit=hdr->vint;

 fi->dim=dim;
 for ( i=0 ; i<dim ; i++ )
  {	fi->naxis[i]=naxis[i];		}

 fi->sx=sx,
 fi->sy=sy,
 fi->bit=bit;

 fits_image_get_scale(header,fi,&fi->read.bscale,&fi->read.bzero);
 fi->curr.bscale=fi->read.bscale;
 fi->curr.bzero =fi->read.bzero;

 return(0);
}

/*****************************************************************************/

int fits_image_rescale(fitsimage *fi)
{
 int		i,ntot;
 double		bscale,bzero,*ddata;

 if ( fi==NULL  || fi->vdata==NULL )	return(0);
 if ( fi->sx==0 || fi->sy==0 )		return(0);

 bscale=fi->curr.bscale;
 bzero =fi->curr.bzero;

 if ( bscale==1.0 && bzero==0.0 )	return(0);	/* Identity... */

 ddata=fits_image_first_pixel(fi->vdata,fi->dim);
 ntot=fits_image_total_pixels(fi->dim,fi->naxis);
 for ( i=0 ; i<ntot ; i++ )
  {	ddata[i]=bzero+bscale*ddata[i];		}

 fi->curr.bscale=1.0;
 fi->curr.bzero =0.0;

 return(0);
}

int fits_image_backscale(fitsimage *fi,double ibscale,double ibzero)
{
 int		i,ntot;
 double		bscale,bzero,*ddata;

 if ( fi==NULL  || fi->vdata==NULL )	return(0);
 if ( fi->sx==0 || fi->sy==0 )		return(0);
 if ( ibscale==0.0 )			return(1);	/* invalid scale */

 bscale=fi->curr.bscale;
 bzero =fi->curr.bzero;

 if ( bscale==1.0 && bzero==0.0 && ibscale==1.0 && ibzero==0.0 ) return(0);

 ddata=fits_image_first_pixel(fi->vdata,fi->dim);
 ntot=fits_image_total_pixels(fi->dim,fi->naxis);

 if ( bscale != 1.0 || bzero != 0.0 )
  {	for ( i=0 ; i<ntot ; i++ )
	 {	ddata[i]=(bzero+bscale*ddata[i]-ibzero)/ibscale;	}
  }
 else
  {	for ( i=0 ; i<ntot ; i++ )
	 {	ddata[i]=(ddata[i]-ibzero)/ibscale;			}
  }

 fi->curr.bscale=ibscale;
 fi->curr.bzero =ibzero;

 return(0);
}

/*****************************************************************************/

int fits_image_quantize(fitsimage *fi,int nquantizebit)
{
 int		i,ntot;
 double		*ddata,x,bm,mm;
 long		l,m;

 if ( fi==NULL  || fi->vdata==NULL )	return(0);
 if ( fi->sx==0 || fi->sy==0 )		return(0);

 if ( nquantizebit>=sizeof(long double)*8 )	return(0);

 ddata=fits_image_first_pixel(fi->vdata,fi->dim);
 ntot=fits_image_total_pixels(fi->dim,fi->naxis);

 if ( nquantizebit==0 )
  {	for ( i=0 ; i<ntot ; i++ )
	 {	ddata[i]=floor(ddata[i]);		}
  }
 else if ( nquantizebit<0 )
  {	m=~((1<<(-nquantizebit))-1);
	for ( i=0 ; i<ntot ; i++ )
	 {	x=floor(ddata[i]);
		l=(long)x;
		ddata[i]=(double)(l&m);
	 }
  }
 else
  {	m=(1<<nquantizebit);
	bm=1.0/(double)(m);
	mm=(double)m;
	for ( i=0 ; i<ntot ; i++ )
	 {	ddata[i]=bm*floor(ddata[i]*mm);		}
  }

 return(0);
}

/*****************************************************************************/

int fits_image_set_params(fitsheaderset *header,fitsimage *fi)
{
 int		k=0,i,nn,nc;
 char		naxishdr[16];

 if ( fi==NULL )	return(1);

 k|=fits_headerset_set_integer(header,headers[HDR_BITPIX],FITS_SH_FIRST,fi->bit,NULL);

 k|=fits_headerset_set_integer(header,headers[HDR_NAXIS],FITS_SH_FIRST,fi->dim,NULL);
 for ( i=0 ; i<fi->dim ; i++ )
  {	sprintf(naxishdr,"%s%d",headers[HDR_NAXIS],i+1);
	k|=fits_headerset_set_integer(header,naxishdr,FITS_SH_FIRST,fi->naxis[i],NULL);
  }
 for ( i=fi->dim ; i<=999 ; i++ )
  {	sprintf(naxishdr,"%s%d",headers[HDR_NAXIS],i+1);
	if ( !  fits_headerset_get_count(header,naxishdr) )	break;
	else	fits_headerset_delete_all(header,naxishdr);
  }

 nn=fits_headerset_get_id(header,headers[HDR_NAXIS],0);
 if ( nn>=0 )	nn++;
 for ( i=0 ; i<fi->dim && nn>=0 ; i++ )
  {	sprintf(naxishdr,"%s%d",headers[HDR_NAXIS],i+1);
	nc=fits_headerset_get_id(header,naxishdr,0);
	/* fprintf(stderr,"nn=%d nc=%d\n",nn,nc); */
	if ( nc>nn )
	 {	fitsheader	th;
		th=header->hdrs[nc];
		memmove(header->hdrs+nn+1,header->hdrs+nn,sizeof(fitsheader)*(nc-nn));
		header->hdrs[nn]=th;
		nn++;
	 }
	else if ( nc==nn )
		nn++;
  }

 k|=fits_headerset_set_double (header,headers[HDR_BSCALE],FITS_SH_FIRST,fi->read.bscale,NULL);
 k|=fits_headerset_set_double (header,headers[HDR_BZERO],FITS_SH_FIRST ,fi->read.bzero,NULL);

 return(k);
}

/*****************************************************************************/

int fits_image_read_cb(int (*cb_read)(void *,void *,int),void *param,fitsimage *fi)
{
 int	i,cread,ntot,nline,sx;
 double	*ddata;

 ddata=fits_image_first_pixel(fi->vdata,fi->dim);
 ntot=fits_image_total_pixels(fi->dim,fi->naxis);
 sx=fi->naxis[0];
 nline=ntot/sx;

 for ( i=0,cread=0 ; i<nline ; i++ )
	cread=(cread+fits_image_read_line_cb(cb_read,param,sx,fi->bit,ddata+i*sx))%2880;

 if ( cread>0 )	cb_read(param,NULL,2880-cread);
 
 return(0);
}
int fits_image_read(FILE *fr,fitsimage *fi)
{
 return(fits_image_read_cb(fits_cb_read,(void *)fr,fi));
}

int fits_image_skip_cb(int (*cb_read)(void *,void *,int),void *param,fitsimage *fi)
{
 int	ntot,bs,btot,bseek;

 bs=fi->bit;
 if ( bs<0 )	bs=-bs;
 bs=bs/8;
 if ( bs<=0 )	return(1);

 ntot=fits_image_total_pixels(fi->dim,fi->naxis);
 if ( ntot<=0 )	return(1);
 btot=ntot*bs;
 bseek=(btot+FITS_TAPE_BLOCKSIZE-1)/FITS_TAPE_BLOCKSIZE;
 cb_read(param,NULL,bseek*FITS_TAPE_BLOCKSIZE);
 
 return(0);
}
int fits_image_skip(FILE *fr,fitsimage *fi)
{
 return(fits_image_skip_cb(fits_cb_read,(void *)fr,fi));
}

/*****************************************************************************/

int fits_image_write_line_cb(int (*cb_write)(void *,void *,int),void *param,int sx,int bit,double *line)
{
 int		j,bs;
 unsigned char	*wbuff,*w;

 bs=bit;
 if ( bs<0 )	bs=-bs;
 bs=bs/8;

 wbuff=(unsigned char *)malloc((unsigned)(sx*bs));
 if ( wbuff==NULL )	return(0);	/* Allocation error, nothing written... */

 switch ( bit )
  {  case 8:
	for ( j=0,w=wbuff ; j<sx ; j++,w++ )
	 {	*(unsigned char *)w=(unsigned char)line[j];	}
	break;
     case 16:
	for ( j=0,w=wbuff ; j<sx ; j++,w+=2 )
	 {	*(short *)w=(short)line[j];			}
	break;
     case 32:
	if ( sizeof(long)==4 )		/* i386 */
	 {	for ( j=0,w=wbuff ; j<sx ; j++,w+=4 )
		 {	*(long *)w=(long)line[j];		}
	 }
	else if ( sizeof(int)==4 )	/* ia64 */
	 {	for ( j=0,w=wbuff ; j<sx ; j++,w+=4 )
		 {	*(int *)w =(int)line[j];		}
	 }
	else				/* ???? */
	 {	for ( j=0,w=wbuff ; j<sx ; j++,w+=4 )
		 {	*w=*(w+1)=*(w+2)=*(w+3)=0;		}
	 }
	break;
     case -32:
	for ( j=0,w=wbuff ; j<sx ; j++,w+=4 )
	 {	*(float *)w=(float)line[j];		}
	break;
     case -64:
	for ( j=0,w=wbuff ; j<sx ; j++,w+=8 )
	 {	*(double *)w=(double)line[j];		}
	break;
  }

 fits_swap_line_bytes(wbuff,bs,sx);

 cb_write(param,wbuff,bs*sx);

 free(wbuff);

 return(bs*sx);		/* total bytes written */
}
int fits_image_write_line(FILE *fw,int sx,int bit,double *line)
{
 return(fits_image_write_line_cb(fits_cb_write,(void *)fw,sx,bit,line));
}

int fits_image_write_cb(int (*cb_write)(void *,void *,int),void *param,fitsimage *fi,int is_pad)
{
 int	i,wr,ntot,nline,sx,psize;
 double	*ddata;
 char	*padding;

 ddata=fits_image_first_pixel(fi->vdata,fi->dim);
 ntot=fits_image_total_pixels(fi->dim,fi->naxis);
 sx=fi->naxis[0];
 nline=ntot/sx;

 wr=0;
 for ( i=0 ; i<nline ; i++ )
  {	wr+=fits_image_write_line_cb(cb_write,param,sx,fi->bit,ddata+i*sx);		}

 if ( is_pad && wr%FITS_TAPE_BLOCKSIZE )
  {	psize=FITS_TAPE_BLOCKSIZE-(wr%FITS_TAPE_BLOCKSIZE);
	padding=(char *)malloc(psize);
	memset(padding,0,psize);
	cb_write(param,padding,psize);
	wr+=psize;
	free(padding);
  }

 return(wr);		/* total bytes written */
}
int fits_image_write(FILE *fw,fitsimage *fi,int is_pad)
{
 return(fits_image_write_cb(fits_cb_write,(void *)fw,fi,is_pad));
}

/*****************************************************************************/

double * fits_image_dump(fitsimage *fi)
{
 int	ntot;
 double	*ddata,*raw;

 if ( fi==NULL || fi->vdata==NULL )	return(NULL);

 ddata=fits_image_first_pixel(fi->vdata,fi->dim);
 ntot=fits_image_total_pixels(fi->dim,fi->naxis);

 raw=(double *)malloc(sizeof(double)*ntot);
 if ( raw==NULL )	return(NULL);
 memcpy(raw,ddata,sizeof(double)*ntot);

 return(raw);
}

/*****************************************************************************/

int fits_image_duplicate(fitsimage *ret,fitsimage *img,int flag)
{
 double	*idata,*rdata;
 int	ntot;

 if ( img->vdata != NULL )
  {	ret->bit=img->bit;
	memcpy(&ret->read,&img->read,sizeof(fitsscale));
	memcpy(&ret->curr,&img->curr,sizeof(fitsscale));

	/* this also sets ret->sx, ret->sy and ret->data properly */
	/*
	fprintf(stderr,"%d [%dx%d]\n",img->dim,img->naxis[0],img->naxis[1]);
	*/
	fits_image_alloc_gen(ret,img->dim,img->naxis);
	ntot=fits_image_total_pixels(img->dim,img->naxis);
	idata=fits_image_first_pixel(img->vdata,img->dim);
	rdata=fits_image_first_pixel(ret->vdata,ret->dim);
	/*
	fprintf(stderr,"rdata=%p\n",(void *)rdata);
	fprintf(stderr,"ret->vdata=%p\n",(void *)ret->vdata);
	*/
        if ( flag )     memcpy(rdata,idata,sizeof(double)*ntot);
        else            memset(rdata,0,sizeof(double)*ntot);
  }
 else
  {	ret->sx=ret->sy=0;
	ret->bit=0;
	ret->data=NULL;
	ret->dim=0;
	ret->vdata=NULL;
	ret->allocdata=NULL;
	ret->read.bscale=1.0,ret->read.bzero=0.0;
	ret->curr.bscale=1.0,ret->curr.bzero=0.0;
  }

 return(0);
}

/*****************************************************************************/

char * fits_image_bitpix_cname(int bitpix)
{
 char	*ret;

 switch ( bitpix )
  {	case   8: ret="byte";  break;
	case  16: ret="short"; break;
	case  32: ret="long";  break;
	case -32: ret="float"; break;
	case -64: ret="double";break;
	default : ret=NULL;break;
  }
 
 return(ret);
}

/*****************************************************************************/

static int fits_image_expand(double ***rdata,int sx,int sy,int nx,int ny,double fill)
{
 double	**data;
 size_t	ns;
 int	i,j;

 if ( nx<sx || ny<sy )
	return(-1);
 else if ( nx==sx && ny==sy )
	return(0);

 data=*rdata;

 ns=sizeof(double *)*(size_t)ny+(size_t)nx*(size_t)ny*sizeof(double);
 data=realloc(data,ns);

 for ( i=sy-1 ; i>=0 ; i-- )
  {	memmove(((double *)(data+ny))+i*nx,((double *)(data+sy))+i*sx,sx*sizeof(double));
	for ( j=sx ; j<nx ; j++ )
	 {	*((double *)(data+ny)+i*nx+j)=fill;		}
  }
 for ( i=sy ; i<ny ; i++ )
  {	for ( j=0 ; j<nx ; j++ )
	 {	*((double *)(data+ny)+i*nx+j)=fill;		}
  }

 for ( i=0 ; i<ny ; i++ )
	data[i]=(double *)(data+ny)+i*nx;

 *rdata=data;

 return(0); 
}

static int fits_image_shrink(double ***rdata,int sx,int sy,int nx,int ny)
{
 double	**data;
 size_t	ns;
 int	i;

 if ( nx>sx || ny>sy )
	return(-1);

 else if ( nx==sx && ny==sy )
	return(0);

 data=*rdata;

 for ( i=0 ; i<ny ; i++ )
  {	memmove(((double *)(data+ny))+i*nx,((double *)(data+sy))+i*sx,nx*sizeof(double));
	data[i]=(double *)(data+ny)+i*nx;
  }

 ns=sizeof(double *)*(size_t)ny+(size_t)nx*(size_t)ny*sizeof(double);
 data=realloc(data,ns);

 *rdata=data;
 
 return(0);
}

int fits_image_copy_line(double **data,int tx,int ty,int x0,int y0,double *line,double outer)
{
 int	j;

 if ( y0<0 || y0>=ty )
  {	for ( j=0 ; j<tx ; j++ )
	 {	line[j]=outer;		}
	return(0);
  }
 else
  {	int	l,c;
	l=tx;
	while ( l>0 )
	 {	if ( x0<0 )
		 {	c=-x0;
			if ( c>l )	c=l;
			for ( j=0 ; j<c ; j++ )
			 {	line[j]=outer;		}
			line+=c;
			l-=c;
			x0+=c;
		 }
		else if ( x0>=tx )
		 {	c=l;
			for ( j=0 ; j<c ; j++ )
			 {	line[j]=outer;		}
			line+=c;
			l-=c;
			x0+=c;
		 }
		else
		 {	c=tx-x0;
			if ( c>l )	c=l;
			for ( j=0 ; j<c ; j++ )
			 {	line[j]=data[y0][x0+j];		}
			line+=c;
			l-=c;
			x0+=c;
		 }
	 }
	return(0);
  }
}

int fits_image_trim(fitsimage *img,int x0,int y0,int nx,int ny,double outer)
{
 double	*line;
 int	tx,ty,i,sx,sy;

 if ( img->allocdata != (void *)img->data || img->vdata != (void *)img->data || img->data==NULL )
	return(-2);
 else if ( img->dim != 2 )
	return(-2);
 
 sx=img->sx;
 sy=img->sy;

 if ( sx<=0 || sy<=0 )
	return(-2);

 tx=(nx>sx?nx:sx);
 ty=(ny>sy?ny:sy);

 fits_image_expand(&img->data,sx,sy,tx,ty,outer);

 line=(double *)malloc(tx*sizeof(double));

 if ( x0 != 0 || y0 != 0 )
  {	if ( y0>0 )
	 {	for ( i=0 ; i<ty ; i++ )
		 {	fits_image_copy_line(img->data,tx,ty,x0,i+y0,line,outer);
			memcpy(img->data[i],line,tx*sizeof(double));
		 }
	 }
	else
	 {	for ( i=ty-1 ; i>=0 ; i-- )
		 {	fits_image_copy_line(img->data,tx,ty,x0,i+y0,line,outer);
			memcpy(img->data[i],line,tx*sizeof(double));
		 }
	 }
  }

 free(line);

 fits_image_shrink(&img->data,tx,ty,nx,ny); 

 img->sx=nx;
 img->sy=ny;
 img->naxis[0]=nx;
 img->naxis[1]=ny;

 img->vdata=img->allocdata=(void *)img->data;

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
                                                           
                         
