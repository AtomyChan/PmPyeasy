/*****************************************************************************/
/* fiinfo-pnm.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line user interface to get some statistics from FITS data:	     */
/* Dumping 2D FITS images (fitsimage: dim=2) into PGM or PPM format.         */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>

#include <fits/fits.h>

#include "fi.h"

#include "fitsmask.h"
#include "tensor.h"
#include "common.h"
#include "statistics.h"

#include "fiinfo.h"

/*********************************************************/ /* Old histeq */ /*
			f=fmin;n=fmax-fmin;
			while ( n )
			 {	if ( w>rawdata[f+n/2] )
					f+=n/2+1,n-=n/2+1;
				else
					n=n/2;
			 };
			d=(double)(f-fmin)/(double)(fmax-fmin);
			contrast=0.1;
			slope=4.0;
			d=2*d-1.0;
			d=(((contrast-1.5+0.5*slope)*d*d+(-2.0*contrast+2.5-
				0.5*slope))*d*d+contrast)*d;
			d=(d+1.0)/2.0;
******************************************************************************/

#define		HISTEQBIN	256

typedef struct
 {	double	cnt;
	double	lower,upper,width,area;
	double	finalupper;
 } histodata;

int fitsimage_dump_pnm(fitsimage *img,char **mask,FILE *fw,pnmparam *pp)
{
 double		*rawdata,dmin,dmax,d,w,am[2][2],bv[2],x,y,
		pmin,pmax,det,slope,icept,contrast;
 int		sx,sy,i,j,n,f,np,is_first,lp,p,lsize,fsp,nh,wi,wj;
 unsigned char	*line;
 char		**fsmask;
 histodata	*hd;

 if ( img==NULL || img->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;
 if ( sx<=0 || sy<=0 )			return(1);
 if ( ! pp->is_color )	fprintf(fw,"P5\n");
 else			fprintf(fw,"P6\n");
 fprintf(fw,"%d %d\n",sx,sy);
 fprintf(fw,"255\n");

 np=0;is_first=1;dmin=dmax=0.0;
 if ( pp->scalemethod==1 || pp->minmaxmethod==4 )
	fsp=sx*sy; /* all pixels are needed by 'histeq' && 'percentage' */

 else				fsp=40*40;
 if  ( fsp>sx*sy )		fsp=sx*sy;
 fsmask=fits_mask_create_floyd(sx,sy,10,10*sx*sy/fsp,1); 
 rawdata=malloc(sizeof(double)*fsp);
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( mask != NULL && mask[i][j] )	continue;
		d=img->data[i][j];
		if ( is_first )		dmin=dmax=d,is_first=0;
		else if ( d<dmin )	dmin=d;
		else if ( d>dmax )	dmax=d;
		if ( ! fsmask[i][j] && np<fsp )	rawdata[np]=d,np++;	
	 }
  }
 fits_mask_free(fsmask);
 median(rawdata,np);

 switch ( pp->minmaxmethod )
  {  case MM_MINMAX:
	pmin=dmin,
	pmax=dmax;
	break;
     case MM_MANUAL:
	if ( pp->mmin_set )	pmin=pp->manmin;
	else			pmin=dmin;
	if ( pp->mmax_set )	pmax=pp->manmax;
	else			pmax=dmax;
	break;
     case MM_ZSCALE: case MM_ZMIN: case MM_ZMAX:
	am[0][0]=am[0][1]=am[1][0]=am[1][1]=0.0;
	bv[0]=bv[1]=0.0;
	for ( i=np/4,lp=0 ; i<3*np/4 ; i++ )
	 {	x=(double)(i-np/2);
		y=rawdata[i];
 		am[0][0]+=x*x;
		am[0][1]+=x;
		am[1][1]+=1.0;
		bv[0]+=x*y;
		bv[1]+=y;
		lp++;
	 }
	am[1][0]=am[0][1];
	det=am[0][0]*am[1][1]-am[1][0]*am[0][1];
	slope=(+am[1][1]*bv[0]-am[0][1]*bv[1])/det;
	icept=(-am[1][0]*bv[0]+am[0][0]*bv[1])/det;
	if ( pp->zcontrast>0.0 )	contrast=pp->zcontrast;
	else				contrast=0.25;
	pmin=rawdata[np/2]-(slope/contrast)*(double)(np/2);
	pmax=rawdata[np/2]+(slope/contrast)*(double)(np/2); 
	if ( pp->minmaxmethod==MM_ZMIN )	pmin=dmin;	/* no, zmin */
	if ( pp->minmaxmethod==MM_ZMAX )	pmax=dmax;	/* no, zmax */
	break;
     case MM_PERCENTAGE: /* percentage */
	i=(int)((0.5-pp->percentage/200.0)*(double)np);
	if ( i<0 )	i=0;
	if ( i>=np )	i=np-1;
	j=(int)((0.5+pp->percentage/200.0)*(double)np);
	if ( j<0 )	j=0;
	if ( j>=np )	j=np-1;
	pmin=rawdata[i],
	pmax=rawdata[j];
	break;
     default:
	pmin=dmin,
	pmax=dmax;
	break;
  }

 if ( pmax<=pmin )	pmax=pmin+1.0;

 if ( pp->scalemethod==SCALE_HISTEQU )
  {	int	il,ih;
	double	tcnt,tavg,w,tw,tt,hmin,hmax;
	nh=HISTEQBIN;
	hd=(histodata *)malloc(sizeof(histodata)*nh);
	tcnt=0.0;hmin=rawdata[0];
	for ( i=0 ; i<nh ; i++ )
	 {	il=i*(np-1)/nh;
		ih=(i+1)*(np-1)/nh;
		hd[i].lower=rawdata[il],
		hmax=hd[i].upper=rawdata[ih];
		hd[i].area=(double)(ih-il);
		w=hd[i].upper-hd[i].lower;
		hd[i].cnt=w,tcnt+=w;

	 }
	tavg=100*tcnt/(double)nh;
	tw=0.0;
	for ( i=0 ; i<nh ; i++ )
	 {	w=hd[i].cnt;
		if ( w>tavg )	w=tavg;
		hd[i].width=w;
		tw+=hd[i].width;
	 }
	tt=0.0;
	for ( i=0 ; i<nh ; i++ )
	 {	tt+=hd[i].width;
		hd[i].finalupper=hmin+(hmax-hmin)*tt/tw;
	 }
  }
 else	hd=NULL,nh=0;

 line=(unsigned char *)malloc(sizeof(short)*sx*3);
 for ( i=0 ; i<sy ; i++ )
  {	lsize=0;
	if ( pp->is_flip )	wi=i;
	else			wi=sy-1-i;
	for ( j=0 ; j<sx ; j++ )
	 {	if ( pp->is_mirror )	wj=sx-1-j;
		else			wj=j;
		w=img->data[wi][wj];
		switch ( pp->scalemethod )
		 {   case SCALE_LINEAR:		/* linear */
			d=(w-pmin)/(pmax-pmin);
			break;
		     case SCALE_HISTEQU:	/* histeq */
			f=0;n=nh;
			while ( n>0 )
			 {	if ( w>hd[f+n/2].finalupper )
					f+=n/2+1,n-=n/2+1;
				else
					n=n/2;
			 };
			d=(double)f/(double)nh;
			break;
		     case SCALE_LOG:		/* log */
			if ( 0.0 < pmin && 0.0 < w )
			 {	d=(log(w)-log(pmin))/(log(pmax)-log(pmin));	}
			else if ( w<=0.0 )
				d=0.0;
			else
			 {	d=(w-pmin)/(pmax-pmin);				}
			break;
		     case SCALE_SQRT:		/* sqrt */
			if ( 0.0 < pmin && 0.0 < w )
			 {	d=(sqrt(w)-sqrt(pmin))/(sqrt(pmax)-sqrt(pmin));	}
			else if ( w<=0.0 )
				d=0.0;
			else
			 {	d=(w-pmin)/(pmax-pmin);				}
			break;
		     case SCALE_SQUARED:	/* squared */
			if ( 0.0 < pmin && 0.0 < w )
			 {	d=(w*w-pmin*pmin)/(pmax*pmax-pmin*pmin);	}
			else if ( w<=0.0 )
				d=0.0;
			else
			 {	d=(w-pmin)/(pmax-pmin);				}
			break;
		     default:
			d=(w-pmin)/(pmax-pmin);
			break;
		 }

		d=0.5+(d-1.0+pp->brightness)*pp->contrast;

		if ( pp->is_invert )	d=1.0-d;
		if ( d<0.0 )		d=0.0;
		else if ( d>1.0 )	d=1.0;

		if ( pp->palette==NULL || pp->ncol<=0 )
		 {	p=(int)(d*256.0);
			if ( p>=256 )	p=255;
			if ( p<0 )	p=0;
			if ( ! pp->is_color )
			 {	line[lsize]=p;
				lsize++;
			 }
			else
			 {	line[lsize]=line[lsize+1]=line[lsize+2]=p;
				lsize+=3;
			 }
		 }
		else
		 {	int		bc,r,g,b;
			gradient	*pl;
			bc=(int)(d*(double)pp->ncol);
			if ( bc<0 )		bc=0;
			if ( bc>=pp->ncol )	bc=pp->ncol-1;
			d=d*(double)pp->ncol-(double)bc;
			if ( d<0.0 )	d=0.0;
			if ( d>1.0 )	d=1.0;
			p=(int)(d*256.0);
			if ( p>=256 )	p=255;
			if ( p<0 )	p=0;
			pl=&pp->palette[bc];
			r=(pl->beg.r*(255-p)+pl->end.r*p)/255;
			g=(pl->beg.g*(255-p)+pl->end.g*p)/255;
			b=(pl->beg.b*(255-p)+pl->end.b*p)/255;
			if ( ! pp->is_color )
			 {	line[lsize]=(r+g+b)/3;
				lsize++;
			 }
			else
			 {	line[lsize+0]=r,
				line[lsize+1]=g,
				line[lsize+2]=b;
				lsize+=3;
			 }
		 }
	 }
	fwrite(line,1,lsize,fw);
  }

 if ( hd != NULL )	free(hd);
 if ( rawdata != NULL )	free(rawdata);
 if ( line != NULL )	free(line);
 return(0);
}

/*****************************************************************************/

static int hex_digit(int c)
{
 if ( '0' <= c && c <= '9' )		return(c-'0');
 else if ( 'a' <= c && c <= 'f' )	return(c-'a'+10);
 else if ( 'A' <= c && c <= 'F' )	return(c-'A'+10);
 else					return(-1);
}
int parse_palette(char *pstr,gradient **rpal,int *rncol)
{
 int		ncol,cd,hx[6],sep,i,j,cl;
 gradient	*pal,*pp;

 hx[0]=hx[1]=hx[2]=hx[3]=hx[4]=hx[5]=0;

 ncol=0;pal=NULL;
 while ( *pstr && *pstr != ',' )
  {	cd=0;
	while ( *pstr && *pstr != ',' && hex_digit(*pstr)>=0 && cd<6 )
	 {	hx[cd]=hex_digit(*pstr);
		pstr++;cd++;
	 };
	if ( ! *pstr || *pstr==',' || *pstr==':' )	sep=0;
	else if  ( *pstr=='/' || *pstr=='-' )		sep=1;
	else
	 {	if ( pal != NULL )	free(pal);
		return(1);
	 }
	if ( *pstr && *pstr != ',' )	pstr++;
	if ( ! ( cd==1 || cd==2 || cd==3 || cd==6 ) )
	 {	if ( pal != NULL )	free(pal);
		return(1);
	 }

	pal=(gradient *)realloc(pal,sizeof(int)*6*(ncol+1));
	pp=&pal[ncol];
	if ( cd==1 )
	 {	pp->beg.r=pp->beg.g=pp->beg.b=hx[0]*17;		}
	else if ( cd==2 )
	 {	pp->beg.r=pp->beg.g=pp->beg.b=hx[0]*16+hx[1];	}
	else if ( cd==3 )
	 {	pp->beg.r=hx[0]*17,
		pp->beg.g=hx[1]*17,
		pp->beg.b=hx[2]*17;
	 }
	else if ( cd==6 )
	 {	pp->beg.r=hx[0]*16+hx[1],
		pp->beg.g=hx[2]*16+hx[3],
		pp->beg.b=hx[4]*16+hx[5];
	 }
	if ( ! sep )	pp->end.r=0;
	else		pp->end.r=1;
	pp->end.g=pp->end.b=0;
	ncol++;
  };
 if ( ncol>0 )	pal[ncol-1].end.r=-1;
 else 		
  {	if ( pal != NULL )	free(pal);
	return(1);
  }

 for ( i=0,j=0,cl=0 ; i<ncol ; i++ )
  {	sep=pal[i].end.r;
	if ( cl==0 && sep != 0 )
	 {	pal[j].beg.r=pal[j].end.r=pal[i].beg.r;
		pal[j].beg.g=pal[j].end.g=pal[i].beg.g;
		pal[j].beg.b=pal[j].end.b=pal[i].beg.b;
		j++;
	 }
	else if ( sep<0 )	break;
	else if ( ! sep )
	 {	pal[j].beg.r=pal[i].beg.r,pal[j].end.r=pal[i+1].beg.r;
		pal[j].beg.g=pal[i].beg.g,pal[j].end.g=pal[i+1].beg.g;
		pal[j].beg.b=pal[i].beg.b,pal[j].end.b=pal[i+1].beg.b;
		j++;cl++;
	 }
	else cl=0;
  }
 ncol=j;

 for ( i=0 ; i<ncol ; i++ )
  {	fprintf(stderr,"%.2x%.2x%.2x-%.2x%.2x%.2x\n",
		pal[i].beg.r,pal[i].beg.g,pal[i].beg.b,
		pal[i].end.r,pal[i].end.g,pal[i].end.b);
  } 

 *rpal=pal,
 *rncol=ncol;

 return(0);
}

/*****************************************************************************/
