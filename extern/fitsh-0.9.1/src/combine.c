/*****************************************************************************/
/* combine.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions related to combining more FITS images.			     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <fits/fits.h>

#include "fi.h"

#include "fitsmask.h"
#include "io/scanarg.h"
#include "statistics.h"
#include "tensor.h"
#include "combine.h"

/*****************************************************************************/

int combine_parse_mode(char *modstr,compar *cp)
{
 int	i;
 double	sigma;

 if ( modstr==NULL )	return(0);
 sigma=-1.0;

 i=scanpar(modstr,SCANPAR_DEFAULT,
	"avg|mean|average:"	SNf(COM_MODE_AVG),&cp->mode,
	"med|median:"		SNf(COM_MODE_MED),&cp->mode,
	"rejavg:"		SNf(COM_MODE_REJ_AVG),&cp->mode,
	"rej|rejmed|rejection:"	SNf(COM_MODE_REJ_MED),&cp->mode,
	"add|sum:"		SNf(COM_MODE_SUM),&cp->mode,
	"sqs|squaresum:"	SNf(COM_MODE_SQSUM),&cp->mode,
	"sct|scatter|stddev:"	SNf(COM_MODE_SCT),&cp->mode,
	"min|minimum:"		SNf(COM_MODE_MIN),&cp->mode,
	"max|maximum:"		SNf(COM_MODE_MAX),&cp->mode,
	"iterations:%d",&cp->niter,
	"lower:%g",&cp->lower,
	"upper:%g",&cp->upper,
	"sigma:%g",&sigma,
	"or:%SN0f",&cp->logicalmethod,
	"and:%SN1f",&cp->logicalmethod,
	"ignorenegative:%N1f",&cp->ignore_flag,
	"ignorezero:%N2f",&cp->ignore_flag,
	"ignorepositive:%N4f",&cp->ignore_flag,
	NULL);

 if ( sigma>0.0 )
	cp->lower=cp->upper=sigma;

 if ( i )	return(1);
 else		return(0);
}

/*****************************************************************************/

double do_ordered_rejection(double *points,int n,double m,double w1,double w2,
int niter,double lower,double upper,int is_median)
{
 double	s,s2,lim;
 int	r;

 if ( n<=0 || points==NULL )
	return(0.0);

 for ( ; niter > 0 ; niter-- )
  {	
	s2=w2+(double)n*m*(m-2.0*w1);
	if ( s2 <= 0.0 )
		return(m);
	else
		s=sqrt(s2/(double)n);

	lim=m-s*lower;
	r=0;
	while ( n>0 && *points < lim )
	 {	w1-=*points,
		w2-=(*points)*(*points);
		points++,n--;
		r++;
	 }
	lim=m+s*upper;
	while ( n>0 && points[n-1] > lim )
	 {	w1-=points[n-1],
		w2-=points[n-1]*points[n-1];
		n--;
		r++;
	 }
	if ( r<=0 )
		return(m);
	else if ( n<=2 )
		return(m);
	else if ( ! is_median )
		m=w1/(double)n;
	else if ( n%2==1 )
		m=points[n/2];
	else
		m=0.5*(points[n/2-1]+points[n/2]);

  }

 return(m);
}

double combine_points(double *points,int n,compar *cp)
{
 int	i;
 double w,w1,w2,m;
 
 median(points,n);

 if ( n==0 )				/* All pixels to combine are false.  */
	return(0.0);		

 else if ( n==1 )			/* One pixel.			     */
	return(points[0]);
 else if ( cp->mode==COM_MODE_MIN )
  {	w=points[0];
	for ( i=1 ; i<n ; i++ )
	 {	if ( points[i]<w )	w=points[i];	}
	return(w);
  }
 else if ( cp->mode==COM_MODE_MAX )
  {	w=points[0];
	for ( i=1 ; i<n ; i++ )
	 {	if ( w<points[i] )	w=points[i];	}
	return(w);
  }
 
 else if ( n==2 )			/* Two, it's the mean/median of them */
	return(0.5*(points[0]+points[1]));

 else if ( cp->mode==COM_MODE_AVG )
  {	w=0.0;
	for ( i=0 ; i<n ; i++ )
	 {	w+=points[i];		}
	w=w/(double)n;
	return(w);
  }
 else if ( cp->mode==COM_MODE_SUM )
  {	w=0.0;
	for ( i=0 ; i<n ; i++ )
	 {	w+=points[i];		}
	return(w);
  }
 else if ( cp->mode==COM_MODE_SQSUM )
  {	w=0.0;
	for ( i=0 ; i<n ; i++ )
	 {	w+=points[i]*points[i];		}
	return(w);
  }
 else if ( cp->mode==COM_MODE_SCT )
  {	double	s,s2;
	s=s2=0.0;
	for ( i=0 ; i<n ; i++ )
	 {	w=points[i];
		s+=w,s2+=w*w;
	 }
	s/=(double)n,s2/=(double)n;
	s2=s2-s*s;
	if ( s2>0.0 )	return(sqrt(s2));
	else		return(0.0);
  }
 else if ( cp->mode==COM_MODE_MED )
  {	if ( n%2==1 )	w=points[n/2];
	else		w=0.5*(points[n/2-1]+points[n/2]);
	return(w);
  }
 else if ( cp->mode==COM_MODE_REJ_DEPRECATED )	/* deprecated... see below   */
  {	if ( n==3 )			/* Three, also...		     */
		return(points[1]);
	else if ( n<=5 )		/* Average of all except the first   */
	 {	w=0.0;			/* and the last one.		     */
		for ( i=1 ; i<n-1 ; i++ )
		 {	w+=points[i];		}
		w=w/(double)(n-2);
		return(w);
	 }
	else				/* Average all of them except the    */
	 {	w=0.0;			/* first two and last two ones.	     */
		for ( i=2 ; i<n-2 ; i++ )
		 {	w+=points[i];		}
		w=w/(double)(n-4);
		return(w);
	 }
  }
 else if ( cp->mode==COM_MODE_REJ_AVG )
  {	w1=0.0;
	w2=0.0;
	for ( i=0 ; i<n ; i++ )
	 {	w1+=points[i];
		w2+=points[i]*points[i];
	 }
	w=do_ordered_rejection(points,n,w1/(double)n,w1,w2,cp->niter,cp->lower,cp->upper,0);
	return(w);
  }
 else if ( cp->mode==COM_MODE_REJ_MED )
  {	w1=0.0;
	w2=0.0;
	for ( i=0 ; i<n ; i++ )
	 {	w1+=points[i];
		w2+=points[i]*points[i];
	 }
	if ( n%2==1 )	m=points[n/2];
	else		m=0.5*(points[n/2-1]+points[n/2]);
	w=do_ordered_rejection(points,n,m,w1,w2,cp->niter,cp->lower,cp->upper,1);
	return(w);
  }
 else					/* Unimplemented method...?!	     */
	return(0.0);
}

int combine_lines(double **lines,int n,int sx,double *out,
	compar *cp,char **wmask,char *outmask)
{
 int	i,j,k,mw,ow;
 double	*warr,w;

 warr=(double *)tensor_alloc_1d(double,n);

 for ( i=0 ; i<sx ; i++ )
  {	k=0;
	ow=0;
	if ( wmask != NULL )
	 {	for ( j=0 ; j<n ; j++ )
		 {	w =lines[j][i];
			mw=wmask[i][j];
			ow|=mw;
			if ( (cp->ignore_flag&COM_IGNORE_NEGATIVE) && w<0.0 )
				mw|=MASK_FAULT;
			if ( ! mw )
				warr[k]=w,k++;
		 }
	 }
	else
	 {	for ( j=0 ; j<n ; j++ )
		 {	w =lines[j][i];
			mw=0;
			if ( (cp->ignore_flag&COM_IGNORE_NEGATIVE) && w<0.0 )
				mw|=MASK_FAULT;
			if ( ! mw )
				warr[k]=w,k++;
		 }
	 }
	if ( outmask != NULL )
	 {	if ( ! cp->logicalmethod && k>0 )		/* --logical-or	 */
			out[i]=combine_points(warr,k,cp),outmask[i]=0;
		else if ( cp->logicalmethod && k==n )		/* --logical-and */
			out[i]=combine_points(warr,k,cp),outmask[i]=0;
		else
		 {	out[i]=0.0;
			outmask[i]=(ow?ow:MASK_FAULT);
		 }
	 }
	else
	 {	if ( ! cp->logicalmethod && k>0 )		/* --logical-or	 */
			out[i]=combine_points(warr,k,cp);
		else if ( cp->logicalmethod && k==n )		/* --logical-and */
			out[i]=combine_points(warr,k,cp);
		else
			out[i]=0.0;
	 }
  }

 tensor_free(warr);
 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int combine_images_from_files_v1(comimg *inputs,int n,fits *outimg,
	compar *cp,char **inmask,char **outmask,
	presubdata *presubs,int npresub,size_t maxmem)
{
 int		i,j,k,l,p,sx,sy,ny,py,cny,nreal,nz,cz;
 double		**lines;
 char		***linemasks,**wmask;
 presubdata	*wp;
 fitsimage	*img;
 
 sx=outimg->i.sx,
 sy=outimg->i.sy;

 for ( i=0,wp=presubs ; i<npresub ; i++,wp++ )
  {	if ( wp->img==NULL || wp->img->data==NULL || wp->img->sx != sx || wp->img->sy != sy )
		npresub=i;
  }

 nreal=0;
 for ( j=0 ; j<n ; j++ )	/* tbd: more sanity checks here? */
  {	img=&inputs[j].img->i;
	if ( img->dim==2 )	nreal++;
	else			nreal+=img->naxis[2];
  }

 if ( maxmem>0 )
 	ny=(int)(maxmem/(sizeof(double)*(size_t)sx*nreal));
 else
	ny=(sizeof(double)*sy)/nreal;

/*
 fprintf(stderr,"maxmem=%d nreal=%d ny=%d\n",(int)(maxmem/1048576),nreal,ny);
*/

 if ( ny<=0 )	ny=1;
 if ( ny>=sy )	ny=sy;

 lines=(double **)tensor_alloc_2d(double,sx,nreal);
 linemasks=(char ***)tensor_alloc_3d(char,sx,ny,nreal); 
 wmask=(char **)tensor_alloc_2d(char,nreal,sx);

 py=-1,cny=0;

 for ( i=0 ; i<sy ; i++ )
  {	if ( py<0 )
	 {	py=i;cny=ny;
		if ( py+cny>sy )	cny=sy-py;
		cz=0;
		for ( j=0 ; j<n ; j++ )
		 {	img=&inputs[j].img->i;
			if ( img->dim==2 )	nz=1;
			else			nz=img->naxis[2];
			for ( l=0 ; l<nz ; l++ )
			 {	for ( k=0 ; k<cny ; k++ )
				 {	memset(linemasks[cz+l][k],0,sx);	}
				fits_mask_mask_more_line(linemasks[cz+l],&inputs[j].img->header,sx,sy,py,cny,NULL);
			 }
			cz+=nz;
		 }
	 }
	cz=0;
	for ( j=0 ; j<n ; j++ )
	 {	img=&inputs[j].img->i;
		if ( img->dim==2 )	nz=1;
		else			nz=img->naxis[2];

		/* currently read image lines must be rescaled... */
		if ( img->vdata == NULL )
		 {	fits_read_image_line(inputs[j].fr,sx,img->bit,lines[cz]);
		 	for ( k=0 ; k<sx ; k++ )
			 {	lines[cz][k]=lines[cz][k]*img->curr.bscale+img->curr.bzero;
				wmask[k][cz]=linemasks[cz][i-py][k];
			 }
		 }
		/* pre-read images are assumed to have been rescaled! */
		else if ( img->dim==2 )
		 {	for ( k=0 ; k<sx ; k++ )
			 {	lines[cz][k]=img->data[i][k];
				wmask[k][cz]=linemasks[cz][i-py][k];
			 }
		 }
		else	/* assuming dim == 3 ... ?! */
		 {	double	***zdata;
			zdata=(double ***)img->vdata;
			for ( l=0 ; l<nz ; l++ )
			 {	for ( k=0 ; k<sx ; k++ )
				 {	lines[cz+l][k]=zdata[l][i][k];
					wmask[k][cz+l]=linemasks[cz][i-py][k];
				 }
			 }
		 }

		cz+=nz;
	 }
	for ( p=0,wp=presubs ; p<npresub ; p++,wp++ )
	 {	for ( j=0 ; j<nreal ; j++ )
		 {	for ( k=0 ; k<sx ; k++ )
			 {	lines[j][k]-=wp->img->data[i][k]*wp->scale;		}
		 }
	 }

	if ( inmask != NULL )
	 {	for ( k=0 ; k<sx ; k++ )
		 {	if ( inmask[i][k] )
			 {	for ( j=0 ; j<nreal ; j++ )
				 {	wmask[k][j] |= inmask[i][k];	}
			 }
		 }
	 }
	if ( outmask != NULL )	combine_lines(lines,nreal,sx,outimg->i.data[i],cp,wmask,outmask[i]);
	else			combine_lines(lines,nreal,sx,outimg->i.data[i],cp,wmask,NULL);
	if ( py+cny>=i+1 )	py=-1;
  }

 tensor_free(wmask);
 tensor_free(linemasks);
 tensor_free(lines);
 
 return(0);
}

int combine_images_from_files_v2(comimg *inputs,int n,fits *outimg,
	compar *cp,char **inmask,char **outmask,
	presubdata *presubs,int npresub,size_t maxmem)
{
 int		i,j,k,l,p,sx,sy,ny,nreal,nz,cz,y0,dy;
 double		***clines;
 char		***cmasks,**wmask;
 presubdata	*wp;
 fitsimage	*img;
 
 sx=outimg->i.sx,
 sy=outimg->i.sy;

 for ( i=0,wp=presubs ; i<npresub ; i++,wp++ )
  {	if ( wp->img==NULL || wp->img->data==NULL || wp->img->sx != sx || wp->img->sy != sy )
		npresub=i;
  }

 nreal=0;
 for ( j=0 ; j<n ; j++ )	/* tbd: more sanity checks here? */
  {	img=&inputs[j].img->i;
	if ( img->dim==2 )	nreal++;
	else			nreal+=img->naxis[2];
  }

 if ( maxmem>0 )
 	ny=(int)(maxmem/((sizeof(double)+1)*(size_t)sx*nreal));
 else
	ny=sy/nreal;

/*
 fprintf(stderr,"maxmem=%d nreal=%d ny=%d\n",(int)(maxmem/1048576),nreal,ny);
*/

 if ( ny<=0 )	ny=1;
 if ( ny>=sy )	ny=sy;

 clines=(double ***)tensor_alloc_3d(double,sx,nreal,ny);
 cmasks=(char ***  )tensor_alloc_3d(char  ,sx,ny,nreal); 
 wmask =(char **   )tensor_alloc_2d(char  ,nreal,sx);

 for ( y0=0 ; y0<sy ;  )
  {	dy=sy-y0;
	if ( dy>ny )	dy=ny;

	cz=0;
	for ( j=0 ; j<n ; j++ )
	 {	img=&inputs[j].img->i;
		if ( img->dim==2 )	nz=1;
		else			nz=img->naxis[2];
		for ( l=0 ; l<nz ; l++ )
		 {	for ( k=0 ; k<dy ; k++ )
			 {	memset(cmasks[cz+l][k],0,sx);	}
			fits_mask_mask_more_line(cmasks[cz+l],&inputs[j].img->header,sx,sy,y0,dy,NULL);
		 }
		cz+=nz;
	 }

	cz=0;
	for ( j=0 ; j<n ; j++ )
	 {	img=&inputs[j].img->i;
		if ( img->dim==2 )	nz=1;
		else			nz=img->naxis[2];
		if ( img->vdata==NULL )
		 {	for ( k=0 ; k<dy ; k++ )
			 {	fits_read_image_line(inputs[j].fr,sx,img->bit,clines[k][cz]);
				if ( img->curr.bscale != 1.0 || img->curr.bzero != 0.0 )
				 {	for ( l=0 ; l<sx ; l++ )
					 {	clines[k][cz][l]=clines[k][cz][l]*img->curr.bscale+img->curr.bzero;	}
				 }
			 }
		 }
		else if ( img->dim==2 )
		 {	for ( k=0 ; k<dy ; k++ )
			 {	memcpy(clines[k][cz],img->data[y0+k],sx*sizeof(double));	}
		 }
		else
		 {	double	***zdata;
			zdata=(double ***)img->vdata;
			for ( l=0 ; l<nz ; l++ )
			 {	for ( k=0 ; k<dy ; k++ )
				 {	memcpy(clines[k][cz+l],zdata[l][y0+k],sx*sizeof(double));	}
			 }
		 }
		cz+=nz;
	 }

	for ( p=0,wp=presubs ; p<npresub ; p++,wp++ )
	 {	for ( j=0 ; j<nreal ; j++ )
		 {	for ( k=0 ; k<dy ; k++ )
			 {	for ( l=0 ; l<sx ; l++ )
				 {	clines[k][j][l]-=wp->img->data[y0+k][l]*wp->scale;	}
			 }
		 }
	 }

	for ( k=0 ; k<dy ; k++ )
	 {	if ( inmask != NULL )
		 {	for ( j=0 ; j<nreal ; j++ )
			 {	for ( l=0 ; l<sx ; l++ )
				 {	wmask[l][j]=cmasks[j][k][l] | inmask[y0+k][l];	}
			 }
		 }
		else
		 {	for ( j=0 ; j<nreal ; j++ )
			 {	for ( l=0 ; l<sx ; l++ )
				 {	wmask[l][j]=cmasks[j][k][l];	}
			 }

		 }
		if ( outmask != NULL )
			combine_lines(clines[k],nreal,sx,outimg->i.data[y0+k],cp,wmask,outmask[y0+k]);
		else
			combine_lines(clines[k],nreal,sx,outimg->i.data[y0+k],cp,wmask,NULL);
	 }		

	y0+=dy;
  }

 tensor_free(wmask);
 tensor_free(cmasks);
 tensor_free(clines);
 
 return(0);
}

int combine_images_from_files(comimg *inputs,int n,fits *outimg,
	compar *cp,char **inmask,char **outmask,
	presubdata *presubs,int npresub,size_t maxmem)
{
 int	r;
 r=combine_images_from_files_v2(inputs,n,outimg,cp,inmask,outmask,presubs,npresub,maxmem);
 return(r);
}

/*****************************************************************************/

int combine_cleanup(comimg *inputs,int ninput)
{
 int    i;
 comimg	*ci;
 
 for ( i=ninput-1 ; i>=0 ; i-- )
  {	ci=&inputs[i];
	if ( ci->fr  != NULL )	fclose(ci->fr);
	if ( ci->img != NULL )	fits_free(ci->img);
  }

 return(0);
}

/*****************************************************************************/
                                
