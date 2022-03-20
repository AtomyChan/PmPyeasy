/*****************************************************************************/
/* transform.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Library for parsing and writing transformation files and perform some     */
/* related calculations (concerning to inverse transformations, Jacobi's     */
/* determinant and so on...).						     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "math/poly.h"
#include "io/iof.h"
#include "io/tokenize.h"

#include "transform.h"

static char	*dxfitstr="dxfit",
		*dyfitstr="dyfit";

/*****************************************************************************/

int transformation_free(transformation *tf)
{
 int	i;
 if ( tf->vfits != NULL )
  {	for ( i=0 ; i<tf->nval ; i++ )
	 {	if ( tf->vfits[i] != NULL )	free(tf->vfits[i]);	}
	free(tf->vfits);
  }
 return(0);
}

/*****************************************************************************/

static int transformation_add_vfit(transformation *tf,int v,char *args)
{
 double	*fit;
 int	i,nvar;

 nvar=(tf->order+1)*(tf->order+2)/2;

 if ( tf->vfits==NULL )
  {	tf->vfits=(double **)malloc((v+1)*sizeof(double *));
	for ( i=0 ; i<v+1 ; i++ )
	 {	tf->vfits[i]=NULL;		}
	tf->nval=v+1;
  }
 else if ( v>=tf->nval )
  {	tf->vfits=(double **)realloc(tf->vfits,(v+1)*sizeof(double *));
	for ( i=tf->nval ; i<v+1 ; i++ )
	 {	tf->vfits[i]=NULL;		}
	tf->nval=v+1;
  }

 if ( tf->vfits[v] != NULL )	return(1);

 fit=(double *)malloc(sizeof(double)*nvar);

 for ( i=0 ; i<nvar ; i++ )
  {	if ( sscanf(args,"%lg",&fit[i])<1 )
	 {	free(fit);return(1);		}
	if ( ! isfinite(fit[i]) )
	 {	free(fit);return(1);		}
	if ( i<nvar-1 )
	 {	while ( *args && *args != ',' )	args++;
		if ( *args != ',' )
		 {	free(fit);return(1);		}
		else	args++;
	 }
  }
	
 tf->vfits[v]=fit;

 return(0);
}

static int transformation_check_vfit(transformation *tf)
{
 int	i;
 if ( tf->vfits==NULL || tf->nval<=0 )	return(1);
 for ( i=0 ; i<tf->nval ; i++ )
  {	if ( tf->vfits[i]==NULL )	return(1);	}
 return(0);
}

int transformation_read_data(FILE *fr,transformation *tf)
{
 char	*rbuff,*cmd[4];
 int	n,type,order,nvar,v;
 type=0;nvar=-1;
 rbuff=NULL;

 tf->bshx=tf->bshy=0.0;	/* these things are global for the transformation... */

 while ( ! feof(fr) )
  {	if ( rbuff != NULL )	free(rbuff);
	rbuff=freadline(fr);
	if ( rbuff==NULL )	break;
 	remove_spaces_and_comments(rbuff);
	if ( strlen(rbuff)==0 )	continue;
	n=tokenize_char(rbuff,cmd,'=',2);
	if ( n<2 )
 	 {	free(rbuff);return(1);		}
	if ( strcmp(cmd[0],"type")==0 || ! type )
	 {	if ( strcmp(cmd[1],"polynomial")==0 )
		 {	tf->type=type=TRANS_POLYNOMIAL;
			tf->ox=tf->oy=0.0;	/* These parameters are optional	*/
			tf->scale=1.0;		/* default values: ox=oy=0, scale=1	*/
			tf->vfits=NULL;
			tf->nval=0;
			nvar=-1;
		 }
		else
		 {	free(rbuff);return(1);	}
	 }
	else if ( ! type )
	 {	free(rbuff);return(5);	}

	else if ( strcmp(cmd[0],"order")==0 )
	 {	if ( sscanf(cmd[1],"%d",&order)<1 )
		 {	free(rbuff);return(5);	}
		tf->order=order;
		nvar=(order+1)*(order+2)/2;
	 }
	else if ( strcmp(cmd[0],"offset")==0 )
	 {	if ( sscanf(cmd[1],"%lg,%lg",&tf->ox,&tf->oy)<2 )
		 {	free(rbuff);return(5);		}
		if ( ! isfinite(tf->ox) || ! isfinite(tf->oy) )
		 {	free(rbuff);return(5);		}
	 }
	else if ( strcmp(cmd[0],"scale")==0 )
	 {	if ( sscanf(cmd[1],"%lg",&tf->scale)<1 )
		 {	free(rbuff);return(5);		}
		if ( ! isfinite(tf->scale) )
		 {	free(rbuff);return(5);		}
	 }
	else if ( strcmp(cmd[0],"basisshift")==0 )
	 {	if ( sscanf(cmd[1],"%lg,%lg",&tf->bshx,&tf->bshy)<2 )
		 {	free(rbuff);return(5);		}
		if ( ! isfinite(tf->bshx) || ! isfinite(tf->bshy) )
		 {	free(rbuff);return(5);		}
	 }

	else if ( strcmp(cmd[0],dxfitstr)==0 || strcmp(cmd[0],dyfitstr)==0 )
	 {	if ( strcmp(cmd[0],dxfitstr)==0 )	v=0;
		else					v=1;
		if ( nvar<=0 )
		 {	free(rbuff);return(5);		}
		if ( transformation_add_vfit(tf,v,cmd[1]) )
		 {	free(rbuff);return(5);		}
	 }
	else if ( sscanf(cmd[0],"vfit:%d",&v)==1 )
	 {	if ( nvar<=0 || v<=0 )
		 {	free(rbuff);return(5);		}
		v--;
		if ( transformation_add_vfit(tf,v,cmd[1]) )
		 {	free(rbuff);return(5);		}
	 }
	else 
	 {	free(rbuff);return(5);	}
  }

 if ( transformation_check_vfit(tf) )	return(1);
 return(0);
}

int strlcmp(char *b1,char *b2)
{
 int	l;
 l=strlen(b2);
 if ( memcmp(b1,b2,l) )		return(1);
 else				return(0);
}

int transformation_parse_params(char *params,transformation *tf)
{
 int	l,ac,polytype,v;
 char	*sp,*args;

 tf->type=0;tf->order=-1;
 tf->bshx=tf->bshy=0.0;
 tf->ox=tf->oy=0.0;tf->scale=1.0;
 tf->vfits=NULL;tf->nval=0;
 polytype=-1;

 while ( *params )
  {	for ( l=0,sp=params ; sp[l] && sp[l] != ',' && sp[l] != '=' ; )	l++;
	if ( sp[l]==',' )	ac=0,args=NULL,params+=l+1;
	else if ( ! sp[l] )	ac=0,args=NULL,params+=l;
	else
	 {	params+=l+1;
		args=params;ac=1;
		while ( *params && ! isalpha(*params) )
		 {	if ( *params==',' )	ac++;
			params++;
		 }
		if ( isalpha(*params) )	ac--;
	 }

	v=-1;
	if ( strlcmp(sp,"polynomial")==0 && ac==0 )
	 {	tf->type=TRANS_POLYNOMIAL;
		tf->vfits=NULL;
		tf->nval=0;
	 }
	else if ( strlcmp(sp,"order")==0 && ac==1 )
	 {	if ( sscanf(args,"%d",&tf->order)<1 )	return(1);
		if ( tf->order<0 )			return(1);
	 }
	else if ( strlcmp(sp,"scale")==0 && ac==1 )
	 {	if ( sscanf(args,"%lg",&tf->scale)<1 )	return(1);
		if ( ! isfinite(tf->scale) )		return(1);
	 }
	else if ( strlcmp(sp,"offset")==0 && ac==2 )
	 {	if ( sscanf(args,"%lg,%lg",&tf->ox,&tf->oy)<2 )	return(1);
		if ( ! isfinite(tf->ox) || ! isfinite(tf->oy) )	return(1);
	 }
	else if ( strlcmp(sp,"basisshift")==0 && ac==2 )
	 {	if ( sscanf(args,"%lg,%lg",&tf->bshx,&tf->bshy)<2 )	return(1);
		if ( ! isfinite(tf->bshx) || ! isfinite(tf->bshy) )	return(1);
	 }

	else if ( strlcmp(sp,dxfitstr)==0 || strlcmp(sp,dyfitstr)==0 || sscanf(sp,"vfit:%d",&v)==1 )
	 {	if ( strlcmp(sp,dxfitstr)==0 )		v=0;
		else if ( strlcmp(sp,dyfitstr)==0 )	v=1;
		else					v--;
		if ( v<0 )	return(1);
		if ( transformation_add_vfit(tf,v,args) )	return(1);
	 }
	else	
		return(1);
  };

 if ( transformation_check_vfit(tf) )	return(1);

 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static int transformation_write_coeffs(FILE *fw,double *fit,int nvar)
{
 int	i;
 for ( i=0 ; i<nvar ; i++ )
  {	fprintf(fw,"%12g",fit[i]);
	if ( i<nvar-1 )	fprintf(fw,", ");
  }
 fprintf(fw,"\n");
 return(0);
}
static int arch_is_swapped(void)
{
 short  a;
 a=1;
 if ( *(char *)(&a) == 1 )	return(1);
 else				return(0);
}

static int sprint_ieee_32(char *buff,double x)
{
 float	f;
 int	i;
 char	*b;

 f=(float)x;
 b=(char *)(&f);
 if ( arch_is_swapped() )
  {	for ( i=0 ; i<4 ; i++ )
	 {	sprintf(buff+i*2,"%.2X",*(unsigned char *)(b+(3-i)));	}
  }
 else
  {	for ( i=0 ; i<4 ; i++ )
	 {	sprintf(buff+i*2,"%.2X",*(unsigned char *)(b+i));	}
  }
 return(0);
}
static int sprint_ieee_64(char *buff,double x)
{
 double	f;
 int	i;
 char	*b;

 f=(double)x;
 b=(char *)(&f);
 if ( arch_is_swapped() )
  {	for ( i=0 ; i<8 ; i++ )
	 {	sprintf(buff+i*2,"%.2X",*(unsigned char *)(b+(7-i)));	}
  }
 else
  {	for ( i=0 ; i<8 ; i++ )
	 {	sprintf(buff+i*2,"%.2X",*(unsigned char *)(b+i));	}
  }
 return(0);
}

static int transformation_write_coeffs_ieee(FILE *fw,double *fit,int nvar,int flag)
{
 int	i;
 char	ibuff[32];
 for ( i=0 ; i<nvar ; i++ )
  {	if ( flag & TRANS_WR_IEEE_32 )
	 {	sprint_ieee_32(ibuff,fit[i]);
  		fprintf(fw,"%s",ibuff);
	 }
	else if ( flag & TRANS_WR_IEEE_64 )
	 {	sprint_ieee_64(ibuff,fit[i]);
  		fprintf(fw,"%s",ibuff);
	 }
	if ( i<nvar-1 )	fprintf(fw,", ");
  }
 fprintf(fw,"\n");
 return(0);
}

int transformation_write_data(FILE *fw,transformation *tf,int flags)
{
 int	j,nvar;

 if ( tf->type==TRANS_POLYNOMIAL )
  {	nvar=(tf->order+1)*(tf->order+2)/2;
	if ( flags&TRANS_WR_COMMENT )
		fprintf(fw,"# Type: polynomial of order=%d (number of coefficients: %d)\n",tf->order,nvar);
	fprintf(fw,"type = polynomial\n");
	fprintf(fw,"order = %d\n",tf->order);
	if ( flags&TRANS_WR_COMMENT )
		fprintf(fw,"# Initial transformation of (x_img,y_img):\n");
	fprintf(fw,"offset = %g, %g\n",tf->ox,tf->oy);
	fprintf(fw,"scale = %g\n",tf->scale);
	fprintf(fw,"basisshift = %g, %g\n",tf->bshx,tf->bshy);

	if ( tf->nval==2 && (flags&TRANS_WR_DXDY) )
	 {	if ( flags&TRANS_WR_COMMENT )
			fprintf(fw,"# Coefficients of the x fit: \n");
		fprintf(fw,"%s= ",dxfitstr);
		transformation_write_coeffs(fw,tf->vfits[0],nvar);
		if ( flags&(TRANS_WR_IEEE_32|TRANS_WR_IEEE_64) )
		 {	fprintf(fw,"%s-ieee= ",dxfitstr);
			transformation_write_coeffs_ieee(fw,tf->vfits[0],nvar,flags);
		 }
		if ( flags&TRANS_WR_COMMENT )
			fprintf(fw,"# Coefficients of the y fit: \n");
		fprintf(fw,"%s= ",dyfitstr);
		transformation_write_coeffs(fw,tf->vfits[1],nvar);
		if ( flags&(TRANS_WR_IEEE_32|TRANS_WR_IEEE_64) )
		 {	fprintf(fw,"%s-ieee= ",dyfitstr);
			transformation_write_coeffs_ieee(fw,tf->vfits[1],nvar,flags);
		 }
	 }
	else
	 {	if ( flags&TRANS_WR_COMMENT )
			fprintf(fw,"# Coefficients: \n");
		for ( j=0 ; j<tf->nval ; j++ )
		 {	fprintf(fw,"vfit:%d= ",j+1);
			transformation_write_coeffs(fw,tf->vfits[j],nvar);
			if ( flags&(TRANS_WR_IEEE_32|TRANS_WR_IEEE_64) )
			 {	fprintf(fw,"vfit:%d-ieee= ",j+1);
				transformation_write_coeffs_ieee(fw,tf->vfits[j],nvar,flags);
			 }
		 }
	 }
  }

 return(0);
}

/*****************************************************************************/

int transformation_check_if_null(transformation *tf)
{
 int	nvar,i,j;

 if ( tf->type == TRANS_POLYNOMIAL )
  {	if ( tf->vfits==NULL )	return(-1);
	nvar=(tf->order+1)*(tf->order+2)/2;
	for ( i=0 ; i<tf->nval ; i++ )
	 {	if ( tf->vfits[i]==NULL )	return(-1);
		for ( j=0 ; j<nvar ; j++ )
		 {	if ( ! isfinite(tf->vfits[i][j]) )
				return(-1);
			if ( tf->vfits[i][j] != 0.0 )	
				return(0);
		 }
	 }
	return(1);
  }
 else
	return(0);

}

/*****************************************************************************/

int transformation_get_jacobi(transformation *tf,double **rjxx,double **rjxy,double **rjyx,double **rjyy)
{
 int	jnvar,i,j,k;
 double	scale,*jxx,*jxy,*jyx,*jyy;

 jnvar=(tf->order+0)*(tf->order+1)/2;
 scale=tf->scale;

 jxx=(double *)malloc(sizeof(double)*jnvar);
 jxy=(double *)malloc(sizeof(double)*jnvar);
 jyx=(double *)malloc(sizeof(double)*jnvar);
 jyy=(double *)malloc(sizeof(double)*jnvar);

 for ( i=0 ; i<jnvar ; i++ )
  {	jxx[i]=jxy[i]=jyx[i]=jyy[i]=0.0;		}
 for ( i=0,k=0 ; i<=tf->order-1 ; i++ )
  {	for ( j=0 ; j<=i ; j++,k++ )
         {	jxx[k]+=tf->vfits[0][k+i+1]/scale;
		jxy[k]+=tf->vfits[0][k+i+2]/scale;
		jyx[k]+=tf->vfits[1][k+i+1]/scale;
		jyy[k]+=tf->vfits[1][k+i+2]/scale;
	 }
  }

 *rjxx=jxx,
 *rjxy=jxy,
 *rjyx=jyx,
 *rjyy=jyy;

 return(0);
}


int transformation_eval_normal_2d(double x,double y,transformation *tf,double *rx,double *ry)
{
 *rx  =eval_2d_poly(x,y,tf->order,tf->vfits[0],tf->ox,tf->oy,tf->scale);
 *ry  =eval_2d_poly(x,y,tf->order,tf->vfits[1],tf->ox,tf->oy,tf->scale);
 return(0);
}

#define		EPS	1e-10

int transformation_eval_invert_2d(double x,double y,transformation *tf,double *rx,double *ry,double *jxx,double *jxy,double *jyx,double *jyy)
{
 double	wx,wy,mxx,mxy,myx,myy,det,imxx,imxy,imyx,imyy,x0,y0,px0,py0,dx,dy;
 double	*px,*py;
 int	i,n;

 if ( tf->order < 1 )	return(1);

 px=tf->vfits[0],
 py=tf->vfits[1];
 
 wx=x-px[0];
 wy=y-py[0];
 mxx=px[1],mxy=px[2],
 myx=py[1],myy=py[2];
 det=1/(mxx*myy-mxy*myx);
 imxx=+myy*det;
 imxy=-mxy*det;
 imyx=-myx*det;
 imyy=+mxx*det;
 x0=imxx*wx+imxy*wy;
 y0=imyx*wx+imyy*wy;
 
 n=100;px0=x0;py0=y0;
 for ( i=0 ; i<n && tf->order>=2 ; i++ )
  {	mxx=eval_2d_poly(x0,y0,tf->order-1,jxx,tf->ox,tf->oy,tf->scale);
	mxy=eval_2d_poly(x0,y0,tf->order-1,jxy,tf->ox,tf->oy,tf->scale);
	myx=eval_2d_poly(x0,y0,tf->order-1,jyx,tf->ox,tf->oy,tf->scale);
	myy=eval_2d_poly(x0,y0,tf->order-1,jyy,tf->ox,tf->oy,tf->scale);
	det=1/(mxx*myy-mxy*myx);
	imxx=+myy*det;
	imxy=-mxy*det;
	imyx=-myx*det;
	imyy=+mxx*det;
	wx=eval_2d_poly(x0,y0,tf->order,px,tf->ox,tf->oy,tf->scale)-x;
	wy=eval_2d_poly(x0,y0,tf->order,py,tf->ox,tf->oy,tf->scale)-y;
	dx=imxx*wx+imxy*wy,
	dy=imyx*wx+imyy*wy;
	x0-=dx;
	y0-=dy;

	if (	fabs(x0-px0)<(fabs(x0)+fabs(px0))*EPS && 
		fabs(y0-py0)<(fabs(y0)+fabs(py0))*EPS )	break;

	px0=x0,
	py0=y0;
  }

/* fprintf(stderr,"transform.c: eval_invert_2d(): iterations=%d\n",i); */
 *rx=x0;
 *ry=y0;

 return(0);
}

