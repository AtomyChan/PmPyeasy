/*****************************************************************************/
/* optcalc.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (C) 2016, Pal, Andras <apal@szofi.net>				     */
/*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

#include "math/matrixvector.h"
#include "io/tokenize.h"
#include "optcalc.h"

int optcalc_raytrace_reset(struct raytrace *rt)
{
 rt->rt_npoint=0;
 rt->rt_points=NULL;
 return(0);
}

int optcalc_raytrace_free(struct raytrace *rt)
{
 if ( rt->rt_points != NULL )	free(rt->rt_points);
 optcalc_raytrace_reset(rt);
 return(0);
}

static int optcalc_raytrace_append(struct raytrace *rt,vector v)
{
 if ( rt==NULL )
	return(0);
 rt->rt_points=realloc(rt->rt_points,sizeof(vector)*(rt->rt_npoint+1));
 rt->rt_points[rt->rt_npoint][0]=v[0];
 rt->rt_points[rt->rt_npoint][1]=v[1];
 rt->rt_points[rt->rt_npoint][2]=v[2];
 rt->rt_npoint++;
 return(0);
}

double optcalc_get_refraction_index(struct glass *g,double lambda)
{
 if ( 0<g->g_nncoeff )
  {	double	x2,sum,n;
	int	i;
	x2=1.0/(lambda*lambda);
	sum=0.0;
	for ( i=0 ; i<g->g_nncoeff ; i++ )
	 {	sum+=g->g_ncoeffs[2*i+0]/(1-g->g_ncoeffs[2*i+1]*x2);	}
	n=sqrt(1+sum);
	return(n);
  }
 else if ( 0.0<g->g_n )
	return(g->g_n);
 else
	return(1.0);
}

int optcalc_glass_refraction_precompute(struct optics *opt,double lambda)
{
 int	i;
 for ( i=0 ; i<opt->opt_nlens ; i++ )
  {	struct	lens	*l;
	l=&opt->opt_lenses[i];
	if ( 0.0<lambda )
	 {	double	n_refr;
		struct	glass	*g;
		if ( l->l_index_glass<0 )
			n_refr=1;
		else
		 {	g=&opt->opt_glasses[l->l_index_glass];
			n_refr=optcalc_get_refraction_index(g,lambda);
		 }
		l->l_n_lambda=n_refr;
	 }
	else
		l->l_n_lambda=0.0;
  }
 return(0);
}

int optcalc_snell_descartes(vector norm,vector s1,vector s2,double n1,double n2)
{
 double	s1norm,d,nr;
 vector	lnorm;

 lnorm[0]=norm[0];
 lnorm[1]=norm[1];
 lnorm[2]=norm[2];
 s1norm=vector_vector_mul(s1,lnorm);
 if ( s1norm<0.0 )
  {	s1norm=-s1norm;
	lnorm[0]=-lnorm[0];
	lnorm[1]=-lnorm[1];
	lnorm[2]=-lnorm[2];
  }
 nr=n1/n2;
 d=1-nr*nr*(1-s1norm*s1norm);
 /* this is total internal reflection: */
 if ( d<0.0 )	
  {	s2[0]=s1[0]-2*lnorm[0]*s1norm;
	s2[1]=s1[1]-2*lnorm[1]*s1norm;
	s2[2]=s1[2]-2*lnorm[2]*s1norm;
	return(1);
  }
 /* this is normal refraction: */
 else
  {	d=sqrt(d);
	s2[0]=nr*(s1[0]-lnorm[0]*s1norm)+lnorm[0]*d;
	s2[1]=nr*(s1[1]-lnorm[1]*s1norm)+lnorm[1]*d;
	s2[2]=nr*(s1[2]-lnorm[2]*s1norm)+lnorm[2]*d;
	return(0);
  }
}

int optcalc_reflection(vector norm,vector s1,vector s2)
{
 double	s1norm;
 vector	lnorm;

 lnorm[0]=norm[0];
 lnorm[1]=norm[1];
 lnorm[2]=norm[2];
 s1norm=vector_vector_mul(s1,lnorm);
 if ( s1norm<0.0 )
  {	s1norm=-s1norm;
	lnorm[0]=-lnorm[0];
	lnorm[1]=-lnorm[1];
	lnorm[2]=-lnorm[2];
  }

 s2[0]=s1[0]-2*lnorm[0]*s1norm;
 s2[1]=s1[1]-2*lnorm[1]*s1norm;
 s2[2]=s1[2]-2*lnorm[2]*s1norm;
 return(0);
}

int optcalc_surface_aspheric_eval(double c,double k,int nalpha,double *alphas,double r2,double *rz)
{
 double	z,d,rp;
 int	i;
 
 d=1-(1+k)*c*c*r2;
 if ( d<0.0 )	return(1);
 z=c*r2/(1+sqrt(d));
 for ( i=0,rp=r2 ; i<+nalpha && alphas != NULL ; i++ )
  {	z+=alphas[i]*rp;
	rp=rp*r2;
  }
 for ( i=0,rp=r2 ; i<-nalpha && alphas != NULL ; i++ )
  {	z-=alphas[i]*rp;
	rp=rp*r2;
  }
 if ( rz != NULL )	*rz=z;
 return(0);
}

int optcalc_surface_aspheric_diff(double c,double k,int nalpha,double *alphas,double r2,double *rdz)
{
 double	dz,d,rp;
 int	i;
 
 d=1-(1+k)*c*c*r2;
 if ( d<0.0 )	return(1);
 d=sqrt(d);
 dz=c/(1+d)*(1+c*c*r2*(1+k)/(2*d*(1+d)));
 for ( i=0,rp=1.0 ; i<+nalpha && alphas != NULL ; i++ )
  {	dz+=(double)(i+1)*alphas[i]*rp;
	rp=rp*r2;
  }
 for ( i=0,rp=1.0 ; i<-nalpha && alphas != NULL ; i++ )
  {	dz-=(double)(i+1)*alphas[i]*rp;
	rp=rp*r2;
  }
 if ( rdz != NULL )	*rdz=dz;
 return(0);
}

/* optcalc_surface_aspheric_ray_trace():
   computes the intersection of the aspheric surface described by the 
   parameter set of (c,k,alphas[nalpha]) and the ray going from x0[] to the
   direction n0[]. 0 is returned if there is an intersection and x_int[]
   and n_int[] is set appropriaterly (intersection point and surface normal
   vector, respectively). Here the aspheric surface is directed to tangent the 
   (x,y) plane at (0,0,0) while positive curvature referes to positive z 
   values.
*/
int optcalc_surface_aspheric_ray_trace(double c,double k,int nalpha,double *alphas,vector x0,vector n0,vector x_int,vector n_int)
{
 double	t,dz,r2;
 vector	x,n;
 int	i;

 x[0]=x0[0];
 x[1]=x0[1];
 x[2]=x0[2];
 n[0]=n0[0];
 n[1]=n0[1];
 n[2]=n0[2];

 t=-x[2]/n[2];
 x[0]+=t*n[0];
 x[1]+=t*n[1];
 x[2]+=t*n[2];

 t=0.0;
 for ( i=0 ; i<15 ; i++ )
  {	double	xi,yi,zi,z,dz,r2,tnext;
	xi=x[0]+t*n[0];
	yi=x[1]+t*n[1];
	zi=x[2]+t*n[2];
	r2=xi*xi+yi*yi;
	
	if ( optcalc_surface_aspheric_eval(c,k,nalpha,alphas,r2,&z) )
		return(1);
	else if ( optcalc_surface_aspheric_diff(c,k,nalpha,alphas,r2,&dz) )
		return(1);
	
	tnext=t-(z-zi)/(dz*(2*(xi*n[0]+yi*n[1])+2*(n[0]*n[0]+n[1]*n[1])*t)-n[2]);
	if ( t==tnext )	break;
	t=tnext;
  }

 x_int[0]=x[0]+t*n[0];
 x_int[1]=x[1]+t*n[1];
 x_int[2]=x[2]+t*n[2];
 r2=x_int[0]*x_int[0]+x_int[1]*x_int[1];
 optcalc_surface_aspheric_diff(c,k,nalpha,alphas,r2,&dz);
 n_int[0]=2*x_int[0]*dz;
 n_int[1]=2*x_int[1]*dz;
 n_int[2]=-1;
 vector_normalize(n_int);

 return(0);
}

int optcalc_refract_quadratic_ray_trace(matrix A,vector B,double C,vector x0,vector n0,double rn1,double rn2)
{
 double	qa,qb,qc;
 double	t;
 vector	norm;

 qa=0.5*vector_matrix_vector_mul(n0,A,n0);
 qb=vector_matrix_vector_mul(x0,A,n0)+vector_vector_mul(B,n0);
 qc=0.5*vector_matrix_vector_mul(x0,A,x0)+vector_vector_mul(B,x0)+C;
 if ( qa != 0.0 )
  {	double	t1,t2,dc;
	dc=qb*qb-4*qa*qc;
	if ( dc<=0.0 )	return(1);
	dc=sqrt(dc);
	t1=(-qb-dc)/(2*qa);
	t2=(-qb+dc)/(2*qa);
	if ( t2<t1 )	t=t1,t1=t2,t2=t;
	if ( t2<0.0 )	return(1);
	else if ( t1<0 )	t=t2;
	else			t=t1;
  }
 else
	t=-qc/qb;

 x0[0]=x0[0]+t*n0[0];
 x0[1]=x0[1]+t*n0[1];
 x0[2]=x0[2]+t*n0[2];
 vector_mul(norm,A,x0);
 norm[0]+=B[0];
 norm[1]+=B[1];
 norm[2]+=B[2];
 vector_normalize(norm);
 optcalc_snell_descartes(norm,n0,n0,rn1,rn2);

 return(0);
}

int optcalc_refract_aspheric_ray_trace(double c,double k,int nalpha,double *alphas,double z0,vector x0,vector n0,double rn1,double rn2)
{
 vector	norm;

 x0[2]-=z0;
 if ( optcalc_surface_aspheric_ray_trace(c,k,nalpha,alphas,x0,n0,x0,norm) )
	return(1);
 x0[2]+=z0;
 optcalc_snell_descartes(norm,n0,n0,rn1,rn2); 

 return(0);
}

int optcalc_surface_reset(struct surface *s)
{
 s->s_curvature=0.0;
 s->s_conic=0.0;
 s->s_nalpha=0;
 return(0);
}
int optcalc_surface_is_quadratic(struct surface *s)
{
 if ( s->s_conic != 0.0 || 0<s->s_nalpha )
	return(0);
 else
	return(1);
}

int optcalc_lens_ray_trace(struct lens *l,double n_refr,vector x0,vector n0,struct raytrace *rt)
{
 matrix	A;
 vector	B,xnext,nnext;
 double	C,cc,z0;

 xnext[0]=x0[0];
 xnext[1]=x0[1];
 xnext[2]=x0[2];
 nnext[0]=n0[0];
 nnext[1]=n0[1];
 nnext[2]=n0[2];

 optcalc_raytrace_append(rt,xnext);
/* fprintf(stdout,"%12g %12g %12g\n",xnext[0],xnext[1],xnext[2]); */

 /* first surface: */
 if ( optcalc_surface_is_quadratic(&l->l_s1) )
  {	cc=+l->l_s1.s_curvature;
	z0=l->l_offset-l->l_thickness/2.0;
	A[0][0]=A[1][1]=A[2][2]=cc;
	A[0][1]=A[1][2]=A[2][0]=0.0;
	A[1][0]=A[2][1]=A[0][2]=0.0;
	B[0]=0;
	B[1]=0;
	B[2]=-1-z0*cc;
	C=z0*(1+0.5*z0*cc);
	if ( optcalc_refract_quadratic_ray_trace(A,B,C,xnext,nnext,1.0,n_refr) )
		return(1);
  }
 else
  {	double	z0;
	z0=l->l_offset-l->l_thickness/2.0;
	if ( optcalc_refract_aspheric_ray_trace(+l->l_s1.s_curvature,l->l_s1.s_conic,+l->l_s1.s_nalpha,l->l_s1.s_alphas,z0,xnext,nnext,1.0,n_refr) )
		return(1);
  }

 if ( l->l_radius1*l->l_radius1 <= xnext[0]*xnext[0]+xnext[1]*xnext[1]  )
	return(1);

 optcalc_raytrace_append(rt,xnext);
/* fprintf(stdout,"%12g %12g %12g\n",xnext[0],xnext[1],xnext[2]); */

 /* second surface: */
 if ( optcalc_surface_is_quadratic(&l->l_s2) )
  {	cc=-l->l_s2.s_curvature;
	z0=l->l_offset+l->l_thickness/2.0;
	A[0][0]=A[1][1]=A[2][2]=cc;
	A[0][1]=A[1][2]=A[2][0]=0.0;
	A[1][0]=A[2][1]=A[0][2]=0.0;
	B[0]=0;
	B[1]=0;
	B[2]=-1-z0*cc;
	C=z0*(1+0.5*z0*cc);
	if ( optcalc_refract_quadratic_ray_trace(A,B,C,xnext,nnext,n_refr,1.0) )
		return(1);
  }
 else
  {	double	z0;
	z0=l->l_offset+l->l_thickness/2.0;
	if ( optcalc_refract_aspheric_ray_trace(-l->l_s2.s_curvature,l->l_s2.s_conic,-l->l_s2.s_nalpha,l->l_s2.s_alphas,z0,xnext,nnext,n_refr,1.0) )
		return(1);
  }

 if ( l->l_radius2*l->l_radius2 <= xnext[0]*xnext[0]+xnext[1]*xnext[1] )
	return(1);

/* fprintf(stdout,"%12g %12g %12g\n",xnext[0],xnext[1],xnext[2]); */

 x0[0]=xnext[0];
 x0[1]=xnext[1]; 
 x0[2]=xnext[2]; 
 n0[0]=nnext[0];
 n0[1]=nnext[1]; 
 n0[2]=nnext[2]; 

 return(0);
}

int optcalc_ray_trace(struct optics *opt,double lambda,vector v0,vector n0,struct raytrace *rt)
{
 int	i;
 double	t;

 vector_normalize(n0);

 for ( i=0 ; i<opt->opt_nlens ; i++ )
  {	double	n_refr;
	struct	lens	*l;
	l=&opt->opt_lenses[i];
	if ( 0<l->l_n_lambda )
		n_refr=l->l_n_lambda;
	else
	 {	struct	glass	*g;
		if ( l->l_index_glass<0 )
			n_refr=1;
		else
		 {	g=&opt->opt_glasses[l->l_index_glass];
			n_refr=optcalc_get_refraction_index(g,lambda);
		 }
	 }
	if ( optcalc_lens_ray_trace(l,n_refr,v0,n0,rt) )
		return(1);
  }

 optcalc_raytrace_append(rt,v0);

 t=(opt->opt_z_focal-v0[2])/n0[2];
 v0[0]+=t*n0[0];
 v0[1]+=t*n0[1];
 v0[2]+=t*n0[2];

 optcalc_raytrace_append(rt,v0);

 return(0);
}

int optcalc_transfer_matrix_compose(transfer_matrix m,transfer_matrix next)
{
 transfer_matrix o;

 o[0][0]=next[0][0]*m[0][0]+next[0][1]*m[1][0];
 o[0][1]=next[0][0]*m[0][1]+next[0][1]*m[1][1];
 o[1][0]=next[1][0]*m[0][0]+next[1][1]*m[1][0];
 o[1][1]=next[1][0]*m[0][1]+next[1][1]*m[1][1];

 m[0][0]=o[0][0];
 m[0][1]=o[0][1];
 m[1][0]=o[1][0];
 m[1][1]=o[1][1];

 return(0); 
}

int optcalc_compute_transfer_matrix(struct optics *opt,double lambda,transfer_matrix m,double z0,double *rzend)
{
 int	i;
 double	distance;

 m[0][0]=m[1][1]=1.0;
 m[0][1]=m[1][0]=0.0;

 if ( 1 )
  {	struct	lens	*l;
	transfer_matrix	n;

	l=&opt->opt_lenses[0];
 	distance=l->l_offset-l->l_thickness/2.0-z0;
	n[0][0]=1.0;
	n[0][1]=distance;
	n[1][0]=0.0;
	n[1][1]=1.0;
	optcalc_transfer_matrix_compose(m,n);
  }
	
 for ( i=0 ; i<opt->opt_nlens ; i++ )
  {	double	n_refr;
	struct	glass	*g;
	struct	lens	*l;
	transfer_matrix	n;

	l=&opt->opt_lenses[i];

	if ( l->l_index_glass<0 )
		n_refr=1;
	else
	 {	g=&opt->opt_glasses[l->l_index_glass];
		n_refr=optcalc_get_refraction_index(g,lambda);
	 }

	n[0][0]=1.0;
	n[0][1]=0.0;
	n[1][0]=+l->l_s1.s_curvature*(1.0-n_refr)/(n_refr);
	n[1][1]=1.0/n_refr;
	optcalc_transfer_matrix_compose(m,n);
	n[0][0]=1.0;
	n[0][1]=l->l_thickness;
	n[1][0]=0.0;
	n[1][1]=1.0;
	optcalc_transfer_matrix_compose(m,n);
	n[0][0]=1.0;
	n[0][1]=0.0;
	n[1][0]=-l->l_s2.s_curvature*(n_refr-1.0)/1.0;
	n[1][1]=n_refr/1.0;
	optcalc_transfer_matrix_compose(m,n);

	if ( i<opt->opt_nlens-1 )
	 {	struct	lens	*lnext;
		lnext=&opt->opt_lenses[i+1];
		distance=(lnext->l_offset-lnext->l_thickness/2.0)-(l->l_offset+l->l_thickness/2.0);
		n[0][0]=1.0;
		n[0][1]=distance;
		n[1][0]=0.0;
		n[1][1]=1.0;
		optcalc_transfer_matrix_compose(m,n);
	 }
	if ( rzend != NULL )
	 	*rzend=l->l_offset+l->l_thickness/2.0;
  }

 return(0);
}

int optcalc_reset(struct optics *opt)
{
 opt->opt_glasses=NULL;
 opt->opt_nglass=0;
 opt->opt_lenses=NULL;
 opt->opt_nlens=0;
 opt->opt_z_focal=0.0;
 optcalc_surface_reset(&opt->opt_s_focal);
 return(0);
}

int optcalc_read_curvature(char *str,double *ret)
{
 double	w;
 if ( sscanf(str,"+1/%lg",&w)==1 )
  {	*ret=+1.0/w;
	return(0);
  }
 else if ( sscanf(str,"-1/%lg",&w)==1 )
  {	*ret=-1.0/w;
	return(0);
  }
 else if ( sscanf(str,"%lg",&w)==1 )
  {	*ret=w;
	return(0);
  }
 else
	return(1);
}

int optcalc_read_aspheric_data(char *str,struct surface *s)
{
 char	*tokens[16];
 int	i,n;

 if ( strcmp(str,"-")==0 )
  {	s->s_nalpha=0;
	return(0);
  }

 n=tokenize_char(str,tokens,':',16);
 for ( i=0 ; i<n ; i++ )
  {	if ( sscanf(tokens[i],"%lg",&s->s_alphas[i])<1 )
		return(1);
  }
 s->s_nalpha=n;
 return(0);
}

int optcalc_read_radii(char *str,double *rr1,double *rr2)
{
 int	n;
 double	r1,r2;
 n=sscanf(str,"%lg:%lg",&r1,&r2);
 if ( n==2 )
  {	*rr1=r1;
	*rr2=r2;
	return(0);
  }
 else if ( n==1 )
  {	*rr1=*rr2=r1;
	return(0);
  }
 else
	return(1);
}


int optcalc_read_optics(FILE *fr,struct optics *opt)
{
 char	buff[4096],**cmd;
 int	i,n;
 int	nline;

 nline=0;
 while ( ! feof(fr) )
  {	if ( fgets(buff,4096,fr)==NULL )
		break;
	nline++;
	buff[4095]=0;
	for ( i=0 ; buff[i] ; i++ )
	 {	if ( buff[i]=='#' )
		 {	buff[i]=0;
			break;
		 }
	 }
	cmd=tokenize_spaces_dyn(buff);
	for ( n=0 ; cmd != NULL && cmd[n] != NULL ; )	n++;
	if ( cmd==NULL || n<=0 )
	 {	free(cmd);
		continue;
	 }
	if ( strcmp(cmd[0],"glass")==0 && n==3 )
	 {	struct	glass	*g;
		char		*ncs[16];
		int		i,n;
		opt->opt_glasses=realloc(opt->opt_glasses,sizeof(struct glass)*(opt->opt_nglass+1));
		if ( 31<strlen(cmd[1]) )
			return(nline);
		g=&opt->opt_glasses[opt->opt_nglass];
		strcpy(g->g_name,cmd[1]);
		n=tokenize_char(cmd[2],ncs,':',16);
		if ( n==1 )
		 {	if ( sscanf(ncs[0],"%lg",&g->g_n)<1 )
				return(nline);
			g->g_nncoeff=0;
		 }
		else if ( n%2==0 && n<=8 )
		 {	for ( i=0 ; i<n ; i++ )
			 {	if ( sscanf(ncs[i],"%lg",&g->g_ncoeffs[i])<1 )
					return(nline);
			 }
			g->g_n=0.0;
			g->g_nncoeff=n/2;
		 }
		else
			return(nline);

		opt->opt_nglass++;
	 }
	else if ( strcmp(cmd[0],"lens")==0 && 7<=n )
	 {	struct	lens	*l;
		int		i;
		opt->opt_lenses=realloc(opt->opt_lenses,sizeof(struct lens)*(opt->opt_nlens+1));
		l=&opt->opt_lenses[opt->opt_nlens];
		optcalc_surface_reset(&l->l_s1);
		optcalc_surface_reset(&l->l_s2);
		if ( sscanf(cmd[1],"%lg",&l->l_offset)<1 )
			return(nline);
		else if ( sscanf(cmd[2],"%lg",&l->l_thickness)<1 )
			return(nline);
		else if ( optcalc_read_radii(cmd[3],&l->l_radius1,&l->l_radius2) )
			return(nline);
		else if ( optcalc_read_curvature(cmd[4],&l->l_s1.s_curvature) )
			return(nline);
		else if ( optcalc_read_curvature(cmd[5],&l->l_s2.s_curvature) )
			return(nline);
		for ( i=0 ; i<opt->opt_nglass ; i++ )
		 {	if ( strcmp(cmd[6],opt->opt_glasses[i].g_name)==0 )
				break;
		 }
		if ( i<opt->opt_nglass )
			l->l_index_glass=i;
		else if ( strcmp(cmd[6],"-")==0 )
			l->l_index_glass=-1;
		else
			return(nline);

		if ( 9<=n )
		 {	if ( optcalc_read_aspheric_data(cmd[7],&l->l_s1) )
				return(nline);
			if ( optcalc_read_aspheric_data(cmd[8],&l->l_s2) )
				return(nline);
		 }
		l->l_n_lambda=0.0;

		opt->opt_nlens++;
	 }
	else if ( strcmp(cmd[0],"focal")==0 && 2<=n )
	 {	if ( sscanf(cmd[1],"%lg",&opt->opt_z_focal)<1 )
			return(nline);
	 }
	else
		return(nline);
  }
 return(0);
}

