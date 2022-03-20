/*****************************************************************************/
/* firandom-eval.c							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* User-interface functions for firandom: parsing and evaluating of command  */
/* line expressions (for background /-m/ and starlist /-l/		     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <stdarg.h>

#include <psn/psn.h>
#include <fits/fits.h>

#include "fi.h"

#include "io/iof.h"
#include "io/scanarg.h"
#include "io/tokenize.h"
#include "math/fit/lmfit.h"
#include "math/expint/expint.h"
#include "math/point.h"
#include "psn/psn-general.h"

#include "basis.h"
#include "tensor.h"
#include "common.h"
#include "stars.h"

#include "firandom.h"

psnsym defvarnames[]=
 {	{ T_VAR, VAR_X     , "x" },	/* scaled x coordinate	             */
	{ T_VAR, VAR_Y     , "y" },	/* scaled y coordinate	             */
	{ T_VAR, VAR_AX    , "X" },	/* abs. x coordinate ( 0<=X<SX )     */
	{ T_VAR, VAR_AY    , "Y" },	/* abs. y coordinate ( 0<=Y<SY )     */
	{ T_VAR, VAR_I     , "i" },	/* intensity (in ADUs)	     	     */
	{ T_VAR, VAR_M     , "m" },	/* magnitude (see --mag-flux)        */
	{ T_VAR, VAR_SH_F  , "f" },	/* FWHM (in pixels)		     */
	{ T_VAR, VAR_SH_E  , "e" },	/* ellipticity ( 0 <= e < 1 )        */
	{ T_VAR, VAR_SH_P  , "p" },	/* position angle (PA, in degrees)   */
	{ T_VAR, VAR_SN_S  , "s" },	/* Gauss standard deviat.(sigma)     */
	{ T_VAR, VAR_SN_D  , "d" },	/* deviation momentum 1. (delta)     */
	{ T_VAR, VAR_SN_K  , "k" },	/* deviation momentum 2. (kappa)     */
	{ T_VAR, VAR_SI_S  , "S" },	/* S				     */
	{ T_VAR, VAR_SI_D  , "D" },	/* D				     */
	{ T_VAR, VAR_SI_K  , "K" },	/* K				     */
	{ T_VAR, PARAM_CNTR, "n" },	/* seq. number of the star (int.)    */
	{ T_VAR, RND_1     , "t" },	/* random number in [0,1] 1.	     */
	{ T_VAR, RND_2     , "u" },	/* random number in [0,1] 2.	     */
	{ T_VAR, RND_3     , "v" },	/* random number in [0,1] 3.	     */
	{ T_VAR, RND_4     , "w" },	/* random number in [0,1] 4.	     */
	{ T_VAR, PARAM_SX  , "SX" },	/* X Size of the image (integer)     */
	{ T_VAR, PARAM_SY  , "SY" },	/* Y Size of the image (integer)     */
	{ 0, 0, NULL }
 };

psnsym defcoordnames[]=
 {	{ T_VAR, 0, "x" },		/* scaled x coordinate		     */
	{ T_VAR, 1, "y" },		/* scaled y coordinate		     */
	{ T_VAR, 2, "X" },		/* abs. x coordinate ( 0<=X<SX )     */
	{ T_VAR, 3, "Y" },		/* abs. y coordinate ( 0<=Y<SY )     */
	{ 0, 0, NULL }
 };

/*****************************************************************************/

double get_gaussian(double mean,double stddev)
{
 double u,v,r;
 u=sqrt(-2.0*log(drand48()));
 v=2*M_PI*drand48();
 r=u*cos(v);
 return(r*stddev+mean);
}
int get_gaussian_2d(double x0,double y0,double is,double id,double ik,double *rx,double *ry)
{
 double u,v,x,y;
 u=sqrt(-2.0*log(drand48()));
 v=2*M_PI*drand48();
 x=u*cos(v);
 y=u*sin(v);
 *rx=x0+(is+id)*x+( ik  )*y;
 *ry=y0+( ik  )*x+(is-id)*y;
 return(0);
}

/*****************************************************************************/

psn *get_background_psn(char *bgarg)
{
 psnsym	*mysyms[6];
 psn	*in,*bg;

 mysyms[0]=psn_general_op;
 mysyms[1]=psn_general_fn;
 mysyms[2]=psn_general_fn_rnd;
 mysyms[3]=defcoordnames;
 mysyms[4]=NULL;

 in=psn_conv_string(bgarg,mysyms);
 if ( in==NULL )	return(NULL);
 psn_init(in,psn_general_prop);
 bg=psn_conv(in,psn_general_prop);
 psn_free(in);
 if ( bg==NULL )	return(NULL);
 return(bg);
}

int create_background(fitsimage *img,char *bgarg,double stddev,double ox,double oy,double scale,int zoom)
{
 int	sx,sy;
 int	i,j;
 psn	*bg;

 double	vars[16],w,n,izoom,zintm;

 if ( img==NULL )	return(1);
 if ( img->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;

 if ( bgarg==NULL )
  {	for ( i=0 ; i<sy ; i++ )
	 { for ( j=0 ; j<sx ; j++ )
	    {	img->data[i][j]=0.0;		}
	 }
	return(0);
  }
 else
  {	bg=get_background_psn(bgarg);
	if ( bg==NULL )	return(1);

	izoom=1.0/(double)zoom;
	zintm=izoom*izoom;

	for ( i=0 ; i<sy ; i++ )
	 { for ( j=0 ; j<sx ; j++ )
	    {	vars[0]=izoom*((double)j-ox)/scale,
		vars[1]=izoom*((double)i-oy)/scale;
		vars[2]=izoom*(double)j,
		vars[3]=izoom*(double)i;
		psn_double_calc(bg,psn_general_funct,&w,vars);
		if ( stddev>0.0 )	n=get_gaussian(0.0,stddev); /* extra noise */
		else			n=0.0;
		img->data[i][j]=(w+n)*zintm;
	    }
	 }
	return(0);
  }
}

/*****************************************************************************/

int fep_to_sdk(double f,double e,double p,double *s,double *d,double *k)
{
 *s=f/SIG_FWHM;
 *d=(*s)*e/(2-e)*cos(2.0*p*M_PI/180.0);
 *k=(*s)*e/(2-e)*sin(2.0*p*M_PI/180.0);
 return(0);
}
int sdk_to_fep(double s,double d,double k,double *f,double *e,double *p)
{
 double	m;
 *f=s*SIG_FWHM;
 m=sqrt(d*d+k*k);
 *e=1-(s-m)/(s+m);
 *p=0.5*180.0/M_PI*atan2(k,d);
 return(0);
}
int sdk_to_isdk(double s,double d,double k,double *is,double *id,double *ik)
{
 double	idet2,s2,m2;
 s2=s*s;
 m2=d*d+k*k;
 idet2=1.0/((s2-m2)*(s2-m2));
 *is=(s+m2)*idet2;
 *id=-2*s*d*idet2;
 *ik=-2*s*k*idet2;
 return(0);
}
int isdk_to_sdk(double is,double id,double ik,double *s,double *d,double *k)
{
 double	det;
 det=is*is-id*id-ik*ik;
 *s=sqrt(0.5*(is+sqrt(det))/det);
 *d=-id/(2*(*s)*det);
 *k=-ik/(2*(*s)*det);
 return(0);
}

/*****************************************************************************/

int replace_limiters(char *buff)
{
 int	i,in_sl,in_pa,n;
 in_sl=in_pa=0;
 remove_spaces(buff);
 n=strlen(buff);
 for ( i=0 ; i<n ; i++ )
  {	if ( buff[i]=='[' )
	 {	if ( in_sl>=1 )	return(1);
		else		in_sl++;
	 }
	else if ( buff[i]==']' )
	 {	if ( in_sl<=0 )	return(1);
		else		in_sl--;
	 }
	else if ( buff[i]=='(' )
				in_pa++;
	else if ( buff[i]==')' )
	 {	if ( in_pa<=0 )	return(1);
		else		in_pa--;
	 }
	else if ( buff[i]==',' || buff[i]==';' )
	 {	if ( ! in_sl && ! in_pa )	buff[i]='|';
		else if ( in_sl && ! in_pa )	buff[i]=';';
		else				buff[i]=',';
	 }
  }	
 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define		LTYPE_MASK_INT		0x03
#define		LTYPE_INT_MAG		0x01
#define		LTYPE_INT_FLUX		0x02
#define		LTYPE_ABS_X		0x04
#define		LTYPE_ABS_Y		0x08
#define		LTYPE_ABS_XY		(LTYPE_ABS_X|LTYPE_ABS_Y)
#define		LTYPE_MASK_SHAPE	0x30
#define		LTYPE_SHAPE_FEP		0x10
#define		LTYPE_SHAPE_SIG		0x20
#define		LTYPE_SHAPE_SDK		0x30

int create_input_list(char *buff,starlistparam *lp,star **rstars,int *rnstar,int basistype,int sseed)
{
 star	*stars,*ws;
 int	nstar;
 char	**terms;
 int	nterm;
 int	i,j,k,l,n,ltype,last_set;

 psnsym	*mysyms[6];
 psn	*in;

 double	gls,gld,glk,glf,gle,glp,glis,glid,glik,bshx,bshy;
 double	vars[32]; /* x,y,X,Y, i,m, f,e,p, s,d,k, S, D, K, n, t,u,v,w, SX,SY  */
 psn	*exs[16]; /* x,y,X,Y, i,m, f,e,p, s,d,k, S, D, K		     */

 mysyms[0]=psn_general_op;
 mysyms[1]=psn_general_fn;
 mysyms[2]=psn_general_fn_rnd;
 mysyms[3]=defvarnames;
 mysyms[4]=NULL;

 *rstars=stars=NULL;
 *rnstar=nstar=0;

 n=strlen(buff);nterm=1;
 for ( i=0 ; i<n ; i++ )
  {	if ( buff[i]=='|' )	nterm++;	}
 terms=(char **)malloc(sizeof(char *)*(nterm+1));
 tokenize_char(buff,terms,'|',nterm);

 gls=1.0,gld=glk=0.0;
 sdk_to_fep(gls,gld,glk,&glf,&gle,&glp);
 sdk_to_isdk(gls,gld,glk,&glis,&glid,&glik);
 last_set=-1;

 if ( basistype>0 )		bshx=bshy=+0.5;
 else if ( basistype<0 )	bshx=bshy=-0.5;
 else				bshx=bshy=0.0;

 srand48((long int)sseed);

 for ( i=0 ; i<nterm ; i++ )
  { if ( isdigit(terms[i][0]) )
     {	int	nstg,ntag;
	char	**cmd,*buff;

	if ( last_set==0 )
	 {	fep_to_sdk (glf,gle,glp,&gls ,&gld ,&glk );
		sdk_to_isdk(gls,gld,glk,&glis,&glid,&glik);
	 }
	else if ( last_set==1 )
	 {	sdk_to_fep (gls,gld,glk,&glf ,&gle ,&glp );
		sdk_to_isdk(gls,gld,glk,&glis,&glid,&glik);
	 }
	else if ( last_set==2 )
	 {	isdk_to_sdk(glis,glid,glik,&gls,&gld,&glk);
		sdk_to_fep (gls ,gld ,glk ,&glf,&gle,&glp);
	 }

	nstg=0;
	sscanf(terms[i],"%d",&nstg);

	if ( nstg<=0 )	continue;
	for ( j=0 ; j<strlen(terms[i]) ; j++ )
	 {	if ( terms[i][j]=='[' )	break;	}
	if ( j==strlen(terms[i]) )	continue;
	buff=terms[i]+j+1;
	ntag=1;
	for ( j=0 ; j<strlen(buff) ; j++ )
	 {	if ( buff[j]==']' )		buff[j]=0;
		else if ( buff[j]==';' )	ntag++;
	 }

	cmd=(char **)malloc(sizeof(char *)*(ntag+1));
	tokenize_char(buff,cmd,';',ntag);

	for ( l=0 ; l<NUM_SET ; l++ )	exs[l]=NULL;
	for ( j=0 ; j<ntag ; j++ )
	 {	for ( k=0 ; k<NUM_SET ; k++ )
		 {	char	*name;
			name=defvarnames[k].name;
			l=strlen(name);
			if ( memcmp(cmd[j],name,l)==0 && cmd[j][l]=='=' )	break;
		 }
		if ( k==NUM_SET )
		 {	free(cmd);free(terms);
			for ( l=0 ; l<NUM_SET ; l++ )
			 {	if ( exs[l] != NULL )	psn_free(exs[l]);	}
			return(1);
		 }
		in=psn_conv_string(cmd[j]+l+1,mysyms);
		if ( in==NULL )
		 {	free(cmd);free(terms);
			for ( l=0 ; l<NUM_SET ; l++ )
			 {	if ( exs[l] != NULL )	psn_free(exs[l]);	}
			return(1);
		 }
		psn_init(in,psn_general_prop);
		exs[k]=psn_conv(in,psn_general_prop);
		if ( exs[k]==NULL )
		 {	free(cmd);free(terms);
			for ( l=0 ; l<NUM_SET ; l++ )
			 {	if ( exs[l] != NULL )	psn_free(exs[l]);	}
			return(1);
		 }
		psn_free(in);
	 }
	for ( l=0,k=0 ; l<NUM_SET ; l++ )
	 {	if ( exs[l] != NULL )	k=k|(1<<l);		}

	ltype=0;
	if ( (k&MSK_INTENSITY)==MSK_I )			 ltype|=LTYPE_INT_FLUX;
	else if ( (k&MSK_INTENSITY)==MSK_M )		 ltype|=LTYPE_INT_MAG;
	else						 ltype=-1;
	if ( ltype>=0 )
	 {	     if ((k&MSK_COORD)==(MSK_X |MSK_Y )) ltype|=0;
		else if ((k&MSK_COORD)==(MSK_AX|MSK_Y )) ltype|=LTYPE_ABS_X;
		else if ((k&MSK_COORD)==(MSK_X |MSK_AY)) ltype|=LTYPE_ABS_Y;
		else if ((k&MSK_COORD)==(MSK_AX|MSK_AY)) ltype|=LTYPE_ABS_XY;
		else					 ltype=-1;
	 }
	if ( ltype>=0 )
	 {	int m;
		m=0;
		if ( k&MSK_SH )		ltype|=LTYPE_SHAPE_FEP,m++;
		else if ( k&MSK_SN )	ltype|=LTYPE_SHAPE_SIG,m++;
		else if ( k&MSK_SI )	ltype|=LTYPE_SHAPE_SDK,m++;
		if ( m>1 )		ltype=-1;
	 }

	if ( ltype<0 )
	 {	free(cmd);free(terms);
		for ( l=0 ; l<NUM_SET ; l++ )
		 {	if ( exs[l] != NULL )	psn_free(exs[l]);	}
		return(1);
	 }
	for ( j=0 ; j<nstg ; j++ )
	 {	double	lx,ly,li,ls,ld,lk;

		vars[VAR_X]=vars[VAR_Y]=vars[VAR_AX]=vars[VAR_AY]=0.0;
		vars[VAR_I]=vars[VAR_M]=0.0;
		
		vars[VAR_SH_F]=glf ,vars[VAR_SH_E]=gle ,vars[VAR_SH_P]=glp;
		vars[VAR_SN_S]=gls ,vars[VAR_SN_D]=gld ,vars[VAR_SN_K]=glk;
		vars[VAR_SI_S]=glis,vars[VAR_SI_D]=glid,vars[VAR_SI_K]=glik;

		vars[RND_1]=drand48(),
		vars[RND_2]=drand48(),
		vars[RND_3]=drand48(),
		vars[RND_4]=drand48();
		vars[PARAM_SX]=(double)lp->sx,
		vars[PARAM_SY]=(double)lp->sy;
		vars[PARAM_CNTR]=(double)j;
		
		for ( l=0 ; l<NUM_SET ; l++ )
		 {	double	x;
			if ( exs[l]==NULL )	continue;
			psn_double_calc(exs[l],psn_general_funct,&x,vars);
			vars[l]=x;	
		 }
		if ( ltype & LTYPE_ABS_X )	
			lx=vars[VAR_AX],vars[VAR_X]=(lx-lp->ox)/lp->scale;
		else	lx=vars[VAR_AX]=vars[VAR_X]*lp->scale+lp->ox;
		if ( ltype & LTYPE_ABS_Y )
			ly=vars[VAR_AY],vars[VAR_Y]=(ly-lp->oy)/lp->scale;
		else	ly=vars[VAR_AY]=vars[VAR_Y]*lp->scale+lp->oy;

		if ( ! ( ltype & LTYPE_INT_MAG ) )
			li=vars[VAR_I];
		else
			li=mag_to_flux(vars[VAR_M],&lp->mf0);

		if ( (ltype&LTYPE_MASK_SHAPE)==LTYPE_SHAPE_FEP )
			fep_to_sdk (vars[VAR_SH_F],vars[VAR_SH_E],vars[VAR_SH_P],&ls,&ld,&lk);
		else if ( (ltype&LTYPE_MASK_SHAPE)==LTYPE_SHAPE_SDK )
			isdk_to_sdk(vars[VAR_SI_S],vars[VAR_SI_D],vars[VAR_SI_K],&ls,&ld,&lk);
		else
			ls=vars[VAR_SN_S],ld=vars[VAR_SN_D],lk=vars[VAR_SN_K];

		stars=(star *)realloc(stars,sizeof(star)*(nstar+1));

		ws=&stars[nstar];

		ws->location.gcx=lx-bshx,
		ws->location.gcy=ly-bshy;
		ws->flux=li;

		ws->gsig=ls;
		ws->gdel=ld;
		ws->gkap=lk;
		ws->shape.model=SHAPE_ELLIPTIC;

		nstar++;
	 }
	for ( l=0 ; l<NUM_SET ; l++ )
	 {	if ( exs[l] != NULL )	psn_free(exs[l]);	}
	free(cmd);
     }

    else if ( strcmp(terms[i],"reset")==0 )
     {	srand48((long int)sseed);		}

    else 
     {	char wbuff[16];
	for ( l=VAR_SH_F ; l<NUM_SET ; l++ )
	 {	double	x;
		sprintf(wbuff,"%s=%%lg",defvarnames[l].name);
		x=0.0;
		if ( sscanf(terms[i],wbuff,&x)<1 )	continue;
		if ( ! isfinite(x) ) 			continue;
		switch ( l )
		 {	case VAR_SH_F : glf=x ;last_set=0;break;
			case VAR_SH_E : gle=x ;last_set=0;break;
			case VAR_SH_P : glp=x ;last_set=0;break;
			case VAR_SN_S : gls=x ;last_set=1;break;
			case VAR_SN_D : gld=x ;last_set=1;break;
			case VAR_SN_K : glk=x ;last_set=1;break;
			case VAR_SI_S : glis=x;last_set=2;break;
			case VAR_SI_D : glid=x;last_set=2;break;
			case VAR_SI_K : glik=x;last_set=2;break;
		 }
		break;
	 }
	if ( l==NUM_SET )
	 {	free(terms);
		return(1);
	 }
     }
  }
 free(terms);

 *rstars=stars;
 *rnstar=nstar;
 return(0);
}

