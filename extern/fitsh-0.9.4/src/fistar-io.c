/*****************************************************************************/
/* fistar-io.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* I/O module for fistar						     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>

#include <fits/fits.h>

#include "fitsh.h"

#include "fitsmask.h"
#include "math/spline/biquad.h"
#include "math/spline/biquad-isc.h"
#include "math/fit/lmfit.h"
#include "math/poly.h"
#include "statistics.h"
#include "io/iof.h"
#include "io/scanarg.h"
#include "io/tokenize.h"

#include "magnitude.h"
#include "imgtrans.h"
#include "tensor.h"
#include "common.h"

#include "stars.h"

#include "fistar.h"

/*****************************************************************************/

static char *wrninappr="Warning: inappropriate content in line %d, skipped.\n";

int read_star_candidates(FILE *fr,colinfo *col,initcand **ricands,int *rnicand,magflux *mf)
{
 initcand	*icands,*wic;
 int		nicand,n,ln;
 char		*rbuff,**cmd;
 double		rx,ry,rs,rd,rk,rf,rb;

 icands=NULL;
 nicand=0;

 rbuff=NULL;cmd=NULL;ln=0;

 while ( ! feof(fr) ) 
  {	if ( rbuff != NULL )
	 {	free(rbuff);rbuff=NULL;		}
	if ( cmd != NULL )
	 {	free(cmd);  cmd=NULL;		}
	rbuff=freadline(fr);
	if ( rbuff==NULL )	break;
	ln++;
	remove_newlines_and_comments(rbuff);
	if ( rbuff[0]==0 )	break;
	cmd=tokenize_spaces_dyn(rbuff);
	if ( cmd==NULL )	break;
	for ( n=0 ; cmd[n] != NULL ; )	n++;
	if ( col->x>=n || col->y>=n )
	 {	fprintf(stderr,wrninappr,ln);continue;		}
	if ( sscanf(cmd[col->x],"%lg",&rx)<1 || sscanf(cmd[col->y],"%lg",&ry)<1 )
	 {	fprintf(stderr,wrninappr,ln);continue;		}
	if ( ! isfinite(rx) || ! isfinite(ry) )	
	 {	fprintf(stderr,wrninappr,ln);continue;		}

	rs=rd=rk=0.0;
	rf=0.0;
	rb=0.0;

	if ( 0<=col->s )
	 {	if ( col->s>=n || sscanf(cmd[col->s],"%lg",&rs)<1 )
		 {	fprintf(stderr,wrninappr,ln);continue;		}
	 }
	if ( 0<=col->d && 0<=col->k )
	 {	if ( col->d>=n || col->k>=n )
		 {	fprintf(stderr,wrninappr,ln);continue;		}
		if ( sscanf(cmd[col->d],"%lg",&rd)<1 || sscanf(cmd[col->k],"%lg",&rk)<1 )
		 {	fprintf(stderr,wrninappr,ln);continue;		}
	 }
	if ( ! isfinite(rs) || ! isfinite(rd) || ! isfinite(rk) )
	 {	fprintf(stderr,wrninappr,ln);continue;		}

	if ( 0<=col->flux )
	 {	if ( col->flux>=n || sscanf(cmd[col->flux],"%lg",&rf)<1 )
		 {	fprintf(stderr,wrninappr,ln);continue;		}
	 }
	else if ( 0<=col->mag ) 
	 {	if ( col->mag>=n || sscanf(cmd[col->mag],"%lg",&rf)<1 )
		 {	fprintf(stderr,wrninappr,ln);continue;		}
		rf=mag_to_flux(rf,mf);
	 }
	if ( 0<=col->bg )
	 {	if ( col->bg>=n || sscanf(cmd[col->bg],"%lg",&rb)<1 )
		 {	fprintf(stderr,wrninappr,ln);continue;		}
	 }

	icands=(initcand *)realloc(icands,sizeof(initcand)*(nicand+1));
	wic=&icands[nicand];
	wic->x=rx,wic->y=ry;
	wic->s=rs,
	wic->d=rd,
	wic->k=rk;
	wic->bg=rb;
	wic->flux=rf;
	nicand++;
  };
 if ( rbuff != NULL )	free(rbuff);
 if ( cmd != NULL )	free(cmd);

 if ( ricands != NULL )	*ricands=icands;
 if ( rnicand != NULL )	*rnicand=nicand;

 return(0);
}

/*****************************************************************************/

int read_input_position_list(FILE *fr,colinfo *col,iposition **rips,int *rnip)
{
 iposition	*ips;
 int		nip,ln,n;
 char		*rbuff,**cmd;
 double		dx,dy;

 rbuff=NULL;
 cmd=NULL;

 ips=NULL;
 nip=0;

 ln=0;

 while ( ! feof(fr) ) 
  {	if ( rbuff != NULL )
	 {	free(rbuff);rbuff=NULL;		}
	if ( cmd != NULL )
	 {	free(cmd);  cmd=NULL;		}

	rbuff=freadline(fr);
	if ( rbuff==NULL )	break;
	ln++;
	remove_newlines_and_comments(rbuff);
	if ( rbuff[0]==0 )	
		continue;
	cmd=tokenize_spaces_dyn(rbuff);
	if ( cmd==NULL )
		continue;
	for ( n=0 ; cmd[n] != NULL ; )	n++;
	if ( col->x>=n || col->y>=n )
		continue;
	if ( sscanf(cmd[col->x],"%lg",&dx)<1 || sscanf(cmd[col->y],"%lg",&dy)<1 )
		continue;

	ips=(iposition *)realloc(ips,sizeof(iposition)*(nip+1));
	if ( col->id<n )
		ips[nip].name=strdup(cmd[col->id]);
	else
		ips[nip].name=NULL;
	ips[nip].x=dx;
	ips[nip].y=dy;

	nip++;

  }

 if ( rbuff != NULL )	free(rbuff);
 if ( cmd != NULL )	free(cmd);

 if ( rips != NULL )	*rips=ips;
 if ( rnip != NULL )	*rnip=nip;

 return(0);
}

/*****************************************************************************/

int draw_mark_star_pixels(fitsimage *img,star *stars,int nstar)
{
 int		n,i,ix,iy;
 candidate	*wc;

 for ( n=0 ; n<nstar ; n++ )
  {	wc=stars[n].cand;
	if ( wc==NULL )	continue;
	if ( wc->nipoint==0 || wc->ipoints==NULL )	continue;
	for ( i=0 ; i<wc->nipoint ; i++ )
	 {	ix=wc->ipoints[i].x,
		iy=wc->ipoints[i].y;
		img->data[iy][ix]=0.0;
	 }
  }
 return(0);
}

int draw_mark_stars(fitsimage *img,star *stars,int nstar,int msym,int msize)
{
 int	i,ix,iy;
 double	color;
 color=0.0;
 for ( i=0 ; i<nstar ; i++ )
  {	if ( stars[i].cand==NULL )	continue;
	ix=stars[i].cand->ix,iy=stars[i].cand->iy;
	switch ( msym )
	 {   case MARK_SYM_DOTS:
		fits_image_draw_pixel(img,ix,iy,color);
		break;
	     case MARK_SYM_SQUARES:
		fits_image_draw_line(img,ix-msize,iy-msize,ix+msize,iy-msize,color,-1);
		fits_image_draw_line(img,ix+msize,iy-msize,ix+msize,iy+msize,color,-1);
		fits_image_draw_line(img,ix+msize,iy+msize,ix-msize,iy+msize,color,-1);
		fits_image_draw_line(img,ix-msize,iy+msize,ix-msize,iy-msize,color,-1);
		break;
	     case MARK_SYM_CIRCLES:
		fits_image_draw_circle(img,ix,iy,msize,color);
		break;
	 }
  }
 return(0);
}

/*****************************************************************************/

typedef struct
 {	int	n;	
	magflux	*mf;
 } formatoutparam;

typedef struct
 {	int	id;
	int	width;
	char	*comment;
	int	(*fprint)(FILE *fw,star *s,formatoutparam *fop);
	char	*tlist[7];
 } formatname;

static int fprint_star_id(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%6d ",fop->n+1);
 return(0);
}
static int fprint_star_ix(FILE *fw,star *s,formatoutparam *fop)
{
 if ( s->cand==NULL )	fprintf(fw,"%4s ","-");
 else			fprintf(fw,"%4d ",s->cand->ix+1);
 return(0);
}
static int fprint_star_iy(FILE *fw,star *s,formatoutparam *fop)
{
 if ( s->cand==NULL )	fprintf(fw,"%4s ","-");
 else			fprintf(fw,"%4d ",s->cand->iy+1);
 return(0);
}
static int fprint_star_cx(FILE *fw,star *s,formatoutparam *fop)
{
 if ( s->cand==NULL )	fprintf(fw,"%8s ","-");
 else			fprintf(fw,"%8.3f ",s->cand->cx);
 return(0);
}
static int fprint_star_cy(FILE *fw,star *s,formatoutparam *fop)
{
 if ( s->cand==NULL )	fprintf(fw,"%8s ","-");
 else			fprintf(fw,"%8.3f ",s->cand->cy);
 return(0);
}
static int fprint_star_cbg(FILE *fw,star *s,formatoutparam *fop)
{
 if ( s->cand==NULL )	fprintf(fw,"%8s ","-");
 else			fprintf(fw,"%8.2f ",s->cand->bg);
 return(0);
}
static int fprint_star_camp(FILE *fw,star *s,formatoutparam *fop)
{
 if ( s->cand==NULL )	fprintf(fw,"%8s ","-");
 else			fprintf(fw,"%8.2f ",s->cand->amp);
 return(0);
}
static int fprint_star_cmax(FILE *fw,star *s,formatoutparam *fop)
{
 if ( s->cand==NULL )	fprintf(fw,"%8s ","-");
 else			fprintf(fw,"%8.2f ",s->cand->peak);
 return(0);
}
static int fprint_star_cdevs(FILE *fw,star *s,formatoutparam *fop)
{
 if ( s->cand==NULL )	fprintf(fw,"%6s ","-");
 else			fprintf(fw,"%6.3f ",0.5*(s->cand->sxx+s->cand->syy));
 return(0);
}
static int fprint_star_cdevd(FILE *fw,star *s,formatoutparam *fop)
{
 if ( s->cand==NULL )	fprintf(fw,"%6s ","-");
 else			fprintf(fw,"%6.3f ",0.5*(s->cand->sxx-s->cand->syy));
 return(0);
}
static int fprint_star_cdevk(FILE *fw,star *s,formatoutparam *fop)
{
 if ( s->cand==NULL )	fprintf(fw,"%6s ","-");
 else			fprintf(fw,"%6.3f ",s->cand->sxy);
 return(0);
}
static int fprint_star_x(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%8.3f ",s->location.gcx);
 return(0);
}
static int fprint_star_y(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%8.3f ",s->location.gcy);
 return(0);
}
static int fprint_star_bg(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%8.2f ",s->location.gbg);
 return(0);
}
static int fprint_star_amp(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%8.2f ",s->location.gamp);
 return(0);
}
static int fprint_star_dsig(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%6.3f ",s->gsig);
 return(0);
}
static int fprint_star_ddel(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%6.3f ",s->gdel);
 return(0);
}
static int fprint_star_dkap(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%6.3f ",s->gkap);
 return(0);
}
static int fprint_star_fwhm(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%5.3f ",s->gfwhm);
 return(0);
}
static int fprint_star_ellip(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%5.3f ",s->gellip);
 return(0);
}
static int fprint_star_pa(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%5.1f ",s->gpa);
 return(0);
}
static int fprint_star_devs(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%6.3f ",s->shape.gs);
 return(0);
}
static int fprint_star_devd(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%6.3f ",s->shape.gd);
 return(0);
}
static int fprint_star_devk(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%6.3f ",s->shape.gk);
 return(0);
}
static int fprint_star_devl(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%6.3f ",s->shape.gl);
 return(0);
}

static int fprint_star_flux(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%10.2f ",s->flux);
 return(0);
}
static int fprint_star_mag(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%9.3f ",flux_to_mag(s->flux,fop->mf));
 return(0);
}
static int fprint_star_noise(FILE *fw,star *s,formatoutparam *fop)
{
 if ( s->cand==NULL )	fprintf(fw,"%8s ","-");
 else			fprintf(fw,"%8.2f ",s->cand->noise);
 return(0);
}
static int fprint_star_sn(FILE *fw,star *s,formatoutparam *fop)
{
 if ( s->cand==NULL )		fprintf(fw,"%8s ","-");
 else if ( s->cand->noise>0.0 )	fprintf(fw,"%8.1f ",s->flux/s->cand->noise);
 else				fprintf(fw,"%8s ","-");
 return(0);
}

static int fprint_star_cpix(FILE *fw,star *s,formatoutparam *fop)
{
 if ( s->cand==NULL )	fprintf(fw,"%5s ","-");
 else 			fprintf(fw,"%5d ",s->cand->nipoint);
 return(0);
}
static int fprint_star_mom(FILE *fw,star *s,formatoutparam *fop)
{
 int	i;
 double	m;
 for ( i=3 ; i<(s->shape.order+1)*(s->shape.order+2)/2 ; i++ )
  {	m=s->shape.mom[i-3];
	if ( i==3 )	fprintf(fw,"%+.4f",m);
	else		fprintf(fw,",%+.4f",m);
  }
 return(0);
}

static int fprint_star_px(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%8.3f ",s->psf.gcx);
 return(0);
}
static int fprint_star_py(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%8.3f ",s->psf.gcy);
 return(0);
}
static int fprint_star_pbg(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%8.2f ",s->psf.gbg);
 return(0);
}
static int fprint_star_pamp(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%8.2f ",s->psf.gamp);
 return(0);
}
static int fprint_star_pmoms(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%6.3f ",1.0);
 return(0);
}
static int fprint_star_pmomx(FILE *fw,star *s,formatoutparam *fop)
{
 fprintf(fw,"%6.3f ",0.0);
 return(0);
}

static formatname formatnames[] =
{
 { I_EMPTY,	 1, "-",	NULL,			{ "-", NULL } },	
 { I_ID,	 6, "Ident",	fprint_star_id,		{ "id", NULL } },
	
 { I_IX,	 4, "IX",	fprint_star_ix,		{ "ix", NULL } },
 { I_IY,	 4, "IY",	fprint_star_iy,		{ "iy", NULL } },

 { I_CX,	 8, "C.X",	fprint_star_cx,		{ "cx", NULL } },
 { I_CY,	 8, "C.Y",	fprint_star_cy,		{ "cy", NULL } },
 { I_CBG,	 8, "C.Bg",	fprint_star_cbg,	{ "cbg", NULL } },
 { I_CAMP,	 8, "C.Amp",	fprint_star_camp,	{ "camp", NULL } },
 { I_CMAX,	 8, "C.Max",	fprint_star_cmax,	{ "cmax", NULL } },
 { I_CPIX,	 5, "NPix",	fprint_star_cpix,	{ "cpix", "pix", "npix", NULL } },
 { I_CDEVS,	 6, "C.S",	fprint_star_cdevs,	{ "cs", NULL } },
 { I_CDEVD,	 6, "C.D",	fprint_star_cdevd,	{ "cd", NULL } },
 { I_CDEVK,	 6, "C.K",	fprint_star_cdevk,	{ "ck", NULL } },
 
 { I_X,		 8, "X",	fprint_star_x,		{ "x", NULL } },
 { I_Y,		 8, "Y",	fprint_star_y,		{ "y", NULL } },
 { I_BG,	 8, "Bg",	fprint_star_bg,		{ "b", "bg", "background", NULL } },
 { I_AMP,	 8, "Amp",	fprint_star_amp,	{ "a", "amp", "amplitude",  NULL } },

 { I_DEVS,	 6, "S",	fprint_star_devs,	{ "s", "S", NULL } },
 { I_DEVD,	 6, "D",	fprint_star_devd,	{ "d", "D", NULL } },
 { I_DEVK,	 6, "K",	fprint_star_devk,	{ "k", "K", NULL } },
 { I_MOM,	 8, "Coeff",	fprint_star_mom,	{ "mom", "moments", NULL } },

 { I_DEVL,	 6, "L",	fprint_star_devl,	{ "l", NULL } },

 { I_DSIG,	 6, "sigma",	fprint_star_dsig,	{ "sig", "sigma", NULL } },
 { I_DDEL,	 6, "delta",	fprint_star_ddel,	{ "del", "delta", NULL } },
 { I_DKAP,	 6, "kappa",	fprint_star_dkap,	{ "kap", "kappa", NULL } },
 { I_FWHM,	 5, "FWHM",	fprint_star_fwhm,	{ "f", "fwhm", NULL } },
 { I_ELLIP,	 5, "Ellip",	fprint_star_ellip,	{ "e", "ellip", NULL } },
 { I_PA,	 5, "P.A.",	fprint_star_pa,		{ "p", "pa", NULL } },

 { I_FLUX,	10, "Flux",	fprint_star_flux,	{ "i", "flux", NULL } },
 { I_MAG,	 9, "Magnitude",fprint_star_mag,	{ "m", "mag", "magnitude", NULL } },
 { I_NOISE,	 8, "Noise",	fprint_star_noise,	{ "noise", NULL } },
 { I_SN,	 8, "S/N",	fprint_star_sn,		{ "sn", "spern", "s/n", NULL } },

 { I_PX,	 8, "PSF.X",	fprint_star_px,		{ "px", NULL } },
 { I_PY,	 8, "PSF.Y",	fprint_star_py,		{ "py", NULL } },
 { I_PBG,	 8, "PSF.Bg",	fprint_star_pbg,	{ "pbg", NULL } },
 { I_PAMP,	 8, "PSF.Amp",	fprint_star_pamp,	{ "pamp", NULL } },
 { I_PMOMS,	 6, "PSF.S",	fprint_star_pmoms,	{ "ps", NULL } },
 { I_PMOMD,	 6, "PSF.D",	fprint_star_pmomx,	{ "pd", NULL } },
 { I_PMOMK,	 6, "PSF.K",	fprint_star_pmomx,	{ "pk", NULL } },
 { I_PMOML,	 6, "PSF.L",	fprint_star_pmomx,	{ "pl", NULL } },

 { -1, 		 0, NULL,	NULL,			{ NULL } }
};
	
int *output_format_stars(char *iformat)
{
 int	*oformat,len,alen,t,i,j,id;
 char	*p,token[16];

 alen=128;
 oformat=(int *)malloc(sizeof(int)*alen);
 len=0;

 while ( *iformat )
  {	p=iformat;t=0;
	while ( *p && *p != ',' && t<15 )
	 {	token[t]=*p,p++,t++;		};
	token[t]=0;
	while ( *p && *p != ',' )	p++;
	if ( *p==',' )	p++;
	iformat=p;

	if ( len>=alen-16 )
	 {	alen+=128;
		oformat=(int *)realloc(oformat,sizeof(int)*alen);
	 }

	for ( i=0,id=-1 ; formatnames[i].id>=0 && id<0 ; i++ )
	 {	for ( j=0 ; formatnames[i].tlist[j] != NULL && id<0 ; j++ )
		 {	if ( strcmp(formatnames[i].tlist[j],token)==0 )
				id=formatnames[i].id;
		 }
	 }
	if ( id<0 )
	 {	free(oformat);
		return(NULL);
	 }

	oformat[len]=id;
	len++;
  }

 oformat[len]=-1;

 return(oformat);
}

int write_stars(FILE *fw,star *stars,int nstar,
	int *indx,int is_comment,int *oformat,magflux *mf)
{
 char		*opnt,fbuff[32];
 star		*s;
 int		n,i,j,w,order,ndev,*of,f;
 int		formatlookup[64];
 formatname	*wf;
 formatoutparam	fop;

 for ( i=0 ; i<64 ; i++ )
	formatlookup[i]=-1;
 for ( i=0 ; formatnames[i].id>=0 ; i++ )
  {	j=formatnames[i].id;
	formatlookup[j]=i;
  }

 order=0;
 for ( j=0 ; j<nstar ; j++ )
  {	if ( indx != NULL )	n=indx[j];
	else			n=j;
	s=&stars[n];
	if ( s->shape.order>order )	order=s->shape.order;
  }
 if ( order>0 )	ndev=(order+1)*(order+2)/2-3;
 else		ndev=0;

 if ( is_comment )
  {	fprintf(fw,"#");
	for ( of=oformat ; *of>=0 ; of++ )
	 {	f=*of;
		if ( f<0 || f>=64 )	f=-1;
		else			f=formatlookup[f];

		if ( f>=0 )
		 {	wf=&formatnames[f];
			w=wf->width;
			opnt=wf->comment;
			if ( wf->id==I_MOM )	w=8*ndev-1;
		 }
		else
		 {	w=0;
			opnt=NULL;
		 }
	
		if ( opnt != NULL )
		 {	sprintf(fbuff,"%%%ds ",w);
			fprintf(fw,fbuff,opnt);
		 }
	 }
	fprintf(fw,"\n");

	fprintf(fw,"#");
	for ( of=oformat,n=0 ; *of>=0 ; of++,n++ )
	 {	f=*of;
		if ( f<0 || f>=64 )	f=-1;
		else			f=formatlookup[f];
		if ( f>=0 )
		 {	wf=&formatnames[f];
			if ( wf->id==I_MOM )	w=wf->width*ndev-1;
			else			w=wf->width;
		 }
		else	w=0;

		if ( w>0 )
		 {	sprintf(fbuff,"%%%ds[%2d] ",w-4,n+1);
			fprintf(fw,fbuff,"");
		 }
	 }
	fprintf(fw,"\n");
   }

 for ( j=0 ; j<nstar ; j++ )
  {	if ( indx != NULL )	n=indx[j];
	else			n=j;
	s=&stars[n];

	fop.n=n;
	fop.mf=mf;	

	fprintf(fw," ");
	for ( of=oformat ; *of>=0 ; of++ )
	 {	f=*of;
		if ( f<0 || f>=64 )	f=-1;
		else			f=formatlookup[f];
		if ( f>=0 )	wf=&formatnames[f];
		else		wf=NULL;
		if ( wf==NULL )	continue;
		wf->fprint(fw,s,&fop);
	 }
	fprintf(fw,"\n");
  }
 return(0);
}

