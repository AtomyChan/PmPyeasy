/*****************************************************************************/
/* fistar.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line tool for searching stars and get star parameters on an image.*/
/*****************************************************************************/
#define	FITSH_FISTAR_VERSION		"1.0rc5"
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>

#include <fits/fits.h>

#include "longhelp.h"
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

#include "imgtrans.h"
#include "tensor.h"
#include "index/sort.h"
#include "common.h"

#include "magnitude.h"
#include "stars.h"
#include "psf.h"
#include "psf-base.h"
#include "psf-determine.h"
#include "psf-io.h"

#include "fistar.h"

#ifdef  HAVE_NO_CC_EXTENSION 
#define __extension__ 
#endif 

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int	is_verbose,is_comment;
char	*progbasename;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int fprint_error(char *expr,...)
{
 va_list	ap;
 fprintf(stderr,"%s: error: ",progbasename);
 va_start(ap,expr);
 vfprintf(stderr,expr,ap);
 va_end(ap);
 fprintf(stderr,"\n");
 return(0);
}

int fprint_warning(char *expr,...)
{
 va_list	ap;
 fprintf(stderr,"%s: warning: ",progbasename);
 va_start(ap,expr);
 vfprintf(stderr,expr,ap);
 va_end(ap);
 fprintf(stderr,"\n");
 return(0);
}

int fprint_info(char *expr,...)
{
 va_list	ap;
 fprintf(stderr,"%s: ",progbasename);
 va_start(ap,expr);
 vfprintf(stderr,expr,ap);
 va_end(ap);
 fprintf(stderr,"\n");
 return(0);
}

/*****************************************************************************/

#define		ALG_PP		1
#define		ALG_TRB		2
#define		ALG_LNK		3
#define		ALG_BIQ		4

#define		COLL_REFINE_POSITION	1	/* 0x01 */
#define		COLL_REFINE_SHAPE	2	/* 0x02 */

typedef struct
 { 	int	niter;
	int	refinelevel;
	double	bhsize;
 } collectivefit;
				/* basic properties of star searching 	*/
typedef struct			/* and star model fitting...		*/
 {	int		is_fit_model;
	int		is_determine_psf;
	starfit		sfp;

	starmodelfit	smfs[16];
	int		nsmf;

	range		srcrange;
	spatial		imgbg;
	int		algorithm;

	candidate	*cands;
	int		ncand;
	initcand	*icands;
	int		nicand;
	star		*stars;
	int		nstar;

	collectivefit	cf;
	

 } starsearch;

/*****************************************************************************/

int make_star_candidates(initcand *icands,int nicand,double rad,
	fitsimage *img,char **mask,candidate **rcands,int *rncand)
{
 candidate	*cands,*wc;
 int		ncand,i,ix1,ix2,iy1,iy2,ix,iy,irad,k,l,t0,t1;
 int		sx,sy;
 double		dx,dy,rad2,dy2,f0;
 initcand	*ic;
 ipoint		*ipoints;
 int		nipoint;
 starshape	shp;

 if ( img==NULL || img->data==NULL )	return(-1);
 sx=img->sx,sy=img->sy;
 if ( sx<=0 || sy<=0 )			return(-1);

 cands=NULL;
 ncand=0;
 
 irad=2+(int)rad;
 rad2=rad*rad;

 for ( i=0 ; i<nicand && icands != NULL ; i++ )
  {	ic=&icands[i];
	ix=(int)floor(ic->x);
	iy=(int)floor(ic->y);
	ix1=ix-irad,ix2=ix+irad;
	iy1=iy-irad,iy2=iy+irad;
	ipoints=(ipoint *)malloc(sizeof(ipoint)*((2*irad+1)*(2*irad+1)));
	nipoint=0,
	t0=t1=0;
	for ( k=iy1 ; k<=iy2 ; k++ )
	 {	
		dy=(ic->y)-((double)k+0.5);
		dy2=dy*dy;
		for ( l=ix1 ; l<=ix2 ; l++ )
		 {	dx=(ic->x)-((double)l+0.5);
			if ( dx*dx+dy2>rad2 )	continue;

			t0++;
			if ( k<0 || sy<=k )	continue;
			if ( l<0 || sx<=l )	continue;
			if ( mask != NULL && mask[k][l] )	continue;
			t1++;

			ipoints[nipoint].x=l,
			ipoints[nipoint].y=k;
			nipoint++;
		 }
	 }
	if ( nipoint==0 || t1<t0 )
	 {	free(ipoints);
		continue;
	 }
	cands=(candidate *)realloc(cands,sizeof(candidate)*(ncand+1));
	wc=&cands[ncand];

	wc->ix=ix,wc->cx=ic->x,
	wc->iy=iy,wc->cy=ic->y;
	
	wc->bg=ic->bg;

	wc->sxx=ic->s+ic->d;
	wc->syy=ic->s-ic->d;
	wc->sxy=ic->k;

	shp.model=SHAPE_ELLIPTIC;
	shp.gs=ic->s,
	shp.gd=ic->d,
	shp.gk=ic->k;
	f0=star_get_unity_flux(&shp);
	if ( f0>0.0 )	wc->peak=wc->amp=ic->flux/f0;
	else		wc->peak=wc->amp=0.0;
	wc->flux=ic->flux;

	wc->ipoints=ipoints;
	wc->nipoint=nipoint;

	wc->area=0.0;

	wc->flags=0;
	wc->marked=0;
	ncand++;
  }

 if ( rcands != NULL )	*rcands=cands;
 if ( rncand != NULL )	*rncand=ncand;

 return(0);
}

/*****************************************************************************/

int compare_x(int i,int j,void *p)
{
 if ( ((star *)p)[i].location.gcx < ((star *)p)[j].location.gcx ) return(1);
 if ( ((star *)p)[i].location.gcx > ((star *)p)[j].location.gcx ) return(-1);
 return(0);
}
int compare_y(int i,int j,void *p)
{
 if ( ((star *)p)[i].location.gcy < ((star *)p)[j].location.gcy ) return(1);
 if ( ((star *)p)[i].location.gcy > ((star *)p)[j].location.gcy ) return(-1);
 return(0);
}
int compare_peak(int i,int j,void *p)
{
 if ( ((star *)p)[i].cand->peak < ((star *)p)[j].cand->peak )	return(1);
 if ( ((star *)p)[i].cand->peak > ((star *)p)[j].cand->peak )	return(-1);
 return(0);
}
int compare_fwhm(int i,int j,void *p)
{
 if ( ((star *)p)[i].gfwhm < ((star *)p)[j].gfwhm )	return(1);
 if ( ((star *)p)[i].gfwhm > ((star *)p)[j].gfwhm )	return(-1);
 return(0);
}
int compare_amp(int i,int j,void *p)
{
 if ( ((star *)p)[i].location.gamp < ((star *)p)[j].location.gamp ) return(1);
 if ( ((star *)p)[i].location.gamp > ((star *)p)[j].location.gamp ) return(-1);
 return(0);
}
int compare_flux(int i,int j,void *p)
{
 if ( ((star *)p)[i].flux < ((star *)p)[j].flux )	return(1);
 if ( ((star *)p)[i].flux > ((star *)p)[j].flux )	return(-1);
 return(0);
}
int compare_noise(int i,int j,void *p)
{
 if ( ((star *)p)[i].cand==NULL || ((star *)p)[j].cand==NULL )	return(0);
 if ( ((star *)p)[i].cand->noise < ((star *)p)[j].cand->noise )	return(1);
 if ( ((star *)p)[i].cand->noise > ((star *)p)[j].cand->noise )	return(-1);
 return(0);
}
int compare_sn(int i,int j,void *p)
{
 double	 sn1,sn2;
 if ( ((star *)p)[i].cand==NULL || ((star *)p)[j].cand==NULL )
	return(0);
 else if ( ((star *)p)[i].cand->noise<=0.0 || ((star *)p)[j].cand->noise<=0.0 )
	return(0);
 else
  {	sn1=((star *)p)[i].flux/((star *)p)[i].cand->noise;
	sn2=((star *)p)[j].flux/((star *)p)[j].cand->noise;
	if ( sn1<sn2 )		return(1);
	else if ( sn1>sn2 )	return(-1);
	else			return(0);
  }
}

/*****************************************************************************/

int fistar_determine_psf(fitsimage *img,char **mask,starsearch *ss,psfdetermine *pdparam,psf *tpd)
{
 psfcandidate	*pcands,*wc;
 int		npcand,n,i,nipoint,ii,niter;
 star		*ws;

 if ( is_verbose )
  {	i=pdparam->hsize*2+1;
	fprintf(stderr,"Determination of PSF [%dx%d]x(%dx%d), "
		"spatial order: %d ... ",
		i,i,pdparam->grid,pdparam->grid,pdparam->order);
	fflush(stderr);
  }

 npcand=ss->nstar;
 pcands=(psfcandidate *)malloc(sizeof(psfcandidate)*npcand);

 for ( n=0 ; n<npcand ; n++ )
  {	wc=&pcands[n];
	ws=&ss->stars[n];

	wc->nipoint=nipoint=ws->cand->nipoint;
	wc->ipoints=ws->cand->ipoints;
	wc->yvals=(double *)malloc(sizeof(double)*nipoint);
	memset(wc->yvals,0,sizeof(double)*nipoint);
	drawback_model(wc->ipoints,nipoint,wc->yvals,
		&ws->location,&ws->shape,+1.0);
	wc->bg =ws->location.gbg;
	wc->amp=ws->location.gamp;
	wc->cx =ws->location.gcx;
	wc->cy =ws->location.gcy;
  }

 niter=0;

 for ( ii=0 ; ii<=niter ; ii++ )
  {	psf_determine(img,mask,pcands,npcand,1,pdparam,tpd);
	psf_bgamp_fit(img,mask,pcands,npcand,1,tpd);
	if ( ii<niter )
	 {	for ( n=0 ; n<npcand ; n++ )
		 {	wc=&pcands[n];
			memset(wc->yvals,0,sizeof(double)*wc->nipoint);
			drawback_psf(wc->ipoints,wc->nipoint,wc->yvals,
				wc->cx,wc->cy,wc->amp,tpd,+1.0);
		 }
	 }
  }

 if ( is_verbose )
  {	fprintf(stderr,"done.\n");		}

 for ( n=0 ; n<npcand ; n++ )
  {	wc=&pcands[n];
	ws=&ss->stars[n];
	ws->psf.gcx=wc->cx,
	ws->psf.gcy=wc->cy;
	ws->psf.gamp=wc->amp;
	ws->psf.gbg =wc->bg ;
  }

 for ( n=npcand-1 ; n>=0 ; n-- )
  {	wc=&pcands[n];
	if ( wc->yvals != NULL )	free(wc->yvals);
  }

 free(pcands);

 return(0);
}

/*****************************************************************************/

int fprint_fistar_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\tfistar [-h|--help|--long-help|--wiki-help] [--version[-short]]\n"
"\t[[-i|--input]<input.fits>[<F>] [--frame <F>]] [-o|--output <out.dat>]\n"
"\t[-V|--verbose] [-C|-C|--comment]\n");
 fprintf(fw,
"General parameters:\n"
"\t[--skynoise|-d <skysigma>]\n"
"\t[--[flux-]threshold|-t|-f <threshold>] [-p|--prominence <crit.prom>]\n"
"\t[-s|--sort {x|y|peak|fwhm|amp|flux|noise|s/n}] [--mag-flux <mag>,<flux>]\n"
"\t[-M|--input-mask <mask.hdr>] [--output-background <bg.fits>]\n"
"\t[-r|--shrink <pre-shrink factor>]\n"
"\t[-F|--format <comma-separated-format>] [--section <x1>:<x2>,<y1>:<y2>]\n"
"\t[--output-mark <mark> [--mark {dot|square|circle} --mark-size <size>]]\n");
 fprintf(fw,
"Format arguments (after --format, separated by commas):\n"
"\t- id,ix,iy (ID, integer part of candidate coordinates)\n"
"\t- cx,cy,cbg,camp,cmax,npix,cs,cd,ck (candidate parameters)\n"
"\t- x,y,bg,amp,s,d,k,mom,l,sigma,delta,kappa,fwhm,ellip,pa,\n"
"\t  flux,noise,s/n,magnitude\n"
"\t  (star parameters: location, various shape parameters, and est'd flux)\n"
"\t- px,py,pbg,pamp,ps,pd,pk,pl (PSF parameters)\n");
 fprintf(fw,
"Input star candidates (conflicts with candidate searching):\n"
"\t[-C|--input-candidates <file> [--col-xy <>,<>] [--col-shape <...>]]\n"
"\t[-R|--candidate-radius <radius>]\n");
 fprintf(fw,
"Input position list (also conflicts with candidate searching):\n"
"\t[-P|--input-positions <file> [--col-id <>] [--col-xy <>,<>]]\n"
"\t[-g|--gain <gain>]\n");
 fprintf(fw,
"Candidate search parameters (conflicts with candidates read from file):\n"
"\t[--algorithm {parabolapeak|uplink}] [--only-candidates]\n"
"Modelling and fit tuning parameters:\n"
"\t[--model {gauss,elliptic,deviated[,order=<N>]}\n"
"\t         [--iterations {symmetric=<n>,general=<n>}] ]\n"
"\t[--collective-fit {iterations=<n>,[position],[shape]|blockhalfsize=<h>}]\n");
 fprintf(fw,
"Parameters of PSF determination:\n"
"\t[--psf grid=<g>,symmetrize,halfsize=<h>,order=<o>\n"
"\t   native[,spline]|integral[,kappa=<k>]|circle[,~width=<w>,~order=<o>]]\n"
"\t[--output-psf <output-psf.fits>]\n");
 fprintf(fw,
"Default parameters are:\n"
"\t--algorithm uplink --model elliptic --iterations symmetric=4,general=2\n"
"\t--threshold 100 --mark dot --mark-size 2 --mag-flux 10,10000\n");
 return(0);
}

longhelp_entry fistar_long_help[]=
{
 LONGHELP_OPTIONS,

 { "General options:", NULL },
 { "-h, --help",
	"Gives general summary about the command line options." },
 { "--long-help, --help-long",
        "Gives a detailed list of command line options." },
 { "--wiki-help, --help-wiki, --mediawiki-help, --help-mediawiki",
        "Gives a detailed list of command line options in Mediawiki format." },
 { "--version, --version-short, --short-version",
 	"Gives some version information about the program." },
 { "-i, --input <image file>",
	"Name of input file from which the sources are extracted." },
 { "-o, --output <data file>",
	"Name of the data file where the list of the extracted sources "
	"and their respective characteristics are written." },

 { "Basic source detection and characterization options:", NULL },
 { "-t, --threshold <threshold>",
	"Detection peak threshold, in ADUs." },
 { "-f, --flux-threshold <flux threshold>",
	"Detection flux threshold, in ADUs." },
 { "--algorithm <algorithm>",
	"Source cadidate extraction algorithm. It can be \"uplink\" (default) "
	"or \"parabolapeak\"." },
 { "-p, --prominence <prominence>",
	"Critical relative prominence parameter used in the \"uplink\" "
	"algorithm. The default is to use no prominence based pixel grouping." },
 { "-r, --shrink <shrink factor>",
	"Shrink factor applied before star candidate detection. The image is "
	"shrunk by this specific factor and after the detection, the "
	"candidate coordinates are multiplied by this value." },
 { "--only-candidates",
	"Do not involve any analitic model funcion during the derivation of "
	"the centroid coordinates and shape parameters, but derive these "
	"from some sort of simple pixel statistics." },
 { "--model <model>[order=<order>]",
	"Analytic model function used in the second stage of source extraction. "
	"This can be \"gauss\" (symmetic Gaussian profile), \"elliptic\" "
	"(assymetric Gaussian profile) and \"deviated\" (Gaussian profile "
	"multiplied by a polynomial up to the specified order)." },
 { "--iterations symmetric=<n>,general=<n>",
	"The number of Levenberg-Marquardt iterations during the analytic "
	"model fit. The fit is done in two substeps, first the symmetric "
	"profile parameters are derived only, which step is followed by "
	"the fit for all of the shape parameters. " },

 { "-s, --sort {x|y|peak|fwhm|amp|flux|noise|s/n}",
	"Sort the output list of extracted sources by X or Y coordinate, "
	"peak flux (no background level subtracted), profile FWHM, peak "
	"intensity (background level is subtracted), total flux, noise level "
	"or signal-to-noise ratio; respectively. " },
 { "--mag-flux <mag>,<flux>",
	"Magnitude - flux conversion level. The specified magnitude will "
	"be equivalent to the specified flux level." },

 { "-M, --input-mask <image file>",
	"Input mask file to co-add to the mask of the input image. Useful for "
	"marking pixels to be excluded from source extraction process "
	"beyond the ones which are previously marked in the input image." },

 { "-F, --format <format>",
	"Comma separated list of format tags, for formatting the "
	"output list of extracted sources. Each row represents an "
	"extracted source while the format specified here defines the "
	"quantities listed in each row of the output file. See \"Format "
	"tags\" for more details. " },

 { "Format tags:", NULL },
 { "id", 
	"An unique identifier for the source (an integer number, in fact)." },
 { "ix",
	"Integer X coordinate for the centroid pixel" },
 { "iy",
	"Integer Y coordinate for the centroid pixel" },
 { "x",
	"Centroid X coordinate in native coordinate convention." },
 { "y",
	"Centroid Y coordinate in native coordinate convention." },
 { "bg",
	"Background level" },
 { "amp",
	"Peak amplitude" },
 { "S", 
	"Gaussian momenum for the stellar profile (S=1/sigma^2)" },
 { "D", 
	"plus-shaped momentum for the stellar profile" },
 { "K", 
	"cross-shaped momentum for the stellar profile" },
 { "sigma", 
	"sigma parameter for the stellar profile (FWHM is roughly 2.35 * sigma)" },
 { "delta", 
	"delta (plus-shaped deviance) parameter for the stellar profile" },
 { "kappa", 
	"kappa (cross-shaped deviance) parameter for the stellar profile" },
 { "fwhm", 
	"full width at half magnitude (FWHM) of the stellar profile" },
 { "ellip", 
	"ellipticity of the stellar profile" },
 { "pa", 
	"position angle of the stellar profile" },
 { "flux",
	"Total flux of the source" },
 { "nosie",
	"Noise level of the source" },
 { "s/n",
	"Signal-to-noise ratio of the detection" },
 { "magnitude",
	"Brightness of the source in magnitudes (see also --mag-flux)" },
 { "cx",
	"Candidate centroid X coordinate (derived from pixel flux statistics)." },
 { "cy",
	"Candidate centroid Y coordinate (derived from pixel flux statistics)." },
 { "cbg",
	"Background level for the source candidate" },
 { "camp",
	"Peak amplitude of the source candidate" },
 { "cmax",
	"Maximum intensity on the source cadidate" },
 { "cs", 
	"Gaussian momenum for the source cadidate, derived from pixel flux statistics" },
 { "cd", 
	"plus-shaped momentum for the source cadidate, derived from pixel flux statistics" },
 { "ck", 
	"cross-shaped momentum for the source cadidate, derived from pixel flux statistics" },
 { "npix",
	"number of pixels assigned to the detetcted source" },

 { "Obtaining the point-spread function (PSF):", NULL },
 { "--psf <parameters>",
	"Comma-separated list of parameters related to the PSF fitting:" },
 { "grid=<grid>",
	"subgrid size for the PSF" },
 { "halfsize=<half size>",
	"half-size of the PSF stamp, the full size of the stamp will always be "
	"2*<half size>+1 and the PSF itself is centered at the center of the "
	"central pixel." },
 { "order=<order>",
	"order of spatial variations in the PSF" },
 { "symmetrize",
	"symmetrize the resulting PSF (i.e. fit a symmetric PSF instead of a normal one)" },
 { "spline",
	"use a spline interpolation during the determination of the PSF" },
 { "--output-psf <output PSF FITS file>",
	"Name of the output file where the PSF is saved in FITS format. "
	"The PSF is stored in a 3 dimensional (a.k.a. \"data cube\") structure "
	"where the z-axis is for the various polynomial coefficients describing "
	"the spatially varied PSF. " },

 { "Alternate sources for object candidates:", NULL },
 { "-C, --input-cadidates <cadidate list file>",
	__extension__ 
	"Name of input cadidate list file. If such a file is defined in the "
	"command line, the cadidates are not searched by the build-in "
	"algorithms. Instead, the centroids are read from this file and the "
	"pixels for each object are defined within a certain radius from "
	"this centroid (see -r|--cadidate-radius also). The role of this "
	"option is twofold. First, it is suitable if only some of the "
	"sources have to be modelled with an analytic function; second, "
	"PSF determination can be done only a previously cleaned list of "
	"sources, in the case when there might be contaminating sources too." },
 { "-r, --candidate-radius <radius>",	
	"This option defines a distance, where pixels within this are "
	"assigned to the candidate. " },
 { "--col-xy <colx>,<coly>",
	"Column indices for X and Y centroid coordinates. "},

 { NULL, NULL }
};

int fprint_fistar_long_help(FILE *fw,int is_wiki)
{
 char	*synopsis=
	"fistar [options] <input> [-o|--output <output>]";
 char	*description=
	__extension__
	"The main purpose of this program is to detect and extract sources (i.e. "
	"star-like objects) from astronomical images. The source detection and "
	"extraction are based on three major steps. First, pixel groups are derived "
	"which are possibly belong to the sources (these preliminary detections are "
	"callad source \"candidates\"). Second, these candidates are modelled with "
	"some sort of analytic model funcion, in order to derive more precise "
	"centroid coordinates and shape parameters. The last step is to extract the "
	"point-spread function (PSF) for the image, based on the detected and "
	"modelled sources. Basically, the input for this program must be an "
	"astronomical image while the output is the list of detected and extracted "
	"sources and their respective characteristics.";

 fprint_generic_long_help(fw,is_wiki,fistar_long_help,synopsis,description);

 return(0);
}

static char *default_format="id,ix,iy,cx,cy,cmax,x,y,b,a,fwhm,ellip,pa,s,d,k,flux";

int main(int argc,char *argv[])
{
 fits		*img;
 char		**mask;

 FILE		*fr,*fw;
 int		i,sx,sy,*indx;
 double		skysigma,threshold,fluxthreshold,critical_prominence,candradius;
 int		shrinkfactor;
 char 		*outname,*inimg,*markimgname,*areaimgname,*incandname,*inposname,
		*iterpar,*algpar,*modelpar,*psfoutname,*psfpar,*collfitpar;
 char		**inmasklist,*formatpar,*basename;
 int		sort,*oformat;
 int		mark_symbol,mark_size;
 int		is_help,frameno;
 colinfo	col;
 starsearch	ss;
 psf		tpd;
 psfdetermine	pdparam;
 magflux	mf0;
 double		gain;

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 inimg=outname=NULL;
 skysigma=0.0;
 threshold=100.0;fluxthreshold=0.0;
 critical_prominence=0.0;

 ss.sfp.iter_symmetric=4;
 ss.sfp.iter_general=2;

 is_verbose=is_help=is_comment=0;
 markimgname=NULL;mark_symbol=MARK_SYM_DOTS;mark_size=2;
 psfoutname=psfpar=NULL;
 collfitpar=NULL;
 incandname=inposname=areaimgname=NULL;inmasklist=NULL;formatpar=NULL;
 frameno=-1;
 sort=-1;

 candradius=0.0;
 col.id=3;
 col.x=1,col.y=2;
 col.s=3,col.d=-1,col.k=-1;
 col.flux=-1,col.mag=-1,col.bg=-1;

 ss.srcrange.xmin=ss.srcrange.xmax=ss.srcrange.ymin=ss.srcrange.ymax=0;

 ss.is_fit_model=1;
 ss.is_determine_psf=0;
 ss.algorithm=ALG_LNK;	    /* default search method is the uplink algorithm */
 ss.smfs[0].model=SHAPE_ELLIPTIC;    /* default model is elliptical gaussian */
 ss.smfs[0].order=2;
 ss.nsmf=1;

 algpar=modelpar=iterpar=NULL;

 ss.stars=NULL;ss.nstar=0;
 ss.cands=NULL;ss.ncand=0;

 ss.cf.niter=-1;
 ss.cf.refinelevel=0;
 ss.cf.bhsize=10.0;

 mf0.intensity=10000.0;
 mf0.magnitude=10.0;
 gain=1.0;
 shrinkfactor=1;

 i=scanarg(argc,argv,SCANARG_ALLOW_FLAGS,
	"--version:%NS-1f%q",&is_help,
	"--version-short|--short-version:%NS-2f%q",&is_help,
	"-i|--input:%s",&inimg,
	"--frame:%d",&frameno,
	"--section:%d:%d,%d:%d",&ss.srcrange.xmin,&ss.srcrange.xmax,
				&ss.srcrange.ymin,&ss.srcrange.ymax,
	"--algorithm:%s",&algpar,
	"-p|--prominence:%g",&critical_prominence,
	"--model:%s",&modelpar,
	"-r|--shrink:%d",&shrinkfactor,
	"-h|--help:%f",&is_help,
        "--long-help|--help-long:%SN2f%q",&is_help,
	"--mediawiki-help|--help-mediawiki|--wiki-help|--help-wiki:%SN3f%q",&is_help,
	"-o|--output:%s",&outname,
	"-M|--input-mask:%t",&inmasklist,

	"-C|--input-candidates:%s",&incandname,
	"-P|--input-positions:%s",&inposname,
	"-R|--candidate-radius:%g",&candradius,	/* equiv to 'hsize'... */
	"--col-id:%Pd",&col.id,
	"--col-xy:%Pd,%Pd",&col.x,&col.y,
	"--col-gauss|--col-shape:%Pd",&col.s,
	"--col-elliptic|--col-shape-dev:%Pd,%Pd,%Pd",&col.s,&col.d,&col.k,
	"--col-flux:%Pd",&col.flux,
	"--col-bg:%Pd",&col.bg,
	"--col-mag:%Pd",&col.mag,

	"-d|--skysigma|--skynoise:%g",&skysigma,
	"-t|--threshold:%g",&threshold,
	"-f|--flux-threshold:%g",&fluxthreshold,
	"--only-candidates:%NS0f",&ss.is_fit_model,
	"--iterations:%s",&iterpar,
	"--collective-fit:%s",&collfitpar,
	"-s|--sort:%(x,y,peak,fwhm,amp,flux,noise,s/n)",&sort,
	"--mag-flux:%g,%g",&mf0.magnitude,&mf0.intensity,
	"-g|--gain:%g",&gain,
	"-F|--format:%s",&formatpar,
	"--output-marked|--mark-output:%s",&markimgname,
	"--mark=dot:%NS0f",&mark_symbol,
	"--mark=square:%NS1f",&mark_symbol,
	"--mark=circle:%NS2f",&mark_symbol,
	"--mark:%(dot,square,circle)",&mark_symbol,
	"--mark-size:%d",&mark_size,
	"--output-area:%s",&areaimgname,
	"--output-psf|--psf-output:%s",&psfoutname,
	"--psf:%f%s",&ss.is_determine_psf,&psfpar,
	"--comment:%i",&is_comment,"(C):%i",&is_comment,
	"--verbose:%i",&is_verbose,"(V):%i",&is_verbose,
 	"-:%w",&inimg,
	"-*|+*:%e",
	"*:%w",&inimg,
	NULL);

 if ( i )		
  {	fprint_error("invalid command line argument near '%s'",argv[i]);
	return(1);
  }
 else if ( is_help<0 )
  {	fprint_generic_version(stdout,argv[0],"fistar",FITSH_FISTAR_VERSION,is_help);
	return(0);
  }
 else if ( 1<is_help )
  {	fprint_fistar_long_help(stdout,2<is_help);
	return(0);
  }
 else if ( is_help )
  {	fprint_fistar_usage(stdout);
	return(0);
  }

 col.x--,col.y--;
 col.s--,col.id--;
 if ( col.d>0 )		col.d--;
 if ( col.k>0 )		col.k--;
 if ( col.flux>0 )	col.flux--;
 if ( col.mag>0 )	col.mag--;
 if ( col.bg>0 )	col.bg--;

 if ( iterpar != NULL )
  {	ss.sfp.iter_general=0;
	i=scanpar(iterpar,SCANPAR_DEFAULT,
		"symmetric:%d",&ss.sfp.iter_symmetric,
		"general:%d",&ss.sfp.iter_general,
		NULL);
	if ( i )	
	 {	fprint_error("invalid fit iteration parameter '%s'",iterpar);
		return(1);
	 }
  }
 if ( modelpar != NULL )
  {	int	model,order,i,n,k;
	char	*lmodelpar,*cmd[16];

	lmodelpar=strdup(modelpar);
	n=tokenize_char(lmodelpar,cmd,'+',15);

	ss.nsmf=0;
	for ( i=0 ; i<n ; i++ )
	 {	model=0;order=2;
		k=scanpar(cmd[i],SCANPAR_DEFAULT,
			"gauss:   " SNf(SHAPE_GAUSS),	&model,
			"elliptic:" SNf(SHAPE_ELLIPTIC),&model,
			"deviated:" SNf(SHAPE_DEVIATED),&model,
			"order:%d",&order,
			NULL);
		if ( k )	
		 {	fprint_error("invalid model specification '%s'",cmd[i]);
			return(1);
		 }
		switch ( model )
		 {	case 1: ss.smfs[i].model=SHAPE_GAUSS;break;
			case 2: ss.smfs[i].model=SHAPE_ELLIPTIC;break;
			case 3: ss.smfs[i].model=SHAPE_DEVIATED;break;
			default: 
				fprint_error("invalid model specification '%s'",cmd[i]);
				return(1);	
		 }

		if ( order<2 )
			order=2;
		else if ( order>MAX_DEVIATION_ORDER )
			order=MAX_DEVIATION_ORDER;

		ss.smfs[i].order=order;
		ss.nsmf++;
	 }
	free(lmodelpar);
  }

 if ( algpar != NULL )
  {	if ( strcmp(algpar,"parabolapeak")==0 )		ss.algorithm=ALG_PP;
	else if ( strcmp(algpar,"triangulation")==0 )	ss.algorithm=ALG_TRB;
	else if ( strcmp(algpar,"uplink")==0 )		ss.algorithm=ALG_LNK;
	else if ( strcmp(algpar,"biquad")==0 )		ss.algorithm=ALG_BIQ;
	else	
	 {	fprint_error("invalid star detection algorithm '%s'",algpar);
		return(1);
	 }
  }
 if ( skysigma <= 0.0 && ss.algorithm==ALG_PP )
  {	fprint_error("unknown sky (background) sigma, it is required by the 'parabolapeak' algorithm");
	return(1);
  }

 pdparam.hsize=4;
 pdparam.grid=4;
 pdparam.order=0;
 pdparam.type=PSF_DET_NATIVE;
 memset(&pdparam.param,0,sizeof(pdparam.param));
 pdparam.is_symmetrize=0;

 if ( psfpar != NULL )
  {	i=scanpar(psfpar,SCANPAR_DEFAULT,
		"halfsize:%d",&pdparam.hsize,
		"grid:%d",&pdparam.grid,
		"order:%d",&pdparam.order,
		"native:"  SNf(PSF_DET_NATIVE)	,&pdparam.type,
		"integral:"SNf(PSF_DET_INTEGRAL),&pdparam.type,
		"circle:"  SNf(PSF_DET_CIRCLE)	,&pdparam.type,
		"spline|biquad:"SNf(PSF_DET_NATIVE)"%f"	 ,&pdparam.type,&pdparam.param.native.use_biquad,
		"kappa:"        SNf(PSF_DET_INTEGRAL)"%g",&pdparam.type,&pdparam.param.integral.kappa,
		"circlewidth:"  SNf(PSF_DET_CIRCLE)"%g"	 ,&pdparam.type,&pdparam.param.circle.width,
		"circleorder:"  SNf(PSF_DET_CIRCLE)"%d"	 ,&pdparam.type,&pdparam.param.circle.order,
		"symmetrize:%f",&pdparam.type,&pdparam.is_symmetrize,
		NULL);
	if ( i )	
	 {	fprint_error("invalid PSF specification in '%s'",psfpar);
		return(1);
	 }
  }

 if ( collfitpar != NULL )
  {	int	niter;
	ss.cf.niter=0;
	niter=-1;
	i=scanpar(collfitpar,SCANPAR_DEFAULT,
		"iterations:%d",&niter,
		"position:%SN1f",&ss.cf.refinelevel,
		"shape:%SN2f",&ss.cf.refinelevel,
		"bhsize|blockhalfsize:%g",&ss.cf.bhsize,
		NULL);
	if ( i )	
	 {	fprint_error("invalid collective fit parameter in '%s'",collfitpar);
		return(1);
	 }
	if ( niter>=0 )	ss.cf.niter=niter+1;
  }

 if ( formatpar != NULL )	oformat=output_format_stars(formatpar);	
 else			 	oformat=output_format_stars(default_format);
 if ( oformat==NULL )		
  {	fprint_error("invalid format parameter '%s'",formatpar);
	return(1);
  }

 if ( inimg==NULL )
  {	fr=stdin;frameno=0;	}
 else
  {	basename=fits_basename(inimg,&frameno);
	if ( (fr=fopenread(basename))==NULL )
	 {	fprint_error("unable to open input file '%s'.",basename);
		return(1);
	 }
  }
 img=fits_read_frame_to_image(fr,frameno);
 fcloseread(fr);
 if ( img==NULL )
  {	fprint_error("unable to interpret input data as a FITS image.");
	return(1);
  }
 else if ( img->i.dim != 2 )
  {	fprint_error("image dimension is differ from 2.");
	return(1);
  }
 fits_rescale(img);

 sx=img->i.sx,
 sy=img->i.sy;
 mask=fits_mask_read_from_header(&img->header,sx,sy,NULL);
 if ( inmasklist != NULL )
  {	if ( join_masks_from_files(mask,sx,sy,inmasklist) )
	 {	fprint_error("unable to read one of the input mask files");
		return(1);
	 }	
  }
 fits_mask_mark_nans(&img->i,mask,MASK_NAN);
	
 if ( is_verbose )	fprintf(stderr,"[%dx%d]\n",sx,sy);

 if ( incandname != NULL && candradius <= 0.0 )	candradius=2.0;

 if ( incandname == NULL )
  {	range	*sr;
	sr=&ss.srcrange;
	if ( sr->xmax<=sr->xmin || sr->ymax<=sr->ymin )	sr=NULL;

	switch ( ss.algorithm )
	 {   /* parabolapeak */
	     case ALG_PP:		
		determine_background(&img->i,&ss.imgbg,3,3,2);
		search_star_candidates(&img->i,mask,&ss.cands,&ss.ncand,sr,threshold,&ss.imgbg,skysigma);
		markout_candidates(&img->i,mask,ss.cands,ss.ncand);
		cleanup_candlist(&ss.cands,&ss.ncand);
		break;
	     /* uplink */
	     case ALG_LNK:
		if ( fluxthreshold>0.0 )	/* use fluxthreshold instead */
			threshold=0.0; 
		search_star_candidates_link(&img->i,mask,&ss.cands,&ss.ncand,sr,
			threshold,fluxthreshold,critical_prominence);
		refine_candidate_params(&img->i,ss.cands,ss.ncand);
		break;
	     /* triangulation-based */
	     /*case ALG_TRB:
		search_star_candidates_trb(img,mask,&ss.cands,&ss.ncand,sr,threshold);
		break;*/
	     default:
		fprint_error("desired candidate search algorihm has not been implemented yet.");
		return(1);
		break;
	 }

	if ( is_verbose )	fprint_info("candidates (found): ncand=%d\n",ss.ncand);
  }
 else
  {	FILE	*fr;
	fr=fopenread(incandname);
	if ( fr==NULL )	
	 {	fprint_error("unable to open input list file '%s'",incandname);
		return(1);
	 }
	read_star_candidates(fr,&col,&ss.icands,&ss.nicand,&mf0);
	if ( is_verbose )	fprint_info("candidates (read): nicand=%d\n",ss.nicand);
	make_star_candidates(ss.icands,ss.nicand,candradius,&img->i,mask,&ss.cands,&ss.ncand);
	fclose(fr);
  }

 if ( inposname != NULL )
  {	iposition	*ips;
	int		nip;
	FILE		*fr,*fw;
	int		**refarr;
	int		x,y,i,j;
	double		dx,dy;
	double		magn,merr,flux,ferr;
	char		*name;
	candidate	*wc;

	fr=fopenread(inposname);
	if ( fr==NULL )	
	 {	fprint_error("unable to open input position list '%s'",inposname);
		return(1);
	 }
	read_input_position_list(fr,&col,&ips,&nip);
	fcloseread(fr);

	refarr=(int **)tensor_alloc_2d(int,sx,sy);
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	refarr[i][j]=-1;		}
	 }
	for ( i=0 ; i<ss.ncand ; i++ )
	 {	for ( j=0 ; j<ss.cands[i].nipoint ; j++ )
		 {	x=ss.cands[i].ipoints[j].x;
			y=ss.cands[i].ipoints[j].y;
			refarr[y][x]=i;
		 }
	 }

	if ( outname == NULL )	fw=stdout;
	else			fw=fopenwrite(outname);
	if ( fw==NULL )	
	 {	fprint_error("unable to create output file '%s'",outname);
		return(1);
	 }

	for ( i=0 ; i<nip ; i++ )
	 {	dx=ips[i].x;x=(int)(floor(dx));
		dy=ips[i].y;y=(int)(floor(dy));
		name=ips[i].name;
		if ( name==NULL )	name="-";
		if ( x<0 || x>=sx || y<0 || y>=sy || ( mask != NULL && mask[y][x] ) || refarr[y][x]<0 )
		 {	wc=NULL;
			flux=ferr=0.0;
		 }
		else
		 {	wc=&ss.cands[refarr[y][x]];
			flux=wc->flux;
			if ( flux<=0.0 )
				flux=ferr=0.0;
			else
				ferr=sqrt(flux/gain+wc->noise*wc->nipoint);
		 }
		
		if ( flux>0.0 )
		 {	flux_to_mag_magerr(flux,ferr,&mf0,&magn,&merr);
			fprintf(fw,"%20s %10.3f %10.3f G--- %8.4f %8.4f\n",name,dx,dy,magn,merr);
		 }
		else
		 {	fprintf(fw,"%20s %10.3f %10.3f X--- %8s %8s\n",name,dx,dy,"-","-");	}
 	 }

	fcloseread(fw);
	fits_free(img);
	return(0);
  }
 
 if ( ss.cf.niter<0 )
  {	if ( ss.is_fit_model )
	 {	ss.sfp.fit_flags=FIT_XY|FIT_AB|FIT_WIDTH|FIT_DEVIATION;
 		fit_star_single_model(&img->i,mask,ss.cands,ss.ncand,&ss.stars,&ss.nstar,&ss.sfp,ss.smfs[0].model,ss.smfs[0].order);
		if ( is_verbose )	fprintf(stderr,"nstar=%d\n",ss.nstar);
	 }
  }
 else
  {	ipointlist	*ipls;
	int		nstar,n;
	convert_candidates(ss.cands,ss.ncand,&ss.stars,&ss.nstar);
	nstar=ss.nstar;
	ipls=(ipointlist *)malloc(sizeof(ipointlist)*nstar);
	for ( n=0 ; n<nstar ; n++ )
	 {	ipls[n].ipoints=ss.cands[n].ipoints;
		ipls[n].nipoint=ss.cands[n].nipoint;
	 }
	
	if ( ss.cf.niter>0 )
	 {	collective_fit_star_single_model_iterative(&img->i,mask,
			ss.stars,nstar,ipls,&ss.sfp,ss.cf.refinelevel,ss.cf.niter-1);
	 }
	else
	 {	collective_fit_star_single_model_blocked(&img->i,mask,
			ss.stars,nstar,ipls,ss.cf.bhsize);
	 }

	cleanup_starlist(&ss.stars,&ss.nstar);
	free(ipls);
  }

 if ( ss.is_determine_psf )
  {	fistar_determine_psf(&img->i,mask,&ss,&pdparam,&tpd);		}

 if ( ss.stars==NULL )
  {	convert_candidates(ss.cands,ss.ncand,&ss.stars,&ss.nstar);	}

 if ( ss.nstar>0 )
  {	indx=(int *)malloc(sizeof(int)*ss.nstar);
	for ( i=0 ; i<ss.nstar ; i++ )	indx[i]=i;
	switch ( sort )
	 {	case 0:index_qsort(indx,ss.nstar,compare_x    ,ss.stars);break;
		case 1:index_qsort(indx,ss.nstar,compare_y    ,ss.stars);break;
		case 2:index_qsort(indx,ss.nstar,compare_peak ,ss.stars);break;
		case 3:index_qsort(indx,ss.nstar,compare_fwhm ,ss.stars);break;
		case 4:index_qsort(indx,ss.nstar,compare_amp  ,ss.stars);break;
		case 5:index_qsort(indx,ss.nstar,compare_flux ,ss.stars);break;
		case 6:index_qsort(indx,ss.nstar,compare_noise,ss.stars);break;
		case 7:index_qsort(indx,ss.nstar,compare_sn   ,ss.stars);break;
	 }
  }
 else	indx=NULL;


 if ( outname == NULL )	fw=stdout;
 else			fw=fopenwrite(outname);
 if ( fw==NULL )	
  {	fprint_error("unable to create output file '%s'",outname);
	return(1); 
  }
 if ( is_comment )
  {	fprintf(fw,"# Created by fistar %s (fi: %s)\n",FITSH_FISTAR_VERSION,FITSH_VERSION);
	fprintf(fw,"# Invoked command:");
	for ( i=0 ; i<argc ; i++ )
	 {	if ( is_any_nasty_char(argv[i]) )
			fprintf(fw," \"%s\"",argv[i]);
		else
			fprintf(fw," %s",argv[i]);
	 }
	fprintf(fw,"\n");
  }

 write_stars(fw,ss.stars,ss.nstar,indx,is_comment,oformat,&mf0);
 fclosewrite(fw);

 if ( indx != NULL )	free(indx);

 if ( psfoutname != NULL && ss.is_determine_psf )
  {	FILE	*fw;
	fw=fopenwrite(psfoutname);
	psf_write_fits(fw,&tpd);
	fclosewrite(fw);
  }

 if ( markimgname != NULL )
  {	FILE	*fw;
	fw=fopenwrite(markimgname);
	if ( fw != NULL )
	 {	draw_mark_stars(&img->i,ss.stars,ss.nstar,mark_symbol,mark_size);
		fits_write(fw,img);
		fclosewrite(fw);
	 }
  }

 if ( areaimgname != NULL )
  {	FILE	*fw;
	fw=fopenwrite(areaimgname);
	if ( fw != NULL )
	 {	int	i,j,ix,iy;
		for ( i=0 ; i<ss.nstar ; i++ )
		 {	if ( ss.stars[i].cand==NULL )	continue;
			for ( j=0 ; j<ss.stars[i].cand->nipoint ; j++ )
			 {	ix=ss.stars[i].cand->ipoints[j].x,
				iy=ss.stars[i].cand->ipoints[j].y;
				img->i.data[iy][ix]=0.0;
			 }
		 } 
		fits_write(fw,img);
		fclosewrite(fw);
	 }
  }

 fits_free(img);

 return(0);
}

/*****************************************************************************/

/* development branch, star finding based on biquad surfaces	*/
/*
 search_star_candidates_biquad(img,mask,&stars,&nstar);
 if ( is_verbose ) fprintf(stderr,"nstar0=%d\n",nstar);
*/
/* to be continued soon...					*/
