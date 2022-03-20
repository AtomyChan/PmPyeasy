/*****************************************************************************/
/* grtrans.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line utility to fit general two-dimensional transformations.	     */
/*****************************************************************************/
#define	FITSH_GRTRANS_VERSION	"1.0"
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>

#include "longhelp.h"
#include "fitsh.h"

#include "io/iof.h"
#include "io/scanarg.h"
#include "io/tokenize.h"
#include "math/fit/lmfit.h"
#include "math/poly.h"
#include "math/polyfit.h"
#include "projection.h"

#include "transform.h"
#include "common.h"

#ifdef  HAVE_NO_CC_EXTENSION 
#define __extension__ 
#endif 

typedef struct
 {      double                  ra0,de0,roll0;
	double			qr0,qi0,qj0,qk0;
	double			max_distance;
        int                     type,order;
        matrix                  mproj;
        double                  zfactor;
        projection_distort      distort;
 } projection_data;

typedef struct
 {	projection_data init;
	double  *pixpoly[2];
	double  *prjpoly[2];
	double  crpix1,crpix2;
	double  cd11,cd12,cd21,cd22;
	double  res_prj,res_pix;
 } wcs_data;


/*****************************************************************************/

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



/*****************************************************************************/

typedef struct
 {	int	refx,refy;
	int	nval;
	int	*fitc;
	int	*outc;
	int	is,id,ik,
		os,od,ok;
	int	weight,wgh_is_mag;
	double	wgh_power;
 } colinfo;

typedef struct
 {	double	x,y;
	double	weight,delta;
	double	*vals;
 } fitpoint;

/*****************************************************************************/

int grtrans_wcs_perform_fit(fitpoint *fps,int nfp,projection_data *wi,wcs_data *wcs)
{
 int	i,k,m,nvar,order;
 point  *p_pix;
 point  *p_prj;
 point  *ps;
 double	*vfits[4],cd11,cd12,cd21,cd22,crpix1,crpix2,det,sig2;
 double	ra,de,x,y,z;

 if ( wi==NULL || wcs==NULL )	return(-1);
 memmove(&wcs->init,wi,sizeof(projection_data));
 
 order=wcs->init.order;
 if ( order<=0 )		return(-1);

 nvar=(order+1)*(order+2)/2;

 for ( k=0 ; k<4 ; k++ )
  {	vfits[k]=(double *)malloc(sizeof(double)*nvar);		}

 p_pix=(point *)malloc(sizeof(point)*nfp);
 p_prj=(point *)malloc(sizeof(point)*nfp);
 ps   =(point *)malloc(sizeof(point)*nfp);

 /* WCS fit stage 1: prepare projection: */
 for ( i=0 ; i<nfp ; i++ )
  {	ra=fps[i].x,
	de=fps[i].y;
	/* Here the roll angle is zero by definition, it is handled differently: */
	projection_do_coord(ra,de,wcs->init.ra0,wcs->init.de0,0.0,&x,&y,&z);
	projection_do_distortion(wcs->init.type,NULL,&x,&y,&z);
        p_prj[i].x=x*M_R2D;
        p_prj[i].y=y*M_R2D;
        p_pix[i].x=fps[i].vals[0];
        p_pix[i].y=fps[i].vals[1];
  }

 /* WCS fit stage 2: get an initial hint for CRPIX[12] and CD[12]_[12]: */
 for ( i=0 ; i<nfp ; i++ )
  {	ps[i].x=p_pix[i].x;
	ps[i].y=p_pix[i].y;
	ps[i].weight=1.0;
  }
 for ( i=0 ; i<nfp ; i++ )
  {	ps[i].value=p_prj[i].x;		}
 fit_2d_poly(ps,nfp,1,vfits[0],0,0,1);
 for ( i=0 ; i<nfp ; i++ )
  {	ps[i].value=p_prj[i].y;		}
 fit_2d_poly(ps,nfp,1,vfits[1],0,0,1);
 cd11=vfits[0][1],cd12=vfits[0][2];
 cd21=vfits[1][1],cd22=vfits[1][2];
 det=(cd11*cd22-cd12*cd21);
 crpix1=-(+cd22*vfits[0][0]-cd12*vfits[1][0])/det;
 crpix2=-(-cd21*vfits[0][0]+cd11*vfits[1][0])/det;

 /* WCS fit stage 3: refine parameters by considering higher order distortions: */
 if ( 2 <= order )
  {	/* re-iterate CRPIX1 and CRPIX2: */
	for ( k=0 ; k<3 ; k++ )
	 {	for ( i=0 ; i<nfp ; i++ )
		 {	ps[i].x=p_pix[i].x-crpix1;
			ps[i].y=p_pix[i].y-crpix2;
			ps[i].weight=1.0;
		 }
		for ( i=0 ; i<nfp ; i++ )
		 {	ps[i].value=p_prj[i].x;		}
		fit_2d_poly(ps,nfp,order,vfits[0],0,0,1);
		for ( i=0 ; i<nfp ; i++ )
		 {	ps[i].value=p_prj[i].y;		}
		fit_2d_poly(ps,nfp,order,vfits[1],0,0,1);
		cd11=vfits[0][1],cd12=vfits[0][2];
		cd21=vfits[1][1],cd22=vfits[1][2];
		det=(cd11*cd22-cd12*cd21);
		crpix1 += -(+cd22*vfits[0][0]-cd12*vfits[1][0])/det;
		crpix2 += -(-cd21*vfits[0][0]+cd11*vfits[1][0])/det;
	  }
	
	/* re-iterate CD1_1, CD1_2, CD2_1 and CD2_2: */
	 for ( k=0 ; k<1 ; k++ )
	  {	double	det,cx11,cx12,cx21,cx22;
		det=(cd11*cd22-cd12*cd21);
		for ( i=0 ; i<nfp ; i++ )
		 {	ps[i].x=p_pix[i].x-crpix1;
			ps[i].y=p_pix[i].y-crpix2;
			ps[i].weight=1.0;
		 }
		for ( i=0 ; i<nfp ; i++ )
		 {	ps[i].value=(+cd22*p_prj[i].x-cd12*p_prj[i].y)/det;	}
		fit_2d_poly(ps,nfp,order,vfits[0],0,0,1);
		for ( i=0 ; i<nfp ; i++ )
		 {	ps[i].value=(-cd21*p_prj[i].x+cd11*p_prj[i].y)/det;	}
		fit_2d_poly(ps,nfp,order,vfits[1],0,0,1);
		
		cx11=cd11*vfits[0][1]+cd12*vfits[1][1];
		cx12=cd11*vfits[0][2]+cd12*vfits[1][2];
		cx21=cd21*vfits[0][1]+cd22*vfits[1][1];
		cx22=cd21*vfits[0][2]+cd22*vfits[1][2];
		cd11=cx11;
		cd12=cx12;
		cd21=cx21;
		cd22=cx22;
	  }
   }

 /* force these polynomial terms to be zero and unity, according to the standard: */
 vfits[0][0]=0.0;
 vfits[0][1]=1.0;
 vfits[0][2]=0.0;
 vfits[1][0]=0.0;
 vfits[1][1]=0.0;
 vfits[1][2]=1.0;

 sig2=0.0;
 for ( i=0 ; i<nfp ; i++ )
  {	double	x,y,px,py;
 	px=p_pix[i].x-crpix1;
	py=p_pix[i].y-crpix2;
	x=eval_2d_poly(px,py,order,vfits[0],0,0,1);
	y=eval_2d_poly(px,py,order,vfits[1],0,0,1);
	px=cd11*x+cd12*y;
	py=cd21*x+cd22*y;
	px-=p_prj[i].x;
	py-=p_prj[i].y;
	sig2+=px*px+py*py;
  };
 if ( 0<nfp )	sig2=sqrt(sig2/(double)nfp);
 wcs->res_prj=sig2;

 /* derive the CD^(-1)*(\xi,\eta) => (x-CRPIX1,y-CRPIX2) coefficients */
 det=(cd11*cd22-cd12*cd21);
 for ( i=0 ; i<nfp ; i++ )
  {	ps[i].x=(+cd22*p_prj[i].x-cd12*p_prj[i].y)/det;
	ps[i].y=(-cd21*p_prj[i].x+cd11*p_prj[i].y)/det;	
	ps[i].weight=1.0;
  }
 for ( i=0 ; i<nfp ; i++ )
  {	ps[i].value=p_pix[i].x-crpix1;		}
 fit_2d_poly(ps,nfp,order,vfits[2],0,0,1);
 for ( i=0 ; i<nfp ; i++ )
  {	ps[i].value=p_pix[i].y-crpix2;		}
 fit_2d_poly(ps,nfp,order,vfits[3],0,0,1);
 vfits[2][0]=0.0;
 vfits[3][0]=0.0;

 /* compute residual of the CD^(-1)*(\xi,\eta) => (x-CRPIX1,y-CRPIX2) fit: */
 det=(cd11*cd22-cd12*cd21);
 sig2=0.0;
 for ( i=0 ; i<nfp ; i++ )
  {	double	x,y,px,py;
	px=(+cd22*p_prj[i].x-cd12*p_prj[i].y)/det;
	py=(-cd21*p_prj[i].x+cd11*p_prj[i].y)/det;	
	x=eval_2d_poly(px,py,order,vfits[2],0,0,1);
	y=eval_2d_poly(px,py,order,vfits[3],0,0,1);
	px=x-(p_pix[i].x-crpix1);
	py=y-(p_pix[i].y-crpix2);
	sig2+=px*px+py*py;
  };
 if ( 0<nfp )	sig2=sqrt(sig2/(double)nfp);
 wcs->res_pix=sig2;

 /* save coefficients: */
 wcs->crpix1=crpix1;
 wcs->crpix2=crpix2;
 wcs->cd11=cd11;
 wcs->cd12=cd12;
 wcs->cd21=cd21;
 wcs->cd22=cd22;

 /* re-scale polynomial terms: divide A_ij, B_ij, AP_ij and BP_ij by (i!j!): */
 k=0;
 m=1;
 for ( i=0 ; i<=order ; i++ )
  {	int	j,w;
	w=m;
	for ( j=0 ; j<=i ; j++ )
	 {	vfits[0][k] /= (double)w;
		vfits[1][k] /= (double)w;
		vfits[2][k] /= (double)w;
		vfits[3][k] /= (double)w;
		if ( j<i )	w=(w/(i-j))*(j+1);
		k++;
	 }
	m=m*(i+1);
  }

 /* subtract unity matrix from the linear part of the reverse transformation: */
 vfits[2][1]-=1.0;
 vfits[3][2]-=1.0;

 /* save polynomial terms: */ 
 wcs->pixpoly[0]=vfits[0];
 wcs->pixpoly[1]=vfits[1];
 wcs->prjpoly[0]=vfits[2];
 wcs->prjpoly[1]=vfits[3];

 free(ps);
 free(p_prj);
 free(p_pix);
 
 return(0);
}

static int grtrans_wcs_dump_coeff_headers(FILE *fw,char *type,int o0,int order,double *coeff)
{
 int	i,j,k;

 fprintf(fw,"%s_ORDER=%d,\n",type,order);

 for ( i=0,k=0 ; i<o0 ; i++ )
	k+=(i+1);

 for ( i=o0 ; i<=order ; i++ )
  {	for ( j=0 ; j<=i ; j++,k++ )
	 {	fprintf(fw,"%s_%d_%d=%.12g,\n",type,i-j,j,coeff[k]);	}
  }
 return(0);
}
int grtrans_wcs_dump_headers(FILE *fw,wcs_data *wcs)
{
 char	*ttype;

 fprintf(fw,"WCSAXES=%d,\n",2);

 switch ( wcs->init.type )
  {  case PROJECTION_SIN:
	ttype="SIN";
	break;
     case PROJECTION_TAN:
	ttype="TAN";
	break;
     case PROJECTION_ARC:
	ttype="ARC";
	break;
     default:
	ttype="???";
	break;
  }

 if ( 2<=wcs->init.order )
  {	fprintf(fw,"CTYPE1='RA---%s-SIP',\n",ttype);
	fprintf(fw,"CTYPE2='DEC--%s-SIP',\n",ttype);
  }
 else
  {	fprintf(fw,"CTYPE1='RA---%s',\n",ttype);
	fprintf(fw,"CTYPE2='DEC--%s',\n",ttype);
  }

 fprintf(fw,"CRVAL1=%.4f,\n",wcs->init.ra0);
 fprintf(fw,"CRVAL2=%.4f,\n",wcs->init.de0);
 fprintf(fw,"CD1_1=%.12g,\n",wcs->cd11);
 fprintf(fw,"CD1_2=%.12g,\n",wcs->cd12);
 fprintf(fw,"CD2_1=%.12g,\n",wcs->cd21);
 fprintf(fw,"CD2_2=%.12g,\n",wcs->cd22);
 fprintf(fw,"CRPIX1=%.3f,\n",wcs->crpix1+0.5);
 fprintf(fw,"CRPIX2=%.3f,\n",wcs->crpix2+0.5);

 if ( 2<=wcs->init.order )
  {	grtrans_wcs_dump_coeff_headers(fw,"A" ,2,wcs->init.order,wcs->pixpoly[0]);
	grtrans_wcs_dump_coeff_headers(fw,"B" ,2,wcs->init.order,wcs->pixpoly[1]);
	grtrans_wcs_dump_coeff_headers(fw,"AP",1,wcs->init.order,wcs->prjpoly[0]);
	grtrans_wcs_dump_coeff_headers(fw,"BP",1,wcs->init.order,wcs->prjpoly[1]);
  }

 fprintf(fw,"WCSR_PRJ=%.5g,\n",wcs->res_prj);
 fprintf(fw,"WCSR_PIX=%.5g\n",wcs->res_pix);

 return(0);
}

int grtrans_wcs_data_free(wcs_data *wcs)
{
 if ( wcs->pixpoly[0] != NULL )	free(wcs->pixpoly[0]);
 if ( wcs->pixpoly[1] != NULL )	free(wcs->pixpoly[1]);
 if ( wcs->prjpoly[0] != NULL )	free(wcs->prjpoly[0]);
 if ( wcs->prjpoly[1] != NULL )	free(wcs->prjpoly[1]);
 return(0);
}

/*****************************************************************************/

static char *wrninappr="Warning: inappropriate content in line %d, skipped.\n";

static int grtrans_read_data_points(FILE *fr,colinfo *col,fitpoint **rfps,int *rnfp)
{
 char		*rbuff,**cmd;
 fitpoint	*fps;
 int		nfp,n,k,i,j,ln;
 double		xr,yr,wr,w,*vals;

 fps=NULL;
 nfp=0;

 ln=0;
 while ( ! feof(fr) )
  {	rbuff=freadline(fr);
	if ( rbuff==NULL )	break;
	ln++;
	remove_newlines_and_comments(rbuff);
	cmd=tokenize_spaces_dyn(rbuff);
	if ( cmd==NULL )
	 {	free(rbuff);continue;			}
	for ( n=0 ; cmd[n] != NULL ; )	n++;
	if ( n==0 )
	 {	free(rbuff);free(cmd);continue;		}

	k=0;
	if ( col->refx<n )	k+=sscanf(cmd[col->refx],"%lg",&xr);
	if ( col->refy<n )	k+=sscanf(cmd[col->refy],"%lg",&yr);
	if ( k<2 || ! ( isfinite(xr) && isfinite(yr) ) )
	 {	if ( is_verbose )	fprintf(stderr,wrninappr,ln);
		free(rbuff);free(cmd);continue;
	 }
	if ( 0<=col->weight && col->weight<n )
	 {	k=sscanf(cmd[col->weight],"%lg",&wr);
		if ( k<1 || ! isfinite(wr) )
		 {	if ( is_verbose )	fprintf(stderr,wrninappr,ln);
			free(rbuff);free(cmd);continue;
		 }
		if ( col->wgh_is_mag )
		 {	wr=exp(-0.4*M_LN10*(wr-20.0));			}
		else if ( wr<0.0 )
		 {	if ( is_verbose )	fprintf(stderr,wrninappr,ln);
			free(rbuff);free(cmd);continue;
		 }
		if ( col->wgh_power != 1.0 && wr>0.0 )
		 {	wr=pow(wr,col->wgh_power);			}
	 }
	else	wr=1.0;

	if ( nfp==0 )	fps=(fitpoint *)malloc(sizeof(fitpoint));
	else		fps=(fitpoint *)realloc(fps,sizeof(fitpoint)*(nfp+1));
	vals=(double *)malloc(sizeof(double)*col->nval);

	k=0;
	for ( i=0 ; i<col->nval ; i++ )
	 {	if ( col->fitc[i]<n )
		 {	j=sscanf(cmd[col->fitc[i]],"%lg",&w);
			if ( j && ! isfinite(w) )	j=0;
			k+=j;
			vals[i]=w;
		 }
	 }
	if ( k<col->nval )
	 {	if ( is_verbose )	fprintf(stderr,wrninappr,ln);
		free(vals);
		free(rbuff);free(cmd);continue;
	 }

	fps[nfp].x=xr;
	fps[nfp].y=yr;
	fps[nfp].weight=wr;
	fps[nfp].vals=vals;
	nfp++;
  };
 *rfps=fps;
 *rnfp=nfp;
 return(0);
 
}

/*****************************************************************************/

static int grtrans_parse_number_list(char *p,int *rn,int **rfitc)
{
 int	w,n,*fitc;

 fitc=NULL;n=0;
 while ( *p )
  {	while ( isspace((int)*p) )	p++;
	if ( sscanf(p,"%d",&w)<1 )
	 {	if ( fitc != NULL )	free(fitc);
		return(1);
	 }
	fitc=(int *)realloc(fitc,sizeof(int)*(n+1));
	fitc[n]=w;
	n++;
	while ( isdigit((int)*p) || *p=='-' )	p++;
	while ( isspace((int)*p) || *p==',' )	p++;
  }

 *rn=n;
 *rfitc=fitc;
 return(0);
}

/*****************************************************************************/

int fprint_grtrans_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\tgrtrans [-h|--help|--long-help|--wiki-help] [--version[-short]]\n"
"\t[-C|--comment] [-V|--verbose]\n"
"\t[<in>|-i <in>|--input <in>] [-o|--output <out>]\n");
 fprintf(fw,
"Fit parameters:\n"
"\t[-a|--order order] [-f|--offset <ox>,<oy>] [--scale <scale>]\n"
"\t[-n|--iterations <iter>] [-r|--rejection-level <level-in-sigmas>]\n"
"\t[--output-transformation|-ot <output-transformation>]\n"
"\t[--weight|-w [magnitude],[power=<power>]\n");
 fprintf(fw,
"Transformation parameters:\n"
"\t[--reverse] [-T|--input-transformation <input-transformation>]\n");
 fprintf(fw,
"Input/output format parameters:\n"
"\t[--col-ref|--col-xy <>,<>] [--col-fit <>,...] [--col-out <>,...]\n"
"\t[--col-weight <>] [--col-sdk <>,<>,<> [--col-out-sdk <>,<>,<>]]\n"
"\t[--preserve-comments]\n");
 fprintf(fw,
"Spherical projections:\n"
"\t[--proj [sin|arc|tan],ra=<RA>,dec=<DEC>,[roll=<ROLL>],[max=<max>],\n"
"\t        [degrees|radians|scale=<scale>],\n"
"\t        [k1=<K1>,k2=<K2>,...,p1=<P1>,p2=<P2>,...]\n"
"\t[--proj [...],qr=<quat_r>,qi=<quat_i>,qj=<quat_j>,qk=<quat_k>]\n"
"\t[--col-radec <>,<>|--col-pixel <>,<>]\n");
 fprintf(fw,
"Parameters related to fitting WCS-related FITS keywords:\n"
"\t[--wcs  [sin|arc|tan],order=<distort.poly.order>,ra=<RA>,dec=<DEC>\n"
"\t --col-radec <>,<> --col-pixel <>,<>]\n");

 return(0);
}

longhelp_entry grtrans_long_help[]= 
{
 LONGHELP_OPTIONS,

 { "General administrative options:", NULL },
 { "-h, --help",
	"Give general summary about the command line options." },
 { "--long-help, --help-long",
	"Gives a detailed list of command line options." },
 { "--wiki-help, --help-wiki, --mediawiki-help, --help-mediawiki",
	"Gives a detailed list of command line options in Mediawiki format." },
 { "--version, --version-short, --short-version",
	"Give some version information about the program." },
 { "-C, --comment",
	"Comment the output." },

 { "Options for input/output specification:", NULL },
 { "<inputfile>, -i <inputile>, --input <inputfile>",
	"Name of the input file. If this switch is omitted, the input  "
	"is read from stdin (specifying some input is mandatory)." },
 { "-o <output>, --output <output>, --output-matched <output>",
	"Name  of the output file, if the program was used for transforming "
	"coordinate lists." },
 { "-T, --input-transformation <output-transformation-file>",
	"Name of the  input  transformation  file  (see  also  the  "
	"notes below)." },
 { "--output-transformation <output-transformation-file>",
	"Name of the output file containing the fitted geometrical "
	"transformation, in human-readable format (see also the notes below)." },

 { "In  all  of the above input/output file specifications, the replacement "
	"of the file name by \"-\" (a single minus sign) forces the  "
	"reading  from stdin  or  writing to stdout. Note that all parts "
	"of the any line after \"#\" (hashmark) are treated as a comment, "
	"therefore ignored.", NULL },

 { 	__extension__
	"Note that there is no explicit switch for  distinguishing  between  the "
	"fitting  and  the evaluating purposes of the program. If "
	"--input-transformation has been specified, the program implies "
	"that the  user  wants to  evaluate a function described by this "
	"existing transformation file. On the other hand, if "
	"--output-transformation has been  specified,  the program  fits  "
	"the  parameters  of the function and stores the resulted "
	"transformation file as it specified by the argument of "
	"this option.  In other  words,  if  no  WCS  or spherical "
	"(de)projection declared by the directives of --wcs, one of "
	"these two switches should be given  in  the "
	"command line.", NULL }, { "", NULL },

 { "General options for fitting polynomial coefficients:", NULL },
 { "--col-xy <x>,<y>",
	"Column indices for the independent values. In the current "
	"implementation, grtrans can only fit polynomial functions of "
	"exactly 2  independent  variables. Lines where these columns do "
	"not contain valid real numbers are excluded." },
 { "--col-fit <>[,<>[,<>]...]",
	"Column indices for the dependent values. In the  current  "
	"implementation,  grtrans  can only fit 2 dimensional polynomial "
	"functions to arbitrary dimensional data. The dimension of the "
	"fitted data  is  specified  indirectly, by the number of column "
	"indices specified after this switch." },
 { "-a <order>, --order <order>",
	"Order of the fitted polynomial function. It can be any "
	"positive integer." },
 { "-n <N>, --iterations <N>; -r <S>, --rejection-level <S>",
	"These  switches specify the total number of rejection iterations "
	"of outlyer points and the rejection level  in  sigma  units.  "
	"By default,  no rejection is applied, therefore all valid "
	"lines are used. " },
 { "--col-weight <w>",
	"Column index for optional weights.  If  specified,  this  column "
	"should contain a valid non-negative real number which is used as a "
	"weight during the least-squares fit. " },
 { "--weight [magnitude],[power=<P>]",
	__extension__
	"These directives specify the weights which are used  during  "
	"the fit  of  the functions or transformations. For example, "
	"in practice it is useful in the following situation. We  try  "
	"to  match star  lists,  then the fainter stars are believed to "
	"have higher astrometrical errors, therefore they should have "
	"smaller  influence  in  the  fit.  We can take the weights "
	"from a given colum, specified by the index  after  --col-weight  "
	"(see  above).   The weights  can  be derived from stellar "
	"magnitudes, if so, specify \"magnitude\" to convert the read "
	"values in magnitude to flux. The real  weights  then  is  "
	"the  \"power\"th  power  of the flux. The default value of "
	"the \"power\" is 1,  however,  for  the  maximum-likelihood  "
	"estimation  of an assumed Gaussian distribution, the weights "
	"should be the second power of the fluxes." },
 { "-ot, --output-transformation <output-transformation-file>",
	"Name of the output transformation file containing the result of "
	"the fit (see above)." },

 { "General options for transformation or function evaluation:", NULL },
 { "--input-transformation <input-transformation-file>",
	"Name  of  the  input  transformation file containing the desired "
	"transformation (see above)." },
 { "--reverse, --inverse",
	"Perform inverse transformation. This is a valid option  only  if "
	"the  dimension  of the fitted function is the same as the dimension "
	"of the independent variables (namely,  2,  because  in  the actual "
	"implementation the latter can only be 2)." },
 { "--col-xy <x>,<y>",
	"Column indices for the independent values. In the current "
	"implementation, `grtrans` can only fit polynomial functions of "
	"exactly 2 independent  variables. Lines where these columns "
	"do not contain valid real numbers are excluded." },
 { "--col-out <>[,<>[,<>]...]",
	"Column indices for the evaluated output variables. The number of "
	"indices listed here should be the same as the number of "
	"independent functions stored in the input transformation file." },
 { "General options for transformation composition:", NULL },
 { "--input-transformation <input-transformation-file>",
	"Name of the input transformation file." },
 { "--scale <s>",
	"Scale factor." },
 { "--offset <dx>,<dy>",
	"Shift. The affine transformation with which the input "
	"transformation is composed is going to be: x'=dx+s*x, y'=dy+s*y." },
 { "--output-transformation <output-transformation-file>",
	"Name  of the output transformation file containing the "
	"result of the composition." },

 { "Options for spherical projection and deprojection:", NULL },
 { "--proj, --project [sin|arc|tan],ra=<R>,dec=<D>,[roll=<F>],[max=<max_distance>],[qr=<R>,qi=<I>,qj=<J>,qk=<K>],[degrees|radians|scale=<S>",
	__extension__
	"This set of directives specify  the  common  parameters  of  "
	"the spherical  projection  or  deprojection.  The  \"sin\",  \"arc\" "
	"and \"tan\" directives set the type of projection to orthographic, "
	"arc  and gnomonic,  respectively.  The values after \"ra\" and "
	"\"dec\" (<R> and <D>) specify the center of the  projection  "
	"(right  ascension  and declination, respectively, in degrees). "
	"The optional parameter <F> is for the roll angle (in degrees). "
	"The \"degrees\", \"radians\" or the \"scale=<S>\" directives "
	"specify the scaling of  the  output. The directive \"degrees\" "
	"is equivalent to set \"scale=57.29577951308232087721\" "
	"(180  over \\Pi),  this  is  the default. The directive "
	"\"radians\" is equivalent to set \"scale=1\". " 
	"Alternatively to the right  ascension, declination and roll "
	"parameterization, one can use quaternion-based pointing here "
	"by specifying qr, qi, qj, qk for the real, and the three "
	"imaginary components, respectively." 
	"By default, this projection tries to compute the projected coordinates "
	"for all of the positions presented in the input list. However, it "
	"might yield unexpected results in certain cases (e.g. an input list "
	"on a nearly full sphere and orthographic projection). One can use the "
	"\"max=<max_distance>\" parameter to filter objects "
	"beyond this maximum distance. " },
 { "--col-radec <r>,<d>", 
	"Column indices for RA and DEC values. This option implies "
	"projection." },
 { "--col-pixel <x>,<y>",
	"Column indices for X and Y projected values. This option "
	"implies deprojection." },
 { "--col-out <>,<>",
	"Column indices for the output values (which are X and Y for "
	"projection and RA, DEC for deprojection)." },

 { "Options for fitting WCS information:", NULL },
 { "--wcs [sin|arc|tan],ra=<R>,dec=<D>,order=<order>",
	__extension__ 
	"This set of directives specify the common parameters of WCS "
	"fitting.  The projection can be orthographic, arc or gnomonic "
	"(however,  there  are  dosens of projections implemented in the "
	"FITS WCS system, but for practical purposes, such projections "
	"seem to be  more than enough).  The center of the fit is set "
	"by \"ra\" and \"dec\" (RA and DEC, in degrees). The distortion "
	"order  is  specified by  order.  Note  that  the RA and DEC values "
	"specified here can only be an assumed values, however, the larger "
	"the optical distortions are, the better to have RA and DEC values "
	"close to the optical axis. " },
 { "--col-ref <>,<>",
	"Column indices for the RA and DEC values." },
 { "--col-fit <>,<>",
	"Column indices for the pixel values." },
 { "Note that in this case, the set of the appropriate  FITS  keyword=value "
   "pairs  are  written  directly  stdout, not in the file specified by "
   "the options --output or --output-transformation." }, { "", NULL },

 LONGHELP_EXAMPLES,

 { __extension__ 
   "Here is an example for a complete astrometry problem which demonstrates "
   "the  proper  usage  of the programs grmatch and grtrans.  Let us assume "
   "that we have 1/ a reference star catalogue, named \"catalog.dat\", a "
   "file with four columns: the first is the identifier of the star, the "
   "second and third are the celestial coordinates (RA and DEC, in  degrees), "
   "and the last is the magnitude of the stars; 2/ an  astronomical  image,  "
   "named  \"img.fits\"  (not  crucial for the astrometry itself, it is "
   "required only by the  demonstration  of the export of FITS WCS headers); "
   "and 3/ a list of decected stars (from \"img.fits\"), named \"img.star\", "
   "a file with three columns: the first two are the pixel coordinates and "
   "the third is an estimation of the flux (in ADUs, not in "
   "magnitudes).", NULL }, { "", NULL },
 
 { __extension__
	"Let us also denote the celestial coordinates of the center of the "
	"image by  R  and  D,  the  RA and DEC values, in degrees and, "
	"for example let R=220 and D=25, a field in the Bootes.  "
	"Let us  also  assume  that  the size  of our field (both the "
	"catalog and the list of the deceted stars) is 3 degrees and "
	"there are approximately 4000-4000 stars  both  in  the reference  "
	"catalog  and  in  the list of the detected stars. Because we have "
	"such a large amount of stars, one can use only a fraction of  "
	"them for triangulation. " }, { "", NULL } ,

 { __extension__
	"The  first  step  is to make a projection from the sky, "
	"centered around the center of our image:", NULL }, { "", NULL },

 { "", 	"grtrans --input catalog.dat --proj tan,degrees,ra=220,dec=25 "
	"--col-radec 2,3 --col-out 5,6 --output img.proj"  },

 { "The second step is the point matching:", NULL}, {"", NULL },

 { "", 	"grmatch --reference img.proj --col-ref 5,6 --col-ref-ordering -4 "
	"--input img.star --col-inp 1,2 --col-inp-ordering +3 "
	"--match-points --order 4 "
	"--triangulation auto,unitarity=0.01,maxnumber=1000,conformable "
	"--max-distance 1 --weight reference,column=4,magnitude,power=2 "
	"--comment --output-transformation img.trans" },

 {	__extension__ 
	"This grmach invocation matches the stars from projected reference "
	"catalog, \"img.proj \" and  the  detected  stars. The \"--order 4\" "
	"specifies a fourth-order polynomial fit, which is, in practice, "
	"good  for  a  field with  the  size  of  3 degrees. The directives "
	"after \"--weight\" makes the magnitudes taken from the reference file "
	"to be used  as  a  weight  for fitting.  This  invocation  yields "
	"one new file, \"img.trans\" which stores the fitted 4th-order "
	"polynomial  transformation  which  transforms  the projected "
	"coordinates to the system of the image." , NULL }, { "", NULL },

 {	"The  next step is the astrometrical transformation, we create a "
	"\"local\" catalog, which is the original catalog extended with "
	"the proper X and Y plate coordinates:", NULL }, { "", NULL },

 { "",	"grtrans --input img.proj --col-xy 5,6 "
	"--input-transformation img.trans --col-out 7,8 --output img.cat" },
 
 {	"This  invocation yields an other new file, \"img.cat\" which has 8 "
	"columns. The first six columns are the same as it was in \"img.proj\" "
	"(identifier, RA, DEC, magnitude and projected X, Y coordinates), "
	"the last two colums are the fitted plate coordinates. Then, the "
	"proper WCS headers  can  be determined by the "
	"following call:", NULL }, { "", NULL },

 { "",	"grtrans --input img.cat --col-ref 2,3 --col-fit 7,8 "
	"--wcs tan,order=4,ra=220,dec=25 >img.wcs" },

 {	"The newly created file, img.wcs contains the FITS "
	"\"keyword\"=\"value\" pairs, which can be exported to \"img.fits\" "
	"to have a standard header extended by the  WCS  information.  "
	"For  exporting, the program fiheader(1) can be used:" }, { "", NULL },
 { "",	"fiheader img.fits --rewrite --update \"$(cat img.wcs)\"" },
 
 {	"Note that the last two grtrans  calls  can be replaced by a single "
	"pipeline, when the file img.cat is not created:", NULL}, { "", NULL },
 { "",	"grtrans --input img.proj --col-xy 5,6 --input-transformation "
	"img.trans --col-out 7,8 --output - | grtrans --input - "
	"--col-ref 2,3 --col-fit 7,8 "
	"--wcs tan,order=4,ra=220,dec=25 >img.wcs" },

 { NULL, NULL },

};

int fprint_grtrans_long_help(FILE *fw,int is_wiki)
{
 char	*synopsis=
	"grtransh [options] <input> [-o <output>]";
 char	*description=
	__extension__
	"The main purpose of the program `grtrans` is to transform coordinate "
	"lists and fit transfomrations to input data. The transformation can be "
	"one of the following methods. 1. Evaluate polynomial functions (with "
	"arbitrary order)of two independent values. 2. Two dimensional spherical "
	"projection (converting to RA/DEC or longitude/latitude values to projected "
	"coordinates on a given tangent plane. 3. Two dimensional spherical "
	"de-projection (converting tangent planar coordinates to RA/DEC or "
	"longitude/latitude values). 4. Compose arbitrary polynomial functions of "
	"two independent values with arbitrary two-dimensional affine (or linear) "
	"transformations. "
	"The  program  also  can  be  used  to  fit  functions,  namely  fit the "
	"coefficients of arbitrary-order polynomial functions of two independent "
	"values or fit WCS distortion parameters. ";

 fprint_generic_long_help(fw,is_wiki,grtrans_long_help,synopsis,description);

 return(0);
}

/*****************************************************************************/

#define		PROJ_RADEC_TO_PIXEL	
#define		PROJ_PIXEL_TO_RADEC

int main(int argc,char *argv[])
{
 FILE		*fr,*fw;
 int		i,k,order,nvar,is_help,is_invert,is_ieee;
 char		*infile,*outfile,*outtransfile,*intransfile,
		*colfitstr,*coloutstr,*colsdkstr,*colsdkoutstr,
		*projection_param,*weightparam;
 wcs_data	wcs_data_static,*wcs;
 int		projection_direction;
 colinfo 	col;
 int		is_preserve_comments;
 fitpoint	*fps;
 point		*ps;
 int		nfp,niter,nt,padsize;
 double		rejlevel;
 transformation	tf_data,*tf=&tf_data;
 double		ox,oy,scale; 

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 infile=outfile=intransfile=outtransfile=NULL;
 is_comment=is_verbose=is_help=is_ieee=0;
 colfitstr=coloutstr=colsdkstr=colsdkoutstr=NULL;

 projection_param=NULL,projection_direction=0;
 order=1;

 ox=oy=0.0,scale=1.0;

 col.refx=1,col.refy=2;

 col.nval=2;
 col.fitc=(int *)malloc(sizeof(int)*col.nval);
 col.fitc[0]=3,
 col.fitc[1]=4;

 col.outc=NULL;
 col.is=col.id=col.ik=-1;
 col.os=col.od=col.ok=-1;
 col.weight=0;
 weightparam=NULL;
 col.wgh_is_mag=0;
 col.wgh_power=1;

 niter=0;rejlevel=3.0;
 padsize=12;is_invert=0;

 is_preserve_comments=0;
 
 i=scanarg(argc,argv,SCANARG_ALLOW_FLAGS,
	"--version:%NS-1f%q",&is_help,
	"--version-short|--short-version:%NS-2f%q",&is_help,
	"-h|--help|--short-help|--help-short:%f%q",&is_help,
	"--long-help|--help-long:%SN2f%q",&is_help,
	"--mediawiki-help|--help-mediawiki|--wiki-help|--help-wiki:%SN3f%q",&is_help,
	"-i|--input:%s",&infile,
	"-o|--output:%s",&outfile,
	"-T|--input-transformation:%s",&intransfile,
	"--output-transformation:%s",&outtransfile,
	"-n|--iterations:%d",&niter,
	"-r|--sigma|--rejection|--rejection-level:%g",&rejlevel,
	"-w|--weight:%s",&weightparam,
	"--wcs|--prj|--proj|--project|--projection:%s",&projection_param,
	"-a|--order:%d",&order,
	"-f|--offset:%g,%g",&ox,&oy,
	"--scale:%g",&scale,
	"--colr|--col-ref|--col-xy:%d,%d",&col.refx,&col.refy,
	"--col-radec:" SNf(1) "%d,%d",&projection_direction,&col.refx,&col.refy,
	"--col-pixel:" SNf(2) "%d,%d",&projection_direction,&col.refx,&col.refy,
	"--coli|--col-fit:%s",&colfitstr,
	"--colo|--col-out:%s",&coloutstr,
	"--col-sdk:%s",&colsdkstr,
	"--col-out-sdk|--col-sdk-out:%s",&colsdkoutstr,
	"--col-weight:%d",&col.weight,
	"--preserve-comments:%f",&is_preserve_comments,
	"--padding:%d",&padsize,
	"--reverse|--inverse:%f",&is_invert,
	"--comment:%f",&is_comment,"(C):%f",&is_comment,
	/* "--ieee:%f",&is_ieee, */
	"--verbose:%i",&is_verbose,"(V):%i",&is_verbose,
	"-:%w",&infile,
	"-*|+*:%e",
	"*:%w",&infile,
	NULL);

 if ( i )
  {	fprint_error("invalid command line argument near '%s'",argv[i]);
	return(1);
  }
 else if ( is_help<0 )
  {	fprint_generic_version(stdout,argv[0],"grtrans",FITSH_GRTRANS_VERSION,is_help);
	return(0);
  }
 else if ( 1<is_help )
  {	fprint_grtrans_long_help(stdout,2<is_help);
	return(0);
  }
 else if ( is_help )
  {	fprint_grtrans_usage(stdout);
	return(0);
  }

 if ( intransfile == NULL && outtransfile == NULL && projection_param == NULL )
  {	fprint_error("invalid combination of command line arguments");
	return(1);
  }
 else if ( colfitstr != NULL && coloutstr != NULL )
  {	fprint_error("invalid combination of command line arguments");
	return(1);
  }

 col.refx--,col.refy--;
 if ( col.refx<0 || col.refy<0 )
  {	fprint_error("invalid column indices for coordinates");
	return(1);
  }
 if ( col.weight>0 )	col.weight--;
 else			col.weight=-1;

 if ( col.weight>=0 && weightparam != NULL )
  {	i=scanpar(weightparam,SCANPAR_DEFAULT,
		"magnitude:%f",&col.wgh_is_mag,
		"power:%g",&col.wgh_power,
		NULL);
	if ( i )	
	 {	fprint_error("invalid weighting parameter in '%s'",weightparam);
		return(1);
	 }
  }

 if ( colfitstr != NULL )
  {	free(col.fitc);
	i=grtrans_parse_number_list(colfitstr,&col.nval,&col.fitc);
	if ( i || col.nval<1 )	
	 {	fprint_error("invalid column specifications");
		return(1);
	 }
	for ( i=0 ; i<col.nval ; i++ )
	 {	col.fitc[i]--;
		if ( col.fitc[i]<0 )	
		 {	fprint_error("invalid column specifications");
			return(1);
		 }
	 }
  }
 if ( projection_param != NULL )
  {	int	is_deg,is_rad,w;
	double	radscale,degscale;

	wcs=&wcs_data_static;
	wcs->init.ra0=0.0;
	wcs->init.de0=0.0;
	wcs->init.roll0=0.0;
	wcs->init.max_distance=0.0;
	wcs->init.qr0=0.0;
	wcs->init.qi0=0.0;
	wcs->init.qj0=0.0;
	wcs->init.qk0=0.0;
	wcs->init.type=PROJECTION_TAN;
	wcs->init.order=1;
	wcs->init.zfactor=M_R2D;	/* 180 \over \pi */
	wcs->init.distort.brown_conrady_radial_ncoeff=0;
	wcs->init.distort.brown_conrady_radial_coeffs[0]=0.0;
	wcs->init.distort.brown_conrady_radial_coeffs[1]=0.0;
	wcs->init.distort.brown_conrady_radial_coeffs[2]=0.0;
	wcs->init.distort.brown_conrady_radial_coeffs[3]=0.0;

	wcs->init.distort.brown_conrady_tangential_ncoeff=0;
	wcs->init.distort.brown_conrady_tangential_coeffs[0]=0.0;
	wcs->init.distort.brown_conrady_tangential_coeffs[1]=0.0;
	wcs->init.distort.brown_conrady_tangential_coeffs[2]=0.0;
	wcs->init.distort.brown_conrady_tangential_coeffs[3]=0.0;

	is_deg=is_rad=0;
	radscale=degscale=0.0;

	i=scanpar(projection_param,SCANPAR_DEFAULT,
		"ra:%g" ,&wcs->init.ra0,
		"dec:%g",&wcs->init.de0,
		"roll:%g",&wcs->init.roll0,
		"max:%g",&wcs->init.max_distance,
		"qr:%g",&wcs->init.qr0,
		"qi:%g",&wcs->init.qi0,
		"qj:%g",&wcs->init.qj0,
		"qk:%g",&wcs->init.qk0,
		"tan:" SNf(PROJECTION_TAN),&wcs->init.type,
		"sin:" SNf(PROJECTION_SIN),&wcs->init.type,
		"arc:" SNf(PROJECTION_ARC),&wcs->init.type,
		"order:%d",&wcs->init.order,
		"scale:%g",&radscale,
		"radscale:%g",&radscale,
		"degscale:%g",&degscale,
		"k1:%g",&wcs->init.distort.brown_conrady_radial_coeffs[0],
		"k2:%g",&wcs->init.distort.brown_conrady_radial_coeffs[1],
		"k3:%g",&wcs->init.distort.brown_conrady_radial_coeffs[2],
		"k4:%g",&wcs->init.distort.brown_conrady_radial_coeffs[3],
		"p1:%g",&wcs->init.distort.brown_conrady_tangential_coeffs[0],
		"p2:%g",&wcs->init.distort.brown_conrady_tangential_coeffs[1],
		"p3:%g",&wcs->init.distort.brown_conrady_tangential_coeffs[2],
		"p4:%g",&wcs->init.distort.brown_conrady_tangential_coeffs[3],
		"deg|degree|degrees:%f",&is_deg,
		"rad|radian|radians:%f",&is_rad,
		NULL);
	if ( i )	
	 {	fprint_error("invalid projection paramter in '%s'",projection_param);
		return(1);
	 }
	else if ( radscale<0.0 || degscale<0.0 )
	 {	fprint_error("invalid (negative) scale value");
		return(1);
	 }
	w=0;
	if ( 0.0<radscale )	w++;
	if ( 0.0<degscale )	w++;
	if ( is_rad )		w++;
	if ( is_deg )		w++;
	if ( 2<=w )
	 {	fprint_error("there is a scale ambiguity in the projection type");
		return(1);
	 }

	else if ( 0.0<radscale )	wcs->init.zfactor=radscale;
	else if ( 0.0<degscale )	wcs->init.zfactor=M_R2D*degscale;
	else if ( is_rad )		wcs->init.zfactor=1.0;
	else /* if ( is_deg ) */	wcs->init.zfactor=M_R2D;

	for ( i=3 ; 0<=i ; i-- )
	 {	if ( wcs->init.distort.brown_conrady_radial_coeffs[i] != 0.0 )
		 {	wcs->init.distort.brown_conrady_radial_ncoeff=i+1;
			break;
		 }
	 }

	for ( i=3 ; 0<=i ; i-- )
	 {	if ( wcs->init.distort.brown_conrady_tangential_coeffs[i] != 0.0 )
		 {	wcs->init.distort.brown_conrady_tangential_ncoeff=i+1;
			break;
		 }
	 }
  }
 else
	wcs=NULL;
		
 if ( intransfile != NULL || projection_direction )
  {	if ( coloutstr != NULL )
	 {	i=grtrans_parse_number_list(coloutstr,&col.nval,&col.outc);
		if ( i || col.nval<1 )	
		 {	fprint_error("invalid output column specification");
			return(1);
	 	 }
		for ( i=0 ; i<col.nval ; i++ )
		 {	col.outc[i]--;
			if ( col.outc[i]<0 )
			 {	fprint_error("invalid output column specification");
				return(1);
	 		 }
		 }
	 }
	else 
	 {	col.nval=2;
		col.outc=(int *)malloc(sizeof(int)*2);
		col.outc[0]=col.refx;	
		col.outc[1]=col.refy;
	 }
  }
 if ( colsdkstr != NULL )
  {	i=sscanf(colsdkstr,"%d,%d,%d",&col.is,&col.id,&col.ik);
	if ( i<3 )	
	 {	fprint_error("invalid shape column specification in '%s'",colsdkstr);
		return(1);
	 }
	col.is--,col.id--,col.ik--;
	if ( col.is<0 || col.id<0 || col.ik<0 )
	 {	fprint_error("invalid shape column specification in '%s'",colsdkstr);
		return(1);
	 }
	if ( colsdkoutstr != NULL )
	 {	i=sscanf(colsdkoutstr,"%d,%d,%d",&col.os,&col.od,&col.ok);
		if ( i<3 )	
		 {	fprint_error("invalid output shape column specification in '%s'",colsdkstr);
			return(1);
		 }
		col.os--,col.od--,col.ok--;
		if ( col.os<0 || col.od<0 || col.ok<0 )
		 {	fprint_error("invalid output shape column specification in '%s'",colsdkstr);
			return(1);
		 }
	 }
	else	
		col.os=col.is,
		col.od=col.id,
		col.ok=col.ik;
  }
 else if ( colsdkoutstr != NULL )
  {	fprint_error("invalid combination of command line arguments");
	return(1);
  }
 else
	col.is=col.id=col.ik=-1;
 
 if ( niter<=0 )	niter=0;

 if ( intransfile != NULL && outtransfile != NULL )
  {	double	xlin[3],ylin[3];
	double	*tcoeff;
	fr=fopenread(intransfile);
	if ( fr==NULL )	
	 {	fprint_error("unable to open input transformation file '%s'",intransfile);
		return(1);
	 }
	i=transformation_read_data(fr,tf);
	fcloseread(fr);

	order=tf->order;
	nvar=(order+1)*(order+2)/2;
	xlin[0]=-ox/scale;
	xlin[1]=1.0/scale;
	xlin[2]=0.0;
	ylin[0]=-oy/scale;
	ylin[1]=0.0;
	ylin[2]=1.0/scale;
	tcoeff=(double *)malloc(sizeof(double)*nvar);
	for ( k=0 ; k<tf->nval ; k++ )
	 {	compose_2d_poly_with_affine(tf->vfits[k],order,xlin,ylin,tcoeff);
		for ( i=0 ; i<nvar ; i++ )
		 {	tf->vfits[k][i]=tcoeff[i];	}
	 }
	free(tcoeff);

	fw=fopenwrite(outtransfile);
	if ( fw==NULL )	
	 {	fprint_error("unable to create output transformation file '%s'",outtransfile);
		return(1);
	 }
	transformation_write_data(fw,tf,(is_comment?TRANS_WR_COMMENT:0)|
		TRANS_WR_DXDY|(is_ieee?TRANS_WR_IEEE_64:0));
	fclosewrite(fw);
  }

 else if ( outtransfile != NULL )
  {	double		**vfits;
	int		nval,nrej;
	double		fd,dd,s,sd,sdd,sig,w;
	fitpoint	*wf;

	nval=col.nval;

	if ( infile==NULL )	fr=stdin;
	else			fr=fopenread(infile);
	if ( fr==NULL )	
	 {	fprint_error("unable to open input data file '%s'",infile);
		return(1);
	 }

	grtrans_read_data_points(fr,&col,&fps,&nfp);
	fclose(fr);

	if ( order<0 )		
	 {	fprint_error("invalid fit order");
		return(1);
	 }

	nvar=(order+1)*(order+2)/2;
	vfits=(double **)malloc(sizeof(double *)*nval);
	for ( k=0 ; k<nval ; k++ )
	 {	vfits[k]=(double *)malloc(sizeof(double)*nvar);		}

	ps=(point *)malloc(sizeof(point)*nfp);

	for ( i=0 ; i<nfp ; i++ )
	 {	ps[i].x=fps[i].x,
		ps[i].y=fps[i].y;
		ps[i].weight=fps[i].weight;
	 }

	nrej=1;
	for ( nt=niter ; nt>=0 && nrej ; nt-- )
	 {	int	fit_failed;

		for ( k=0,fit_failed=0 ; k<nval && ! fit_failed ; k++ )
		 {	for ( i=0 ; i<nfp ; i++ )
			 {	ps[i].value=fps[i].vals[k];	  }
			if ( fit_2d_poly(ps,nfp,order,vfits[k],ox,oy,scale) )
				fit_failed=1;
		 }
		if ( fit_failed )
		 {	fprint_error("polynomial fit cannot be obtained.");
			exit(0);
		 }
		
		if ( nt>0 )
		 {	s=sdd=0.0;
			for ( i=0 ; i<nfp ; i++ )
			 {	w=ps[i].weight;
				wf=&fps[i];
				dd=0.0;
				for ( k=0 ; k<nval ; k++ )
				 {	fd=eval_2d_poly(wf->x,wf->y,order,vfits[k],ox,oy,scale)-wf->vals[k];
					dd+=fd*fd;
				 }
				wf->delta=sqrt(dd);
				s+=w;
				sdd+=dd*w;
			 }
			sdd/=s;
			sig=sqrt(sdd);
			nrej=0;
			for ( i=0 ; i<nfp ; i++ )
			 {	wf=&fps[i];
				dd=wf->delta;
				if ( dd>rejlevel*sig )	
				 {	ps[i].weight=0.0;
					nrej++;
				 }
			 }	
		 }
		else
			nrej=0;
	 }	

	if ( outtransfile==NULL )	fw=stdout;
	else				fw=fopenwrite(outtransfile);
	if ( fw==NULL )	
	 {	fprint_error("unable to create output transformation file '%s'",outtransfile);
		return(1);
	 }

	if ( is_comment )
	 {	fprintf(fw,"# Generated by grtrans\n");		}

	tf->type=TRANS_POLYNOMIAL;
	tf->order=order,
	tf->ox=ox,tf->oy=oy,
	tf->scale=scale;
	tf->vfits=vfits;
	tf->nval=nval;

	transformation_write_data(fw,tf,(is_comment?TRANS_WR_COMMENT:0)|
		TRANS_WR_DXDY|(is_ieee?TRANS_WR_IEEE_64:0));

	if ( is_comment )
	 {	s=sd=sdd=0.0;
		for ( i=0 ; i<nfp ; i++ )
		 {	w=ps[i].weight;
			dd=0.0;
			for ( k=0 ; k<nval ; k++ )
			 {	fd=eval_2d_poly(fps[i].x,fps[i].y,order,vfits[k],ox,oy,scale)-fps[i].vals[k];
				dd+=fd*fd;
			 }
			s+=w;
			sdd+=dd*w;
		 }
		sdd/=s;
		sig=sqrt(sdd);
		fprintf(fw,"# Residuals: sigma = %g\n",sig);
	 }
  }

 else if ( intransfile != NULL )
  {	char	*rbuff,**cmd,*nstr,**npnt,*cpnt,outformat[32];
	int	n,nval,mxc;
	double	x,y,nw,is,id,ik,nx,ny;
	double	*jxx,*jxy,*jyx,*jyy;

	fr=fopenread(intransfile);
	if ( fr==NULL )	
	 {	fprint_error("unable to open input transformation file '%s'",intransfile);
		return(1);
	 }
	i=transformation_read_data(fr,tf);
	fcloseread(fr);
	if ( i )	
	 {	fprint_error("unable to parse input transformation file '%s",intransfile);
		return(1);
	 }
	nval=tf->nval;
	if ( col.nval != nval )		
	 {	fprint_error("column/dimension mismatch");
		return(1);
	 }
	if ( is_invert && nval != 2 )
	 {	fprint_error("only 2D -> 2D transformations can be inverted");
		return(1);
	 }

	if ( infile==NULL )	fr=stdin;
	else			fr=fopenread(infile);
	if ( fr==NULL )	
	 {	fprint_error("unable to open input data file '%s'",infile);
		return(1);
	 }

	if ( outfile==NULL )	fw=stdout;
	else			fw=fopenwrite(outfile);
	if ( fw==NULL )	
	 {	fprint_error("unable to create output data file '%s'",outfile);
		return(1);
	 }

	jxx=jxy=jyx=jyy=NULL;
	if ( col.is>=0 )
	 {	if ( nval != 2 )	col.is=col.id=col.ik=-1;
		else
		 {	transformation_get_jacobi(tf,&jxx,&jxy,&jyx,&jyy);	}
	 }
	if ( is_invert && jxx==NULL )
		transformation_get_jacobi(tf,&jxx,&jxy,&jyx,&jyy);

	if ( is_invert && transformation_check_if_null(tf) )
	 {	fprint_error("the transformation is definitely nil, cannot be inverted");
		return(1);
	 }	

	npnt=(char **)malloc((nval+3)*sizeof(char *));
	nstr=(char *)malloc((nval+3)*32);
	for ( k=0 ; k<nval+3 ; k++ )	npnt[k]=nstr+k*32;

	rbuff=NULL;cmd=NULL;
	sprintf(outformat,"%%%ds ",padsize);
	while ( ! feof(fr) )
	 {	if ( rbuff != NULL )	free(rbuff);
		if ( cmd != NULL )	{ free(cmd);cmd=NULL; }
		rbuff=freadline(fr);
		if ( rbuff==NULL )
			break;
		else if ( rbuff[0]=='#' && is_preserve_comments )
		 {	fprintf(fw,"%s",rbuff);
			continue;
		 }
		remove_newlines_and_comments(rbuff);
		cmd=tokenize_spaces_dyn(rbuff);
		if ( cmd==NULL )	continue;
		for ( n=0 ; cmd[n] != NULL ; )	n++;
		if ( n<=col.refx || n<=col.refy )	continue;
		if ( sscanf(cmd[col.refx],"%lg",&x)<1 )	continue;
		if ( sscanf(cmd[col.refy],"%lg",&y)<1 )	continue;
		if ( ! isfinite(x) || ! isfinite(y) )	continue;

		mxc=0;
		if ( is_invert )
		 {	transformation_eval_invert_2d(x,y,tf,&nx,&ny,jxx,jxy,jyx,jyy);
			sprintf(npnt[0],"%15.9g",nx);
			sprintf(npnt[1],"%15.9g",ny);
			if ( col.outc[0]>mxc )	mxc=col.outc[0];
			if ( col.outc[1]>mxc )	mxc=col.outc[1];
		 }
		else
		 {	for ( k=0 ; k<nval ; k++ )
			 {	nw=eval_2d_poly(x,y,tf->order,tf->vfits[k],tf->ox,tf->oy,tf->scale);
				sprintf(npnt[k],"%15.9g",nw);
				if ( col.outc[k]>mxc )	mxc=col.outc[k];
			 }
		 }
		i=0;is=id=ik=0.0;
		if ( col.is>=0 && col.is<n && col.id<n && col.ik<n )
		 {	i+=sscanf(cmd[col.is],"%lg",&is);
			i+=sscanf(cmd[col.id],"%lg",&id);
			i+=sscanf(cmd[col.ik],"%lg",&ik);
		 }
		if ( i==3 )
		 {	double	cjxx,cjxy,cjyx,cjyy,cjdet;
			double	sa,sb,sc,j11,j12,j21,j22,ma,mb,mc,ns,nd,nk;
			sa=is+id,sb=ik,sc=is-id;
			cjxx=eval_2d_poly(x,y,tf->order-1,jxx,tf->ox,tf->oy,tf->scale);
			cjxy=eval_2d_poly(x,y,tf->order-1,jxy,tf->ox,tf->oy,tf->scale);
			cjyx=eval_2d_poly(x,y,tf->order-1,jyx,tf->ox,tf->oy,tf->scale);
			cjyy=eval_2d_poly(x,y,tf->order-1,jyy,tf->ox,tf->oy,tf->scale);
			if ( ! is_invert )
			 {	cjdet=cjxx*cjyy-cjxy*cjyx;
				j11=+cjyy/cjdet;
				j12=+cjxy/cjdet;
				j21=+cjyx/cjdet;
				j22=+cjxx/cjdet;
			 }
			else
			 {	j11=+cjxx;
				j12=-cjxy;
				j21=-cjyx;
				j22=+cjyy;
			 }
			ma=j11*j11*sa+2*j11*j12*sb+j12*j12*sc;
			mb=j21*j11*sa+j21*j12*sb+j22*j11*sb+j22*j12*sc;
			mc=j21*j21*sa+2*j21*j22*sb+j22*j22*sc;
			ns=0.5*(ma+mc);
			nd=0.5*(ma-mc);
			nk=mb;
			sprintf(npnt[nval+0],"%12g",ns);
			sprintf(npnt[nval+1],"%12g",nd);
			sprintf(npnt[nval+2],"%12g",nk);
			if ( col.os>mxc )	mxc=col.os;
			if ( col.od>mxc )	mxc=col.od;
			if ( col.ok>mxc )	mxc=col.ok;
			/*fprintf(stderr,"mxc=%d\n",mxc);*/
		 }

		mxc++;

		for ( i=0 ; i<n || i<mxc ; i++ )
		 {	for ( k=0 ; k<nval ; k++ )
			 {	if ( i==col.outc[k] )	break;	}
			if ( k<nval )		cpnt=npnt[k];
			else if ( i==col.os )	cpnt=npnt[nval+0];
			else if ( i==col.od )	cpnt=npnt[nval+1];
			else if ( i==col.ok )	cpnt=npnt[nval+2];
			else if ( i<n )		cpnt=cmd[i];
			else			cpnt="-";
			fprintf(fw,outformat,cpnt);
		 }
		fprintf(fw,"\n");
		fflush(fw);
	 };
	if ( cmd != NULL )	free(cmd);
	if ( rbuff != NULL )	free(rbuff);
	free(nstr);
	free(npnt);
	if ( jyy != NULL )	free(jyy);
	if ( jyx != NULL )	free(jyx);
	if ( jxy != NULL )	free(jxy);
	if ( jxx != NULL )	free(jxx);

	fclosewrite(fw);
	fcloseread(fr);
  }

 else if ( wcs != NULL )
  {	
	if ( projection_direction )	/* wcs: ra,dec -> x,y or x,y -> ra,dec */
	 {	double	ra,dec,zmin;
		char	*rbuff,**cmd,*nstr,**npnt,*cpnt,outformat[32];
		int	n,nval,mxc;
		double	x,y,z;
	
		if ( wcs->init.qr0 != 0.0 || wcs->init.qi0 != 0.0 || wcs->init.qj0 != 0.0 || wcs->init.qk0 != 0.0 )
			projection_get_matrix_quat(wcs->init.qr0,wcs->init.qi0,wcs->init.qj0,wcs->init.qk0,wcs->init.mproj);
		else
			projection_get_matrix_rdr(wcs->init.ra0,wcs->init.de0,wcs->init.roll0,wcs->init.mproj);

		if ( infile==NULL )	fr=stdin;
		else			fr=fopenread(infile);
		if ( fr==NULL )		
		 {	fprint_error("unable to open input data file '%s'",infile);
			return(1);
		 }
		
		nval=2;
		if ( col.nval != nval )	
		 {	fprint_error("column number must be 2 for WCS information extraction");
			return(1);
		 }

		if ( outfile==NULL )	fw=stdout;
		else			fw=fopenwrite(outfile);
		if ( fw==NULL )		
		 {	fprint_error("unable to create output file '%s'",outfile);
			return(1);
		 }

		npnt=(char **)malloc(nval*sizeof(char *));
		nstr=(char *)malloc(nval*32);
		for ( k=0 ; k<nval ; k++ )	npnt[k]=nstr+k*32;

		rbuff=NULL;cmd=NULL;
		sprintf(outformat,"%%%ds ",padsize);

		if ( 0.0<wcs->init.max_distance )
			zmin=-cos(wcs->init.max_distance*M_PI/180.0);
		else
			zmin=1.1;

		while ( ! feof(fr) )
		 {	
			if ( rbuff != NULL )	free(rbuff);
			if ( cmd != NULL )	{ free(cmd);cmd=NULL; }

			rbuff=freadline(fr);
			if ( rbuff==NULL )	break;
			remove_newlines_and_comments(rbuff);
			cmd=tokenize_spaces_dyn(rbuff);

			if ( cmd==NULL )	continue;
			for ( n=0 ; cmd[n] != NULL ; )	n++;
			if ( n<=col.refx || n<=col.refy )	continue;
			if ( sscanf(cmd[col.refx],"%lg",&x)<1 )	continue;
			if ( sscanf(cmd[col.refy],"%lg",&y)<1 )	continue;
			if ( ! isfinite(x) || ! isfinite(y) )	continue;

			if ( ! ( projection_direction-1 ) )	/* ra,dec -> x,y */
			 {	ra=x,dec=y;
				projection_do_matrix_coord(wcs->init.mproj,ra,dec,&x,&y,&z);
				if ( zmin < z )
					continue;
				projection_do_distortion(wcs->init.type,&wcs->init.distort,&x,&y,&z);
				/* note the '-' sign here: */
				sprintf(npnt[0],"%15.9g",-x*wcs->init.zfactor);
				sprintf(npnt[1],"%15.9g",+y*wcs->init.zfactor);
			 }
			else				/* x,y -> ra,dec */
			 {	/* note the '-' sign here: */
				x=-x/wcs->init.zfactor;
				y=+y/wcs->init.zfactor;
				z=0.0;
				ra=dec=0.0;
				projection_do_inverse_distortion(wcs->init.type,&wcs->init.distort,&x,&y,&z);
				projection_do_inverse_matrix_coord(wcs->init.mproj,x,y,&ra,&dec);
				sprintf(npnt[0],"%15.9g",ra);
				sprintf(npnt[1],"%15.9g",dec);
			 }
			mxc=0;
			if ( col.outc[0]>mxc )	mxc=col.outc[0];
			if ( col.outc[1]>mxc )	mxc=col.outc[1];

			mxc++;
			for ( i=0 ; i<n || i<mxc ; i++ )
			 {	for ( k=0 ; k<nval ; k++ )
				 {	if ( i==col.outc[k] )	break;	}
				if ( k<nval )		cpnt=npnt[k];
				else if ( i<n )		cpnt=cmd[i];
				else			cpnt="-";
				fprintf(fw,outformat,cpnt);
			 }
			fprintf(fw,"\n");
			fflush(fw);
		 };
		if ( cmd != NULL )	free(cmd);
		if ( rbuff != NULL )	free(rbuff);
		free(nstr);
		free(npnt);
		fclosewrite(fw);
		fcloseread(fr);
	 }
	else			/* fit wcs */
	 {     	if ( col.nval != 2 )	
		 {	fprint_error("exactly two input columns are required for WCS determination");
			return(1);
		 }
		if ( infile==NULL )	fr=stdin;
		else			fr=fopenread(infile);
		if ( fr==NULL )		
		 {	fprint_error("unable to open input data file '%s'",infile);
			return(1);
		 }
		grtrans_read_data_points(fr,&col,&fps,&nfp);
		fclose(fr);
		grtrans_wcs_perform_fit(fps,nfp,&wcs->init,wcs);
		fw=stdout;
		grtrans_wcs_dump_headers(fw,wcs);
		grtrans_wcs_data_free(wcs);
	 }
  }

 return(0);
}

