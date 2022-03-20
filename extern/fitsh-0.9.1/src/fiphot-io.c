/*****************************************************************************/
/* fiphot-io.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line tool for performing photometry on FITS images:		     */
/* input/output handling.						     */
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
#include "io/iof.h"
#include "io/scanarg.h"
#include "io/tokenize.h"
#include "magnitude.h"
#include "math/poly.h"

#include "basis.h"
#include "common.h"
#include "kernel.h"

#include "fiphot.h"

/*****************************************************************************/

int read_input_star_list(FILE *fr,photstar **rps,int *rnp,colread *col,int zoom)
{
 photstar	*ps;
 int		np,n,i,k,cr0,cra,cda;
 char		*buff,*cmd[64],*id;
 double		x,y,r0,ra,da,bshx,bshy,ref_mag,ref_col,ref_err;
 apgeom		*inaps;
 int		ninap,use_ref;

 ps=NULL;np=0;

 buff=NULL;
 bshx=bshy=0.0;

 if ( zoom<1 )	zoom=1;

 while ( ! feof(fr) )
  {	if ( buff != NULL )	free(buff);
	buff=freadline(fr);
	if ( buff==NULL )
		break;
        else if ( basis_read_comment(buff,&bshx,&bshy) )
		continue;
	remove_newlines_and_comments(buff);
	if ( strlen(buff)==0 )	continue;
	n=tokenize_spaces(buff,cmd,63);
	if ( col->colx>=n || col->coly>=n )	continue;

	sscanf(cmd[col->colx],"%lg",&x);
	sscanf(cmd[col->coly],"%lg",&y);
	if ( ! ( isfinite(x) && isfinite(y) ) )	continue;

	inaps=NULL;
	ninap=0;
	for ( i=0 ; col->colap[i]>=0 ; i+=3 )
	 {	cr0=col->colap[i],
		cra=col->colap[i+1],
		cda=col->colap[i+2];
		if ( cr0>=n )	continue;
		r0=ra=da=0.0;
		k=sscanf(cmd[cr0],"%lg:%lg:%lg",&r0,&ra,&da);
		if ( k<1 )
			continue;
		if ( ! isfinite(r0) || ! isfinite(ra) || ! isfinite(da) )
			continue;
		inaps=(apgeom *)realloc(inaps,sizeof(apgeom)*(ninap+1));
		inaps[ninap].r0=r0;
		if ( cra>=0 && cda>=0 )
		 {	if ( cra<n && sscanf(cmd[cra],"%lg",&ra)==1 && isfinite(ra) )
				inaps[ninap].ra=ra;
			else
				inaps[ninap].ra=0.0;
			if ( cda<n && sscanf(cmd[cda],"%lg",&da)==1 && isfinite(da) )
				inaps[ninap].da=da;
			else
				inaps[ninap].da=0.0;
		 }
		else
		 {	if ( k>=2 )	inaps[ninap].ra=ra;
			else		inaps[ninap].ra=0.0;
			if ( k>=3 )	inaps[ninap].da=da;
			else		inaps[ninap].da=0.0;
		 }

		inaps[ninap].r0 *= (double)zoom;
		inaps[ninap].ra *= (double)zoom;
		inaps[ninap].da *= (double)zoom;

		ninap++;
	 }
		
	if ( 0<=col->colid && col->colid<n )	id=cmd[col->colid];
	else					id=NULL;

	if ( 0<=col->colmag && col->colmag<n )
		use_ref=sscanf(cmd[col->colmag],"%lg",&ref_mag);
	else	
		ref_mag=0.0,use_ref=0;
	if ( 0<=col->colcol && col->colcol<n )
		sscanf(cmd[col->colcol],"%lg",&ref_col);
	else	
		ref_col=0.0;
	if ( 0<=col->colerr && col->colerr<n )
		sscanf(cmd[col->colerr],"%lg",&ref_err);
	else
		ref_err=0.0;

	ps=(photstar *)realloc(ps,sizeof(photstar)*(np+1));

	ps[np].x=(double)zoom*(x-bshx);
	ps[np].y=(double)zoom*(y-bshy);
	ps[np].inaps=inaps;
	ps[np].ninap=ninap;

	ps[np].use_ref=use_ref;
	ps[np].ref_mag=ref_mag;
	ps[np].ref_col=ref_col;
	ps[np].ref_err=ref_err;

	ps[np].fluxes=NULL;
	ps[np].rfflux=NULL;

	if ( id==NULL )	ps[np].id=NULL;
	else		ps[np].id=strdup(id);

	np++;
  };

 *rps=ps,
 *rnp=np;
 return(0);
}

int read_raw_photometry(FILE *fr,photstar **rps,int *rnp)
{
 photstar	*ps;
 photflux	*rfflux;
 int		np,t,n,i,r,flag;
 char		*buff,**cmd,*id;
 double		x,y,bshx,bshy,ref_mag,ref_col;
 apgeom		*inaps;

 ps=NULL;np=0;

 buff=NULL;cmd=NULL;
 bshx=bshy=0.0;
 while ( ! feof(fr) )
  {	if ( cmd != NULL )
	 {	free(cmd);cmd=NULL;	}
	if ( buff != NULL )
	 {	free(buff);buff=NULL;	}
	buff=freadline(fr);
	if ( buff==NULL )
		break;
        else if ( basis_read_comment(buff,&bshx,&bshy) )
		continue;
	remove_newlines_and_comments(buff);
	if ( strlen(buff)==0 )	continue;
	cmd=tokenize_spaces_dyn(buff);
	for ( t=0 ; cmd[t] != NULL ; )	t++;

	if ( t<12 )	continue;
	id=cmd[0];
	if ( sscanf(cmd[1],"%lg",&x)<1  )	continue;
	if ( sscanf(cmd[2],"%lg",&y)<1  )	continue;
 	if ( ! ( isfinite(x) && isfinite(y) ) )	continue;
	if ( sscanf(cmd[3],"%d",&n)<1 )		continue;
	if ( t<6+n*6 )				continue;

	ps=(photstar *)realloc(ps,sizeof(photstar)*(np+1));

	ps[np].x=x-bshx,
	ps[np].y=y-bshy;
	ps[np].n=ps[np].ninap=n;

	if ( sscanf(cmd[4],"%lg",&ref_mag)==1 )
	 {	ps[np].use_ref=1;
		ps[np].ref_mag=ref_mag;
	 }
	else
	 {	ps[np].use_ref=0;
		ps[np].ref_mag=0.0;
	 }
	if ( sscanf(cmd[5],"%lg",&ref_col)==1 )
	 	ps[np].ref_col=ref_col;
	else
		ps[np].ref_col=0.0;

	inaps=(apgeom *)malloc(sizeof(apgeom)*n);
	for ( i=0 ; i<n ; i++ )
	 {	sscanf(cmd[6+i*6+0],"%lg",&inaps[i].r0);
		sscanf(cmd[6+i*6+1],"%lg",&inaps[i].ra);
		sscanf(cmd[6+i*6+2],"%lg",&inaps[i].da);
	 }
	ps[np].inaps=inaps;

	rfflux=(photflux *)malloc(sizeof(photflux)*n);
	for ( i=0 ; i<n ; i++ )
	 {	r=0;
		r+=sscanf(cmd[6+i*6+3],"%lg",&rfflux[i].flux);
		r+=sscanf(cmd[6+i*6+4],"%lg",&rfflux[i].fluxerr);
		r+=sscanf(cmd[6+i*6+5],"%i",&flag);
		if ( r<3 )
		 {	rfflux[i].flag=MASK_BAD;
			rfflux[i].flux=0.0;
			rfflux[i].fluxerr=0.0;
		 }
		else	
			rfflux[i].flag=flag;

		rfflux[i].bgarea=rfflux[i].bgflux=0.0;
		rfflux[i].bgmedian=0.0;rfflux[i].bgsigma=0.0;
		rfflux[i].mag=rfflux[i].magerr=0;
		rfflux[i].rtot=rfflux[i].rbad=rfflux[i].rign=0;
		rfflux[i].atot=rfflux[i].abad=0;
	 }
	ps[np].rfflux=rfflux;
	ps[np].fluxes=NULL;

	if ( id==NULL )	ps[np].id=NULL;	
	else		ps[np].id=strdup(id);

	np++;
  };

 if ( cmd != NULL )	free(cmd);
 if ( buff != NULL )	free(buff);

 *rps=ps,
 *rnp=np;
 return(0);
}

/*****************************************************************************/

int get_id_format_parameters(photstar *ps,int np,int *rid_len,char *idf,char *ndf)
{
 int	i,j,id_len;

 id_len=7;
 for ( i=0 ; i<np ; i++ )
  {	if ( ps[i].id==NULL )	continue;
	j=strlen(ps[i].id);
	if ( j>id_len )	id_len=j;
  }
 sprintf(idf,"%%-%ds ",id_len);
 sprintf(ndf,"%%%dd ",id_len);
 *rid_len=id_len;
 return(0);
}

int write_output_star_list(FILE *fw,photstar *ps,int np)
{
 int	i,id_len;
 char	id_format[8],nd_format[8];

 if ( np==0 )	return(0);
 get_id_format_parameters(ps,np,&id_len,id_format,nd_format);

 for ( i=0 ; i<np ; i++ )
  {
	if ( ps[i].id==NULL )	fprintf(fw,nd_format,i+1);
	else			fprintf(fw,id_format,ps[i].id);
	
	fprintf(fw,"%11.3f %11.3f",ps[i].x+1.0,ps[i].y+1.0);
	if ( ps[i].optimal.r0>0.0 )
	 {	fprintf(fw," %6.3f",ps[i].optimal.r0);
		if ( ps[i].optimal.ra>0.0 )
		 {	fprintf(fw," %6.3f %6.3f",ps[i].optimal.ra,ps[i].optimal.da);	}
	 }
	fprintf(fw,"\n");
  }
 return(0);
}

/*****************************************************************************/

static int photometry_status_letter(photflux *pf)
{
 if ( ( pf->rbad>0 && pf->rbad>pf->rign ) || pf->atot==0 || pf->flux <= 0.0 )
	return('X');
 else if ( pf->rign>0 )
	return('C');
 else
	return('G');
}

int photometry_status(char *buff,photflux *pf)
{
 	   /* "X-bfhcsio-CCC/BBB/TTT-GGG/TTT", length: 25 */
 sprintf(buff,"----------%.3d/%.3d/%.3d-%.3d/%.3d",pf->rign,pf->rbad,pf->rtot,pf->abad,pf->atot);
 if ( pf->flag & MASK_NOBACKGROUND )	buff[2]='b';
 if ( pf->flag & MASK_FAULT )		buff[3]='f';
 if ( pf->flag & MASK_HOT )		buff[4]='h';
 if ( pf->flag & MASK_COSMIC )		buff[5]='c';
 if ( pf->flag & MASK_SATURATED )	buff[6]='s';
 if ( pf->flag & MASK_INTERPOLATED )	buff[7]='i';
 if ( pf->flag & MASK_OUTER )		buff[8]='o';

 buff[0]=photometry_status_letter(pf);

 return(0);
}
int photometry_short_status(char *buff,photflux *pf)
{
 buff[0]=photometry_status_letter(pf);
 buff[1]=0;
 
 return(0);
}

int write_photometry(FILE *fw,photstar *ps,int np,char *ofxy,char *ofph,
	spatialgain *sg,int basistype,char *nanstring,char *serialstring,int sx,int sy)
{
 int		i,j,id_len,serial_len,fchar,ok_magn,ok_flux,ok_back;
 char		id_format[8],nd_format[8],buff[128],*f;
 double		bshx,bshy,gain;
 photflux *	pf;

 if ( np==0 )	return(0);
 get_id_format_parameters(ps,np,&id_len,id_format,nd_format);

 basis_set_shift_coords(basistype,&bshx,&bshy);
 basis_write_comment(fw,bshx,bshy);

 if ( serialstring==NULL )
  {	serial_len=1;
	serialstring="-";
  }
 else
	serial_len=strlen(serialstring);

 if ( is_comment )
  {	fchar='#';
	for ( f=ofxy ; *f && *f != ',' ; f++ )
	 {	switch ( *f )
		 {	case 'X': fprintf(fw,"%c  X Coord. ",fchar);break;
			case 'Y': fprintf(fw,"%c  Y Coord. ",fchar);break;
		        case 'S': fprintf(fw,"%c Serial ",fchar);
				  for ( i=8 ; i<serial_len ; i++ )	fprintf(fw," ");
				  break;
			case 'I': fprintf(fw,"%c  ID ",fchar);
				  for ( i=5 ; i<id_len ; i++ ) fprintf(fw," ");
				  break;
			case '-': fprintf(fw,"%c - ",fchar);break;
		 }
		fchar=' ';
	 }
	for ( j=0 ; j<ps[0].n ; j++ )
	 {	for ( f=ofph ; *f ; f++ )
		 {	switch ( *f )
			 {	case 'S': fprintf(fw,"%cStatus flags                 ",fchar);break;
				case 's': fprintf(fw,"%cS",fchar);break;
				case 'M': fprintf(fw,"%c Magnitude ",fchar);break;
				case 'm': fprintf(fw,"%cMagn.error ",fchar);break;
				case 'B': fprintf(fw,"%cBackground ",fchar);break;
				case 'b': fprintf(fw,"%cBg.scatter ",fchar);break;
				case 'F': 
				case 'A': fprintf(fw,"%c Flux/ADU/ ",fchar);break;
				case 'f': 
				case 'a': fprintf(fw,"%c ADU error ",fchar);break;
				case 'E': fprintf(fw,"%c Flux/elec/",fchar);break;
				case 'e': fprintf(fw,"%cElec.error ",fchar);break;
				case 'X': fprintf(fw,"%cCenter X   ",fchar);break;
				case 'x': fprintf(fw,"%cCntr X err.",fchar);break;
				case 'Y': fprintf(fw,"%cCenter Y   ",fchar);break;
				case 'y': fprintf(fw,"%cCntr Y err.",fchar);break;
				case 'W': fprintf(fw,"%cProfile shp",fchar);break;
				case 'w': fprintf(fw,"%cPr.shp.err ",fchar);break;
				case 'D': fprintf(fw,"%cPr.dev. D  ",fchar);break;
				case 'K': fprintf(fw,"%cPr.dev. K  ",fchar);break;
				case '-': fprintf(fw,"%c - ",fchar);break;
			 }
			fchar=' ';
		 }
	 }
	fprintf(fw,"\n");
  }

 for ( i=0 ; i<np ; i++ )
  {	gain=eval_2d_poly(ps[i].x,ps[i].y,sg->order,sg->coeff,
		0.5*(double)sx,0.5*(double)sy,0.5*(double)sx);
	if ( sg->vmin > 0 && gain < sg->vmin )	gain=sg->vmin;
	else if ( gain < 0.0 )			gain=0.0;

	for ( f=ofxy ; *f && *f != ',' ; f++ )
	 {	switch ( *f )
		 {   case 'X': 
			fprintf(fw,"%11.3f ",ps[i].x+bshx);
			break;
		     case 'Y':
			fprintf(fw,"%11.3f ",ps[i].y+bshy);
			break;
		     case 'S':
			fprintf(fw,"%8s ",serialstring);
			break;
		     case 'I':
			if ( ps[i].id==NULL )	fprintf(fw,nd_format,i+1);
			else			fprintf(fw,id_format,ps[i].id);
			break;
		     case '-': fprintf(fw,"  - ");break;
		 }
	 }
	for ( j=0 ; j<ps[i].n ; j++ )
	 {	
		pf=&ps[i].fluxes[j];

		for ( f=ofph ; *f ; f++ )
		 {	ok_back=!(pf->flag & MASK_NOBACKGROUND);
			ok_flux=ok_back;
			ok_magn=(ok_back && pf->flux>0.0 );
			if ( *f=='-' )
				fprintf(fw,"  - ");
			else if ( *f=='B' || *f=='b' )
			 {	if ( ! ok_back )
					fprintf(fw,"%11s ",nanstring);
				else if ( *f=='B' )
					fprintf(fw,"%11.2f ",pf->bgmedian);
				else	
					fprintf(fw,"%11.2f ",pf->bgsigma);
			 }
			else if ( *f=='S' )
			 {	photometry_status(buff,pf);
				fprintf(fw,"%s ",buff);
			 }
			else if ( *f=='s' )
			 {	photometry_short_status(buff,pf);
				fprintf(fw,"%s ",buff);
			 }
			else if ( pf->rign<pf->rbad && pf->rbad==pf->rtot )
				fprintf(fw,"%11s ",nanstring);	
			else
			 {	switch ( *f )
				 {   case 'M':	
					if ( ok_magn )
						fprintf(fw,"%11.5f ",pf->mag);
					else
						fprintf(fw,"%11s ",nanstring);
					break;
				     case 'm':	
					if ( ok_magn )
						fprintf(fw,"%11.5f ",pf->magerr);
					else
						fprintf(fw,"%11s ",nanstring);
					break;
				     case 'F': case 'A':
					if ( ok_flux )
						fprintf(fw,"%11.2f ",pf->flux);
					else
						fprintf(fw,"%11s ",nanstring);
					break;
				     case 'f': case 'a':
					if ( ok_flux )
						fprintf(fw,"%11.2f ",pf->fluxerr);
					else
						fprintf(fw,"%11s ",nanstring);
					break;
				     case 'X': 
					if ( ok_flux )
						fprintf(fw,"%11.3f ",pf->cntr_x);
					else
						fprintf(fw,"%11s ",nanstring);
					break;
				     case 'x': 
					if ( ok_flux )
						fprintf(fw,"%11.3f ",pf->cntr_x_err);
					else
						fprintf(fw,"%11s ",nanstring);
					break;
				     case 'Y': 
					if ( ok_flux )
						fprintf(fw,"%11.3f ",pf->cntr_y);
					else
						fprintf(fw,"%11s ",nanstring);
					break;
				     case 'y': 
					if ( ok_flux )
						fprintf(fw,"%11.3f ",pf->cntr_y_err);
					else
						fprintf(fw,"%11s ",nanstring);
					break;
				     case 'W': 
					if ( ok_flux && 0.0<pf->cntr_width )
						fprintf(fw,"%11.3f ",pf->cntr_width);
					else
						fprintf(fw,"%11s ",nanstring);
					break;
				     case 'w': 
					if ( ok_flux && 0.0<pf->cntr_width )
						fprintf(fw,"%11.3f ",pf->cntr_w_err);
					else
						fprintf(fw,"%11s ",nanstring);
					break;
				     case 'D': 
					if ( ok_flux && 0.0<pf->cntr_width )
						fprintf(fw,"%11.3f ",pf->cntr_w_d);
					else
						fprintf(fw,"%11s ",nanstring);
					break;
				     case 'K': 
					if ( ok_flux && 0.0<pf->cntr_width )
						fprintf(fw,"%11.3f ",pf->cntr_w_k);
					else
						fprintf(fw,"%11s ",nanstring);
					break;
				     case 'E': 
					if ( ok_flux )
						fprintf(fw,"%11.2f ",pf->flux*gain);
					else
						fprintf(fw,"%11s ",nanstring);
					break;
				     case 'e': 
					if ( ok_flux )
						fprintf(fw,"%11.2f ",pf->fluxerr*gain);
					else
						fprintf(fw,"%11s ",nanstring);
					break;
				     default:
					fprintf(fw,"%11s ",nanstring);
					break;
				 }
			 }
		 }
	 }
	fprintf(fw,"\n");
  }

 return(0); 
}

int write_raw_photometry(FILE *fw,photstar *ps,int np,int basistype,char *nanstring)
{
 int		i,j,id_len,n;
 char		id_format[8],nd_format[8];
 photflux	*pf;
 double		x,y,bshx,bshy;

 if ( np==0 )	return(0);
 get_id_format_parameters(ps,np,&id_len,id_format,nd_format);

 if ( is_comment )
  {	fprintf(fw,"# Raw photometry file, summarize all information which may required by \n");
	fprintf(fw,"# photometry on subtracted images. Format is fixed, cannot be changed.\n");
	fprintf(fw,"# Do not edit unless you exactly know what are you doing!\n");
  }
 basis_set_shift_coords(basistype,&bshx,&bshy);
 basis_write_comment(fw,bshx,bshy);

 for ( i=0 ; i<np ; i++ )
  {	if ( ps[i].id==NULL )	fprintf(fw,nd_format,i+1);
	else			fprintf(fw,id_format,ps[i].id);
	x=ps[i].x+bshx,
	y=ps[i].y+bshy;
	fprintf(fw," %11.3f %11.3f",x,y);

	n=ps[i].n;
	fprintf(fw," %3d",n);

	if ( ps[i].use_ref )	
		fprintf(fw," %10.5f %10.5f",ps[i].ref_mag,ps[i].ref_col);
	else
		fprintf(fw," %10s %10s","-","-");

	for ( j=0 ; j<n ; j++ )
	 {	pf=&ps[i].fluxes[j];
		fprintf(fw," %5.2f %5.2f %5.2f",pf->ag.r0,pf->ag.ra,pf->ag.da);
		if ( pf->flag & 0x0100 ) 
			fprintf(fw," %11g %11g 0x%04x",0.0,0.0,pf->flag);	
		else
			fprintf(fw," %11g %11g 0x%04x",pf->flux,pf->fluxerr,pf->flag);	
	 }
	fprintf(fw,"\n");
  }

 return(0); 
}

/*****************************************************************************/

char *output_format_xy(char *oformat)
{
 char 	*ofxy;
 int	c;
 for ( ofxy=oformat ;*oformat && *oformat != ',' ; oformat++ )
  {	c=*oformat;
	if ( ! ( c=='X' || c=='Y' || c=='I' || c=='S' || c=='-' ) )
		return(NULL);
  }
 return(ofxy);
}
char *output_format_ph(char *oformat)
{
 char 	*ofph,*achr;
 int	c;
 while ( *oformat && *oformat != ',' )	oformat++;
 if ( *oformat != ',' )	return(NULL);
 oformat++;
 achr="MmBbFfEeAaXxYyWwDKSs-";
 for ( ofph=oformat ; *oformat && *oformat != ',' ; oformat++ )
  {	c=*oformat;
	if ( strchr(achr,c)==NULL )
		return(NULL);
  }
 return(ofph);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

char *subtracted_format_xy(char *sformat)
{
 char 	*sfxy;
 int	c;
 for ( sfxy=sformat ;*sformat && *sformat != ',' ; sformat++ )
  {	c=*sformat;
	if ( ! ( c=='N' || c=='H' ) )	return(NULL);
  }
 return(sfxy);
}
char *subtracted_format_ph(char *sformat)
{
 char 	*sfph;
 int	c;
 while ( *sformat && *sformat != ',' )	sformat++;
 if ( *sformat != ',' )	return(NULL);
 sformat++;
 for ( sfph=sformat ;*sformat && *sformat != ',' ; sformat++ )
  {	c=*sformat;
	if ( ! ( c=='M' || c=='m' || c=='F' || c=='f' || c=='E' || c=='e' || c=='A' || c=='a' ) )
		return(NULL);
  }
 return(sfph);
}

/*****************************************************************************/

int * create_col_ap_data(char *apcolpar)
{ 
 int	*colap;
 char	*acp,*cmd[64];
 int	i,n,k,cr0,cra,cda,ncc;

 colap=NULL;ncc=0;

 if ( apcolpar==NULL )
  {	colap=(int *)malloc(sizeof(int)*1);
	colap[0]=-1;
	return(colap);
  }
 else
  {	acp=strdup(apcolpar);
	n=tokenize_char(acp,cmd,',',63);
	for ( i=0 ; i<n ; i++ )
	 {	cr0=cra=cda=-1;
		k=sscanf(cmd[i],"%d:%d:%d",&cr0,&cra,&cda);
		if ( k<1 || cr0<=0 )	continue;
		colap=(int *)realloc(colap,(3*ncc+4)*sizeof(int));
		colap[3*ncc]=cr0-1;
		if ( cra>0 )	colap[3*ncc+1]=cra-1;
		else		colap[3*ncc+1]=-1;
		if ( cda>0 )	colap[3*ncc+2]=cda-1;
		else		colap[3*ncc+2]=-1;
		ncc++;
	 }
	colap[3*ncc]=-1;
	free(acp);
 	return(colap);
  }

}

int create_input_ap_data(char *appar,apgeom **rinaps,int *rninap,int zoom)
{
 apgeom	*inaps;
 int	ninap,n,i,k;
 char	*ap,*cmd[64];
 double	r0,ra,da;

 if ( zoom<1 )	zoom=1;

 if ( appar==NULL )
  {	inaps=NULL;
	ninap=0;
  }
 else
  {	ap=strdup(appar);
	n=tokenize_char(ap,cmd,',',63);
	inaps=NULL;ninap=0;
	for ( i=0 ; i<n ; i++ )
	 {	r0=ra=da=0.0;
		k=sscanf(cmd[i],"%lg:%lg:%lg",&r0,&ra,&da);
		if ( k<1 || r0<=0.0 )	continue;
		inaps=(apgeom *)realloc(inaps,(ninap+1)*sizeof(apgeom));
		inaps[ninap].r0=r0;
		if ( k<2 )   ra=2.0*r0;
		if ( k<3 )   da=ra;
		inaps[ninap].ra=ra;
		inaps[ninap].da=da;
		if ( r0<=0.0 || ra<r0 || da<=0.0 )
		 {	free(inaps);
			return(1);
		 }
		ninap++;
	 }
	free(ap);
  }
 for ( i=0 ; i<ninap ; i++ )
  {	inaps[i].r0 *= (double)zoom;
	inaps[i].ra *= (double)zoom;
	inaps[i].da *= (double)zoom;
  }
 *rinaps=inaps;
 *rninap=ninap;
 return(0);
}

int check_apertures(photstar *ps,int np)
{
 int	i,j,ninap;
 double	r0,ra,da;
 if ( ps==NULL || np<=0 )	return(0);
 ninap=ps[0].ninap;
 for ( i=1 ; i<np ; i++ )
  {	if ( ninap != ps[i].ninap )	return(-1);	}
 if ( ninap==0 )	return(0);
 for ( i=0 ; i<np ; i++ )
  {	for ( j=0 ; j<ninap ; j++ )
	 {	r0=ps[i].inaps[j].r0;
		ra=ps[i].inaps[j].ra;
		da=ps[i].inaps[j].da;
		if ( r0<=0.0 || ra<r0 || da<=0.0 )	return(-1);
	 }
  }
 return(ninap);
}

/*****************************************************************************/
