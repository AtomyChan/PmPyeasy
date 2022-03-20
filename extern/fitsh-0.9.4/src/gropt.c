/*****************************************************************************/
/* gropt.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line tool to perform calculations related to geometrical optics.  */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (C) 2016, Pal, Andras <apal@szofi.net>				     */
/*****************************************************************************/
#define FITSH_GROPT_VERSION	"0.1.0"
/*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <sys/time.h>

#include "math/matrixvector.h"
#include "math/polygon.h"
#include "io/tokenize.h"
#include "longhelp.h"
#include "list.h"
#include "optcalc.h"
#include "fitsh.h"

typedef	struct	triangle triangle;

typedef struct 
 {	double	x;
	double	y;
 } triangle_point;

struct triangle
 {	triangle	*prev,*next;
	triangle_point	t_vertex[3];
	int		t_flag;
	double		t_area,t_proj;
 };

int fprint_gropt_usage(FILE *fw,char *argv0)
{
 fprintf(fw,"Usage:\t%s [-h|--help|--long-help|--wiki-help] [--version[-short]]\n",argv0);
 fprintf(fw,"\t[[-i|--input] <input.opt>]\n");
 fprintf(fw,"Testing and benchmarking:\n");
 fprintf(fw,"\t[--speed-test <time/seconds>]\n");
 fprintf(fw,"Spot diagram analysis:\n");
 fprintf(fw,"\t[-s|--spot-aperture <radius>,<number-of-rings>[,<position>]]\n");
 fprintf(fw,"\t[-l|--lambda|--wavelength <wavelength/microns>]\n");
 fprintf(fw,"\t[-f|--focus <focal-plane-position> [-x|--scale <pixel-scale>]]\n");
 fprintf(fw,"\t[-a|--angle <incident-angle/radians>|<normal_x>,<normal_y>]\n");
 fprintf(fw,"\t[-o|--output|--output-spot <spot-diagram-output>]\n");
 fprintf(fw,"Ray transfer matrix analysis:\n");
 fprintf(fw,"\t[-t|--transfer]\n");
 fprintf(fw,"Exporting geometry:\n");
 fprintf(fw,"\t[-d|--output-scad <openscad-file>]\n");
 fprintf(fw,"\t[-d|--output-eps <encapsulated-postscript-file>]\n");
 fprintf(fw,"Determination of point-spread function:\n");
 fprintf(fw,"\t[-z|--psf-half-size <psf-half-size>\n");
 fprintf(fw,"\t[-p|--output-psf <PSF-file>\n");
 return(0);
}

longhelp_entry gropt_long_help[] = 
{
 LONGHELP_OPTIONS,

 { "General options:", NULL },
 { "-h, --help",
	"Give general summary about the command line options." },
 { "--long-help, --help-long",
        "Gives a detailed list of command line options." },
 { "--wiki-help, --help-wiki, --mediawiki-help, --help-mediawiki",
        "Gives a detailed list of command line options in Mediawiki format." },
 { "--version, --version-short, --short-version",
	"Give some version information about the program." },

 { "<input>, -i <input>, --input <input>",
	"Name of the input file describing the optical system. Note that all "
	"of the length dimensions, including offsets, curvatures, curvature "
	"radii and higher order aspherical constants are needed to be defined "
	"in the units of millimeters. "
	"Reading from standard input can be forced using a single dash \"-\" "
	"as input file name. " },

 { "Testing and benchmarking:", NULL },
 { "--speed-test <time/seconds>",
	"In this mode, gropt runs a speed test to figure out the number of "
	"rays that can be traced during the given interval of time. Note that " 
	"this speed depends on the optical system itself, so one should "
	"specify a valid input (see -i, --input)." },

 { "Common options for the various analysis modes:", NULL },
 { "-l, --lambda, --wavelength <wavelength/microns>", 
	"The wavelength (in microns) of the ray set of the actual analysis. " },
 { "-f, --focus <focal-plane-position>",
	"The position of the focal plane. It overrides the \"focal\" keyword "
	"in the input optical system description (see -i, --input for more). " },
 { "-x, --scale <pixel-scale>",
	"The pixel scale (i.e. the pixel size), also in millimeters. " },
 { "-a, --angle <incident-angle/radians>|<normal_x>,<normal_y>",
	"The incident angle of the incoming parallel rays. If a single parameter "
	"is given, then it is going to be interpreted as an angle in radians. "
	"If two parameters are given after -a or --angle then these are treated "
	"as the x and y components of the ray normal vector. By conventions, "
	"a single paramater of \"a\" is equivalent to \"0,sin(a)\". " },

 { "Spot diagram analysis:", NULL },
 { "-s, --spot-aperture <radius>[,<number-of-rings>[,<position>]]", 
	"The radius, the number of the rings and the offset of spot aperture. "
	"If omitted, the number of the rings is going to be one and the "
	"offset of the spot aperture is going to be zero. Note that "
	"the option -o, --output or --output-spot is needed to be specified "
	"in order to perform a spot diagram analysis. " },
 { "-o, --output, --output-spot <spot-diagram-output>",
	"The name of the file to which the spot diagram output is intended to "
	"be written. To write the output to the standard output, use a dash (-) "
	"as a file name. To write the output to the file named as a dash, "
	"use -o ./- or something equivalent. " }, 

 { "Ray transfer matrix analysis: ", NULL },
 { "-t, --transfer",
	"Perform a ray transfer matrix analysis. The results (namely, the "
	"computed focal plane offset and the effective focal length) are "
	"written to the standard output. " },

 { "Exporting geometry:", NULL },
 { "-d, --output-scad <openscad-file>",
	"Exports the geometry of the input optical setup to an "
	"OpenSCAD file. This 3D model of the setup can then be viewed by running "
	"`openscad` whose input is the output this process. " },
 { "-e, --output-eps <encapsulated-postscript-file>",
	"Exports the geometry of the input optical setup to an "
	"encapsulated PostScript file. This PostScript image can be considered "
	"as a planar diagram of the lens system. If -s or --spot-aperture "
	"is defined, then rays are also drawn according to the specifications "
	"following these command-line switches. " }, 

 { "Determination of point-spread function:", NULL },
 { "-z, --psf-half-size <half-size>",
	"The half size of the output point-spread function. If this parameter "
	"is, for example, 3, then the resolution of the output point-spread "
	"function will be 7x7. In general, if the parameter is H, then the "
	"output point-spread function will have a resolution of (2H+1)x(2H+1). "
	"Note that the pixel scale of this computed point-spread function "
	"is defined by -x or --scale. " },

 { NULL, NULL }
};

int fprint_gropt_long_help(FILE *fw,int is_wiki)
{ 
 char	*synopsis=
	"gropt [options] <input> [...] [-o <output>]";
 char	*description=
	__extension__
	"The purpose of the task `gropt` is to perform various computations "
	"in the framework of geometrical optics. The main input of the "
	"task is a descriptor file quantifying the various optical elements "
	"(lenses, mirrors, aperture stops, etc. and the respective materials) "
	"as well as the alignment geometry of these pieces. The output of "
	"`gropt` includes well-known products of optical analysis such as "
	"ray transfer matrices, analysis for multiple wavelengths, "
	"ray diagrams, 3D models for visualization, "
	"spot diagrams, models for spatially quantized point-spread "
	"functions, plate solutions and vignetting. The underlying ray trace "
	"library (see ./src/optcalc.[ch]) is also included in the program "
	"`firandom` in order to generate both precise and accurate simulated "
	"images that can be acquired with the given optical setup. ";

 fprint_generic_long_help(fw,is_wiki,gropt_long_help,synopsis,description);

 return(0);
}



int gropt_test(void)
{
/*
 if ( 1 )
  {	vector	x0,n0;
	x0[0]=0;
	x0[1]=-8;
	x0[2]=0;
	n0[0]=0;
	n0[1]=0.7;
	n0[2]=0.7;
	optcalc_surface_aspheric_ray_trace(0.1,0.0,0,NULL,x0,n0,x0,n0);
	return(0);
  }
*/
/*
 if ( 1 )
  {	double	cc,z0;
	matrix	A;
	vector	B,x0,n0;
	double	C;

	z0=3.0;
	cc=0.01;

	A[0][0]=A[1][1]=A[2][2]=cc;
	A[0][1]=A[1][2]=A[2][0]=0.0;
	A[1][0]=A[2][1]=A[0][2]=0.0;
	B[0]=0;
	B[1]=0;
	B[2]=-1-z0*cc;
	C=z0*(1+0.5*z0*cc);

	x0[0]=0;x0[1]=60;x0[2]=0;n0[0]=0;n0[1]=0;n0[2]=1;
	optcalc_refract_aspheric_ray_trace(cc,0.0,0,NULL,z0,x0,n0,1.0,2.0);
	fprintf(stderr,"%12g %12g %12g\n",x0[0],x0[1],x0[2]);
	fprintf(stderr,"%12g %12g %12g\n",n0[0],n0[1],n0[2]);

	x0[0]=0;x0[1]=60;x0[2]=0;n0[0]=0;n0[1]=0;n0[2]=1;
	optcalc_refract_quadratic_ray_trace(A,B,C,x0,n0,1.0,2.0);
	fprintf(stderr,"%12g %12g %12g\n",x0[0],x0[1],x0[2]);
	fprintf(stderr,"%12g %12g %12g\n",n0[0],n0[1],n0[2]);
	
	return(0);
  }
*/
 return(0);
}

int main(int argc,char *argv[])
{
 int		i;
 char		*infile;
 FILE		*fr;
 struct optics	opt;
 double		lambda,zfoc,zstart,pixel_scale,test_time;
 double		angle,angle_nx,angle_ny,angle_nz;
 int		is_zfoc_set,psf_hsize,test_mode;
 double		aperture_radius;
 int		aperture_nring;
 int		is_calc_transfer;
 char		*output_spot;
 char		*output_scad,*output_eps,*output_psf;

 infile=NULL;
 lambda=0.6;
 aperture_radius=0.0;
 aperture_nring=1;
 angle=angle_nx=angle_ny=0.0;
 zstart=0.0;
 pixel_scale=1.0;
 psf_hsize=1;

 is_calc_transfer=0;

 output_scad=NULL;
 output_eps=NULL;
 output_spot=NULL;
 output_psf=NULL;

 zfoc=0.0;
 is_zfoc_set=0;

 test_mode=0;
 test_time=1.0;

 /* gropt_test(); */

 for ( i=1 ; i<argc; i++ )
  {	if ( strcmp(argv[i],"-h")==0 || strcmp(argv[i],"--help")==0 )
	 {	fprint_gropt_usage(stdout,argv[0]);
		return(0);
	 }
	else if ( strcmp(argv[i],"--long-help")==0 || strcmp(argv[i],"--help-long")==0 )
	 {	fprint_gropt_long_help(stdout,0);
		return(0);
	 }
	else if ( strcmp(argv[i],"--wiki-help")==0 || strcmp(argv[i],"--help-wiki")==0 || strcmp(argv[i],"--mediawiki-help")==0 || strcmp(argv[i],"--help-mediawiki")==0 )
	 {	fprint_gropt_long_help(stdout,1);
		return(0);
	 }
	else if ( strcmp(argv[i],"--version")==0 )
	 {	fprint_generic_version(stdout,argv[0],"gropt",FITSH_GROPT_VERSION,-1);
		return(0);
	 }
	else if ( strcmp(argv[i],"--version-short")==0 || strcmp(argv[i],"--short-version")==0 )
	 {	fprint_generic_version(stdout,argv[0],"gropt",FITSH_GROPT_VERSION,-2);
		return(0);
	 }
	else if ( ( strcmp(argv[i],"-i")==0 || strcmp(argv[i],"--input")==0 ) && i<argc-1 )
		i++,infile=argv[i];
	else if ( strcmp(argv[i],"--speed-test")==0 && i<argc-1 && sscanf(argv[i+1],"%lg",&test_time)==1 )
		i++,test_mode=1;
	else if ( ( strcmp(argv[i],"-o")==0 || strcmp(argv[i],"--output")==0 || strcmp(argv[i],"--output-spot")==0 ) && i<argc-1 )
		i++,output_spot=argv[i];
	else if ( ( strcmp(argv[i],"-d")==0 || strcmp(argv[i],"--output-scad")==0 ) && i<argc-1 )
		i++,output_scad=argv[i];
	else if ( ( strcmp(argv[i],"-e")==0 || strcmp(argv[i],"--output-eps")==0 ) && i<argc-1 )
		i++,output_eps=argv[i];
	else if ( ( strcmp(argv[i],"-p")==0 || strcmp(argv[i],"--output-psf")==0 ) && i<argc-1 )
		i++,output_psf=argv[i];
	else if ( ( strcmp(argv[i],"-z")==0 || strcmp(argv[i],"--psf-half-size")==0 ) && i<argc-1 && sscanf(argv[i+1],"%d",&psf_hsize)==1 )
		i++;
	else if ( ( strcmp(argv[i],"-s")==0 || strcmp(argv[i],"--spot-aperture")==0 ) && i<argc-1 && sscanf(argv[i+1],"%lg,%d,%lg",&aperture_radius,&aperture_nring,&zstart) )
		i++;
	else if ( ( strcmp(argv[i],"-l")==0 || strcmp(argv[i],"--lambda")==0 || strcmp(argv[i],"--wavelength")==0) && i<argc-1 && sscanf(argv[i+1],"%lg",&lambda)==1 )
		i++;
	else if ( ( strcmp(argv[i],"-x")==0 || strcmp(argv[i],"--scale")==0 || strcmp(argv[i],"--wavelength")==0) && i<argc-1 && sscanf(argv[i+1],"%lg",&pixel_scale)==1 )
		i++;
	else if ( ( strcmp(argv[i],"-f")==0 || strcmp(argv[i],"--focus")==0 ) && i<argc-1 && sscanf(argv[i+1],"%lg",&zfoc)==1 )
		i++,is_zfoc_set=1;
	else if ( ( strcmp(argv[i],"-a")==0 || strcmp(argv[i],"--angle")==0 ) && i<argc-1 && sscanf(argv[i+1],"%lg,%lg",&angle_nx,&angle_ny)==2 )
		i++;
	else if ( ( strcmp(argv[i],"-a")==0 || strcmp(argv[i],"--angle")==0 ) && i<argc-1 && sscanf(argv[i+1],"%lg",&angle)==1 )
		i++;
	else if ( strcmp(argv[i],"-t")==0 || strcmp(argv[i],"--transfer")==0 )
		is_calc_transfer=1;
	else if ( strcmp(argv[i],"-")==0 )
		infile=argv[i];
	else if ( argv[i][0]=='-' )
	 {	fprintf(stderr,"%s: error: invalid command line argument near '%s'.\n",argv[0],argv[i]);
		return(1);
	 }
	else
		infile=argv[i];
  }

 if ( 0.0<angle )
  {	angle_nx=0.0;
	angle_ny=sin(angle);
	angle_nz=cos(angle);
  }
 else
  {	double	w;
	w=1.0-(angle_nx*angle_nx+angle_ny*angle_ny);
	if ( w<0.0 )
	 {	fprintf(stderr,"%s: error: unexpected (i.e. too large) incident angle vector.\n",argv[0]);
		return(1);
	 }
	angle_nz=sqrt(w);
  }

 if ( infile==NULL || strcmp(infile,"-")==0 )
	fr=stdin;
 else if ( (fr=fopen(infile,"rb"))==NULL )
  {	fprintf(stderr,"%s: error: unable to open input file '%s'.\n",argv[0],infile);
	return(1);
  }

 optcalc_reset(&opt);

 if ( 0<(i=optcalc_read_optics(fr,&opt)) )
  {	fprintf(stderr,"%s: error: unable to parse input file (line %d).\n",argv[0],i);
	fclose(fr);
	return(1);
  }

 if ( is_zfoc_set )
	opt.opt_z_focal=zfoc;

 fclose(fr);

 if ( test_mode==1 )
  {	vector	v0,n0;
	struct	timeval	tv0,tv1;
	int	ntot;
	v0[0]=0.0;
	v0[1]=0.0;
	v0[2]=zstart-0.0001;
	n0[0]=angle_nx;
	n0[1]=angle_ny;
	n0[2]=angle_nz;

	gettimeofday(&tv0,NULL);
	ntot=0;
	optcalc_glass_refraction_precompute(&opt,lambda);
	while ( 1 )
	 {	int	cnt;
		double	dt;
		for ( cnt=0; cnt<1000; cnt++ )
		 {	optcalc_ray_trace(&opt,lambda,v0,n0,NULL);
		 }
		ntot+=cnt;
		gettimeofday(&tv1,NULL);
		tv1.tv_usec -= tv0.tv_usec;
		tv1.tv_sec  -= tv0.tv_sec;
		if ( tv1.tv_usec<0 )
	 	 {	tv1.tv_usec+=1000000;
			tv1.tv_sec--;
		 }
		dt=(double)tv1.tv_sec+(double)tv1.tv_usec/1000000.0;
		if ( test_time <= dt )
			break;
	 }
	fprintf(stderr,"# %d\n",ntot);
  }

 if ( output_spot != NULL && 0<aperture_nring && 0.0<aperture_radius )
  {	int	i,j,nsec;
	double	sx0,sy0;
	FILE	*fw;
	vector	v0,n0;

	if ( strcmp(output_spot,"-")==0 )
 		fw=stdout;
	else if ( (fw=fopen(output_spot,"wb"))==NULL )
	 {	fprintf(stderr,"%s: error: unable to open or create output spot file '%s': %s\n",argv[0],output_spot,strerror(errno));
		return(1);
	 }

	optcalc_glass_refraction_precompute(&opt,lambda);

	v0[0]=0.0;
	v0[1]=0.0;
	v0[2]=zstart-0.0001;
	n0[0]=angle_nx;
	n0[1]=angle_ny;
	n0[2]=angle_nz;
	optcalc_ray_trace(&opt,lambda,v0,n0,NULL);
	sx0=v0[0];
	sy0=v0[1];

	for ( i=1 ; i<=aperture_nring ; i++ )
	 {	if ( i==0 )	nsec=1;
		else		nsec=i*6;
		for ( j=0 ; j<nsec ; j++ )
		 {	double	x,y,r,w;
			r=aperture_radius*(double)i/(double)aperture_nring;
			w=2*M_PI*j/(double)nsec;
			x=r*cos(w);
			y=r*sin(w);
 			v0[0]=x;
			v0[1]=y;
			v0[2]=zstart-0.0001;
			n0[0]=angle_nx;
			n0[1]=angle_ny;
			n0[2]=angle_nz;
			optcalc_ray_trace(&opt,lambda,v0,n0,NULL);
			fprintf(fw,"%12g %12g\n",(v0[0]-sx0)/pixel_scale,(v0[1]-sy0)/pixel_scale);
		 }
	 }
	fprintf(fw,"# %12g %12g\n",sx0,sy0);

	fclose(fw);

  }

 if ( is_calc_transfer )
  {	transfer_matrix	m;
	double		zend,distance,eff_focus,z_focus;
	FILE		*fw;

	fw=stdout;

	optcalc_compute_transfer_matrix(&opt,lambda,m,0.0,&zend);
	distance=-m[0][0]/m[1][0];
	eff_focus=m[0][1]+distance*m[1][1];
	z_focus=zend+distance;
	fprintf(fw,"focal_plane:\t\t%12g\n",z_focus);
	fprintf(fw,"effective_focus:\t%12g\n",eff_focus);

	fclose(fw);
  }


 if ( output_scad != NULL )
  {	int	i;
	FILE	*fw;

	static	char *openscad_color[]={"red","green","blue","magenta","darkcyan","brown","grey"};

	if ( strcmp(output_scad,"-")==0 )
	 	fw=stdout;
	else if ( (fw=fopen(output_scad,"wb"))==NULL )
	 {	fprintf(stderr,"%s: error: unable to open or create output SCAD file '%s': %s\n",argv[0],output_scad,strerror(errno));
		return(1);
	 }

	fprintf(fw,"/* generated by %s */\n",argv[0]);
	for ( i=0 ; i<opt.opt_nlens ; i++ )
	 {	struct	lens	*l;
		struct	surface	*s1,*s2;
		int		j,jmax;
		double		z0,z1;
		l=&opt.opt_lenses[i];
		z0=l->l_offset-l->l_thickness/2.0;
		z1=l->l_offset+l->l_thickness/2.0;
		jmax=32;
		fprintf(fw,"color(\"%s\")rotate_extrude($fn=200)polygon([[%g,%g]",openscad_color[i%7],0.0,z0);
		s1=&l->l_s1;
		s2=&l->l_s2;
		for ( j=1 ; j<=jmax ; j++ )
		 {	double	r,z;
			r=l->l_radius1*(double)j/(double)jmax;
			optcalc_surface_aspheric_eval(s1->s_curvature,s1->s_conic,s1->s_nalpha,s1->s_alphas,r*r,&z);
			fprintf(fw,",[%g,%g]",r,z0+z);
		 }
		if ( l->l_radius1<l->l_radius2 )
		 {	double	r,z;
			r=l->l_radius1;
			optcalc_surface_aspheric_eval(s1->s_curvature,s1->s_conic,s1->s_nalpha,s1->s_alphas,r*r,&z);
			fprintf(fw,",[%g,%g]",l->l_radius2,z0+z);
		 }
		else if ( l->l_radius2<l->l_radius1 )
		 {	double	r,z;
			r=l->l_radius2;
			optcalc_surface_aspheric_eval(s2->s_curvature,s2->s_conic,s2->s_nalpha,s2->s_alphas,r*r,&z);
			fprintf(fw,",[%g,%g]",l->l_radius1,z1-z);
		 }
		for ( j=jmax ; 0<=j ; j-- )
		 {	double	r,z;
			r=l->l_radius2*(double)j/(double)jmax;
			optcalc_surface_aspheric_eval(s2->s_curvature,s2->s_conic,s2->s_nalpha,s2->s_alphas,r*r,&z);
			fprintf(fw,",[%g,%g]",r,z1-z);
		 }
		fprintf(fw,"]);\n");
		 				
	 }

	fclose(fw);
  }

 if ( output_eps != NULL )
  {	int	i;
	double	ox,oy,scale,rmax,zmin,zmax,sx,sy,bndx,bndy;
	FILE	*fw;

	if ( strcmp(output_eps,"-")==0 )
	 	fw=stdout;
	else if ( (fw=fopen(output_eps,"wb"))==NULL )
	 {	fprintf(stderr,"%s: error: unable to open or create output EPS file '%s': %s\n",argv[0],output_scad,strerror(errno));
		return(1);
	 }

	rmax=0.0;
	zmin=zmax=opt.opt_z_focal;
	for ( i=0 ; i<opt.opt_nlens ; i++ )
	 {	struct	lens	*l;
		double	z0,z1;
		l=&opt.opt_lenses[i];
		if ( rmax<l->l_radius1 )	rmax=l->l_radius1;
		if ( rmax<l->l_radius2 )	rmax=l->l_radius2;
		z0=l->l_offset-l->l_thickness/2.0;
		z1=l->l_offset+l->l_thickness/2.0;
		if ( z0<zmin )	zmin=z0;
		if ( zmax<z1 )	zmax=z1;
	 }

	zmin=floor(zmin);
	zmax=1+floor(zmax);
	rmax=1+floor(rmax);

	scale=2.0;

	bndx=10.0;
	bndy=10.0;

	ox=scale*(-zmin+bndx);
	sx=scale*(zmax-zmin+2*bndx);
	oy=scale*(rmax+bndy);
	sy=2*oy;

	fprintf(fw,"%%!PS-Adobe-2.0\n");
	fprintf(fw,"%%%%Title: %s\n",output_eps);
	fprintf(fw,"%%%%Creator: %s\n",argv[0]);
	fprintf(fw,"%%%%BoundingBox: 0 0 %g %g\n",sx,sy);

	for ( i=0 ; i<opt.opt_nlens ; i++ )
	 {	struct	lens	*l;
		struct	surface	*s0,*s1;
		double		z0,z1,r0,r1;
		double		c0,dz0,zh0,zf0;
		double		c1,dz1,zh1,zf1;
		l=&opt.opt_lenses[i];
		z0=l->l_offset-l->l_thickness/2.0;
		z1=l->l_offset+l->l_thickness/2.0;
		fprintf(fw,"%g %g moveto\n",ox+scale*z0,oy);
		s0=&l->l_s1;
		s1=&l->l_s2;
		r0=l->l_radius1;
		r1=l->l_radius2;
		optcalc_surface_aspheric_eval(s0->s_curvature,s0->s_conic,s0->s_nalpha,s0->s_alphas,r0*r0/4,&zh0);
		optcalc_surface_aspheric_eval(s0->s_curvature,s0->s_conic,s0->s_nalpha,s0->s_alphas,r0*r0  ,&zf0);
		optcalc_surface_aspheric_diff(s0->s_curvature,s0->s_conic,s0->s_nalpha,s0->s_alphas,r0*r0  ,&dz0);
		dz0=2*r0*dz0;
		if ( r0*dz0 != 0 )	c0=8.0*(0.5*zf0-zh0)/(3.0*r0*dz0);
		else			c0=1.0/3.0;
		optcalc_surface_aspheric_eval(s1->s_curvature,s1->s_conic,s1->s_nalpha,s1->s_alphas,r1*r1/4,&zh1);
		optcalc_surface_aspheric_eval(s1->s_curvature,s1->s_conic,s1->s_nalpha,s1->s_alphas,r1*r1  ,&zf1);
		optcalc_surface_aspheric_diff(s1->s_curvature,s1->s_conic,s1->s_nalpha,s1->s_alphas,r1*r1  ,&dz1);
		dz1=2*r1*dz1;
		if ( r1*dz1 != 0 )	c1=8.0*(0.5*zf1-zh1)/(3.0*r1*dz1);
		else			c1=1.0/3.0;
		fprintf(fw,"%g %g %g %g %g %g curveto\n",
			ox+scale*z0,oy+scale*(c0*r0),
			ox+scale*(z0+zf0-c0*r0*dz0),oy+scale*(1-c0)*r0,
			ox+scale*(z0+zf0),oy+scale*r0);

		if ( r0<r1 )
			fprintf(fw,"%g %g lineto %g %g lineto\n",
				ox+scale*(z0+zf0),oy+scale*r1,
				ox+scale*(z1-zf1),oy+scale*r1);
		else if ( r1<r0 )
		 	fprintf(fw,"%g %g lineto %g %g lineto\n",
				ox+scale*(z1-zf1),oy+scale*r0,
				ox+scale*(z1-zf1),oy+scale*r1);
		else
			fprintf(fw,"%g %g lineto\n",ox+scale*(z1-zf1),oy+scale*r0);

		fprintf(fw,"%g %g %g %g %g %g curveto\n",
			ox+scale*(z1-zf1+c1*r1*dz1),oy+scale*(1-c1)*r1,
			ox+scale*z1,oy+scale*(c1*r1),
			ox+scale*z1,oy);

		fprintf(fw,"%g %g %g %g %g %g curveto\n",
			ox+scale*z1,oy-scale*(c1*r1),
			ox+scale*(z1-zf1+c1*r1*dz1),oy-scale*(1-c1)*r1,
			ox+scale*(z1-zf1),oy-scale*r1);

		if ( r0<r1 )
		 	fprintf(fw,"%g %g lineto %g %g lineto\n",
				ox+scale*(z0+zf0),oy-scale*r1,
				ox+scale*(z0+zf0),oy-scale*r0);
		else if ( r1<r0 )
		 	fprintf(fw,"%g %g lineto %g %g lineto\n",
				ox+scale*(z1-zf1),oy-scale*r0,
				ox+scale*(z0+zf0),oy-scale*r0);
		else
			fprintf(fw,"%g %g lineto\n",
				ox+scale*(z0+zf0),oy-scale*r0);

		fprintf(fw,"%g %g %g %g %g %g curveto\n",
			ox+scale*(z0+zf0-c0*r0*dz0),oy-scale*(1-c0)*r0,
			ox+scale*z0,oy-scale*(c0*r0),
			ox+scale*z0,oy);

		fprintf(fw,"gsave 0.9 setgray fill grestore\n");
		fprintf(fw,"1 setlinewidth stroke\n");
	 }

	for ( i=-aperture_nring ; i<=aperture_nring ; i++ )
	 {	double	r,x,y;
		vector	v0,n0;
		struct	raytrace	rt;
		int	j;

		if ( aperture_nring )
			r=aperture_radius*(double)i/(double)aperture_nring;
		else
			r=0.0;
		x=0.0;
		y=r;
 		v0[0]=x;
		v0[1]=y;
		v0[2]=zstart-0.0001;
		n0[0]=angle_nx;
		n0[1]=angle_ny;
		n0[2]=angle_nz;
		optcalc_raytrace_reset(&rt);
		optcalc_ray_trace(&opt,lambda,v0,n0,&rt);
		for ( j=0 ; j<rt.rt_npoint ; j++ )
		 {	fprintf(fw,"%g %g %s\n",
				ox+scale*rt.rt_points[j][2],oy+scale*rt.rt_points[j][1],
				j?"lineto":"moveto");
		 }
		fprintf(fw,"0.5 setlinewidth stroke\n");
		optcalc_raytrace_free(&rt);
	 }


	fclose(fw);
  }

 if ( output_psf != NULL )
  {	triangle *triangles,*t;
	int	i,iteration_depth,max_iteration;
	double	*polygon,sx0,sy0;
	int	px,py;
	vector	v0,n0;

	FILE	*fw;

	if ( strcmp(output_psf,"-")==0 )
	 	fw=stdout;
	else if ( (fw=fopen(output_psf,"wb"))==NULL )
	 {	fprintf(stderr,"%s: error: unable to open or create output PSF file '%s': %s\n",argv[0],output_scad,strerror(errno));
		return(1);
	 }

	optcalc_glass_refraction_precompute(&opt,lambda);

	if ( 1 )
	 {	vector	v0,n0;
		v0[0]=0.0;
		v0[1]=0.0;
		v0[2]=zstart-0.0001;
		n0[0]=angle_nx;
		n0[1]=angle_ny;
		n0[2]=angle_nz;
		optcalc_ray_trace(&opt,lambda,v0,n0,NULL);
		sx0=v0[0];
		sy0=v0[1];
	 }

	triangles=NULL;
	for ( i=0 ; i<6 ; i++ )
	 {	double	w1,w2;
		t=list_new(triangle);
		t->t_vertex[0].x=0.0;
		t->t_vertex[0].y=0.0;
		w1=(double)(i+0)*M_PI/3.0;
		w2=(double)(i+1)*M_PI/3.0;
		t->t_vertex[1].x=aperture_radius*cos(w1);
		t->t_vertex[1].y=aperture_radius*sin(w1);
		t->t_vertex[2].x=aperture_radius*cos(w2);
		t->t_vertex[2].y=aperture_radius*sin(w2);
		t->t_flag=1;
		list_insert_first(triangles,t);
	 }

	max_iteration=7;
	for ( iteration_depth=0 ; iteration_depth<max_iteration ; iteration_depth++ )
	 {	for ( t=triangles ; t != NULL ; )
		 {	int	i,cnt;
			
			if ( ! t->t_flag )	
			 {	t=t->next;
				continue;
			 }
			cnt=0;
			for ( i=0 ; i<3 ; i++ )
			 {	v0[0]=t->t_vertex[i].x;
				v0[1]=t->t_vertex[i].y;
				v0[2]=zstart-0.0001;
				n0[0]=angle_nx;
				n0[1]=angle_ny;
				n0[2]=angle_nz;
				if ( ! optcalc_ray_trace(&opt,lambda,v0,n0,NULL) )
					cnt++;
			 }
			if ( cnt==3 )
			 {	t->t_flag=0;
				t=t->next;
				continue;
			 }
			else if ( cnt==0 )
			 {	triangle	*n;
				n=t->next;
				list_remove(triangles,t);
				t=n;
			 }
			else
			 {	triangle	*n,*t1,*t2,*t3,*t4;
				triangle_point	v01,v12,v20;
				n=t->next;
				t1=list_new(triangle);
				t2=list_new(triangle);
				t3=list_new(triangle);
				t4=list_new(triangle);
				v01.x=(t->t_vertex[0].x+t->t_vertex[1].x)/2.0;
				v01.y=(t->t_vertex[0].y+t->t_vertex[1].y)/2.0;
				v12.x=(t->t_vertex[1].x+t->t_vertex[2].x)/2.0;
				v12.y=(t->t_vertex[1].y+t->t_vertex[2].y)/2.0;
				v20.x=(t->t_vertex[2].x+t->t_vertex[0].x)/2.0;
				v20.y=(t->t_vertex[2].y+t->t_vertex[0].y)/2.0;
				t1->t_vertex[0]=t->t_vertex[0];
				t1->t_vertex[1]=v01;
				t1->t_vertex[2]=v20;
				t1->t_flag=1;
				t2->t_vertex[0]=t->t_vertex[1];
				t2->t_vertex[1]=v12;
				t2->t_vertex[2]=v01;
				t2->t_flag=1;
				t3->t_vertex[0]=t->t_vertex[2];
				t3->t_vertex[1]=v20;
				t3->t_vertex[2]=v12;
				t3->t_flag=1;
				t4->t_vertex[0]=v01;
				t4->t_vertex[1]=v12;
				t4->t_vertex[2]=v20;
				t4->t_flag=1;
				list_remove(triangles,t);
				list_insert_first(triangles,t1);
				list_insert_first(triangles,t2);
				list_insert_first(triangles,t3);
				list_insert_first(triangles,t4);
				t=n;
			 }
		 }	
	 }	

	/*
	for ( t=triangles ; t != NULL ; t=t->next )
	 {	if ( t->t_flag )	continue;
		fprintf(fw,"%12g %12g\n",t->t_vertex[0].x,t->t_vertex[0].y);
		fprintf(fw,"%12g %12g\n",t->t_vertex[1].x,t->t_vertex[1].y);
		fprintf(fw,"%12g %12g\n",t->t_vertex[2].x,t->t_vertex[2].y);
		fprintf(fw,"%12g %12g\n",t->t_vertex[0].x,t->t_vertex[0].y);
		fprintf(fw,"\n");
 	 }
	*/

		
	for ( t=triangles ; t != NULL ; t=t->next )
	 {	double	dx1,dy1,dx2,dy2;
		int	i;
		if ( t->t_flag )	continue;
		dx1=t->t_vertex[1].x-t->t_vertex[0].x;
		dy1=t->t_vertex[1].y-t->t_vertex[0].y;
		dx2=t->t_vertex[2].x-t->t_vertex[0].x;
		dy2=t->t_vertex[2].y-t->t_vertex[0].y;
		t->t_area=fabs(dx1*dy2-dx2*dy1)/2.0;
		for ( i=0 ; i<3 ; i++ )
		 {	v0[0]=t->t_vertex[i].x;
			v0[1]=t->t_vertex[i].y;
			v0[2]=zstart-0.0001;
			n0[0]=angle_nx;
			n0[1]=angle_ny;
			n0[2]=angle_nz;
			optcalc_ray_trace(&opt,lambda,v0,n0,NULL);
			t->t_vertex[i].x=(v0[0]-sx0)/pixel_scale;
			t->t_vertex[i].y=(v0[1]-sy0)/pixel_scale;
		 }
		dx1=t->t_vertex[1].x-t->t_vertex[0].x;
		dy1=t->t_vertex[1].y-t->t_vertex[0].y;
		dx2=t->t_vertex[2].x-t->t_vertex[0].x;
		dy2=t->t_vertex[2].y-t->t_vertex[0].y;
		t->t_proj=fabs(dx1*dy2-dx2*dy1)/2.0;
		/* fprintf(stderr,"%12g => %12g\n",t->t_area,t->t_proj); */
	 }

	polygon=malloc(sizeof(double)*2*16);
	for ( py=-psf_hsize ; py<=psf_hsize ; py++ )
 	 {	for ( px=-psf_hsize ; px<=psf_hsize ; px++ )
		 {	double		a;
			triangle	*t;
			a=0.0;
			for ( t=triangles ; t != NULL ; t=t->next )
			 {	int	j,n;
				double	w;
				if ( t->t_flag )	continue;
				for ( j=0 ; j<3 ; j++ )
				 {	polygon[2*j+0]=t->t_vertex[j].x;
					polygon[2*j+1]=t->t_vertex[j].y;
				 }
				n=polygon_intersection_square(polygon,3,(double)px-0.5,(double)py-0.5,1.0,1.0);
				if ( n<=0 )	continue;
				w=polygon_area(polygon,n);
				a+=t->t_area*w/t->t_proj;
			 }
			fprintf(fw,"%4d %4d %12g\n",psf_hsize+px,psf_hsize+py,a);
		 }
	 }
	free(polygon);

	fclose(fw);
	
	optcalc_glass_refraction_precompute(&opt,0.0);

  }

 return(0);
}
