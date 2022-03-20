/*****************************************************************************/
/* grselect.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line tool to select points from a list using the methods:	     */
/*   -	selecting points from a given box more or less homogeneously;	     */
/*   -	(de)selecting close pairs from a given list.			     */
/*****************************************************************************/
#define	FI_GRSELECT_VERSION	"1.0rc4"
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include "fi.h"

#include "io/iof.h"
#include "io/scanarg.h"
#include "io/tokenize.h"
#include "math/point.h"
#include "math/tpoint.h"
#include "pselect.h"

/*****************************************************************************/

#define	SELECT_NOSPECIAL	0
#define	SELECT_BY_ID		1
#define	SELECT_HOMOGENEOUS	2
#define	SELECT_PAIRS		3

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
 {	char	*line;
 } linedata;

typedef struct
 {	linedata	*lines;
	point		*points;
	int		nline;
 } dataset;

typedef struct
 {	char		*id;
 } iddata;

typedef struct
 {	int		number;
	int		level;
	int		is_mixed;
 } homdata;

typedef struct
 {	int		is_reject;
	double		distance;
	double		weightdiff;
 } pairdata;

/*****************************************************************************/

static int mix_dataset(dataset *ds)
{
 long		lrand48(void);
 int		i,j;
 linedata	wl;
 point		wp;

 for ( i=0 ; i<ds->nline ; i++ )
  {	j=(int)lrand48();
	j=j%(ds->nline-i)+i;
	if ( j==i )	continue;
	wl=ds->lines [i],ds->lines [i]=ds->lines [j],ds->lines [j]=wl;
	wp=ds->points[i],ds->points[i]=ds->points[j],ds->points[j]=wp;
  }
 return(0);
}

/*****************************************************************************/

int read_input_file(FILE *fr,dataset *ds,int colx,int coly,int colw,int neg_weight)
{
 char		*rbuff,*wbuff,*cmd[64];
 int		nline,nt;
 linedata	*lines;
 point		*points;
 double		x,y,w;

 lines=NULL;
 points=NULL;
 nline=0;

 while ( ! feof(fr) )
  {	rbuff=freadline(fr);
	if ( rbuff == NULL )	break;
	wbuff=(char *)malloc(strlen(rbuff)+1);

	strcpy(wbuff,rbuff);
	remove_newlines_and_comments(wbuff);
	nt=tokenize_spaces(wbuff,cmd,63);
	if ( colx>=nt || coly>=nt )
	 {	free(wbuff);free(rbuff);
		continue;
	 }
	sscanf(cmd[colx],"%lg",&x);
	sscanf(cmd[coly],"%lg",&y);
	if ( ! isfinite(x) || ! isfinite(y) )
		continue;
	if ( colw<nt && colw>=0 )
	 {	sscanf(cmd[colw],"%lg",&w);
		if ( ! isfinite(w) )	continue;
		if ( neg_weight )	w=-w;
	 }
	else
		w=0.0;
	free(wbuff);
	
	lines =(linedata *)realloc(lines ,(nline+1)*sizeof(linedata)),
	points=(point    *)realloc(points,(nline+1)*sizeof(point));

	points[nline].x=x,
	points[nline].y=y;
	points[nline].weight=w;

	lines[nline].line=rbuff;
	nline++;
  };
	
 ds->points=points;
 ds->lines =lines;
 ds->nline =nline;

 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int write_output_file(FILE *fw,dataset *ds,char *mask)
{
 int		i,k;
 linedata	*l;
 char		*buff;

 for ( i=0 ; i<ds->nline ; i++ )
  {	if ( mask != NULL )
	 {	if ( ! mask[i] )	continue;		}
	l=&ds->lines[i];
	if ( l==NULL )	continue;

	buff=l->line;

	k=strlen(buff);
	if ( k==0 || buff[k-1] != 10 )
		fprintf(fw,"%s\n",buff);
	else
		fprintf(fw,"%s",buff);
  };

 return(0);
}

/*****************************************************************************/

int select_from_section(point *points,int npoint,char *mask,double x1,double y1,double x2,double y2)
{
 int	ns;
 point	*p;
 double	w;

 if ( x2<x1 )	w=x1,x1=x2,x2=w;
 if ( y2<y1 )	w=y1,y1=y2,y2=w;

 ns=0;
 for ( p=points ; npoint ; npoint--,p++,mask++ )
  {	if ( x1<=p->x && p->x<=x2 && y1<=p->y && p->y<=y2 )
		*mask=1,ns++;
	else
		*mask=0;
  }
 
 return(ns);
}

/*****************************************************************************/

int self_check(double diff,double difflimit)
{
 if ( diff>-difflimit )	return(1);
 else			return(0);
}
int pair_check(double diff,double difflimit)
{
 if ( diff<difflimit )	return(1);
 else			return(0);
}

int select_pairs(point *points,int npoint,char *mask,pairdata *pd)
{
 tpoint	*tps;
 int	i,j,ns,ii,jj;
 double	dist,d2,x0,y0,dx,dy,diff;

 dist=pd->distance;
 d2=dist*dist;

 tps=(tpoint *)malloc(sizeof(tpoint)*npoint);
 for ( i=0 ; i<npoint ; i++ )
  {	tps[i].id=i;
	tps[i].xcoord=points[i].x,
	tps[i].ycoord=points[i].y;
	mask[i]=1;
  }

 qsort((void *)tps,npoint,sizeof(tpoint),tpoint_sortx);

 for ( i=0 ; i<npoint ; i++ )
  {	
	ii=tps[i].id;
	/*if ( ! mask[ii] )	continue;*/

	x0=tps[i].xcoord,
	y0=tps[i].ycoord;

	for ( j=i+1 ; j<npoint ; j++ )
	 {	jj=tps[j].id;
		/*if ( ! mask[jj] )		continue;*/
		if ( tps[j].xcoord>x0+dist )	break;
		dx=tps[j].xcoord-x0,
		dy=tps[j].ycoord-y0;
		if ( dx*dx+dy*dy > d2 )		continue;
		if ( pd->weightdiff<=0.0 )
		 {	mask[ii]=0;
			mask[jj]=0;
	 	 }
		else
		 {	/* diff=tps[j].prop-tps[i].prop; */
			diff=points[tps[j].id].weight-points[tps[i].id].weight;
			if ( self_check(diff,pd->weightdiff) )	mask[ii]=0;
			if ( pair_check(diff,pd->weightdiff) )	mask[jj]=0;
		 }
	 }

	for ( j=i-1 ; j>=0 ; j-- )
	 {	jj=tps[j].id;
		/*if ( ! mask[jj] )		continue;*/
		if ( tps[j].xcoord<x0-dist )	break;
		dx=tps[j].xcoord-x0,
		dy=tps[j].ycoord-y0;
		if ( dx*dx+dy*dy > d2 )		continue;	
		if ( pd->weightdiff<=0.0 )
		 {	mask[ii]=0;
			mask[jj]=0;
	 	 }
		else
		 {	/* diff=tps[j].prop-tps[i].prop; */
			diff=points[tps[j].id].weight-points[tps[i].id].weight;
			if ( self_check(diff,pd->weightdiff) )	mask[ii]=0;
			if ( pair_check(diff,pd->weightdiff) )	mask[jj]=0;
		 }
	 }
  }

 free(tps);

 ns=0;
 if ( ! pd->is_reject )
  {	for ( i=0 ; i<npoint ; i++ )
	 {	mask[i] = ! mask[i];		
		if ( mask[i] )	ns++;
	 }
  }
 else
  {	for ( i=0 ; i<npoint ; i++ )
	 {	if ( mask[i] )	ns++;	 	}
  }

 return(ns);
}

/*****************************************************************************/

int fprint_grselect_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\tgrselect [in] [-o out] [-h|--help]  [-V|--verbose]\n"
"\t[--offset <x1>,<y1> --size <sx>,<sy>|--section <x1>,<y1>,<x2>,<y2>]\n"
"\t[--col-xy <>,<>] [--col-weight [-]<>] [--col-id <>]\n");
 fprintf(fw,
"Selection by identifiers:\n"
"\t[--id <id>]\n");
 fprintf(fw,
"Homogenous selection:\n"
"\t[--homogeneous number=<N>,level=<L>,mix --selection <...>]\n");
 fprintf(fw,
"Selection by pairing criteria:\n"
"\t[--pairs [reject],distance=<distance>,weightdiff=<weight.diff>]\n");
 return(0);
}
int fprint_grselect_long_help(FILE *fw)
{ 
 return(0);
}

/*****************************************************************************/

int main(int argc,char *argv[])
{
 FILE		*fr,*fw;
 int		i,is_help,nline;
 int		colx,coly,colw,neg_weight,size_set,box_set;
 char		*infile,*outfile;
 char		*hdparam,*pdparam;

 dataset	ds_data,*ds=&ds_data;
 char		*mask;
 double		x1,y1,x2,y2,sx,sy;

 int		selectmethod;
 iddata		id;
 homdata	hd;
 pairdata	pd;

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 infile=outfile=NULL;

 is_comment=0;is_verbose=is_help=0;
 colx=1,coly=2;colw=-1,neg_weight=0;

 size_set=box_set=0;x1=y1=0.0;

 selectmethod=SELECT_NOSPECIAL;
 id.id=NULL; 
 hdparam=NULL;
	 hd.number=-1;
	 hd.level=-1;
	 hd.is_mixed=0;
 pdparam=NULL;
	 pd.is_reject=0;
	 pd.distance=0.0;
	 pd.weightdiff=0.0;

 i=scanarg(argc,argv,SCANARG_ALLOW_FLAGS,
	"--version:%NS-1f%q",&is_help,
        "--version-short|--short-version:%NS-2f%q",&is_help,
	"-h|--help|--help-short|--short-help:%f%q",&is_help,
	"--long-help|--help-long:%SN2f%q",&is_help,
	"-o|--output:%s",&outfile,

	"--col-xy:%d,%d",&colx,&coly,
	"--col-weight:%d",&colw,
	"--offset:%f%g,%g",&box_set,&x1,&y1,
	"--size:%f%f%g,%g",&box_set,&size_set,&sx,&sy,
	"--section:%f%g,%g,%g,%g",&box_set,&x1,&y1,&x2,&y2,

	"--id:" SNf(SELECT_BY_ID) "%s",&selectmethod,&id.id,
	"--homogeneous:" SNf(SELECT_HOMOGENEOUS) "%s",&selectmethod,&hdparam,
	"--pairs:" SNf(SELECT_PAIRS) "%s",&selectmethod,&pdparam,

	"--comment:%f",&is_comment,"(C):%f",&is_comment,
	"--verbose:%i",&is_verbose,"(V):%i",&is_verbose,
	"-*|+*:%e",
	"*:%w",&infile,
	NULL);

 if ( i )		
  {	fprint_error("invalid command line argument near '%s'",argv[i]);
	return(1);
  }
 else if ( is_help<0 )
  {	fprint_generic_version(stdout,argv[0],"grselect",FI_GRSELECT_VERSION,is_help);
	return(0);
  }
 else if ( is_help>1 )
  {	fprint_grselect_long_help(stdout);
	return(0);
  }
 else if ( is_help )
  {	fprint_grselect_usage(stdout);
	return(0);
  }

 if ( selectmethod==SELECT_HOMOGENEOUS && hd.number<0 )
  {	fprint_error("number of points to select is unknown, use the switch `-n` to specify");
	return(1);
  }

 colx--,coly--;
 if ( colx<0 || coly<0 || colw==0 )
  {	fprint_error("invalid column specification(s)");
	return(1);
  }
 if ( colw>0 )	colw--,neg_weight=0;
 else		colw=-colw-1,neg_weight=1;

 if ( infile==NULL )	fr=stdin;
 else			fr=fopenread(infile);
 if ( fr==NULL )	
  {	fprint_error("unable to open input file '%s'",infile);
	return(1);
  }
 read_input_file(fr,ds,colx,coly,colw,neg_weight);
 fcloseread(fr);

 if ( ! box_set )
  {	if ( ds->nline>0 )
	 {	x1=x2=ds->points[i].x,
		y1=y2=ds->points[i].y;
	 }
	for ( i=1 ; i<ds->nline ; i++ )
	 {	if ( ds->points[i].x>x2 )	x2=ds->points[i].x;
		if ( ds->points[i].y>y2 )	y2=ds->points[i].y;
		if ( ds->points[i].x<x1 )	x1=ds->points[i].x;
		if ( ds->points[i].y<y1 )	y1=ds->points[i].y;
	 }
  }
 else if ( ! size_set )
  {	x2=x1;y2=y1;
	for ( i=0 ; i<ds->nline ; i++ )
	 {	if ( ds->points[i].x>x2 )	x2=ds->points[i].x;
		if ( ds->points[i].y>y2 )	y2=ds->points[i].y;
	 }
  }
 else
  {	x2=x1+sx;
	y2=y1+sy;
  }

 if ( hdparam != NULL )
  {	i=scanpar(hdparam,SCANPAR_DEFAULT,
		"number:%d",&hd.number,
		"level:%d",&hd.level,
		"mix:%f",&hd.is_mixed,
		NULL);
	if ( i )
	 {	fprint_error("invalid homogenity parameter in '%s'",hdparam);
		return(1);
	 }
  }

 if ( pdparam != NULL )
  {	i=scanpar(pdparam,SCANPAR_DEFAULT,
		"distance:%g",&pd.distance,
		"reject:%f",&pd.is_reject,
		"weightdiff:%g",&pd.weightdiff,
		NULL);
	if ( i || pd.distance < 0.0 )
	 {	fprint_error("invalid pairing parameter in '%s'",pdparam);
		return(1);
	 }
  }

 if ( ds->nline > 0 )
  {	mask=(char *)malloc(sizeof(char)*ds->nline);
	memset(mask,0,ds->nline);

	switch ( selectmethod )
	 {  case SELECT_HOMOGENEOUS:
		if ( is_verbose )
			fprintf(stderr,"Number of points to select (homogenous): %d from %d\n",hd.number,ds->nline);
		if ( hd.is_mixed )	mix_dataset(ds);
		nline=point_select(ds->points,ds->nline,mask,hd.number,x1,y1,x2,y2,hd.level);
		break;
	     case SELECT_PAIRS:
		nline=select_pairs(ds->points,ds->nline,mask,&pd);
		break;
	     default:
		nline=select_from_section(ds->points,ds->nline,mask,x1,y1,x2,y2);
		break;
	  }
  }
 else
  {	nline=0;
	mask=NULL;
  }
	
 if ( is_verbose )
 fprintf(stderr,"Number of selected points: %d (from %d)\n",nline,ds->nline);
 
 if ( outfile==NULL )	fw=stdout;
 else			fw=fopenwrite(outfile);
 if ( fw==NULL )	
  {	fprint_error("unable to create output file '%s'",outfile);
	return(1);
  }
 if ( ds->nline>0 && nline>0 )
	write_output_file(fw,ds,mask);
 fclosewrite(fw);

 return(0);
}
