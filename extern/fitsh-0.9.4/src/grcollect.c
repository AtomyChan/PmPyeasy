/*****************************************************************************/
/* grcollect.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line tool to:						     */
/*  - transpose data read from multiple files;				     */
/*  - do various kind of statistics.					     */
/*****************************************************************************/
#define	FITSH_GRCOLLECT_VERSION	"1.0pre1"
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>

#include "io/iof.h"
#include "io/scanarg.h"
#include "io/tokenize.h"
#include "cache.h"
#include "statistics.h"
#include "str.h"

#include "longhelp.h"
#include "fitsh.h"

#ifdef  HAVE_NO_CC_EXTENSION 
#define __extension__ 
#endif 

/*****************************************************************************/

#define		STAT_NONE		0
#define		STAT_COUNT		1
#define		STAT_MEAN		2
#define		STAT_MEANDEV		3
#define		STAT_MIN		4
#define		STAT_MAX		5
#define		STAT_SUM		6
#define		STAT_SUM2		7

#define		STAT_PURE_CUMULATIVE	7

#define		STAT_MEDIAN		8
#define		STAT_MEDIANDEV		9
#define		STAT_MEANMEDDEV		10
#define		STAT_MEDIANMEDDEV	11
#define		STAT_MODE		12
#define		STAT_MODEDEV		13
#define		STAT_MODEMEDDEV		14

#define		STAT_RCOUNT		21
#define		STAT_RMEAN		22
#define		STAT_RMEANDEV		23
#define		STAT_RMIN		24
#define		STAT_RMAX		25
#define		STAT_RSUM		26
#define		STAT_RSUM2		27
#define		STAT_RMEDIAN		28
#define		STAT_RMEDIANDEV		29
#define		STAT_RMEANMEDDEV	30
#define		STAT_RMEDIANMEDDEV	31
#define		STAT_RMODE		32
#define		STAT_RMODEDEV		33
#define		STAT_RMODEMEDDEV	34

typedef struct
 {	char    *name;
	int     code;
 } dataname;

static dataname datanamelist[]=
 {	{ "count",		STAT_COUNT		},
	{ "mean",		STAT_MEAN		},
	{ "average",		STAT_MEAN		},
	{ "avg",		STAT_MEAN		},
	{ "stddev",		STAT_MEANDEV		},
	{ "sigma",		STAT_MEANDEV		},
	{ "meandev",		STAT_MEANDEV		},
	{ "meanstddev",		STAT_MEANDEV		},
	{ "min",		STAT_MIN		},
	{ "max",		STAT_MAX		},
	{ "sum",		STAT_SUM		},
	{ "sum2",		STAT_SUM2		},

	{ "median",		STAT_MEDIAN		},
	{ "mediandev",		STAT_MEDIANDEV		},
	{ "medianstddev",	STAT_MEDIANDEV		},
	{ "meanmeddev",		STAT_MEANMEDDEV		},
	{ "medianmeddev",	STAT_MEDIANMEDDEV	},
	{ "mode",		STAT_MODE		},
	{ "modedev",		STAT_MODEDEV		},
	{ "modestddev",		STAT_MODEDEV		},
	{ "modemeddev",		STAT_MODEDEV		},

	{ "rcount",		STAT_RCOUNT		},
	{ "rmean",		STAT_RMEAN		},
	{ "raverage",		STAT_RMEAN		},
	{ "ravg",		STAT_RMEAN		},
	{ "rstddev",		STAT_RMEANDEV		},
	{ "rsigma",		STAT_RMEANDEV		},
	{ "rmeandev",		STAT_RMEANDEV		},
	{ "rmeanstddev",	STAT_RMEANDEV		},
	{ "rmin",		STAT_RMIN		},
	{ "rmax",		STAT_RMAX		},
	{ "rsum",		STAT_RSUM		},
	{ "rsum2",		STAT_RSUM2		},
	{ "rmedian",		STAT_RMEDIAN		},
	{ "rmediandev",		STAT_RMEDIANDEV		},
	{ "rmedianstddev",	STAT_RMEDIANDEV		},
	{ "rmeanmeddev",	STAT_RMEANMEDDEV	},
	{ "rmedianmeddev",	STAT_RMEDIANMEDDEV	},
	{ "rmode",		STAT_RMODE		},
	{ "rmodedev",		STAT_RMODEDEV		},
	{ "rmodestddev",	STAT_RMODEDEV		},
	{ "rmodemeddev",	STAT_RMODEDEV		},
	{ NULL,			-1			}
 };

static int	statlist_default[] = { STAT_COUNT,STAT_MEAN,STAT_MEANDEV,-1 };

/*****************************************************************************/

typedef struct
 {	char	*key;
	int	id;
 } key;

typedef struct lookupkey lookupkey;

struct lookupkey
 {	lookupkey	*lookups;
	key		value;
 };

/*****************************************************************************/

typedef struct
 {	size_t	offset;		/* offset in collectbuffer->buffer	*/
	int	key;		/* position of the base column		*/
	int	file;		/* infiles[] array index		*/
	char	*keyptr;	/* pointer to base column		*/
 } line;

typedef struct
 {	lookupkey	lk;
	char		*buffer;
	size_t		size,asize;
	line		*lines;
	int		nline;
	int		is_comment;
	char		**comments;
	int		argc;
	char		**argv;
 } collectbuffer;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define			REJECTION_CENTER	0x03
#define			REJECTION_C_NONE	0
#define			REJECTION_C_MEAN	1
#define			REJECTION_C_MEDIAN	2
#define			REJECTION_C_MODE	3
#define			REJECTION_MAGNITUDE	0x0C
#define			REJECTION_M_STDDEV	0
#define			REJECTION_M_MEDDEV	4
#define			REJECTION_M_ABSOLUTE	8

typedef struct
 {	int		center;
	double		magnitude;
	int		iterations;
 } rejection;

typedef struct
 {	int		index;
	rejection	r;
 } column;

typedef struct
 {	int	id;
	int	flag;
 } tmpvaluecommon;

typedef struct
 {	char	*key;
	int	id;
	int	count;
 } statfieldcommon;

#define		N_CENTER		3	/* mean, median, mode	     */
#define		C_MEAN			0
#define		C_MEDIAN		1
#define		C_MODE			2
#define		N_DEVIATION		3	/* std, med, abs	     */
#define		D_STD			0
#define		D_MED			1
#define		D_ABS			2

typedef struct
 {	double	center[N_CENTER];
	double	deviation[N_CENTER][N_DEVIATION];
 } cendev;
typedef struct
 {	int	count;
	double	sum,sum2,min,max;
 } bcumul;

typedef struct
  {	double	sum,sum2,min,max;		/* cumulative		     */
	double	median,mediandev;		/* normal non-cumulative 1)  */
	double	meanmeddev,medianmeddev;	/* normal non-cumulative 2)  */
	double	mode,modedev,modemeddev;	/* normal non-cumulative 3)  */
	int	rcount;				/* rej'd  non-cumulative 1)  */
	double	rsum,rsum2,rmin,rmax;		/* rej'd  non-cumulative 2)  */
	double	rmedian,rmediandev;		/* rej'd  non-cumulative 3)  */
	double	rmeanmeddev,rmedianmeddev;	/* rej'd  non-cumulative 4)  */
	double	rmode,rmodedev,rmodemeddev;	/* rej'd  non-cumulative 5)  */
 } statfielddata;

typedef struct
 {	bcumul	nb;
	cendev	nc;
	bcumul	rb;
	cendev	rc;
 } statfielddata2;

typedef struct
 {	lookupkey	lk;		/* hash tree for cumulative stat's   */
	void		*fields;	/* keys				     */
	int		nfield,afield;	/* number of keys (+ alloc'ed)	     */
	char		*tmptemplate;	/* template name of temporary file   */
	int		fh;		/* file handle of the temporary file */
	int		nstat;		/* number of columns for making stat */
	column		*colstats;	/* column spec info (for rejection)  */
	void		*tmps;		/* list of temporary values, cached  */
	size_t		ntmp;		/* number of temporary values.       */
	size_t		atmp;		/* number of allocated entries (tmp) */
	off_t		wtmp;		/* num of written temp entries	     */
	size_t		mxtmp;		/* max number of temporary values    */
	size_t		maxmem;		/* max memory in bytes		     */
 } collectstat;

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

int	*get_stat_list(char *statstr)
{
 char	*wstr,*cmd[16];
 int	i,n,*ret,j;

 wstr=strdup(statstr);
 n=tokenize_char(statstr,cmd,',',15);

 if ( n <= 0 )
  {	free(wstr);
	return(NULL);
  }
 ret=(int *)malloc(sizeof(int)*(n+1));

 for ( i=0 ; i<n ; i++ )
  {	for ( j=0 ; datanamelist[j].name != NULL ; j++ )
	 {	if ( strcmp(cmd[i],datanamelist[j].name)==0 )	break;	}
	if ( datanamelist[j].name==NULL || datanamelist[j].code<0 )
	 {	free(ret);
		free(wstr);
	 	return(NULL);
	 }
	else	ret[i]=datanamelist[j].code;
  }
 ret[n]=-1;

 free(wstr);

 return(ret);
}

/*****************************************************************************/

#define		HASH_BITS		2

int lookup_get_hash_bit(char *kstr,int offset)
{
 int	k;
#if HASH_BITS == 8
 k=(unsigned char)kstr[offset];
#elif HASH_BITS == 4
 int	w;
 w=(unsigned char)kstr[offset/2];
 if ( offset%2==0 )	k=(w>>4) & 0x0f;
 else			k=(w>>0) & 0x0f;
#elif HASH_BITS == 2 
 int	w;
 w=(unsigned char)kstr[offset/4];
 k=(w>>(2*(3-(offset%4)))) & 0x03;
#elif HASH_BITS == 1 
 int	w;
 w=(unsigned char)kstr[offset/8];
 k=(w>>(7-(offset%8))) & 0x01;
#else
#error "Invalid hash bit value"
#endif
 return(k);
}

int lookup_get_end(char *kstr,int offset)
{
 int	k;
#if HASH_BITS == 8
 k=(unsigned char)kstr[offset];
#elif HASH_BITS == 4
 k=(unsigned char)kstr[offset/2];
#elif HASH_BITS == 2 
 k=(unsigned char)kstr[offset/4];
#elif HASH_BITS == 1 
 k=(unsigned char)kstr[offset/8];
#else
#error "Invalid hash bit value"
#endif
 return(k);
}

key * lookup_search_or_add_key(lookupkey *lk,char *kstr,int offset)
{
 int		k,w,o;
 key		*ret,tmp;
 char		*wkey;
 lookupkey	*wl;

 k=lookup_get_hash_bit(kstr,offset);

 if ( lk->lookups==NULL )		/* empty set */
  {	lk->lookups=(lookupkey *)malloc(sizeof(lookupkey)*(1<<HASH_BITS));
	memset(lk->lookups,0,sizeof(lookupkey)*(1<<HASH_BITS));
	lk->lookups[k].lookups=NULL;
	lk->lookups[k].value.key=strdup(kstr);
	ret=&lk->lookups[k].value;
  }
 else if ( lk->lookups[k].lookups != NULL )
  {	ret=lookup_search_or_add_key(&lk->lookups[k],kstr,offset+1);	}
 else if ( lk->lookups[k].value.key==NULL )
  {	lk->lookups[k].value.key=strdup(kstr);
	ret=&lk->lookups[k].value;
  }
 else if ( strcmp(lk->lookups[k].value.key,kstr)==0 )
  {	ret=&lk->lookups[k].value;		}
 else
  {	wkey=lk->lookups[k].value.key;
	memcpy(&tmp,&lk->lookups[k].value,sizeof(key));
	memset(&lk->lookups[k].value,0,sizeof(key));
	wl=&lk->lookups[k];
	for ( o=offset+1,ret=NULL ; lookup_get_end(wkey,o) || lookup_get_end(kstr,o) ; o++ )
	 {	wl->lookups=(lookupkey *)malloc(sizeof(lookupkey)*(1<<HASH_BITS));
		memset(wl->lookups,0,sizeof(lookupkey)*(1<<HASH_BITS));
		k=lookup_get_hash_bit(kstr,o);
		w=lookup_get_hash_bit(wkey,o);
		if ( k==w )
			wl=&wl->lookups[k];
		else
		 {	memcpy(&wl->lookups[w].value,&tmp,sizeof(key));
			wl->lookups[k].value.key=strdup(kstr);
			ret=&wl->lookups[k].value;
			break;
		 }
	 }
  }
 return(ret);
}

int lookup_free_node(lookupkey *lk)
{
 int		i;
 lookupkey	*wl;
 if ( lk==NULL || lk->lookups==NULL )
	return(0);
 for ( i=0 ; i<(1<<HASH_BITS) ; i++ )
  {	wl=&lk->lookups[i];
	if ( wl->lookups != NULL )
		lookup_free_node(wl);
	else if ( wl->value.key != 0 )
		free(wl->value.key);
  }
 free(lk->lookups);
 return(0);
}

/******************************************************************************
int lookup_dump(FILE *fw,lookupkey *lk,int depth)
{
 int	i,j;

 if ( lk==NULL || lk->lookups==NULL )
	return(0);
 for ( i=0 ; i<256 ; i++ )
  {	if ( lk->lookups[i].value.key != NULL )
	 {	for ( j=0 ; j<depth ; j++ )
		 {	fprintf(fw,"+ ");		}
		fprintf(fw,"%s\n",lk->lookups[i].value.key);
	 }
	else if ( lk->lookups[i].lookups != NULL )
		lookup_dump(fw,&lk->lookups[i],depth+1);
  }
 return(0);
}

int lookup_test(void)
{
 lookupkey	lk;

 memset(&lk,0,sizeof(lookupkey));

 lookup_search_or_add_key(&lk,"B1",0);
 lookup_search_or_add_key(&lk,"B2",0);
 lookup_search_or_add_key(&lk,"A2",0);
 lookup_search_or_add_key(&lk,"B2c",0);

 lookup_dump(stderr,&lk,0);

 return(0);
}
******************************************************************************/

/*****************************************************************************/

int fprint_grcollect_usage(FILE *fw)
{
 fprintf(fw,
"Usage:\tgrcollect [-h|--help|--long-help|--wiki-help] [--version[-short]]\n"
"\t[-V|--verbose]\n"
"Options for data collection and transposition:\n"
"\t<input> [<input> ...] -c|--col-base <> [-t|--tmpdir <tmp-dir>]\n"
"\t{-b|--basename <base-%%b-name>|-x|--extension <xtn> -p|--prefix <pre>}\n"
"\t[--rejection column=<col>,{mean|median|mode},\n"
"\t             {stddev|meddev|absolute}=<magnitude>,[comment|exclude]]\n"
"\t[-a|--pre-allocate <size>[,<nblock>]\n"
"\t[-C|--comment] [-S|--additional-comment <...> [-S ...]]\n");
 fprintf(fw,
"Options for cumulative statistics:\n"
"\t<input> [...] -c|--col-base <> -d|--col-stat <>[,...] -o|--output <out>\n"
"\t[--stat [r]count,[r]{mean|median|mode}[[med|std]dev],\n"
"\t        [r]sum,[r]sum2,[r]min,[r]max]\n"
"\t[--rejection column=<col>,iterations=<n>,{mean|median|mode},\n"
"\t             {stddev|meddev|absolute}=<magnitude>] [--rejection ...]\n");
 fprintf(fw,
"Other options:\n"
"\t[-t|--tmpdir <directory-for-temporary-file>]\n"
"\t[-l|--line-size <maximum-line-size>] [-m|--max-memory <max-memory>]\n");
 fprintf(fw,
"Notes:\n"
" - the default directory for tempoary files is always the current one (which\n"
"   is equivalent to define `--tmpdir ./`);\n"
" - the suffixes 'k', 'm' or 'g' can be added after the amount of max. memory,\n"
"   however, this amount is only a hint, the actual usage depends on the data\n"
"   (the default value of max. memory is %d megabytes)\n",DEFAULT_MAX_MEMORY);

 return(0);
}

longhelp_entry grcollect_long_help[] = 
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

 { "<input> [,<input>, ...]",
	"Name of the input file. At least, one file should be specified. "
	"Reading from standard input can be forced using a single dash \"-\" "
	"as input file name. More dashes are silently ignored." },
 { "-c, --col-base <key column index>",
	"Column index for the key." },

 { "Data transposition specific options:", NULL },
 { "-b, --basename <base-%b-name>",
	"Base name of the output files. The base name string should conatain "
	"at least one \"%b\" tag, which is replaced by the respective key "
	"string on the creation of the file." },
 { "-x, --extension <extension>, -p, --prefix <prefix>",
	__extension__
	"Equivalent to \"-b|--basename <prefix>%b.<extension>\". Note that "
	"in practice, <prefix> might be some sort of directory name and extension is "
	"a regular file extension, but the above substitution is done literally. "
	"Therefore, the \"dot\" between the key and the <extension> is always "
	"inserted in the final name of the output files but a trailing "
	"slash is required at the end of <prefix> if the files are to be "
	"created in that particular directory. Note also that this case, the target "
	"directory must exist before the invocation of `grcollect`, otherwise "
	"the output files cannot be created. " },
 { "-C, --comment",
	"Insert a commented line (starting with \"#\") containing information "
	"about the version and command line invocation syntax of `grcollect` "
	"to the beginning of the transposed files." },
 { "-S, --additional-comment <...>",
	"Insert an additional commented lines (starting with \"#\") "
	"to the beginning of the transposed files." },

 { "Options for cumulative statistics:", NULL },
 { "-d, --col-stat <>[,...]",
	"Comma-separated list of column indices on which the statistics are "
	"to be calculated. Columns with non-numerical contents are ignored." 
	"Note that this option imply the cumulative statistics mode of "
	"`grcollect`." },
 { "-o, --output <filename>",
	"The name of the output file to which the output statistics are "
	"written. The total number of columns in this file will be "
	"1+C*N, where C is the number of columns (see -d|--col-stat) "
	"on which the statistics are calculated and N is the number of "
	"statistic quantities (see --stat). The first column in the output "	
	"file is the key, which is followed by the per-column list of statistics, "
	"in the same order as the user defined after -d|--col-stat and "
	"--stat." },
 { "-s, --stat <list of statistics>",
	"Comma-separated list of statistics to be estimated on the input data. "
	"These can be one or more of the following:" },
 { "count",
	"Total number of records, for the given key." },
 { "rcount",
	"The number of records after rejecting outliers (i.e. it is always "
	"the same as the \"count\" value if no \"--rejection\" was used). " },
 { "mean, median, mode",
	"Mean, median or mode statistics of the data." },
 { "rmean, rmedian, rmode",
	"Mean, median or mode, after rejecting outliers." },
 { "{mean|median|mode}stddev, {mean|median|mode}meddev, stddev",
	"Scatter of the data around the mean, median or mode. The scatter can "
	"either be standard deviation (stddev) or median deviance (meddev). "
	"The literal \"stddev\" is the classic standard deviation, "
	"equivalent to \"meanstddev\". " },
 { "r{mean|median|mode}stddev, r{mean|median|mode}meddev, rstddev",
	"The same scatters as above but after rejecting outliers." },
 { "sum, rsum",
	"Sum of the data, esp. total sum and sum after rejecting outliers." },
 { "sum2, rsum2",
	"Sum of the squares, total and after rejecting outliers." },
 { "min, max",
	"Minimal and maximal data values." },
 { "rmin, rmax",
	"Minimal and maximal data values after the rejection of outliers." },

 { "-r, --rejection column=<index>,<rejection parameters>",
	"Comma-separated directives for outlier rejection for the "
	"specified column. The rejection parameters are:" },
 { "iterations=<n>",
	"Maximum number of iterations to reject outliers." },
 { "mean, median, mode",
	"Use the mean, median or mode for the center of the rejection." },
 { "stddev, meddev, absolute=<limit>",
	"Use the standard deviation or median deviance  "
	"for rejection limit units or define an absolute limit for "
	"rejection level." },

 { "Note that each column can have different kind of rejection method, thus "
   "more than one \"--rejection ...\" command line option can be used at "
   "the invocation of `grcollect`.", NULL }, {"", NULL },
 
 { "Other options:", NULL },
 { "-m, --max-memory <memory>[kmg]",
	__extension__
	"Maximum amount of memory available for `grcollect`. The prefixes "
	"\"k\", \"m\" or \"g\" can be used for kilobytes, megabytes and "
	"gigabytes, respectively. On 32bit systems, the maximum memory is "
	"limited to 3gigabytes. Note that `grcollect` does not use any "
	"kind of operating system specific methods to determine the maximum "
	"amount of memory, it always should be set by the user. The default "
	"value of 8 megabytes is somewhat small, so upon massive data "
	"transposition (tens or hundreds of gigabytes), this limit is "
	"worth to be set accordingly to the physical memory available. " },
 { "-t, --tmpdir <directory>",
	"Directory for temporary file storage. Note that the default "
	"temporary directory is always the current one (which is "
	"is equivalent to define \"--tmpdir ./\"), since in a usual "
	"configuration the /tmp directory is small, moreover, it can be "
	"some sort of \"tmpfs\", temporary file system mount on the "
	"physical memory itself. " },

 { NULL, NULL }
};

int fprint_grcollect_long_help(FILE *fw,int is_wiki)
{ 
 char	*synopsis=
	"grcollect [options] <input> [...] [-o <output>|-b <basename>]";
 char	*description=
	__extension__
	"The main purpose of the program `grcollect` is twofold. First, it is intended "
	"to do data transposition on the input data, i.e. the input (which is read "
	"from files or standard input) is sorted and splitted to separate files where "
	"the splitting is based on a respective key. These keys are taken from "
	"the input data. In such a case where the input is from more files and each "
	"key is unique in a given file, this process is called data transposition "
	"(since it is similar when a 2 dimensional data matrix is stored in the form "
	"as each row is in a separate file, and one intends to transpose the matrix, "
	"i.e. store each column in a separate file). "
	"The other feature of `grcollect` is to do some sort of statistics on data "
	"associated to different keys. These statistics include average (mean, median, "
	"mode) and scatter (standard deviation or median deviance) estimations with "
	"the optional deselection of outlier points, summation, count statistics and "
	"so on.";

 fprint_generic_long_help(fw,is_wiki,grcollect_long_help,synopsis,description);

 return(0);
}

/*****************************************************************************/

int collect_buffer_reset(collectbuffer *cb)
{

 cb->buffer=NULL;
 cb->size=(size_t)0;
 cb->asize=(size_t)0;
 cb->lines=NULL;
 cb->nline=0;
 return(0);
}
int collect_buffer_free(collectbuffer *cb)
{
 if ( cb->buffer != NULL )	free(cb->buffer);
 if ( cb->lines != NULL )	free(cb->lines);
 collect_buffer_reset(cb);
 return(0);
}

int get_column_position(char *line,int col)
{
 int	pos;

 pos=0;
 while ( isspace((int)*line) )	line++,pos++;
 if ( ! (*line) )		return(-1);
 while ( col>0 )
  {	while ( ! isspace((int)*line) && *line )	line++,pos++;
	while ( isspace((int)*line) )			line++,pos++;
	if ( ! (*line) )				return(-1);
	col--;
  }
 return(pos);
}
int compare_line_by_keyptr(const void *vl1,const void *vl2)
{
 line	*l1=(line *)vl1;
 line	*l2=(line *)vl2;
 return(strcmp(l1->keyptr,l2->keyptr));
}

#define		BUFFSIZE	256

int append_string(char **rbuff,int *rsize,int *rasize,char *str,int len)
{
 if ( str==NULL )	return(0);

 if ( *rsize+len > *rasize )
  {	while ( *rsize+len > *rasize )	(*rasize)+=BUFFSIZE;
	*rbuff=(char *)realloc(*rbuff,*rasize);
  }
 memcpy((*rbuff)+(*rsize),str,len);
 (*rsize)+=len;
 return(len);
}

char *expand_basename(char *basename,char *key)
{
 char	*out;
 int	size,asize,klen,i;

 asize=size=0;
 out=NULL;

 if ( key != NULL )	klen=strlen(key);
 else			klen=0;

 if ( klen>0 )		/* ugly, very ugly, suxx */
  {	for ( i=0 ; i<klen ; i++ )
	 {	if ( isspace((int)key[i]) )
			klen=i;
	 }
  }

 while ( *basename )
  {	if ( *basename=='%' )
	 {	basename++;
		switch ( *basename )
		 {   case 'b':
			append_string(&out,&size,&asize,key,klen);
			break;
		     default:
			append_string(&out,&size,&asize,basename,1);
		 }
		basename++;
	 }
	else
	 {	append_string(&out,&size,&asize,basename,1);
		basename++;
	 }
  }
 append_string(&out,&size,&asize,"",1);
 return(out);
}

/*
int collect_flush(collectbuffer *cb,char *basename,off_t preallocate)
{
 int	i;
 line	*wl;
 char	*filename,*prevfile;
 FILE	*fw;
 key	*wk;

 for ( i=0 ; i<cb->nline ; i++ )
  {	wl=&cb->lines[i];
	wl->keyptr=cb->buffer+wl->offset+(size_t)wl->key;
  }
 qsort(cb->lines,cb->nline,sizeof(line),compare_line_by_keyptr);

 prevfile=NULL;
 fw=NULL;
 for ( i=0 ; i<cb->nline ; i++ )
  {	wl=&cb->lines[i];
	filename=expand_basename(basename,wl->keyptr);
	if ( ( prevfile != NULL && strcmp(prevfile,filename) ) || prevfile==NULL || fw==NULL )
	 {	if ( fw != NULL )	fclose(fw);
		wk=lookup_search_or_add_key(&cb->lk,filename,0);
		if ( wk->id <= 0 )
		 {	fw=fopen(filename,"wb");
			wk->id=1;
		 }
		else
			fw=fopen(filename,"ab");
	 }
	if ( fw != NULL )	fprintf(fw,"%s\n",cb->buffer+wl->offset);

	if ( prevfile != NULL )		free(prevfile);
	prevfile=filename;
  }

 if ( prevfile != NULL )	free(prevfile);
 if ( fw != NULL )		fclose(fw);

 collect_buffer_free(cb);

 return(0);
}
*/

/*****************************************************************************/

#define	HPRINT_BUFFER_SIZE		256

int hprintf(int handle,char *msg,...)
{
 char		buff[HPRINT_BUFFER_SIZE],*tbuff;
 va_list	ap;
 int		n;

 va_start(ap,msg);
 n=vsnprintf(buff,HPRINT_BUFFER_SIZE,msg,ap);
 va_end(ap);
 if ( n<HPRINT_BUFFER_SIZE )
  {	(void)write(handle,buff,n);
	return(n);
  }
 else
  {	tbuff=NULL;
	va_start(ap,msg);
	vstrappendf(&tbuff,msg,ap);
	va_end(ap);
	if ( tbuff != NULL )
	 {	n=strlen(tbuff);
		(void)write(handle,tbuff,n);
		free(tbuff);
	 }
	else
		n=0;

	return(n);
  }

}

/*****************************************************************************/

int collect_flush(collectbuffer *cb,char *basename,off_t preallocate)
{
 int	i,j,*pincr;
 line	*wl;
 char	*filename,*prevfile,*cline;
 int	fd;
 key	*wk;

 for ( i=0 ; i<cb->nline ; i++ )
  {	wl=&cb->lines[i];
	wl->keyptr=cb->buffer+wl->offset+(size_t)wl->key;
  }
 qsort(cb->lines,cb->nline,sizeof(line),compare_line_by_keyptr);

 prevfile=NULL;
 fd=-1;
 pincr=NULL;

 for ( i=0 ; i<cb->nline ; i++ )
  {	
	wl=&cb->lines[i];
	filename=expand_basename(basename,wl->keyptr);

	if ( ( prevfile != NULL && strcmp(prevfile,filename) ) || prevfile==NULL || fd<0 )
	 {	if ( fd >= 0 )
		 {	close(fd);
			fd=-1;
		 }

		wk=lookup_search_or_add_key(&cb->lk,filename,0);

		if ( wk->id <= 0 )
		 {	fd=open(filename,O_CREAT|O_TRUNC|O_RDWR,0666);
			if ( cb->is_comment>=1 )
			 {	hprintf(fd,"# Created by grcollect %s (fi: %s)\n",FITSH_GRCOLLECT_VERSION,FITSH_VERSION);	}
			if ( cb->is_comment>=2 )
			 {	hprintf(fd,"# Invoked command:");
				for ( j=0 ; j<cb->argc ; j++ )
				 {	if ( is_any_nasty_char(cb->argv[j]) )
						hprintf(fd," \"%s\"",cb->argv[j]);
					else
						hprintf(fd," %s",cb->argv[j]);
				 }
				hprintf(fd,"\n");
			 }
			for ( j=0 ; cb->comments != NULL && cb->comments[j] != NULL ; j++ )
			 {	cline=cb->comments[j];
				hprintf(fd,"# %s\n",cline);
			 }			
			wk->id=0;
			if ( preallocate>0 )
			 {	char	*buff;
				int	i,pagesize;
				off_t	size;
#ifdef	HOST_WIN32
				pagesize=4096;
#else
				pagesize=getpagesize();
#endif
				buff=(char *)malloc((size_t)pagesize);
				for ( i=0 ; i<pagesize ; i++ )
				 {	buff[i]=i%256;			}
				size=preallocate;
				while ( size>0 )
				 {	if ( size>(off_t)pagesize )
						i=pagesize;
					else
						i=(int)size;
					(void)write(fd,buff,i);
					size-=(off_t)i;
				 }
				free(buff);
				lseek(fd,(off_t)0,SEEK_SET);
			 }
		 }

		else
		 {	fd=open(filename,O_RDWR);
			lseek(fd,(off_t)wk->id,SEEK_SET);
		 }

		pincr=&wk->id;
	 }

	if ( fd >= 0 )
	 {	int	len;
		char	*wrbuff;
		wrbuff=cb->buffer+wl->offset;
		len=strlen(wrbuff);
		(void)write(fd,wrbuff,(size_t)len);
		(void)write(fd,"\n",(size_t)1);
		if ( pincr != NULL )
			(*pincr)+=len+1;
	 }

	if ( prevfile != NULL )		free(prevfile);

	prevfile=filename;
  }

 if ( prevfile != NULL )	free(prevfile);

 if ( fd>=0 )
	close(fd);

 collect_buffer_free(cb);

 return(0);
}

#define		FRAGSIZE	16384

int collect_read_file(collectbuffer *cb,int fid,FILE *fr,int colbase,
	char *basename,size_t maxmem,off_t preallocate)
{
 char	*buff;
 int	key,len;
 size_t	n;
 line	*wl;

 while ( ! feof(fr) )
  {	buff=freadline(fr);
	if ( buff==NULL )	break;
	remove_newlines_and_comments(buff);
	key=get_column_position(buff,colbase);
	if ( key<0 )
	 {	free(buff);
		continue;
	 }
	cb->lines=(line *)realloc(cb->lines,sizeof(line)*(cb->nline+1));
	wl=&cb->lines[cb->nline];
	cb->nline++;

	len=strlen(buff)+1;	
	if ( cb->size+(size_t)len > cb->asize )
	 {	n=(cb->size+(size_t)(len+FRAGSIZE-1))/(size_t)FRAGSIZE;
		cb->asize=(size_t)n*(size_t)FRAGSIZE;
		cb->buffer=(char *)realloc(cb->buffer,cb->asize);
	 }
	memcpy(cb->buffer+(size_t)cb->size,buff,(size_t)len);

	free(buff);

	wl->offset=cb->size;
	wl->key=key;
	wl->file=fid;
	cb->size+=(size_t)len;

	if ( cb->size > maxmem )
	 {	collect_flush(cb,basename,preallocate);		}

  };

 return(0); 
}

int collect_transposition_truncate(collectbuffer *cb,lookupkey *lk,off_t preallocate)
{
 int			i;
 lookupkey		*wl;

 if ( lk==NULL || lk->lookups==NULL )	return(0);

 for ( i=0 ; i<(1<<HASH_BITS) ; i++ )
  {	wl=&lk->lookups[i];
	if ( wl->lookups != NULL )
		collect_transposition_truncate(cb,wl,preallocate);	
#ifndef	HOST_WIN32
	else if ( wl->value.key != NULL )
	 {	key	*wk;
		wk=&wl->value;
		if ( (off_t)wk->id < preallocate )
			(void)truncate(wk->key,(off_t)wk->id);
	 }
#endif
  }

 return(0);
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int do_transposition_normal(char **infiles,int ninfile,int colbase,
	char *basename,size_t maxmem,off_t preallocate,
	char **additionalcomments,int argc,char **argv)
{
 collectbuffer	cb; 
 FILE		*fr;
 int		i;
 
 collect_buffer_reset(&cb);
 cb.lk.lookups=NULL;
 cb.lk.value.key=NULL;
 cb.lk.value.id=0;

 cb.is_comment=is_comment;
 cb.comments=additionalcomments;
 cb.argc=argc;
 cb.argv=argv;

 for ( i=0 ; i<ninfile ; i++ )
  {	fr=fopenread(infiles[i]);
	if ( fr==NULL )
	 {	fprint_warning("unable to open file '%s', skipped",infiles[i]);
		continue;
	 }
	collect_read_file(&cb,i,fr,colbase,basename,maxmem,preallocate);
	fcloseread(fr);
  }

 collect_flush(&cb,basename,preallocate);

 if ( preallocate>0 )
	collect_transposition_truncate(&cb,&cb.lk,preallocate);

 lookup_free_node(&cb.lk);

 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int compare_line_by_keyoffset(const void *vl1,const void *vl2)
{
 char	*l1=(char *)vl1;
 char	*l2=(char *)vl2;
 int	k1=(*(int *)l1)+sizeof(int);
 int	k2=(*(int *)l2)+sizeof(int);
 return(strcmp(l1+k1,l2+k2));
}

int do_transposition_cached(char **infiles,int ninfile,int colbase,
	char *basename,size_t maxmem,char *tmpdir,int linesize)
{
 FILE		*fr,*fw;
 int		i,k,l,ln,key;
 char		*tmptemplate;
 int		fh;
 char		*buff,*memcache,*p;
 int		mcline,ncline,ttline,ccline;
 char		*filename,*prevfile;

 tmptemplate=NULL; 
 strappendf(&tmptemplate,"%s/grcollect.XXXXXX",tmpdir);

/* fprintf(stderr,"do_transposition_cached(): template=%s\n",tmptemplate); */

 fh=-1;

 memcache=NULL;
 mcline=(int)(maxmem/linesize);
 ncline=0;
 ttline=0;
 
 for ( i=0 ; i<ninfile ; i++ )
  {	fr=fopenread(infiles[i]);
	if ( fr==NULL )
	 {	fprint_warning("unable to open file '%s', skipped",infiles[i]);
		continue;
	 }
	ln=0;
	while ( ! feof(fr) )
	 {	ln++;
		buff=freadline(fr);
		if ( buff==NULL )	break;
		remove_newlines_and_comments(buff);
		l=strlen(buff);
		if ( l+1>linesize-sizeof(int) )
		 {	fprint_warning("line %d in file '%s' is too long (%d rather than %d), line is truncated",ln,infiles[i],l,(int)(linesize-sizeof(int)));
			l=linesize-sizeof(int)-1;
			buff[l]=0;
		 }
		key=get_column_position(buff,colbase);
		if ( key<0 )
		 {	free(buff);
			continue;
		 }
		memcache=(char *)realloc(memcache,(size_t)(ncline+1)*(size_t)linesize);
		p=memcache+(size_t)ncline*(size_t)linesize;
		*(int *)p=key;
		memcpy(p+sizeof(int),buff,linesize-sizeof(int));
		free(buff);
		ncline++;
		ttline++;
		if ( ncline >= mcline )
		 {	if ( fh<0 )
			 {	
#ifndef	HOST_WIN32
				fh=mkstemp(tmptemplate);
#else
				mktemp(tmptemplate);
				fh=open(tmptemplate,O_CREAT|O_TRUNC|O_RDWR);
#endif
				if ( fh<0 )
				 {	fprint_error("unable to create temporary file using '%s'",tmptemplate);
					exit(1);
				 }
			 }
			qsort(memcache,ncline,linesize,compare_line_by_keyoffset);
			for ( k=0 ; k<ncline ; k++ )
			 {	(void)write(fh,memcache+(size_t)k*(size_t)linesize,(size_t)linesize);
			 }
			free(memcache);
			ncline=0;
			memcache=NULL;
		 }
	 }
 	fcloseread(fr);
  }

 if ( fh>=0 && ncline>0 )
  {	qsort(memcache,ncline,linesize,compare_line_by_keyoffset);
	for ( k=0 ; k<ncline ; k++ )
	 {	(void)write(fh,memcache+(size_t)k*(size_t)linesize,(size_t)linesize);
 	 }
	free(memcache);
	memcache=NULL;
	ncline=0;
  }

 /* all data is in memory: */
 if ( ttline>0 )
  {	if ( fh<0 )
	 	qsort(memcache,ttline,linesize,compare_line_by_keyoffset);
	else 
	 {	cache	ch;
		size_t	blocksize,membcount;
		int	multip;
	
		blocksize=cache_blocksize(linesize);
		membcount=maxmem/blocksize;
		if ( membcount<4 )	membcount=4;
		multip=1;
		while ( membcount>16384 )
		 {	multip *= 2;
			membcount /= 2;
		 }

		cache_init(&ch,linesize,(off_t)ttline,multip,(int)membcount,fh,CACHE_READ_AND_WRITE);
		cache_sort(&ch,compare_line_by_keyoffset);
		cache_finalize(&ch);
	 }
  }

 prevfile=NULL;
 fw=NULL;
 ncline=0;
 ccline=0;
 for ( i=0 ; i<ttline ; i++ )
  {	if ( fh<0 )
	 {	p=memcache+(size_t)i*(size_t)linesize;		}
	else
	 {	if ( ncline<=0 )
		 {	ncline=(mcline<ttline-i?mcline:ttline-i);
			memcache=realloc(memcache,(size_t)ncline*(size_t)linesize);
			for ( k=0 ; k<ncline ; k++ )
			 {	(void)read(fh,memcache+(size_t)k*(size_t)linesize,(size_t)linesize);
			 }
			ccline=0;
		 }
		p=memcache+(size_t)ccline*(size_t)linesize;
		ccline++;
	 }

	key=*(int *)p;
	filename=expand_basename(basename,p+sizeof(int)+key);

	if ( ( prevfile != NULL && strcmp(prevfile,filename) ) || prevfile==NULL || fw==NULL )
	 {	if ( fw != NULL )	fclose(fw);
		fw=fopen(filename,"wb");
	 }
	if ( fw != NULL )	fprintf(fw,"%s\n",p+sizeof(int));

	if ( prevfile != NULL )		free(prevfile);
	prevfile=filename;

	if ( fh>=0 && ccline>=ncline )
		ncline=0;

  }

 if ( prevfile != NULL )	free(prevfile);
 if ( fw != NULL )		fclose(fw);

 if ( memcache != NULL )
  {	free(memcache);
	memcache=NULL;
	ncline=0;
  }

 if ( fh>=0 )
  {	close(fh);
	unlink(tmptemplate);
	free(tmptemplate);
  }

 return(0);
}


/*****************************************************************************/

int collect_stat_reset(collectstat *cs)
{
 memset(&cs->lk,0,sizeof(lookupkey));
 cs->nfield=0;
 cs->afield=0;
 cs->fields=NULL;
 return(0);
}
int collect_stat_flush_tempfile(collectstat *cs)
{
 size_t	cw,cp,cc;

 cw=cs->ntmp*(size_t)(sizeof(tmpvaluecommon)+sizeof(double)*cs->nstat);
 for ( cp=0 ; cw>0 ; )
  {	cc=cw;
	if ( cc>16*1024*1024 )	cc=16*1024*1024;
	(void)write(cs->fh,(void *)((char *)cs->tmps+cp),cc);
	cp+=cc;
	cw-=cc;
  }

 cs->wtmp += (off_t)cs->ntmp;
 free(cs->tmps);

 cs->atmp=(size_t)0;
 cs->ntmp=(size_t)0;
 cs->tmps=NULL;

 return(0);
}

int collect_stat_read_file(collectstat *cs,FILE *fr,int colbase,column *colstats)
{
 int			mxc,i,n;
 char			*rbuff,**cmd,*kstr;
 double			*values;
 key			*wk;
 statfieldcommon	*wsc;
 statfielddata		*wsd;
 tmpvaluecommon		*wtc;
 double			*wtd;

 mxc=colbase;
 for ( i=0 ; i<cs->nstat ; i++ )
  {	if ( colstats[i].index>mxc )
		mxc=colstats[i].index;
  }
 mxc++;

 values=(double *)malloc(sizeof(double)*cs->nstat);

 while ( ! feof(fr) )
  {	rbuff=freadline(fr);
	if ( rbuff==NULL )	break;
	remove_newlines_and_comments(rbuff);
	cmd=tokenize_spaces_dyn(rbuff);
	if ( cmd==NULL || cmd[0]==NULL )
	 {	if ( cmd != NULL )	free(cmd);
		free(rbuff);
		continue;
	 }
	for ( n=0 ; cmd[n] != NULL ; )	n++;
	if ( n<mxc )
	 {	free(cmd);
		free(rbuff);
		continue;
	 }
	for ( i=0 ; i<cs->nstat ; i++ )
	 {	if ( sscanf(cmd[colstats[i].index],"%lg",&values[i])<1 )
			break;
	 }
	if ( i<cs->nstat )
	 {	free(cmd);
		free(rbuff);
		continue;
	 }

	kstr=cmd[colbase];
	wk=lookup_search_or_add_key(&cs->lk,kstr,0);
	if ( wk->id<=0 )
	 {	cs->nfield++;
		wk->id=cs->nfield;
		if ( cs->nfield>cs->afield )
		 {	cs->afield+=64;
			cs->fields=realloc(cs->fields,(sizeof(statfieldcommon)+sizeof(statfielddata)*cs->nstat)*cs->afield);
		 }
		wsc=(statfieldcommon *)((char *)cs->fields+(cs->nfield-1)*(sizeof(statfieldcommon)+sizeof(statfielddata)*cs->nstat));
		wsd=(statfielddata *)((char *)wsc+sizeof(statfieldcommon));
		wsc->id=wk->id;
		wsc->key=wk->key;
		wsc->count=1;
		for ( i=0 ; i<cs->nstat ; i++ )
		 {	wsd[i].min=wsd[i].max=values[i];
			wsd[i].sum=values[i];
			wsd[i].sum2=values[i]*values[i];
		 }
	 }
	else
	 {	wsc=(statfieldcommon *)((char *)cs->fields+(wk->id-1)*(sizeof(statfieldcommon)+sizeof(statfielddata)*cs->nstat));
		wsd=(statfielddata *)((char *)wsc+sizeof(statfieldcommon));
		for ( i=0 ; i<cs->nstat ; i++ )
		 {	if ( values[i]<wsd[i].min )	wsd[i].min=values[i];
			if ( values[i]>wsd[i].max )	wsd[i].max=values[i];
			wsd[i].sum+=values[i];
			wsd[i].sum2+=values[i]*values[i];
		 }
		wsc->count++;
	 }

	if ( cs->tmptemplate != NULL )
	 {	if ( cs->ntmp>=cs->atmp )
		 {	cs->atmp+=(size_t)64;
			cs->tmps=realloc(cs->tmps,(size_t)(sizeof(tmpvaluecommon)+sizeof(double)*cs->nstat)*cs->atmp);
		 }

		wtc=(tmpvaluecommon *)((char *)cs->tmps+(size_t)(sizeof(tmpvaluecommon)+sizeof(double)*cs->nstat)*cs->ntmp);
		wtd=(double *)((char *)wtc+sizeof(tmpvaluecommon));
		cs->ntmp++;
		wtc->id=wk->id;
		wtc->flag=0;
		for ( i=0 ; i<cs->nstat ; i++ )
		 {	wtd[i]=values[i];		}
		if ( cs->ntmp>=cs->mxtmp )
		 {	if ( cs->fh<0 )
			 {	
#ifndef	HOST_WIN32
				cs->fh=mkstemp(cs->tmptemplate);
#else
				mktemp(cs->tmptemplate);
				cs->fh=open(cs->tmptemplate,O_CREAT|O_TRUNC|O_RDWR);
#endif
				if ( cs->fh<0 )
				 {	fprint_error("unable to create temporary file using '%s'",cs->tmptemplate);
					exit(1);
				 }
			 }
			collect_stat_flush_tempfile(cs);
		 }
	 }

	free(cmd);
	free(rbuff);
  };

 return(0);
}

int collect_stat_temp_compare(const void *vt1,const void *vt2)
{
 tmpvaluecommon	*t1=(tmpvaluecommon *)vt1,
		*t2=(tmpvaluecommon *)vt2;
 if ( t1->id < t2->id )					return(-1);
 else if ( t1->id > t2->id )				return(+1);
 else if ( *(double *)(&t1[1]) < *(double *)(&t2[1]) )	return(-1);
 else if ( *(double *)(&t1[1]) > *(double *)(&t2[1]) )	return(+1);
 else							return(0);
}

int collect_stat_order_mem(collectstat *cs)
{
 int	recordsize;
 recordsize=sizeof(tmpvaluecommon)+sizeof(double)*cs->nstat;
 qsort(cs->tmps,cs->ntmp,recordsize,collect_stat_temp_compare);
 return(0);
}

int collect_stat_order_tmp(collectstat *cs)
{
 cache	ch;
 int	recordsize,multip;
 size_t	blocksize,membcount;

 recordsize=sizeof(tmpvaluecommon)+sizeof(double)*cs->nstat;

 blocksize=cache_blocksize(recordsize);
 membcount=cs->maxmem/blocksize;
 if ( membcount<4 )	membcount=4;
 multip=1;
 while ( membcount>16384 )
  {	multip *= 2;
	membcount /= 2;
  }

 cache_init(&ch,recordsize,cs->wtmp,multip,(int)membcount,cs->fh,CACHE_READ_AND_WRITE);
 cache_sort(&ch,collect_stat_temp_compare);
 cache_finalize(&ch);

 return(0);
}

int collect_stat_order(collectstat *cs)
{
 if ( cs->wtmp<=0 )	/* all data in memory, nothing written to tmpfile */
  {	collect_stat_order_mem(cs);	}
 else			/* sort tmpfile instead				  */
  {	collect_stat_order_tmp(cs);	}
 return(0);
}

double get_std_deviation(double *arr,int n,double m)
{
 int	i;
 double	s2,dev;

 if ( n<=0 || arr==NULL )	return(-1.0);
 
 s2=0.0;
 for ( i=0 ; i<n ; i++ )
  {	s2 += (arr[i]-m)*(arr[i]-m);		}
 s2/=(double)n;
 if ( s2>0.0 )	dev=sqrt(s2);
 else		dev=0.0;
 return(dev); 
}

double get_abs_deviation(double *arr,int n,double m)
{
 int	i;
 double	dev;

 if ( n<=0 || arr==NULL )	return(-1.0);
 
 dev=0.0;
 for ( i=0 ; i<n ; i++ )
  {	dev += fabs(arr[i]-m);			}
 dev/=(double)n;
 return(dev); 
}

int get_value_index(double *arr,int n,double m)
{
 int    f,k;
 f=0;
 while ( n )
  {	k=f+n/2;
	if ( arr[k]<m )	f+=n/2+1,n-=n/2+1;
	else            n=n/2;
  };
 return(f);
}

double get_med_deviation(double *arr,int n,double m)
{
 int	i;
 double	meddev,*wrr;

 wrr=(double *)malloc(sizeof(double)*n);
 for ( i=0 ; i<n ; i++ )
  {	wrr[i]=fabs(arr[i]-m);		}
 meddev=median(wrr,n);
 free(wrr);

 return(meddev);
}

int get_noncumulative_statistics(double *arr,int n,statfielddata *wsd,rejection *r)
{
 int	i,lo,nn,ilo,inn,ii;
 double	mean,mlo,mhi,mct,mdv,rmean;

 wsd->median=0.5*(arr[(n-1)/2]+arr[n/2]);

 mean=wsd->sum/(double)n;

 wsd->mediandev=get_std_deviation(arr,n,wsd->median);
 wsd->medianmeddev=get_med_deviation(arr,n,wsd->median);
 wsd->meanmeddev=get_med_deviation(arr,n,mean);
 wsd->mode=2.5*wsd->median-1.5*mean;
 wsd->modedev=get_std_deviation(arr,n,wsd->mode);
 wsd->modemeddev=get_med_deviation(arr,n,wsd->mode);

 wsd->rcount=n;
 wsd->rsum=wsd->sum;
 wsd->rsum2=wsd->sum2;
 wsd->rmin=wsd->min;
 wsd->rmax=wsd->max;
 wsd->rmedian=wsd->median;
 wsd->rmediandev=wsd->mediandev;
 wsd->rmeanmeddev=wsd->meanmeddev;
 wsd->rmedianmeddev=wsd->medianmeddev;
 wsd->rmode=wsd->mode;
 wsd->rmodedev=wsd->modedev;
 wsd->rmodemeddev=wsd->modemeddev;

 lo=0;
 nn=n;

 for ( ii=0 ; r != NULL && (r->center & REJECTION_CENTER) && 
 (ii<r->iterations || r->iterations<0) ; ii++ )
  {	rmean=wsd->rsum/(double)wsd->rcount;
	switch ( r->center & REJECTION_CENTER )
	 {   case REJECTION_C_MEAN:
		mct=rmean;
		switch ( r->center & REJECTION_MAGNITUDE )
		 {   case REJECTION_M_STDDEV:
			mdv=wsd->rsum2/(double)n-rmean*rmean;
			if ( mdv>0 )	mdv=sqrt(mdv);
			else		mdv=0.0;
			break;
		     case REJECTION_M_MEDDEV:
			mdv=wsd->rmeanmeddev;
			break;
		     case REJECTION_M_ABSOLUTE:
			mdv=1.0;
			break;
		     default:
			mdv=0.0;
			break;
		 }
		break;
	     case REJECTION_C_MEDIAN:
		mct=wsd->rmedian;
		switch ( r->center & REJECTION_MAGNITUDE )
		 {   case REJECTION_M_STDDEV:
			mdv=wsd->rmediandev;
			break;
		     case REJECTION_M_MEDDEV:
			mdv=wsd->rmedianmeddev;
			break;
		     case REJECTION_M_ABSOLUTE:
			mdv=1.0;
			break;
		     default:
			mdv=0.0;
			break;
		 }
		break;
	     case REJECTION_C_MODE:
		mct=wsd->rmode;
		switch ( r->center & REJECTION_MAGNITUDE )
		 {   case REJECTION_M_STDDEV:
			mdv=wsd->rmodedev;
			break;
		     case REJECTION_M_MEDDEV:
			mdv=wsd->rmodemeddev;
			break;
		     case REJECTION_M_ABSOLUTE:
			mdv=1.0;
			break;
		     default:
			mdv=0.0;
			break;
		 }
		break;
	     default:
		mct=0.0;
		mdv=0.0;
		break;
	 }
	if ( mdv>0.0 )
	 {	mlo=mct-mdv*r->magnitude;
		mhi=mct+mdv*r->magnitude;
		ilo=get_value_index(arr,n,mlo);
		inn=get_value_index(arr,n,mhi)-ilo;
	 }
	else
	 {	ilo=lo;
		inn=nn;
	 }
	/*
	fprintf(stderr,"mct=%12g mdv=%12g ilo=%3d inn=%3d\n",mct,mdv,ilo,inn);
	*/
	if ( inn <= 0 )
	 {	wsd->rcount=0;
		wsd->rmin=wsd->rmax=0.0;
		wsd->rsum=wsd->rsum2=0.0;
		wsd->rmedian=0.0;
		wsd->rmediandev=0.0;
		wsd->rmedianmeddev=0.0;
		wsd->rmeanmeddev=0.0;
		wsd->rmode=0.0;
		wsd->rmodedev=0.0;
		wsd->rmodemeddev=0.0;
		break;
	 }
	else if ( ilo>lo || inn<nn )
	 {	lo=ilo;
		nn=inn;

		wsd->rcount=nn;
		wsd->rmin=arr[lo];
		wsd->rmax=arr[lo+nn-1];
		wsd->rsum=wsd->rsum2=0.0;
		for ( i=0 ; i<nn ; i++ )
		 {	wsd->rsum+=arr[lo+i];
			wsd->rsum2+=arr[lo+i]*arr[lo+i];
		 }
		rmean=wsd->rsum/(double)wsd->rcount;
		wsd->rmedian=0.5*(arr[lo+(nn-1)/2]+arr[lo+nn/2]);

		wsd->rmediandev=get_std_deviation(arr+lo,nn,wsd->rmedian);
		wsd->rmedianmeddev=get_med_deviation(arr+lo,nn,wsd->rmedian);
		wsd->rmeanmeddev=get_med_deviation(arr+lo,nn,rmean);
		wsd->rmode=2.5*wsd->rmedian-1.5*rmean;
		wsd->rmodedev=get_std_deviation(arr+lo,nn,wsd->rmode);
		wsd->rmodemeddev=get_med_deviation(arr+lo,nn,wsd->rmode);
	 }
	else
		break;
  }

 return(0);
}

int collect_stat_non_cumulative(collectstat *cs)
{
 int			i,j,k,aa,n;
 size_t			offset;
 statfieldcommon	*wsc;
 statfielddata		*wsd;
 tmpvaluecommon		*wtc;
 double			*wtd;
 double			*arr;
 cache			ch;
 int			recordsize;

 aa=0;
 arr=NULL;

 offset=(size_t)0;

 recordsize=sizeof(tmpvaluecommon)+sizeof(double)*cs->nstat;
 if ( cs->fh>=0 )
  {	lseek(cs->fh,(off_t)0,SEEK_SET);
	cache_init(&ch,recordsize,cs->wtmp,1,4,cs->fh,CACHE_READONLY);
  }

 for ( k=0 ; k<cs->nfield ; k++ )
  {	wsc=(statfieldcommon *)((char *)cs->fields+(size_t)k*(size_t)(sizeof(statfieldcommon)+sizeof(statfielddata)*cs->nstat));
	wsd=(statfielddata *)((char *)wsc+sizeof(statfieldcommon));

	n=wsc->count;
	if ( n>aa )
	 {	aa=n;
		arr=(double *)realloc(arr,sizeof(double)*aa);
	 }

	for ( i=0 ; i<cs->nstat ; i++ )
	 {	if ( cs->wtmp<=0 )
		 {	for ( j=0 ; j<n ; j++ )
			 {	wtc=(tmpvaluecommon *)((char *)cs->tmps+(offset+(size_t)j)*(size_t)recordsize);
				wtd=(double *)((char *)wtc+sizeof(tmpvaluecommon));
				arr[j]=wtd[i];
			 }
		 }
		else
		 {	for ( j=0 ; j<n ; j++ )
			 {	wtc=cache_read_record(&ch,(off_t)(offset+(size_t)j));
				wtd=(double *)((char *)wtc+sizeof(tmpvaluecommon));
				arr[j]=wtd[i];
			 }
		 }
		if ( i>0 )	median(arr,n);	/* only the first column has been ordered */

		get_noncumulative_statistics(arr,n,&wsd[i],&cs->colstats[i].r);

	 }
	offset+=(size_t)wsc->count;
  }

 if ( arr != NULL )	free(arr);

 if ( cs->fh>=0 )
	cache_finalize(&ch);

 return(0);
}

int collect_stat_write_stat_one(FILE *fw,statfieldcommon *wsc,statfielddata *wsd,int *statlist)
{
 double	mean,m2,sigma,rmean,rsigma;
 int	j,s;

 mean=wsd->sum /(double)wsc->count;
 m2  =wsd->sum2/(double)wsc->count;
 sigma=m2-mean*mean;
 if ( sigma>0.0 )	sigma=sqrt(sigma);
 else			sigma=0.0;

 if ( wsd->rcount>0 )
  {	rmean=wsd->rsum/(double)wsd->rcount;
 	m2=wsd->rsum2/(double)wsd->rcount;
	rsigma=m2-rmean*rmean;
	if ( rsigma>0.0 )	rsigma=sqrt(rsigma);
	else			rsigma=0.0;
  }
 else
  {	rmean=0.0;
	rsigma=0.0;
  }
	
 for ( j=0 ; statlist[j]>0 ; j++ )
  {	s=statlist[j];
	switch ( s )
	 {   case STAT_COUNT:
		fprintf(fw,"%6d ",wsc->count);
		break;
	     case STAT_MEAN:
		fprintf(fw,"%12g ",mean);
		break;
	     case STAT_MEANDEV:
		fprintf(fw,"%12g ",sigma);
		break;
	     case STAT_MIN:
		fprintf(fw,"%12g ",wsd->min);
		break;
	     case STAT_MAX:
		fprintf(fw,"%12g ",wsd->max);
		break;
	     case STAT_SUM:
		fprintf(fw,"%12g ",wsd->sum);
		break;
	     case STAT_SUM2:
		fprintf(fw,"%12g ",wsd->sum2);
		break;

	     case STAT_MEDIAN:
		fprintf(fw,"%12g ",wsd->median);
		break;
	     case STAT_MEDIANDEV:
		fprintf(fw,"%12g ",wsd->mediandev);
		break;
	     case STAT_MEDIANMEDDEV:
		fprintf(fw,"%12g ",wsd->medianmeddev);
		break;
	     case STAT_MEANMEDDEV:
		fprintf(fw,"%12g ",wsd->meanmeddev);
		break;
	     case STAT_MODE:
		fprintf(fw,"%12g ",wsd->mode);
		break;
	     case STAT_MODEDEV:
		fprintf(fw,"%12g ",wsd->modedev);
		break;
	     case STAT_MODEMEDDEV:
		fprintf(fw,"%12g ",wsd->modemeddev);
		break;

	     case STAT_RCOUNT:
		fprintf(fw,"%6d ",wsd->rcount);
		break;
	     case STAT_RMEAN:
		fprintf(fw,"%12g ",rmean);
		break;
	     case STAT_RMEANDEV:
		fprintf(fw,"%12g ",rsigma);
		break;
	     case STAT_RMIN:
		fprintf(fw,"%12g ",wsd->rmin);
		break;
	     case STAT_RMAX:
		fprintf(fw,"%12g ",wsd->rmax);
		break;
	     case STAT_RSUM:
		fprintf(fw,"%12g ",wsd->rsum);
		break;
	     case STAT_RSUM2:
		fprintf(fw,"%12g ",wsd->rsum2);
		break;
	     case STAT_RMEDIAN:
		fprintf(fw,"%12g ",wsd->rmedian);
		break;
	     case STAT_RMEDIANDEV:
		fprintf(fw,"%12g ",wsd->rmediandev);
		break;
	     case STAT_RMEDIANMEDDEV:
		fprintf(fw,"%12g ",wsd->rmedianmeddev);
		break;
	     case STAT_RMEANMEDDEV:
		fprintf(fw,"%12g ",wsd->rmeanmeddev);
		break;
	     case STAT_RMODE:
		fprintf(fw,"%12g ",wsd->rmode);
		break;
	     case STAT_RMODEDEV:
		fprintf(fw,"%12g ",wsd->rmodedev);
		break;
	     case STAT_RMODEMEDDEV:
		fprintf(fw,"%12g ",wsd->rmodemeddev);
		break;

	 }
  }
 return(0);
}

int collect_stat_write_stat_line(FILE *fw,
	statfieldcommon *wsc,statfielddata *wsd,int nstat,int *statlist)
{
 int	i;
 fprintf(fw,"%s\t",wsc->key);
 for ( i=0 ; i<nstat ; i++ )
  {	collect_stat_write_stat_one(fw,wsc,&wsd[i],statlist);		}
 fprintf(fw,"\n");
 return(0);
}

int collect_stat_write_node(FILE *fw,collectstat *cs,lookupkey *lk,int *statlist)
{
 int			i;
 key			*wk;
 lookupkey		*wl;
 statfieldcommon	*wsc;
 statfielddata		*wsd;

 if ( lk==NULL || lk->lookups==NULL )	return(0);

 for ( i=0 ; i<(1<<HASH_BITS) ; i++ )
  {	wl=&lk->lookups[i];
	if ( wl->lookups != NULL )
		collect_stat_write_node(fw,cs,wl,statlist);	
	else if ( wl->value.key != NULL )
	 {	wk=&wl->value;
		wsc=(statfieldcommon *)((char *)cs->fields+(wk->id-1)*(sizeof(statfieldcommon)+sizeof(statfielddata)*cs->nstat));
		wsd=(statfielddata *)((char *)wsc+sizeof(statfieldcommon));
		collect_stat_write_stat_line(fw,wsc,wsd,cs->nstat,statlist);
	 }
  }

 return(0);
}

int collect_stat_write(FILE *fw,collectstat *cs,int *statlist)
{
 collect_stat_write_node(fw,cs,&cs->lk,statlist);
 return(0);
}

int collect_stat_free(collectstat *cs)
{
 lookup_free_node(&cs->lk);
 return(0);
}


int do_collective_statistics(char **infiles,int ninfile,FILE *fw,
	int colbase,column *colstats,int ncolstat,
	int *statlist,size_t maxmem,char *tmpdir)
{
 int		i;
 int		is_non_cumulative;
 FILE		*fr;
 collectstat	cs;
 size_t		mx;

 if ( statlist==NULL )
	return(1);

 is_non_cumulative=0;
 for ( i=0 ; statlist[i]>0 && ! is_non_cumulative ; i++ )
  {	if ( statlist[i]>STAT_PURE_CUMULATIVE )
		is_non_cumulative=1;
  }

 collect_stat_reset(&cs);

 if ( is_non_cumulative )
  {	if ( tmpdir==NULL || strlen(tmpdir) <= 0 )
		cs.tmptemplate=strdup("./grcollect.XXXXXX");
	else
	 {	i=strlen(tmpdir);
		cs.tmptemplate=(char *)malloc(i+32);
		strcpy(cs.tmptemplate,tmpdir);
		if ( tmpdir[i-1] !='/' )	strcat(cs.tmptemplate,"/");
		strcat(cs.tmptemplate,"grcollect.XXXXXX");
	 }
  }
 else
	cs.tmptemplate=NULL;

 cs.tmps=NULL;
 cs.ntmp=(size_t)0;
 cs.atmp=(size_t)0;
 cs.wtmp=(off_t)0;
 mx=(maxmem/(size_t)(sizeof(tmpvaluecommon)+ncolstat*sizeof(double))+(size_t)63)/(size_t)64;
 cs.mxtmp=(size_t)64*mx;
 cs.fh=-1;
 cs.maxmem=maxmem;

 cs.nstat=ncolstat;
 cs.colstats=colstats;
 
 for ( i=0 ; i<ninfile ; i++ )
  {	fr=fopenread(infiles[i]);
	if ( fr==NULL )
	 {	fprint_warning("unable to open file '%s', skipped",infiles[i]);
		continue;
	 }
	collect_stat_read_file(&cs,fr,colbase,colstats);
	fcloseread(fr);
  }

 if ( is_non_cumulative )
  {	if ( cs.wtmp>(size_t)0 && cs.fh>=0 )
	 {	collect_stat_flush_tempfile(&cs);		}
	collect_stat_order(&cs);
	collect_stat_non_cumulative(&cs);
	if ( cs.fh>=0 )
	 {	close(cs.fh);
		unlink(cs.tmptemplate);
	 }
  }

 collect_stat_write(fw,&cs,statlist);
 collect_stat_free(&cs);
 
 return(0);
}

/*****************************************************************************/

int main(int argc,char *argv[])
{
 int		i,is_help;
 int		colbase;
 char		**infiles,*basename,*extension,*tmpdir,
		*prefix,*maxmemstr,*outfile,*statstr,*colstatstr,*preallocarg,
		**rejectionlist,**additionalcomments;
 int		ninfile,*statlist;
 column		*colstats;
 int		ncolstat;
 int		linesize;
 off_t		preallocate;

 size_t		maxmem;

 progbasename=strrchr(argv[0],'/');
 if ( progbasename != NULL )	progbasename++;
 else				progbasename=argv[0];

 is_comment=0;is_verbose=is_help=0;
 additionalcomments=NULL;
 rejectionlist=NULL;

 infiles=NULL;
 outfile=basename=extension=prefix=NULL;
 statstr=maxmemstr=NULL;
 tmpdir=NULL;

 colbase=1;
 colstatstr="2";

 linesize=240;
 preallocarg=NULL;

 i=scanarg(argc,argv,SCANARG_ALLOW_FLAGS,
	"--version:%NS-1f%q",&is_help,
        "--version-short|--short-version:%NS-2f%q",&is_help,
	"-C|--comment:%i",&is_comment,
	"-S|--additional-comment:%Dt",&additionalcomments,
	"-h|--help|--short-help|--help-short:%f%q",&is_help,
	"--long-help|--help-long:%SN2f%q",&is_help,
	"--mediawiki-help|--help-mediawiki|--wiki-help|--help-wiki:%SN3f%q",&is_help,
	"-b|--basename:%s",&basename,
	"-o|--output:%s",&outfile,
	"-x|--extension:%s",&extension,
	"-p|--prefix:%s",&prefix,
	"-a|--pre-allocate:%s",&preallocarg,
	"-m|--max-mem|--max-memory:%s",&maxmemstr,
	"-s|--stat:%s",&statstr,
	"-t|--tmpdir|--tempdir|--temporary-directory:%s",&tmpdir,
	"-r|--rejection:%Dt",&rejectionlist,
	"-l|--linesize|--line-size:%d",&linesize,
	"-c|--col-base:%d",&colbase,
	"-d|--col-stat:%s",&colstatstr,
	"--comment:%f",&is_comment,"(C):%f",&is_comment,
	"--verbose:%i",&is_verbose,"(V):%i",&is_verbose,
	"-:%Dl",&infiles,
	"-*|+*:%e",
	"*:%Dl",&infiles,
	NULL);

 if ( i )		
  {	fprint_error("invalid command line argument near '%s'",argv[i]);
	return(1);
  }
 else if ( is_help<0 )
  {	fprint_generic_version(stdout,argv[0],"grcollect",FITSH_GRCOLLECT_VERSION,is_help);
	return(0);
  }
 else if ( 1<is_help )
  {	fprint_grcollect_long_help(stdout,2<is_help);
	return(0);
  }
 else if ( is_help )
  {	fprint_grcollect_usage(stdout);
	return(0);
  }

 if ( maxmemstr != NULL )
  {     maxmem=parse_max_memory_string(maxmemstr);
        if ( maxmem<=1 )
         {      fprint_error("invalid maximum memory specification '%s'",maxmemstr);
                return(1);
         }
  }
 else
	maxmem=0;

 if ( maxmem <= 0 )
	maxmem=(size_t)(DEFAULT_MAX_MEMORY*1024*1024);

 if ( statstr != NULL )
  {	statlist=get_stat_list(statstr);
	if ( statlist==NULL )
	 {	fprint_error("invalid list of statistics");
		return(1);
	 }
  }
 else
	statlist=statlist_default;

 /* fprintf(stderr,"maxmem=%zd\n",maxmem); */

 colbase--;
 if ( colbase<0 )	
  {	fprint_error("invalid column index specification");
	return(1);
  }

 if ( colstatstr != NULL )
  {	char	**cmd;
	int	c;
	cmd=tokenize_char_dyn_wwt(colstatstr,',',0);
	colstats=NULL;
	for ( ncolstat=0 ; cmd[ncolstat] != NULL && cmd != NULL ; )
	 {	colstats=(column *)realloc(colstats,sizeof(column)*(ncolstat+1));
		if ( sscanf(cmd[ncolstat],"%d",&c)<1 || c<=0 )
		 {	fprint_error("invalid column index specification");
			return(1);
		 }
		colstats[ncolstat].index=c-1;
		colstats[ncolstat].r.center=0;
		colstats[ncolstat].r.magnitude=0.0;
		ncolstat++;
	 }
	if ( cmd != NULL )	free(cmd);
  }
 else
  {	ncolstat=0;
	colstats=NULL;
  }

 if ( rejectionlist != NULL )
  {	int		index,magtype,j,k;
	rejection	r;
	
	for ( i=0 ; rejectionlist[i] != NULL ; i++ )
	 {	index=magtype=0;
		r.center=0;
		r.magnitude=0.0;
		r.iterations=1;
		k=scanpar(rejectionlist[i],SCANPAR_DEFAULT,
			"column:%d",	&index,
			"mean:"		SNf(REJECTION_C_MEAN),    &r.center,
			"median:"	SNf(REJECTION_C_MEDIAN),  &r.center,
			"mode:"		SNf(REJECTION_C_MEAN),	  &r.center,
			"stddev:%g"	SNf(REJECTION_M_STDDEV),  &r.magnitude,&magtype,
			"meddev:%g"	SNf(REJECTION_M_MEDDEV),  &r.magnitude,&magtype,
			"absolute:%g"	SNf(REJECTION_M_ABSOLUTE),&r.magnitude,&magtype,
			"iterations:%d",&r.iterations,
			NULL);
		for ( j=0 ; j<ncolstat ; j++ )
		 {	if ( colstats[j].index==index-1 )
				break;
		 }
		if ( r.magnitude <= 0.0 )
		 {	r.center=magtype=0;
			r.magnitude=0;
		 }
		if ( k || ( j>=ncolstat ) )
		 {	fprint_error("invalid rejection content (index=%d)",index);
			return(1);
		 }
		colstats[j].r.center=r.center|magtype;
		colstats[j].r.magnitude=r.magnitude;
		colstats[j].r.iterations=r.iterations;
	 }
  }
			

 if ( basename != NULL && (extension != NULL || prefix != NULL) )
  {	fprint_error("invalid combination of command line arguments");
	return(1);
  }
 else if ( extension != NULL || prefix != NULL )
  {	int	len;
	len=3;
	if ( extension != NULL )	len+=strlen(extension)+1;
	if ( prefix != NULL )		len+=strlen(prefix)+1;
	basename=malloc(len);
	if ( extension != NULL && prefix != NULL )
		sprintf(basename,"%s%%b.%s",prefix,extension);
	else if ( extension==NULL && prefix != NULL )
		sprintf(basename,"%s%%b",prefix);
	else if ( extension != NULL && prefix==NULL )
		sprintf(basename,"%%b.%s",extension);
  }

 if ( basename==NULL && outfile==NULL )
  {	fprint_error("neither basename nor output file name has been specified");
	return(1);
  }
 else if ( basename != NULL && outfile != NULL )
  {	fprint_error("both basename and output file name have been specified");
	return(1);
  }

 if ( infiles==NULL )	return(0);	/* no input files */
 for ( ninfile=0 ; infiles[ninfile] != NULL ; )	ninfile++;
 if ( ninfile<=0 )	return(0);	/* no input files */

 if ( linesize<0 )	linesize=240;
 linesize=(linesize+8+16)&(~0xf);
/* fprintf(stderr,"linesize=%d\n",linesize); */

 if ( basename != NULL )		/* do transposition */
  {	/* fprintf(stderr,"tmpdir=%s\n",tmpdir); */

	if ( preallocarg==NULL )
		preallocate=0;

	else
	 {	int	s1,s2;
		off_t	pagesize;

		i=sscanf(preallocarg,"%d,%d",&s1,&s2);
		if ( i<1 )
		 {	fprint_error("unexpected pre-allocation size definition '%s'",preallocarg);
			return(1);
		 }
		else if ( i<2 )
			preallocate=(off_t)s1;
		else
		 {	s1=(s1+15)&(~15);
			s2=(s2+15)&(~15);
			preallocate=((off_t)s1)*((off_t)s2);
		 }
#ifdef	HOST_WIN32
		pagesize=(off_t)4096;
#else
		pagesize=(off_t)getpagesize();
#endif
		preallocate=(preallocate+pagesize-1)&(~(pagesize-1));
	 }

	/*
	if ( tmpdir == NULL )
	*/

	do_transposition_normal(infiles,ninfile,
		colbase,basename,maxmem,preallocate,
		additionalcomments,argc,argv);

	/*
	else
		do_transposition_cached(infiles,ninfile,colbase,basename,maxmem,tmpdir,linesize);
	*/

	return(0);
  }

 if ( outfile != NULL )
  {	FILE	*fw;
	fw=fopenwrite(outfile);
	if ( fw==NULL )
	 {	fprint_error("unable to create output file '%s'",outfile);
		return(1);
	 }
	do_collective_statistics(infiles,ninfile,fw,colbase,
		colstats,ncolstat,statlist,maxmem,tmpdir);
	fclosewrite(fw);
	return(0);
  }

 return(0);
}

/*****************************************************************************/

