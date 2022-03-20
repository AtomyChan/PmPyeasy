/*****************************************************************************/
/* scanarg.h								     */
/*****************************************************************************/

#ifndef	__SCANARG_H_INCLUDED
#define	__SCANARG_H_INCLUDED	1

#define	SCANARG_DEFAULT		0x00
#define	SCANARG_DASHDASH_ARG	0x01
#define	SCANARG_DASHDASH_EXP	0x02
#define	SCANARG_DASHDASH	(SCANARG_DASHDASH_EXP|SCANARG_DASHDASH_ARG)
#define	SCANARG_WILDCARD_DASH	0x04
#define	SCANARG_NO_SKIP_FIRST	0x08
#define	SCANARG_DYNAMIC_LISTS	0x00	/* this is the default for %l and %t */
#define	SCANARG_STATIC_LISTS	0x10
#define	SCANARG_USE_STRCMP	0x20
#define	SCANARG_ALLOW_FLAGS	0x40

#define SNf(i)			SNf_str(i)
#define SNf_str(i)		"%SN" #i "f"

/* scanarg():
   Scans the (basically) command-line argument list (*argc,argv[]) by the
   rules listed in 'format's. The results are stored in the variables
   listed after the format string, like the scanf()-like functions. All
   patterns to be matched should be separated by the (list of) pointer(s)
   where the result(s) are going to be stored.				     */

int	scanarg(int argc,char **argv,int flags,...);

/******************************************************************************
The syntax of the 'format' strings:
	arg[|arg2[|...]]:{%[<optional-formats>][e|q|d|g|s|m|r|w|l|t|f|i|(...)]}
Here 'arg' can be any fnmatch-like expression or a string like '(x)' if
SCANARG_ALLOW_FLAGS is set.

    %e:	
	Error (stops argument scan processing and returns an error).
    %q:	
	Quit (stops argument processing and returns successfully).
    %c:
	Continue (do nothing).

    %[PZN]d: 
	Integer. P allows positive, Z allows zero and N allows negative
	number. If none of them are specified, all integer values are 
	accepted. For example, to accept only natural (non-negative 
	integer) numbers, write "%PZd". 

    %g:	
	Real (stored as double)

    %s:	
	String (save the pointer of the argument)

    %m:	
	String (duplicates the argument with malloc() and save this pointer)

    %[{C|S|M|L}}r:
	String. 
	(concates after the dynamical string, separate the strings optionally 
	with comma[C], whitespace[S], semicolon[M], colon[L] or without any 
	separator[simply %r])
    %w:
	The switch as a string.
	(save the pointer of the switch)
    %[{C|S|M|L}}a:
	The switch as a string. 
	(concates after the dynamical string, separate the strings optionally 
	with comma[C], whitespace[S], semicolon[M], colon[L] or without any 
	separator[simply %a])

    %[{D|S[<max>]}]l:
	List of the matched arguments.
	(D: dynamic, S: static with maximum <max> elements.)

    %[{D|S[<max>]}][M]t:
	List of the arguments of the switch.
	(D: dynamic, S: static with maximum <max> elements, M: more than
	one argument after the swith is accepted, until another switch...)

    %[N][S][<n>f]:
	Flag.
	(N: threat <n> as a number otherwise a bit index - then 0 is the LSB,
	S: set the value instead of or'ing)

    %[<n>]i:
	Incrementation (by <n>).

    %(<arg1>[,<arg2>[,...]]):
	Selection.

Examples:
    Usage: prog	[-n|--number <n>] [-k] [-o|--output <output>] 	=>
    scanarg(argc,argv,0,
	"-n|--number:%d",&n,
	"-k:%f",&k_flag,
	"-o|--output:%s",&outputname,
	NULL);

    Usage: prog [-i <input-files>] [-x <xmin>:<xmax>] [-y <ymin>:<ymax>] 
		[-o <output>]	=>
    char **list=NULL;
    scanarg(argc,argv,0,
	"-i:%Dl",&list,
	"-x:%g:%g",&xmin,&xmax,
	"-y:%g:%g",&ymin,&ymax,
	"-o:%s",&outputname,
	NULL);

    Usage: prog	[-l <expr1>[,<expr2>]] [-l <expr3>,...] [...]
		[-v|--verbose] [--opt-*] [-m|--method {alpha|betha|gamma}]  =>
    char *explist=NULL,**opt_args=NULL;
    int	 verbose_level=0,method=-1;
    scanarg(argc,argv,0,
	"-l:%Cr",&explist,
	"-v|--verbose:%1i",&verbose_level,
	"--opt-*:%t",&opt_args,
	"-m|--method:%(alpha,beta,gamma)",&method,
	NULL);

******************************************************************************/

#define	SCANPAR_DEFAULT		0x00
#define	SCANPAR_USE_STRCMP	0x20

int	scanpar(char *par,int flags,...);

/*****************************************************************************/

#define	SCANFLAG_DEFAULT	0x00
#define	SCANFLAG_ALLOW_NEGATE	0x01
#define	SCANFLAG_ALLOW_RESET	0x02

int	scanflag(char *par,int flags,...);

/*****************************************************************************/

#endif
                     
