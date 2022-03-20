/*****************************************************************************/
/* fitsh.h								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Global configuration stuff for FITSH package...			     */
/*****************************************************************************/

#ifndef	__FITSH_H_INCLUDED
#define	__FITSH_H_INCLUDED	1

/*****************************************************************************/

#include "longhelp.h"
#include "../config.h"

/*****************************************************************************/

/* It's the exact value: */
#define		SIG_FWHM		(2.354820045030949382)
/* in the previous releases (<0.9c0) it was 2.35. */

/*****************************************************************************/

#define		MASK_OK			0x00	/* everything is ok...	     */
#define		MASK_CLEAR		MASK_OK	/* alias		     */

#define		MASK_FAULT		0x01	/* faulty pixel		     */
#define		MASK_HOT		0x02	/* hot (nonlinear) pixel     */
#define		MASK_COSMIC		0x04	/* hit by a cosmic particle  */
#define		MASK_OUTER		0x08	/* outer pixel		     */
#define		MASK_OVERSATURATED	0x10	/* oversaturated	     */
#define		MASK_LEAKED		0x20	/* leaked (during readout)   */
#define		MASK_SATURATED		(MASK_OVERSATURATED|MASK_LEAKED)
#define		MASK_INTERPOLATED	0x40	/* interpolated (not real)   */

#define		MASK_BAD		0x7F	/* any error		     */
#define		MASK_ALL		0x7F	/* alias		     */

#define		MASK_NOBACKGROUND	0x0100	/* indefinite background, not*/
						/* used in a (char **) mask  */

#define		MASK_NAN		MASK_FAULT

/*****************************************************************************/

typedef struct
 {	int	code;	
	char	*message;
 } errormessage;

/*****************************************************************************/

#define	DEFAULT_MAX_MEMORY		8	/* in megabytes	*/

/*****************************************************************************/

/* print_error(), print_error_cmdline():
   Prints the error with the specified code 'code' to the stream 'fw'.	     */
/*
int	print_error(FILE *fw,int code,...);
int	print_error_cmdline(FILE *fw,char *arg);
*/

/* print_warning():
   Prints the warning message with the code 'code' to the stream 'fw'.	     */
/*
int	print_warning(FILE *fw,int code,...);
*/

char	*get_progname(char *argv0);

/* fprint_generic_version():
   Print generic version information to the stream 'fw'. */
int	fprint_generic_version(FILE *fw,char *arg0,char *name,char *pv,int type);

int	fprint_generic_long_help(FILE *fw,int is_wiki,longhelp_entry *help,char *synopsis,char *description);

/* logmsg():
   Writes the sprintf-formatted message 'msg' to stderr, if 'flag' is true.  */
int	logmsg(int flag,char *msg,...);

/*****************************************************************************/

int	parse_mask_flags(char *list);

typedef struct
 {	int	bitpix;
	int	is_scale;
	double	bscale;
	double	bzero;
	int	nquantizebit;
 } fitsdataparam;

/* parse_fits_data_param():
   Parses the parameter string 'param' and sets the values in 'fdp'. If
   'param' is NULL or some of the parameters are not defined, the default
   parameters (bitpix=0, bscale=1.0 and bzero=0.0) are set.		     */
int	parse_fits_data_param(char *param,fitsdataparam *fdp);

/*****************************************************************************/

/* is_nasty_char(), is_any_nasty_char():
   These functions check whether is there any nasty characters (which are 
   preprocessed by the standard unix command shells) in 't' and 'buff'.      */
int	is_nasty_char(int t);
int	is_any_nasty_char(char *buff);

/*****************************************************************************/

size_t	parse_max_memory_string(char *maxmemstr);

/*****************************************************************************/

#endif

/*****************************************************************************/
                             
