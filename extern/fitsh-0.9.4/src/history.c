/*****************************************************************************/
/* history.c								     */
/*****************************************************************************/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include <fits/fits.h>

#include "fitsh.h"
#include "history.h"

/*****************************************************************************/

int fits_history_export_command_line(fits *img,char *prg,char *vrs,int argc,char *argv[])
{
 char	buff[256];
 time_t	t;
 struct	tm	*tm;

 time(&t);
 tm=localtime(&t);

 sprintf(buff,"%s: %s [%s@%s] %.4d.%.2d.%.2d %.2d:%.2d:%.2d (%s)",
	prg,vrs,FITSH_VERSION,FITSH_RELEASE,
	tm->tm_year+1900,tm->tm_mon+1,tm->tm_mday,
	tm->tm_hour,tm->tm_min,tm->tm_sec,
	(!daylight?tzname[0]:tzname[1]));
 fits_header_export_command_line(img,"FI_HSTRY",buff,"> ",argc,argv);

 return(0);
}

/*****************************************************************************/
                                                                    
                             
