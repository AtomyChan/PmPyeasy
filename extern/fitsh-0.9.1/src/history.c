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

#include "fi.h"
#include "history.h"

/*****************************************************************************/

int fits_history_export_command_line(fits *img,char *prg,char *vrs,int argc,char *argv[])
{
 char	buff[64];
 time_t	t;
 struct	tm	tm;

 time(&t);
 localtime_r(&t,&tm);

 sprintf(buff,"%s: %s [%s@%s] %.4d.%.2d.%.2d %.2d:%.2d:%.2d (%s)",
	prg,vrs,FI_VERSION,FI_RELEASE_DATE,
	tm.tm_year+1900,tm.tm_mon+1,tm.tm_mday,
	tm.tm_hour,tm.tm_min,tm.tm_sec,
	(!daylight?tzname[0]:tzname[1]));
 fits_header_export_command_line(img,"FI_HSTRY",buff,"> ",argc,argv);

 return(0);
}

/*****************************************************************************/
                                            
