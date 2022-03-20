/*****************************************************************************/
/* linear.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Demonstartion and/or test module for `lfit` dynamic extensions.	     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2008; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <lfit/lfit.h>

/*****************************************************************************/

/* these functions can be declared as static ones since the wrapper inside
 `lfit` reads only the exported `lfitfunction` arrays (here: "linear") */
static int function_line(double *vars,double *idvs,double *ret,double *diff);

lfitfunction linear [] =
 {	{ "line", LFITFUNCTION_LINEAR,	2, 1, function_line },
	{ NULL, 0, 0, 0, NULL }
 };

/*****************************************************************************/

static int function_line(double *vars,double *idvs,double *ret,double *diff)
{

 *ret = vars[0]*idvs[0]+vars[1];
 if ( diff != NULL )
  {	diff[0]=idvs[0];
	diff[1]=1.0;
  }
/*
 *ret = vars[0]*vars[2]+vars[1];
 if ( diff != NULL )
  {	diff[0]=vars[2];
	diff[1]=1.0;
	diff[2]=vars[0];
  }
*/
 return(0);
}

/*****************************************************************************/
                                  
