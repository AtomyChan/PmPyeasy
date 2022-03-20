/*****************************************************************************/
/* lfit.h								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* This file should be included in all external (dynamically linked) library */
/* extensions designed for the `lfit` program. This file provides the object */
/* `lfitfunction` which shall be used to bind the external functions to the  */
/* main `lfit` code. The modules must be compiled using the command	     */
/* gcc ... -fPIC -o module.so module.c ... -shared			     */
/* (with the additional compilation flags and object modules)		     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 1996, 2002, 2004-2005, 2006, 2007-2008, 2009;			     */
/* Pal, Andras (apal@szofi.net)						     */
/*****************************************************************************/

#ifndef	__LFIT_H_INCLUDED

#define	__LFIT_H_INCLUDED	1

/*****************************************************************************/

#define		LFITFUNCTION_NONDIFFERENTIABLE		0x00
#define		LFITFUNCTION_DIFFERENTIABLE		0x02
#define		LFITFUNCTION_LINEAR			0x03

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

typedef struct
 {	char	*name;			/* name of the given function	     */
	int	flags;			/* see LFITFUNCTION_* above	     */
	int	nvar;			/* number of adjustable variables    */
	int	nidv;			/* number of independent parameters  */
	int	(*function)		/* the callback function itself:     */
		(double *vars,		/*  - array of the adj. variables    */
		 double	*idvs,		/*  - array of the indep. param's    */
		 double	*ret,		/*  - return value		     */
		 double *diff);		/*  - partial derivatives	     */
 } lfitfunction;

/*****************************************************************************/

typedef struct
 {	char	*description;
 } psnfunctinfo;

typedef struct
 {	char		*name;
	int		type,major,minor;
	int		argnum,precedency,associativity;
	int		(*funct)(double *);
	short		*diff;
	char		*string;
	int		strength,affixation;
	psnfunctinfo	*info;
 } psnlfit;

/*****************************************************************************/

#endif

/*****************************************************************************/
                    
