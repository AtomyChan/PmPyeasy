/*****************************************************************************/
/* spmatrix.h								     */
/*****************************************************************************/

/* spmatrix.h : a 2D sparse matrix implementation (domsa@konkoly.hu (2002)) */
/* written for the HAT Project */
/* spmatrix.h,v 5.5 2003/05/12 18:18:11 domsa Exp */

#ifndef __SPMATRIX_H_INCLUDED
#define __SPMATRIX_H_INCLUDED	1

typedef enum
 {	SpNoErr,
	SpMem,
	SpInvSub,
 	SpNull
 } sperr;

typedef struct
 {	void  *matrix;
 } spmatrix;

typedef void *  spnode;

/* spm_create(), spm_destroy():
   create and destroy a sparse matrix					     */
spmatrix *spm_create(void);
void    spm_destroy(spmatrix *);

/* spm_geterr(), spm_numofelems():
   read error status or number of elements				     */
sperr	spm_geterr(spmatrix *);
int	spm_numofelems(spmatrix *);

/* spm_setnode(), spm_incrnode(), spm_addval(), spm_getval():
   set and get the value of a node					     */
void	spm_setnode(spmatrix *,int row,int col,int value);
void	spm_incrnode(spmatrix *,int row,int col);
void	spm_addval(spmatrix *,int val,int row,int col);
int     spm_getval(spmatrix *,int row,int col);

/* get next non-zero node value (set spnode to NULL to start from the beg) */
spnode	spm_getnode(spmatrix *,int row,int col);
spnode	spm_getnextnode(spmatrix *,spnode last,spnode *witness);

/* get values from a node */
int	spm_getcol(spnode);
int	spm_getrow(spnode);
int	spm_getnodeval(spnode);
void	spm_setnodeval(spnode, int value);

#endif
                                                                   
                     
