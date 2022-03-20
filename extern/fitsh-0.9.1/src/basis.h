/*****************************************************************************/
/* basis.h								     */
/*****************************************************************************/

#ifndef	__BASIS_H_INCLUDED
#define	__BASIS_H_INCLUDED	1	

extern char *basis_comment_string;


int	basis_set_shift_coords(int basistype,double *rbshx,double *rbshy);
int	basis_write_comment(FILE *fw,double bshx,double bshy);
int	basis_read_comment(char *buff,double *rbshx,double *rbshy);

#endif
               
