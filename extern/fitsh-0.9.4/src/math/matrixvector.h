/*****************************************************************************/
/* matrixvector.h							     */
/*****************************************************************************/

#ifndef	__MATRIXVECTOR_H_INCLUDED
#define	__MATRIXVECTOR_H_INCLUDED	1

/*****************************************************************************/

typedef double  matrix[3][3];
typedef double  vector[3];

int	matrix_mul(matrix a,matrix b,matrix c);
int	vector_mul(vector a,matrix b,vector c);
double	vector_matrix_vector_mul(vector a,matrix b,vector c);
double	vector_vector_mul(vector a,vector c);
int	matrix_add(matrix a,matrix b,matrix c);
int	matrix_diad(matrix a,vector b,vector c);
int	matrix_sub(matrix a,matrix b,matrix c);
int	matrix_copy(matrix a,matrix b);
int	matrix_real_mul(matrix a,matrix b,double d);
int	matrix_unity(matrix a);
int	matrix_zero(matrix a);
int	matrix_antisym(matrix m,double a,double b,double c);
int	matrix_is_diff(matrix a,matrix b);
int	matrix_exp(matrix e,matrix a);
int	matrix_rot(matrix o,double a,double b,double c);
double	vector_length(vector v);
int	vector_normalize(vector v);

/*****************************************************************************/

#endif

/*****************************************************************************/
           
