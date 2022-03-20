/*****************************************************************************/
/* pbfft.h								     */
/*****************************************************************************/

#ifndef	__PBFFT_H_INCLUDED
#define	__PBFFT_H_INCLUDED	1

typedef struct _complex_struct
{	double re;
 	double im;
} complex;

/* pbfft_conv():
   Calculates the discrete Fourier transfomational of the complex array 'dt'
   (which has 'n' points) using the prime-based fast Fourier transformation 
   algorithm. The DFT of 'dt' is stored in 'dt', thus, if the original array
   is also wanted to be used after the call, it should be preserved before
   the call. If 'type' is false (zero), the function calculates the 
   DFT, if 'type' is true (non-zero), it calculates the inverse DFT.
   The function uses internal static and dynamic temporary arraies, and
   subsequential calls of the functions with different values of 'n' or 'type'
   can be slower than if 'n' and 'type' is not varied between the calls.     */
int	pbfft_conv(complex *dt,int n,int type);

#endif
                                          
