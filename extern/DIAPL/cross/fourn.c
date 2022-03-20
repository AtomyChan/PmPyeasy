
/* This is Press's fourn() tidied up and converted to 0-indexed arrays */

#include <math.h>

#define PI2 6.28318530717959
#define SWAP(a,b) tempr = (a);(a) = (b);(b) = tempr

void fourn(float data[], int nn[], int ndim, int isign)
{
        int     i1, i2, i3, i2rev, i3rev, ip1, ip2, ifp1, ifp2;
        int     ibit, idim, k1, k2, n, nprev, ntot;
        float   tempi, tempr;
        double  theta, wi, wpi, wpr, wr, wtemp;

  ntot = 1;
  for (idim=0; idim < ndim; idim++)    ntot *= nn[idim];

/* loop over the dimensions */
  nprev = 1;
  for (idim=ndim-1; idim>=0; idim--)
  {
    n = nn[idim];
    ip1 = 2*nprev;
    ip2 = ip1*n;

/* Bit reversal */
    i2rev = 0;
    for (i2=0; i2<ip2; i2+=ip1)
    {
/* current axis */
      if (i2 < i2rev)
      {
        for (i1=i2; i1<i2+ip1-1; i1+=2)
        {
/* previous axes */
          for (i3=i1; i3<2*ntot; i3+=ip2)
          {
/* future dimensions */
            i3rev = i2rev + i3 - i2;
            SWAP(data[i3], data[i3rev]);
            SWAP(data[i3+1], data[i3rev+1]);
          }
        }
      }

      ibit = ip2/2;
      while (ibit >= ip1 && i2rev >= ibit)
      {
        i2rev -= ibit;
        ibit /= 2;
      }
      i2rev += ibit;
    }

/* And the Danielson-Lanczos recursion */
    ifp1 = ip1;
    while (ifp1 < ip2)
    {
      ifp2 = 2*ifp1;
      theta = isign*PI2/(ifp2/ip1);
      wtemp = sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi = sin(theta);
      wr = 1.0;
      wi = 0.0;

      for (i2=0; i2<ifp1; i2+=ip1)
      {
        for (i1=i2; i1<i2+ip1-1; i1+=2)
        {
          for (i3=i1; i3<2*ntot; i3+=ifp2)
          {
            k1 = i3;
            k2 = k1 + ifp1;
            tempr = wr*data[k2] - wi*data[k2 + 1];
            tempi = wr*data[k2 + 1] + wi*data[k2];
            data[k2] = data[k1] - tempr;
            data[k2 + 1] = data[k1 + 1] - tempi;
            data[k1] += tempr;
            data[k1 + 1] += tempi;
          }
        }

        wr = (wtemp = wr)*wpr - wi*wpi + wr;
        wi = wi*wpr + wtemp*wpi + wi;
      }
      ifp1 = ifp2;
    }
    nprev *= n;
  }

  return;
}
/* (C) Copr. 1986-92 Numerical Recipes Software 1!57. */
