! Subroutine MATHSUB
! This file contains subroutines that are not I/O related, but rather
! involve arithmetic or character operations that my be somewhat 
! machine-specific-- it may be necessary or desirable to make some
! changes to optimize the code to run on your computer.
!
! Current contents:
!
! function VALUE     computes at a given point the value of a point-spread
!                      function that consists of the sum of a two-dimensional
!                      integral under a bivariate Gaussian function and a 
!                      correction obtained by interpolation within a lookup
!                      table with half-pixel resolution.
! function RINTER    does the two-dimensional interpolation within the lookup 
!                      table.
! function CUBINT    given known values at four equally-spaced points,computes
!                      a value at an aribitary point between the second and 
!                      third points with known values, by cubic interpolation.
! function ERFF      computes the integral under a one-dimensional Gaussian 
!                      by Simpson's 1/3 rule.
!! subroutine INVERS  inverts a square matrix.
!! subroutine VMUL    multiplies a square matrix (on the left) by a column 
!                      vector (on the right) yielding a new column vector.
! subrotine QUICK    does a quicksort on a vector of data.
!! subroutine SHOW    produces what is effectively an eleven-level gray scale 
!!                      plot of a rectangular data array on the terminal.
!! function SWITCH    peels the filename-extension off of a filename, and 
!!                      replaces it with another.
! function ICNVRT    takes a character string and returns an integer 
!                      representing the first two characters in the string.
! function ICNVRT3    takes a character string and returns an integer 
!                      representing the first two characters in the string.
!! subroutine ASKFIL  asks (neatly) for a filename, offering a default.
!! subroutine INQUIRE asks (neatly) for input.
!
! subroutine MMM
!
!-----------------------------------------------------------------------
!
! Function VALUE
! This function returns the value of a point-spread function at a
! given point.  The value of the point-spread function is the sum of a
! two-dimensional integral under a bivariate Gaussian function, and 
! a value obtained by interpolation in a look-up table.
!
!     XX,YY    are the real coordinates of the desired point relative 
!                  to the centroid of the point-spread function.
!   GAUSS(1)   is the peak height of the best-fitting Gaussian profile.
!   GAUSS(2,3) are x and y offsets from the centroid of the point-spread
!                  function to the center of the best-fitting Gaussian.
!   GAUSS(4,5) are the x and y sigmas of the best-fitting Gaussian.
!      PSF     is an array containing the look-up table.
!   DVDX,DVDY  are the first derivatives of the composite point-spread
!                  function with respect to x and y.
!
! Please note that, although the arguments XX,YY of the function VALUE are
! relative to the centroid of the PSF, the function RINTER which VALUE calls
! requires coordinates relative to the corner of the array (see below).

function value(xx, yy, gauss, psf, npsf, dvdx, dvdy)

  integer,parameter :: maxpsf = 69
  real :: psf(maxpsf, maxpsf), gauss(5)

  half = float(npsf + 1) / 2   !This is not standard f90?
  value = 0.
  dvdx = 0.
  dvdy = 0.
  x = (2. * xx) + half
  y = (2. * yy) + half
  if ((((x .lt. 2.) .or. (x .gt. (npsf - 1.))) .or. (y .lt. 2.))&
       .or. (y .gt. (npsf - 1.))) return 
  erfx = erff(xx,gauss(2),gauss(4),dedxc,dummy)
  erfy = erff(yy,gauss(3),gauss(5),dedyc,dummy)

!   Since the lookup table has a grid size of one-half pixel in each
!   coordinate, the derivatives must be multiplied by two to yield
!   the derivatives in units of ADU/pixel in the big frame.

  value = ((gauss(1) * erfx) * erfy) + rinter(psf,maxpsf,maxpsf,x,y,&
       dfdx,dfdy,ist)
  dvdx = (2. * dfdx) - ((gauss(1) * dedxc) * erfy)
  dvdy = (2. * dfdy) - ((gauss(1) * dedyc) * erfx)
  
  return 
  
end function value

!-----------------------------------------------------------------------

! function RINTER
! This function interpolates in a two-dimensional look-up table, F,
! to estimate a value at the position X,Y.  It also provides estimates
! of the first derivative of F with respect to x and y at the position
! X,Y.  F is dimensioned (1:NX,1:NY). X=1.0, Y=1.0 refers to the center
! of the first pixel; X=FLOAT(NX), Y=FLOAT(NY) refers to the center of
! the last pixel.  Values for RINTER, DFDX, and DFDY are returned only
! when 2.0 <= X <= FLOAT(NX-1), 2.0 <= Y <= FLOAT(NY-1); in this case
! ISTAT is set to zero.  Otherwise, RINTER=DFDX=DFDY=0.0 and ISTAT=1.
!
! The method used is cubic interpolation: first, at integral values of
! y, the values of F at integral values of x are interpolated to x=X to
! give a value for the function and its first derivative with respect 
! to x.  Then these are interpolated to the point x,y=X,Y to give a 
! functional value and its first derivative with respect to y.  
! Finally, the first derivative of F with respect to x at x,y = X,Y is 
! obtained by interpolating among the values of dF/dx obtained at the 
! integral values of y.  It so happens that the values obtained for the
! function and its derivatives would be the same if we had interpolated
! in y first, rather than in x.

function rinter(f, nx, ny, x, y, dfdx, dfdy, istat)
  
  real :: f(nx, ny), g(-1:2), dgdx(-1:2)
  
  if ((((x .lt. 2.) .or. (x .gt. float(nx - 1))) .or. (y .lt. 2.))&
       .or. (y .gt. float(ny - 1))) then
     istat = 1
     rinter = 0.
     dfdx = 0.
     dfdy = 0.
     return 
  end if
  
  ix = int(x)
  iy = int(y)
  dx = x - ix
  dy = y - iy
  do i = -1, 2
     j = iy + i
     g(i) = cubint(f(ix,j),dx,dgdx(i))
  end do
  rinter = cubint(g(0),dy,dfdy)
  dfdx = cubint(dgdx(0),dy,dummy)
  return 
  
end function rinter

!-----------------------------------------------------------------------

! Function CUBINT
! Given values for a function F(x) at positions x=-1,0,1,2, this 
! function estimates a value for the function at a position X, using a
! cubic interpolating polynomial.  The specific polynomial used has the
! properties CUBINT(0)=F(0), CUBINT(1)=F(1), and the first derivative 
! of CUBINT with respect to x is continuous at x=0 and x=1.
!
! Input arguments:
!     F is the 0-element of the series F(-1), F(0), F(+1), F(+2).
!     X is the REAL*4 distance between the position of the 0-element
!          and the desired point:  0.0 <= X < 1.0.
!
! Output arguments:
!     CUBINT is the interpolated value of the function F at position X.
!     DFDX is the estimate of dF/dx at position X.

function cubint(f, x, dfdx)

  real :: f(-1:2)

  c1 = 0.5 * (f(1) - f(-1))
  c2 = ((2. * f(1)) + f(-1)) - (0.5 * ((5. * f(0)) + f(2)))
  c3 = 0.5 * (((3. * (f(0) - f(1))) + f(2)) - f(-1))
  cubint = (x * ((x * ((x * c3) + c2)) + c1)) + f(0)
  dfdx = (x * (((x * c3) * 3.) + (2. * c2))) + c1
  return 
  
end function cubint

!-----------------------------------------------------------------------

! Function ERFF
! Numerically integrate a Gaussian function 
!
!         { F = EXP[-0.5*(DELTAX/BETA)**2] },
!
! from XIN-0.5 to XIN+0.5 using Simpson's 1/3 rule.  Also provide first
! derivative of the integral with respect to Xo and BETA.
!
! The number of intervals required to end up with an error less than
! ALPHA is greater than 
!
!   fourth root of [(fourth derivative of F w.r.t. x) / (180 * ALPHA)].
!
! Here I am using ALPHA = 0.00005.   N must be an even number.

function erff(xin, xo, beta, dfdxo, dfdbet)

  betasq = beta ** 2
  x = ((xin - xo) / beta) ** 2
  f = exp(- (0.5 * x))
  n = max(2,ifix((3.247*((f*abs((x*(x-6.))+3.))**0.25))/beta)+1)
  if (mod(n,2) .ne. 0) n = n + 1

!   Start with the lower end point (weight = 1).

  dx = 1. / float(n)
  deltax = (xin - xo) - 0.5
  dxsq = deltax ** 2
  f = exp(- ((0.5 * dxsq) / betasq))
  erff = f
  dfdxo = f * deltax

!   Now include the end points of each subinterval except the last one. 
!   If it is an odd-numbered subinterval, weight = 4.  If even, weight = 2.

  dfdbet = f * dxsq
  do i = 1, n - 1
     deltax = deltax + dx
     dxsq = deltax ** 2
     f = exp(- ((0.5 * dxsq) / betasq))
     fwt = (f * 2.) * (1. + mod(i,2))
     erff = erff + fwt
     dfdxo = dfdxo + (deltax * fwt)     
!   Now add the upper end point (weight = 1) and multiply by DX/3.
     dfdbet = dfdbet + (dxsq * fwt)
  end do
  deltax = deltax + dx
  dxsq = deltax ** 2
  f = exp(- ((0.5 * dxsq) / betasq))
  dx = dx / 3.
  erff = dx * (erff + f)
  dfdxo = (dx * (dfdxo + (deltax * f))) / betasq
  dfdbet = (dx * (dfdbet + (f * dxsq))) / (betasq * beta)
  return 
  
end function erff

!-----------------------------------------------------------------------
!
!! Subroutine INVERS
!! Although it seems counter-intuitive, the tests that I have run
!! so far suggest that the 180 x 180 matrices that NSTAR needs can
!! be inverted with sufficient accuracy if the elements are REAL*4
!! rather than REAL*8.
!
!subroutine invers(a, max, n, iflag)
!
!  real(8) :: a(max, max)
!  
!  iflag = 0
!  i = 1
!300 if (a(i,i) .eq. 0.0) goto 999
!  a(i,i) = 1.0 / a(i,i)
!  j = 1
!301 if (j .eq. i) goto 304
!  a(j,i) = - (a(j,i) * a(i,i))
!  k = 1
!302 if (k .eq. i) goto 303
!  a(j,k) = a(j,k) + (a(j,i) * a(i,k))
!303 if (k .eq. n) goto 304
!  k = k + 1
!  goto 302
!304 if (j .eq. n) goto 305
!  j = j + 1
!  goto 301
!305 if (k .eq. i) goto 306
!  a(i,k) = a(i,k) * a(i,i)
!306 if (k .eq. 1) goto 307
!  k = k - 1
!  goto 305
!307 if (i .eq. n) return 
!  i = i + 1
!  goto 300
!999 iflag = 1
!  return 
!
!end subroutine invers
!
!-----------------------------------------------------------------------
!
!subroutine vmul(a, max, n, v, x)
!
!  real(8) :: sum, a(max, max), v(max), x(max)
!
!  i = 1
!200 sum = 0.0d0
!  j = 1
!201 sum = sum + (dble(a(i,j)) * dble(v(j)))
!  if (j .eq. n) goto 203
!  j = j + 1
!  goto 201
!203 x(i) = sngl(sum)
!  if (i .eq. n) return 
!  i = i + 1
!  goto 200
!  
!end subroutine vmul
!
!-----------------------------------------------------------------------

! Subroutine QUICK
! A quick-sorting algorithm suggested by the discussion on pages 114-119
! of THE ART OF COMPUTER PROGRAMMING, Vol. 3, SORTING AND SEARCHING, by
! D.E. Knuth, which was referenced in Don Wells' subroutine QUIK.  This
! is my own attempt at encoding a quicksort-- PBS.
!
! The array DATUM contains randomly ordered data. The integer vector
! index1(i) will return to the calling program where the i-th element
! in the array AFTER sorting had been BEFORE the array was sorted.
! The limiting stack length of 14 will limit this quicksort subroutine
! to vectors of maximum length of order 32,768 (= 2**15).

subroutine quick(datum, n, index1)

  integer,parameter :: maxstak = 14
  real :: datum(n)
  integer(4) :: index1(n), stklo(maxstak), stkhi(maxstak), hi

!   Initialize INDEX1.
  do i = 1, n
     index1(i) = i
  end do
!    Initialize the pointers.
  nstak = 0
  limlo = 1

  limhi = n
  100 dkey = datum(limlo)

! Compare all elements in the sub-vector between LIMLO and LIMHI with
! the current key datum.

  ikey = index1(limlo)
  lo = limlo
  hi = limhi
101 continue
  if (lo .eq. hi) goto 200
  if (datum(hi) .le. dkey) goto 109

! The pointer HI is to be left pointing at a datum SMALLER than the
! key, which is intended to be overwritten.

  hi = hi - 1

  goto 101
109 datum(lo) = datum(hi)
  index1(lo) = index1(hi)
  lo = lo + 1
110 continue
  if (lo .eq. hi) goto 200
  if (datum(lo) .ge. dkey) goto 119
  lo = lo + 1

  goto 110
119 datum(hi) = datum(lo)
  index1(hi) = index1(lo)

! The pointer LO is to be left pointing at a datum LARGER than the
! key, which is intended to be overwritten.

  hi = hi - 1

  goto 101

! LO and HI are equal, and point at a value which is intended to
! be overwritten.  Since all values below this point are less than
! the key and all values above this point are greater than the key,
! this is where we stick the key back into the vector.

200 continue
  datum(lo) = dkey

! At this point in the subroutine, all data between LIMLO and LO-1, 
! inclusive, are less than DATUM(LO), and all data between LO+1 and 
! LIMHI are larger than DATUM(LO).

! If both subarrays contain no more than one element, then take the most
! recent interval from the stack (if the stack is empty, we're done).
! If the larger of the two subarrays contains more than one element, and
! if the shorter subarray contains one or no elements, then forget the 
! shorter one and reduce the other subarray.  If the shorter subarray
! contains two or more elements, then place the larger subarray on the
! stack and process the subarray.

  index1(lo) = ikey

! Case 1:  the lower subarray is longer.  If it contains one or no 
! elements then take the most recent interval from the stack and go 
! back and operate on it.

  if ((limhi - lo) .gt. (lo - limlo)) goto 300

! If the upper (shorter) subinterval contains one or no elements, then
! process the lower (longer) one, but if the upper subinterval contains
! more than one element, then place the lower (longer) subinterval on
! the stack and process the upper one.

  if ((lo - limlo) .le. 1) goto 400

! Case 1a:  the upper (shorter) subinterval contains no or one elements,
! so we go back and operate on the lower (longer) subinterval.

  if ((limhi - lo) .ge. 2) goto 250
  limhi = lo - 1

  goto 100

! Case 1b:  the upper (shorter) subinterval contains at least two 
! elements, so we place the lower (longer) subinterval on the stack and
! then go back and operate on the upper subinterval.
 
250 continue
  nstak = nstak + 1
  stklo(nstak) = limlo
  stkhi(nstak) = lo - 1
!     DO 3666 I=1,NSTAK
!3666 TYPE *,'STACK: ',STKLO(I),STKHI(I)
  limlo = lo + 1

  goto 100

! Case 2:  the upper subarray is longer.  If it contains one or no 
! elements then take the most recent interval from the stack and 
! operate on it.

300 continue

! If the lower (shorter) subinterval contains one or no elements, then
! process the upper (longer) one, but if the lower subinterval contains
! more than one element, then place the upper (longer) subinterval on
! the stack and process the lower one.

  if ((limhi - lo) .le. 1) goto 400

! Case 2a:  the lower (shorter) subinterval contains no or one elements,
! so we go back and operate on the upper (longer) subinterval.

  if ((lo - limlo) .ge. 2) goto 350
  limlo = lo + 1

  goto 100

! Case 2b:  the lower (shorter) subinterval contains at least two 
! elements, so we place the upper (longer) subinterval on the stack and
! then go back and operate on the lower subinterval.
 
350 continue
  nstak = nstak + 1
  stklo(nstak) = lo + 1
  stkhi(nstak) = limhi
  limhi = lo - 1

  goto 100

! Take the most recent interval from the stack.  If the stack happens 
! to be empty, we are done.

400 continue
  if (nstak .le. 0) return 
  limlo = stklo(nstak)
  limhi = stkhi(nstak)
  nstak = nstak - 1

  goto 100

end subroutine quick

!-----------------------------------------------------------------------
!
!subroutine show(f, fmax, fzero, nx, ny, maxdim)
!
!  integer,parameter :: nchar = 11
!  real :: f(maxdim, maxdim)
!!      integer*2 char(nchar)
!  character(2) :: char1(nchar)
!  data char1 /'  ','- ','--','::','==','ll','II','%%','00','HH','##'/
!
!  low = (40 - nx) - 1
!  s = sqrt(max(float(nchar),fmax - fzero))   !float is an extension to f90?
!  write(unit=6, fmt=601) 
!601 format(<low+2>x,<nx>(2h--))
!  do 801 iy = 1, ny
!801  write(unit=6, fmt=602) (char1(min(nchar,&
!          ifix((nchar * sqrt(max(0.,f(ix,iy) - fzero))) / s) + 1)),ix = 1, nx)
!602  format(1x,<low>x,1h|,<nx>a2,1h|)
!     write(unit=6, fmt=601) 
!     return 
!
!end subroutine show

!-----------------------------------------------------------------------
!
!function switch(nafile, addend)
!
!  character(*) ::  nafile, addend,switch
!  character c
!!  character, allocatable :: switch(:)
!
!  nnafile=len_trim(nafile)
!  naddend=len_trim(addend)
!  do i = nnafile, 1, -1
!     c = nafile(i:i)
!     if (c .eq. '.') exit
!  end do
!!  allocate (switch(i+naddend))
!  switch = nafile(1:i) // addend
!  return 
!
!end function switch
!
!-----------------------------------------------------------------------

! This little function is supposed to take two ASCII characters and
! express them as an integer in the range 0-(32**2-1) without 
! distinguishing upper and lower case:
! AA = Aa = aA = aa = 0, AB = Ab = aB = ab = 1, BA = Ba = bA = ba = 32,etc.

function icnvrt(astring)

  character(2) :: astring

  icnvrt = (32 * mod(ichar(astring(1:1))-1,32)) + mod(ichar(astring(2:2))-1,32)

end function icnvrt

!-----------------------------------------------------------------------

function icnvrt3(astring)

  character(3) :: astring

  icnvrt3 = ((1024 * mod(ichar(astring(1:1)) - 1,32)) + (32 * &
       mod(ichar(astring(2:2)) - 1,32))) + mod(ichar(astring(3:3)) - 1,32)
  
end function icnvrt3

!-----------------------------------------------------------------------
!
!subroutine askfil(prompt, file)
!
! If input file name includes a version number, wipe it out.
!
!  character prompt*40, file*30, infile*30
!  do 780 i = 1, 30
!     if (file(i:i) .eq. ';') goto 781
!780  continue
!     goto 783
!781  do 782 j = i, 30
!782     file(j:j) = ' '
!783     continue
!     do 700 i = 1, 40
!        n = ichar(prompt(i:i))
!        if (n .eq. 58) goto 800
!700     continue
!800   i = i - 1
!      do 701 j = 30, 1, -1
!         n = ichar(file(j:j))
!         if ((n .ne. 0) .and. (n .ne. 32)) goto 801
!701      continue
!      j = max(1,49 - i)
!702   write(unit=6, fmt=600) prompt(1:i)
!  600 format(<j>x,a<i>,2h: ,$)
!      read(unit=5, fmt=500, end=999, err=702) file
!500   format(a32)
!      k=lenc(file)
!      if (k .le. 0) goto 702
!      return 
!801   k = max(1,(38 - i) - j)
!703   write(unit=6, fmt=601) prompt(1:i), file(1:j)
!601   format(<k>x,a<i>,10h (default ,a<j>,3h): ,$)
!      read(unit=5, fmt=500, end=998, err=703) k, infile
!      if (k .ne. 0) file = infile
!      return 
!998   file = ' '
!999   return 
!   
!end subroutine askfil
!
!-----------------------------------------------------------------------
!
!subroutine inquire(prompt)
!  character prompt*50
!  do 700 i = 1, 50
!     n = ichar(prompt(i:i))
!     if ((n .eq. 58) .or. (n .eq. 63)) goto 800
!700  continue
!800  k = max(1,50 - i)
!     write(unit=6, fmt=600) prompt(1:i)
!600  format(<k>x,a<i>,1h ,$)
!     return 
!
!end subroutine inquire
!
!=======================================================================
!
! This version of MMM (modified by PBS 1984.IV.10ff) assumes that
! the sky brightnesses in the one-dimensional array SKY are already
! sorted on entering this routine, and that pixels outside the "bad"
! limits have already been eliminated.
!
! This particular version of MMM also takes cognizance of the fact that,
! pixels falling below the LOWBAD threshold already having been 
! eliminated, the contaminated sky pixels values overwhelmingly display
! POSITIVE departures from the true value.
!
! If for some reason it is impossible to obtain the mode of the sky
! distribution, this will be flagged by setting SIGMA = -1.0.
!
! Arguments
!
!     SKY (INPUT) is a real vector containing actual sorted sky values.
!    NSKY (INPUT) is the number of defined elements in SKY.
! SKYMODE (OUTPUT) is the estimated mode of the sky values.
!   SIGMA (OUTPUT) is the computed standard deviation of the peak in
!         the sky histogram.
!    SKEW (OUTPUT) is the computed skewness of the peak in the sky
!         histogram.
!
!=======================================================================

subroutine mmm(sky, nsky, highbad, skymode, sigma, skew)

  real(8) :: sum, sumsq
  real :: sky(nsky)
  logical :: redo

!-----------------------------------------------------------------------

! SECTION 1

  data maxiter / 20 /
  data minsky / 20 /

!   An error condition has been detected.
  if (nsky .le. 0) then 
     sigma = -1.0
     skew = 0.0
     return
  endif

!   SKYMID is the median value for the whole ensemble of sky pixels.
!   Notice that if NSKY is odd, then (NSKY+1)/2 and (NSKY/2)+1 will be the
!   same number, whereas if NSKY is even, they will be different numbers.
!   This same trick will be used again later.
  skymid = 0.5 * (sky((nsky + 1) / 2) + sky((nsky / 2) + 1))

!   Initialize the variables for accumulating the mean and standard
!   deviation, and initialize the rejection limits.
  sum = 0.
  sumsq = 0.
  minimum = 0
  maximum = 0

!   For the first pass we will consider only pixels in a symmetric 
!   interval of brightness values about the median value.  This exploits
!   the assumption that all the bad pixels are already rejected from the
!   lower end of the brightness range.
  cut1 = min(skymid - sky(1),sky(nsky) - skymid,highbad - skymid)
  cut2 = skymid + cut1

  cut1 = skymid - cut1
  do i = 1, nsky
     if (sky(i) .lt. cut1) then
        minimum = i
        cycle
     end if
     if (sky(i) .gt. cut2) exit
     delta = sky(i) - skymid
     sum = sum + delta
     sumsq = sumsq + (delta ** 2)
     maximum = i
  end do
!   Henceforth in this subroutine, MINIMUM will point to the highest value
!   rejected at the lower end of the vector, and MAXIMUM will point to the
!   highest value accepted at the upper end of the vector.
!   MAXIMUM-MINIMUM is the number of pixels within the acceptance range.

!   Compute mean and sigma (from the first pass).

  skymed = 0.5 * (sky(((minimum + maximum) + 1) / 2) + sky(((minimum&
       + maximum) / 2) + 1))
  skymean = sum / dble(maximum - minimum)
  sigma = dsqrt((sumsq / dble(maximum - minimum)) - (skymean ** 2))

!   The middle sky value, SKYMID, was subtracted off up above and added 
!   back in down here to reduce the truncation error in the computation 
!   of SIGMA.
!   Note that this definition of SIGMA is incorrect by a factor of
!   SQRT [NSKY/(NSKY-1.)], but for all but pathological cases (where none
!   of this can be trusted anyway), it's close enough.

  skymean = skymean + skymid
  skymode = skymean

!  If the mean is less than the mode, that means the contamination is
!  slight, and the mean value is what we really want.  Note that this
!  introduces a slight bias toward underestimating the sky when
!  the scatter in the sky is caused by random fluctuations rather than
!  by contamination, but I think this bias is negligible compared to the
!  problem of contamination.

!-----------------------------------------------------------------------

! SECTION 2

!   Rejection and recomputation loop:

  if (skymed .lt. skymean) skymode = (3. * skymed) - (2. * skymean)
  niter = 0
  redo=.true.

  do while (redo)
     niter = niter + 1

!   Compute Chauvenet rejection criterion.

     if ((niter .gt. maxiter) .or. ((maximum - minimum) .lt. minsky)) then 
        sigma = -1.0
        skew = 0.0
        return
     endif
     r = alog10(float(maximum - minimum))
     
!   Compute rejection limits (symmetric about the current mode).

     r = max(2.,(((- (.1042 * r)) + 1.1695) * r) + .8895)
     cut = (r * sigma) + (0.5 * abs(skymean - skymode))
     cut = max(1.5,cut)
     cut1 = skymode - cut

!   Recompute mean and sigma by adding and/or subtracting sky values
!   at both ends of the interval of acceptable values.

!   At each end of the interval, ISTEP will show the direction we have to 
!   step through the vector to go from the old partition to the new one.
!   Pixels are added or subtracted depending upon whether the limit is 
!   moving toward or away from the mode.

     cut2 = skymode + cut

!   Is CUT1 above or below the minimum currently-accepted value?

     redo = .false.
     istep = sign(1.0001,cut1 - sky(minimum + 1))

!   If ISTEP = +1, JSTEP = 1.  If ISTEP = -1, JSTEP=0.  If ISTEP = +1, 
!   then we know that at least one pixel must be deleted at the low end.

     jstep = (istep + 1) / 2
     if (istep .gt. 0) then
        delta = sky(minimum + jstep) - skymid
        sum = sum - (istep * delta)
        sumsq = sumsq - (istep * (delta ** 2))
        minimum = minimum + istep
!   A change has occured
        redo = .true.
     endif
     do
!   Quit when SKY(MINIMUM) < CUT1 <= SKY(MINIMUM+1)
        if ((istep .lt. 0) .and. (minimum .le. 0)) exit
!   If ISTEP is positive, subtract out the sky value at MINIMUM+1; if 
!   ISTEP is negative, add in the sky value at MINIMUM.
        if ((sky(minimum) .le. cut1) .and. (sky(minimum + 1) .ge. cut1)) exit
        delta = sky(minimum + jstep) - skymid
        sum = sum - (istep * delta)
        sumsq = sumsq - (istep * (delta ** 2))
        minimum = minimum + istep
!   A change has occured
        redo = .true.
     end do

!   Is CUT2 above or below the current maximum?

     istep = sign(1.0001,cut2 - sky(maximum))

!   If ISTEP = +1, JSTEP = 1.  If ISTEP = -1, JSTEP=0.  If ISTEP = -1, then
!   we know that we must subtract at least one pixel from the high end.

     jstep = (istep + 1) / 2
     if (istep .lt. 0) then
        delta = sky(maximum + jstep) - skymid
        sum = sum + (istep * delta)
        sumsq = sumsq + (istep * (delta ** 2))
        maximum = maximum + istep
!   A change has occured
        redo = .true.
     end if

     do
!   Quit when SKY(MAXIMUM) <= CUT2 < SKY(MAXIMUM+1)
        if((istep .gt. 0) .and. (maximum .ge. nsky)) exit
!   If ISTEP is positive, add in the sky value at MAXIMUM+1; if ISTEP is 
!   negative, subtract off the sky value at MAXIMUM.
        if ((sky(maximum) .le. cut2) .and. (sky(maximum + 1) .ge. cut2)) exit 
        delta = sky(maximum + jstep) - skymid
        sum = sum + (istep * delta)
        sumsq = sumsq + (istep * (delta ** 2))
        maximum = maximum + istep
!   A change has occured
        redo = .true.

     end do

!   Compute mean and sigma (from this pass).

     skymean = sum / dble(maximum - minimum)
     sigma = dsqrt((sumsq / dble(maximum - minimum)) - (skymean ** 2))

!   Obtain the median.  To first approximation, the median would be the
!   value of the sky in the middle pixel in the sorted data (if the
!   total number is odd) or the mean of the two pixels straddling
!   the middle (if the total number of pixels is even).
!
!     SKYMED=0.5*(SKY((MINIMUM+MAXIMUM+1)/2)+SKY((MINIMUM+MAXIMUM)/2+1))
!
!   However, this is not good enough.  If you look at the estimator for
!   the mode, you will note that a tiny change in the list of sky pixels,
!   just sufficient to alter the median value of the sky brightness by
!   one unit, will change the estimator of the mode by three units.  We
!   really want something more robust than this.  As a first attempt
!   at a more robust median estimator, I propose to estimate the median
!   of the distribution by the mean of the central ten or eleven sky
!   values:

     skymean = skymean + skymid
     skymed = 0.0
     
     x = 0.0
     do  i=((minimum + maximum)-8)/2,((minimum + maximum)+11)/2
        skymed = skymed + sky(i)
        x = x + 1.
     end do
     skymed = skymed / x
     skymode = skymean

!   If the mean is less than the mode, that means the contamination is
!   slight, and the mean value is what we really want.  Note that this
!   introduces a slight bias toward underestimating the sky when
!   the scatter in the sky is caused by random fluctuations rather than
!   by contamination, but I think this bias is negligible compared to the
!   problem of contamination.
!
!   If the limits have not yet stopped moving, try again.

     if (skymed .lt. skymean) skymode = (3. * skymed) - (2. * skymean)
     
  end do
   
  skew = (skymean - skymode) / max(1.,sigma)
  nsky = maximum - minimum

  return 

end subroutine mmm
