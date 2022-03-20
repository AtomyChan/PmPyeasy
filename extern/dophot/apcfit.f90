! Subroutine apcfit(x,y,mag,apin,apout,nap,files,nff,niter,lverb,&
!    npapmag,npapxy, > apparmag,apparxy,abig)
! Calculates the coefficients of the aperture correction functions
! as function of mag and xy, in addition to calculating the outer aperture
! correction value apbig 

subroutine apcfit(x,y,mag,apin,apout,nap,files,nff,niter,lverb,&
     npapmag,npapxy,apparmag,apparxy,abig)

  real :: x(nap),y(nap),mag(nap),apin(nap),apout(nap)
  real :: xy(2,nap),rmag(2,nap)
  real :: xyref(2,1),rmagref(2,1)
  real :: apparmag(npapmag),dumarrayf(npapmag)
  real :: apparxy(npapxy),dumarray(npapxy)
  real :: apcxy(nap),apcmag(nap),ee(nap)
  real :: resid(nap),apcfitt(nap),apcbig(nap)
  integer :: indexx(nap)
  character(*) :: files(nff)
  external surface
  external surface1

! Initialize the variables
  apparmag=0.
  apparxy=0.
  xyref=0.
  do i=1,nap
     xy(1,i) = x(i)
     xy(2,i) = y(i)
     rmag(1,i) = mag(i)
     rmag(2,i) = 0.0 !rmag(2,k)=0 only used cause neccessary for chisq
     ee(i) = 10.0**((mag(i)-10.0)/3.0)/5.0
     apcmag(i) = apin(i)
  end do

  do it=1,3
!    Gives the coeffs for the aperture correction as a function of magnitude
     call fixcoeff(surface1,rmag,apcmag,ee,nap,lverb,niter,npapmag,apparmag)
!    New values for apcxy
     do i=1,nap
        call surface1(rmag(1,i),apparmag,npapmag,dum,dum,dumarrayf,fit)
        apcxy(i) = apin(i) - fit
     end do
!    Gives the coeffs for the aperture correction as a function of xy
     call fixcoeff(surface,xy,apcxy,ee,nap,lverb,niter,npapxy,apparxy)
!    New values for apcmag
     do i=1,nap
        call surface(xy(1,i),apparxy,npapxy,dum,dum,dumarray,fit)
        apcmag(i) = apin(i) - fit
     end do
  end do

! Calculates outer aperture correction

  do k=1,nap
     call surface1(rmag(1,k),apparmag,npapmag,dum,dum,dumarrayf,fitmag)
     call surface(xy(1,k),apparxy,npapxy,dum,dum,dumarray,fitxy)
     apcfitt(k) = fitmag+fitxy
     resid(k) = apin(k) - apcfitt(k)
     apcbig(k) = apout(k) - apcfitt(k)
  end do
  call getapcbig(apcbig,ee,nap,abig,abigsig)

! Writes report file ordered by magnitude
  call quick(mag,nap,indexx)
  call write_report(files(12),nap,indexx,ee,x,y,mag,resid,apcfitt,apcbig,&
     apparmag,npapmag,apparxy,npapxy,abig,abigsig)

end subroutine apcfit

!-----------------------------------------------------------------------

! Subroutine fixcoeff(onestar,vals,apcvals,ee,nap,lverb,nit,npap,
!  > coeffs)
! Gives the coeffs for the aperture correction as a function of magnitude 

subroutine fixcoeff(onestar,vals,apcvals,ee,nap,lverb,nit,npap,coeffs)
     
  real :: vals(2,nap),apcvals(nap),ee(nap)
  real :: resid(nap)
  real :: coeffs(npap),err(npap)
  real :: wacc(npap), wlim(npap)
  external onestar
  
! Initialize the variables
  coeffs=0.
  err = 0.
  wacc = npap*0.001
  wlim = 0.

! Calculate the coefficients of the apcor function 
  call firstit(nap,apcvals,npap,coeffs)
  call chisq(onestar,dum,dum,vals,apcvals,ee,nap,coeffs,npap,npap,wacc,wlim,&
       nit,rmagic,lverb,err,chi)

! Sigma clip to improve fit
  call find_sigma(onestar,vals,apcvals,ee,nap,coeffs,npap,resid,sigma)
  clipfact=3.0
  call clip_stars(nap,resid,sigma,clipfact,ee)
  call chisq(onestar,dum,dum,vals,apcvals,ee,nap,coeffs,npap,npap,wacc,wlim,&
       nit,rmagic,lverb,err,chi)

! Repeat the whole operation
  call find_sigma(onestar,vals,apcvals,ee,nap,coeffs,npap,resid,sigma)
  clipfact=2.0
  call clip_stars(nap,resid,sigma,clipfact,ee)
  call chisq(onestar,dum,dum,vals,apcvals,ee,nap,coeffs,npap,npap,wacc,wlim,&
       nit,rmagic,lverb,err,chi)
  
end subroutine fixcoeff

!-----------------------------------------------------------------------

! Subroutine firstit(nap,apc,npap, > a)
! Calculates the coefficients of the funtion on the first iteration
! In this first iteration the best guess is the mean for a1 
! and 0 for the other coefficients
 
subroutine firstit(nap,apc,npap,a)

  real :: apc(nap),a(npap)

  asum=0.
  do i=1,nap
     asum = asum + apc(i)
  end do
  a(1) = asum/float(nap)
  do i=2,npap
     a(i) = 0.0
  end do
  
end subroutine firstit

!-----------------------------------------------------------------------

! Subroutine find_sigma(vals,apcvals,ee,nap,coeffs,npap, > resid,sigmaval)
! Work out residuals from the fit for each star
! and calculates the deviation

subroutine find_sigma(onestar,vals,apcvals,ee,nap,coeffs,npap,resid,sigmaval)

  real :: vals(2,nap),apcvals(nap),ee(nap),resid(nap)
  real :: coeffs(npap),dumarray(npap)
  external onestar
  
  sum = 0.
  nt = 0
  do k=1,nap
     call onestar(vals,coeffs,npap,dum,dum,dumarray,fit)
     resid(k) = apcvals(k) - fit
     if(ee(k).lt.10000.0) then
        sum = sum + resid(k)**2
        nt = nt + 1
     end if
  end do
  sigmaval = sqrt(sum/float(nt))
     
end subroutine find_sigma

!-----------------------------------------------------------------------

subroutine clip_stars(nap,resid,sigma,clipfact,ee)

  real :: resid(nap),ee(nap)

  do k=1,nap
     if(abs(resid(k)).gt.clipfact*sigma) then
        ee(k) = 10000.0
     end if
  end do

end subroutine clip_stars

!-----------------------------------------------------------------------

! Subroutine getapcbig(apcbig,ee,nap, > apcb,apcbsig)
! Calculates outer aperture correction

subroutine getapcbig(apcbig,ee,nap,apcb,apcbsig)

  real :: apcbig(nap),ee(nap),aphold(nap)

  ib = 1
  do i=1,nap
     if(ee(i).lt.10000.0) then
        aphold(ib) = apcbig(i)
        ib = ib + 1
     end if
  end do
  nb = ib - 1 
  apcb = 0.0
  apcbsig = 0.0
  if(nb.ge.1) then
     call median(aphold,nb,amed)
     call mad(aphold,nb,amed,amad)
     call biweightloc(aphold,nb,amed,amad,apcb)
     if(nb.gt.1) then
        call biweightscale(aphold,nb,amed,amad,apcbsig)
     end if
  end if

end subroutine getapcbig

!-----------------------------------------------------------------------

subroutine write_report(file,nap,indexx,ee,x,y,mag,resid,apcfitt,apcbig,&
     apparmag,npapmag,apparxy,npapxy,abig,abigsig)

  character(*) :: file
  real :: apparxy(npapxy),apparmag(npapmag)
  integer indexx(nap)
  real ee(nap),x(nap),y(nap),mag(nap),resid(nap),apcfitt(nap),apcbig(nap)
      
  nfiles = len_trim(file)
  call opena(12,file(1:nfiles)//'2',2,ierr)
  
  write(12,'(a,2f9.4)') 'abig,asigbig = ',abig,abigsig

  write(12,*) 'Aperture correction coefficients as a function of mag'
  do i=1,npapmag
     write(12,'(2x,e16.8)') apparmag(i)
  end do

  write(12,*) 'Aperture correction coefficients as a function of x and y'
  do i=1,npapxy
     write(12,'(2x,e16.8)') apparxy(i)
  end do
  
  do i=1,nap
     ii = indexx(i)
     if(ee(ii).lt.10000.0) then
        write(12,'(1x,2f8.2,4f7.3,a)') x(ii),y(ii),mag(i),resid(ii),&
             apcfitt(ii),apcbig(ii),' = x,y,mag,resid,apcxy+apcmag,apcbig'
     end if
  end do
  close(unit=12)
  
end subroutine write_report

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
