! Subroutine aperpar(big,noise,nfast,nslow,starpar,npstar,rmagic,&
!   apmax,apskymin,apskymax,nbadleft,nbadright,nbadtop,nbadbot,napertures,npap,&
!   > APPAR,APCORR,SKIP)
! This subroutine calculates the aperture photometry, and also the aperture 
! correction comparing the fit and aperture magnitudes.
! The aperture information is kept in an array of 4 elements, 3 of which are
! calculated here
! apar(1) keeps the aperture flux from the object
! apar(2) keeps the aperture flux from the sky
! apar(3) is the aperture correction magnitude (fit minus aperture)

subroutine aperpar(big,noise,nfast,nslow,starpar,npstar,rmagic,&
    apmax,apskymin,apskymax,nbadleft,nbadright,nbadtop,nbadbot,napertures,npap,&
    appar,apcorr,skip)

  real :: big(nfast,nslow), noise(nfast,nslow)
  real :: starpar(npstar), appar(npap)
  real :: apcorr(napertures),apval(napertures),aparea(napertures)
  logical :: badpix,skip

  appar(1) = 0.
  appar(2) = 0.
  appar(3) = 0.
  badpix = .false.
  skip = .false.
  skywt = 0.0
  skyval = 0.0
  aparea = 0.0
  apval = 0.0
  apcorr = 99.999

  xcenter = starpar(3)
  ycenter = starpar(4)
  ixlow = int(xcenter-apskymax)
  if (ixlow.le.nbadleft) skip = .true.
  ixhigh = int(xcenter+apskymax)+1
  if(ixhigh.gt.(nfast-nbadright)) skip = .true.
  iylow = int(ycenter-apskymax)
  if(iylow.le.nbadbot) skip = .true.
  iyhigh = int(ycenter+apskymax)+1
  if(iyhigh.gt.(nslow-nbadtop)) skip = .true.
  if (skip) return

  do jy=iylow,iyhigh
     do jx=ixlow,ixhigh
        badpix = noise(jx,jy).ge.rmagic
        x = float(jx)
        y = float(jy)
!       Distance from star center.
        dist = sqrt(((xcenter-x)**2)+((ycenter-y)**2))  
        do k=1,napertures
           app = (float(k)/float(napertures))*apmax
!          Fraction of the pixel inside the aperture circle
           fractnstar = amax1(0.0,amin1(1.0,app-dist+0.5))
           if((fractnstar.gt.0.0).and.badpix) then
              skip = .true.
              return
           end if
           aparea(k) = aparea(k)+fractnstar
           apval(k) = apval(k)+(fractnstar*big(jx,jy))
        end do
        fractnout = amax1(0.0,amin1(1.0,apskymax-dist+0.5))
        fractnin = amax1(0.0,amin1(1.0,dist-apskymin+0.5))
        fractnsky = amin1(fractnin,fractnout)
        if((fractnsky.gt.0.0).and.badpix) cycle
        if(noise(jx,jy).le.0) cycle
        skyval = skyval+big(jx,jy)*fractnsky/noise(jx,jy)
        skywt = skywt+fractnsky/noise(jx,jy)
     end do
  end do

  starlum = starpar(2)*ELLIPAREA(starpar(5),starpar(6),starpar(7))
  do ii=1,napertures
     appar(1) = apval(ii)-skyval*aparea(ii)/skywt
     if(appar(1).le.0.0) then
        appar(1) = 0.0
        apcorr(ii) = 99.999
     else
        appar(2) = skyval/skywt
!       This is an independent way of calculating the difference between
!       aperture and fit mags. Large values for abs(appar(3)) for star 
!       objects will indicate problems.  This value equals the difference 
!       between the sky and fit mags.
        ratlum = starlum/appar(1)
        if (ratlum.gt.0) appar(3) = 2.5*alog10(ratlum)
        apcorr(ii) = 2.5*alog10(ratlum)
     end if
  end do

end subroutine aperpar
