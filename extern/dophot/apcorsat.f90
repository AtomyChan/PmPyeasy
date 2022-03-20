! Subroutine apcorr2(big,noise,nfast,nslow,starpar,itype,npstar,nsmax,&
!     appar,npap,onestar,par1,par2,nstot,flags,files,nff,niter,&
!     mapparmag,mapparxy,napsurffit,apfaintmag,&
!     napcmin,norder,nbadleft,nbadright,nbadtop,nbadbot,iadd,isub,&
!     fsub,xpnd,fac,rmagic,lverb,eperdn,rnoise,&
!     napertures,aprad,apskymin,apskymax,aperrmax,enuffap)
! This subroutine determines the aperture corrections and curves of growth.

subroutine apcorsat(big,noise,nfast,nslow,starpar,itype,npstar,nsmax,&
     appar,npap,onestar,par1,par2,nstot,flags,nff,&
     nbadleft,nbadright,nbadtop,nbadbot,iadd,isub,&
     fsub,xpnd,fac,rmagic,lverb,eperdn,rnoise,&
     aprad,apskymin,apskymax,aperrmax,niter,mapparmag)
  
  integer :: itype(nsmax)
  real :: starpar(npstar,nsmax),appar(npap,nsmax),dummy(npstar)
  real :: big(nfast,nslow),noise(nfast,nslow),xy(2),rrmag(2)
  real :: rmag(2,nsmax),apcmag(nsmax),ee(nsmax)
  real,allocatable :: apparmagn(:),dummagn(:)
!  real :: apparmagn(2),dummagn(2)
  character(*) :: flags(nff)
  external onestar
  external surface1



! Based on apcorr2.f90
! Loop through the stars.
! Deal only with type 9 with formal errors better than some predetermined limit

  iap = 0

  do i=1,nstot
     
     if (flags(13).eq.'NO') then
        if(.not.((itype(i).eq.9).or.(itype(i).eq.19))) cycle
        if((appar(4,i).gt.aperrmax).or.(appar(4,i).le.0.)) cycle
     end if

     call addstar(iadd,starpar(1,i),npstar,onestar,par1,par2,&
          fsub,xpnd,fac,rmagic,nfast,nslow,big,noise)
     
     xcenter = starpar(3,i)
     ycenter = starpar(4,i)
     ixlow = max0(int(xcenter-apskymax),nbadleft)
     ixhigh = min0(int(xcenter+apskymax)+1,nfast-nbadright)
     iylow = max0(int(ycenter-apskymax),nbadbot)
     iyhigh = min0(int(ycenter+apskymax)+1,nslow-nbadtop)
     aparea=0.
     apval=0.
     starlum1=0.
     starlum2=0.
     skyval=0.
     skywt=0.
     do jy=iylow,iyhigh
        do jx=ixlow,ixhigh
           x = float(jx)
           y = float(jy)
           xy(1) =  x
           xy(2) =  y
!          Distance from star center.
           dist = sqrt(((xcenter-x)**2)+((ycenter-y)**2))  
!          Fraction of the pixel inside the aperture circle
           fractnstar = amax1(0.0,amin1(1.0,aprad-dist+0.5))
           if (fractnstar.gt.0.) then
              aparea = aparea+fractnstar
              if((noise(jx,jy).le.0).or.(noise(jx,jy).ge.rmagic)) then 
                 call onestar(xy,starpar(1,i),npstar,par1,par2,dummy,fitval)
                 apval=apval+fitval*fractnstar
              else
                 apval=apval+(big(jx,jy)*fractnstar)
              end if
           end if
           fractnout = amax1(0.0,amin1(1.0,apskymax-dist+0.5))
           fractnin = amax1(0.0,amin1(1.0,dist-apskymin+0.5))
           fractnsky = amin1(fractnin,fractnout)
           skyval = skyval+big(jx,jy)*fractnsky/noise(jx,jy)
           skywt = skywt+fractnsky/noise(jx,jy)
        end do
     end do

     call addstar(isub,starpar(1,i),npstar,onestar,par1,par2,&
          fsub,xpnd,fac,rmagic,nfast,nslow,big,noise)

!    Based on apcfit.f90
     starlum = starpar(2,i)*ELLIPAREA(starpar(5,i),starpar(6,i),starpar(7,i))
     appar(1,i) = apval-skyval*aparea/skywt
     appar(2,i) = skyval/skywt
     ratlum = starlum/appar(1,i)
!print*,apval,skyval,aparea,skywt,skyval/skywt
     if (ratlum.le.0.) cycle
     iap=iap+1
     apcmag(iap) =  (2.5*alog10(ratlum)) - appar (3,i)
apcmag(iap) =  (2.5*alog10(ratlum))
!apcmag(iap) =  2.5*alog10(ratlum1) - appar (3,i)
!appar(3,i) =  (2.5*alog10(ratlum))
!print*,2.5*alog10(ratlum),apcmag(iap)
     rmag(1,iap) = SHAPEMAG(starpar(1,i),npstar)
     rmag(2,iap) = 0.0 !rmag(2,k)=0 only used cause neccessary for chisq
     ee(iap) = 10.0**((rmag(1,iap)-10.0)/3.0)/5.0
!     ee(iap) = 1/appar(4,i)**2
  end do

  napparmag = 0
  do i=1,mapparmag+1
     napparmag = napparmag + 1
  end do
  allocate (apparmagn(napparmag),dummagn(napparmag))
  apparmagn=0.
  dummagn=0.
! Gives the coeffs for the aperture correction as a function of magnitude
  call fixcoeff(surface1,rmag,apcmag,ee,iap,lverb,niter,napparmag,apparmagn)
! Based on apapply.f90
  do i=1,nstot

     if (flags(13).eq.'NO') then
        if(.not.((itype(i).eq.9).or.(itype(i).eq.19))) cycle
     end if

     rrmag(1) = SHAPEMAG(starpar(1,i), npstar)
     rrmag(2) = 0. !rmag(2,k)=0 only used cause neccessary for surface1
     call surface1(rrmag,apparmagn,napparmag,dum,dum,dummagn,amag)
     appar(3,i) = amag+appar(3,i)
appar(3,i) = amag

  end do

end subroutine apcorsat
