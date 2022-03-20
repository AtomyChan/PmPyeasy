! Subroutine apcorr2(big,noise,nfast,nslow,starpar,itype,npstar,nsmax,&
!     appar,npap,onestar,par1,par2,nstot,flags,files,nff,niter,&
!     mapparmag,mapparxy,napsurffit,apfaintmag,&
!     napcmin,norder,nbadleft,nbadright,nbadtop,nbadbot,iadd,isub,&
!     fsub,xpnd,fac,rmagic,lverb,eperdn,rnoise,&
!     napertures,aprad,apskymin,apskymax,aperrmax,enuffap)
! This subroutine determines the aperture corrections and curves of growth.

subroutine apcorr2(big,noise,nfast,nslow,starpar,itype,npstar,nsmax,&
     appar,npap,onestar,par1,par2,nstot,flags,files,nff,niter,&
     mapparmag,mapparxy,napsurffit,apfaintmag,&
     napcmin,norder,nbadleft,nbadright,nbadtop,nbadbot,iadd,isub,&
     fsub,xpnd,fac,rmagic,lverb,eperdn,rnoise,&
     napertures,aprad,apskymin,apskymax,aperrmax,enuffap)

  real,parameter :: pi=3.14159
  integer :: itype(nsmax)
  real :: starpar(npstar,nsmax),appar(npap,nsmax)
  real :: big(nfast,nslow),noise(nfast,nslow)
  real :: apcorr(napertures),apsum(napertures),aparea(napertures)
  real :: apx(nstot),apy(nstot),apamag(nstot)
  real :: apinner(nstot),apouter(nstot)
  real,allocatable :: apparxy(:),apparmag(:)
  logical :: skip,enuffap
  character(2) :: fmtvar
  character(120) :: fmt
  character(*) :: flags(nff),files(nff)
  external onestar

  napparmag = 0
  do i=1,mapparmag+1
     napparmag = napparmag + 1
  end do
  napparxy = MGETMFIT(mapparxy)

  call opena(12,files(12),2,ierr)

! Loop through the stars.
! Deal only with type 1 with formal errors better than some predetermined limit
! and not to close to the borders of the image
  iii = 0
  iap = 1

  do i=1,nstot

     if (flags(13).eq.'NO') then
        if(.not.((itype(i).eq.1).or.(itype(i).eq.11))) cycle
        if(appar(4,i).gt.aperrmax) cycle
     end if

     call addstar(iadd,starpar(1,i),npstar,onestar,par1,par2,&
          fsub,xpnd,fac,rmagic,nfast,nslow,big,noise)

     call aperpar(big,noise,nfast,nslow,starpar(1,i),npstar,rmagic,&
          aprad,apskymin,apskymax,nbadleft,nbadright,nbadtop,nbadbot,&
          napertures,npap,appar(1,i),apcorr,skip)

!    In case of doing just aperture photometry
     if (flags(13).eq.'YES') then
        appar(3,i) = apcorr(napsurffit)
        if (appar(3,i).lt.99.) then
           fitstar = starpar(2,i)*ELLIPAREA(starpar(5,i),starpar(6,i),starpar(7,i))
           apstar = fitstar/10**(0.4*appar(3,i))
           varstar = apstar*eperdn
           app = (float(napsurffit)/float(napertures))*aprad
           varsky = (pi*app**2)*(appar(2,i)*eperdn+rnoise**2)
           errelec = sqrt(varstar+varsky)
           errdn = errelec/eperdn
           appar(4,i) = 1.086*errdn/apstar   !dm=2.5d(logf)=2.5/ln(10)/f*df
        end if
     end if

     iii = iii + 1

     call addstar(isub,starpar(1,i),npstar,onestar,par1,par2,&
          fsub,xpnd,fac,rmagic,nfast,nslow,big,noise)

     if (skip) cycle

     amag = SHAPEMAG(starpar(1,i),npstar)

     write(fmtvar,'(i2)') napertures+1
     fmt = "(2i6,2f8.2,f7.1,"//fmtvar//"f7.3)"
     write(12,fmt) iii,itype(i),starpar(3,i),starpar(4,i),appar(2,i),amag,&
          (apcorr(k),k=1,napertures)
     
     if ((amag.le.apfaintmag).and.&
          (apcorr(napsurffit).lt.99).and.(apcorr(napertures).lt.99)) then
        apinner(iap) = apcorr(napsurffit)
        apouter(iap) = apcorr(napertures)
        apamag(iap) = amag
        apx(iap) = starpar(3,i)
        apy(iap) = starpar(4,i)
        iap = iap + 1
     end if

  end do
  close(unit=12) 

  nap = iap - 1
  if (lverb.gt.10) write(13,*) 'nap,napcmin = ',nap,napcmin

  allocate (apparmag(napparmag))
  allocate (apparxy(napparxy))
 
  enuffap = nap.ge.napcmin

  if(enuffap.and.(flags(13).eq.'NO')) then
     call apcfit(apx,apy,apamag,apinner,apouter,nap,files,nff,niter,lverb,&
          napparmag,napparxy,apparmag,apparxy,abig)
     call apapply(starpar,itype,npstar,nsmax,nstot,apparmag,napparmag,&
          apparxy,napparxy,abig,npap,appar)
  end if

end subroutine apcorr2
