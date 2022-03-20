! Subroutine apapply(starpar,itype,npstar,nsmax,nstot,apparmag,npapmag,&
!     apparxy,npapxy,abig,npap,appar)
! This subroutine calculate the aperture correction for every star

subroutine apapply(starpar,itype,npstar,nsmax,nstot,apparmag,npapmag,&
     apparxy,npapxy,abig,npap,appar)

  real :: starpar(npstar,nsmax),appar(npap,nsmax)
  integer :: itype(nsmax)
  real :: xy(2),rmag(2)
  real :: apparmag(npapmag),dummag(npapmag)
  real :: apparxy(npapxy), dumxy(npapxy)
  logical :: dumlog

  do i=1,nstot
     if((itype(i).eq.4).or.(itype(i).eq.6).or.(itype(i).eq.8)) cycle
     if((itype(i).eq.14).or.(itype(i).eq.16).or.(itype(i).eq.18)) cycle
     xy(1) = starpar(3,i)
     xy(2) = starpar(4,i)
     rmag(1) = shapemag(starpar(1,i), npstar)
     rmag(2) = 0. !rmag(2,k)=0 only used cause neccessary for surface1
     call surface1(rmag,apparmag,npapmag,dum,dum,dummag,amag)
     call surface(xy,apparxy,npapxy,dum,dum,dumxy,axy)
     appar(3,i) = amag+axy+abig
     if((itype(i).eq.9).or.(itype(i).eq.19)) appar(3,i) = axy
  end do

end subroutine apapply

