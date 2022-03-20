! Subroutine vpsf (nsperf,nstot,starpar,itype,shadow,shaderr,ava,npstar,nsmax,&
!     files,nff,npsffit,npsfmin,rmagic,lverb, > ENUFFVPSF,A5,A6,A7)
! This subroutine calculates the variable psf coeffs.


subroutine vpsf (nsperf,nstot,starpar,itype,shadow,shaderr,ava,npstar,nsmax,&
     files,nff,npsffit,npsfmin,rmagic,lverb,enuffvpsf,a5,a6,a7)

  real :: starpar(npstar,nsmax),shadow(npstar,nsmax),shaderr(npstar,nsmax)
  integer :: itype(nsmax)
  real ::  xy(2,nsmax)
  real :: sp5(nsmax),sp6(nsmax),sp7(nsmax)
  real :: sperr5(nsmax),sperr6(nsmax),sperr7(nsmax)
  real :: err5(npsffit),err6(npsffit),err7(npsffit)
  real :: accvpsf(npsffit),alimvpsf(npsffit)
  real :: a5(npsffit),a6(npsffit),a7(npsffit)
  real :: ava(npstar)
  character(*) :: files(nff)
  character(120) :: fmt
  logical :: enuffvpsf
  external surface

! First determine if there are enough stars to justify a variable psf.
  if(nsperf.lt.npsfmin) then
     enuffvpsf = .false.
     return
  end if
  enuffvpsf = .true.

  nsperf = 0
  do i=1,nstot
     if (itype(i).eq.1) then
        nsperf = nsperf+1
        xy(1,nsperf) = starpar(3,i)
        xy(2,nsperf) = starpar(4,i)
        sp5(nsperf) = shadow(5,i)
        sp6(nsperf) = shadow(6,i)
        sp7(nsperf) = shadow(7,i)
        sperr5(nsperf) = shaderr(5,i)
        sperr6(nsperf) = shaderr(6,i)
        sperr7(nsperf) = shaderr(7,i)
     end if
  end do

! Initialize parameters for subroutine chisq 
  nnit = 50
  accvpsf = 0.10
  alimvpsf = 0.0
  if(a5(1).eq.0.0) then
     a5 = 0.
     a5(1) = ava(5)
  end if
  if(a6(1).eq.0.0) then
     a6 = 0.
     a6(1) = ava(6)
  end if
  if(a7(1).eq.0.0) then
     a7 = 0.
     a7(1) = ava(7)
  end if
  
  call chisq(surface,dum,dum,xy,sp5,sperr5,nsperf,a5,npsffit,&
       npsffit,accvpsf,alimvpsf,nnit,rmagic,lverb,err5,chi5)
  call chisq(surface,dum,dum,xy,sp6,sperr6,nsperf,a6,npsffit,&
       npsffit,accvpsf,alimvpsf,nnit,rmagic,lverb,err6,chi6)
  call chisq(surface,dum,dum,xy,sp7,sperr7,nsperf,a7,npsffit,&
       npsffit,accvpsf,alimvpsf,nnit,rmagic,lverb,err7,chi7)

  if(lverb.gt.10) write(13,*) 'VARIABLE PSF: CHI5 = ',chi5,&
       '   A5 = ',(a5(k),k=1,npsffit)
  if(lverb.gt.10) write(13,*) 'VARIABLE PSF: CHI6 = ',chi6,&
       '   A6 = ',(a6(k),k=1,npsffit)
  if(lverb.gt.10)  write(13,*) 'VARIABLE PSF: CHI7 = ',chi7,&
       '   A7 = ',(a7(k),k=1,npsffit)

  if(files(11)(1:1).ne.' ') then
     call opena(11,files(11),2,ierr)
     fmt="(a,i2,a,3(e12.4,e11.4))"
     do i=1,npsffit
        write(11,fmt) ' Coefficient(',i,') = ',a5(i),err5(i),a6(i),err6(i),&
             a7(i),err7(i)
     end do
     close(unit=11)     
  end if

end subroutine vpsf
