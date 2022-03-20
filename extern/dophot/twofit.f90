! Subroutine twofit(starpar,a,npstar,twostar,par1,par2,&
!     flags,nff,enuffvpsf,a5,a6,a7,npsffit,ava,&
!     rmagic,lverb,fitrad,xy,z,ze,maxfil,npt,nfit2,acc,alim,nit,&
!     > B1,B2,CHI2)
! This subroutine tests if the lights from the star fits a model of a blend
! of two close stars

subroutine twofit(starpar,a,npstar,twostar,par1,par2,&
     flags,nff,enuffvpsf,a5,a6,a7,npsffit,ava,&
     rmagic,lverb,fitrad,xy,z,ze,maxfil,npt,nfit2,acc,alim,nit,&
     b1,b2,chi2)

  real :: starpar(npstar),a(npstar)
  real :: acc(npstar),alim(npstar),bacc(npstar+3),blim(npstar+3)
  real :: a5(npsffit),a6(npsffit),a7(npsffit),ava(npstar)
  real :: xy(2,maxfil),z(maxfil),ze(maxfil)
  real :: b(npstar+3),berr(npstar+3),b1(npstar),b2(npstar)
  character(*) :: flags(nff)
  logical :: enuffvpsf,badfit
  external twostar

!...sense is that if x is long, angle is small.
!...sense is that if x&y are positively correlated, angle is positive.

  badfit = .false.

  b=0.
  call stpar567(a(3),a(4),flags,nff,enuffvpsf,a5,a6,a7,npsffit,ava,npstar,b)
  call ellipse(a(5),a(6),a(7),garea,dum,dum,angle)
  sarea = ELLIPAREA(b(5),b(6),b(7))
  dx2 = a(5)-b(5)
  dy2 = a(7)-b(7)
  dx = sqrt(amax1(dx2,0.))
  dy = sqrt(amax1(dy2,0.))
  badfit = amax1(dx,dy).eq.0
  dx = sign(dx,cos(angle))              
  dy = sign(dy,sin(angle))

  dx74 = starpar(3) - a(3)              
  dy74 = starpar(4) - a(4)              

  dot = dx74*dx + dy74*dy                       
  if (dot.gt.0) then                  
     fac1 = 0.3333333
     fac2 = 0.6666666
  else
     fac1 = 0.6666666
     fac2 = 0.3333333
  end if

  b(8) = b(5)
  b(9) = b(6)
  b(10) = b(7)
  b(1) = a(1)
  b(2) = log(a(2)*garea/sarea*fac2)              
  b(3) = a(3) - dx*fac1
  b(4) = a(4) - dy*fac1
  b(5) = log(a(2)*garea/sarea*fac1)              
  b(6) = a(3) + dx*fac2
  b(7) = a(4) + dy*fac2
  do i = 1,4
     bacc(i+3) = acc(i)
     bacc(i) = acc(i)
     blim(i+3) = alim(i)
     blim(i) = alim(i)
  end do
  bacc(2) = -0.01
  bacc(5) = -0.01
  blim(2) = -10.0
  blim(5) = -10.0     

  dxmax = amax1(abs(b(3)-a(3)),abs(b(6)-a(3)))
  dymax = amax1(abs(b(4)-a(4)),abs(b(7)-a(4)))
  badfit = badfit.or.(dxmax.gt.fitrad)
  badfit = badfit.or.(dymax.gt.fitrad)
  if (.not.badfit) then 
     call chisq(twostar,par1,par2,xy,z,ze,npt,b,npstar+3,&
          nfit2,bacc,blim,2*nit,rmagic,lverb,berr,chi2)
  end if
  dxmax = amax1(abs(b(3)-a(3)),abs(b(6)-a(3)))
  dymax = amax1(abs(b(4)-a(4)),abs(b(7)-a(4)))
  badfit = badfit.or.(dxmax.gt.fitrad)
  badfit = badfit.or.(dymax.gt.fitrad)

  b1(1) = b(1)
  b1(2) = exp(b(2))
  b1(3) = b(3)
  b1(4) = b(4)
  b1(5) = b(8)
  b1(6) = b(9)
  b1(7) = b(10)
  
  b2(1) = b(1)
  b2(2) = exp(b(5))
  b2(3) = b(6)
  b2(4) = b(7)
  b2(5) = b(8)
  b2(6) = b(9)
  b2(7) = b(10)

  if (badfit) chi2 = rmagic

end subroutine twofit

