! Subroutine addstar(iadd,starpar,npstar,onestar,par1,par2,fsub,xpnd,fac,rmagic,
!     nfast,nslow, > big,noise)
! Subroutine adds (or substract) star model in the image array 
! and substracts (or adds) noise to the noise array.

subroutine addstar(iadd,starpar,npstar,onestar,par1,par2,fsub,xpnd,fac,rmagic,&
     nfast,nslow,big,noise)

  real,parameter :: sig2fw=2.35482
  real :: big(nfast,nslow),noise(nfast,nslow)
  real :: starpar(npstar),b(npstar),dummy(npstar),xy(2)
  integer :: jrect(4)
  external onestar


  sky = starpar(1)
  b = starpar
  if(xpnd.le.0.0) print *,'divide by zero 906.'
  b(5) = starpar(5)*xpnd**2
  b(6) = starpar(6)/xpnd**2                       
  b(7) = starpar(7)*xpnd**2

  call ellipse(starpar(5),starpar(6),starpar(7),dum,axmaj,axmin,dum)
  radi=fsub*axmaj/sig2fw
  ixmin = max0(nint(starpar(3)-radi),1)
  ixmax = min0(nint(starpar(3)+radi),nfast)
  iymin = max0(nint(starpar(4)-radi),1)
  iymax = min0(nint(starpar(4)+radi),nslow)
  do j = iymin,iymax
     xy(2) = float(j)
     do i = ixmin,ixmax
        if (noise(i,j).lt.rmagic) then
           xy(1) = float(i)
           dist = sqrt(((starpar(3)-xy(1))**2)+((starpar(4)-xy(2))**2))  
!          Fraction of the pixel inside the circle
           fractn = amax1(0.0,amin1(1.0,radi-dist+0.5))
           if (fractn.gt.0.0) then
              call onestar(xy,starpar,npstar,par1,par2,dummy,starval)
              big(i,j) = big(i,j)+iadd*(starval-sky)
              call onestar(xy,b,npstar,par1,par2,dummy,starval)
              rnoisehold = noise(i,j)-iadd*(fac*(starval-sky))**2
              if(rnoisehold.le.0.) then
                 big(i,j) = big(i,j)-iadd*(starval-sky)
              else
                 noise(i,j) = rnoisehold
              end if
           end if
        end if
     end do
  end do

end subroutine addstar
