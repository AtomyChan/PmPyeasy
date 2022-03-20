! Subroutine fillerup(ixin,iyin,big,noise,nfast,nslow,snlim,fitrad,rmagic,
!     maxfil,npstar > XY,Z,ZE,NPT,NPTTOT,A)
! This subroutine populates the arrays xy, z, and ze with npt elements.
! These arrays are later used by chisq.f90
! This subroutine also assigns integer values to the analytic stellar function 
! parameters 1 and 2 (see explanation in stpar.f90)


subroutine fillerup(x,y,big,noise,nfast,nslow,snlim,fitrad,rmagic,&
     maxfil,npstar,xy,z,ze,npt,areafill,areatot,a)

  real :: big(nfast,nslow),noise(nfast,nslow)
  real :: xy(2,maxfil),z(maxfil),ze(maxfil)
  real :: a(npstar)

  snlim2 = snlim**2
  areatot = 0.
  areafill = 0.
  npt = 0
  xy = 0.
  z = 0.
  ze = 0.
  sum0 = 0.
  sum1 = 0.

  ixmin = max0(int(x-fitrad),1)
  ixmax = min0(int(x+fitrad)+1,nfast)
  iymin = max0(int(y-fitrad),1)
  iymax = min0(int(y+fitrad)+1,nslow)
  do j = iymin,iymax
     do i = ixmin,ixmax
        dist = sqrt((x-i)**2+(y-j)**2)  
!       Fraction of the pixel inside the circle
        fractn = amax1(0.0,amin1(1.0,fitrad-dist+0.5))
        areatot = areatot+fractn
        if ((big(i,j)**2.lt.snlim2*noise(i,j)).or.(noise(i,j).ge.rmagic)&
             .or.(fractn.eq.0.0)) cycle
        areafill = areafill+fractn
        npt = npt+1
        xy(1,npt) = i
        xy(2,npt) = j
        z(npt) = big(i,j)
        ze(npt) = noise(i,j)
        if(noise(i,j).eq.0.) then
           temp = 0.0
        else
           temp = 1/noise(i,j)
        end if
        sum0 = sum0 + temp
        sum1 = sum1 + big(i,j)*temp
     end do
  end do

  a(1) = 0.
  a(2) = 0.
  if (sum0.ne.0) then
     a(1) = sum1/sum0
     a(2) = big(int(x),int(y))-a(1)
  end if
  a(3) = x
  a(4) = y

end subroutine fillerup
