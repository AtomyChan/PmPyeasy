! Functions toobright(starpar,npstar,big,noise,nfast,nslow,top,fitrad,cmax,icrit,&
!     rmagic,lverb)
! This function finds out if star is too bright and saturate

function toobright(starpar,npstar,big,noise,nfast,nslow,top,fitrad,cmax,icrit,&
     rmagic,lverb)

  real :: starpar(npstar)
  real :: big(nfast,nslow), noise(nfast,nslow)
  logical :: toobright
   

  toobright = .false.

  if (starpar(1)+starpar(2).lt.top/2) return

  nsat = 0
  ixmin = max0(nint(starpar(3)-fitrad),1)
  ixmax = min0(nint(starpar(3)+fitrad),nfast)
  iymin = max0(nint(starpar(4)-fitrad),1)
  iymax = min0(nint(starpar(4)+fitrad),nslow)
  do j = iymin,iymax
     do i = ixmin,ixmax
        dist = sqrt(((starpar(3)-i)**2)+((starpar(4)-j)**2))  
!       Fraction of the pixel inside the circle
        fractn = amax1(0.0,amin1(1.0,fitrad-dist+0.5))
!       noise.ge.rmagic when big.gt.top or big.lt.bot (the case of the centers 
!       of badly saturated stars in some images)
        if ((fractn.gt.0.0).and.(noise(i,j).ge.rmagic)) nsat = nsat+1
     end do
  end do

  toobright = (starpar(1)+starpar(2).gt.cmax).or.(nsat.ge.icrit)
  if (toobright.and.(lverb.gt.20)) &
       write(13,*) nsat, ' saturated pixels in object at:', ix,iy 
      
end function toobright
