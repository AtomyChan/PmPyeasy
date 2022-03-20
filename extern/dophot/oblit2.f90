! Subroutine oblit(starpar,npstar,topsat,fitrad,lverb,rmagic,nfast,nslow,&
!     > BIG,NOISE)
! This subroutine is based on subroutine oblit.f90
! but it only obliterates a small inner region of a saturating star
! not useful to calculate its parameters

subroutine oblit2(starpar,npstar,topsat,fitrad,lverb,rmagic,nfast,nslow,&
     big,noise)

  real :: big(nfast,nslow),noise(nfast,nslow)
  real :: starpar(npstar)
  
  ixmin = max0(int(starpar(3)-fitrad),1)
  ixmax = min0(int(starpar(3)+fitrad)+1,nfast)
  iymin = max0(int(starpar(4)-fitrad),1)
  iymax = min0(int(starpar(4)+fitrad)+1,nslow)
  do j = iymin,iymax
     do i = ixmin,ixmax
        dist = sqrt(((starpar(3)-i)**2)+((starpar(4)-j)**2)) 
!       Fraction of the pixel inside the obliteration circle
        fractn = amax1(0.0,amin1(1.0,fitrad-dist+0.5)) 
        if ((fractn.gt.0).and.(big(i,j).gt.topsat)) then
           noise(i,j)=rmagic
        end if
     end do
  end do
  
end subroutine oblit2
