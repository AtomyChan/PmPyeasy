! Subroutine oblit(starpar,npstar,obltype,fobl,cosm,lverb,rmagic,nfast,nslow,&
!     > BIG,NOISE)
! This subroutine oblitarates a region (circular or rectangular).
! Stars will not be found in the obliterated region.

subroutine oblit(starpar,npstar,obltype,fobl,cosm,lverb,rmagic,nfast,nslow,&
     big,noise)

  real,parameter :: sig2fw=2.35482
  real :: big(nfast,nslow),noise(nfast,nslow)
  real :: starpar(npstar)
  logical :: cosm
  character(1) :: obltype    ! obltype = 'r' (rectangular) or 'c' (circle)
  
  x = starpar(3)
  y = starpar(4)

  if ((cosm).or.(obltype.eq.'r')) then
     roblx = starpar(5)
     robly = starpar(7)
  else
     call ellipse(starpar(5),starpar(6),starpar(7),dum,axmaj,axmin,dum)
     roblx = fobl*axmaj/sig2fw
     robly = roblx
  end if

  ixhi = min0(nint(x+roblx),nfast)
  ixlo = max0(nint(x-roblx),1)       
  iyhi = min0(nint(y+robly),nslow)   
  iylo = max0(nint(y-robly),1)
  if(lverb.gt.30) then
     if (obltype.eq.'c') &
          write(13,*) ' Obliterating a circle in the following region :'
     if (obltype.eq.'r') &
          write(13,*) ' Obliterating following region :'
     write(13,*) 'IXLO, IXHI, IYLO, IYHI = ',ixlo, ixhi, iylo, iyhi
  end if
  do jy = iylo,iyhi
     do jx = ixlo,ixhi
        dist = sqrt(((x-jx)**2)+((y-jy)**2))  
!       Fraction of the pixel inside the obliteration circle
        fractn = amax1(0.0,amin1(1.0,roblx-dist+0.5))
        if (((fractn.gt.0).and.(obltype.eq.'c')).or.(obltype.eq.'r')) then
           big(jx,jy) = -rmagic
           noise(jx,jy) = rmagic
        end if
     end do
  end do
  if(lverb.gt.40) then
     write(13,*) 'OBLITERATION: X & Y = ',x,y
     if (obltype.eq.'c') &
          write(13,*) 'OBLITERATION: R_OBL = ',robl
     if (obltype.eq.'r') &
          write(13,*) 'OBLITERATION: WX & WY = ',roblx,robly
  end if
  
end subroutine oblit
