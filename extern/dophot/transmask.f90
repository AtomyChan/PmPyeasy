! Function transmask(ix,iy,sky,starmask,maskrad,big,noise,nfast,nslow,&
!     crit,rmagic,lverb)
! This function tests the flux though a mask representing the current model 
! of a stellar PSF.

function transmask(x,y,sky,starmask,maskrad,big,noise,nfast,nslow,&
     crit,rmagic,lverb)

  real :: big(nfast,nslow),noise(nfast,nslow)
  real :: starmask(-maskrad:maskrad,-maskrad:maskrad)
  logical :: transmask

  transmask = .false.
  crit2 = crit**2
  sum1 = 0                        
  sum2 = 0
  do j = -maskrad, maskrad
     jj = int(y)+j
     if ((jj.lt.1).or.(jj.gt.nslow)) cycle
     do i = -maskrad, maskrad
        ii = int(x)+i
        if ((ii.lt.1).or.(ii.gt.nfast)) cycle
        dist = sqrt((x-ii)**2+(y-jj)**2)  
!       Fraction of the pixel inside the circle
        fractn = amax1(0.0,amin1(1.0,maskrad-dist+0.5))
        if ((fractn.le.0.0).or.(noise(ii,jj).ge.rmagic)) cycle
        if(noise(ii,jj).eq.0) then
           temp = 0.0
        else
           temp = starmask(i,j)/noise(ii,jj)
        end if
        sum1 = sum1+starmask(i,j)*temp          
        sum2 = sum2+temp*(big(ii,jj)-sky)
     end do
  end do
  if (sum2.gt.0) then
     if(sum1.eq.0) then
        rat = 0.0
     else
        rat = sum2**2/sum1
     end if
     if (rat.gt.crit2) then
        transmask = .true.
        if ((rat.lt.1.1*crit2).and.(lverb.gt.30)) &
             write(13,*) 'MARGINAL: (S/N)**2 through Mask = ', rat
     end if
  end if
    
end function transmask
