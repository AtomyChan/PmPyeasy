! This subroutine finds 'volcano' stars (saturated stars 
! where the counts decrease in the central part) 
! and increase the counts and noise of the central region
! so other subroutines do not have problems with these stars

subroutine volcano(big,noise,nfast,nslow,rmagic,maxfil,fitrad,top,topsat)

  real :: big(nfast,nslow),noise(nfast,nslow)
  integer :: xsat(maxfil),ysat(maxfil)

  topmax=top*2

  do iy = 1,nslow
     do ix = 1,nfast
        if (big(ix,iy).lt.topsat) cycle
        xcen = float(ix)
        ycen = float(iy)
        do it=1,3
           ixmin = max0(int(xcen-fitrad),1)
           ixmax = min0(int(xcen+fitrad)+1,nfast)
           iymin = max0(int(ycen-fitrad),1)
           iymax = min0(int(ycen+fitrad)+1,nslow)
           bigmax= maxval(big(ixmin:ixmax,iymin:iymax))
           if (bigmax.ge.topmax) exit
           k=0
           do jy = iymin,iymax
              do jx = ixmin,ixmax
                 if (big(jx,jy).gt.topsat) then
                    k=k+1
                    xsat(k)=jx
                    ysat(k)=jy
                 end if
              end do
           end do
           idistkcurrent=-1
           do k2=1,k
              do k3=1,k
                 idistk=(xsat(k2)-xsat(k3))**2+(ysat(k2)-ysat(k3))**2
                 if (idistk.gt.idistkcurrent) then
                    idistkcurrent=idistk
                    iextreme1=k2
                    iextreme2=k3
                 end if
              end do
           end do
           xcen=xsat(iextreme1)+(xsat(iextreme2)-xsat(iextreme1))/2.
           ycen=ysat(iextreme1)+(ysat(iextreme2)-ysat(iextreme1))/2.
           distobl=sqrt(float(idistkcurrent))/2.
        end do
        if (bigmax.ge.topmax) cycle
        ixmin = max0(int(xcen-distobl),1)
        ixmax = min0(int(xcen+distobl)+1,nfast)
        iymin = max0(int(ycen-distobl),1)
        iymax = min0(int(ycen+distobl)+1,nslow)
        distcurrent1=0.
        do jy = iymin,iymax
           do jx = ixmin,ixmax
              dist=sqrt((xcen-jx)**2+(ycen-jy)**2)
              if ((noise(jx,jy).ge.rmagic).and.(dist.gt.distcurrent1)) &
                distcurrent1=dist
           end do
        end do
        distcurrent2=0.
        bigmax= maxval(big(ixmin:ixmax,iymin:iymax))
        do jy = iymin,iymax
           do jx = ixmin,ixmax
              dist = sqrt((xcen-jx)**2+(ycen-jy)**2) 
              if (big(jx,jy).eq.bigmax) distcurrent2 = dist
           end do
        end do
        distobl2=max(distcurrent1,distcurrent2)
        ixmin = max0(int(xcen-distobl2),1)
        ixmax = min0(int(xcen+distobl2)+1,nfast)
        iymin = max0(int(ycen-distobl2),1)
        iymax = min0(int(ycen+distobl2)+1,nslow)
        do jy = iymin,iymax
           do jx = ixmin,ixmax
              dist = sqrt((xcen-jx)**2+(ycen-jy)**2)
!             Fraction of the pixel inside the obliteration circle
              fractn = amax1(0.0,amin1(1.0,distobl2-dist+0.5))
              if (fractn.ge.1.) big(jx,jy)=topmax
              if (fractn.ge.1.) noise(jx,jy)=rmagic
           end do
        end do
        
     end do
  end do

end subroutine volcano
