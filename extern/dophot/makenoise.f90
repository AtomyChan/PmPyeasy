!Subroutine makenoise(big,nfast,nslow,rnoise,eperdn,top,bot,&
!     nbadtop,nbadbot,nbadleft,nbadright,rmagic, > noise)

subroutine makenoise(big,nfast,nslow,rnoise,eperdn,top,bot,&
     nbadtop,nbadbot,nbadleft,nbadright,rmagic,noise)

  real ::  big(nfast,nslow),noise(nfast,nslow)
  logical :: badline, badpix

  temp1 = (rnoise/eperdn)**2
  temp2 = 1/eperdn
  do i = 1, nslow
     badline = ((i.le.nbadbot).or.(i.gt.nslow-nbadtop))
     do j = 1, nfast
        if (badline) then
           noise(j,i) = rmagic
        else
           badpix = ((j.le.nbadleft).or.(j.gt.nfast-nbadright))
           if ((badpix).or.(big(j,i).ge.top).or.(big(j,i).le.bot)) then
              noise(j,i) = rmagic
           else
              if (big(j,i).gt.0) then
                 noise(j,i) = big(j,i)*temp2 + temp1
              else
                 noise(j,i) = temp1
              end if
           end if
        end if
     end do
  end do

end subroutine makenoise
