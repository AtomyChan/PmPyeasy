! Makes certain statistical calculus
!
! Includes these subroutines:
! AVERAGE
! WAVERAGE (Weighted Average)
! MEDIAN
! MAD (Median Absolute Deviation)
! BIWEIGHTLOC
! BIWEIGHTSCALE
!
!-----------------------------------------------------------------------

subroutine average(array,narray,avgout,sigout,sabsout)

  real :: array(narray)
  real(8) :: sum,sum2,suma,avg,sig,sabs

  sum = 0.0
  sum2 = 0.0
  suma = 0.0

  do i=1,narray
     sum = sum + dble(array(i))
  end do
  an = dble(float(narray))
  avg = sum/an

  do i=1,narray
     sum2 = sum2 + (dble(array(i))-avg)**2
     suma = suma + dabs(dble(array(i))-avg)
  enddo

  if(narray.gt.1) then
     sig = dsqrt(sum2/(an-1.0))
     sabs = suma/(an-1.0)
  else 
     sig = -10.0
     sabs = -10.0
  end if

  avgout = sngl(avg)
  sigout = sngl(sig)
  sabsout = sngl(sabs)

  return
end subroutine average

!-----------------------------------------------------------------------

! Subroutine to calculate the weighted averaged of a clipped set of values
! If you don't want to clip the data, just make clip=0

subroutine waverage(array,sigarray,narray,clip,wavgout,wsigout1,wsigout2,inotclip)

  real :: array(narray),sigarray(narray)
  real(8) :: sum1,sum2,sum3,sum4,weight(narray)
  real(8) :: wavg,wsig1,wsig2,wavgold,wsigold1,wsigold2

  itest=0
  if (clip.eq.0) itest=1
  sum1 = 0d0
  sum2 = 0d0
  sum3 = 0d0
  sum4 = 0d0
  inotclip = narray
  do i=1,narray
     if (sigarray(i).eq.0) cycle
     weight(i) = 1.0/dble(sigarray(i))**2
     sum1 = sum1 + weight(i)*dble(array(i))
     sum2 = sum2 + weight(i)
     sum3 = sum3 + weight(i)**2
  end do
  wavg = sum1/sum2
  do i=1,narray
     sum4 = sum4 + weight(i)*((array(i)-wavg)**2)
  end do
  wsig1 = dsqrt(1.0/sum2)
  if (narray.eq.1) then
     wsig2=(array(1)-wavg)**2
  else
     wsig2 = dsqrt(sum2*sum4/(sum2**2-sum3))
  end if
  wavgout = sngl(wavg)
  wsigout1 = sngl(wsig1)
  wsigout2 = sngl(wsig2)
  do while (itest.eq.0)
     wavgold = wavg
     wsigold1 = wsig1
     wsigold2 = wsig2
     sum1 = 0d0
     sum2 = 0d0
     sum3 = 0d0
     sum4 = 0d0
     inotclip = 0
     do i=1,narray
        weight(i)=0d0 
        if (sigarray(i).eq.0) cycle
        if ((array(i).gt.wavg-clip*sigarray(i)).and.&
             (array(i).lt.wavg+clip*sigarray(i))) then
           weight(i) = 1.0/dble(sigarray(i))**2
           inotclip = inotclip+1
        end if
        sum1 = sum1 + weight(i)*dble(array(i))
        sum2 = sum2 + weight(i)
        sum3 = sum3 + weight(i)**2
      end do
     if (inotclip.eq.0) then
        print*,'Clipping eliminates all the elements'
        wavgout=0.0
        wsigout1=0.0
        wsigout2=0.0
        return
     end if
     wavg = sum1/sum2 
     do i=1,narray
        sum4 = sum4 + weight(i)*((array(i)-wavg)**2)
     end do
     wsig1 = dsqrt(1.0/sum2)
     if (inotclip.eq.1) then
        wsig2=(array(1)-wavg)**2
     else
        wsig2 = dsqrt(sum2*sum4/(sum2**2-sum3))
     end if
     wavgout = sngl(wavg)
     wsigout1 = sngl(wsig1)
     wsigout2 = sngl(wsig2)
     if ((wavg.eq.wavgold).and.(wsig1.eq.wsigold1).and.(wsig2.eq.wsigold2)) itest=1
  end do

end subroutine waverage

!-----------------------------------------------------------------------

! Subroutine to calculate the median of a set of values.

subroutine median(values,nvalues,amedian)

  real :: values(nvalues)
  integer :: index(nvalues)
  real, allocatable :: holdval(:)

  allocate (holdval(nvalues))
  do i=1,nvalues
     holdval(i) = values(i)
  enddo
  call quick(holdval,nvalues,index)

  if(mod(nvalues,2).eq.0) then
     nmid = nvalues/2
     amedian = 0.5*(holdval(nmid)+holdval(nmid+1))
  else
     nmid = nvalues/2 + 1
     amedian = holdval(nmid)
  end if

  return
  
end subroutine median

!-----------------------------------------------------------------------

! Subroutine to calculate the median absolute deviation (MAD).

subroutine mad(values,nvalues,amedian,amad)

  real :: values(nvalues)
  real,allocatable :: diff(:)

  allocate (diff(nvalues))
  do i=1,nvalues
     diff(i) = abs(values(i)-amedian)
  end do

  call median(diff,nvalues,amad)

  return
  
end subroutine mad

!-----------------------------------------------------------------------

subroutine biweightloc(values,nvalues,amedian,amad,bwloc)

  real,parameter :: tuner=6.0
  real :: values(nvalues),u(nvalues)

  sum1 = 0.0
  sum2 = 0.0

  do i=1,nvalues
     diff = values(i)-amedian
     u(i) =diff/(tuner*amad)
     udiff = (1.0 - u(i)**2)**2
     if(abs(u(i)).lt.1.0) then
        sum1 = sum1 + diff*udiff
        sum2 = sum2 + udiff
     end if
  end do
    
  bwloc = amedian + sum1/sum2      

  return
  
end subroutine biweightloc

!-----------------------------------------------------------------------

subroutine biweightscale(values,nvalues,amedian,amad,bwscale)

  real,parameter :: tuner=9.0
  real values(nvalues),u(nvalues)

  sum1 = 0.0
  sum2 = 0.0

  do i=1,nvalues
     diff = values(i)-amedian
     u(i) =diff/(tuner*amad)
     udiff = 1.0 - u(i)**2
     if(abs(u(i)).lt.1.0) then
        sum1 = sum1 + (diff**2)*(udiff**4)
        sum2 = sum2 + udiff*(1.0 - 5.0*u(i)**2)
     end if
  end do

  bwscale = sqrt(float(nvalues))*(sqrt(sum1)/abs(sum2))

  return
  
end subroutine biweightscale

!-----------------------------------------------------------------------
