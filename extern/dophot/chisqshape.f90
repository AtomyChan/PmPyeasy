! Function chisqshape(a,err,starparrms,npstar,sig,lverb,&
!     flags,nff,enuffvpsf,a5,a6,a7,npsffit,ava)
! This function measures the agreement between the shape parameters of an object
! with those of a "typical" star to check if the object is "big"
! To better understand this function check section 4.4 in the DoPHOT paper 
! (Schechter, Mateo, and Saha 1993 PASP, 105, 1345).

function chisqshape(a,err,starparrms,npstar,sig,lverb,flags,nff,&
     enuffvpsf,a5,a6,a7,npsffit,ava)

  real :: a(npstar),err(npstar),b(npstar),starparrms(npstar) 
  real :: sig(3),tot(3),chi(3)
  real :: a5(npsffit),a6(npsffit),a7(npsffit),ava(npstar)
  character(*) :: flags(nff)
  logical :: enuffvpsf

  if ((a(5).ge.2*sqrt(err(5))).and.(a(6).ge.2*sqrt(err(6))).and.&
       (a(7).ge.2*sqrt(err(7)))) then          
     call stpar567 (a(3),a(4),flags,nff,enuffvpsf,a5,a6,a7,npsffit,ava,npstar,b)

     ! Check if the area of the associated footprint is increased or not
     areaa=ELLIPAREA(a(5),a(6),a(7))
     areab=ELLIPAREA(b(5),b(6),b(7))
     if (areaa.lt.areab) then
        chisqshape=0
        return
     end if

     ! Calculate chisq
     tot(1) = amax1(starparrms(5),(sig(1)*b(5))**2)
     if (b(6).ne.0.) then
        tot(2) = amax1(starparrms(6),sig(2)**2/b(6))
     else
        tot(2) = starparrms(6)
     end if
     tot(3) = amax1(starparrms(7),(sig(3)*b(7))**2)
     chi(1) = (a(5)-b(5))**2/(tot(1)+err(5))
     chi(2) = (a(6)-b(6))**2/(tot(2)+err(6))
     chi(3) = (a(7)-b(7))**2/(tot(3)+err(7))
     chisqshape = chi(1)+chi(2)+chi(3)
     if(lverb.gt.20) then
        write(13,*) ' (galaxy...) object at:', a(2), a(3)  
        write(13,*) ' chisq-x, chisq-xy, chisq-y, chisq-tot: ', chi    
     end if
  end if
    
end function chisqshape
