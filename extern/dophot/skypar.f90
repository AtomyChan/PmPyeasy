! Subroutine skypar(nstot,nfast,nslow,maxfil,nsmax,npstar,&
!     big,noise,starpar,shadow,itype,ava,flags,nff,rmagic,lverb,&
!     npskyp,npskyh,iskymedstep,&
!     nsskypmin,itskyp,accskyp,alimskyp,nsskyhmin,itskyh,accskyh,alimskyh,&
!     > SKYPLANEPAR,SKYHUBPAR,SKYMEDARRAY)
! This subroutine calculates the values of the parameters or arrays neccesary
! to evaluate the sky functions.

subroutine skypar(nstot,nfast,nslow,maxfil,nsmax,npstar,&
     big,noise,starpar,shadow,itype,ava,flags,nff,rmagic,lverb,&
     npskyp,npskyh,iskymedstep,&
     nsskypmin,itskyp,accskyp,alimskyp,nsskyhmin,itskyh,accskyh,alimskyh,&
     skyplanepar,skyhubpar,skymedarray)

  real :: skyplanepar(npskyp),skyhubpar(npskyh),skyhubparinit(npskyh)
  real :: skymedarray(nfast/iskymedstep,nslow/iskymedstep)
  real :: big(nfast,nslow),noise(nfast,nslow)
  real :: starpar(npstar,nsmax),shadow(npstar,nsmax),ava(npstar)
  real :: xy(2,maxfil),z(maxfil),ze(maxfil)
  real :: accskyp(npskyp),alimskyp(npskyp),accskyh(npskyh),alimskyh(npskyh)
  real :: dummyskyp(npskyp,npskyp),dummyskyh(npskyh,npskyh)
  real :: values((nfast*nslow)/(iskymedstep*iskymedstep))
  integer :: itype(nsmax)
  character(*) :: flags(nff)
  logical :: goodstar,conv
  external skyfun_plane     
  external skyfun_hub



  skyplanepar = 0.
  skyhubpar = 0.
  skymedarray = 0.



  if ((flags(2)(1:5).eq.'PLANE').or.(flags(2)(1:5).eq.'HUBBL')) then
 

     ngood = 0
     skymax = -rmagic
     do i = 1,nstot 
        goodstar = (itype(i).eq.1).or.((itype(i).eq.3).and.(shadow(1,i).ne.0))
        if (goodstar) then
           if (ngood.lt.maxfil) then
              if (starpar(1,i).gt.skymax) then
                 skymax = starpar(1,i)
                 xm = starpar(3,i)
                 ym = starpar(4,i)
              end if
              ngood = ngood+1
              xy(1,ngood) = starpar(3,i)
              xy(2,ngood) = starpar(4,i)
              z(ngood) = starpar(1,i)
              ze(ngood) = 1
           end if
        end if
     end do     
     if(lverb.gt.10) then
        write(13,*) '# of stars available for computing model '&
             //'sky (ngood) = ', ngood
     end if
     
     if (flags(2)(1:5).eq.'PLANE') then
        if (ngood.ge.nsskypmin) then
           if (skyplanepar(1).eq.0.) skyplanepar(1) = ava(1)
           call chisq(skyfun_plane,dum,dum,xy,z,ze,ngood,skyplanepar,npskyp,&
                npskyp,accskyp,alimskyp,itskyp,rmagic,lverb,dummyskyp,dum)
        else
           skyplanepar(1) = ava(1)
        end if
        if(lverb.gt.10) write(13,*) ' Skyfun_plane parameters: '
        if(lverb.gt.10) write(13,*) (skyplanepar(i),i=1,npskyp)
     else if(flags(2)(1:5).eq.'HUBBL') then
        if (ngood.ge.nsskyhmin) then
           if (skyhubpar(2) .eq. 0) then
              skyhubpar(1) = skyguess
              skyhubpar(2) = skymax - skyhubpar(1)
              skyhubpar(3) = xm
              skyhubpar(4) = ym
              skyhubpar(5) = nfast
              skyhubpar(6) = 0
              skyhubpar(7) = nslow
           end if
           do i=1,npskyh
              skyhubparinit(i) = skyhubpar(i)
           end do
           call chisq(skyfun_hub,dum,dum,xy,z,ze,ngood,skyhubpar,npskyh,&
                npskyh,accskyh,alimskyh,itskyh,rmagic,lverb,dummyskyh,chisqh)
           conv = chisqh.lt.rmagic 
           if(.not.conv) then
              do i=1,npskyh
                 skyhubpar(i) = skyhubparinit(i) 
              end do
           end if
        else
           skyhubpar(1) = skyguess
           skyhubpar(2) = 0
           skyhubpar(3) = nfast/2
           skyhubpar(4) = nslow/2
           skyhubpar(5) = nfast
           skyhubpar(6) = 0
           skyhubpar(7) = nslow
        end if
        if(lverb.gt.10) write(13,*) ' Skyfun_hub parameters:'
        if(lverb.gt.10) write(13,*) (skyhubpar(i),i=1,npskyh)
     end if


  else if (flags(2)(1:5).eq.'MEDIA') then


     iskymedhalf = iskymedstep/2 
     jm = 0
     do j=iskymedstep,nslow,iskymedstep
        jm = jm + 1
        jjlo = max0(1,j-iskymedhalf)
        jjhi = min0(nslow,j+iskymedhalf)
        im = 0
        do i=iskymedstep,nfast,iskymedstep
           im = im + 1
           iilo = max0(1,i-iskymedhalf)
           iihi = min0(nfast,i+iskymedhalf)
           k = 1
           do jj=jjlo,jjhi
              do ii=iilo,iihi
                 if (noise(ii,jj).ge.rmagic) cycle
                 values(k) = big(ii,jj)
                 k = k + 1
              end do
           end do
           nvalues = k - 1
           if(nvalues.eq.0) then
              skymedarray(im,jm) = -rmagic
              cycle
           end if
           call median(values,nvalues,skymedarray(im,jm))
        end do
     end do
     if(lverb.ge.10) write(13,*) ' Formed median '


  end if

end subroutine skypar
