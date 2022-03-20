! Dophot classifies the stars in 9 types:
!
! Type 1: A "perfect" star. Type 1 objects are used to compute the 7-parameter
!         weighted-mean stellar shape used to extract photometric data 
!         from an image.
! Type 2: The object is not as peaked as a star and has been classified as 
!         a "galaxy".
! Type 3: A blend of 2 close stars, with the data being for one component of 
!         the pair. The final fit uses mean shape parameters for single stars 
!         determined from all type 1 objects on the image.
! Type 4: Fit failed to converge on a 4-parameter fit based on the mean 
!         single-star shape. Photometry is unreliable.
! Type 5: Not enough points with adequate S/N for a full 7-parameter fit, 
!         so a 4-parameter fit attempted and succeeded. Photometry is suspect, 
!         and may indicate image problems or that the object is not a star.
! Type 6: Too few points with adequate S/N for even a 4-parameter fit with 
!         the mean Type 1 star shape. Not subtracted from the image.
! Type 7: Object too faint for a full 7-parameter mean-star fit. No object 
!         classification is possible. The results are for a 4-parameter fit, 
!         and are reasonable only if the object really is a star.
! Type 8: Object was "obliterated" becase it was determined to be saturated, 
!         either above the pre-set threshold or found to be flat-topped.
! Type 9: Center of the star out of the image. Photometry is unreliable.
!         It is also used for saturated stars when not obliterated.
! In addition to this types, stars from a warmstart with positions and/or
! shapes fixed are given types equal to the ones described plus 10.
!
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!
! Subroutine shappe(big,noise,nfast,nslow,onestar,twostar,par1,par2,&
!     fixpos,fixshape,nfix,NSTOT,NSNEW,nsperf,nsrmsmin,nsmax,npstar,lverb,&
!     flags,nff,enuffvpsf,a5,a6,a7,npsffit,ava,acc,alim,nit,&
!     iadd,isub,fsub,xpnd,fac,rmagic,snlim,fitrad,maxfil,maxval,&
!     enuff7,crit7,ixby2,iyby2,sig,xtra,chicrit,nfit2,stograt,&
!     > STARPAR,SHADOW,SHADERR,ITYPE)
! This subroutine classifies the stars according to the types described above.
! It takes care specially of finding stars that are too 'big' and tries to
! discriminate between galaxies and blends of stars (types 2 and 3).
! To better understand this subroutine read section 4.4 and 
! paragraph 5 of section 3 of the DoPHOT paper (Schechter, Mateo, and Saha 1993
! PASP, 105, 1345).

subroutine shappe(big,noise,nfast,nslow,onestar,twostar,par1,par2,&
     fixpos,fixshape,nfix,nstot,nsnew,nsperf,nsrmsmin,nsmax,npstar,lverb,&
     flags,nff,enuffvpsf,a5,a6,a7,npsffit,ava,acc,alim,nit,&
     iadd,isub,fsub,xpnd,fac,rmagic,snlim,fitrad,maxfil,maxval,&
     enuff7,sig,xtra,chicrit,nfit2,stograt,&
     starpar,shadow,shaderr,itype)

  integer :: itype(nsmax)
  real :: starpar(npstar,nsmax),shadow(npstar,nsmax),shaderr(npstar,nsmax)
  real :: big(nfast,nslow),noise(nfast,nslow)
  real :: xy(2,maxfil),z(maxfil), ze(maxfil)
  real :: a(npstar),b1(npstar),b2(npstar),dummy(npstar)
  real :: ava(npstar),acc(npstar),alim(npstar),starparrms(npstar)
  real :: chi(4),sig(3)
  real :: a5(npsffit),a6(npsffit),a7(npsffit)
  character(*) :: flags(nff)
  logical :: fixpos,fixshape,fixxy
  logical :: enuffvpsf,enuffpts
  logical :: verybig,converge,offp,offpic
  external onestar,twostar


! starparrms is used later by chisqshape
  starparrms = 0.
  if (nsperf.ge.nsrmsmin) then
     call stparrms(nstot,nsmax,npstar,starpar,itype,shadow,shaderr,ava,&
       flags,nff,lverb,enuffvpsf,a5,a6,a7,npsffit,starparrms)
  end if

  nsprev = nstot
  do i=1,nsprev
     if(lverb.gt.20) write(13,*) ' Determining SHAPE for object No.', I
      do k = 1,npstar
         if (shadow(1,i).eq.0) then
            a(k) = starpar(k,i)
         else
            a(k) = shadow(k,i)
         end if
      end do
     if ((itype(i).eq.7).and.(lverb.gt.20)) write(13,*) &
          ' TOO FAINT, SKIPPING Object # ',i,' AT: ',(starpar(j,i),j=3,4)
     if ((itype(i).eq.8).and.(lverb.gt.20)) write(13,*) &
          ' TOO BRIGHT, SKIPPING Object # ',i,' AT: ',(starpar(j,i),j=3,4)
     if ((itype(i).eq.9).and.(lverb.gt.20)) write(13,*) &
          ' CENTER OUT IMG., SKIPPING Object #: ',i,' AT: ',(starpar(j,i),j=3,4)
     if (((itype(i).eq.4).or.(itype(i).eq.6)).and.(lverb.gt.20)) write(13,*) &
          ' NONCONVERGENCE, SKIPPING Object #: ',i,' AT: ',(starpar(j,i),j=3,4)
     if((itype(i).le.3).or.(itype(i).eq.5)) then 
        call addstar(iadd,starpar(1,i),npstar,onestar,par1,par2,&
             fsub,xpnd,fac,rmagic,nfast,nslow,big,noise)
        call fillerup (a(3),a(4),big,noise,nfast,nslow,snlim,fitrad,rmagic,&
             maxfil,npstar,xy,z,ze,npt,areafill,areatot,dummy)
        enuffpts = areafill.ge.enuff7*areatot
        if (.not.enuffpts) then           
           if(lverb.gt.20) write(13,*) &
                ' Obj#, NPTS FIT, IX, IY = ',i,npt,ix,iy,&
                '   .... SKIPPIN STAR: not enough pixels for 7-parm fit '
           if (itype(i).ne.2) itype(i) = 5             
           call addstar(isub,starpar(1,i),npstar,onestar,par1,par2,&
                fsub,xpnd,fac,rmagic,nfast,nslow,big,noise)
           cycle
        end if

        if(lverb.gt.20) write(13,*) &
             ' Obj#, NPTS FIT, IX, IY = ',i,npt,ix,iy
        call chisq(onestar,par1,par2,xy,z,ze,npt,a,npstar,&
             nfit2,acc,alim,nit,rmagic,lverb,shaderr(1,i),galchi)
        do k = 1,npstar
           shadow(k,i) = a(k)
        end do
        converge = galchi.lt.rmagic
        offp = OFFPIC(a(3),a(4),nfast,nslow)
        chisqshappe = CHISQSHAPE(shadow(1,i),shaderr(1,i),starparrms,npstar,&
             sig,lverb,flags,nff,enuffvpsf,a5,a6,a7,npsffit,ava)
        verybig = (chisqshappe.ge.chicrit)
        if(itype(i).eq.3) verybig = verybig.and.(chisqshappe.gt.xtra)
        fixxy=((fixpos.or.fixshape).and.(i.le.nfix))

!       If the object position or shape is fixed then we don't test 
!       for duplicity.  
!       If it failed to converge, it is so flagged (using the shadow array),
!       and its original status as type 3 (if possessed) is preserved.
        if ((.not.verybig).or.(.not.converge).or.offp.or.fixxy) then
           if (itype(i).ne.3) itype(i) = 1
           if (verybig.and.fixxy) itype(i) = 2
           if (.not.converge) then
              if (itype(i).ne.3) itype(i) = 4
              shadow(1,i) = 0                       
              if(lverb.gt.20)  write(13,*) &
                   'Obj #',i,' at ',(starpar(j,i),j=3,4),&
                   'FAILED TO CONVERGE! ' 
           else if (offp) then
              itype(i) = 9                            
              shadow(1,i) = 0                       
              if(lverb.gt.20) then
                 write(13,*) & 
                      ' ABSURD SHAPE VALUES for Obj # ',i,' at: ',&
                      (starpar(j,i),j=3,4)
                 write(13,*) &
                      ' Fit center outside fit subraster ... DISCARD SOLUTION!'
              end if
           end if
           call addstar(isub,starpar(1,i),npstar,onestar,par1,par2,&
                fsub,xpnd,fac,rmagic,nfast,nslow,big,noise)
           cycle
        end if
           
        if(lverb.gt.20) then
           write(13,*) ' Obj # ',i,' at: ',(starpar(j,i),j=3,4),' is VERY BIG..'
           write(13,*) ' ... Testing GALAXY vs. DBLE-STAR ...'
        end if
        call twofit(starpar(1,i),a,npstar,twostar,par1,par2,&
             flags,nff,enuffvpsf,a5,a6,a7,npsffit,ava,&
             rmagic,lverb,fitrad,xy,z,ze,maxfil,npt,nfit2,acc,alim,nit,&
             b1,b2,starchi)
        if (starchi/galchi.lt.stograt) then   
           if(lverb.gt.20) write(13,*) ' Result-> A SPLIT STAR:',&
                ' GAL-CHI:', GALCHI, ' STAR-CHI:', STARCHI
           itype(i) = 3
           do k = 1,npstar
              starpar(k,i) = b1(k)
              shadow(k,i) = b1(k)
           end do
           call addstar(isub,starpar(1,i),npstar,onestar,par1,par2,&
                fsub,xpnd,fac,rmagic,nfast,nslow,big,noise)
           nstot = nstot+1
           nsnew = nsnew+1
           itype(nstot) = 3             
           do k = 1,npstar
              starpar(k,nstot) = b2(k)
              shadow(k,nstot) = b2(k)
           end do
           call addstar(isub,starpar(1,nstot),npstar,onestar,par1,par2,&
                fsub,xpnd,fac,rmagic,nfast,nslow,big,noise)
        else
           if(lverb.gt.20) write(13,*) ' Result-> A GALAXY:', &
                ' GAL-CHI:', GALCHI, ' STAR-CHI:', STARCHI
           itype(i) = 2
           do k = 1,npstar
              starpar(k,i) = a(k)
           end do
           call addstar(isub,starpar(1,i),npstar,onestar,par1,par2,&
                fsub,xpnd,fac,rmagic,nfast,nslow,big,noise)
        end if
     end if
  end do
  
end subroutine shappe






