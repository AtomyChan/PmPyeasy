! Subroutine isearch(big,noise,nfast,nslow,thresh,pixthresh,rmagic,&
!     flags,nff,skyplanepar,npskyp,skyhubpar,npskyh,iskymedstep,skymedarray,&
!     enuffvpsf,a5,a6,a7,npsffit,ava,nsmax,npstar,&
!     onestar,par1,par2,maskrad,&
!     crit7,bumpcrit,lverb,&
!     snlim,enuff4,fitrad,maxfil,nfit1,acc,alim,nit,&
!     top,topsat,cmax,icrit,fobl,cmin,widobl,discrim,sncos,&
!     fsub,xpnd,isub,fac,&
!     nstot,starpar,itype,nsnew)
! This subroutine looks for new objects in the image in a way explained in 
! paragraph 4 of section 3 of the DoPHOT paper (Schechter, Mateo, and Saha 1993
! PASP, 105, 1345).

subroutine isearch(big,noise,nfast,nslow,thresh,pixthresh,rmagic,&
     flags,nff,skyplanepar,npskyp,skyhubpar,npskyh,iskymedstep,skymedarray,&
     enuffvpsf,a5,a6,a7,npsffit,ava,nsmax,npstar,&
     onestar,par1,par2,maskrad,&
     crit7,bumpcrit,lverb,&
     snlim,enuff4,fitrad,maxfil,nfit1,acc,alim,nit,&
     top,topsat,cmax,icrit,fobl,cmin,widobl,discrim,sncos,&
     fsub,xpnd,isub,fac,&
     nstot,starpar,itype,nsnew)

  real :: big(nfast,nslow),noise(nfast,nslow)
  real :: skymedarray(nfast/iskymedstep,nslow/iskymedstep)
  real :: ava(npstar),skyplanepar(npskyp),skyhubpar(npskyh)
  real :: starmask(-maskrad:maskrad,-maskrad:maskrad)
  real :: a5(npsffit),a6(npsffit),a7(npsffit)
  real :: a(npstar),err(npstar),acc(npstar),alim(npstar)
  real :: xy(2,maxfil),z(maxfil),ze(maxfil)
  character(*) :: flags(nff)
  logical:: enuffvpsf,transmask,enuffpts,offpic
  logical :: cosm,toobright,toob,wipe,toofaint
  real :: starpar(npstar,nsmax)
  integer :: itype(nsmax)
  external onestar

  nsprev = nstot
  do iy = 1,nslow
     y = float(iy)
     do ix = 1,nfast
        x = float(ix)
        if (noise(ix,iy).ge.amin1(rmagic,(thresh/pixthresh)**2)) cycle
        bestsky=STPAR1(x,y,flags,nff,skyplanepar,npskyp,skyhubpar,npskyh,&
             skymedarray,nfast,nslow,iskymedstep)
        highthresh = bestsky+thresh
        if (big(ix,iy).lt.highthresh) cycle 
        a = 0.
        call makemask(x,y,flags,nff,enuffvpsf,a5,a6,a7,npsffit,ava,npstar,&
             onestar,par1,par2,maskrad,starmask)
        if (.not.(TRANSMASK(x,y,bestsky,starmask,maskrad,&
             big,noise,nfast,nslow,bumpcrit,rmagic,lverb))) cycle
        if(lverb.gt.20) &
             write(13,*) ' Triggering on new object:',ix,iy,big(ix,iy)        
        call fillerup (x,y,big,noise,nfast,nslow,snlim,fitrad,rmagic,&
             maxfil,npstar,xy,z,ze,npt,areafill,areatot,a)
        enuffpts = areafill.ge.enuff4*areatot
        if (.not.enuffpts) then
           if(lverb.gt.20) write(13,*) 'Skipping: npt = ', npt
           cycle
        end if
        if(.not.TRANSMASK(x,y,a(1),starmask,maskrad,&
             big,noise,nfast,nslow,bumpcrit,rmagic,lverb)) then
           if(lverb.gt.20) write(13,*) 'Failed transmask on average sky'
           cycle
        end if
       call stpar567(x,y,flags,nff,enuffvpsf,a5,a6,a7,npsffit,ava,npstar,a)
        err=0.
        call chisq(onestar,par1,par2,xy,z,ze,npt,a,npstar,&
             nfit1,acc,alim,nit,rmagic,lverb,err,chi2)
        if(chi2.ge.rmagic) then
           if(lverb.gt.20) write(13,*) 'Failed to converge.'
           cycle
        end if
       
        nstot = nstot+1
        if(nstot.gt.nsmax) then
           print *
           print *,'The next star would overflow star arrays.'
           print *,'Reductions complete to threshold = ',thresh
           print *,'To continue, resize arrays and do warmstart.'
           call exit(7)     ! subroutine exit is an ext to f90
        end if
        if(lverb.gt.20) write(13,*) 'This is star no. ',nstot
       
        do k = 1,npstar
           starpar(k,nstot) = a(k)
        end do
       
        if (OFFPIC(a(3),a(4),nfast,nslow)) itype(nstot) = 9

        call cosmic(big,noise,nfast,nslow,onestar,par1,par2,npstar,&
             widobl,discrim,sncos,rmagic,lverb,&
             starpar(1,nstot),cosm)
        toob = TOOBRIGHT(starpar(1,nstot),npstar,big,noise,nfast,nslow,&
             top,fitrad,cmax,icrit,rmagic,lverb)
        if (cosm.or.(toob.and.(flags(16)(1:1).eq.'Y'))) then
           itype(nstot) = 8
           call oblit(starpar(1,nstot),npstar,'c',fobl,cosm,lverb,rmagic,&
                nfast,nslow,big,noise)
           cycle
        else if (toob.and.(flags(16)(1:1).eq.'N')) then
           itype(nstot) = 9
           call oblit2(starpar(1,nstot),npstar,topsat,fitrad,lverb,rmagic,&
                nfast,nslow,big,noise)
       end if

        if (TOOFAINT(starpar(1,nstot),err,npstar,cmin,crit7,lverb)) &
             itype(nstot) = 7   

        call addstar(isub,starpar(1,nstot),npstar,onestar,par1,par2,&
             fsub,xpnd,fac,rmagic,nfast,nslow,big,noise)
        
     end do
  end do
  
  nsnew = nstot-nsprev

end subroutine isearch






