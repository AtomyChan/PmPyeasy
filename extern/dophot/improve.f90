! The aperture information is kept in an array of 4 elements
! apar(1) keeps the aperture flux from the object
! apar(2) keeps the aperture flux from the sky
! apar(3) is the aperture correction magnitude (fit minus aperture)
! apar(4) is the fit magnitude variance
!
!---------------------------------------------------------------------------
!
! Subroutine improve(big,noise,nfast,nslow,nstot,fixpos,fixshape,nfix,&
!    nfit0,nfit1,nit,acc,alim,apradius,enuff4,enuff7,apmagmaxerr,eperdn,rnoise,&
!    nsperf,fsub,xpnd,fac,onestar,par1,par2,flags,nff,&
!    ava,a5,a6,a7,npsffit,enuffvpsf,lverb,iadd,isub,&
!    maxfil,snlim,fitrad,rmagic,cmin,crit7,nsmax,npstar,npap,&
!    > STARPAR,ITYPE,APPAR,PROBG)
! This subroutine decides (again) the type of the star (types defined in 
! subroutine shape), this time without altering the shape parameters. 
! In addition to this, the aperture photometry is obtained and 
! the aperture parameters described above are obtained, 
! along with an extendness parameter that can help ex post facto delineation 
! of stars from galaxies. 
! To better understand this subroutine, see second to last paragraph of 
! section 3 in DoPHOT paper (Schechter, Mateo, and Saha 1993 PASP, 105, 1345))
  
subroutine improve(big,noise,nfast,nslow,nstot,fixpos,fixshape,nfix,&
     nfit0,nfit1,nit,acc,alim,enuff4,enuff7,eperdn,rnoise,&
     apmagmaxerr,apradius,apskymin,apskymax,&
     nsperf,fsub,xpnd,fac,onestar,par1,par2,flags,nff,&
     ava,a5,a6,a7,npsffit,enuffvpsf,lverb,iadd,isub,&
     maxfil,snlim,fitrad,rmagic,cmin,crit7,nsmax,npstar,npap,&
     starpar,itype,appar,probg)

  real,parameter :: pi=3.14159
  integer :: itype(nsmax)
  real :: appar(npap,nsmax),probg(nsmax),starpar(npstar,nsmax)
  real :: big(nfast,nslow),noise(nfast,nslow)
  real :: xy(2,maxfil),z(maxfil),ze(maxfil)
  real :: a(npstar),err(npstar),dummy(npstar)
  real :: a5(npsffit),a6(npsffit),a7(npsffit)
  real :: ava(npstar),acc(npstar),alim(npstar)
  logical :: fixpos,fixshape,enuffpts,enuffvpsf,offpic,converge,toofaint,dumlog
  character(*) :: flags(nff)
  external onestar

  nsperf = 0
  do i=1,nstot
     if ((itype(i).eq.0).or.(itype(i).eq.2).or.(itype(i).eq.6).or.&
          (itype(i).eq.8)) cycle
     call addstar(iadd,starpar(1,i),npstar,onestar,par1,par2,&
          fsub,xpnd,fac,rmagic,nfast,nslow,big,noise)
     do k = 1,npstar
        a(k) = starpar(k,i)
     end do
     if (.not.(fixshape.and.(i.le.nfix))) &
          call stpar567 (a(3),a(4),flags,nff,enuffvpsf,a5,a6,a7,npsffit,ava,&
          npstar,a)
     if(lverb.gt.20) write(13,*) ' IMPROVING STAR # ',i,' AT ',ix,iy
 
     call fillerup (a(3),a(4),big,noise,nfast,nslow,snlim,&
          fitrad,rmagic,maxfil,npstar,xy,z,ze,npt,areafill,areatot,dummy)
     if (itype(i).eq.1) then
        enuffpts = areafill.ge.enuff7*areatot
        if (.not.enuffpts) itype(i) = 5
     end if
     enuffpts = areafill.ge.enuff4*areatot
     if (.not.enuffpts) then
        itype(i) = 6        
        if(lverb.gt.20) write(13,*) &
             ' DEACTIVATING STAR # ',i, ' AT IX, IY = ',ix,iy
        cycle
     end if
     nfit = nfit1
     if (fixpos.and.(i.le.nfix)) nfit = nfit0
     call chisq(onestar,par1,par2,xy,z,ze,npt,a,npstar,&
          nfit,acc,alim,nit,rmagic,lverb,err,starchi2)
     converge = starchi2.lt.rmagic
 
     if(.not.converge) then
        if (itype(i).ne.3) itype(i) = 4
     else if (OFFPIC(a(3),a(4),nfast,nslow)) then
        itype(i) = 9
     else
        do k = 1,npstar
           starpar(k,i) = a(k)
        end do
        if ((itype(i).ne.3).and.&
             TOOFAINT(starpar(1,i),err,npstar,cmin,crit7,lverb))&
             itype(i)=7
        totstar = starpar(2,i)*ELLIPAREA(starpar(5,i),starpar(6,i),starpar(7,i))
        varstar = totstar*eperdn
        varsky = (pi*apradius**2)*(starpar(1,i)*eperdn+rnoise**2)
        errelec = sqrt(varstar+varsky)
        errdn = errelec/eperdn
        snpred = 1.086*errdn/totstar
        if (snpred.lt.apmagmaxerr) &
             call aperpar(big,noise,nfast,nslow,starpar(1,i),npstar,rmagic,&
             apradius,apskymin,apskymax,0,0,0,0,1,npap,appar(1,i),dum,dumlog)
        appar(4,i) = 1.086*sqrt(err(2))/a(2)    !dm=2.5d(logf)=2.5/ln(10)/f*df
        if(abs(appar(4,i)).gt.100.0) appar(4,i) = 99.999
        probg(i) = PROBGAL(xy,z,ze,npt,npstar,a,onestar,par1,par2,lverb)
     end if

     if(starpar(5,i).le.0.0.or.starpar(7,i).le.0.0) then
        itype(i) = 6        
        if(lverb.gt.20) write(13,*) &
             ' DEACTIVATING STAR #', i, ' AT IX, IY = ',ix,iy
     else
        call addstar(isub,starpar(1,i),npstar,onestar,par1,par2,&
             fsub,xpnd,fac,rmagic,nfast,nslow,big,noise)
     end if
     if(itype(i).eq.1) nsperf = nsperf+1
  end do
  
end subroutine improve
