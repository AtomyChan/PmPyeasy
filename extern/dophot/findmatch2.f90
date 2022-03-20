! Subroutine findmatch2(starpar,itype,npstar,nsmax,nstot,npap,flags,files,nff,&
!     amagminm,amagmaxm,amagerrmaxm,xoffinit,yoffinit,rmagic,lverb,&
!     radius,dradius,radmin,niter,maxfit,norder,ax,ay,axinv,ayinv,matched)
! This subroutine tries to match two sets of stars, one set from the automatch
! file and the other formed for all the perfect stars found (itype=1)
! It is based on program OFFSET5

subroutine findmatch2(starpar,itype,npstar,nsmax,nstot,npap,flags,files,nff,&
     amagminm,amagmaxm,amagerrmaxm,xoffinit,yoffinit,rmagic,lverb,&
     radius,dradius,radmin,niter,maxfit,norder,ax,ay,axinv,ayinv,matched)

  real, parameter :: factor=3.0,siglimup=0.50,siglimdn=0.20
  real, parameter :: fbadlim=0.01,sigok=0.10
  integer, parameter :: mbadlim=20

  integer :: itype(nsmax)
  real :: starpar(npstar,nsmax)
  real :: matchpar(npstar),apmatchpar(npap)
  real :: xtempl(nsmax),ytempl(nsmax),magtempl(nsmax)
  real :: xtarg(nsmax),ytarg(nsmax),magtarg(nsmax)
  real :: xpostempl(nsmax),ypostempl(nsmax),mmagtempl(nsmax)
  real :: xpostarg(nsmax),ypostarg(nsmax),mmagtarg(nsmax)
  real :: xhold1(nsmax),yhold1(nsmax)
  real :: xhold2(nsmax),yhold2(nsmax)
  real :: dx(nsmax),dy(nsmax)
  real :: ddx(nsmax),ddy(nsmax)
  real :: diffx(nsmax), diffy(nsmax)
  real :: ax(maxfit),ay(maxfit),axinv(maxfit),ayinv(maxfit)
  integer :: indextempl(nsmax),indextarg(nsmax),index3(nsmax)
  integer :: ibad(nsmax)
  character(*) :: flags(nff),files(nff)
  character(80) :: prompt
  logical :: matched,ibadchange

  xmode = 0.0
  ymode = 0.0
  xmodeprev = 0.0
  ymodeprev = 0.0

! Open the matching template file.
  if(flags(10)(1:5).eq.'BINAR') then
     call openab(9,files(9),0,ierr,33)
     if(ierr.eq.1) then
        print *,' Automatch file doesn''t exist on warmstart!'
        prompt = ' Enter automatch file name'
        call opensb(9,prompt,0,33)
     end if
  else
     call opena(9,files(9),0,ierr)
     if(ierr.eq.1) then
        print *,' Automatch file doesn''t exist on warmstart!'
        prompt = ' Enter automatch file name'
        call opens(9,prompt,0)
     end if
  end if

! Now, read the matching template file.
  nstempl = 0
  i1 = 1
  if(flags(10)(1:5).eq.'DAOPH') then
     do i=1,3
        read(9,*)
     end do
  end if
  do
     call choose_in(9,flags(10),idumb,iitype,matchpar,npstar,apmatchpar,npap,&
          amag,probg,ios)
     if (ios.ne.0) exit
     if (.not.((iitype.eq.1).or.(iitype.eq.11)))  cycle
     if (amag.lt.amagminm) cycle
     if (amag.gt.amagmaxm) cycle
     amagerr = apmatchpar(4)
     if (amagerr.gt.amagerrmaxm) cycle
     xtempl(i1) = matchpar(3)+xoffinit
     ytempl(i1) = matchpar(4)+yoffinit
     magtempl(i1) = amag
     i1 = i1+1
     if(i1.gt.nsmax) then
        print *
        print *,' Warning:  Reached AUTOMATCH limit NSMAX.'
        exit
     end if
  end do
  nstempl = i1-1
  if(nstempl.gt.0) then
     print *,nstempl,' stars read from AUTOMATCH file.'
  else
     print *,' NO stars read from AUTOMATCH file; Exiting.'
     call exit(6)      !subroutine exit is an extension to f90
  end if

! Select only perfect stars for the matching from the photometry file.
  i2 = 1
  do k=1,nstot
     if(itype(k).ne.1) cycle
     amag=shapemag(starpar(1,k),npstar)
     if((amag.gt.amagminm).and.(amag.le.amagmaxm)) then
        xtarg(i2) = starpar(3,k)
        ytarg(i2) = starpar(4,k)
        magtarg(i2) = amag
        i2 = i2+1
     end if
  end do
  nstarg = i2-1


  call quick(ytempl,nstempl,indextempl)
  call quick(ytarg,nstarg,indextarg)
  
  dradius = abs(dradius)
  rad = radius+dradius
  do
     rad = amax1(radmin,rad-dradius) 

     k = 1
     jlow2hold = 1
     do i = 1,nstarg
        ylow = ytarg(i)-rad
        yhi  = ytarg(i)+rad
        jlow2 = jlow2hold
        do j = jlow2,nstempl
           if (ytempl(j).gt.yhi) exit
           if (ytempl(j).lt.ylow) then
              jlow2hold = j
              cycle
           end if
           deltax = xtarg(indextarg(i))-xtempl(indextempl(j))
           deltay = ytarg(i)-ytempl(j)
           distance = sqrt((deltax**2)+(deltay**2))
           if (distance.le.rad) then
              diffx(k) = deltax
              diffy(k) = deltay
              xpostarg(k) = xtarg(indextarg(i))
              ypostarg(k) = ytarg(i)
              xpostempl(k) = xtempl(indextempl(j))-xmodeprev
              ypostempl(k) = ytempl(j)-ymodeprev
              mmagtarg(k) = magtarg(indextarg(i))
              mmagtempl(k) = magtempl(indextempl(j))
              k = k+1
           end if
        end do
     end do
     nsmatch = k-1
     if(lverb.gt.10) then
        write(13,*)
        write(13,*)' Radius = ',rad,'  Number of matches found = ',nsmatch
     end if

     call quick(diffx,nsmatch,index3)
     call quick(diffy,nsmatch,index3)
     nmmm = nsmatch
     call mmm(diffx,nmmm,rmagic,xmode,xsigma,xskew)
     call mmm(diffy,nmmm,rmagic,ymode,ysigma,yskew)
     if(lverb.gt.10) then
        write(13,*) 
        write(13,*)'Modal offset in x (targ-templ) = ',xmode,'+/-',xsigma
        write(13,*)'Modal offset in y (targ-templ) = ',ymode,'+/-',ysigma
        write(13,*)
        write(13,*)'Total offset so far in x(targ-templ)=',xmode+xmodeprev
        write(13,*)'Total offset so far in y(targ-templ)=',ymode+ymodeprev
        write(13,*)
     end if

     xmodeprev = xmodeprev+xmode
     ymodeprev = ymodeprev+ymode

     if(rad.eq.radmin) exit

     do i = 1,nstempl
        xtempl(indextempl(i)) = xtempl(indextempl(i))+xmode
        ytempl(i) = ytempl(i)+ymode
     end do
  end do

  iiter = 1
  ibad = 0
  mfit = MGETMFIT(norder)

  do
     ibadchange = .false.
     
     ii = 1
     do i=1,nsmatch
        if(ibad(i).eq.0) then
           xhold1(ii) = xpostarg(i)
           yhold1(ii) = ypostarg(i)
           xhold2(ii) = xpostempl(i)
           yhold2(ii) = ypostempl(i)
           ii = ii+1
        end if
     end do
     ntot = ii-1

     nbadlim = max0(mbadlim,nint(fbadlim*float(ntot)))

     if(mfit.gt.1) then
        call surfpoly(xhold1,xhold2,yhold2,ntot,ax,mfit)
        call surfpoly(yhold1,xhold2,yhold2,ntot,ay,mfit)
        call surfpoly(xhold2,xhold1,yhold1,ntot,axinv,mfit)
        call surfpoly(yhold2,xhold1,yhold1,ntot,ayinv,mfit)
     else
        ax(1) = xmodeprev
        ay(1) = ymodeprev
        axinv(1) = -xmodeprev
        ayinv(1) = -ymodeprev
     end if

     ii = 1
     do i=1,nsmatch
        if(ibad(i).eq.0) then
           dx(i) = xpostarg(i)-SURFEVAL(xpostempl(i),ypostempl(i),ax,mfit)
           dy(i) = ypostarg(i)-SURFEVAL(xpostempl(i),ypostempl(i),ay,mfit)
           ddx(ii) = dx(i)
           ddy(ii) = dy(i)
           ii = ii + 1
        end if
     end do     
     ndd = ii-1
     
     call average(ddx,ndd,dxavg,dxsig,dxsabs)
     call average(ddy,ndd,dyavg,dysig,dysabs)

     if(lverb.gt.10) then
        write(13,*)
        write(13,*) ' Iteration, ndd = ',iiter,ndd
        write(13,*)
        write(13,*) ' x,y standard deviations = ',dxsig,dysig
        write(13,*) ' x,y standard abs deviations = ',dxsabs,dysabs
        write(13,*)
     end if

!    See if any points should be removed.
     if(dxsig.gt.sigok.or.dysig.gt.sigok) then
        if(iiter.lt.niter) then
           xlim = amin1(siglimup,amax1(factor*dxsig,siglimdn))
           ylim = amin1(siglimup,amax1(factor*dysig,siglimdn))
           nbad = 0
           do i=1,nsmatch
              if((ibad(i).eq.0).and.(.not.ibadchange)&
                   .and.((abs(dx(i)).gt.xlim).or.(abs(dy(i)).gt.ylim))) then
                 ibad(i) = 1
                 nbad = nbad+1
                 if (nbad.gt.nbadlim) ibadchange=.true.
              end if
           end do
        else
           return
        end if
     else
        matched = .true.
        exit
     end if
     
     if(ibadchange) then 
        iiter = iiter + 1
        if(mfit.le.1) exit
     end if
  end do

!  Correct coefficients for any initial offsets so I don't have to mess
!  with this later.
  ax(1) = ax(1) + xoffinit
  ay(1) = ay(1) + yoffinit
  axinv(1) = axinv(1) - xoffinit
  ayinv(1) = ayinv(1) - yoffinit
  if(flags(14)(1:3).eq.'YES') then
     call opena(14,files(14),2,ierr)
     write(14,*) ' '
     write(14,*) ' Nsmatch, Norder, Ncoefs = ',nsmatch-nbad,norder,mfit
     write(14,*) ' x,y standard deviations = ',dxsig,dysig
     write(14,*) ' x,y standard abs deviations = ',dxsabs,dysabs
     do j=1,2
       write(14,*) ' ' 
       if(j.eq.1) write(14,*)'Forward transformation; Template to Target (x,y):'
       if(j.eq.2) write(14,*)'Inverse transformation; Target to Template (x,y):'
       write(14,*) ' '
       do i=1,mfit
          if(j.eq.1) write(14,*) ax(i),ay(i)
          if(j.eq.2) write(14,*) axinv(i),ayinv(i)
       end do
     end do
     write(14,*) ' ' 
  end if

end subroutine findmatch2

