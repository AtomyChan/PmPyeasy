! Subroutine warmstart(big,noise,nfast,nslow,flags,files,nff,&
!     ava,npstar,nsmax,npap,nstot,&
!     onestar,par1,par2,fobl,isub,fsub,xpnd,fac,rmagic,lverb,&
!     starpar,itype,shadow,shaderr)
! This subroutine provides DoPHOT with a starting list of objects. 
! The position and shapes of these objects can be fixed 
! or can be allowed to change.

subroutine warmstart(big,noise,nfast,nslow,flags,files,nff,&
     ava,npstar,nsmax,npap,&
     onestar,par1,par2,fobl,isub,fsub,xpnd,fac,rmagic,lverb,&
     nstot,starpar,itype,shadow,shaderr)

  real :: big(nfast,nslow),noise(nfast,nslow)
  real :: starpar(npstar,nsmax),appar(npstar,nsmax)
  real :: shadow(npstar,nsmax),shaderr(npstar,nsmax)
  real :: summ(npstar),ava(npstar)
  integer :: itype(nsmax)
  character(80) :: flags(nff),files(nff),prompt
  external onestar

  if(flags(6)(1:5).eq.'BINAR') then
     call openab(3,files(3),0,ierr,33)
     if(ierr.eq.1) then
        print *,' Input object file doesn''t exist on warmstart!'
        prompt = ' Enter input object file name'
        call opensb(3,prompt,0,33)
     end if
  else
     call opena(3,files(3),0,ierr)
     if(ierr.eq.1) then
        print *,' Input object file doesn''t exist on warmstart!'
        prompt = ' Enter input object file name'
        call opens(3,prompt,0)
     end if
  end if
  
  if(flags(7)(1:1).eq.'Y') then
     call opena(7,files(7),0,ierr)
     if(ierr.eq.1) then
        print *,' Input shadow file doesn''t exist on warmstart!'
        prompt = ' Enter input shadow file name'
        call opens(7,prompt,0)
     end if
  end if
  
  i = 1
  nstperfect = 0
  nstot = 0
  do
     call choose_in(3,flags(6),i1,i2,starpar(1,i),npstar,appar(1,i),npap,&
          dum,probg,ios)
     if (ios.ne.0) exit
     
     if(flags(7)(1:1).eq.'Y') then
        read (7,*,iostat=io)  j1,j2,(shadow(j,i),j=1,npstar)
        if(shadow(2,i).gt.0) then
           temp = (shadow(2,i)+shadow(1,i))/shadow(2,i)**2
        else
           temp = rmagic
        end if
        do j = 1,npstar
           shaderr(j,i) = temp
        end do
        if((j1.ne.i1).or.(j2.ne.i2)) then
           write(13,*) ' Forced exit in WARMSTART:'
           write(13,*) ' Object and Shadow file lengths unequal.'
           write(13,*) ' J1.NE.I1 OR J2.NE.I2 ', J1,I1,J2,I2
           call exit(3)        ! Subroutine exit is an extension to f90
        end if
     end if
     
     itype(i) = i2
     if(itype(i).eq.1) then
        nstperfect = nstperfect+1
        do  j = 1,npstar
           summ(j) = summ(j)+starpar(j,i)
        end do
     end if
     if((itype(i).ne.6).and.(itype(i).ne.8).and.(itype(i).ne.9)) then
        call addstar(isub,starpar(1,nstot),npstar,onestar,par1,par2,&
             fsub,xpnd,fac,rmagic,nfast,nslow,big,noise)
     else if(itype(i).eq.8) then
        call oblit(starpar(1,i),npstar,'r',fobl,.false.,lverb,rmagic,&
             nfast,nslow,big,noise)
     end if
     i = i + 1
  end do
  nstot = i-1
  
  if ((lverb.gt.10).and.(nstot.gt.0)) then
     write(13,*)
     write(13,*) nstot,' objects from WARMSTART'
  end if
  
  if(nstperfect.ne.0) then
     do i = 1,npstar
        ava(i) = summ(i)/float(nstperfect)
     end do
  end if

end subroutine warmstart
