! Subroutine finish(files,flags,nff,ava,npstar,axinv,ayinv,maxfit,itype,nstot,&
!     norder,timediff,elaptime,date1)
! This subroutine saves the information about the program run

subroutine finish(files,flags,nff,ava,npstar,axinv,ayinv,maxfit,itype,nstot,&
     norder,timediff,elaptime,date1)

  character(*) :: date1,files(nff),flags(nff)
  integer :: itype(nstot),itype2(9)
  real :: axinv(maxfit),ayinv(maxfit),ava(npstar)
  character(2) :: fmtvar
  character(120) :: fmt
  logical :: fixpos

  nfile = len_trim(files(15))
  if(nfile.eq.0) return
       
  call opena(15,files(15),3,ierr)

! Count the number of star with equal itype
  itype2=0
  do i=1,nstot
     j = mod(itype(i),10)
     itype2(j) = itype2(j)+1
  end do

  call ellipse(ava(5),ava(6),ava(7),area,fwhmmaj,fwhmmin,angle)
  fwhm = (fwhmmaj+fwhmmin)/2.0

  nimage = len_trim(files(1))
  ndate1 = len_trim(date1)

  write(15,*)
  write(15,'(3a)') files(1)(1:nimage),' - Completed reduction on ',&
       date1(1:ndate1)

  elaptime = elaptime/3600.0
  write(15,'(2a,f5.3,a)') files(1)(1:nimage),&
       ' - Approximate total elapsed time = ',elaptime,' hours'

  timediff = timediff/60.0
  write(15,'(2a,f5.3,a)') files(1)(1:nimage),&
       ' - Approximate total CPU time = ',timediff,' minutes'

  rate = float(nstot)/timediff
  write(15,'(2a,f8.1,a)') files(1)(1:nimage),&
       ' - Approximate reduction rate = ',rate,' stars/minute'

  write(15,'(2a,i7,a,f5.2,a,f8.1)') files(1)(1:nimage),' - NSTOT = ',nstot,&
       ' FWHM = ',fwhm,' <SKY> = ',ava(1)

  write(15,'(2a,9i7)') files(1)(1:nimage),' - CENSUS: ',(itype2(k),k=1,9)

  if(flags(15)(1:1).eq.'Y') then
     mfit = MGETMFIT(norder)
     write(fmtvar,'(i2)') mfit
     fmt = "(2a,i2,"//fmtvar//"e10.2)"
     write(15,fmt) files(1)(1:nimage),'- X COEFS: ',mfit,(axinv(k),k=1,mfit)
     write(15,fmt) files(1)(1:nimage),'- Y COEFS: ',mfit,(ayinv(k),k=1,mfit)
  end if

  write(15,*)
  
end subroutine finish

