! Program fitsubs
! Works with FITS files, using the subroutines from the CFITSIO library
!
! Includes these subroutines:
!  OPENIMAGE
!  CLOSEIMAGE
!
!  GETDIMEN
!  GETHEADER1
!  GETHEADER2
!  PUTHEADER
!  GETITEMR
!  GETITEMI
!  GETITEMC
!!  GETTIMEINFO
!!  GETDATE
!!  GETUT
!
!  GETDATAI2
!  GETDATAI4
!  GETDATAR4
!
!  PUTDATAI2
!  PUTDATAI4 
!  PUTDATAR4
!
! It compiles with the libcfitsio.a library

!-----------------------------------------------------------------------

subroutine openimage(imagename,iunit,status)

  character(*) :: imagename
  integer :: status,iunit,readwrite,blocksize
  
  status = 0
  readwrite = 1
  nimagename=len_trim(imagename)
  call ftopen(iunit,imagename(1:nimagename),readwrite,blocksize,status)
  if(status.ne.0) then
     status = 0 
     call ftinit(iunit,imagename(1:nimagename),blocksize,status)
  end if
  
end subroutine openimage

!-----------------------------------------------------------------------

subroutine closeimage(iunit,status)
  
  integer :: status
  
  status=0
  call ftclos(iunit,status)
  
end subroutine closeimage

!-----------------------------------------------------------------------

subroutine getdimen(iunit,n1,n2) 
  
  integer :: iunit,n1,n2,naxes(2),nfound,status
  
  status=0
  call ftgknj(iunit,'NAXIS',1,2,naxes,nfound,status)
  if(nfound.ne.2) then
     print *,'Problem; naxes is not 2!'
     stop
  end if
  
  n1 = naxes(1)
  n2 = naxes(2)
  
end subroutine getdimen

!-----------------------------------------------------------------------

subroutine getheader1(iunit,nheader)
  
  integer :: iunit
  integer :: status,nheader,nspace
  
  status = 0
  call ftghsp(iunit,nheader,nspace,status)
  
end subroutine getheader1

!-----------------------------------------------------------------------

subroutine getheader2(iunit,nheader,header)
  
  integer :: iunit
  integer :: status,nheader
  character(80) :: record,header(nheader)
  
  status = 0
  do i = 1,nheader
     call ftgrec(iunit,i,record,status)
     header(i) = record
  end do
  
end subroutine getheader2

!-----------------------------------------------------------------------

subroutine putheader(iunit,header,nheader)
  
  integer :: iunit,status
  character(*) :: header(nheader)
  character(80) :: record
  logical :: endpresent
  
! Check if there is already a header on image.
! If there is, exit
  status=0
  call ftghsp(iunit,keysexist,keysadd,status)
  if (keysadd.ne.-1) return

! Add header to image
  status=0
 
  do i=1,nheader
     record = header(i)
     call ftprec(iunit,record,status)
  end do
  
! Check for last record.
  nrecord = len_trim(record)
  endpresent = (record(1:3).eq.'END').and.(nrecord.eq.3)
  if(.not.endpresent) then
     do i=1,80
        record(i:i) = ' '
     end do
     record(1:3) = 'END'
     call ftprec(iunit,record,status)
  end if
  
end subroutine putheader

!-----------------------------------------------------------------------

subroutine getitemr(iunit,keyword,item,comment,status)
  
  integer :: status
  character(*) :: keyword, comment
  real :: item
  
  status=0
  call ftgkye(iunit,keyword,item,comment,status)
  
end subroutine getitemr

!-----------------------------------------------------------------------

subroutine getitemi(iunit,keyword,item,comment,status)
  
  integer :: status,item
  character(*) :: keyword, comment
  
  status=0
  call ftgkyj(iunit,keyword,item,comment,status)
  
end subroutine getitemi

!-----------------------------------------------------------------------

subroutine getitemc(iunit,keyword,item,comment,status)
  
  integer :: status
  character(*) :: keyword, comment, item
  
  status=0
  call ftgkys(iunit,keyword,item,comment,status)
  
end subroutine getitemc
!-----------------------------------------------------------------------
!
!! DANGER: The names of these parameters can change in the header
!! depending on where the images were taken at.
!
!  subroutine gettimeinfo(iunit,date,ut,etime,gain,ron)
!
!    real :: date(1),ut(1),etime,gain,ron,ncombine
!    integer :: iobs,status,ipic
!    character(80) :: comment,utstring,datestring
!
!    status=0
!    call ftgkye(iunit,'EXPTIME',etime,comment,status)
!    if(status.ne.0) then
!       print*,'Problems in the subrotine gettimeinfo; check exptime'
!       stop
!    end if
!
!    call ftgkys(iunit,'TIME-OBS',utstring,comment,status)
!    if(status.ne.0) then
!       status = 0
!       call ftgkys(iunit,'UT-TIME',utstring,comment,status)
!    end if
!    if(status.ne.0) then
!       print*,'Problems in the subrotine gettimeinfo; check time-obs'
!       stop
!    end if
!
!    call ftgkys(iunit,'DATE-OBS',datestring,comment,status)
!    if(status.ne.0) then
!       print*,'Problems in the subrotine gettimeinfo; check date-obs'
!       stop
!    end if
!    
!    call ftgkye(iunit,'ENOISE',ron,comment,status)
!    if(status.ne.0) then
!       print*,'Problems in the subrotine gettimeinfo; check enoise'
!       stop
!    end if
!
!    call ftgkye(iunit,'EGAIN',gain,comment,status)
!    if(status.ne.0) then
!       print*,'Problems in the subrotine gettimeinfo; check egain'
!       stop
!    end if
!
!    call ftgkye(iunit,'NCOMBINE',ncombine,comment,status)
!    if(status.ne.0) then
!       ncombine=1
!    endif
!
!    gain=gain*ncombine
!    ron=ron*sqrt(ncombine)
!
!    call getdate(datestring,date1)
!    call getut(utstring,ut1)
!
!    return
!    
!  end subroutine gettimeinfo
!
!-----------------------------------------------------------------------
!
!  subroutine getdate(datestring,date1)
!
!    character(*) :: datestring
!    real :: date1(3)
!    integer :: jseparator(2)
!    
!    ndatestring = len_trim(datestring)
!    j = 1
!    do i=1,ndatestring
!       if((llt(datestring(i),'0')).or.(lgt(datestring(i),'9'))) then
!          jseparator(j) = i
!          j = j + 1
!       end if
!    end do
!   
!    if (jseparator(1)-1.eq.4) then
!      read(datestring(1:jseparator(1)-1),*) iyear
!      read(datestring(jseparator(1)+1:jseparator(2)-1),*) imonth
!      read(datestring(jseparator(2)+1:ndatestring),*) iday
!    else
!      read(datestring(1:jseparator(1)-1),*) iyear
!      read(datestring(jseparator(1)+1:jseparator(2)-1),*) imonth
!      read(datestring(jseparator(2)+1:ndatestring),*) iday
!    endif
!
!    date1(1) = float(iday)
!    date1(2) = float(imonth)
!    date1(3) = float(iyear)
!    
!    return
!    
!  end subroutine getdate
!
!-----------------------------------------------------------------------
!
!  subroutine getut(utstring,ut1)
!
!    character(*) :: utstring
!    real :: ut1(3)
!    integer :: jseparator(2)
!
!    nutstring = len_trim(utstring)
!    j = 1
!    do i=1,nutstring
!       if((llt(utstring(i),'0')).or.(lgt(utstring(i),'9'))) then
!          jseparator(j) = i
!          j = j + 1
!       end if
!    end do
!
!    read(utstring(1:jseparator(1)-1),*) ihour
!    read(utstring(jseparator(1)+1:jseparator(2)-1),*) iminute
!    read(utstring(jseparator(2)+1:nutstring),*) second
!    
!    ut1(1) = float(ihour)
!    ut1(2) = float(iminute)
!    ut1(3) = second
!
!    return
!
!  end subroutine getut
!
!-----------------------------------------------------------------------

subroutine getdatai2(iunit,array,n1,n2)
  
  integer(2) :: array(n1,n2)
  integer :: group,firstpix,status
  logical :: anynull
  
  group = 1
  nullval = -999
  firstpix = 1
  nelements = n1*n2
  status= 0
  call ftgpvi(iunit,group,firstpix,nelements,nullval,array,anynull,status)
  
end subroutine getdatai2

!-----------------------------------------------------------------------

subroutine getdatai4(iunit,array,n1,n2)
  
  integer(4) :: array(n1,n2)
  integer :: group,firstpix,status
  logical :: anynull

  group = 1
  nullval = -999
  firstpix = 1
  nelements = n1*n2
  status= 0
  call ftgpvj(iunit,group,firstpix,nelements,nullval,array,anynull,status)

end subroutine getdatai4

!-----------------------------------------------------------------------

subroutine getdatar4(iunit,array,n1,n2)
  
  real(4) :: array(n1,n2)
  integer :: group,firstpix,status
  logical :: anynull
  
  group = 1
  nullval = -999
  firstpix = 1
  nelements = n1*n2
  status= 0
  call ftgpve(iunit,group,firstpix,nelements,nullval,array,anynull,status)
  
end subroutine getdatar4

!-----------------------------------------------------------------------

subroutine putdatai2(iunit,array,n1,n2)
  
  real(8) :: bscale,bzero
  integer(2) :: array(n1,n2)
  integer :: group,firstpix,status
  
  bscale = 1.0
  bzero = 0.0
  status = 0
  call ftpscl(iunit,bscale,bzero,status)
  
  nullval = -999
  call ftpnul(iunit,nullval,status)
  
  group = 1
  firstpix = 1
  nelements = n1*n2
  call ftppni(iunit,group,firstpix,nelements,array,nullval,status)
  
end subroutine putdatai2

!-----------------------------------------------------------------------

subroutine putdatai4(iunit,array,n1,n2)
  
  real(8) :: bscale,bzero
  integer(4) :: array(n1,n2)
  integer :: group,firstpix,status
  
  bscale = 1.0
  bzero = 0.0
  status = 0
  call ftpscl(iunit,bscale,bzero,status)
  
  nullval = -999
  call ftpnul(iunit,nullval,status)
  
  group = 1
  firstpix = 1
  nelements = n1*n2
  call ftppnj(iunit,group,firstpix,nelements,array,nullval,status)
  
  return
end subroutine putdatai4

!-----------------------------------------------------------------------

subroutine putdatar4(iunit,array,n1,n2)
  
  real(8) :: bscale,bzero
  real(4) :: array(n1,n2)
  integer :: group,firstpix,status
  
  bscale = 1.0
  bzero = 0.0
  status = 0
  call ftpscl(iunit,bscale,bzero,status)
  
  nullval = -999
  call ftpnul(iunit,nullval,status)
  
  group = 1
  firstpix = 1
  nelements = n1*n2
  call ftppne(iunit,group,firstpix,nelements,array,nullval,status)
  
end subroutine putdatar4

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
