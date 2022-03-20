! Program io_dophot
! Works with the input and output files of dophot
!
! Includes these subroutines:
! QUERY
! QUERYA
! INFO
!
! OPENS
! OPENSB
! OPENA
! OPENAB
!
! READKEY
! READITEM
! READCOMMENT
! WRITEITEMR
! WRITEITEMI
! WRITEITEMC
! WRITEEND
!
! UPPER
! YES_OR_NO
!
!
!---------------------------------------------------------------------------

! query (prompt)
! This subroutine asks for input.

subroutine query(prompt)
    
  character(*) :: prompt
  
  nprompt = len_trim(prompt)        
  write(*,'(1x,a,$)') prompt(1:nprompt)//' '   ! format $ is not standard f90
  
end subroutine query

!---------------------------------------------------------------------------

! querya (prompt,value,variable)
! This subroutine asks for input.

  subroutine querya(prompt,value,variable)
  
    character(*) :: prompt
    character(80) :: line

    call query(prompt)   
    read(*,'(a)') line 
    nline = len_trim(line)
    if(nline.le.0) then
       variable = value
    else
    read(line,*) variable
    end if

  end subroutine querya

!---------------------------------------------------------------------------

! info(prompt)
! This supplies information.

  subroutine info(prompt)
        
    character(*) prompt

    nprompt = len_trim(prompt)
    write(6,'(a)') prompt(1:nprompt)

  end subroutine info

!---------------------------------------------------------------------------

! opens (iunit,prompt,iopen)
! This subroutine opens or closes a sequential file.  
! It also handles basic errors.

subroutine opens(iunit,prompt,iopen)
  
  character(*) :: prompt
  character(80) :: name

10 call query(prompt)

  read(*,'(a)',err=210) name
  nname=len_trim(name)
  
  if(iopen.eq.0) open(unit=iunit,status='old',file=name(1:nname),err=211)
  if(iopen.eq.1) open(unit=iunit,status='new',file=name(1:nname),err=211)
  if(iopen.eq.2) open(unit=iunit,status='unknown',file=name(1:nname),err=211)
  return
  
210 print *,'Error in file name. Try again.'
  go to 10
  
211 print *,'Error opening file. Try again.'
  go to 10
  
end subroutine opens

!---------------------------------------------------------------------------

! opensb(iunit,prompt,iopen,irecl)
! This subroutine opens or closes a binary file.  
! It also handles basic errors.

subroutine opensb(iunit,prompt,iopen,irecl)
  
  character(*) :: prompt
  character(80) :: name
  
  call query(prompt)
  
  read(*,'(a)',err=210) name
  nname=len_trim(name)
  
  if(iopen.eq.0) open(unit=iunit,status='old',access='direct',&
       recl=irecl,file=name(1:nname),err=211)
  if(iopen.eq.1) open(unit=iunit,status='new',access='direct',&
       recl=irecl,file=name(1:nname),err=211)
  if(iopen.eq.2) open(unit=iunit,status='unknown',access='direct',&
       recl=irecl,file=name(1:nname),err=211)
  
  return
  
210 print *,'Error in file name. Exiting.'
  call exit(1)     !subroutine exit is an extension to f90
  
211 print *, 'Error opening file. Exiting.'
  call exit(1)     !subroutine exit is an extension to f90
  
end subroutine opensb

!---------------------------------------------------------------------------

! opena(iunit,name,iopen,ierr)
! This subroutine opens a sequential file given a file name.
! It returns an error code, 0 if ok, 1 if there is a problem.

subroutine opena(iunit,name,iopen,ierr)
  
  character(*) :: name
  
  ierr = 0
  
  if(iopen.eq.0) open(unit=iunit,status='old',file=name,err=210)
  if(iopen.eq.1) open(unit=iunit,status='new',file=name,err=210)
  if(iopen.eq.2) open(unit=iunit,status='unknown',file=name,err=210)
  if(iopen.eq.3) open(unit=iunit,status='unknown',position='append',&
       file=name,err=210)
  
  return
  
210 ierr = 1
  return
  
end subroutine opena

!---------------------------------------------------------------------------

! openab(iunit,name,iopen,ierr,irecl)
! This subroutine opens a binary file given a file name.
! It returns an error code, 0 if ok, 1 if there is a problem.

subroutine openab(iunit,name,iopen,ierr,irecl)
  
  character(*) :: name
  
  ierr = 0
  
  if(iopen.eq.0) open(unit=iunit,status='old',access='direct',&
       recl=irecl,file=name,err=210)
  if(iopen.eq.1) open(unit=iunit,status='new',access='direct',&
       recl=irecl,file=name,err=210)
  if(iopen.eq.2) open(unit=iunit,status='unknown',access='direct',&
       recl=irecl,file=name,err=210)
  
  return
  
210 ierr = 1
  return
  
end subroutine openab

!---------------------------------------------------------------------------

!  readkey(header, > keyword,lenkey)
!  Reads a line and extracts the keyword from a single header line.  
!  LENKEY is the length of the keyword.

subroutine readkey(header,keyword,lenkey)
  
  character(*) :: header,keyword
  
  lenmax = len_trim(header)
  do j=1,lenmax
     if(header(j:j).eq.'=') then
        lenkey = j - 1
        exit
     end if
  end do
  
  if(lenkey.eq.0) return
  
  do j=lenkey-1,1,-1
     if(header(j:j).ne.' ') then
        lenkey = j
        exit
     end if
  end do
  
  keyword = header(1:lenkey)
  call upper(keyword)
  
end subroutine readkey

!---------------------------------------------------------------------------
 
! readitem(n,header,keyword, > nr,itype,item,mline)
! Reads the item with key keyword from a fits header or a fits type header.   
! The arguments are: 
!   n        = number of header lines
!   header   = character variable to hold the header data.
!   keyword  = parameter name (always capitalized)
!   nr       = number of character or digits returned
!   itype    = type of parameter (0=real;1=integer;2=character)
!   item     = value of the parameter
!   mline    = number of the line in the header where the keyword is 
! ITYPE is 0 for a real variable, 1 for an integer (assumed from 
! leading letter of keyword), and 2 if a character variable (assumed 
! because the parameter is inside quotation marks).

subroutine readitem(n,header,keyword,nr,itype,item,mline)
  
  character(*) :: header(n),keyword,item
  character :: first
  
! Locate the item.
  
  itype = 0
  mline = 0
  nr = 0
  do i = 1,n
     lenmx = len_trim(header(i))
     do j=1,lenmx
        if(header(i)(j:j).eq.'=') exit
     end do
     lenhead1 = len_trim(header(i)(1:j - 1))
     lenhead2 = j

!  Find a match to the keyword.

    call upper(header(i)(1:lenhead1))
     if(header(i)(1:lenhead1).eq.keyword) then
        mline = i
        exit
     end if
!    No match? Then return
     if((i.eq.n).and.(header(i)(1:lenhead1).ne.keyword)) then
        nr=0           
        return
     end if
  end do
  
! Locate the beginnings and ends of the item value

  i1 = lenhead2 + 1
  

!   First, locate the beginnings of the string
  is = i1
  do i=is,80
!    Strip leading blanks.
     if(header(mline)(i:i).ne.' ') then

!    Determine if this is a character variable.   
        if(header(mline)(i:i).eq.'''') then
           itype = 2
           i1 = i + 1
           exit
        end if
        
        i1 = i
        exit
     end if
  end do

!  Now, find the end of the string.
  lenmax = len_trim(header(mline))
  i2 = i1
  do i=i1+1,lenmax
!    Take care if this is a character.
     if((itype.eq.2.and.header(mline)(i:i).eq.'''').or.&
          (header(mline)(i:i).eq.' ')) then
        i2 = i - 1
        exit
     else if(i.eq.lenmax) then
        i2 = i
        exit
     end if
     
  end do
 
!  Strip trailing blanks
  nr = len_trim(header(mline)(i1:i2))
    
!  Read the item value
  item = header(mline)(i1:i1+nr-1)

!  Determine if real or integer if not character.
  if(itype.ne.2) then
     itype = 0
     first = keyword(1:1)
     if(first.eq.'I'.or.first.eq.'i') itype = 1
     if(first.eq.'J'.or.first.eq.'j') itype = 1
     if(first.eq.'K'.or.first.eq.'k') itype = 1
     if(first.eq.'L'.or.first.eq.'l') itype = 1
     if(first.eq.'M'.or.first.eq.'m') itype = 1
     if(first.eq.'N'.or.first.eq.'n') itype = 1
  end if

end subroutine readitem

!---------------------------------------------------------------------------

! readcomment(header,> comment,ncomm)
! Extracts the comment, if any exists, from a header line.

subroutine readcomment(header,comment,ncomm)
  
  character(*) :: header,comment
  
  lenmax = len_trim(header)
  
  if(lenmax.le.2) then
     ncomm = 0         ! If lenmax.le.2 then there is no comment here.
     return
  end if
  
  do i=1,lenmax
     if(header(i:i).eq.'=') then
        len1 = i
        exit
     end if
  enddo
    
! Now, find the end of the parameter.  If len1 = 1, then the comment
! starts after the first blank; otherwise it starts after the second one.

  do i=len1+1,lenmax
     if(header(i:i).ne.' ') then
        nz1 = i
        exit
     end if
  end do
  
  if(len1.eq.1) then
     ncomm = lenmax
     comment(1:ncomm) = header(1:lenmax)
     return
  end if
  
  if(nz1.eq.lenmax) then
     ncomm = 0
     return
  end if
  
  do i=nz1+2,lenmax
     if(header(i:i).eq.' ') then
        nc1 = i
        ncomm = lenmax - nc1 + 1
        comment(1:ncomm) = header(nc1:lenmax)
        return
     end if
  end do
  ncomm = 0
    
end subroutine readcomment

!---------------------------------------------------------------------------

! writeitemr(header,keyword,item,formatt,nformatt,iunit,comment,ncomm)
! Writes a real item out to a header line.
! The arguments are: 
!   header   = character variable to hold the header data.
!   keyword  = parameter name (always capitalized)
!   item     = value of the parameter
!   formatt  = format for the data
!   nformatt = number of digits in the format
!   iunit    = number of the output unit
!   comment  = comment at the end of the line 
!   ncomm    = number of characters in the comment 

subroutine writeitemr(header,keyword,item,formatt,nformatt,iunit,comment,ncomm)

  character(*) :: header,keyword,formatt,comment
  real :: item
  character(80) :: form
  
  nf = len_trim(formatt) + 2
  form(1:nf) = '('//formatt(1:nf-2)//')'
  nkeyword = len_trim(keyword)
  nlast = nkeyword+3
  header(1:nlast) = keyword(1:nkeyword)//' = '
  write(header(nlast+1:nlast+nformatt),form(1:nf)) item
  nheader = nlast + nformatt
  if (ncomm.gt.0) then
     nstart = nheader+1
     ncomm2 = min0(ncomm,80-nstart)
     nheader = nstart + ncomm2
     header(nstart:nheader) = comment(1:ncomm2)
  end if
  if (iunit.ge.1) write(iunit,'(a)') header(1:nheader)
  
end subroutine writeitemr

!---------------------------------------------------------------------------

! writeitemi(header,keyword,item,formatt,nformatt,iunit,comment,ncomm)
! Writes a real item out to a header line.
! The arguments are: 
!   header   = character variable to hold the header data.
!   keyword  = parameter name (always capitalized)
!   item     = value of the parameter
!   formatt  = format for the data
!   nformatt = number of digits in the format
!   iunit    = number of the output unit
!   comment  = comment at the end of the line 
!   ncomm    = number of characters in the comment 

subroutine writeitemi(header,keyword,item,formatt,nformatt,iunit,comment,ncomm)

  character(*) :: header,keyword,formatt,comment
  integer :: item
  character(80) :: form

  nf = len_trim(formatt) + 2
  form(1:nf) = '('//formatt(1:nf-2)//')'
  nkeyword = len_trim(keyword)
  nlast = nkeyword+3
  header(1:nlast) = keyword(1:nkeyword)//' = '
  write(header(nlast+1:nlast+nformatt),form(1:nf)) item
  nheader = nlast + nformatt
  if (ncomm.gt.0) then
     nstart = nheader+1
     ncomm2 = min0(ncomm,80-nstart)
     nheader = nstart + ncomm2
     header(nstart:nheader) = comment(1:ncomm2)
  end if
  if (iunit.ge.1) write(iunit,'(a)') header(1:nheader)

end subroutine writeitemi

!---------------------------------------------------------------------------

! writeitemc(header,keyword,item,formatt,nformatt,iunit,comment,ncomm)
! Writes a real item out to a header line.
! The arguments are: 
!   header   = character variable to hold the header data.
!   keyword  = parameter name (always capitalized)
!   item     = value of the parameter
!   formatt  = format for the data
!   nformatt = number of digits in the format
!   iunit    = number of the output unit (0 => not written in output unit)
!   comment  = comment at the end of the line 
!   ncomm    = number of characters in the comment 

subroutine writeitemc(header,keyword,item,iunit,comment,ncomm)

  character(*) :: header,keyword,item,comment
  
  nkeyword = len_trim(keyword)
  nitem = len_trim(item)
  nheader = nkeyword + 3 + nitem
  header(1:nheader) = keyword(1:nkeyword)//' = '//item(1:nitem)
  if (ncomm.gt.0) then
     nstart = nheader+1
     ncomm2 = min0(ncomm,80-nstart)
     nheader = nstart + ncomm2
     header(nstart:nheader) = comment(1:ncomm2)
  end if
  if (iunit.ge.1) write(iunit,'(a)') header(1:nheader)
  
end subroutine writeitemc

!---------------------------------------------------------------------------

! writeend(iunit)
! This subroutine finishes a parameter file
  
subroutine writeend(iunit)
  
  write(iunit,'(a)') 'END'
  
end subroutine writeend

!---------------------------------------------------------------------------

! upper(stringg)
! This subroutine transforms little letters into capital letters

subroutine upper(stringg)
    
  character(*) :: stringg
    
  nstring = len_trim(stringg)
  do i=1,nstring
     nhold = ichar(stringg(i:i))
     if(nhold.ge.97.and.nhold.le.122) nhold = nhold - 32
     stringg(i:i) = char(nhold)
  end do

end subroutine upper
  
!---------------------------------------------------------------------------

! yes_or_no(prompt, > answer)
! This subroutine handles yes or no questions, ensuring that
! y or n are the only allowable answers returned to the calling program.

subroutine yes_or_no(prompt,answer)
  
  character(*) :: prompt
  character(1) :: answer

  nprompt = len_trim(prompt)
  answer = 't'
  do while (answer.ne.'y'.and.answer.ne.'n')
     call query(prompt(1:nprompt))
     read (5,'(a)') answer
     if(answer.eq.'y'.or.answer.eq.'Y') answer = 'y'      
     if(answer.eq.'n'.or.answer.eq.'N') answer = 'n'      
  end do
      
end subroutine yes_or_no

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
