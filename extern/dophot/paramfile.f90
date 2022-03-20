! Subroutine param_in(nparam, > header,nlines)
! This subroutine asks the user for the parameter files.  If the modified file
! does not exist, the default file is transferred directly to the output file.
! A valid default parameter file (paramdefault) is supplied in 
! the DoPHOT directory.
     
subroutine param_in(nparam,header,nlines)

  integer,parameter :: iunitdef=40,iunitmod=41
  character(80) :: modheader(nparam),defheader(nparam)
  character(80) :: prompt,keyword,item,comment,header(nparam)
  character(40) :: fimod
 
!  Ask for the modification parameter file name.
!     or take it from the first command line parameter
  if (iargc().ge.1) then         ! function arg is not standard f90
     call getarg(1,prompt)       ! subroutine getarg is not standard f90
     nprompt = len_trim(prompt)     
     call opena(iunitmod,prompt(1:nprompt),0,ierr)
  else
     call query('Enter modification parameters file name. If none, enter NONE')
     read(*,*) fimod
     call opena(iunitmod,fimod,0,ierr)
  endif

  if (ierr.ne.0) then
     write(*,*) 'Error opening parameter modification file.'
     write(*,*) 'Using only default parameter file.'
     nlmod=1
  else
!    Read the modified parameter (.tuneup) file
     do i=1,nparam
        read(iunitmod,'(a)',err=200,end=200) modheader(i)
     end do
200  nlmod =i-1
  endif

!  Find the name of the default parameter file in the modified parameter file
!  If not found, ask for it
  call readitem(nlmod,modheader,'PARAMS_DEFAULT',nr,itype,item,ll)
  if(nr.gt.0) then
     ierr = 0
     call opena(iunitdef,item,0,ierr)
     if(ierr.eq.1) then
        print *,'Default parameters file not found!'
        nr = 0
     end if     
  end if
  if(nr.eq.0) then
     prompt = 'Enter default parameters file name'
     call opens(iunitdef,prompt,0)
  end if

!  Loop over the default parameter file and look for those keywords in the
!  modified parameter file.
  ii=0
  do i=1,nparam
     read(iunitdef,'(a)',err=210,end=210) defheader(i)
     call readkey(defheader(i),keyword,lenkey)   ! Get keywords from default
     if(lenkey.eq.0) cycle   ! No key=>full-line comment; not written output 
     ii=ii+1
     header(ii)=defheader(i)
!    See if you can find a value for that keyword in modified parameter file.
     call readitem(nlmod,modheader,keyword,nr,itype,item,ll)
!    Did you find it? If not, keep the default header line.
!    If yes, read the full line and write it out to the appropriate output file
     if(nr.gt.0.or.itype.eq.2) then
        call readcomment(modheader(ll),comment,ncomm)
        if(ncomm.eq.0) call readcomment(defheader(i),comment,ncomm)
        call writeitemc(modheader(ll),keyword,item,iunitdump2,comment,ncomm)
        header(ii) = modheader(ll)
     end if
  end do
210 nldef = i - 1
  nlines=ii

  close(iunitmod)
  close(iunitdef)

  return
  
end subroutine param_in

!---------------------------------------------------------------------------

! Subroutine param_out(header,nlines)
! Write current parameters in a file.

subroutine param_out(header,nlines)

  integer,parameter :: iunitparout=42
  character(80) :: header(nlines),prompt,item

!   Checks for an output parameter file. 
!   If none is found, asks for it.
  call readitem(nlines,header,'PARAMS_OUT',nr,itype,item,ll)
  if(nr.gt.0) then     
     call opena(iunitparout,item,2,ierr)
     if(ierr.eq.1) then
        print *,'Error opening output parameters file.'
        nr = 0
     end if
  end if
  if(nr.eq.0) then
     prompt = 'Enter output parameters file name'
     call opens(iunit,prompt,2)
  end if

!   Writes values to the output parameter file
  do i=1,nlines
     write(iunitparout,'(a)') header(i)
  end do

  close(iunitparout)
  
  return

end subroutine param_out
