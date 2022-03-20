! This file contains the subroutines to choose between 
! the different subroutines that reads (write) the formatted input (output).
!
! CHOOSE_IN
! CHOOSE_OUT
!
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!
! Subroutine choose_in(iunit,flag,i,itype,starpar,npstar,appar,npap,&
!     probg,ios)
! This subroutine chooses between the different subroutines 
! that read the formatted input.

subroutine choose_in(iunit,flag,i,itype,starpar,npstar,appar,npap,&
     amag,probg,ios)

  real :: starpar(npstar),appar(npap)
  character(*) :: flag

  starpar = 0.
  appar = 0.
  probg = 0.
  if(flag(1:5).eq.'INTER') then
     call inte_in(i,itype,starpar,npstar,appar,npap,probg,ios,iunit)
     amag=SHAPEMAG(starpar,npstar)
  else if(flag(1:5).eq.'COMPL') then
     call comp_in(i,itype,starpar,npstar,appar,npap,probg,ios,iunit)
     amag=SHAPEMAG(starpar,npstar)
  else if(flag(1:5).eq.'BINAR') then
     call bina_in(i,itype,starpar,npstar,appar,npap,ios,iunit)
     amag=SHAPEMAG(starpar,npstar)
  else if(flag(1:5).eq.'DAOPH') then
     call daoph_in(i,itype,starpar,npstar,appar,npap,amag,ios,iunit)
  end if

end subroutine choose_in

!---------------------------------------------------------------------------
!
! Subroutine choose_out(iunit,flag,i,itype,starpar,npstar,appar,npap,&
!     probg,applyap)
! This subroutine chooses between the different subroutines 
! that write the formatted output.

subroutine choose_out(iunit,flag,i,itype,starpar,npstar,appar,npap,&
     probg,applyap)

  real :: starpar(npstar),appar(npap)
  character(*) :: flag
  logical :: applyap

  if(flag(1:5).eq.'INTER') then
     call inte_out(i,itype,starpar,npstar,appar,npap,probg,iunit)
  else if(flag(1:5).eq.'COMPL') then
     call comp_out(i,itype,starpar,npstar,appar,npap,probg,iunit,applyap)
  else if(flag(1:5).eq.'BINAR') then
     call bina_out(i,itype,starpar,npstar,appar,npap,iunit,applyap)
  else if(flag(1:5).eq.'INCOM') then
     call incomp_out(i,itype,starpar,npstar,appar,npap,probg,iunit,applyap)
  else if(flag(1:5).eq.'DAOPH') then
     call daoph_out(i,itype,starpar,npstar,appar,npap,iunit,applyap)
  end if

end subroutine choose_out

