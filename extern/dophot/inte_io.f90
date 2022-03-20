! This file contains the subroutines to read or write the INTERNAL dophot format
!
! INTE_IN
! INTE_OUT
!
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!
! Subroutine inte_in(i,itype,starpar,npstar,appar,npap,probg,ios,iunit)
! This subroutine reads a file with INTERNAL format.

subroutine inte_in(i,itype,starpar,npstar,appar,npap,probg,ios,iunit)
     
  real :: starpar(npstar),appar(npap)

  read(iunit,*,iostat=ios) i,itype,(starpar(k1),k1=1,npstar),&
       (appar(k2),k2=1,npap),probg

end subroutine inte_in

!---------------------------------------------------------------------------
!
! Subroutine inte_out(i,itype,starpar,npstar,appar,npap,probg,iunit)
! This subroutine writes a file with INTERNAL format.

subroutine inte_out(i,itype,starpar,npstar,appar,npap,probg,iunit)
  
  real ::  starpar(npstar),appar(npap)
  character(2) :: fmtvar
  character(120) :: fmt

  write(fmtvar,'(i2)') npstar-4+npap+1
  fmt = "(i6,i3,2e11.3,2f11.2,"//fmtvar//"e11.3)"
  write(iunit,fmt) i,itype,(starpar(k),k=1,npstar),(appar(kk),kk=1,npap),probg

       
end subroutine inte_out
