! Subroutine shadow_out(i,itype,starpar,npstar,iunit)
! This subroutine writes the output for a shadow file,
! in a short version of the INTERNAL format type (see INTE_IO.f90)

subroutine shadow_out(i,itype,starpar,npstar,iunit)

  real :: starpar(npstar)
  character(2) :: fmtvar
  character(120) :: fmt

  write(fmtvar,'(i2)') npstar-4
  fmt = "(i6,i3,2e11.3,2f11.2,"//fmtvar//"e11.3)"
  write(iunit,fmt) i,itype,(starpar(k),k=1,npstar)

end subroutine shadow_out

