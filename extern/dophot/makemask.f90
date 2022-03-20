!subroutine makemask(x,y,flags,nff,enuffvpsf,a5,a6,a7,npsffit,ava,npstar,&
!     starfunc,par1,par2,imaskrad,starmask)
! This subroutine creates a mask representing the current model of a stellar PSF
! to test (later) if an object exists in a certain region of the sky.

subroutine makemask(x,y,flags,nff,enuffvpsf,a5,a6,a7,npsffit,ava,npstar,&
     starfunc,par1,par2,maskrad,starmask)

  real :: starmask(-maskrad:maskrad,-maskrad:maskrad)
  real :: xy(2)
  real :: a(npstar),ava(npstar),dummy(npstar)
  real :: a5(npsffit),a6(npsffit),a7(npsffit)
  character(*) :: flags(nff)
  logical :: enuffvpsf
  external starfunc

  starmask = 0.

  call stpar567(x,y,flags,nff,enuffvpsf,a5,a6,a7,npsffit,ava,npstar,a)
  a(1) = 0.
  a(2) = 1.                 
  a(3) = x               
  a(4) = y
       
  do j = -maskrad, maskrad
     xy(2) = a(4)+j
     do i = -maskrad, maskrad
        xy(1) = a(3)+i
        call starfunc(xy,a,npstar,par1,par2,dummy,starmask(i,j))
     enddo
  enddo

end subroutine makemask
