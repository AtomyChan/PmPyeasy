! Function probgal(xy,z,ze,npt,npstar,a,onestar,par1,par2,lverb)
!
! This function provides an ex post facto delineation of stars from galaxies
! (see second to last paragraph of section 3 in DoPHOT paper (Schechter, 
! Mateo, and Saha 1993 PASP, 105, 1345))
  
function probgal(xy,z,ze,npt,npstar,a,onestar,par1,par2,lverb)

  real :: a(npstar),fa(npstar)
  real :: b(npstar,npstar),c(npstar,npstar),v(npstar),vv(npstar)
  integer :: indx(npstar)
  real :: xy(2,npt),z(npt),ze(npt)
  external onestar

  b = 0.0
  v = 0.0
  chi = 0
  do i = 1,npt
     call onestar(xy(1,i),a,npstar,par1,par2,fa,ff)
     f = ff-z(i)
     if(ze(i).eq.0.0) then
        chi = 0.0
        exit
     end if
     chi = chi+f**2/ze(i)
     do j = 1,npstar
        vv(j) = vv(j)+f*fa(j)/ze(i)
        do k = 1,j
           b(j,k) = b(j,k)+fa(k)*fa(j)/ze(i)
           b(k,j) = b(j,k)
        end do
     end do
  end do
  if(chi.eq.0.0) then
     probgal = 0.0
     return
  end if

  v=vv
  c=b
  call ludcmp(c,npstar,npstar,indx,d)
  call lubksb(c,npstar,npstar,indx,v)
  do i = 1, npstar
     dchi = dchi+vv(i)*(-v(i))
     do j = 1, npstar
        dchi = dchi+v(i)*b(j,i)*v(j)/2
     end do
  end do

  h5 = a(5)-v(5)
  h6 = a(6)-v(6)
  h7 = a(7)-v(7)
  aarea = ELLIPAREA (a(5),a(6),a(7))
  harea = ELLIPAREA (h5,h6,h7)
  if (dchi.gt.0) then
     if(lverb.gt.20) then
        write(13,*) ' object at ',a(3),a(4),' has -ve delta chi**2 !' 
        write(13,*) ' dchi = ', dchi
     end if
     probgal = 0.0
  else if (aarea.le.0) then
     probgal = 0.0
  else
     chiper = chi/npt                        
     probgal = sign(dchi/chiper,harea-aarea)
  end if

end function probgal










