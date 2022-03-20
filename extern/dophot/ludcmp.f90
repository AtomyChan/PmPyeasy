! Subroutine ludcmp(a,n,np,indx,d)
! This subroutine does an LU decomposition, i.e., a matrix decomposition 
! which writes a matrix as the product of a lower triangular matrix and 
! an upper triangular matrix.
! The method used for the LU decomposition is the Crout's method with partial
! pivoting.
! This subroutine is used in combination with subroutine LUBKSB.
! Note that the original matrix A is substituted by the elements of 
! the L and U triangular matrices.
! The vector INDX records the row permutation effected by the partialpivoting.
! The number d has three possible values: -1 if the number of row interchanges
! was odd; +1 if it was even; and 0 if the matrix was singular.  
! To better understand this subroutine and the LU decomposition, 
! check section 2.3 in Numerical Recipes in Fortran.

subroutine ludcmp(a,n,np,indx,d)

  real :: a(np,np),vv(n)
  integer :: indx(n)

  d=1.
  do i=1,n
     aamax=0.
     do j=1,n
        if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
     end do
     if (aamax.eq.0.) then
        d = 0
        return
     end if
     vv(i)=1./aamax
  end do
  do j=1,n
     do i=1,j-1
        sum=a(i,j)
        do k=1,i-1
           sum=sum-a(i,k)*a(k,j)
        end do
        a(i,j)=sum
     end do
     aamax=0.
     do i=j,n
        sum=a(i,j)
        do k=1,j-1
           sum=sum-a(i,k)*a(k,j)
        end do
        a(i,j)=sum
        dum=vv(i)*abs(sum)
        if (dum.ge.aamax) then
           imax=i
           aamax=dum
        end if
     end do
     if (j.ne.imax)then
        do k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
        end do
        d=-d
        vv(imax)=vv(j)
     endif
     indx(j)=imax
     if (a(j,j).eq.0) then
        d = 0
        return
     end if
     if(j.ne.n)then
        dum=1./a(j,j)
        do i=j+1,n
           a(i,j)=a(i,j)*dum
        end do
     endif
  end do
     
end subroutine ludcmp





