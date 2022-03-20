! Subroutine notifyend(iunit,thresh,nstar,nstot)
! This subroutine print some useful information of the loop we are in.

subroutine notifyend(iunit,thresh,nstar,nstot)

  write(iunit,100) thresh
100 format('  Ending loop at thresh = ',f9.1)
  write(iunit,101) nstar
101 format('  Number of new objects found = ',i7)
  write(iunit,102) nstot
102 format('  Total number of objects found so far = ',i7)
  write(iunit,*) ' '

end subroutine notifyend
