! Subroutine chisq(functn,par1,par2,xy,z,zerr,n,A,npar,m,acc,alim,it,&
!     rmagic,lverb, > ERR,CHISQU)
! This subroutine looks for better parameters of the fit function, 
! using nonlinear least-squares fitting.
! The nonlinear least-squares fitting is done by minimizing chi square, 
! using the Levenberg-Marquardt method.
! The input of these subroutine are 
!  FUNCTN = the function to be fit; must be external in calling program.
!  PAR1,PAR2 = parameters of the function FUNCTN
!  XY = array with independent variables
!  Z = dependent variable
!  ZERR = variance of dependent variable
!  N = number of values to be fit
!  A = parameters of fit function 
!  NPAR = number of parameters of the fit function 
!  M = order of fit
!  ACC = accuracy array; defines convergence
!  ALIM = limits array; defines valid limits for parameters
!  IT = maximum number of iterations allowed
! As an intermediate step FA, the partial derivatives of 
! the fit function with respect to the parametersare generated.
! Also as an intermediate step the covariance matrix C is generated, 
! whose diagonal terms are going to give us the vector ERR with the 
! variance of the parameters. 
! The output are A, ERR, and CHISQ.
! To better understand Levenberg-Marquardt method and this subroutine, 
! check section 15.5 in Numerical Recipes in Fortran, and section 4.11 in
! DoPHOT paper (Schechter, Mateo, and Saha 1993 PASP, 105, 1345)  

subroutine chisq(functn,par1,par2,xy,z,zerr,n,a,npar,m,acc,alim,it,&
     rmagic,lverb,err,chisqu)

  real :: a(npar),err(npar),acc(npar),alim(npar)
  real :: xy(2,n),z(n),zerr(n)
  real :: fa(npar),c(m,m),b(m,m),vv(m),v(m)
  integer :: indx(m)
  logical :: conv,marq,limit
  external functn

  conv = .false.
  limit = .false.
  do j=1,m
!  Check if parameters exceeds limits at the start.  
!  Check only absolute parameters (alim(J) < 0).
     if(alim(j).lt.0.0) then
        limit = abs(a(j)).gt.abs(alim(j))
        if(limit) then
           if(lverb.gt.20) write(13,*) & 
                'SELF-DECEPTION HAS OCCURED: INITIAL LIMITS.'
        end if
     end if
  end do
  
  if(lverb.gt.30) print '(i3,i2,7e10.2)', 0,0,(a(kk),kk=1,m)
  
  i = 1
  ifact = 0
  do while ((i.le.it).and.(.not.conv).and.(.not.limit))
     chi = 0.0
     b = 0.0
     vv = 0.0
     
     do j=1,n
        call functn(xy(1,j),a,npar,par1,par2,fa,ff)
        f = ff-z(j)
        if(zerr(j).ne.0.0) chi = chi+f**2/zerr(j)
        if(zerr(j).gt.0.0) then 
           do kk=1,m
              vv(kk) = vv(kk)+f*fa(kk)/zerr(j)
              do l=1,kk
                 b(kk,l) = b(kk,l)+fa(l)*fa(kk)/zerr(j)      
                 b(l,kk) = b(kk,l)
              end do
           end do
        end if
     end do
     chiold = chi

     k = 1
     fact = 0.0
     marq = .false.
     do while ((k.le.10).and.(.not.marq)) 
        conv = (k.eq.1)
        v = vv
        c = b
        do j=1,m
           c(j,j) = (1+fact)*b(j,j)
        end do

        call ludcmp(c,m,m,indx,d)
!       If d = 0, the matrix was singular; no convergence.
        if(d.eq.0.) then
           chisqu = rmagic
           return
        end if
        call lubksb(c,m,m,indx,v)

        do j=1,m
           a(j) = a(j)-v(j)
           if(a(j).eq.0.0) cycle
!     Check if change in parameters exceeds limits.  If alim(j) > 0, then
!     consider fractional changes.  If alim(j) < 0, consider absolute
!     changes.  If alim(j) = 0, ignore this test.
           if(alim(j).gt.0.0) then
              limit = limit.or.(abs(v(j)/a(j)).gt.alim(j))
              if(limit) then
                 if(lverb.gt.20) write(13,*) &
                      'SELF-DECEPTION HAS OCCURED: FRAC LIMITS.'
              end if
           else if(alim(j).lt.0.0) then
              limit = limit.or.(abs(a(j)).gt.abs(alim(j)))
              if(limit) then
                 if(lverb.gt.20) write(13,*) &
                      'SELF-DECEPTION HAS OCCURED: ABS LIMITS.'
              end if
           end if
           
           if(acc(j).ge.0) then
              conv = conv.and.(abs(v(j)/a(j)).le.acc(j))
           else
              conv = conv.and.(abs(v(j)).le.abs(acc(j)))
           end if
        end do

        if(lverb.gt.30) print'(i3,i2,7e10.2)',i,k,(a(kk),kk=1,m)
        
        if(conv) then            
           marq = .true.
        else if(.not.limit) then
           chi = 0.0
           do j = 1,n 
              call functn(xy(1,j),a,npar,par1,par2,fa,ff)
              f = ff-z(j)
              if (zerr(j).ne.0.0) chi = chi+f**2/zerr(j)
           end do
           if(chi.le.chiold) then
              marq = .true.
           else
              ifact = ifact+1
              fact = 2.0**ifact
              do j=1,m
                 a(j) = a(j)+v(j)     !So a() remain the same as before line 90
              end do
           end if
        end if
        k = k + 1
     end do
     i = i + 1
  end do
    
! Finished with Levengerg-Marquardt after obtaining new a()
! Now looking for values of c(,)
  if(.not.limit) then
     b=0.0
     do j=1,m                  
        b(j,j) = 1
        call lubksb(c,m,m,indx,b(1,j))
     end do
  end if
  
  if(conv.and.(.not.limit)) then
     perdeg = chi/max0(n-m,1)
     
     if(perdeg.le.0.0) then
        chisqu = rmagic
        return
     end if
     
     perdeg = sqrt(perdeg)
     
     do i=1,m
        if(b(i,i).gt.0) then
           savve = sqrt(b(i,i))
        else
           if(lverb.gt.20) write(13,*) & 
                'TROUBLE: NEGATIVE AUTOVARIANCE FOR I, B(I,I) = ',i,b(i,i)
           chisqu = rmagic
           return
        end if
        
        do j=1,m
           if(savve.ne.0.0) then
              c(i,j) = b(i,j)/savve
              c(j,i) = b(j,i)/savve
           end if
        end do
        c(i,i) = savve*perdeg
        err(i) = c(i,i)**2
     end do
     chisqu = chiold
  else
     chisqu = rmagic
  end if

end subroutine chisq
