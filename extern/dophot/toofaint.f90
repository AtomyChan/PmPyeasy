! Function toofaint(starpar,err,npstar,cmin,crit7,lverb)
!
! This function defines which stars are used for a 7 parameter fit,
! i.e., to be PSF stars. 
! Objects with a smaller S/N, i.e., fainter stars, are not used.
   
function toofaint(starpar,err,npstar,cmin,crit7,lverb)

  real :: starpar(npstar),err(npstar)
  logical :: toofaint

  toofaint = starpar(2).lt.cmin

  sig2 = starpar(2)**2/err(2)
  toofaint = toofaint.and.(sig2.lt.crit7**2)
  
  if(lverb.gt.30) write(13,*) 'sig2 = ', sig2 

end function toofaint
