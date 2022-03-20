! Subroutine cosmic(big,noise,nfast,nslow,onestar,par1,par2,npstar,&
!     widobl,discrim,sncos,rmagic,lverb,>starpar,cosmicval)
! This subroutine decides if the found object is a cosmic ray

subroutine cosmic(big,noise,nfast,nslow,onestar,par1,par2,npstar,&
     widobl,discrim,sncos,rmagic,lverb,starpar,cosmicval)

  real :: big(nfast,nslow),noise(nfast,nslow)
  real :: starpar(npstar),aa(npstar),dummy(npstar),xy(2)
  logical :: cosmicval
  external onestar

  sncos2=sncos**2
  cosmicval = .false.
  sky=starpar(1)               
  rmaxob = -rmagic
  chistar = 0.
  chicos  = 0.
  npix = 0
  aa = starpar
  do ky = -1, 1
     do kx = -1, 1
        ix=int(aa(3))+kx
        iy=int(aa(4))+ky
        if ((ix.le.nfast).and.(ix.ge.1).and.(iy.le.nslow).and.(iy.ge.1)) then
           xy(1) = float(ix)
           xy(2) = float(iy)
           if (noise(ix,iy).lt.rmagic) then
              npix = npix+1
              call onestar(xy,aa,npstar,par1,par2,dummy,pred)
              if(noise(ix,iy).eq.0) then
                 if(lverb.gt.20) write(13,*)& 
                      'Noise = ZERO! ix, iy = ',ix,iy
                 cosmicval = .true.
                 return
              end if
              obs = big(ix,iy)
              temp = 1.0/noise(ix,iy)
              chistar = chistar+(obs-pred)**2*temp
              sn2 = (obs-sky)**2*temp
              chicos = chicos+sn2
              if ((obs.gt.rmaxob).and.(sn2.ge.sncos2)) then
                 imax = ix
                 jmax = iy
                 rmaxob = obs
                 tnoise = temp
              end if
           end if
        end if
     end do
  end do
  if ((npix.ge.7).and.(rmaxob.ne.-rmagic)) then        
     chicos = chicos-(rmaxob-sky)**2*tnoise
     cosmicval = chicos.lt.discrim*chistar          
     if(lverb.gt.30) then
        write(13,*) 'location: ', imax,jmax
        write(13,*) 'chi-star & chi-cosmic = ', chistar, chicos
     end if
  end if
  if (cosmicval) then
     if(lverb.gt.20) then
        write(13,*) 'cosmic ray intensity, ix, iy : ',big(imax,jmax),imax,jmax
     end if
     starpar(1) = sky
     starpar(2) = rmaxob
     starpar(3) = imax
     starpar(4) = jmax
     starpar(5) = widobl
     starpar(6) = -1
     starpar(7) = widobl
  end if
     
end subroutine cosmic
