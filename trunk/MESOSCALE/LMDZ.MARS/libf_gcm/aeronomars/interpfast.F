      subroutine interpfast(escout,p,nlayer,escin,pin,nl,limdown,limup)
C
C subroutine to perform linear interpolation in pressure from 1D profile 
C escin(nl) sampled on pressure grid pin(nl) to profile
C escout(nlayer) on pressure grid p(nlayer).
C
      real escout(nlayer),p(nlayer)
      real escin(nl),pin(nl),wm,wp
      real limup,limdown
      integer nl,nlayer,n1,n,nm,np
      nm=1
      do 5 n1=1,nlayer
         if(p(n1) .gt. limup .or. p(n1) .lt. limdown) then
            escout(n1) = 0.0
         else
            do n = nm,nl-1
               if (p(n1).ge.pin(n).and.p(n1).le.pin(n+1)) then
                  nm=n
                  np=n+1
                  wm=abs(pin(nm)-p(n1))/(pin(np)-pin(nm))
                  wp=1.0 - wm
                  goto 33
               endif
            enddo
 33         escout(n1) = escin(np)*wm + escin(nm)*wp
         endif
    5 continue
      return
      end
