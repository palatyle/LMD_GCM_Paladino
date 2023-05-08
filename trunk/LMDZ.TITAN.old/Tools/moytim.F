      SUBROUTINE moytim(iim,jjp1,nl,nt,indefini,x,xmean)
c=======================================================================
c
c
c   Subject:
c   ------
c   Calcul de la moyenne zonale de la variable au point scalaire x
c
c=======================================================================
      IMPLICIT NONE
c-----------------------------------------------------------------------
c   Declararations:
c   ---------------

c   Arguments:
c   ----------

      INTEGER iim,jjp1,nl,nt
      real indefini
      REAL x(iim,jjp1,nl,nt)
      REAL xmean(iim,jjp1,nl)

c   Local:
c   ------

      iNTEGER i,j,l,t , n

c------------------------------------------------------------------------
      do l=1,nl
       do j=1,jjp1
        do i=1,iim
           xmean(i,j,l)=0.
           n = 0
           do t=1,nt
              if (x(i,j,l,t).ne.indefini) then 
                 xmean(i,j,l) = xmean(i,j,l) + x(i,j,l,t)
                 n = n+1
              end if
           end do
           if (n.ne.0) then
              xmean(i,j,l) = xmean(i,j,l)/float(n)
           else
              xmean(i,j,l) = indefini
           end if  
        end do
       end do
      end do 

      RETURN
      END
