      SUBROUTINE moyzon(iim,jjp1,nl,indefini,x,xbar)
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

      INTEGER iim,jjp1,nl
      real indefini
      REAL x(iim+1,jjp1,nl)
      REAL xbar(jjp1,nl)

c   Local:
c   ------

      iNTEGER i,j,l , n

c------------------------------------------------------------------------
      do l=1,nl
        do j=1,jjp1
           xbar(j,l)=0.
           n = 0
           do i=1,iim
              if (x(i,j,l).ne.indefini) then 
                 xbar(j,l) = xbar(j,l) + x(i,j,l)
                 n = n+1
              end if
           end do
           if (n.ne.0) then
              xbar(j,l) = xbar(j,l)/float(n)
           else
              xbar(j,l) = indefini
           end if  
        end do
      end do 

      RETURN
      END
