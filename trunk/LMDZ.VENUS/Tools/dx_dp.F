      SUBROUTINE dx_dp(jjp1,llm,indefini,pniv,x,dxdp)
c=======================================================================
c
c
c   Subject:
c   ------
c   Calcul de la derivee /p d'1 moyenne zonale 
c   EXTRAPOLEE EN COORDONNEE DE PRESSION
c
c=======================================================================
      IMPLICIT NONE
c-----------------------------------------------------------------------
c   Declararations:
c   ---------------

c   Arguments:
c   ----------

      integer jjp1,llm
      real indefini
      REAL pniv(llm)
      REAL x(jjp1,llm)
      REAL dxdp(jjp1,llm)

c   Local:
c   ------

      INTEGER i,j,l , k

c-----------------------------------------------------------------------
      do j=1,jjp1
         if ((x(j,1).ne.indefini).and.(x(j,2).ne.indefini)) then  
           dxdp(j,1) = (x(j,2)-x(j,1))/(pniv(2) - pniv(1))
         else
           dxdp(j,1) = indefini
         end if
        
         do l=2,llm-1
           if ((x(j,l-1).ne.indefini).and.(x(j,l+1).ne.indefini))then
             dxdp(j,l)= (x(j,l+1)-x(j,l-1))/(pniv(l+1) - pniv(l-1))
         else if((x(j,l+1).ne.indefini).and.(x(j,l).ne.indefini))then
             dxdp(j,l)= (x(j,l+1)-x(j,l))  /(pniv(l+1) - pniv(l))
         else if((x(j,l-1).ne.indefini).and.(x(j,l).ne.indefini))then
             dxdp(j,l)= (x(j,l)-x(j,l-1))  /(pniv(l)   - pniv(l-1))
           else 
             dxdp(j,l)= indefini
           end if
         end do
         if ((x(j,llm).ne.indefini).and.(x(j,llm-1).ne.indefini)) then
         dxdp(j,llm)= (x(j,llm)-x(j,llm-1))/(pniv(llm)-pniv(llm-1))
         else
           dxdp(j,llm)= indefini
         end if
      end do

c     Test
      
      do j=1,jjp1
         do l = 1, llm
           if ( (abs(dxdp(j,l)).gt.(abs(indefini)/100.)).and.
     .          (dxdp(j,l).ne.indefini)) then
              write(*,*) '----> j= ', j , ' l= ' ,  l
              write(*,*) 'dxdp(j,l) ' , dxdp(j,l)
              write(*,*) 'x' , (x(j,k),k=1,llm)
              write(*,*) 'pniv' , pniv
              write(*,*)
           end if
         end do
      end do 



      RETURN
      END
