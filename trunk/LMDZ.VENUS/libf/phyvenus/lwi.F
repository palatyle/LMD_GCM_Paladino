      subroutine lwi(nl,netrad,deltapsi,deltap,temp,coolrate)

      use dimphy
      use cpdet_phy_mod, only: cpdet
      implicit none


#include "YOMCST.h"
#include "timerad.h"
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                             -   lwi    -    
C
C     PURPOSE:       Schema semi - implicite 
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c************************************************************************
c
c        0.    Declarations
c              ------------
c
c-------------------------------------------------------------------------
c        0.1   Arguments
c              ---------
c Inputs
c.......

      integer nl
      real    deltapsi(0:nl+1,0:nl+1) ! Dpsi/DT = sum(nu)[ksi DB/DT] 
      real    netrad(0:nl)     ! radiative budget (W/m**2)
      real    deltap(nl) !epaisseur de la couche en pression (Pa) 
! ADAPTATION GCM POUR CP(T)
      real    temp(nl)   !temperature de la couche 
      
c Outputs
c........

      real    coolrate(nl)
      
c-------------------------------------------------------------------------
c        0.2   local arrays
c              ------------
c
      real    di(nl) 
     .      , spec(nl) 
     .      , hi(nl) 
     .      , bi(nl) 

      real    ci(nl) 
     .      , ai(nl) 

      real    deltat

      real    semi,semit,denom

      integer i,jl,band

c************************************************************************
c
c        1.    Initialisations
c              ---------------
c
c-----------------------------------------------------------------------
        deltat = dtimerad
c       print*,'SEMI = ',semi, '(expl:0  semi-implicite:0.5  impl:1)'
c        semi  = 0.5
c       semi = 0.5
        semi = 1.
        semit = semi * deltat

c       print*,'dtimerad,deltat,semit:',dtimerad,deltat,semit

c************************************************************************
c
c        2.    
c              ---------------
c
c-------------------------------------------------------------------------
c        2.1   Calcul des di
c              -------------
c

      do i = 1 , nl-1
        spec(i) = 
     .      ( deltapsi(i,nl+1)
     .      + deltapsi(i,i+1)
     .      + deltapsi(i,i-1) )
      
        di(i) =  1. + semit * (RG/deltap(i)/cpdet(temp(i))) * spec(i)
c     print*,i,' di(i)=',di(i)
      enddo
         

c couche nl
c ------------
c      , on enleve i,i+1 sinon on a 2 fois le cooling2space
 
        spec(nl) = 
     .      ( deltapsi(nl,nl+1)
     .      + deltapsi(nl,nl-1) )
     
        di(nl) =  1. + semit * (RG/deltap(nl)/cpdet(temp(nl)))*spec(nl)

c-------------------------------------------------------------------------
c        2.2   Calcul des hi
c              -------------
c

      do i = 1 , nl-1
      
        spec(i) = deltapsi(i+1,i)

        hi(i) = - semit * (RG/deltap(i)/cpdet(temp(i))) * spec(i)

      enddo
                        
c     print*,'hi(i)',hi(i)

c-------------------------------------------------------------------------
c        2.3   Calcul des bi
c              -------------
c


      do i = 2 , nl
      
        spec(i) = deltapsi(i-1,i)
       
        bi(i) = - semit * (RG/deltap(i)/cpdet(temp(i))) * spec(i)
        
      enddo

c     print*,'bi(i)',bi(i)

c couche 1
c --------
c  tant qu'on a pas un calcul propre de deltab(0) qui tienne compte de 
c    la discontinuite de temperature au sol , on met  b1 = 0

      
        bi(1) = 0 

c-------------------------------------------------------------------------
c        2.4   
c              -------------
c

c couche nl
c ------------

c     ci(nl) = ((RG/CP) * netrad(nl) / deltap(nl))
c    .                   / di(nl)
      ci(nl) = netrad(nl) *RG/cpdet(temp(nl)) / deltap(nl)

      ai(nl) = - bi(nl) / di(nl)

      do i = nl-1 , 1 , -1
        denom = di(i) + hi(i) * ai(i+1)

c       ci(i) = (  (RG/CP) * netrad(i) / deltap(i)
c    .               - hi(i) * ci(i+1)  )  / denom
        ci(i) = netrad(i) *RG/cpdet(temp(i)) / deltap(i)
 
        ai(i) = -bi(i) / denom
      enddo

c-------------------------------------------------------------------------
c        2.5   
c              -------------
c

      coolrate(1) = ci(1)

      do i = 2 , nl
           coolrate(i) = ci(i) + ai(i) * coolrate(i-1)
c       print*,i,' coolrate(i)=',coolrate(i)
      enddo

c-------------------------------------------------------------------------

      RETURN
      END
