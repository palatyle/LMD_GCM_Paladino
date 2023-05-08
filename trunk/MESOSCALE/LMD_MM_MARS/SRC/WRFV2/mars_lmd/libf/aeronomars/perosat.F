      SUBROUTINE perosat(ig, ptimestep,
     $                   pplev, pplay, zt,
     &                   zy, pdqcloud, pdqscloud)
      IMPLICIT NONE

c=======================================================================
c     Treatment of saturation of hydrogen peroxide (H2O2)
c
c     Modif de zq si saturation dans l'atmopshere
c     si zq(ig,l)> zqsat(ig,l) ->    zq(ig,l)=zqsat(ig,l)
c     Le test est effectue de bas en haut. H2O2 condense
c    (si saturation) est remis dans la couche en dessous.
c     H2O2 condense dans la couche du bas est depose a la surface
c
c     WARNING : H2O2 mixing ratio is assumed to be q(i_h2o2)    
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "chimiedata.h"
#include "tracer.h"
#include "conc.h"
c
c   arguments:
c   ----------

      INTEGER ig
      REAL ptimestep                ! pas de temps physique (s)
      REAL pplev(ngridmx,nlayermx+1)! pression aux inter-couches (Pa)
      REAL pplay(ngridmx,nlayermx)  ! pression au milieu des couches (Pa)
      REAL zt(nlayermx)             ! temperature au centre des couches (K)
                                    ! deja mise a jour dans calchim

c   Traceurs :
      real zy(nlayermx,nqmx)        ! traceur (fraction molaire sortie chimie)
      real pdqcloud(ngridmx,nlayermx,nqmx) ! tendance condensation (kg/kg.s-1)
      real pdqscloud(ngridmx,nqmx)         ! flux en surface (kg.m-2.s-1)
      
c   local:
c   ------

      INTEGER l,iq
      INTEGER i_h2o2

      REAL zysat(nlayermx)
      REAL zynew(nlayermx)             ! mole fraction after condensation
      REAL psat_hg                     ! pression saturante (mm Hg)
      REAL psat_hpa                    ! pression saturante (hPa)

c     Pour diagnostique :
c     ~~~~~~~~~~~~~~~~~
      REAL taucond(ngridmx,nlayermx)   ! taux de condensation (kg/kg/s-1)

c-----------------------------------------------------------------------
c    1. initialisation
c    -----------------
c
      do iq=1,nqmx
       if (noms(iq).eq."h2o2") then
          i_h2o2 = iq
c         print *,'H2O2: ',noms(iq),' iq=',iq
c      else
c         print *,noms(iq),' iq=',iq
       endif
      enddo

c    ----------------------------------------------
c    
c       Rapport de melange a saturation dans la couche l :
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c       d'apres Lindner, Planet. Space Sci., 36, 125, 1988.
c       domaine d'application: T < 220 K
c
        do l = 1,nlayermx

c       print *,'ig=',ig,' l=',l,' i_h2o2=',i_h2o2
c       print *,'y=',zy(l,i_h2o2),' T=',zt(l)

             zynew(l) = zy(l,i_h2o2)

             if (zt(l) .le. 220.) then
               psat_hg = 10.**(11.98 - (3422./zt(l)))
               psat_hpa = psat_hg*760./1013.
               zysat(l) = (psat_hpa*100./pplay(ig,l))
             else
               zysat(l) = 1.e+30
             end if

c       print *,'ysat=',zysat(l)

        end do

c       taux de condensation (kg/kg/s-1) dans les differentes couches
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c       (Pour diagnostic seulement !)
c
        do l=1, nlayermx
            taucond(ig,l)= max((zy(l,i_h2o2)-zysat(l))*mmol(i_h2o2)
     $                         /(mmean(ig,l)*ptimestep),0.)
        end do
c
c       Saturation couche nlay a 2 :  
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~
c
        do l=nlayermx,2, -1
           if (zynew(l).gt.zysat(l)) then
              zynew(l-1) =  zynew(l-1) + (zynew(l) - zysat(l))
     &      *(pplev(ig,l)-pplev(ig,l+1))/(pplev(ig,l-1)-pplev(ig,l))

              zynew(l)=zysat(l)
           endif
        enddo
c
c       Saturation couche l=1 
c       ~~~~~~~~~~~~~~~~~~~~~
c
        if (zynew(1).gt.zysat(1)) then
           pdqscloud(ig,i_h2o2)= (zynew(1)-zysat(1))*mmol(i_h2o2)
     $   *(pplev(ig,1)-pplev(ig,2))/(mmean(ig,1)*g*ptimestep)
c
           zynew(1)=zysat(1)
        end if
c
c       Tendance finale
c       ~~~~~~~~~~~~~~~
c
        do l=1, nlayermx
           pdqcloud(ig,l,i_h2o2)=(zynew(l) - zy(l,i_h2o2))*mmol(i_h2o2)
     &                           /(mmean(ig,l)*ptimestep) 
c          print *,'pdqcloud=',pdqcloud(ig,l,i_h2o2)
        end do

      RETURN
      END
