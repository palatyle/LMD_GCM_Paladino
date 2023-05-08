      SUBROUTINE perosat(ngrid,nlayer,nq,ig, ptimestep,
     $                   pplev, pplay, zt,
     &                   zy, pdqcloud, pdqscloud)
      
      use tracer_mod, only: igcm_h2o2, mmol
      use conc_mod, only: mmean
      USE comcstfi_h
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
c     WARNING : H2O2 mixing ratio is assumed to be q(igcm_h2o2)   
c               index igcm_h2o2 is known from tracer_mod
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

c
c   arguments:
c   ----------

      integer,intent(in) :: ngrid   ! number of atmospheric columns
      integer,intent(in) :: nlayer  ! number of atmospheric layers
      integer,intent(in) :: nq      ! number of tracers
      INTEGER ig
      REAL ptimestep                ! pas de temps physique (s)
      REAL pplev(ngrid,nlayer+1)    ! pression aux inter-couches (Pa)
      REAL pplay(ngrid,nlayer)      ! pression au milieu des couches (Pa)
      REAL zt(nlayer)               ! temperature au centre des couches (K)
                                    ! deja mise a jour dans calchim

c   Traceurs :
      real zy(nlayer,nq)        ! traceur (fraction molaire sortie chimie)
      real pdqcloud(ngrid,nlayer,nq) ! tendance condensation (kg/kg.s-1)
      real pdqscloud(ngrid,nq)         ! flux en surface (kg.m-2.s-1)
      
c   local:
c   ------

      INTEGER l,iq

      REAL zysat(nlayer)
      REAL zynew(nlayer)               ! mole fraction after condensation
      REAL psat_hg                     ! pression saturante (mm Hg)
      REAL psat_hpa                    ! pression saturante (hPa)
      logical,save :: firstcall=.true.

c     Pour diagnostique :
c     ~~~~~~~~~~~~~~~~~
      REAL taucond(ngrid,nlayer)   ! taux de condensation (kg/kg/s-1)

c-----------------------------------------------------------------------
c    1. initialisation/verification
c    ------------------------------
c
       if (firstcall) then
         ! check that there is an h2o2 tracer:
         if (igcm_h2o2.eq.0) then
           write(*,*) "perosat: error; no h2o2 tracer !!!!"
           stop
         endif
         firstcall=.false.
       endif

c    ----------------------------------------------
c    
c       Rapport de melange a saturation dans la couche l :
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c       d'apres Lindner, Planet. Space Sci., 36, 125, 1988.
c       domaine d'application: T < 220 K
c
        do l = 1,nlayer

c       print *,'ig=',ig,' l=',l,' igcm_h2o2=',igcm_h2o2
c       print *,'y=',zy(l,igcm_h2o2),' T=',zt(l)

             zynew(l) = zy(l,igcm_h2o2)

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
        do l=1, nlayer
          taucond(ig,l)=max((zy(l,igcm_h2o2)-zysat(l))*mmol(igcm_h2o2)
     $                         /(mmean(ig,l)*ptimestep),0.)
        end do
c
c       Saturation couche nlay a 2 :  
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~
c
        do l=nlayer,2, -1
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
           pdqscloud(ig,igcm_h2o2)= (zynew(1)-zysat(1))*mmol(igcm_h2o2)
     $   *(pplev(ig,1)-pplev(ig,2))/(mmean(ig,1)*g*ptimestep)
c
           zynew(1)=zysat(1)
        else
          pdqscloud(ig,igcm_h2o2)=0
        end if
c
c       Tendance finale
c       ~~~~~~~~~~~~~~~
c
        do l=1, nlayer
          pdqcloud(ig,l,igcm_h2o2)=(zynew(l) - zy(l,igcm_h2o2))
     &                     *mmol(igcm_h2o2)/(mmean(ig,l)*ptimestep) 
c          print *,'pdqcloud=',pdqcloud(ig,l,igcm_h2o2)
        end do

      RETURN
      END
