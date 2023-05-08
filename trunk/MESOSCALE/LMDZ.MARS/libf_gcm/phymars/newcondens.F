      SUBROUTINE newcondens(ngrid,nlayer,nq,ptimestep,
     $                  pcapcal,pplay,pplev,ptsrf,pt,
     $                  pphi,pdt,pdu,pdv,pdtsrf,pu,pv,pq,pdq,
     $                  piceco2,psolaralb,pemisurf,
     $                  pdtc,pdtsrfc,pdpsrf,pduc,pdvc,pdqc,
     $                  fluxsurf_sw)
                                                   
       IMPLICIT NONE
c=======================================================================
c   subject:
c   --------
c   Condensation/sublimation of CO2 ice on the ground and in the
c   atmosphere 
c   (Scheme described in Forget et al., Icarus, 1998)
c
c   author:   Francois Forget     1994-1996
c   ------
c
c   input:
c   ----- 
c    ngrid                 nombre de points de verticales
c                          (toutes les boucles de la physique sont au
c                          moins vectorisees sur ngrid)
c    nlayer                nombre de couches
c    pplay(ngrid,nlayer)   Pressure levels 
c    pplev(ngrid,nlayer+1) Pressure levels 
c    pt(ngrid,nlayer)      temperature (en K)
c    ptsrf(ngrid)          temperature de surface
c 
c                    \
c    pdt(ngrid,nlayermx)   \  derivee temporelle physique avant condensation
c                     /  ou sublimation pour  pt,ptsrf
c    pdtsrf(ngrid)   /
c
c   output:
c   -------
c
c    pdpsrf(ngrid)   \  derivee temporelle physique (contribution de
c    pdtc(ngrid,nlayermx) /  la condensation ou sublimation) pour Ps,pt,ptsrf
c    pdtsrfc(ngrid) /
c
c   Entree/sortie
c   -------------
c    
c    piceco2(ngrid) :      quantite de glace co2 au sol (kg/m2)
c    psolaralb(ngrid,2) :  albedo au sol
c    pemisurf(ngrid)     :  emissivite du sol             

c
c=======================================================================
c
c    0.  Declarations :
c    ------------------
c
#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "surfdat.h"
#include "comgeomfi.h"
#include "comvert.h"
#include "paramet.h"
#include "callkeys.h"
#include "tracer.h"

#include "fisice.h"

c-----------------------------------------------------------------------
c    Arguments :
c    ---------
      INTEGER ngrid, nlayer, nq

      REAL ptimestep 
      REAL pplay(ngrid,nlayer),pplev(ngrid,nlayer+1)
      REAL pcapcal(ngrid)
      REAL pt(ngrid,nlayer)
      REAL ptsrf(ngrid)
      REAL pphi(ngrid,nlayer)
      REAL pdt(ngrid,nlayer),pdtsrf(ngrid),pdtc(ngrid,nlayer)
      REAL pdtsrfc(ngrid),pdpsrf(ngrid)
      REAL piceco2(ngrid),psolaralb(ngrid,2),pemisurf(ngrid)

      REAL pu(ngrid,nlayer) ,  pv(ngrid,nlayer)
      REAL pdu(ngrid,nlayer) , pdv(ngrid,nlayer)
      REAL pduc(ngrid,nlayer) , pdvc(ngrid,nlayer)
      REAL pq(ngridmx,nlayer,nq),pdq(ngrid,nlayer,nq)
      REAL pdqc(ngrid,nlayer,nq)
      REAL fluxsurf_sw(ngrid,2) ! added to calculate flux dependent albedo 
c
c    Local variables :
c    -----------------

c   variables used for albedo parametrization       
c --------------------------------------------      
      INTEGER i,j
      REAL Fluxmean(jjp1)
      INTEGER l,ig,iq,icap,nmix
      LOGICAL transparency, fluxdependent
c   flag transparency if you want to make the co2ice semi-transparent       
      PARAMETER(transparency=.true.)
c   flag fluxdependent if you want the co2ice albedo to be dependent on
c   the incident solar flux      
      PARAMETER(fluxdependent=.false.)
      REAL slopy,alpha,constA,constB,constT,albediceF_new(ngridmx)
      REAL zt(ngridmx,nlayermx)
      REAL zcpi
      REAL ztcond (ngridmx,nlayermx)
      REAL ztcondsol(ngridmx) 
      REAL zdiceco2(ngridmx)
c      REAL zcondicea(ngridmx,nlayermx) ! Already defined in fisice.h 
      REAL zcondices(ngridmx)
      REAL zfallice(ngridmx,nlayermx+1) , zfallheat
      REAL zmflux(nlayermx+1)
      REAL zu(nlayermx),zv(nlayermx)
      REAL zq(nlayermx,nqmx),zq1(nlayermx)
      REAL ztsrf(ngridmx)
      REAL ztc(nlayermx), ztm(nlayermx+1) 
      REAL zum(nlayermx+1) , zvm(nlayermx+1)
      REAL zqm(nlayermx+1,nqmx),zqm1(nlayermx+1)
      REAL masse(nlayermx),w(nlayermx+1)
      REAL Sm(nlayermx),Smq(nlayermx,nqmx),mixmas,qmix
      LOGICAL condsub(ngridmx)

c variable speciale diagnostique
      real tconda1(ngridmx,nlayermx)
      real tconda2(ngridmx,nlayermx)
c     REAL zdiceco2a(ngridmx) ! for diagnostic only
      real zdtsig (ngridmx,nlayermx)
      real zdt (ngridmx,nlayermx)


c   local saved variables
      integer ico2
      real qco2min
      save ico2,qco2min
      REAL emisref(ngridmx)
      REAL latcond,tcond1mb
      REAL acond,bcond,ccond,cpice
      REAL albediceF(ngridmx)
      SAVE emisref
      SAVE latcond,acond,bcond,ccond,cpice
      SAVE albediceF
 

      LOGICAL firstcall,firstcall2
      SAVE firstcall,firstcall2
      REAL SSUM
      EXTERNAL SSUM

c      common/scratch/zt,ztcond,zcondicea,zfallice,tconda1
c     ,    ,tconda2,zdtsig,zdt,zu,zv


      DATA latcond,tcond1mb/5.9e5,136.27/
      DATA cpice /1000./
      DATA firstcall/.true./
      DATA firstcall2/.true./

      integer flag

c----------------------------------------------------------------------

c   Initialisation
c   --------------
c
      IF (firstcall) THEN
         bcond=1./tcond1mb
         ccond=cpp/(g*latcond)
         acond=r/latcond

         firstcall=.false.
         PRINT*,'In newcondens:Tcond(P=1mb)=',tcond1mb,' Lcond=',latcond
         PRINT*,'acond,bcond,ccond',acond,bcond,acond

         ico2=0
         if (tracer) then  
c          Prepare Special treatment if one of the tracer is CO2 gas
           do iq=1,nqmx
             if (noms(iq).eq."co2") ico2=iq
           enddo
c          minimum CO2 mix. ratio below which mixing occur with layer above: 
           qco2min =0.75  
         end if
      ENDIF
      zcpi=1./cpp
c
c======================================================================
c    Calcul of CO2 condensation sublimation 
c    ============================================================
c  
c    Used variable :
c       piceco2(ngrid)   :  amount of co2 ice on the ground (kg/m2)
c       zcondicea(ngrid,l):  condensation rate in layer  l  (kg/m2/s)
c       zcondices(ngrid):  condensation rate on the ground (kg/m2/s)
c       zfallice(ngrid,l):amount of ice falling from layer l (kg/m2/s)
c                           
c       pdtc(ngrid,nlayermx)  : dT/dt due to cond/sub
c
c
c     Tendencies set to 0 (except pdtc)
c     -------------------------------------
      DO l=1,nlayer
         DO ig=1,ngrid
           zcondicea(ig,l) = 0.
           zfallice(ig,l) = 0.
           pduc(ig,l)  = 0
           pdvc(ig,l)  = 0
         END DO
      END DO
         
      DO iq=1,nqmx         
      DO l=1,nlayer
         DO ig=1,ngrid
           pdqc(ig,l,iq)  = 0
         END DO
      END DO
      END DO

      DO ig=1,ngrid
         zfallice(ig,nlayer+1) = 0.
         zcondices(ig) = 0.
         pdtsrfc(ig) = 0.
         pdpsrf(ig) = 0.
         condsub(ig) = .false.
         zdiceco2(ig) = 0.
      ENDDO
      zfallheat=0

c     *************************
c     ATMOSPHERIC CONDENSATION
c     *************************

c     forecast of atmospheric temperature zt and frost temperature ztcond
c     --------------------------------------------------------------------

      DO l=1,nlayer
         DO ig=1,ngrid
            zt(ig,l)=pt(ig,l)+ pdt(ig,l)*ptimestep
            ztcond(ig,l)=1./(bcond-acond*log(.0095*pplay(ig,l)))
            if (pplay(ig,l).lt.1e-4) ztcond(ig,l)=0.0 !mars Monica
         ENDDO
      ENDDO
 
c      Condensation/sublimation in the atmosphere
c      ------------------------------------------
c      (calcul of zcondicea , zfallice and pdtc)
c
      DO l=nlayer , 1, -1
         DO ig=1,ngrid
           pdtc(ig,l)=0.
           IF((zt(ig,l).LT.ztcond(ig,l)).or.(zfallice(ig,l+1).gt.0))THEN
               condsub(ig)=.true.
               IF (zfallice(ig,l+1).gt.0) then  
                 zfallheat=zfallice(ig,l+1)*
     &           (pphi(ig,l+1)-pphi(ig,l) +
     &          cpice*(ztcond(ig,l+1)-ztcond(ig,l)))/latcond
               ELSE
                   zfallheat=0.
               ENDIF
               pdtc(ig,l)=(ztcond(ig,l) - zt(ig,l))/ptimestep
               zcondicea(ig,l)=(pplev(ig,l)-pplev(ig,l+1))
     &                        *ccond*pdtc(ig,l)- zfallheat
c              Case when the ice from above sublimes entirely
c              """""""""""""""""""""""""""""""""""""""""""""""
               IF (zfallice(ig,l+1).lt.- zcondicea(ig,l)) then
                  pdtc(ig,l)=(-zfallice(ig,l+1)+zfallheat)/
     &              (ccond*(pplev(ig,l)-pplev(ig,l+1)))
                  zcondicea(ig,l)= -zfallice(ig,l+1)
               END IF

               zfallice(ig,l) = zcondicea(ig,l)+zfallice(ig,l+1)
            END IF
         ENDDO
      ENDDO

c     *************************
c       SURFACE  CONDENSATION
c     *************************

c     forecast of ground temperature ztsrf and frost temperature ztcondsol
c     --------------------------------------------------------------------
      DO ig=1,ngrid
         ztcondsol(ig)=1./(bcond-acond*log(.0095*pplev(ig,1)))
         ztsrf(ig) = ptsrf(ig) + pdtsrf(ig)*ptimestep
      ENDDO

c
c      Condensation/sublimation on the ground 
c      --------------------------------------
c      (calcul of zcondices , pdtsrfc)
c
      DO ig=1,ngrid
         IF(ig.GT.ngrid/2+1) THEN
            icap=2
         ELSE
            icap=1
         ENDIF
       
c
c        Loop on where we have  condensation/ sublimation 
         IF ((ztsrf(ig) .LT. ztcondsol(ig)) .OR.   ! ground cond 
     $       (zfallice(ig,1).NE.0.) .OR.           ! falling snow
     $      ((ztsrf(ig) .GT. ztcondsol(ig)) .AND.  ! ground sublim.
     $      ((piceco2(ig)+zfallice(ig,1)*ptimestep) .NE. 0.))) THEN
            condsub(ig) = .true. 

            IF (zfallice(ig,1).gt.0) then  
                 zfallheat=zfallice(ig,1)*
     &           (pphi(ig,1)- phisfi(ig) +
     &           cpice*(ztcond(ig,1)-ztcondsol(ig)))/(latcond*ptimestep)
            ELSE
                 zfallheat=0.
            ENDIF

c           condensation or partial sublimation of CO2 ice
c           """""""""""""""""""""""""""""""""""""""""""""""
            zcondices(ig)=pcapcal(ig)*(ztcondsol(ig)-ztsrf(ig)) 
     &      /(latcond*ptimestep)         - zfallheat
            pdtsrfc(ig) = (ztcondsol(ig) - ztsrf(ig))/ptimestep
       
c           If the entire CO_2 ice layer sublimes
c           """"""""""""""""""""""""""""""""""""""""""""""""""""
c           (including what has just condensed in the atmosphere)

            IF((piceco2(ig)/ptimestep+zfallice(ig,1)).LE.
     &          -zcondices(ig))THEN
               zcondices(ig) = -piceco2(ig)/ptimestep - zfallice(ig,1) 
               pdtsrfc(ig)=(latcond/pcapcal(ig))*
     &         (zcondices(ig)+zfallheat)
            END IF

c           Changing CO2 ice amount and pressure :
c           """"""""""""""""""""""""""""""""""""

            zdiceco2(ig) = zcondices(ig) + zfallice(ig,1)
            piceco2(ig) = piceco2(ig) + zdiceco2(ig)*ptimestep
            pdpsrf(ig) = -zdiceco2(ig)*g

            IF(ABS(pdpsrf(ig)*ptimestep).GT.pplev(ig,1)) THEN
               PRINT*,'STOP in condens'
               PRINT*,'condensing more than total mass'
               PRINT*,'Grid point ',ig
               PRINT*,'Ps = ',pplev(ig,1)
               PRINT*,'d Ps = ',pdpsrf(ig)
               STOP
            ENDIF
         END IF
      ENDDO

c  ********************************************************************
c  Surface albedo and emissivity of the surface below the snow (emisref)
c  ********************************************************************
c     Prepare the case where albedo varies with insolation:
c     ----------------------------------------------------
      if (fluxdependent) then

c       Calcul du flux moyen (zonal mean)
        do j=1,jjp1
           Fluxmean(j)=0      
           do i=1,iim
               ig=1+(j-2)*iim +i
               if(j.eq.1) ig=1
               if(j.eq.jjp1) ig=ngrid
               Fluxmean(j)=Fluxmean(j)+fluxsurf_sw(ig,1)
     $            +fluxsurf_sw(ig,2)
           enddo
           Fluxmean(j)=Fluxmean(j)/float(iim)
        enddo

c       const A and B used to calculate the albedo which depends on solar flux 
c       albedice=constA+constB*Flux
c       constT = time step to calculate the solar flux when flux decreases         
         constA=0.26
c        constA=0.33
c        constA=0.186
         constB=0.00187
         constT=10
      endif

c     Calcul de l'albedo
c     ------------------
         do ig =1,ngrid
           IF(ig.GT.ngrid/2+1) THEN
              icap=2
           ELSE
              icap=1
           ENDIF
           IF(firstcall2) THEN
           albediceF(ig)=albedice(icap)
           ENDIF
c       test on the existence of co2ice ccccccccccccccccccccccc
           if(.not.piceco2(ig).ge.0.) THEN
              if(piceco2(ig).le.-5.e-8) print*,
     $         'WARNING newcondens piceco2(',ig,')=', piceco2(ig)
               piceco2(ig)=0.
           endif
c       if there is still co2ice        ccccccccccccccccccccccc
           if (piceco2(ig).gt.0) then
              emisref(ig) = emisice(icap)

c     if flux dependent albedo is used 
c     --------------------------------
             if (fluxdependent) then
                j=INT((ig-2)/iim)+2
                if(ig.eq.1) j=1
                if(ig.eq.ngrid) j=jjp1
c                albediceF_new(ig)=MIN(constA+constB*Fluxmean(j),
c     $          constA+constB*250)
                  albediceF_new(ig)=constA+constB*Fluxmean(j)
                if (albediceF(ig).gt.albediceF_new(ig)) then
                   albediceF(ig)=albediceF(ig)+ ptimestep/(daysec*
     $             constT)*(albediceF_new(ig)-albediceF(ig))
                else
                   albediceF(ig)=albediceF_new(ig)
                endif
c               if part of the ice is transparent
c               slopy = pente de la droite:  alpha = y*co2ice/1620
c               pour une valeur superieur a une epaisseur de glace donnee
c               ici, epaisseur limite = 10cm
                if (transparency) then
                slopy=1/(1620*0.10)
                alpha=MIN(slopy*piceco2(ig),1.)
                psolaralb(ig,1) = alpha*albediceF(ig)
     $                           +(1-alpha)*albedodat(ig) 
                psolaralb(ig,2) = psolaralb(ig,1) 
                else
                psolaralb(ig,1) = albediceF(ig) 
                psolaralb(ig,2) = psolaralb(ig,1) 
                endif
             else
c           transparency set to true and fluxdependent set to false            
                if (transparency) then
                slopy=1/(1620*0.10)
                alpha=MIN(slopy*piceco2(ig),1.)
                psolaralb(ig,1) = alpha*albedice(icap)
     $                           +(1-alpha)*albedodat(ig) 
                psolaralb(ig,2) = psolaralb(ig,1) 
                else
c           simplest case: transparency and flux dependent set to false
                psolaralb(ig,1) = albedice(icap)
                psolaralb(ig,2) = albedice(icap)
                endif
             endif    
c         no more co2ice, albedo = ground albedo
           else
              psolaralb(ig,1) = albedodat(ig)
              psolaralb(ig,2) = albedodat(ig)
              emisref(ig) = emissiv
              pemisurf(ig) = emissiv
           endif
         end do ! end of the ig loop

         firstcall2=.false.
c ***************************************************************
c  Correction to account for redistribution between sigma or hybrid 
c  layers when changing surface pressure (and warming/cooling
c  of the CO2 currently changing phase).
c *************************************************************

      DO ig=1,ngrid
        if (condsub(ig)) then
           do l=1,nlayer
             ztc(l)  =zt(ig,l)   +pdtc(ig,l)  *ptimestep
             zu(l)   =pu(ig,l)   +pdu( ig,l)  *ptimestep
             zv(l)   =pv(ig,l)   +pdv( ig,l)  *ptimestep
            do iq=1,nqmx
             zq(l,iq)=pq(ig,l,iq)+pdq(ig,l,iq)*ptimestep
            enddo
           end do

c  Mass fluxes through the sigma levels (kg.m-2.s-1)  (>0 when up)
c  """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

            zmflux(1) = -zcondices(ig)
            DO l=1,nlayer
             zmflux(l+1) = zmflux(l) -zcondicea(ig,l)
     &        + (bp(l)-bp(l+1))*(zfallice(ig,1)-zmflux(1))
c zmflux set to 0 if very low to avoid: top layer is disappearing in v1ld  
          if (abs(zmflux(l+1)).lt.1E-13.OR.bp(l+1).eq.0.) zmflux(l+1)=0.
            END DO

c Mass of each layer
c ------------------ 
            DO l=1,nlayer
              masse(l)=(pplev(ig,l) - pplev(ig,l+1))/g
            END DO


c  Corresponding fluxes for T,U,V,Q
c  """"""""""""""""""""""""""""""""

c           averaging operator for TRANSPORT  
c           """"""""""""""""""""""""""""""""
c           Value transfert at the surface interface when condensation
c           sublimation:
            ztm(1) = ztsrf(ig) + pdtsrfc(ig)*ptimestep
            zum(1) = 0  
            zvm(1) = 0  
            do iq=1,nqmx
              zqm(1,iq)=0. ! most tracer do not condense !
            enddo
c           Special case if one of the tracer is CO2 gas
            if (ico2.ne.0) zqm(1,ico2)=1. ! flux is 100% CO2

c           Van Leer scheme:
            DO l=1,nlayer+1
                w(l)=-zmflux(l)*ptimestep
            END DO
            call vl1d(ztc,2.,masse,w,ztm) 
            call vl1d(zu ,2.,masse,w,zum) 
            call vl1d(zv ,2.,masse,w,zvm) 
            do iq=1,nqmx
             do l=1,nlayer
              zq1(l)=zq(l,iq)
             enddo
             zqm1(1)=zqm(1,iq)
             call vl1d(zq1,2.,masse,w,zqm1)
             do l=2,nlayer
              zq( l,iq)=zq1(l)
              zqm(l,iq)=zqm1(l)
             enddo
            enddo

c           Surface condensation affects low winds
            if (zmflux(1).lt.0) then 
                zum(1)= zu(1) *  (w(1)/masse(1))
                zvm(1)= zv(1) *  (w(1)/masse(1))
                if (w(1).gt.masse(1)) then ! ensure numerical stability
                  zum(1)= (zu(1)-zum(2))*masse(1)/w(1) + zum(2)
                  zvm(1)= (zv(1)-zvm(2))*masse(1)/w(1) + zvm(2)
                end if
            end if
                    
            ztm(nlayer+1)= ztc(nlayer) ! should not be used, but... 
            zum(nlayer+1)= zu(nlayer)  ! should not be used, but...
            zvm(nlayer+1)= zv(nlayer)  ! should not be used, but...
            do iq=1,nqmx
             zqm(nlayer+1,iq)= zq(nlayer,iq)
            enddo
  
c           Tendencies on T, U, V, Q 
c           """"""""""""""""""""""""
            DO l=1,nlayer
 
c             Tendencies on T
                zdtsig(ig,l) = (1/masse(l)) *
     &        ( zmflux(l)*(ztm(l) - ztc(l)) 
     &        - zmflux(l+1)*(ztm(l+1) - ztc(l))
     &        + zcondicea(ig,l)*(ztcond(ig,l)-ztc(l))  )
                pdtc(ig,l) =  pdtc(ig,l) + zdtsig(ig,l)

c             Tendencies on U
                pduc(ig,l)   = (1/masse(l)) *
     &        ( zmflux(l)*(zum(l) - zu(l))
     &        - zmflux(l+1)*(zum(l+1) - zu(l)) )


c             Tendencies on V
                pdvc(ig,l)   = (1/masse(l)) *
     &        ( zmflux(l)*(zvm(l) - zv(l))
     &        - zmflux(l+1)*(zvm(l+1) - zv(l)) )

            END DO

c           Tendencies on Q
            do iq=1,nqmx
              if (noms(iq).eq.'co2') then 
c               SPECIAL Case when the tracer IS CO2 :
                DO l=1,nlayer
                  pdqc(ig,l,iq)= (1/masse(l)) *
     &           ( zmflux(l)*(zqm(l,iq) - zq(l,iq))
     &           - zmflux(l+1)*(zqm(l+1,iq) - zq(l,iq)) 
     &           + zcondicea(ig,l)*(zq(l,iq)-1.) )
                END DO
              else
                DO l=1,nlayer
                  pdqc(ig,l,iq)= (1/masse(l)) *
     &           ( zmflux(l)*(zqm(l,iq) - zq(l,iq))
     &           - zmflux(l+1)*(zqm(l+1,iq) - zq(l,iq)) 
     &           + zcondicea(ig,l)*zq(l,iq) )
                END DO
              end if
            enddo

c           --------------------------------------------------------
c           Roughly Simulate Molecular mixing when CO2 is too depleted by
c           Surface condensation (mixing starts if qco2 < qco2min )  
c           FF 06/2004
c           WARNING : this is now done in convadj, better (FF 02/2005)
c           --------------------------------------------------------
            flag=0  ! now done in convadj : must be =0
            if (flag.eq.1) then 
            if(ico2.gt.0) then   ! relevant only if one tracer is CO2
             if(pq(ig,1,ico2)+(pdq(ig,1,ico2)+pdqc(ig,1,ico2))*ptimestep
     &       .lt.qco2min) then 
                do iq=1,nqmx
                  zq(1,iq)=pq(ig,1,iq)
     &                     +(pdq(ig,1,iq)+pdqc(ig,1,iq))*ptimestep
                  Smq(1,iq) = masse(1)*zq(1,iq)
                end do
                Sm(1)  = masse(1)
                do l =2,nlayermx
                  do iq=1,nqmx
                     zq(l,iq)=pq(ig,l,iq)
     &                        +(pdq(ig,l,iq)+pdqc(ig,l,iq))*ptimestep
                     smq(l,iq) = smq(l-1,iq) + masse(l)*zq(l,iq)
                  end do 
                  sm(l) = sm(l-1) + masse(l)
                  if(zq(l,ico2).gt.qco2min) then 
c                   mixmas: mass of atmosphere that must be mixed to reach qco2min
                    mixmas = (sm(l-1)*zq(l,ico2)-Smq(l-1,ico2))
     &                      /(zq(l,ico2)-qco2min)
                    if((mixmas.le.sm(l)))then
c                      OK if mixed mass less than mass of layers affected
                       nmix=l   ! number of layer affected by mixing
                       goto 99 
                    end if
                  end if
                 end do
 99              continue 
                 do iq=1,nqmx
                   qmix=zq(nmix,iq) 
     &             +(Smq(nmix-1,iq)-zq(nmix,iq)*Sm(nmix-1))/mixmas
                   do l=1,nmix-1
                      pdqc(ig,l,iq)=
     &                  (qmix-pq(ig,l,iq))/ptimestep - pdq(ig,l,iq)
                   end do 
c                  layer only partly mixed :
                   pdqc(ig,nmix,iq)=(
     &             qmix+(Sm(nmix)-mixmas)*(zq(nmix,iq)-qmix)/masse(nmix)
     &             -pq(ig,nmix,iq))/ptimestep -pdq(ig,nmix,iq)

                 end do 
             end if
            end if

        endif   ! (flag.eq.1)
       end if ! if (condsub)
      END DO  ! loop on ig 

c ***************************************************************
c   CO2 snow / clouds scheme
c ***************************************************************

      call co2snow(ngrid,nlayer,ptimestep,emisref,condsub,pplev,
     &        zcondicea,zcondices,zfallice,pemisurf)

c ***************************************************************
c Ecriture des diagnostiques
c ***************************************************************

c     DO l=1,nlayer
c        DO ig=1,ngrid
c         Taux de cond en kg.m-2.pa-1.s-1
c          tconda1(ig,l)=zcondicea(ig,l)/(pplev(ig,l)-pplev(ig,l+1))
c          Taux de cond en kg.m-3.s-1
c          tconda2(ig,l)=tconda1(ig,l)*pplay(ig,l)*g/(r*pt(ig,l))
c        END DO
c     END DO
c     call WRITEDIAGFI(ngridmx,'tconda1',
c    &'Taux de condensation CO2 atmospherique /Pa',
c    & 'kg.m-2.Pa-1.s-1',3,tconda1)
c     call WRITEDIAGFI(ngridmx,'tconda2',
c    &'Taux de condensation CO2 atmospherique /m',
c    & 'kg.m-3.s-1',3,tconda2)

      return
      end



c *****************************************************************
      SUBROUTINE vl1d(q,pente_max,masse,w,qm)
c
c    
c     Operateur de moyenne inter-couche pour calcul de transport type
c     Van-Leer " pseudo amont " dans la verticale
c    q,w sont des arguments d'entree  pour le s-pg ....
c    masse : masse de la couche Dp/g
c    w : masse d'atm ``transferee'' a chaque pas de temps (kg.m-2)
c    pente_max = 2 conseillee
c
c
c   --------------------------------------------------------------------
      IMPLICIT NONE

#include "dimensions.h"

c
c
c
c   Arguments:
c   ----------
      real masse(llm),pente_max
      REAL q(llm),qm(llm+1)
      REAL w(llm+1)
c
c      Local 
c   ---------
c
      INTEGER l
c
      real dzq(llm),dzqw(llm),adzqw(llm),dzqmax
      real sigw, Mtot, MQtot
      integer m 
c     integer ismax,ismin 


c    On oriente tout dans le sens de la pression 
c     W > 0 WHEN DOWN !!!!!!!!!!!!!

      do l=2,llm
            dzqw(l)=q(l-1)-q(l)
            adzqw(l)=abs(dzqw(l))
      enddo

      do l=2,llm-1
            if(dzqw(l)*dzqw(l+1).gt.0.) then
                dzq(l)=0.5*(dzqw(l)+dzqw(l+1))
            else
                dzq(l)=0.
            endif
            dzqmax=pente_max*min(adzqw(l),adzqw(l+1))
            dzq(l)=sign(min(abs(dzq(l)),dzqmax),dzq(l))
      enddo

         dzq(1)=0.
         dzq(llm)=0.

       do l = 1,llm-1

c         Regular scheme (transfered mass < layer mass)
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(w(l+1).gt.0. .and. w(l+1).le.masse(l+1)) then
             sigw=w(l+1)/masse(l+1)
             qm(l+1)=(q(l+1)+0.5*(1.-sigw)*dzq(l+1))
          else if(w(l+1).le.0. .and. -w(l+1).le.masse(l)) then
             sigw=w(l+1)/masse(l)
             qm(l+1)=(q(l)-0.5*(1.+sigw)*dzq(l))

c         Extended scheme (transfered mass > layer mass)
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          else if(w(l+1).gt.0.) then
             m=l+1
             Mtot = masse(m)
             MQtot = masse(m)*q(m)
             do while ((m.lt.llm).and.(w(l+1).gt.(Mtot+masse(m+1))))
                m=m+1
                Mtot = Mtot + masse(m)
                MQtot = MQtot + masse(m)*q(m)
             end do
             if (m.lt.llm) then
                sigw=(w(l+1)-Mtot)/masse(m+1)
                qm(l+1)= (1/w(l+1))*(MQtot + (w(l+1)-Mtot)*
     &          (q(m+1)+0.5*(1.-sigw)*dzq(m+1)) )
             else
                w(l+1) = Mtot
                qm(l+1) = Mqtot / Mtot
                write(*,*) 'top layer is disapearing !'
                stop
             end if
          else      ! if(w(l+1).lt.0) 
             m = l-1 
             Mtot = masse(m+1)
             MQtot = masse(m+1)*q(m+1)
             do while (m.gt.0)
              do while (-w(l+1).gt.(Mtot+masse(m)))
                m=m-1
                Mtot = Mtot + masse(m+1)
                MQtot = MQtot + masse(m+1)*q(m+1)
              end do
             end do
             if (m.gt.0) then
                sigw=(w(l+1)+Mtot)/masse(m)
                qm(l+1)= (-1/w(l+1))*(MQtot + (-w(l+1)-Mtot)* 
     &          (q(m)-0.5*(1.+sigw)*dzq(m))  )
             else
                qm(l+1)= (-1/w(l+1))*(MQtot + (-w(l+1)-Mtot)*qm(1))
             end if   
          end if
       enddo

c     boundary conditions (not used in newcondens !!)
c         qm(llm+1)=0.
c         if(w(1).gt.0.) then
c            qm(1)=q(1)
c         else 
c           qm(1)=0.
c         end if

      return
      end
