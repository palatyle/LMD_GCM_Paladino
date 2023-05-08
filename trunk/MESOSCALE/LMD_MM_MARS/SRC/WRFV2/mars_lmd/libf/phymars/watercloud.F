      SUBROUTINE watercloud(ngrid,nlay, ptimestep,
     &                pplev,pplay,pdpsrf,pzlev,pzlay,pt,pdt,
     &                pq,pdq,pdqcloud,qsurf,pdqscloud,pdtcloud,
     &                nq,nsize,tau,icount,zls)
      IMPLICIT NONE

c=======================================================================
c     Treatment of saturation of water vapor
c
c
c     Modif de zq si saturation dans l'atmopshere
c     si zq(ig,l)> zqsat(ig,l) ->    zq(ig,l)=zqsat(ig,l)
c     Le test est effectue de bas en haut. L'eau condensee
c    (si saturation) est remise dans la couche en dessous.
c     L'eau condensee dans la couche du bas est deposee a la surface
c
c     WARNING : water vapor mixing ratio is assumed to be 
c     q(nqmx)
c
c    Modification: Franck Montmessin water ice scheme
c
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

#include "dimensions.h"
#include "dimphys.h"
#include "comcstfi.h"
#include "callkeys.h"
#include "tracer.h"

#include "fisice.h"
#include "comgeomfi.h"
c
c   arguments:
c   ----------

      INTEGER ngrid,nlay,nsize
      INTEGER icount
      REAL ptimestep             ! pas de temps physique (s)
      REAL pplev(ngrid,nlay+1)   ! pression aux inter-couches (Pa)
      REAL pplay(ngrid,nlay)     ! pression au milieu des couches (Pa)
      REAL pzlev(ngrid,nlay+1)   ! altitude at layer boundaries
      REAL pzlay(ngrid,nlay)     ! altitude at the middle of the layers
      REAL pdpsrf(ngrid)         ! tendance surf pressure
      REAL pt(ngrid,nlay)        ! temperature au centre des couches (K)
      REAL pdt(ngrid,nlay)       ! tendance temperature des autres param.
      REAL pdtcloud(ngrid,nlay)  ! tendance temperature due a la chaleur latente
      REAL tau(ngridmx,nsize)
      REAL zls

c    Traceurs :
      integer nq                   ! nombre de traceurs 
      real pq(ngrid,nlay,nq)       ! traceur (kg/kg)
      real pdq(ngrid,nlay,nq)      ! tendance avant condensation  (kg/kg.s-1)
      real pdqcloud(ngrid,nlay,nq) ! tendance de la condesation H2O(kg/kg.s-1)
      real qsurf(ngrid,nq)
      real pdqscloud(ngrid,nq)     ! flux en surface (kg.m-2.s-1)
      
c   local:
c   ------

      REAL CBRT
      EXTERNAL CBRT
      INTEGER ig,l

      LOGICAL firstcall,improved
      SAVE firstcall

      REAL zq(ngridmx,nlayermx,nqmx)
      REAL zq0(ngridmx,nlayermx,nqmx)
      REAL zqsat(ngridmx,nlayermx)
      REAL zt(ngridmx,nlayermx)


      REAL masse (ngridmx,nlayermx)
      REAL epaisseur (ngridmx,nlayermx)
      REAL dustcores(ngridmx,nlayermx)  !Dust number density (#/kg) 
      REAL rfinal        ! Ice crystal radius after condensation(m)
      REAL seq           ! Equilibrium saturation ration (accounting for curvature effect)
      REAL dzq           ! masse de glace echangee (kg/kg)
      REAL lw       !Latent heat of sublimation (J.kg-1) 
      REAL To

      REAL Ctot
      REAL*8 ph2o,satu
      REAL gr,Cste,up,dwn,newvap

      DATA firstcall/.true./
cc      DATA improved/.false./      ! Actionne une microphysique plus raffinee
      DATA improved/.true./

c     Reference temperature, T=273,15 K
      data To/273.15/

      SAVE improved

c     POur diagnostique :
c     ~~~~~~~~~~~~~~~~~
      REAL taucond(ngridmx,nlayermx)   ! taux de condensation (kg/kg/s-1)
ccc****WRF
c      REAL mtot(ngridmx)               ! Total mass of water vapor (kg/m2)
c      REAL icetot(ngridmx)             ! Total mass of water ice (kg/m2)
c      REAL rave(ngridmx)               ! Mean crystal radius in a column (m)
ccc****WRF
      REAL cloudmass(ngridmx,nlayermx) ! Integrated mass of water ice in each layer (kg/m2)
      REAL waterhem(2,3)

      REAL icetot2pm(ngridmx)
      REAL mtot2pm(ngridmx)
      REAL rave2pm(ngridmx)
      REAL qsurf2pm(ngridmx)
      INTEGER i,nit
      REAL ecart,loctime
      integer univtime,lon2pm(ngridmx)

      integer ig_vl2
      real latvl2,lonvl2

      logical tesfield  ! output of TES like data
      data tesfield /.false./

c    ** un petit test de coherence
c       --------------------------

      IF (firstcall) THEN
         IF(ngrid.NE.ngridmx) THEN
            PRINT*,'STOP dans watercloud'
            PRINT*,'probleme de dimensions :'
            PRINT*,'ngrid  =',ngrid
            PRINT*,'ngridmx  =',ngridmx
            STOP
         ENDIF
         firstcall=.false.
         if (tesfield) then
           latvl2= 47.57
           lonvl2= 134.26
           ig_vl2= 1+ int( (1.5-(latvl2-90.)*jjm/180.)  -2 )*iim +
     &              int(1.5+(lonvl2+180)*iim/360.)
           open(201,file='VL2_data')
           open(199,file='water_budgets')
           print*,'File water_budgets created'
c          Find the universal time corresponding to 2PM (TES obs)
           print*,'Universal time corresponding to 2PM (TES obs)'
           print*,'of each longitude - Used when iceparty is on.'
           print*,'2PM outputs are written in a specific file.'
           print*,'|   Longitude   |   U.T = 2PM LT   |'
           nit= 48
           do ig=1,ngridmx
             ecart=24.
             do i=1,nit-1
               loctime = i/float(nit)*24. + long(ig)*12./pi
               if (loctime.gt.24.) loctime = loctime - 24.
               if (loctime.lt.0)   loctime = loctime + 24.
               if (abs(14.-loctime).lt.ecart) then
                 ecart      = abs(14.-loctime)
                 lon2pm(ig) = i
               endif
             enddo
             if(int(lati(ig)*180./pi).eq.0) 
     *             print*,long(ig)/pi*180.,'|',lon2pm(ig)/2.
           enddo
           open(200,file='2pm_outputs',form='unformatted')
           print*,'File 2pm_outputs created'
           write(200) long,lati
         endif ! (tesfield)
      ENDIF


c-----------------------------------------------------------------------
c    1. initialisation
c    -----------------

c    On "update" la valeur de q(nqmx) (water vapor) et temperature.
c    On effectue qqes calculs preliminaires sur les couches : 
c    masse (kg.m-2), epaisseur(m).

      do l=1,nlay
          do ig=1,ngrid
            zq(ig,l,nq)=pq(ig,l,nq)+pdq(ig,l,nq)*ptimestep
            zq(ig,l,nq)=max(zq(ig,l,nq),1.E-30) ! FF 12/2004 
            zq0(ig,l,nq)=zq(ig,l,nq)
            zt(ig,l)=pt(ig,l)+ pdt(ig,l)*ptimestep
            masse(ig,l)=(pplev(ig,l) - pplev(ig,l+1)) /g
            epaisseur(ig,l)= pzlev(ig,l+1) - pzlev(ig,l)
            if(iceparty) then
              zq(ig,l,nq-1)=pq(ig,l,nq-1)+pdq(ig,l,nq-1)*ptimestep
              zq(ig,l,nq-1)=max(zq(ig,l,nq-1),0.) ! FF 12/2004 
              zq0(ig,l,nq-1)=zq(ig,l,nq-1)
              rdust(ig,l)= max(.8e-6*exp(-pzlay(ig,l)/18000.),1.e-9)
c         Calcul du rayon moyen des particules de glace.
c         Hypothese : Dans une couche, la glace presente se
c         repartit uniformement autour du nbre de poussieres
c         dont le rayon moyen est prescrit par rdust.
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if ((dustbin.eq.1).and..not.(doubleq)) then
                zq(ig,l,1)=pq(ig,l,1)+pdq(ig,l,1)*ptimestep
                dustcores(ig,l)=max(zq(ig,l,1)/
     &           (rho_dust*4./3.*pi*radius(1)**3.),1.e-9)
              else
                dustcores(ig,l)=( epaisseur(ig,l)/masse(ig,l) ) *
     &       2.e+6/0.1*max(tau(ig,1),0.001)*exp(-pzlay(ig,l)/10000.)

!!!PATCH WATER
c            TEMPORAIRE : réduction du nombre de nuclei FF 04/2008 :
                dustcores(ig,l) = dustcores(ig,l) / 8.
                ! reduction facteur 2

              endif
              rice(ig,l)=CBRT( ( zq(ig,l,nq-1)/rho_ice+
     &         dustcores(ig,l)*(4./3.)*pi*rdust(ig,l)**3. )
     &         / (dustcores(ig,l)*4./3.*pi) )
              rice(ig,l)=max(rice(ig,l),rdust(ig,l))
            endif
          enddo
      enddo

      call zerophys(ngrid*nq,pdqscloud)
      call zerophys(ngrid*nlay*nq,pdqcloud)
      call zerophys(ngrid*nlay,pdtcloud)
      call zerophys(ngrid,mtot)

      call zerophys(ngrid,icetot)
      call zerophys(ngrid,rave)
      call zerophys(ngrid,cloudmass)

c    ----------------------------------------------
c
c
c       Rapport de melange a saturation dans la couche l : -------
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        call watersat(ngridmx*nlayermx,zt,pplay,zqsat)

c       taux de condensation (kg/kg/s-1) dans les differentes couches
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c       (Pour diagnostic seulement !)
        if(.not.iceparty)then
        do l=1, nlay
          do ig=1,ngridmx
            taucond(ig,l)= max((zq(ig,l,nq)-zqsat(ig,l))/ptimestep,0.)
          end do
        end do
        endif


        if(iceparty) then
          do l=1,nlay
            do ig=1,ngrid

            If (improved) then
c           Improved microphysics scheme
c           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

             Ctot = zq(ig,l,nq) + zq(ig,l,nq-1)
             ph2o = zq(ig,l,nq) * 44. / 18. * pplay(ig,l)
             satu = zq(ig,l,nq) / zqsat(ig,l)

             call growthrate(ptimestep,zt(ig,l),pplay(ig,l),
     &        ph2o,ph2o/satu,seq,rice(ig,l),gr)
             Cste = ptimestep * 4. * pi * rice(ig,l)
     *              * rho_ice * dustcores(ig,l)
             up   = zq(ig,l,nq) + Cste * gr * seq
             dwn  =    1.       + Cste * gr / zqsat(ig,l)
             newvap = min(up/dwn,Ctot)

             gr  = gr * ( newvap/zqsat(ig,l) - seq )
             dzq = min( max( Cste * gr,-zq(ig,l,nq-1) ) 
     *             , zq(ig,l,nq) )

c            Nucleation (sat ratio must be larger than a critical value)
c            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             if (satu.gt.1.) then
               if (satu.le.1.4.and.zq(ig,l,nq-1).lt.1.e-8) 
     *            dzq = 0.
             endif

            Else
c           Old version
c           ~~~~~~~~~~~
             if (zq(ig,l,nq).ge.zqsat(ig,l))then  !  Condensation
               dzq=zq(ig,l,nq)-zqsat(ig,l)                
             elseif(zq(ig,l,nq).lt.zqsat(ig,l))then  ! Sublimation
               dzq=-min(zqsat(ig,l)-zq(ig,l,nq),zq(ig,l,nq-1))
             endif

            Endif

c           Water Mass change
c           ~~~~~~~~~~~~~~~~~
            zq(ig,l,nq-1)=zq(ig,l,nq-1)+dzq
            zq(ig,l,nq)=zq(ig,l,nq)-dzq

            rice(ig,l)=max( CBRT ( (zq(ig,l,nq-1)/rho_ice
     &       +dustcores(ig,l)*(4./3.)*pi*rdust(ig,l)**3.)
     &       /(dustcores(ig,l)*4./3.*pi) ), rdust(ig,l)) 

            if(activice.and.pclc(ig).gt.0)
     &        rice(ig,l)=rice(ig,l)/CBRT(pclc(ig))
              
            enddo
          enddo

        else   ! if not iceparty

c         Saturation couche nlay a 2 :
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~
          do l=nlay,2, -1
             do ig=1,ngrid
             if (zq(ig,l,nq).gt.zqsat(ig,l))then 
               zq(ig,l-1,nq) =  zq(ig,l-1,nq)+(zq(ig,l,nq)-zqsat(ig,l))
     &          *(pplev(ig,l)-pplev(ig,l+1))/(pplev(ig,l-1)-pplev(ig,l))
               zq(ig,l,nq)=zqsat(ig,l)
             endif
             enddo
          enddo

c       Saturation couche l=1 si pas iceparty 
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          do ig=1,ngridmx
              if (zq(ig,1,nq).gt.zqsat(ig,1))then
                 pdqscloud(ig,nq)= (zq(ig,1,nq)-zqsat(ig,1))
     &         *(pplev(ig,1)-pplev(ig,2))/(g*ptimestep)
                 zq(ig,1,nq)=zqsat(ig,1)
              endif
          end do

        endif   ! (iceparty)

c       Tendance finale
c       ~~~~~~~~~~~~~~~
        do l=1, nlay
          do ig=1,ngridmx
            pdqcloud(ig,l,nq) = (zq(ig,l,nq) - zq0(ig,l,nq))/ptimestep
            if(iceparty) then
              pdqcloud(ig,l,nq-1) =
     &        (zq(ig,l,nq-1) - zq0(ig,l,nq-1))/ptimestep
            endif
            lw=(2834.3-0.28*(zt(ig,l)-To)-0.004*(zt(ig,l)-To)**2)*1.e+3
            pdtcloud(ig,l)=-pdqcloud(ig,l,nq)*lw/cpp
          end do 
        end do

c       A correction if a lot of subliming CO2 fills the 1st layer FF04/2005
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c       Then that should not affect the ice particle radius
        do ig=1,ngridmx
          if(pdpsrf(ig)*ptimestep.gt.0.9*(pplev(ig,1)-pplev(ig,2)))then
            if(pdpsrf(ig)*ptimestep.gt.0.9*(pplev(ig,1)-pplev(ig,3)))
     &      rice(ig,2)=rice(ig,3) 
            rice(ig,1)=rice(ig,2)
          end if
        end do
        
c**************************************************
c       Output  
c**************************************************
        univtime = icount - int(float(icount)/48.) * 48

c       if "tesfield"
c       ~~~~~~~~~~~~~
        if (tesfield.and.mod(icount,12).eq.0) 
     .          call zerophys(6,waterhem)
c****WRF: ligne trop longue        
         do ig=1,ngridmx
          do l=1 ,nlay
c           masse de vapeur d'eau dans la couche l
            mtot(ig)=mtot(ig)+pq(ig,l,nq)*masse(ig,l)
            if(iceparty) then
c           masse de glace d'eau dans la couche l
              icetot(ig)=icetot(ig)+masse(ig,l)*pq(ig,l,nq-1)
c           rayon moyen des cristaux dans la colonne ig
              rave(ig)=rave(ig)+masse(ig,l)*pq(ig,l,nq-1)*rice(ig,l)
              cloudmass(ig,l)=masse(ig,l)*pq(ig,l,nq-1)
            endif
          enddo
          if(iceparty) then
            rave(ig)=rave(ig)/max(icetot(ig),1.e-30)
            if (icetot(ig)*1000.lt.0.01) rave(ig)=0.
          endif
          if (tesfield) then
            if (univtime.eq.lon2pm(ig)) mtot2pm(ig)   = mtot(ig)
            if (univtime.eq.lon2pm(ig).and.iceparty) 
     *        icetot2pm(ig) = icetot(ig)
            if (univtime.eq.lon2pm(ig).and.iceparty) 
     *        rave2pm(ig)   = rave(ig)
            if (univtime.eq.lon2pm(ig)) qsurf2pm(ig)  = qsurf(ig,nqmx)
        endif

c       Outputs total mass of water and clouds    *
c       ******************************************
          IF (tesfield.and.mod(icount,12).EQ.0) THEN
           if (int(lati(ig)*180/pi).gt.0) then
             waterhem(1,1)=waterhem(1,1)+area(ig)*mtot(ig)
             waterhem(1,2)=waterhem(1,2)+area(ig)*qsurf(ig,nqmx)
            if (iceparty)waterhem(1,3)=waterhem(1,3)+area(ig)*icetot(ig)
           endif
           if (int(lati(ig)*180/pi).lt.0) then
             waterhem(2,1)=waterhem(2,1)+area(ig)*mtot(ig)
             waterhem(2,2)=waterhem(2,2)+area(ig)*qsurf(ig,nqmx)
            if (iceparty)waterhem(2,3)=waterhem(2,3)+area(ig)*icetot(ig)
           endif
           if (int(lati(ig)*180/pi).eq.0) then
             waterhem(1,1)=waterhem(1,1)+area(ig)*mtot(ig)/2.
             waterhem(2,1)=waterhem(2,1)+area(ig)*mtot(ig)/2.
             waterhem(1,2)=waterhem(1,2)+area(ig)*qsurf(ig,nqmx)/2.
             waterhem(2,2)=waterhem(2,2)+area(ig)*qsurf(ig,nqmx)/2.
             if (iceparty) then
               waterhem(1,3)=waterhem(1,3)+area(ig)*icetot(ig)/2.
               waterhem(2,3)=waterhem(2,3)+area(ig)*icetot(ig)/2.
             endif
           endif
          ENDIF ! (total mass)
        enddo   ! (ngridmx)

        if (tesfield.and.mod(icount,12).eq.0) then
          write(199,6000)zls*180./pi,waterhem(1,1),waterhem(1,2),
     &    waterhem(1,3),waterhem(2,1),waterhem(2,2),waterhem(2,3)
!6000      format(f8.3,x,6(e12.5,x))
!!!!!!!!!!!!!!!!g95
6000      format(f8.3,1x,6(e12.5,1x))
          write(201,6001) zls*180./pi,pplev(ig_vl2,1),icetot(ig_vl2)
     &    ,zt(ig_vl2,1)
!6001      format(f8.3,x,3(e12.5,x))
!!!!!!!!!!!!!!!!g95
6001      format(f8.3,1x,3(e12.5,1x))
        endif

***************************************************
*       Write outputs at 2PM LT in a dedicated file
***************************************************
        if (tesfield) then
          if (mod(icount,48).eq.0.and.icount.ne.1.and.iceparty) then
            print*,'ecriture dans 2pm_outputs'
            write(200) zls
            write(200) mtot2pm
            if (iceparty) write(200) icetot2pm
            write(200) qsurf2pm
            if (iceparty) write(200) rave2pm
          endif
        endif

c ******************************************    
c      Output in diagfi
c ******************************************    

      if (ngrid.gt.1) then

c          CALL WRITEDIAGFI(ngridmx,'mtot',
c     &  'total mass of water vapor','kg/m2',2,mtot)
c
c         if (callstats) then
c            call wstats(ngrid,"mtot",
c     .   "total mass of water vapor","kg/m2",2,mtot)
c         endif
c
cc        if(.not.iceparty)
cc     &  call WRITEDIAGFI(ngridmx,'taucond',
cc     &  'taux cond H2O ice','kg/kg/s',3,taucond)
cc
c         if(iceparty) then
c            CALL WRITEDIAGFI(ngridmx,'icetot',
c     &    'total mass of water ice','kg/m2',2,icetot)
c           if (callstats) then
c              call wstats(ngrid,"icetot",
c     .        "total mass of water ice","kg/m2",2,icetot)
c           endif
cc          CALL WRITEDIAGFI(ngridmx,'cloudmass',
cc     &    'mass of water ice/layer','kg/m2',3,cloudmass)
cc          CALL WRITEDIAGFI(ngridmx,'rice','ice radius',
cc     &    'meter',3,rice)
cc          CALL WRITEDIAGFI(ngridmx,'rave','Mean ice radius',
cc     &    'meter',2,rave)
c         endif

      else

c          CALL WRITEG1D(ngridmx,1,mtot,'gas','kg/m2')
c          CALL WRITEG1D(ngridmx,1,icetot,'ice','kg/m2')
c          call WRITEG1D(ngridmx,nlay,rice,'rice','um')
c          call WRITEG1D(ngridmx,nlay,dustcores,'pouss','#/layer')

      endif
 99   continue

      RETURN
      END
