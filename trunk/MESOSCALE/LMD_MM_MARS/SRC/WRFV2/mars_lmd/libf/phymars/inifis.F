      SUBROUTINE inifis(ngrid,nlayer,
     $           wday_ini,wdaysec,
     $           wappel_phys,
     $           plat,plon,parea,
     $           prad,pg,pr,pcpp,
     $           nq, wdt,
     $           womeg,wmugaz,
     $           wyear_day,wperiheli,waphelie,wperi_day,wobliquit, 
     $           wz0,wemin_turb,wlmixmin,
     $           wemissiv,wemissiceN,wemissiceS,walbediceN,walbediceS, 
     $           wiceradiusN,wiceradiusS,wdtemisiceN,wdtemisiceS,
     $           walbedodat,winertiedat,wphisfi, 
     $           wzmea,wzstd,wzsig,wzgam,wzthe,
     $           wtheta,wpsi)


      IMPLICIT NONE
c=======================================================================
c
c	CAREFUL: THIS IS A VERSION TO BE USED WITH WRF !!!
c
c	... CHECK THE ****WRF lines
c
c=======================================================================
c
c   subject:
c   --------
c
c   Initialisation for the physical parametrisations of the LMD 
c   martian atmospheric general circulation modele.
c
c   author: Frederic Hourdin 15 / 10 /93
c   -------
c   modified: Sebastien Lebonnois 11/06/2003 (new callphys.def)
c   adapted to the WRF use - Aymeric Spiga - Jan 2007
c
c   arguments:
c   ----------
c
c   input:
c   ------
c
c    ngrid                 Size of the horizontal grid.
c                          All internal loops are performed on that grid.
c    nlayer                Number of vertical layers.
c    pdayref               Day of reference for the simulation
c    firstcall             True at the first call
c    lastcall              True at the last call
c    pday                  Number of days counted from the North. Spring
c                          equinoxe.
c
c=======================================================================
c
c-----------------------------------------------------------------------
c   declarations:
c   -------------
 
#include "dimensions.h"
#include "dimphys.h"
#include "planete.h"
#include "comcstfi.h"
#include "comsaison.h"
#include "comdiurn.h"
#include "comgeomfi.h"
#include "callkeys.h"
#include "surfdat.h"
#include "slope.h"


      INTEGER ngrid,nlayer,nq

      REAL prad,pg,pr,pcpp,pdaysec
      REAL womeg,wmugaz,wdaysec  
      REAL wyear_day,wperiheli,waphelie,wperi_day,wobliquit
      REAL wz0,wemin_turb,wlmixmin
      REAL wemissiv,wemissiceN,wemissiceS,walbediceN,walbediceS

      REAL wiceradiusN,wiceradiusS,wdtemisiceN,wdtemisiceS     
      REAL walbedodat(ngrid),winertiedat(ngrid),wphisfi(ngrid) 
      REAL wzmea(ngrid),wzstd(ngrid),wzsig(ngrid)
      REAL wzgam(ngrid),wzthe(ngrid)
      REAL wtheta(ngrid),wpsi(ngrid)

      REAL plat(ngrid),plon(ngrid),parea(ngridmx)
      integer wday_ini
      REAL wdt
      INTEGER ig,ierr

      INTEGER wecri_phys
      REAL wappel_phys
 
      EXTERNAL iniorbit,orbite
      EXTERNAL SSUM
      REAL SSUM
 
      CHARACTER ch1*12
      CHARACTER ch80*80

      logical chem, h2o

c ****WRF
c
c------------------------------------------------------
c  Fill some parameters in the 'include' files
c  >> Do part of the job previously done by phyetat0.F
c  >> Complete list of parameters is found in tabfi.F
c------------------------------------------------------
c
c Values are defined in the module_model_constants.F WRF routine
c      
      ! in 'comcstfi.h'
      rad=prad                  
      cpp=pcpp
      g=pg
      r=pr                    
      rcp=r/cpp
      daysec=wdaysec  
      omeg=womeg                
      mugaz=wmugaz  
      print*,"check: rad,cpp,g,r,rcp,daysec,omeg,mugaz"
      print*,rad,cpp,g,r,rcp,daysec,omeg,mugaz
    
      ! in 'planet.h' 
      year_day=wyear_day
      periheli=wperiheli
      aphelie=waphelie
      peri_day=wperi_day
      obliquit=wobliquit
      z0=wz0
      emin_turb=wemin_turb
      lmixmin=wlmixmin
      print*,"check: year_day,periheli,aphelie,peri_day,obliquit"
      print*,year_day,periheli,aphelie,peri_day,obliquit
      print*,"check: z0,emin_turb,lmixmin"
      print*,z0,emin_turb,lmixmin

      ! in 'surfdat.h'
      emissiv=wemissiv
      emisice(1)=wemissiceN
      emisice(2)=wemissiceS
      albedice(1)=walbediceN
      albedice(2)=walbediceS
      iceradius(1)=wiceradiusN
      iceradius(2)=wiceradiusS
      dtemisice(1)=wdtemisiceN
      dtemisice(2)=wdtemisiceS
      print*,"check: emissiv,emisice,albedice,iceradius,dtemisice"
      print*,emissiv,emisice,albedice,iceradius,dtemisice

c
c Values are defined in the WPS processing
c  
        albedodat(:)=walbedodat(:)
        inertiedat(:)=winertiedat(:)
        phisfi(:)=wphisfi(:)
        print*,"check: albedodat(1),inertiedat(1),phisfi(1)"
        print*,albedodat(1),inertiedat(1),phisfi(1)
        print*,"check: albedodat(end),inertiedat(end),phisfi(end)"
        print*,albedodat(ngrid),inertiedat(ngrid),phisfi(ngrid)

        ! NB: usually, gravity wave scheme is useless in mesoscale modeling
        ! NB: we however keep the option for coarse grid case ... 	
        zmea(:)=wzmea(:)
        zstd(:)=wzstd(:)
        zsig(:)=wzsig(:)
        zgam(:)=wzgam(:)
        zthe(:)=wzthe(:)
        print*,"check: gw param"
        print*,zmea(1),zmea(ngrid)
        print*,zstd(1),zstd(ngrid)
        print*,zsig(1),zsig(ngrid)
        print*,zgam(1),zgam(ngrid)
        print*,zthe(1),zthe(ngrid)

        !
        ! in slope.h
        !
        theta_sl(:)=wtheta(:)
        psi_sl(:)=wpsi(:)
        print*,"check: theta_sl(1),psi_sl(1)"
        print*,theta_sl(1),psi_sl(1)
        print*,"check: theta_sl(end),psi_sl(end)"
        print*,theta_sl(ngrid),psi_sl(ngrid)


c ****WRF

 

c --------------------------------------------------------
c     The usual Tests
c     --------------

c ****WRF
c useless here, because it was already done in the WRF driver

      IF (nlayer.NE.nlayermx) THEN
         PRINT*,'STOP in inifis'
         PRINT*,'Probleme de dimensions :'
         PRINT*,'nlayer     = ',nlayer
         PRINT*,'nlayermx   = ',nlayermx
         STOP
      ENDIF

      IF (ngrid.NE.ngridmx) THEN
         PRINT*,'STOP in inifis'
         PRINT*,'Probleme de dimensions :'
         PRINT*,'ngrid     = ',ngrid
         PRINT*,'ngridmx   = ',ngridmx
         STOP
      ENDIF


c --------------------------------------------------------------
c  Reading the "callphys.def" file controlling some key options
c --------------------------------------------------------------

      callrad=.true.
      calldifv=.true.
      calladj=.true.
      callcond=.true.
      callsoil=.true.
      season=.true.
      diurnal=.false.
      lwrite=.false.
      calllott=.true.
      iaervar=2
      iddist=3
      topdustref=55.
      OPEN(99,file='callphys.def',status='old',form='formatted'
     .     ,iostat=ierr)
      IF(ierr.EQ.0) THEN
         !PRINT*
         !PRINT*
         PRINT*,'--------------------------------------------'
         PRINT*,' Parametres pour la physique (callphys.def)'
         PRINT*,'--------------------------------------------'

         READ(99,*)
         READ(99,*)

         READ(99,fmt='(a)') ch1 
         READ(99,*) tracer
         WRITE(*,8000) ch1,tracer

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') diurnal
         WRITE(*,8000) ch1,diurnal

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') season
         WRITE(*,8000) ch1,season

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') lwrite
         WRITE(*,8000) ch1,lwrite

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') callstats
         WRITE(*,8000) ch1,callstats

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') calleofdump
         WRITE(*,8000) ch1,calleofdump

         READ(99,*)
         READ(99,*)

         READ(99,fmt='(a)') ch1
         READ(99,*,iostat=ierr) iaervar
         if(ierr.ne.0) stop'no iaervar in callphys.def (old?)'
c****WRF: ligne trop longue ....
         WRITE(*,8001) ch1,iaervar

         READ(99,fmt='(a)') ch1
         READ(99,*) iddist
         WRITE(*,8001) ch1,iddist

         READ(99,fmt='(a)') ch1
         READ(99,*) topdustref
         WRITE(*,8002) ch1,topdustref

         READ(99,*)
         READ(99,*)

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') callrad
         WRITE(*,8000) ch1,callrad

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') callnlte
         WRITE(*,8000) ch1,callnlte
         
         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') callnirco2
         WRITE(*,8000) ch1,callnirco2

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') calldifv
         WRITE(*,8000) ch1,calldifv

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') calladj
         WRITE(*,8000) ch1,calladj

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') callcond
         WRITE(*,8000) ch1,callcond

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') callsoil
         WRITE(*,8000) ch1,callsoil

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') calllott
         WRITE(*,8000) ch1,calllott

         READ(99,*) 
         READ(99,*)

         READ(99,fmt='(a)') ch1
         READ(99,*) iradia
         WRITE(*,8001) ch1,iradia

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') callg2d
         WRITE(*,8000) ch1,callg2d

         READ(99,fmt='(a)') ch1
         READ(99,*) rayleigh
         WRITE(*,8000) ch1,rayleigh

         READ(99,*) 
         READ(99,*)

c TRACERS:

         READ(99,fmt='(a)') ch1
         READ(99,*) dustbin
         WRITE(*,8001) ch1,dustbin

         READ(99,fmt='(a)') ch1
         READ(99,*) active
         WRITE(*,8000) ch1,active

c Test of incompatibility:
c if active is used, then dustbin should be > 0

         if (active.and.(dustbin.lt.1)) then
           print*,'if active is used, then dustbin should > 0'
           stop
         endif

         READ(99,fmt='(a)') ch1
         READ(99,*) doubleq
         WRITE(*,8000) ch1,doubleq

c Test of incompatibility:
c if doubleq is used, then dustbin should be 1

         if (doubleq.and.(dustbin.ne.1)) then
           print*,'if doubleq is used, then dustbin should be 1'
           stop
         endif

         READ(99,fmt='(a)') ch1
         READ(99,*) lifting
         WRITE(*,8000) ch1,lifting

c Test of incompatibility:
c if lifting is used, then dustbin should be > 0

         if (lifting.and.(dustbin.lt.1)) then
           print*,'if lifting is used, then dustbin should > 0'
           stop
         endif

         READ(99,fmt='(a)') ch1
         READ(99,*) callddevil
         WRITE(*,8000) ch1,callddevil

c Test of incompatibility:
c if dustdevil is used, then dustbin should be > 0

         if (callddevil.and.(dustbin.lt.1)) then
           print*,'if dustdevil is used, then dustbin should > 0'
           stop
         endif

         READ(99,fmt='(a)') ch1
         READ(99,*) scavenging
         WRITE(*,8000) ch1,scavenging

c Test of incompatibility:
c if scavenging is used, then dustbin should be > 0

         if (scavenging.and.(dustbin.lt.1)) then
           print*,'if scavenging is used, then dustbin should > 0'
           stop
         endif

         READ(99,fmt='(a)') ch1
         READ(99,*) sedimentation
         WRITE(*,8000) ch1,sedimentation

         READ(99,fmt='(a)') ch1
         READ(99,*) iceparty
         WRITE(*,8000) ch1,iceparty

         READ(99,fmt='(a)') ch1
         READ(99,*) activice
         WRITE(*,8000) ch1,activice

c Test of incompatibility:
c if activice is used, then iceparty should be used too

         if (activice.and..not.iceparty) then
           print*,'if activice is used, iceparty should be used too'
           stop
         endif

         READ(99,fmt='(a)') ch1
         READ(99,*) water
         WRITE(*,8000) ch1,water

c Test of incompatibility:
c if iceparty is used, then water should be used too

         if (.not.water) then
            iceparty = .false. 
         endif

         READ(99,fmt='(a)') ch1
         READ(99,*) caps
         WRITE(*,8000) ch1,caps

         READ(99,fmt='(a)') ch1
         READ(99,*) photochem
         WRITE(*,8000) ch1,photochem

         READ(99,*)
         READ(99,*)

c THERMOSPHERE

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') callthermos
         WRITE(*,8000) ch1,callthermos

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') thermoswater
         WRITE(*,8000) ch1,thermoswater

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') callconduct
         WRITE(*,8000) ch1,callconduct

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') calleuv
         WRITE(*,8000) ch1,calleuv

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') callmolvis
         WRITE(*,8000) ch1,callmolvis

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') callmoldiff
         WRITE(*,8000) ch1,callmoldiff

         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') thermochem
         WRITE(*,8000) ch1,thermochem

         READ(99,fmt='(a)') ch1
         READ(99,*) solarcondate
         WRITE(*,*) ch1,solarcondate

         if (.not.callthermos) then
                 thermoswater=.false.          
                 callconduct=.false.          
                 calleuv=.false.          
                 callmolvis=.false.          
                 callmoldiff=.false.          
                 thermochem=.false.          
        end if

c Test of incompatibility:
c if photochem is used, then water should be used too

         if (photochem.and..not.water) then
           print*,'if photochem is used, water should be used too'
           stop
         endif

c if callthermos is used, then thermoswater should be used too 
c (if water not used already)

         if (callthermos .and. .not.water) then
           if (callthermos .and. .not.thermoswater) then
             print*,'if callthermos is used, water or thermoswater 
     &               should be used too'
             stop
           endif
         endif

         PRINT*,'--------------------------------------------'
         PRINT*
         PRINT*
      ELSE
         write(*,*)
         write(*,*) 'Cannot read file callphys.def. Is it here ?'
         stop
      ENDIF
      CLOSE(99)

      !!-----------------------------------------
      !! titus cap
      OPEN(99,file='titus.def',status='old',form='formatted'
     .     ,iostat=ierr)
      IF(ierr.EQ.0) THEN
         READ(99,*,iostat=ierr) tituscap
      ELSE
         write(*,*) 'no titus.def files. assuming no Titus cap.'
         tituscap = .false.
      ENDIF
      write(*,*) 'Titus cap is: ',tituscap
      CLOSE(99)
      !!-----------------------------------------


8000  FORMAT(t5,a12,l8)
8001  FORMAT(t5,a12,i8)
8002  FORMAT(t5,a12,f8.1)


c ****WRF
c*****************************************************
c Since it comes from WRF settings, we have to
c fill dtphys in the include file 
c It must be set now, because it is used afterwards
c
c Opportunity is taken to fill ecri_phys as well
c*****************************************************
        dtphys=wdt*wappel_phys
        print*,'Physical timestep (s) ',dtphys
        ecri_phys=10 !8.e18  !! a dummy low frequency
        print*,'Physical frequency for writing ',ecri_phys
c
c ****WRF


      PRINT*
      PRINT*,'daysec',daysec
      PRINT*
      PRINT*,'The radiative transfer is computed:'
      PRINT*,'           each ',iradia,' physical time-step'
      PRINT*,'        or each ',iradia*dtphys,' seconds'
      PRINT*



c --------------------------------------------------------------
c  Managing the Longwave radiative transfer
c --------------------------------------------------------------

c     In most cases, the run just use the following values :
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      callemis=.true.     
c     ilwd=10*int(daysec/dtphys) ! bug before 22/10/01       
      ilwd=10*int(daysec/dtphys) 
      ilwn=2               
      linear=.true.        
      ncouche=3
      alphan=0.4
      ilwb=2
      semi=0

c     BUT people working hard on the LW may want to read them in 'radia.def' 
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      OPEN(99,file='radia.def',status='old',form='formatted'
     .     ,iostat=ierr)
      IF(ierr.EQ.0) THEN
         write(*,*) 'Reading radia.def !!!'
         READ(99,fmt='(a)') ch1
         READ(99,*) callemis
         WRITE(*,8000) ch1,callemis

         READ(99,fmt='(a)') ch1
         READ(99,*) iradia
         WRITE(*,8001) ch1,iradia

         READ(99,fmt='(a)') ch1
         READ(99,*) ilwd
         WRITE(*,8001) ch1,ilwd

         READ(99,fmt='(a)') ch1
         READ(99,*) ilwn
         WRITE(*,8001) ch1,ilwn

         READ(99,fmt='(a)') ch1
         READ(99,*) linear
         WRITE(*,8000) ch1,linear

         READ(99,fmt='(a)') ch1
         READ(99,*) ncouche
         WRITE(*,8001) ch1,ncouche

         READ(99,fmt='(a)') ch1
         READ(99,*) alphan
         WRITE(*,*) ch1,alphan

         READ(99,fmt='(a)') ch1
         READ(99,*) ilwb
         WRITE(*,8001) ch1,ilwb


         READ(99,fmt='(a)') ch1
         READ(99,'(l1)') callg2d
         WRITE(*,8000) ch1,callg2d

         READ(99,fmt='(a)') ch1
         READ(99,*) semi
         WRITE(*,*) ch1,semi
      end if

c-----------------------------------------------------------------------
c     Some more initialization:
c     ------------------------

      ! in 'comgeomfi.h'	
      CALL SCOPY(ngrid,plon,1,long,1)
      CALL SCOPY(ngrid,plat,1,lati,1)
      CALL SCOPY(ngrid,parea,1,area,1)
      totarea=SSUM(ngridmx,area,1)

      ! in 'comdiurn.h'
      DO ig=1,ngrid
         sinlat(ig)=sin(plat(ig))
         coslat(ig)=cos(plat(ig))
         sinlon(ig)=sin(plon(ig))
         coslon(ig)=cos(plon(ig))
      ENDDO

      pi=2.*asin(1.)

c     managing the tracers, and tests:
c     -------------------------------

      if(tracer) then

c          when photochem is used, nqchem_min is the rank
c          of the first chemical species

       nqchem_min = 1
       if (photochem .or. callthermos) then
         chem = .true.
        if (doubleq) then
          nqchem_min = 3
        else
          nqchem_min = dustbin+1
        end if
       end if

       if (water .or. thermoswater) h2o = .true.

c          TESTS

       print*,'TRACERS:'

       if ((doubleq).and.(h2o).and.
     $     (chem).and.(iceparty)) then
         print*,' 1: dust ; 2: dust (doubleq)'
         print*,' 3 to ',nqmx-2,': chemistry'
         print*,nqmx-1,': water ice ; ',nqmx,': water vapor'
       endif

       if ((doubleq).and.(h2o).and.
     $     (chem).and..not.(iceparty)) then
         print*,' 1: dust ; 2: dust (doubleq)'
         print*,' 3 to ',nqmx-1,': chemistry'
         print*,nqmx,': water vapor'
       endif

       if ((doubleq).and.(h2o).and.
     $     .not.(chem).and.(iceparty)) then
         print*,' 1: dust ; 2: dust (doubleq)'
         print*,nqmx-1,': water ice ; ',nqmx,': water vapor'
         if (nqmx.ne.4) then
           print*,'nqmx should be 4 with these options.'
           print*,'(or check callphys.def)'
           stop
         endif
       endif

       if ((doubleq).and.(h2o).and.
     $     .not.(chem).and..not.(iceparty)) then
         print*,' 1: dust ; 2: dust (doubleq)'
         print*,nqmx,': water vapor'
         if (nqmx.ne.3) then
           print*,'nqmx should be 3 with these options...'
           print*,'(or check callphys.def)'
           stop
         endif
       endif

       if ((doubleq).and..not.(h2o)) then
         print*,' 1: dust ; 2: dust (doubleq)'
         if (nqmx.ne.2) then
           print*,'nqmx should be 2 with these options...'
           print*,'(or check callphys.def)'
           stop
         endif
       endif

       if (.not.(doubleq).and.(h2o).and.
     $     (chem).and.(iceparty)) then
         if (dustbin.gt.0) then
           print*,' 1 to ',dustbin,': dust bins'
         endif
         print*,nqchem_min,' to ',nqmx-2,': chemistry'
         print*,nqmx-1,': water ice ; ',nqmx,': water vapor'
       endif
       if (.not.(doubleq).and.(h2o).and.
     $     (chem).and..not.(iceparty)) then
         if (dustbin.gt.0) then
           print*,' 1 to ',dustbin,': dust bins'
         endif
         print*,nqchem_min,' to ',nqmx-1,': chemistry'
         print*,nqmx,': water vapor'
       endif
       if (.not.(doubleq).and.(h2o).and.
     $     .not.(chem).and.(iceparty)) then
         if (dustbin.gt.0) then
           print*,' 1 to ',dustbin,': dust bins'
         endif
         print*,nqmx-1,': water ice ; ',nqmx,': water vapor'
         if (nqmx.ne.(dustbin+2)) then
           print*,'nqmx should be ',(dustbin+2),
     $            ' with these options...'
           print*,'(or check callphys.def)'
           stop
         endif
       endif
       if (.not.(doubleq).and.(h2o).and.
     $     .not.(chem).and..not.(iceparty)) then
         if (dustbin.gt.0) then
           print*,' 1 to ',dustbin,': dust bins'
         endif
         print*,nqmx,': water vapor'
         if (nqmx.ne.(dustbin+1)) then
           print*,'nqmx should be ',(dustbin+1),
     $            ' with these options...'
           print*,'(or check callphys.def)'
           stop
         endif
       endif
       if (.not.(doubleq).and..not.(h2o)) then
         if (dustbin.gt.0) then
           print*,' 1 to ',dustbin,': dust bins'
           if (nqmx.ne.dustbin) then
             print*,'nqmx should be ',dustbin,
     $              ' with these options...'
             print*,'(or check callphys.def)'
             stop
           endif
         else
           print*,'dustbin=',dustbin,
     $            ': tracer should be F with these options...'
     $           ,'UNLESS you just want to move tracers around '
         endif
       endif

      endif

      RETURN
      END
