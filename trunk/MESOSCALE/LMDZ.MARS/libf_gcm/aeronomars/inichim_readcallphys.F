      SUBROUTINE inichim_readcallphys
      IMPLICIT NONE
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
c
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

      character*12 ch1
      integer ierr
      logical chem, h2o


c --------------------------------------------------------------
c  Reading the "callphys.def" file controlling some key options
c --------------------------------------------------------------

      OPEN(99,file='callphys.def',status='old',form='formatted'
     .     ,iostat=ierr)
      IF(ierr.EQ.0) THEN
         PRINT*
         PRINT*
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
         if(ierr.ne.0) stop'Can t read iaervar in callphys.def (old?)'
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

         if (iceparty.and..not.water) then
           print*,'if iceparty is used, then water should be used too'
           stop
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

      pi=2.*asin(1.)

c     managing the tracers, and tests:
c     -------------------------------

      if(tracer) then

c          when photochem is used, nqchem_min is the rank
c          of the first chemical species

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
           print*,'nqmx should be 4 with these options...'
           stop
         endif
       endif
       if ((doubleq).and.(h2o).and.
     $     .not.(chem).and..not.(iceparty)) then
         print*,' 1: dust ; 2: dust (doubleq)'
         print*,nqmx,': water vapor'
         if (nqmx.ne.3) then
           print*,'nqmx should be 3 with these options...'
           stop
         endif
       endif

       if ((doubleq).and..not.(h2o)) then
         print*,' 1: dust ; 2: dust (doubleq)'
         if (nqmx.ne.2) then
           print*,'nqmx should be 2 with these options...'
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
           stop
         endif
       endif
       if (.not.(doubleq).and..not.(h2o)) then
         if (dustbin.gt.0) then
           print*,' 1 to ',dustbin,': dust bins'
           if (nqmx.ne.dustbin) then
             print*,'nqmx should be ',dustbin,
     $              ' with these options...'
             stop
           endif
         else
           print*,'dustbin=',dustbin,
     $            ': tracer should be F with these options...'
           stop
         endif
       endif

      endif

8000  FORMAT(t5,a12,l8)
8001  FORMAT(t5,a12,i8)
8002  FORMAT(t5,a12,f8.1)

      RETURN
      END
