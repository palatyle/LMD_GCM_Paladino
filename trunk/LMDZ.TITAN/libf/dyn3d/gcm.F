      PROGRAM gcm

      use ioipsl_getincom, only: getin
      use infotrac, only: iniadvtrac, nqtot, iadv
      use sponge_mod,only: callsponge,mode_sponge,sponge
      use control_mod, only: nday, day_step, iperiod, iphysiq,
     &                       iconser, dissip_period
!      use comgeomphy, only: initcomgeomphy
      USE mod_const_mpi, ONLY: COMM_LMDZ
      use filtreg_mod, only: inifilr
      USE comvert_mod, ONLY: ap,bp
      USE comconst_mod, ONLY: daysec,dtvr,dtphys,dtdiss,rad,g,r,cpp
      USE logic_mod, ONLY: ecripar,forward,leapf,apphys,statcl,conser,
     .		apdiss,purmats,iflag_phys,apphys
      USE temps_mod, ONLY: day_ini,day_end,itaufin,dt
      USE iniphysiq_mod, ONLY: iniphysiq
      IMPLICIT NONE

c      ......   Version  du 10/01/98    ..........

c             avec  coordonnees  verticales hybrides 
c   avec nouveaux operat. dissipation * ( gradiv2,divgrad2,nxgraro2 )

c=======================================================================
c
c   Auteur:  P. Le Van /L. Fairhead/F.Hourdin
c   -------
c
c   Objet:
c   ------
c
c   GCM LMD nouvelle grille
c
c=======================================================================
c
c  ... Dans inigeom , nouveaux calculs pour les elongations  cu , cv
c      et possibilite d'appeler une fonction f(y)  a derivee tangente
c      hyperbolique a la  place de la fonction a derivee sinusoidale.

c  ... Possibilite de choisir le shema de Van-leer pour l'advection de
c        q  , en faisant iadv = 3  dans   traceur  (29/04/97) .
c
c      Pour Van-Leer + Vapeur d'eau saturee, iadv(1)=4. (F.Codron,10/99)
c
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "dimensions.h"
#include "paramet.h"
#include "comdissnew.h"
#include "comgeom.h"
!#include "control.h"
#include "netcdf.inc"
#include "tracstoke.h"
!#include"advtrac.h"

      INTEGER*4  iday ! jour julien
      REAL       time ! Heure de la journee en fraction d''1 jour
      REAL zdtvr

c   variables dynamiques
      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm) ! vents covariants
      real, dimension(ip1jmp1,llm) :: teta   ! temperature potentielle 
      REAL,allocatable :: q(:,:,:)           ! champs advectes
      REAL ps(ip1jmp1)                       ! pression  au sol
      REAL pext(ip1jmp1)                     ! pression  extensive
      REAL p (ip1jmp1,llmp1  )               ! pression aux interfac.des couches
      REAL pks(ip1jmp1)                      ! exner au  sol
      REAL pk(ip1jmp1,llm)                   ! exner au milieu des couches
      REAL pkf(ip1jmp1,llm)                  ! exner filt.au milieu des couches
      REAL masse(ip1jmp1,llm)                ! masse d''air
      REAL phis(ip1jmp1)                     ! geopotentiel au sol
      REAL phi(ip1jmp1,llm)                  ! geopotentiel
      REAL w(ip1jmp1,llm)                    ! vitesse verticale

c variables dynamiques intermediaire pour le transport
      REAL pbaru(ip1jmp1,llm),pbarv(ip1jm,llm) !flux de masse

c   variables dynamiques au pas -1
      REAL vcovm1(ip1jm,llm),ucovm1(ip1jmp1,llm)

      REAL tetam1(ip1jmp1,llm),psm1(ip1jmp1)
      REAL massem1(ip1jmp1,llm)

c   tendances dynamiques
      REAL dv(ip1jm,llm),du(ip1jmp1,llm)
      REAL dteta(ip1jmp1,llm),dp(ip1jmp1)
      REAL,ALLOCATABLE :: dq(:,:,:)

c   tendances de la dissipation
      REAL dvdis(ip1jm,llm),dudis(ip1jmp1,llm)
      REAL dhdis(ip1jmp1,llm)

c   tendances physiques
      REAL dvfi(ip1jm,llm),dufi(ip1jmp1,llm)
      REAL dhfi(ip1jmp1,llm),dpfi(ip1jmp1)
      REAL,ALLOCATABLE :: dqfi(:,:,:)

c   variables pour le fichier histoire
      REAL dtav      ! intervalle de temps elementaire

      REAL tppn(iim),tpps(iim),tpn,tps
c
!      INTEGER iadv(nqmx) ! indice schema de transport pour le traceur iq

      INTEGER itau,itaufinp1,iav


      EXTERNAL caldyn, traceur
      EXTERNAL dissip,geopot,iniconst
      EXTERNAL integrd,SCOPY
      EXTERNAL inigeom
      EXTERNAL exner_hyb
      EXTERNAL defrun_new, test_period
      REAL  SSUM
      REAL time_0 , finvmaold(ip1jmp1,llm)

      LOGICAL lafin
      INTEGER ij,iq,l,ierr,numvanle,iapp_tracvl

      REAL rdayvrai,rdaym_ini,rday_ecri
!      LOGICAL first
      REAL beta(ip1jmp1,llm)

      LOGICAL offline  ! Controle du stockage ds "fluxmass"
      PARAMETER (offline=.false.)

      character*20 modname
      character*80 abort_message

! flag to run with or without tracer transport (read from run.def)
      logical tracer
      
C Calendrier
      LOGICAL true_calendar
      PARAMETER (true_calendar = .false.)

! flag to set/remove calls to groupeun
      logical callgroupeun
      parameter (callgroupeun = .false.)

!     added by RDW for tests without dynamics
      logical calldyn
      parameter (calldyn = .true.)

!     added by RW for test
!      real pmean,airetot

!     added by FF for dissipation / energy conservation tests
      logical dissip_conservative
      parameter (dissip_conservative = .false.)
      REAL ecin(ip1jmp1,llm),ecin0(ip1jmp1,llm)
!     d (theta)/ d t (pot. temperature) due to the creation 
!     of thermal energy by dissipation
      REAL dtetaecdt(ip1jmp1,llm)
      REAL vcont(ip1jm,llm),ucont(ip1jmp1,llm)

c-----------------------------------------------------------------------
c    variables pour l'initialisation de la physique :
c    ------------------------------------------------
!      INTEGER ngridmx
!      PARAMETER( ngridmx = 2+(jjm-1)*iim - 1/jjm   )
!      REAL zcufi(ngridmx),zcvfi(ngridmx)
!      REAL latfi(ngridmx),lonfi(ngridmx)
!      REAL airefi(ngridmx)
!      SAVE latfi, lonfi, airefi
!      INTEGER i,j

c-----------------------------------------------------------------------
c   Initialisations:
c   ----------------

      modname = 'gcm'
      lafin    = .FALSE.

c-----------------------------------------------------------------------
      CALL defrun_new( 99, .TRUE. )
      WRITE(*,*) ""
      WRITE(*,*) "Run with or without tracer transport ?"
      tracer=.true. ! default value
      call getin("tracer",tracer)
      WRITE(*,*)" tracer = ",tracer

! Initialize tracers
      CALL iniadvtrac(nqtot,numvanle)
! Allocation de la tableau q : champs advectes   
      allocate(q(ip1jmp1,llm,nqtot))
      allocate(dq(ip1jmp1,llm,nqtot))
      allocate(dqfi(ip1jmp1,llm,nqtot))

      CALL dynetat0("start.nc",vcov,ucov,
     .              teta,q,masse,ps,phis, time_0)

c  on recalcule eventuellement le pas de temps

      IF(MOD(day_step,iperiod).NE.0)
     * STOP'Il faut choisir un nb de pas par jour multiple de iperiod'

      IF(MOD(day_step,iphysiq).NE.0)
     * STOP'Il faut choisir un nb de pas par jour multiple de iphysiq'

      zdtvr    = daysec/FLOAT(day_step)

        IF(dtvr.NE.zdtvr) THEN
         PRINT*,'WARNING!!! changement de pas de temps',dtvr,'>',zdtvr
        ENDIF

c  nombre d'etats dans les fichiers demarrage et histoire

      dtvr = zdtvr
      CALL iniconst
      CALL inigeom
      CALL inifilr

c
c   ......    P.Le Van    ( modif  le 29/04/97 )   .........  
c

      CALL inidissip ( lstardis, nitergdiv, nitergrot, niterh   ,
     *                tetagdiv, tetagrot , tetatemp              )
c

c-----------------------------------------------------------------------
c   Initialisation de la physique :
c   -------------------------------

!      IF (call_iniphys.and.(iflag_phys==1.or.iflag_phys>=100)) THEN
!         latfi(1)=rlatu(1)
!         lonfi(1)=0.
!         zcufi(1) = cu(1)
!         zcvfi(1) = cv(1)
!         DO j=2,jjm
!            DO i=1,iim
!               latfi((j-2)*iim+1+i)= rlatu(j)
!               lonfi((j-2)*iim+1+i)= rlonv(i)
!               zcufi((j-2)*iim+1+i) = cu((j-1)*iip1+i)
!               zcvfi((j-2)*iim+1+i) = cv((j-1)*iip1+i)
!            ENDDO
!         ENDDO
!         latfi(ngridmx)= rlatu(jjp1)
!         lonfi(ngridmx)= 0.
!         zcufi(ngridmx) = cu(ip1jm+1)
!         zcvfi(ngridmx) = cv(ip1jm-iim)
!
!         ! build airefi(), mesh area on physics grid
!         CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,aire,airefi)
!         ! Poles are single points on physics grid
!         airefi(1)=airefi(1)*iim
!         airefi(ngridmx)=airefi(ngridmx)*iim

! Initialisation de la physique: pose probleme quand on tourne 
! SANS physique, car iniphysiq.F est dans le repertoire phy[]... 
! Il faut une cle CPP_PHYS
!#ifdef CPP_PHYS
!         CALL iniphysiq(iim,jjm,llm,daysec,day_ini,dtphys/nsplit_phys,
         CALL iniphysiq(iim,jjm,llm,
     &                  (jjm-1)*iim+2,comm_lmdz,
     &                  daysec,day_ini,dtphys,
     &                  rlatu,rlatv,rlonu,rlonv,
     &                  aire,cu,cv,rad,g,r,cpp,
     &                iflag_phys)
!#endif
!         call_iniphys=.false.
!      ENDIF ! of IF (call_iniphys.and.(iflag_phys.eq.1))

      CALL pression ( ip1jmp1, ap, bp, ps, p       )

      call dump2d(iip1,jjp1,ps,'PRESSION SURFACE')
      CALL exner_hyb( ip1jmp1, ps, p,beta, pks, pk, pkf )


c  numero de stockage pour les fichiers de redemarrage:

c-----------------------------------------------------------------------
c   temps de depart et de fin:
c   --------------------------

      itau = 0
      iday = day_ini+itau/day_step
      time = FLOAT(itau-(iday-day_ini)*day_step)/day_step+time_0


         IF(time.GT.1.) THEN
          time = time-1.
          iday = iday+1
         ENDIF
      itaufin   = nday*day_step
c ********************************
c      itaufin = 120   ! temporaire !!
c ********************************
      itaufinp1 = itaufin +1

      day_end = day_ini + nday
      PRINT 300, itau,itaufin,day_ini,day_end
 300  FORMAT('1'/,15x,'run du pas',i7,2x,'au pas',i7,2x,  
     . 'c''est a dire du jour',i7,3x,'au jour',i7//)

      CALL dynredem0("restart.nc",day_end,phis)

      ecripar = .TRUE.

      dtav = iperiod*dtvr/daysec


c   Quelques initialisations pour les traceurs
      call initial0(ijp1llm*nqtot,dq)
c     istdyn=day_step/4     ! stockage toutes les 6h=1jour/4
c     istphy=istdyn/iphysiq

      write(*,*) "gcm: callgroupeun set to:",callgroupeun
c-----------------------------------------------------------------------
c   Debut de l'integration temporelle:
c   ----------------------------------

   1  CONTINUE
c
      IF( MOD( itau, 10* day_step ).EQ.0 )  THEN
        CALL test_period ( ucov,vcov,teta,q,p,phis )
        PRINT *,' ----   Test_period apres continue   OK ! -----', itau
      ENDIF 

      if (callgroupeun) then
        call groupeun(jjp1,llm,ucov,.true.)
        call groupeun(jjm,llm,vcov,.true.)
        call groupeun(jjp1,llm,teta,.true.)
        call groupeun(jjp1,llm,masse,.true.)
        call groupeun(jjp1,1,ps,.false.)
      endif

      CALL SCOPY( ijmllm ,vcov , 1, vcovm1 , 1 )
      CALL SCOPY( ijp1llm,ucov , 1, ucovm1 , 1 )
      CALL SCOPY( ijp1llm,teta , 1, tetam1 , 1 )
      CALL SCOPY( ijp1llm,masse, 1, massem1, 1 )
      CALL SCOPY( ip1jmp1, ps  , 1,   psm1 , 1 )

      forward = .TRUE.
      leapf   = .FALSE.
      dt      =  dtvr

c   ...    P.Le Van .26/04/94  ....

      CALL SCOPY   ( ijp1llm,   masse, 1, finvmaold,     1 )
      CALL filtreg ( finvmaold ,jjp1, llm, -2,2, .TRUE., 1 )


   2  CONTINUE

c-----------------------------------------------------------------------

c   date:
c   -----

!      write(*,*) 'GCM: itau=',itau

c   gestion des appels de la physique et des dissipations:
c   ------------------------------------------------------
c
c   ...    P.Le Van  ( 6/02/95 )  ....

      apphys = .FALSE.
      statcl = .FALSE.
      conser = .FALSE.
      apdiss = .FALSE.

      IF( purmats ) THEN
         IF( MOD(itau,iconser) .EQ.0.AND.  forward    ) conser = .TRUE.
         IF( MOD(itau,dissip_period ).EQ.0
     &                              .AND..NOT.forward ) apdiss = .TRUE.
         IF( MOD(itau,iphysiq ).EQ.0.AND..NOT.forward 
     $                        .AND. (iflag_phys.eq.1) ) apphys = .TRUE.
      ELSE
         IF( MOD(itau   ,iconser) .EQ. 0              ) conser = .TRUE.
         IF( MOD(itau+1,dissip_period)  .EQ. 0        ) apdiss = .TRUE.
         IF( MOD(itau+1,iphysiq).EQ.0
     &                     .AND. (iflag_phys.eq.1)    ) apphys = .TRUE.
      END IF

c-----------------------------------------------------------------------
c   calcul des tendances dynamiques:
c   --------------------------------

      CALL geopot  ( ip1jmp1, teta  , pk , pks,  phis  , phi   )

      if(calldyn)then
         CALL caldyn 
     $        ( itau,ucov,vcov,teta,ps,masse,pk,pkf,phis ,
     $        phi,conser,du,dv,dteta,dp,w,pbaru,pbarv,time+iday-day_ini)
      endif

c-----------------------------------------------------------------------
c   calcul des tendances advection des traceurs (dont l'humidite)
c   -------------------------------------------------------------



      if (tracer) then
       IF( forward. OR . leapf )  THEN

        DO iq = 1, nqtot
c
         IF ( iadv(iq).EQ.1.OR.iadv(iq).EQ.2 )  THEN
            CALL traceur( iq,iadv,q,teta,pk,w, pbaru, pbarv, dq )

         ELSE IF( iq.EQ. nqtot )   THEN
c
            iapp_tracvl = 5
c
cccc     iapp_tracvl est la frequence en pas du groupement des flux
cccc      de masse pour  Van-Leer dans la routine  tracvl  .
c

            CALL vanleer(numvanle,iapp_tracvl,nqtot,q,pbaru,pbarv,
     *                      p, masse, dq,  iadv(1), teta, pk     )

c
c                   ...  Modif  F.Codron  ....
c
         ENDIF
c
        ENDDO
C
c        IF (offline) THEN
C maf stokage du flux de masse pour traceurs OFF-LINE

c           CALL fluxstokenc(pbaru,pbarv,masse,teta,phi,phis,
c    .   time_step, itau)

c        ENDIF
c
      ENDIF
	  END IF   ! tracer


c-----------------------------------------------------------------------
c   integrations dynamique et traceurs:
c   ----------------------------------

          if(calldyn)then
             CALL integrd(2,vcovm1,ucovm1,tetam1,psm1,massem1, 
     $            dv,du,dteta,dq,dp,vcov,ucov,teta,q,ps,masse,phis,
     $            finvmaold)
          else
             print*,'Currently no dynamics in this GCM...'
          endif





c .P.Le Van (26/04/94  ajout de  finvpold dans l'appel d'integrd)
c
c-----------------------------------------------------------------------
c   calcul des tendances physiques:
c   -------------------------------
c    ########   P.Le Van ( Modif le  6/02/95 )   ###########
c
       IF( purmats )  THEN
          IF( itau.EQ.itaufin.AND..NOT.forward ) lafin = .TRUE.
       ELSE
          IF( itau+1. EQ. itaufin )              lafin = .TRUE.
       ENDIF
c
c
       IF( apphys )  THEN
c
c     .......   Ajout   P.Le Van ( 17/04/96 )   ...........
c

         CALL pression (  ip1jmp1, ap, bp, ps,  p      )
         CALL exner_hyb(  ip1jmp1, ps, p,beta,pks, pk, pkf )

           rdaym_ini  = itau * dtvr / daysec
           rdayvrai   = rdaym_ini  + day_ini

! Ehouarn: what was this for ??
!           IF ( ecritphy.LT.1. )  THEN
!             rday_ecri = rdaym_ini
!           ELSE
             rday_ecri = INT(rdaym_ini)+INT(day_ini)
!           ENDIF
c



        CALL calfis( nqtot, lafin ,rdayvrai,rday_ecri,time  ,
     $                 ucov,vcov,teta,q,masse,ps,p,pk,phis,phi ,
     $     du,dv,dteta,dq,w, dufi,dvfi,dhfi,dqfi,dpfi)


c      ajout des tendances physiques:
c      ------------------------------
!        if(1.eq.2)then
          CALL addfi( nqtot, dtphys, leapf, forward   ,
     $                  ucov, vcov, teta , q   ,ps , masse,
     $                 dufi, dvfi, dhfi , dqfi ,dpfi  )
!       else
!          print*,'Currently no physics in this GCM...'
!       endif
c
       ENDIF

       CALL pression ( ip1jmp1, ap, bp, ps, p                  )

       CALL exner_hyb( ip1jmp1, ps, p,beta, pks, pk, pkf )
c   ----------------------------------------------------------

c
c
c   dissipation horizontale et verticale  des petites echelles:
c   ----------------------------------------------------------

      IF(apdiss) THEN

c        Sponge layer
c        ~~~~~~~~~~~~

         IF (callsponge) THEN
           pext(1:ip1jmp1)=ps(1:ip1jmp1)*aire(1:ip1jmp1)
           CALL sponge(ucov,vcov,teta,pext,dtdiss,mode_sponge)
         ENDIF

c        Dissipation horizontale
c        ~~~~~~~~~~~~~~~~~~~~~~~


         if(dissip_conservative) then
!     calculate kinetic energy before dissipation
            call covcont(llm,ucov,vcov,ucont,vcont)
            call enercin(vcov,ucov,vcont,ucont,ecin0)
         end if

         CALL dissip(vcov,ucov,teta,p,dvdis,dudis,dhdis)

         ucov(:,:)=ucov(:,:)+dudis(:,:)
         vcov(:,:)=vcov(:,:)+dvdis(:,:)


         if (dissip_conservative) then
C           On rajoute la tendance due a la transform. Ec -> E therm. cree
C           lors de la dissipation
            call covcont(llm,ucov,vcov,ucont,vcont)
            call enercin(vcov,ucov,vcont,ucont,ecin)
            dtetaecdt(:,:) = (ecin0(:,:)-ecin(:,:))/ pk(:,:)
            dhdis(:,:)     = dhdis(:,:) + dtetaecdt(:,:)
         endif

         teta(:,:)=teta(:,:)+dhdis(:,:)


c    .......        P. Le Van (  ajout  le 17/04/96  )   ...........
c   ...      Calcul de la valeur moyenne, unique de h aux poles  .....
c

        DO l  =  1, llm
          DO ij =  1,iim
           tppn(ij)  = aire(  ij    ) * teta(  ij    ,l)
           tpps(ij)  = aire(ij+ip1jm) * teta(ij+ip1jm,l)
          ENDDO

           tpn  = SSUM(iim,tppn,1)/apoln
           tps  = SSUM(iim,tpps,1)/apols

          DO ij = 1, iip1
           teta(  ij    ,l) = tpn
           teta(ij+ip1jm,l) = tps
          ENDDO
        ENDDO

        DO ij =  1,iim
          tppn(ij)  = aire(  ij    ) * ps (  ij    )
          tpps(ij)  = aire(ij+ip1jm) * ps (ij+ip1jm)
        ENDDO
          tpn  = SSUM(iim,tppn,1)/apoln

          tps  = SSUM(iim,tpps,1)/apols

        DO ij = 1, iip1
          ps(  ij    ) = tpn
          ps(ij+ip1jm) = tps
        ENDDO







      END IF
        
c   ********************************************************************
c   ********************************************************************
c   .... fin de l'integration dynamique  et physique pour le pas itau ..
c   ********************************************************************
c   ********************************************************************

c   preparation du pas d'integration suivant  ......

      IF ( .NOT.purmats ) THEN
c       ........................................................
c       ..............  schema matsuno + leapfrog  ..............
c       ........................................................

            IF(forward. OR. leapf) THEN
              itau= itau + 1
              iday= day_ini+itau/day_step
              time= FLOAT(itau-(iday-day_ini)*day_step)/day_step+time_0
                IF(time.GT.1.) THEN
                  time = time-1.
                  iday = iday+1
                ENDIF
            ENDIF


            IF( itau. EQ. itaufinp1 ) then  
              abort_message = 'Simulation finished'
              call abort_gcm(modname,abort_message,0)
            ENDIF
c-----------------------------------------------------------------------
c   ecriture du fichier histoire moyenne:
c   -------------------------------------

c           IF(MOD(itau,iperiod).EQ.0 .OR. itau.EQ.itaufin) THEN
c              IF(itau.EQ.itaufin) THEN
c                 iav=1
c              ELSE
c                 iav=0
c              ENDIF
c              CALL writedynav(histaveid, nqtot, itau,vcov ,
c    ,                          ucov,teta,pk,phi,q,masse,ps,phis)
c           ENDIF

c-----------------------------------------------------------------------


            IF(itau.EQ.itaufin) THEN


       PRINT *,' Appel test_period avant redem ', itau
       CALL test_period ( ucov,vcov,teta,q,p,phis )
       CALL dynredem1("restart.nc",0.0,
     .                     vcov,ucov,teta,q,masse,ps)

              CLOSE(99)
            ENDIF

c-----------------------------------------------------------------------
c   gestion de l'integration temporelle:
c   ------------------------------------

            IF( MOD(itau,iperiod).EQ.0 )    THEN
                    GO TO 1
            ELSE IF ( MOD(itau-1,iperiod). EQ. 0 ) THEN

                   IF( forward )  THEN
c      fin du pas forward et debut du pas backward

                      forward = .FALSE.
                        leapf = .FALSE.
                           GO TO 2

                   ELSE
c      fin du pas backward et debut du premier pas leapfrog

                        leapf =  .TRUE.
                        dt  =  2.*dtvr
                        GO TO 2
                   END IF
            ELSE

c      ......   pas leapfrog  .....

                 leapf = .TRUE.
                 dt  = 2.*dtvr
                 GO TO 2
            END IF

      ELSE

c       ........................................................
c       ..............       schema  matsuno        ...............
c       ........................................................
            IF( forward )  THEN

             itau =  itau + 1
             iday = day_ini+itau/day_step
             time = FLOAT(itau-(iday-day_ini)*day_step)/day_step+time_0

                  IF(time.GT.1.) THEN
                   time = time-1.
                   iday = iday+1
                  ENDIF

               forward =  .FALSE.
               IF( itau. EQ. itaufinp1 ) then  
                 abort_message = 'Simulation finished'
                 call abort_gcm(modname,abort_message,0)
               ENDIF
               GO TO 2

            ELSE

            IF(MOD(itau,iperiod).EQ.0 .OR. itau.EQ.itaufin) THEN
               IF(itau.EQ.itaufin) THEN
                  iav=1
               ELSE
                  iav=0
               ENDIF
c              CALL writedynav(histaveid, nqtot, itau,vcov ,
c    ,                          ucov,teta,pk,phi,q,masse,ps,phis)

            ENDIF


                 IF(itau.EQ.itaufin)
     . CALL dynredem1("restart.nc",0.0,
     .                     vcov,ucov,teta,q,nqtot,masse,ps)

                 forward = .TRUE.
                 GO TO  1


            ENDIF

      END IF

      STOP
      END
