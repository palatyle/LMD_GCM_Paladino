!
! $Id: leapfrog.F 1446 2010-10-22 09:27:25Z emillour $
!
c
c
      SUBROUTINE leapfrog(ucov,vcov,teta,ps,masse,phis,q,
     &                    time_0)


cIM : pour sortir les param. du modele dans un fis. netcdf 110106
#ifdef CPP_IOIPSL
      use IOIPSL
#endif
      USE infotrac, ONLY: nqtot,ok_iso_verif
      USE guide_mod, ONLY : guide_main
      USE write_field, ONLY: writefield
      USE control_mod, ONLY: planet_type,nday,day_step,iperiod,iphysiq,
     &                       less1day,fractday,ndynstep,iconser,
     &                       dissip_period,offline,ip_ebil_dyn,
     &                       ok_dynzon,periodav,ok_dyn_ave,iecri,
     &                       ok_dyn_ins,output_grads_dyn,ecritstart
      use exner_hyb_m, only: exner_hyb
      use exner_milieu_m, only: exner_milieu
      use cpdet_mod, only: cpdet,tpot2t,t2tpot
      use sponge_mod, only: callsponge,mode_sponge,sponge
       use comuforc_h
      USE comvert_mod, ONLY: ap,bp,pressure_exner,presnivs
      USE comconst_mod, ONLY: daysec,dtvr,dtphys,dtdiss,
     .			cpp,ihf,iflag_top_bound,pi
      USE logic_mod, ONLY: iflag_phys,ok_guide,forward,leapf,apphys,
     .			statcl,conser,apdiss,purmats,tidal,ok_strato
      USE temps_mod, ONLY: jD_ref,jH_ref,itaufin,day_ini,day_ref,
     .			start_time,dt,hour_ini

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

c  ... Possibilite de choisir le shema pour l'advection de
c        q  , en modifiant iadv dans traceur.def  (10/02) .
c
c      Pour Van-Leer + Vapeur d'eau saturee, iadv(1)=4. (F.Codron,10/99)
c      Pour Van-Leer iadv=10 
c
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

#include "dimensions.h"
#include "paramet.h"
#include "comdissnew.h"
#include "comgeom.h"
!#include "com_io_dyn.h"
#include "iniprint.h"
#include "academic.h"

! FH 2008/05/09 On elimine toutes les clefs physiques dans la dynamique
! #include "clesphys.h"

      REAL,INTENT(IN) :: time_0 ! not used

c   dynamical variables:
      REAL,INTENT(INOUT) :: ucov(ip1jmp1,llm)    ! zonal covariant wind
      REAL,INTENT(INOUT) :: vcov(ip1jm,llm)      ! meridional covariant wind
      REAL,INTENT(INOUT) :: teta(ip1jmp1,llm)    ! potential temperature
      REAL,INTENT(INOUT) :: ps(ip1jmp1)          ! surface pressure (Pa)
      REAL,INTENT(INOUT) :: masse(ip1jmp1,llm)   ! air mass
      REAL,INTENT(INOUT) :: phis(ip1jmp1)        ! geopotentiat at the surface
      REAL,INTENT(INOUT) :: q(ip1jmp1,llm,nqtot) ! advected tracers

      REAL p (ip1jmp1,llmp1  )               ! interlayer pressure
      REAL pks(ip1jmp1)                      ! exner at the surface
      REAL pk(ip1jmp1,llm)                   ! exner at mid-layer
      REAL pkf(ip1jmp1,llm)                  ! filtered exner at mid-layer
      REAL phi(ip1jmp1,llm)                  ! geopotential
      REAL w(ip1jmp1,llm)                    ! vertical velocity
! ADAPTATION GCM POUR CP(T)
      REAL temp(ip1jmp1,llm)                 ! temperature  
      REAL tsurpk(ip1jmp1,llm)               ! cpp*T/pk  

!      real zqmin,zqmax

c variables dynamiques intermediaire pour le transport
      REAL pbaru(ip1jmp1,llm),pbarv(ip1jm,llm) !flux de masse

c   variables dynamiques au pas -1
      REAL vcovm1(ip1jm,llm),ucovm1(ip1jmp1,llm)
      REAL tetam1(ip1jmp1,llm),psm1(ip1jmp1)
      REAL massem1(ip1jmp1,llm)

c   tendances dynamiques en */s
      REAL dv(ip1jm,llm),du(ip1jmp1,llm)
      REAL dteta(ip1jmp1,llm),dq(ip1jmp1,llm,nqtot),dp(ip1jmp1)

c   tendances de la dissipation en */s
      REAL dvdis(ip1jm,llm),dudis(ip1jmp1,llm)
      REAL dtetadis(ip1jmp1,llm)

c   tendances de la couche superieure */s
c      REAL dvtop(ip1jm,llm)
      REAL dutop(ip1jmp1,llm)
c      REAL dtetatop(ip1jmp1,llm)
c      REAL dqtop(ip1jmp1,llm,nqtot),dptop(ip1jmp1)

c   TITAN : tendances due au forces de marees */s
      REAL dvtidal(ip1jm,llm),dutidal(ip1jmp1,llm)

c   tendances physiques */s
      REAL dvfi(ip1jm,llm),dufi(ip1jmp1,llm)
      REAL dtetafi(ip1jmp1,llm),dqfi(ip1jmp1,llm,nqtot),dpfi(ip1jmp1)

c   variables pour le fichier histoire
!      REAL dtav      ! intervalle de temps elementaire
      LOGICAL lrestart

      REAL tppn(iim),tpps(iim),tpn,tps
c
      INTEGER itau,itaufinp1,iav
!      INTEGER  iday ! jour julien
      REAL       time 

      REAL  SSUM
!     REAL finvmaold(ip1jmp1,llm)

cym      LOGICAL  lafin
      LOGICAL :: lafin=.false.
      INTEGER ij,iq,l
!      INTEGER ik

!      real time_step, t_wrt, t_ops

      REAL rdaym_ini
! jD_cur: jour julien courant
! jH_cur: heure julienne courante
      REAL :: jD_cur, jH_cur
!      INTEGER :: an, mois, jour
!      REAL :: secondes

      LOGICAL first,callinigrads
cIM : pour sortir les param. du modele dans un fis. netcdf 110106
      save first
      data first/.true./
!      real dt_cum
!      character*10 infile
!      integer zan, tau0, thoriid
!      integer nid_ctesGCM
!      save nid_ctesGCM
!      real degres
!      real rlong(iip1), rlatg(jjp1)
!      real zx_tmp_2d(iip1,jjp1)
!      integer ndex2d(iip1*jjp1)
      logical ok_sync
      parameter (ok_sync = .true.) 
      logical physics

      data callinigrads/.true./
      character*10 string10

!      REAL alpha(ip1jmp1,llm),beta(ip1jmp1,llm)
      REAL :: flxw(ip1jmp1,llm)  ! flux de masse verticale

c+jld variables test conservation energie
      REAL ecin(ip1jmp1,llm),ecin0(ip1jmp1,llm)
C     Tendance de la temp. potentiel d (theta)/ d t due a la 
C     tansformation d'energie cinetique en energie thermique
C     cree par la dissipation
      REAL dtetaecdt(ip1jmp1,llm)
      REAL vcont(ip1jm,llm),ucont(ip1jmp1,llm)
      REAL vnat(ip1jm,llm),unat(ip1jmp1,llm)
!      REAL      d_h_vcol, d_qt, d_qw, d_ql, d_ec
      CHARACTER*15 ztit
!IM   INTEGER   ip_ebil_dyn  ! PRINT level for energy conserv. diag.
!IM   SAVE      ip_ebil_dyn
!IM   DATA      ip_ebil_dyn/0/
c-jld 

!      integer :: itau_w ! for write_paramLMDZ_dyn.h

!      character*80 dynhist_file, dynhistave_file
      character(len=*),parameter :: modname="leapfrog"
      character*80 abort_message

      logical dissip_conservative
      save dissip_conservative
      data dissip_conservative/.true./

      INTEGER testita
      PARAMETER (testita = 9)

      logical , parameter :: flag_verif = .false.
      
      ! for CP(T)
      real :: dtec
      real :: ztetaec(ip1jmp1,llm)

      if (nday>=0) then
         itaufin   = nday*day_step
      else
         ! to run a given (-nday) number of dynamical steps
         itaufin   = -nday
      endif
      if (less1day) then
c MODIF VENUS: to run less than one day:
        itaufin   = int(fractday*day_step)
      endif
      if (ndynstep.gt.0) then
        ! running a given number of dynamical steps
        itaufin=ndynstep
      endif
      itaufinp1 = itaufin +1
      
c INITIALISATIONS
        dudis(:,:)   =0.
        dvdis(:,:)   =0.
        dtetadis(:,:)=0.
        dutop(:,:)   =0.
c        dvtop(:,:)   =0.
c        dtetatop(:,:)=0.
c        dqtop(:,:,:) =0.
c        dptop(:)     =0.
        dufi(:,:)   =0.
        dvfi(:,:)   =0.
        dtetafi(:,:)=0.
        dqfi(:,:,:) =0.
        dpfi(:)     =0.

      itau = 0
      physics=.true.
      if (iflag_phys==0.or.iflag_phys==2) physics=.false.

c      iday = day_ini+itau/day_step
c      time = REAL(itau-(iday-day_ini)*day_step)/day_step+time_0
c         IF(time.GT.1.) THEN
c          time = time-1.
c          iday = iday+1
c         ENDIF


c-----------------------------------------------------------------------
c   On initialise la pression et la fonction d'Exner :
c   --------------------------------------------------

      dq(:,:,:)=0.
      CALL pression ( ip1jmp1, ap, bp, ps, p       )
      if (pressure_exner) then
        CALL exner_hyb( ip1jmp1, ps, p, pks, pk, pkf )
      else
        CALL exner_milieu( ip1jmp1, ps, p, pks, pk, pkf )
      endif

c------------------
c TEST PK MONOTONE
c------------------
      write(*,*) "Test PK"
      do ij=1,ip1jmp1
        do l=2,llm
          if(pk(ij,l).gt.pk(ij,l-1)) then
c           write(*,*) ij,l,pk(ij,l)
            abort_message = 'PK non strictement decroissante'
            call abort_gcm(modname,abort_message,1)
c           write(*,*) "ATTENTION, Test PK deconnecté..."
          endif
        enddo
      enddo
      write(*,*) "Fin Test PK"
c     stop
c------------------

c-----------------------------------------------------------------------
c   Debut de l'integration temporelle:
c   ----------------------------------

c     RMBY: check that hour_ini and start_time are not both non-zero
      if ((hour_ini.ne.0.0).and.(start_time.ne.0.0)) then
        write(*,*) "ERROR: hour_ini = ", hour_ini, 
     &             "start_time = ", start_time
        abort_message = 'hour_ini and start_time both nonzero'
        call abort_gcm(modname,abort_message,1)
      endif

   1  CONTINUE ! Matsuno Forward step begins here

c   date: (NB: date remains unchanged for Backward step)
c   -----

      jD_cur = jD_ref + day_ini - day_ref +                             &
     &          (itau+1)/day_step 
      IF (planet_type .eq. "mars") THEN
        jH_cur = jH_ref + hour_ini +                                    &
     &           mod(itau+1,day_step)/float(day_step) 
      ELSE
        jH_cur = jH_ref + start_time +                                  &
     &           mod(itau+1,day_step)/float(day_step)
      ENDIF
      jD_cur = jD_cur + int(jH_cur)
      jH_cur = jH_cur - int(jH_cur)

        if (ok_iso_verif) then
           call check_isotopes_seq(q,ip1jmp1,'leapfrog 321')
        endif !if (ok_iso_verif) then

#ifdef CPP_IOIPSL
      IF (planet_type.eq."earth") THEN
      if (ok_guide) then
        call guide_main(itau,ucov,vcov,teta,q,masse,ps)
      endif
      ENDIF
#endif


c
c     IF( MOD( itau, 10* day_step ).EQ.0 )  THEN
c       CALL  test_period ( ucov,vcov,teta,q,p,phis )
c       PRINT *,' ----   Test_period apres continue   OK ! -----', itau
c     ENDIF 
c

! Save fields obtained at previous time step as '...m1'
      CALL SCOPY( ijmllm ,vcov , 1, vcovm1 , 1 )
      CALL SCOPY( ijp1llm,ucov , 1, ucovm1 , 1 )
      CALL SCOPY( ijp1llm,teta , 1, tetam1 , 1 )
      CALL SCOPY( ijp1llm,masse, 1, massem1, 1 )
      CALL SCOPY( ip1jmp1, ps  , 1,   psm1 , 1 )

      forward = .TRUE.
      leapf   = .FALSE.
      dt      =  dtvr

c   ...    P.Le Van .26/04/94  ....
! Ehouarn: finvmaold is actually not used
!      CALL SCOPY   ( ijp1llm,   masse, 1, finvmaold,     1 )
!      CALL filtreg ( finvmaold ,jjp1, llm, -2,2, .TRUE., 1 )

        if (ok_iso_verif) then
           call check_isotopes_seq(q,ip1jmp1,'leapfrog 400')
        endif !if (ok_iso_verif) then

   2  CONTINUE ! Matsuno backward or leapfrog step begins here

c-----------------------------------------------------------------------

c   date: (NB: only leapfrog step requires recomputing date)
c   -----

      IF (leapf) THEN
        jD_cur = jD_ref + day_ini - day_ref +
     &            (itau+1)/day_step
        IF (planet_type .eq. "mars") THEN
          jH_cur = jH_ref + hour_ini + 
     &             mod(itau+1,day_step)/float(day_step) 
        ELSE
          jH_cur = jH_ref + start_time +
     &             mod(itau+1,day_step)/float(day_step)
        ENDIF
        jD_cur = jD_cur + int(jH_cur)
        jH_cur = jH_cur - int(jH_cur)
      ENDIF


c   gestion des appels de la physique et des dissipations:
c   ------------------------------------------------------
c
c   ...    P.Le Van  ( 6/02/95 )  ....

      apphys = .FALSE.
      statcl = .FALSE.
      conser = .FALSE.
      apdiss = .FALSE.

      IF( purmats ) THEN
      ! Purely Matsuno time stepping
         IF( MOD(itau,iconser) .EQ.0.AND.  forward    ) conser = .TRUE.
         IF( MOD(itau,dissip_period ).EQ.0.AND..NOT.forward ) 
     s        apdiss = .TRUE.
         IF( MOD(itau,iphysiq ).EQ.0.AND..NOT.forward 
     s          .and. physics                        ) apphys = .TRUE.
      ELSE
      ! Leapfrog/Matsuno time stepping 
         IF( MOD(itau   ,iconser) .EQ. 0              ) conser = .TRUE.
         IF( MOD(itau+1,dissip_period).EQ.0 .AND. .NOT. forward )
     s        apdiss = .TRUE.
         IF( MOD(itau+1,iphysiq).EQ.0.AND.physics       ) apphys=.TRUE.
      END IF

! Ehouarn: for Shallow Water case (ie: 1 vertical layer),
!          supress dissipation step
      if (llm.eq.1) then
        apdiss=.false.
      endif

#ifdef NODYN
      apdiss=.false.
#endif

        if (ok_iso_verif) then
           call check_isotopes_seq(q,ip1jmp1,'leapfrog 589')
        endif !if (ok_iso_verif) then

c-----------------------------------------------------------------------
c   calcul des tendances dynamiques:
c   --------------------------------

! ADAPTATION GCM POUR CP(T)
      call tpot2t(ijp1llm,teta,temp,pk)
      tsurpk = cpp*temp/pk
      ! compute geopotential phi()
      CALL geopot  ( ip1jmp1, tsurpk  , pk , pks,  phis  , phi   )

           rdaym_ini  = itau * dtvr / daysec

      time = jD_cur + jH_cur

#ifdef NODYN
      dv(:,:) = 0.D+0
      du(:,:) = 0.D+0
      dteta(:,:) = 0.D+0
      dq(:,:,:) = 0.D+0
      dp(:) = 0.D+0
#else
      CALL caldyn 
     $  ( itau,ucov,vcov,teta,ps,masse,pk,pkf,tsurpk,phis ,
     $    phi,conser,du,dv,dteta,dp,w, pbaru,pbarv, time )

      ! Simple zonal wind nudging for generic planetary model
      ! AS 09/2013
      ! ---------------------------------------------------
      if (planet_type.eq."generic") then
       if (ok_guide) then
         du(:,:) = du(:,:) + ((uforc(:,:)-ucov(:,:)) / facwind)
       endif
      endif

c-----------------------------------------------------------------------
c   calcul des tendances advection des traceurs (dont l'humidite)
c   -------------------------------------------------------------

        if (ok_iso_verif) then
           call check_isotopes_seq(q,ip1jmp1,
     &           'leapfrog 686: avant caladvtrac')
        endif !if (ok_iso_verif) then

      IF( forward. OR . leapf )  THEN
! Ehouarn: NB: fields sent to advtrac are those at the beginning of the time step
         CALL caladvtrac(q,pbaru,pbarv,
     *        p, masse, dq,  teta,
     .        flxw, pk)
         
         IF (offline) THEN
Cmaf stokage du flux de masse pour traceurs OFF-LINE

#ifdef CPP_IOIPSL
           CALL fluxstokenc(pbaru,pbarv,masse,teta,phi,phis,
     .   dtvr, itau)
#endif


         ENDIF ! of IF (offline)
c
      ENDIF ! of IF( forward. OR . leapf )


c-----------------------------------------------------------------------
c   integrations dynamique et traceurs:
c   ----------------------------------

        if (ok_iso_verif) then
           write(*,*) 'leapfrog 720' 
           call check_isotopes_seq(q,ip1jmp1,'leapfrog 756')
        endif !if (ok_iso_verif) then
        
       CALL integrd ( nqtot,vcovm1,ucovm1,tetam1,psm1,massem1 ,
     $         dv,du,dteta,dq,dp,vcov,ucov,teta,q,ps,masse,phis )
!     $              finvmaold                                    )

       if (ok_iso_verif) then
          write(*,*) 'leapfrog 724'
           call check_isotopes_seq(q,ip1jmp1,'leapfrog 762')
        endif !if (ok_iso_verif) then

       IF ((planet_type.eq."titan").and.(tidal)) then
c-----------------------------------------------------------------------
c   Marées gravitationnelles causées par Saturne
c   B. Charnay (28/10/2010)
c   ----------------------------------------------------------
            CALL tidal_forces(rdaym_ini, dutidal, dvtidal)
            ucov=ucov+dutidal*dt
            vcov=vcov+dvtidal*dt
       ENDIF

! NODYN precompiling flag
#endif

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
         if (pressure_exner) then
           CALL exner_hyb(  ip1jmp1, ps, p,pks, pk, pkf )
         else
           CALL exner_milieu( ip1jmp1, ps, p, pks, pk, pkf )
         endif

! Compute geopotential (physics might need it)
! GEOP CORRECTION
! ADAPTATION GCM POUR CP(T)
!        CALL geopot  ( ip1jmp1, teta  , pk , pks,  phis  , phi   )
         call tpot2t(ijp1llm,teta,temp,pk)
         tsurpk = cpp*temp/pk
         CALL geopot  ( ip1jmp1, tsurpk, pk, pks, phis, phi )

           jD_cur = jD_ref + day_ini - day_ref +                        &
     &          (itau+1)/day_step 

           IF ((planet_type .eq."generic").or.
     &         (planet_type .eq."mars")) THEN
              ! AS: we make jD_cur to be pday
              jD_cur = int(day_ini + itau/day_step)
           ENDIF

           IF (planet_type .eq. "mars") THEN
             jH_cur = jH_ref + hour_ini +                                 &
     &                mod(itau,day_step)/float(day_step)
           ELSE IF (planet_type .eq. "generic") THEN
             jH_cur = jH_ref + start_time +                               &
     &                mod(itau,day_step)/float(day_step)
           ELSE
             jH_cur = jH_ref + start_time +                               &
     &                mod(itau+1,day_step)/float(day_step)
           ENDIF
           jD_cur = jD_cur + int(jH_cur)
           jH_cur = jH_cur - int(jH_cur)
!         write(lunout,*)'itau, jD_cur = ', itau, jD_cur, jH_cur
!         call ju2ymds(jD_cur+jH_cur, an, mois, jour, secondes)
!         write(lunout,*)'current date = ',an, mois, jour, secondes 

c rajout debug
c       lafin = .true.


c   Interface avec les routines de phylmd (phymars ... )
c   -----------------------------------------------------

c+jld

c  Diagnostique de conservation de l'énergie : initialisation
         IF (ip_ebil_dyn.ge.1 ) THEN 
          ztit='bil dyn'
! Ehouarn: be careful, diagedyn is Earth-specific!
           IF (planet_type.eq."earth") THEN
            CALL diagedyn(ztit,2,1,1,dtphys
     &    , ucov    , vcov , ps, p ,pk , teta , q(:,:,1), q(:,:,2))
           ENDIF
         ENDIF ! of IF (ip_ebil_dyn.ge.1 )
c-jld
#ifdef CPP_IOIPSL
cIM decommenter les 6 lignes suivantes pour sortir quelques parametres dynamiques de LMDZ 
cIM uncomment next 6 lines to get some parameters for LMDZ dynamics
c        IF (first) THEN
c         first=.false.
c#include "ini_paramLMDZ_dyn.h"
c        ENDIF
c
c#include "write_paramLMDZ_dyn.h"
c
#endif
! #endif of #ifdef CPP_IOIPSL

c          call WriteField('pfi',reshape(p,(/iip1,jmp1,llmp1/)))

         CALL calfis( lafin , jD_cur, jH_cur,
     $               ucov,vcov,teta,q,masse,ps,p,pk,phis,phi ,
     $               du,dv,dteta,dq,
     $               flxw,
     $               dufi,dvfi,dtetafi,dqfi,dpfi  )

c          call WriteField('dufi',reshape(dufi,(/iip1,jmp1,llm/)))
c          call WriteField('dvfi',reshape(dvfi,(/iip1,jjm,llm/)))
c          call WriteField('dtetafi',reshape(dtetafi,(/iip1,jmp1,llm/)))

c      ajout des tendances physiques:
c      ------------------------------
          CALL addfi( dtphys, leapf, forward   ,
     $                  ucov, vcov, teta , q   ,ps ,
     $                 dufi, dvfi, dtetafi , dqfi ,dpfi  )
          ! since addfi updates ps(), also update p(), masse() and pk()
          CALL pression (ip1jmp1,ap,bp,ps,p)
          CALL massdair(p,masse)
          if (pressure_exner) then
            CALL exner_hyb( ip1jmp1, ps, p, pks, pk, pkf )
          else
            CALL exner_milieu( ip1jmp1, ps, p, pks, pk, pkf )
          endif
          
c      Couche superieure :
c      -------------------
         IF (iflag_top_bound > 0) THEN
           CALL top_bound(vcov,ucov,teta,masse,dtphys,dutop)
           dutop(:,:)=dutop(:,:)/dtphys   ! convert to a tendency in (m/s)/s
         ENDIF

c  Diagnostique de conservation de l'énergie : difference
         IF (ip_ebil_dyn.ge.1 ) THEN 
          ztit='bil phys'
          IF (planet_type.eq."earth") THEN
           CALL diagedyn(ztit,2,1,1,dtphys
     &     , ucov    , vcov , ps, p ,pk , teta , q(:,:,1), q(:,:,2))
          ENDIF
         ENDIF ! of IF (ip_ebil_dyn.ge.1 )

       ENDIF ! of IF( apphys )

      IF(iflag_phys.EQ.2) THEN ! "Newtonian" case
!   Academic case : Simple friction and Newtonan relaxation 
!   -------------------------------------------------------
        DO l=1,llm   
          DO ij=1,ip1jmp1
           teta(ij,l)=teta(ij,l)-dtvr*
     &      (teta(ij,l)-tetarappel(ij,l))*(knewt_g+knewt_t(l)*clat4(ij))
          ENDDO
        ENDDO ! of DO l=1,llm 
    
        if (planet_type.eq."giant") then
          ! add an intrinsic heat flux at the base of the atmosphere
          teta(:,1)=teta(:,1)+dtvr*aire(:)*ihf/cpp/masse(:,1)
        endif

        call friction(ucov,vcov,dtvr)

    
        ! Sponge layer (if any)
        IF (ok_strato) THEN
           CALL top_bound(vcov,ucov,teta,masse,dtvr,dutop)
           dutop(:,:)=dutop(:,:)/dtvr   ! convert to a tendency in (m/s)/s
        ENDIF ! of IF (ok_strato) 
      ENDIF ! of IF (iflag_phys.EQ.2)


c-jld

        CALL pression ( ip1jmp1, ap, bp, ps, p                  )
        if (pressure_exner) then
          CALL exner_hyb( ip1jmp1, ps, p, pks, pk, pkf )
        else
          CALL exner_milieu( ip1jmp1, ps, p, pks, pk, pkf )
        endif
        CALL massdair(p,masse)

        if (ok_iso_verif) then
           call check_isotopes_seq(q,ip1jmp1,'leapfrog 1196')
        endif !if (ok_iso_verif) then

c-----------------------------------------------------------------------
c   dissipation horizontale et verticale  des petites echelles:
c   ----------------------------------------------------------

      IF(apdiss) THEN

        ! sponge layer
        if (callsponge) then
          CALL sponge(ucov,vcov,teta,ps,dtdiss,mode_sponge)
        endif

c   calcul de l'energie cinetique avant dissipation
        call covcont(llm,ucov,vcov,ucont,vcont)
        call enercin(vcov,ucov,vcont,ucont,ecin0)

c   dissipation
! ADAPTATION GCM POUR CP(T)
        call tpot2t(ijp1llm,teta,temp,pk)

        CALL dissip(vcov,ucov,teta,p,dvdis,dudis,dtetadis)
        ucov=ucov+dudis
        vcov=vcov+dvdis
        dudis=dudis/dtdiss   ! passage en (m/s)/s
        dvdis=dvdis/dtdiss   ! passage en (m/s)/s

c------------------------------------------------------------------------
        if (dissip_conservative) then
C       On rajoute la tendance due a la transform. Ec -> E therm. cree
C       lors de la dissipation
            call covcont(llm,ucov,vcov,ucont,vcont)
            call enercin(vcov,ucov,vcont,ucont,ecin)
! ADAPTATION GCM POUR CP(T)
            do ij=1,ip1jmp1
              do l=1,llm
                dtec = (ecin0(ij,l)-ecin(ij,l))/cpdet(temp(ij,l))
                temp(ij,l) = temp(ij,l) + dtec
              enddo
            enddo
            call t2tpot(ijp1llm,temp,ztetaec,pk)
            dtetaecdt=ztetaec-teta
            dtetadis=dtetadis+dtetaecdt
        endif
        teta=teta+dtetadis
        dtetadis=dtetadis/dtdiss   ! passage en K/s
c------------------------------------------------------------------------


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

        if (1 == 0) then
!!! Ehouarn: lines here 1) kill 1+1=2 in the dynamics
!!!                     2) should probably not be here anyway
!!! but are kept for those who would want to revert to previous behaviour
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
        endif ! of if (1 == 0)

      END IF ! of IF(apdiss)

c ajout debug
c              IF( lafin ) then  
c                abort_message = 'Simulation finished'
c                call abort_gcm(modname,abort_message,0)
c              ENDIF
        
c   ********************************************************************
c   ********************************************************************
c   .... fin de l'integration dynamique  et physique pour le pas itau ..
c   ********************************************************************
c   ********************************************************************

c   preparation du pas d'integration suivant  ......

        if (ok_iso_verif) then
           call check_isotopes_seq(q,ip1jmp1,'leapfrog 1509')
        endif !if (ok_iso_verif) then

      IF ( .NOT.purmats ) THEN
c       ........................................................
c       ..............  schema matsuno + leapfrog  ..............
c       ........................................................

            IF(forward. OR. leapf) THEN
              itau= itau + 1
c              iday= day_ini+itau/day_step
c              time= REAL(itau-(iday-day_ini)*day_step)/day_step+time_0
c                IF(time.GT.1.) THEN
c                  time = time-1.
c                  iday = iday+1
c                ENDIF
            ENDIF


            IF( itau. EQ. itaufinp1 ) then  
              if (flag_verif) then
                write(79,*) 'ucov',ucov
                write(80,*) 'vcov',vcov
                write(81,*) 'teta',teta
                write(82,*) 'ps',ps
                write(83,*) 'q',q
                WRITE(85,*) 'q1 = ',q(:,:,1)
                WRITE(86,*) 'q3 = ',q(:,:,3)
              endif

              abort_message = 'Simulation finished'

              call abort_gcm(modname,abort_message,0)
            ENDIF
c-----------------------------------------------------------------------
c   ecriture du fichier histoire moyenne:
c   -------------------------------------

            IF(MOD(itau,iperiod).EQ.0 .OR. itau.EQ.itaufin) THEN
               IF(itau.EQ.itaufin) THEN
                  iav=1
               ELSE
                  iav=0
               ENDIF
               
!              ! Ehouarn: re-compute geopotential for outputs
! GEOP CORRECTION
! ADAPTATION GCM POUR CP(T)
!              CALL geopot(ip1jmp1,teta,pk,pks,phis,phi)
               call tpot2t(ijp1llm,teta,temp,pk)
               tsurpk = cpp*temp/pk
               CALL geopot(ip1jmp1,tsurpk,pk,pks,phis,phi)

               IF (ok_dynzon) THEN
#ifdef CPP_IOIPSL
c les traceurs ne sont pas sortis, trop lourd. 
c Peut changer eventuellement si besoin.
                 CALL bilan_dyn(dtvr*iperiod,dtvr*day_step*periodav,
     &                 ps,masse,pk,pbaru,pbarv,teta,phi,ucov,vcov,
     &                 du,dudis,dutop,dufi)
#endif
               END IF
               IF (ok_dyn_ave) THEN
#ifdef CPP_IOIPSL
                 CALL writedynav(itau,vcov,
     &                 ucov,teta,pk,phi,q,masse,ps,phis)
#endif
               ENDIF

            ENDIF ! of IF((MOD(itau,iperiod).EQ.0).OR.(itau.EQ.itaufin))

        if (ok_iso_verif) then
           call check_isotopes_seq(q,ip1jmp1,'leapfrog 1584')
        endif !if (ok_iso_verif) then

c-----------------------------------------------------------------------
c   ecriture de la bande histoire:
c   ------------------------------

            IF( MOD(itau,iecri).EQ.0) THEN
             ! Ehouarn: output only during LF or Backward Matsuno
	     if (leapf.or.(.not.leapf.and.(.not.forward))) then
! ADAPTATION GCM POUR CP(T)
              call tpot2t(ijp1llm,teta,temp,pk)
              tsurpk = cpp*temp/pk
              CALL geopot(ip1jmp1,tsurpk,pk,pks,phis,phi)
              unat=0.
              do l=1,llm
                unat(iip2:ip1jm,l)=ucov(iip2:ip1jm,l)/cu(iip2:ip1jm)
                vnat(:,l)=vcov(:,l)/cv(:)
              enddo
#ifdef CPP_IOIPSL
              if (ok_dyn_ins) then
!               write(lunout,*) "leapfrog: call writehist, itau=",itau
	       CALL writehist(itau,vcov,ucov,teta,phi,q,masse,ps,phis)
!               call WriteField('ucov',reshape(ucov,(/iip1,jmp1,llm/)))
!               call WriteField('vcov',reshape(vcov,(/iip1,jjm,llm/)))
!              call WriteField('teta',reshape(teta,(/iip1,jmp1,llm/)))
!               call WriteField('ps',reshape(ps,(/iip1,jmp1/)))
!               call WriteField('masse',reshape(masse,(/iip1,jmp1,llm/)))
              endif ! of if (ok_dyn_ins)
#endif
! For some Grads outputs of fields
              if (output_grads_dyn) then
#include "write_grads_dyn.h"
              endif
             endif ! of if (leapf.or.(.not.leapf.and.(.not.forward)))
            ENDIF ! of IF(MOD(itau,iecri).EQ.0)

c           Determine whether to write to the restart.nc file
c           Decision can't be made in one IF statement as if
c           ecritstart==0 there will be a divide-by-zero error
            lrestart = .false.
            IF (itau.EQ.itaufin) THEN
              lrestart = .true.
            ELSE IF (ecritstart.GT.0) THEN
              IF (MOD(itau,ecritstart).EQ.0) lrestart  = .true.
            ENDIF

c           Write to the restart.nc file
            IF (lrestart) THEN
              if (planet_type=="mars") then
                CALL dynredem1("restart.nc",REAL(itau)/REAL(day_step),
     &                         vcov,ucov,teta,q,masse,ps)
              else
                CALL dynredem1("restart.nc",start_time,
     &                         vcov,ucov,teta,q,masse,ps)
              endif
              CLOSE(99)
              !!! Ehouarn: Why not stop here and now?
            ENDIF ! of IF (lrestart)

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
                   END IF ! of IF (forward)
            ELSE

c      ......   pas leapfrog  .....

                 leapf = .TRUE.
                 dt  = 2.*dtvr
                 GO TO 2
            END IF ! of IF (MOD(itau,iperiod).EQ.0)
                   !    ELSEIF (MOD(itau-1,iperiod).EQ.0)

      ELSE ! of IF (.not.purmats)

        if (ok_iso_verif) then
           call check_isotopes_seq(q,ip1jmp1,'leapfrog 1664')
        endif !if (ok_iso_verif) then

c       ........................................................
c       ..............       schema  matsuno        ...............
c       ........................................................
            IF( forward )  THEN

             itau =  itau + 1
c             iday = day_ini+itau/day_step
c             time = REAL(itau-(iday-day_ini)*day_step)/day_step+time_0
c
c                  IF(time.GT.1.) THEN
c                   time = time-1.
c                   iday = iday+1
c                  ENDIF

               forward =  .FALSE.
               IF( itau. EQ. itaufinp1 ) then  
                 abort_message = 'Simulation finished'
                 call abort_gcm(modname,abort_message,0)
               ENDIF
               GO TO 2

            ELSE ! of IF(forward) i.e. backward step

        if (ok_iso_verif) then
           call check_isotopes_seq(q,ip1jmp1,'leapfrog 1698')
        endif !if (ok_iso_verif) then  

              IF(MOD(itau,iperiod).EQ.0 .OR. itau.EQ.itaufin) THEN
               IF(itau.EQ.itaufin) THEN
                  iav=1
               ELSE
                  iav=0
               ENDIF

!              ! Ehouarn: re-compute geopotential for outputs
! GEOP CORRECTION
! ADAPTATION GCM POUR CP(T)
!              CALL geopot(ip1jmp1,teta,pk,pks,phis,phi)
               call tpot2t(ijp1llm,teta,temp,pk)
               tsurpk = cpp*temp/pk
               CALL geopot(ip1jmp1,tsurpk,pk,pks,phis,phi)

               IF (ok_dynzon) THEN 
#ifdef CPP_IOIPSL
c les traceurs ne sont pas sortis, trop lourd. 
c Peut changer eventuellement si besoin.
                 CALL bilan_dyn(dtvr*iperiod,dtvr*day_step*periodav,
     &                 ps,masse,pk,pbaru,pbarv,teta,phi,ucov,vcov,
     &                 du,dudis,dutop,dufi)
#endif
               ENDIF
               IF (ok_dyn_ave) THEN
#ifdef CPP_IOIPSL
                 CALL writedynav(itau,vcov,
     &                 ucov,teta,pk,phi,q,masse,ps,phis)
#endif
               ENDIF

              ENDIF ! of IF(MOD(itau,iperiod).EQ.0 .OR. itau.EQ.itaufin)

              IF(MOD(itau,iecri         ).EQ.0) THEN
c              IF(MOD(itau,iecri*day_step).EQ.0) THEN
! ADAPTATION GCM POUR CP(T)
                call tpot2t(ijp1llm,teta,temp,pk)
                tsurpk = cpp*temp/pk
                CALL geopot(ip1jmp1,tsurpk,pk,pks,phis,phi)
                unat=0.
                do l=1,llm
                  unat(iip2:ip1jm,l)=ucov(iip2:ip1jm,l)/cu(iip2:ip1jm)
                  vnat(:,l)=vcov(:,l)/cv(:)
                enddo
#ifdef CPP_IOIPSL
              if (ok_dyn_ins) then
!                write(lunout,*) "leapfrog: call writehist (b)",
!     &                        itau,iecri
		CALL writehist(itau,vcov,ucov,teta,phi,q,masse,ps,phis)
              endif ! of if (ok_dyn_ins)
#endif
! For some Grads outputs
                if (output_grads_dyn) then
#include "write_grads_dyn.h"
                endif

              ENDIF ! of IF(MOD(itau,iecri         ).EQ.0) 

c             Determine whether to write to the restart.nc file
c             Decision can't be made in one IF statement as if
c             ecritstart==0 there will be a divide-by-zero error
              lrestart = .false.
              IF (itau.EQ.itaufin) THEN
                lrestart = .true.
              ELSE IF (ecritstart.GT.0) THEN
                IF (MOD(itau,ecritstart).EQ.0) lrestart = .true.
              ENDIF

c             Write to the restart.nc file
              IF (lrestart) THEN
                if (planet_type=="mars") then
                  CALL dynredem1("restart.nc",REAL(itau)/REAL(day_step),
     &                         vcov,ucov,teta,q,masse,ps)
                else
                  CALL dynredem1("restart.nc",start_time,
     &                         vcov,ucov,teta,q,masse,ps)
                endif
              ENDIF ! of IF (lrestart)

              forward = .TRUE.
              GO TO  1

            ENDIF ! of IF (forward)

      END IF ! of IF(.not.purmats)

      STOP
      END
