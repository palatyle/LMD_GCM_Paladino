!
! $Id: leapfrog.F 1446 2010-10-22 09:27:25Z emillour $
!
c
c
      SUBROUTINE leapfrog_nogcm(ucov,vcov,teta,ps,masse,phis,q,
     &                    time_0)


cIM : pour sortir les param. du modele dans un fis. netcdf 110106
#ifdef CPP_IOIPSL
      use IOIPSL
#endif
      USE infotrac, ONLY: nqtot,ok_iso_verif,tname
      USE guide_mod, ONLY : guide_main
      USE write_field, ONLY: writefield
      USE control_mod, ONLY: planet_type,nday,day_step,iperiod,iphysiq,
     &                       less1day,fractday,ndynstep,iconser,
     &                       dissip_period,offline,ip_ebil_dyn,
     &                       ok_dynzon,periodav,ok_dyn_ave,iecri,
     &                       ok_dyn_ins,output_grads_dyn
      use exner_hyb_m, only: exner_hyb
      use exner_milieu_m, only: exner_milieu
      use cpdet_mod, only: cpdet,tpot2t,t2tpot
      use sponge_mod, only: callsponge,mode_sponge,sponge
       use comuforc_h
      USE comvert_mod, ONLY: ap,bp,pressure_exner,presnivs,
     &                   aps,bps,presnivs,pseudoalt,preff,scaleheight
      USE comconst_mod, ONLY: daysec,dtvr,dtphys,dtdiss,
     .			cpp,ihf,iflag_top_bound,pi,kappa,r
      USE logic_mod, ONLY: iflag_phys,ok_guide,forward,leapf,apphys,
     .			statcl,conser,purmats,tidal,ok_strato
      USE temps_mod, ONLY: jD_ref,jH_ref,itaufin,day_ini,day_ref,
     .			start_time,dt




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
      REAL kpd(ip1jmp1)                       ! TB15 = exp (-z/H)

! ADAPTATION GCM POUR CP(T)
      REAL temp(ip1jmp1,llm)                 ! temperature  
      REAL tsurpk(ip1jmp1,llm)               ! cpp*T/pk  

!      real zqmin,zqmax


!       ED18 nogcm
      REAL tau_ps
      REAL tau_co2
      REAL tau_teta
      REAL tetadpmean
      REAL tetadp(ip1jmp1,llm)
      REAL dpmean
      REAL dp(ip1jmp1,llm)
      REAL tetamean
      real globaverage2d
!      LOGICAL firstcall_globaverage2d


c variables dynamiques intermediaire pour le transport
      REAL pbaru(ip1jmp1,llm),pbarv(ip1jm,llm) !flux de masse

c   variables dynamiques au pas -1
      REAL vcovm1(ip1jm,llm),ucovm1(ip1jmp1,llm)
      REAL tetam1(ip1jmp1,llm),psm1(ip1jmp1)
      REAL massem1(ip1jmp1,llm)

c   tendances dynamiques en */s
      REAL dv(ip1jm,llm),du(ip1jmp1,llm)
      REAL dteta(ip1jmp1,llm),dq(ip1jmp1,llm,nqtot)

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

      REAL psmean                    ! pression moyenne 
      REAL pqco2mean                 ! moyenne globale ps*qco
      REAL p0                        ! pression de reference
      REAL p00d                        ! globalaverage(kpd)
      REAL qmean_co2,qmean_co2_vert ! mass mean mixing ratio vap co2
      REAL pqco2(ip1jmp1)           ! average co2 mass index : ps*q_co2
      REAL oldps(ip1jmp1)           ! saving old pressure ps to calculate qch4

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

      ! TEMP : diagnostic mass
      real :: co2mass(iip1,jjp1)
      real :: co2ice_ij(iip1,jjp1)
      integer,save :: igcm_co2=0 ! index of CO2 tracer (if any)
      integer :: i,j,ig
      integer, parameter :: ngrid = 2+(jjm-1)*iim
      ! local mass
      real :: mq(ip1jmp1,llm)

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

! ED18 nogcm
!      firstcall_globaverage2d = .true.

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
!          PRINT*,'pk(ij,l) = ',pk(ij,l)
!          PRINT*,'pk(ij,l-1) = ',pk(ij,l-1)
          if(pk(ij,l).gt.pk(ij,l-1)) then
c           write(*,*) ij,l,pk(ij,l)
            abort_message = 'PK non strictement decroissante'
            call abort_gcm(modname,abort_message,1)
c           write(*,*) "ATTENTION, Test PK deconnecte..."
          endif
        enddo
      enddo
      write(*,*) "Fin Test PK"
c     stop
c------------------
c     Preparing mixing of pressure and tracers in nogcm 
c     kpd=exp(-z/H)  Needed to define a reference pressure :
c     P0=pmean/globalaverage(kpd)  
c     thus P(i) = P0*exp(-z(i)/H) = P0*kpd(i)
c     It is checked that globalaverage2d(Pi)=pmean
      DO ij=1,ip1jmp1
             kpd(ij) = exp(-phis(ij)/(r*200.))
      ENDDO
      p00d=globaverage2d(kpd) ! mean pres at ref level
      tau_ps = 1. ! constante de rappel for pressure  (s) 
      tau_co2 = 1.E5 !E5 ! constante de rappel for mix ratio qco2 (s) 
      tau_teta = 1.E7 !constante de rappel for potentiel temperature

! ED18 TEST
!      PRINT*,'igcm_co2 = ',igcm_co2
! Locate tracer "co2" and set igcm_co2:
      do iq=1,nqtot
        if (tname(iq)=="co2") then
          igcm_co2=iq
          exit
        endif
      enddo

c-----------------------------------------------------------------------
c   Debut de l'integration temporelle:
c   ----------------------------------

   1  CONTINUE ! Matsuno Forward step begins here

c   date: (NB: date remains unchanged for Backward step)
c   -----

      jD_cur = jD_ref + day_ini - day_ref +                             &
     &          (itau+1)/day_step 
      jH_cur = jH_ref + start_time +                                    &
     &          mod(itau+1,day_step)/float(day_step) 
      jD_cur = jD_cur + int(jH_cur)
      jH_cur = jH_cur - int(jH_cur)

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

   2  CONTINUE ! Matsuno backward or leapfrog step begins here

c-----------------------------------------------------------------------

c   date: (NB: only leapfrog step requires recomputing date)
c   -----

      IF (leapf) THEN
        jD_cur = jD_ref + day_ini - day_ref +
     &            (itau+1)/day_step
        jH_cur = jH_ref + start_time +
     &            mod(itau+1,day_step)/float(day_step) 
        jD_cur = jD_cur + int(jH_cur)
        jH_cur = jH_cur - int(jH_cur)
      ENDIF


c   gestion des appels de la physique et des dissipations:
c   ------------------------------------------------------
c
c   ...    P.Le Van  ( 6/02/95 )  ....

! ED18: suppression des mentions de la variable apdiss dans le cas
! 'nogcm'

      apphys = .FALSE.
      statcl = .FALSE.
      conser = .FALSE.
     

      IF( purmats ) THEN
      ! Purely Matsuno time stepping
         IF( MOD(itau,iconser) .EQ.0.AND.  forward    ) conser = .TRUE.
    
   
         IF( MOD(itau,iphysiq ).EQ.0.AND..NOT.forward 
     s          .and. physics        ) apphys = .TRUE.
      ELSE
      ! Leapfrog/Matsuno time stepping 
        IF( MOD(itau   ,iconser) .EQ. 0  ) conser = .TRUE.
  
 
        IF( MOD(itau+1,iphysiq).EQ.0.AND.physics) apphys=.TRUE.
      END IF

! Ehouarn: for Shallow Water case (ie: 1 vertical layer),
!          supress dissipation step



c-----------------------------------------------------------------------
c   calcul des tendances dynamiques:
c   --------------------------------

! ED18: suppression de l'onglet pour le cas nogcm


         dv(:,:) = 0.
         du(:,:) = 0.
         dteta(:,:) = 0.
         dq(:,:,:) = 0.
          


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
         CALL geopot  ( ip1jmp1, teta  , pk , pks,  phis  , phi   )

           jD_cur = jD_ref + day_ini - day_ref +                        &
     &          (itau+1)/day_step 

           IF ((planet_type .eq."generic").or.
     &         (planet_type .eq."mars")) THEN
              ! AS: we make jD_cur to be pday
              jD_cur = int(day_ini + itau/day_step)
           ENDIF

!           print*,'itau =',itau
!           print*,'day_step =',day_step
!           print*,'itau/day_step =',itau/day_step


           jH_cur = jH_ref + start_time +                               &
     &          mod(itau+1,day_step)/float(day_step) 
           IF ((planet_type .eq."generic").or.
     &         (planet_type .eq."mars")) THEN
             jH_cur = jH_ref + start_time +                               &
     &          mod(itau,day_step)/float(day_step)
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

c  Diagnostique de conservation de l'Energie : initialisation
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
!         print*,'---------LEAPFROG---------------'
!         print*,''
!         print*,'> AVANT CALFIS'
!         print*,''
!         print*,'teta(3049,:) = ',teta(3049,:)
!         print*,''
!         print*,'dteta(3049,:) = ',dteta(3049,:)
!         print*,''
!         print*,'dtetafi(3049,:) = ',dtetafi(3049,:)
!         print*,''


         CALL calfis( lafin , jD_cur, jH_cur,
     $               ucov,vcov,teta,q,masse,ps,p,pk,phis,phi ,
     $               du,dv,dteta,dq,
     $               flxw,
     $               dufi,dvfi,dtetafi,dqfi,dpfi  )

c          call WriteField('dufi',reshape(dufi,(/iip1,jmp1,llm/)))
c          call WriteField('dvfi',reshape(dvfi,(/iip1,jjm,llm/)))
c          call WriteField('dtetafi',reshape(dtetafi,(/iip1,jmp1,llm/)))
!         print*,'> APRES CALFIS (AVANT ADDFI)'
!         print*,''
!         print*,'teta(3049,:) = ',teta(3049,:)
!         print*,''
!         print*,'dteta(3049,:) = ',dteta(3049,:)
!         print*,''
!         print*,'dtetafi(3049,:) = ',dtetafi(3049,:)
!         print*,''


c      ajout des tendances physiques 
c      ------------------------------

          CALL addfi( dtphys, leapf, forward   ,
     $                  ucov, vcov, teta , q   ,ps ,
     $                 dufi, dvfi, dtetafi , dqfi ,dpfi  )

!         print*,'> APRES ADDFI'
!         print*,''
!         print*,'teta(3049,:) = ',teta(3049,:)
!         print*,''
!         print*,'dteta(3049,:) = ',dteta(3049,:)
!         print*,''
!         print*,'dtetafi(3049,:) = ',dtetafi(3049,:)
!         print*,''

          CALL pression (ip1jmp1,ap,bp,ps,p)
          CALL massdair(p,masse)
   
         ! Mass of tracers in each mess (kg) before Horizontal mixing of pressure
c        -------------------------------------------------------------------
          
         DO l=1, llm
           DO ij=1,ip1jmp1
              mq(ij,l) = masse(ij,l)*q(ij,l,igcm_co2)
           ENDDO
         ENDDO

c        Horizontal mixing of pressure
c        ------------------------------
         ! Rappel newtonien vers psmean
           psmean= globaverage2d(ps)  ! mean pressure
           p0=psmean/p00d
           DO ij=1,ip1jmp1
                oldps(ij)=ps(ij)
                ps(ij)=ps(ij) +(p0*kpd(ij)
     &                 -ps(ij))*(1-exp(-dtphys/tau_ps))
           ENDDO

           write(72,*) 'psmean = ',psmean  ! mean pressure redistributed

         ! security check: test pressure negative
           DO ij=1,ip1jmp1
             IF (ps(ij).le.0.) then
                 !write(*,*) 'Negative pressure :'
                 !write(*,*) 'ps = ',ps(ij),' ig = ',ij
                 !write(*,*) 'pmean = ',p0*kpd(ij)
                 !write(*,*) 'tau_ps = ',tau_ps
                 !STOP
                 ps(ij)=0.0000001*kpd(ij)/p00d
             ENDIF
           ENDDO
!***********************
!          Correction on teta due to surface pressure changes
           DO l = 1,llm
            DO ij = 1,ip1jmp1
              teta(ij,l)= teta(ij,l)*(ps(ij)/oldps(ij))**kappa
            ENDDO
           ENDDO
!***********************



!        ! update pressure and update p and pk
!         DO ij=1,ip1jmp1
!            ps(ij) = ps(ij) + dpfi(ij)*dtphys
!         ENDDO
          CALL pression (ip1jmp1,ap,bp,ps,p)
          CALL massdair(p,masse)
          if (pressure_exner) then
            CALL exner_hyb( ip1jmp1, ps, p, pks, pk, pkf )
          else
            CALL exner_milieu( ip1jmp1, ps, p, pks, pk, pkf )
          endif

         ! Update tracers after Horizontal mixing of pressure to ! conserve  tracer mass
c        -------------------------------------------------------------------
         DO l=1, llm
           DO ij=1,ip1jmp1
              q(ij,l,igcm_co2) = mq(ij,l)/ masse(ij,l)
           ENDDO
         ENDDO
          


c        Horizontal mixing of pressure
c        ------------------------------
         ! Rappel newtonien vers psmean
           psmean= globaverage2d(ps)  ! mean pressure
!        ! increment q_co2  with physical tendancy
!          IF (igcm_co2.ne.0) then
!            DO l=1, llm 
!               DO ij=1,ip1jmp1
!                q(ij,l,igcm_co2)=q(ij,l,igcm_co2)+
!    &                    dqfi(ij,l,igcm_co2)*dtphys
!               ENDDO
!            ENDDO
!          ENDIF

c          Mixing CO2 vertically
c          --------------------------
           if (igcm_co2.ne.0) then
            DO ij=1,ip1jmp1
               qmean_co2_vert=0.
               DO l=1, llm 
                 qmean_co2_vert= qmean_co2_vert
     &          + q(ij,l,igcm_co2)*( p(ij,l) - p(ij,l+1))
               END DO
               qmean_co2_vert= qmean_co2_vert/ps(ij)
               DO l=1, llm 
                 q(ij,l,igcm_co2)= qmean_co2_vert
               END DO
            END DO
           end if



c        Horizontal mixing of pressure
c        ------------------------------
         ! Rappel newtonien vers psmean
           psmean= globaverage2d(ps)  ! mean pressure

c        Horizontal mixing tracers, and temperature (replace dynamics in nogcm)
c        ---------------------------------------------------------------------------------  

         ! Simulate redistribution by dynamics for qco2 
           if (igcm_co2.ne.0) then

              DO ij=1,ip1jmp1
                 pqco2(ij)= ps(ij) * q(ij,1,igcm_co2)
              ENDDO
              pqco2mean=globaverage2d(pqco2) 

         !    Rappel newtonien vers qco2_mean
              qmean_co2= pqco2mean / psmean

              DO ij=1,ip1jmp1
                  q(ij,1,igcm_co2)=q(ij,1,igcm_co2)+
     &                  (qmean_co2-q(ij,1,igcm_co2))*
     &                  (1.-exp(-dtphys/tau_co2))
              ENDDO

              DO l=2, llm 
                 DO ij=1,ip1jmp1
                     q(ij,l,igcm_co2)=q(ij,1,igcm_co2)
                 END DO
              END DO

!             TEMPORAIRE (ED)
!             PRINT*,'psmean = ',psmean
!             PRINT*,'qmean_co2 = ',qmean_co2
!             PRINT*,'pqco2mean = ',pqco2mean
!             PRINT*,'q(50,1,igcm_co2) = ',q(50,1,igcm_co2)
!             PRINT*,'q(50,2,igcm_co2) = ',q(50,2,igcm_co2)
!             PRINT*,'q(50,3,igcm_co2) = ',q(50,3,igcm_co2)

           endif ! igcm_co2.ne.0


!       ********************************************************

c        Horizontal mixing of pressure
c        ------------------------------
         ! Rappel newtonien vers psmean
           psmean= globaverage2d(ps)  ! mean pressure


c          Mixing Temperature horizontally
c          -------------------------------
!          initialize variables that will be averaged
           DO l=1,llm
             DO ij=1,ip1jmp1
               dp(ij,l) = p(ij,l) - p(ij,l+1)
               tetadp(ij,l) = teta(ij,l)*dp(ij,l)
             ENDDO
           ENDDO 

           DO l=1,llm
             tetadpmean = globaverage2d(tetadp(:,l))
             dpmean = globaverage2d(dp(:,l))
             tetamean = tetadpmean / dpmean
             DO ij=1,ip1jmp1
               teta(ij,l) = teta(ij,l) + (tetamean - teta(ij,l)) *
     &                      (1 - exp(-dtphys/tau_teta))
             ENDDO
           ENDDO
           

       ENDIF ! of IF( apphys )

        
c   ********************************************************************
c   ********************************************************************
c   .... fin de l'integration dynamique  et physique pour le pas itau ..
c   ********************************************************************
c   ********************************************************************


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
               CALL geopot(ip1jmp1,teta,pk,pks,phis,phi)

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

            IF(itau.EQ.itaufin) THEN

              if (planet_type=="mars") then
                CALL dynredem1("restart.nc",REAL(itau)/REAL(day_step),
     &                         vcov,ucov,teta,q,masse,ps)
              else
                CALL dynredem1("restart.nc",start_time,
     &                         vcov,ucov,teta,q,masse,ps)
              endif
              CLOSE(99)
              !!! Ehouarn: Why not stop here and now?
            ENDIF ! of IF (itau.EQ.itaufin)

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
               CALL geopot(ip1jmp1,teta,pk,pks,phis,phi)

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

              IF(itau.EQ.itaufin) THEN
                if (planet_type=="mars") then
                  CALL dynredem1("restart.nc",REAL(itau)/REAL(day_step),
     &                         vcov,ucov,teta,q,masse,ps)
                else
                  CALL dynredem1("restart.nc",start_time,
     &                         vcov,ucov,teta,q,masse,ps)
                endif
              ENDIF ! of IF(itau.EQ.itaufin)

              forward = .TRUE.
              GO TO  1

            ENDIF ! of IF (forward)

      END IF ! of IF(.not.purmats)

      STOP
      END

! ************************************************************
!     globalaverage VERSION PLuto
! ************************************************************
      FUNCTION globaverage2d(var)
! ************************************************************
      IMPLICIT NONE
#include "dimensions.h"
#include "paramet.h"
#include "comgeom2.h"
       !ip1jmp1 called in paramet.h = iip1 x jjp1
       REAL var(iip1,jjp1), globaverage2d
       integer i,j
       REAL airetot
       save airetot
       logical firstcall
       data firstcall/.true./
       save firstcall

       if (firstcall) then
         firstcall=.false.
         airetot =0.
         DO j=2,jjp1-1      ! lat
           DO i = 1, iim    ! lon
             airetot = airetot + aire(i,j)
           END DO
         END DO
         DO i=1,iim
           airetot=airetot + aire(i,1)+aire(i,jjp1)
         ENDDO
       end if

       globaverage2d = 0.
       DO j=2,jjp1-1
         DO i = 1, iim
           globaverage2d = globaverage2d + var(i,j)*aire(i,j)
         END DO
       END DO
       
       DO i=1,iim
         globaverage2d = globaverage2d + var(i,1)*aire(i,1)
         globaverage2d = globaverage2d + var(i,jjp1)*aire(i,jjp1)
       ENDDO

       globaverage2d = globaverage2d/airetot
      return
      end

