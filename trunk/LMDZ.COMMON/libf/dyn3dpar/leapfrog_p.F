! 
! $Id: leapfrog_p.F 1446 2010-10-22 09:27:25Z emillour $
!
c
c

      SUBROUTINE leapfrog_p(ucov,vcov,teta,ps,masse,phis,q,
     &                    time_0)

      use exner_hyb_m, only: exner_hyb
      use exner_milieu_m, only: exner_milieu
      use exner_hyb_p_m, only: exner_hyb_p
      use exner_milieu_p_m, only: exner_milieu_p
       USE misc_mod
       USE parallel_lmdz
       USE times
       USE mod_hallo
       USE Bands
       USE Write_Field
       USE Write_Field_p
       USE vampir
       USE timer_filtre, ONLY : print_filtre_timer
       USE infotrac, ONLY: nqtot
       USE guide_p_mod, ONLY : guide_main
       USE getparam
       USE control_mod, only: planet_type,nday,day_step,iperiod,iphysiq,
     &                       less1day,fractday,ndynstep,iconser,
     &                       dissip_period,offline,ip_ebil_dyn,
     &                       ok_dynzon,periodav,ok_dyn_ave,iecri,
     &                       ok_dyn_ins,output_grads_dyn,
     &                       iapp_tracvl,ecritstart
       use cpdet_mod, only: cpdet,tpot2t_glo_p,t2tpot_glo_p
       use sponge_mod_p, only: callsponge,mode_sponge,sponge_p
       use comuforc_h
       USE comvert_mod, ONLY: ap,bp,pressure_exner,presnivs
       USE comconst_mod, ONLY: jmp1,daysec,dtvr,dtphys,dtdiss,
     .			cpp,ihf,iflag_top_bound,pi
       USE logic_mod, ONLY: iflag_phys,ok_guide,forward,leapf,apphys,
     .			statcl,conser,apdiss,purmats,tidal,ok_strato
       USE temps_mod, ONLY: itaufin,jD_ref,jH_ref,day_ini,
     .			day_ref,start_time,dt,hour_ini


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
      
      REAL,INTENT(IN) :: time_0 ! not used

c   dynamical variables:
      REAL,INTENT(INOUT) :: ucov(ip1jmp1,llm)    ! zonal covariant wind
      REAL,INTENT(INOUT) :: vcov(ip1jm,llm)      ! meridional covariant wind
      REAL,INTENT(INOUT) :: teta(ip1jmp1,llm)    ! potential temperature
      REAL,INTENT(INOUT) :: ps(ip1jmp1)          ! surface pressure (Pa)
      REAL,INTENT(INOUT) :: masse(ip1jmp1,llm)   ! air mass
      REAL,INTENT(INOUT) :: phis(ip1jmp1)        ! geopotentiat at the surface
      REAL,INTENT(INOUT) :: q(ip1jmp1,llm,nqtot) ! advected tracers
      
      REAL,SAVE :: p (ip1jmp1,llmp1  )       ! interlayer pressure
      REAL,SAVE :: pks(ip1jmp1)              ! exner at the surface
      REAL,SAVE :: pk(ip1jmp1,llm)           ! exner at mid-layer
      REAL,SAVE :: pkf(ip1jmp1,llm)          ! filtered exner at mid-layer
      REAL,SAVE :: phi(ip1jmp1,llm)          ! geopotential
      REAL,SAVE :: w(ip1jmp1,llm)            ! vertical velocity
! ADAPTATION GCM POUR CP(T)
      REAL,SAVE :: temp(ip1jmp1,llm)                 ! temperature  
      REAL,SAVE :: tsurpk(ip1jmp1,llm)               ! cpp*T/pk  

      real zqmin,zqmax

c variables dynamiques intermediaire pour le transport
      REAL,SAVE :: pbaru(ip1jmp1,llm),pbarv(ip1jm,llm) !flux de masse

c   variables dynamiques au pas -1
      REAL,SAVE :: vcovm1(ip1jm,llm),ucovm1(ip1jmp1,llm)
      REAL,SAVE :: tetam1(ip1jmp1,llm),psm1(ip1jmp1)
      REAL,SAVE :: massem1(ip1jmp1,llm)

c   tendances dynamiques
      REAL,SAVE :: dv(ip1jm,llm),du(ip1jmp1,llm)
      REAL,SAVE :: dteta(ip1jmp1,llm),dp(ip1jmp1)
      REAL,DIMENSION(:,:,:), ALLOCATABLE, SAVE :: dq

c   tendances de la dissipation
      REAL,SAVE :: dvdis(ip1jm,llm),dudis(ip1jmp1,llm)
      REAL,SAVE :: dtetadis(ip1jmp1,llm)

c   tendances physiques
      REAL,SAVE :: dvfi(ip1jm,llm),dufi(ip1jmp1,llm)
      REAL,SAVE :: dtetafi(ip1jmp1,llm)
      REAL,SAVE :: dpfi(ip1jmp1)
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: dqfi

c   tendances top_bound (sponge layer)
c      REAL,SAVE :: dvtop(ip1jm,llm)
      REAL,SAVE :: dutop(ip1jmp1,llm)
c      REAL,SAVE :: dtetatop(ip1jmp1,llm)
c      REAL,SAVE :: dptop(ip1jmp1)
c      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: dqtop

c   TITAN : tendances due au forces de marees */s
      REAL,SAVE :: dvtidal(ip1jm,llm),dutidal(ip1jmp1,llm)

c   variables pour le fichier histoire
      REAL dtav      ! intervalle de temps elementaire
      LOGICAL lrestart

      REAL tppn(iim),tpps(iim),tpn,tps
c
      INTEGER itau,itaufinp1,iav
!      INTEGER  iday ! jour julien
      REAL       time 

      REAL  SSUM 
!      REAL,SAVE :: finvmaold(ip1jmp1,llm)

cym      LOGICAL  lafin
      LOGICAL :: lafin
      INTEGER ij,iq,l
      INTEGER ik

      real time_step, t_wrt, t_ops

! jD_cur: jour julien courant
! jH_cur: heure julienne courante
      REAL :: jD_cur, jH_cur
      INTEGER :: an, mois, jour
      REAL :: secondes
      real :: rdaym_ini
      logical :: physics
      LOGICAL first,callinigrads

      data callinigrads/.true./
      character*10 string10

      REAL,SAVE :: flxw(ip1jmp1,llm) ! flux de masse verticale

c+jld variables test conservation energie
      REAL,SAVE :: ecin(ip1jmp1,llm),ecin0(ip1jmp1,llm)
C     Tendance de la temp. potentiel d (theta)/ d t due a la 
C     tansformation d'energie cinetique en energie thermique
C     cree par la dissipation
      REAL,SAVE :: dtetaecdt(ip1jmp1,llm)
      REAL,SAVE :: vcont(ip1jm,llm),ucont(ip1jmp1,llm)
      REAL,SAVE :: vnat(ip1jm,llm),unat(ip1jmp1,llm)
      REAL      d_h_vcol, d_qt, d_qw, d_ql, d_ec
      CHARACTER*15 ztit
!      INTEGER   ip_ebil_dyn  ! PRINT level for energy conserv. diag.
!      SAVE      ip_ebil_dyn
!      DATA      ip_ebil_dyn/0/
c-jld 

      character*80 dynhist_file, dynhistave_file
      character(len=*),parameter :: modname="leapfrog"
      character*80 abort_message


      logical,PARAMETER :: dissip_conservative=.TRUE.
 
      INTEGER testita
      PARAMETER (testita = 9)

      logical , parameter :: flag_verif = .false.

      ! for CP(T)  -- Aymeric
      real :: dtec
      real,save :: ztetaec(ip1jmp1,llm)  !!SAVE ???

c declaration liees au parallelisme
      INTEGER :: ierr
      LOGICAL :: FirstCaldyn
      LOGICAL :: FirstPhysic
      INTEGER :: ijb,ije,j,i
      type(Request) :: TestRequest
      type(Request) :: Request_Dissip
      type(Request) :: Request_physic
      REAL,SAVE :: dvfi_tmp(iip1,llm),dufi_tmp(iip1,llm)
      REAL,SAVE :: dtetafi_tmp(iip1,llm)
      REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: dqfi_tmp
      REAL,SAVE :: dpfi_tmp(iip1)

      INTEGER :: true_itau
      INTEGER :: iapptrac
      INTEGER :: AdjustCount
!      INTEGER :: var_time
      LOGICAL :: ok_start_timer=.FALSE.
      LOGICAL, SAVE :: firstcall=.TRUE.

c$OMP MASTER
      ItCount=0
c$OMP END MASTER      
      true_itau=0
      FirstCaldyn=.TRUE.
      FirstPhysic=.TRUE.
      iapptrac=0
      AdjustCount = 0
      lafin=.false.
      
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

      itau = 0
      physics=.true.
      if (iflag_phys==0.or.iflag_phys==2) physics=.false.
!      iday = day_ini+itau/day_step
!      time = REAL(itau-(iday-day_ini)*day_step)/day_step+time_0
!         IF(time.GT.1.) THEN
!          time = time-1.
!          iday = iday+1
!         ENDIF

c Allocate variables depending on dynamic variable nqtot
c$OMP MASTER
         IF (firstcall) THEN
            firstcall=.FALSE.
            ALLOCATE(dq(ip1jmp1,llm,nqtot))
            ALLOCATE(dqfi(ip1jmp1,llm,nqtot))
            ALLOCATE(dqfi_tmp(iip1,llm,nqtot))
c            ALLOCATE(dqtop(ip1jmp1,llm,nqtot))
         END IF
c$OMP END MASTER      
c$OMP BARRIER

c-----------------------------------------------------------------------
c   On initialise la pression et la fonction d'Exner :
c   --------------------------------------------------

c$OMP MASTER
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
        dq(:,:,:)   =0.
        dp(:)=0
        pbaru(:,:)=0
        pbarv(:,:)=0

      CALL pression ( ip1jmp1, ap, bp, ps, p       )
      if (pressure_exner) then
        CALL exner_hyb( ip1jmp1, ps, p, pks, pk, pkf )
      else
        CALL exner_milieu( ip1jmp1, ps, p, pks, pk, pkf )
      endif
c$OMP END MASTER
c-----------------------------------------------------------------------
c   Debut de l'integration temporelle:
c   ----------------------------------
c et du parallelisme !!

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
      if (jH_cur > 1.0 ) then
        jD_cur = jD_cur +1.
        jH_cur = jH_cur -1.
      endif


#ifdef CPP_IOIPSL
      IF (planet_type.eq."earth") THEN
      if (ok_guide) then
!$OMP MASTER
        call guide_main(itau,ucov,vcov,teta,q,masse,ps)
!$OMP END MASTER
!$OMP BARRIER
      endif
      ENDIF
#endif

c
c     IF( MOD( itau, 10* day_step ).EQ.0 )  THEN
c       CALL  test_period ( ucov,vcov,teta,q,p,phis )
c       PRINT *,' ----   Test_period apres continue   OK ! -----', itau
c     ENDIF 
c
cym      CALL SCOPY( ijmllm ,vcov , 1, vcovm1 , 1 )
cym      CALL SCOPY( ijp1llm,ucov , 1, ucovm1 , 1 )
cym      CALL SCOPY( ijp1llm,teta , 1, tetam1 , 1 )
cym      CALL SCOPY( ijp1llm,masse, 1, massem1, 1 )
cym      CALL SCOPY( ip1jmp1, ps  , 1,   psm1 , 1 )

       if (FirstCaldyn) then
c$OMP MASTER
         ucovm1=ucov
         vcovm1=vcov
         tetam1= teta
         massem1= masse
         psm1= ps
         
! Ehouarn: finvmaold is actually not used       
!         finvmaold = masse
!         CALL filtreg ( finvmaold ,jjp1, llm, -2,2, .TRUE., 1 )
c$OMP END MASTER
c$OMP BARRIER
       else
! Save fields obtained at previous time step as '...m1'
         ijb=ij_begin
         ije=ij_end

c$OMP MASTER           
         psm1     (ijb:ije) = ps    (ijb:ije)
c$OMP END MASTER

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)         
         DO l=1,llm      
           ije=ij_end
           ucovm1   (ijb:ije,l) = ucov  (ijb:ije,l)
           tetam1   (ijb:ije,l) = teta  (ijb:ije,l)
           massem1  (ijb:ije,l) = masse (ijb:ije,l)
!           finvmaold(ijb:ije,l)=masse(ijb:ije,l)
                 
           if (pole_sud) ije=ij_end-iip1
           vcovm1(ijb:ije,l) = vcov  (ijb:ije,l)
       

         ENDDO
c$OMP ENDDO  


! Ehouarn: finvmaold not used
!          CALL filtreg_p ( finvmaold ,jj_begin,jj_end,jjp1, 
!     .                    llm, -2,2, .TRUE., 1 )

       endif ! of if (FirstCaldyn)
       
      forward = .TRUE.
      leapf   = .FALSE.
      dt      =  dtvr

c   ...    P.Le Van .26/04/94  ....

cym      CALL SCOPY   ( ijp1llm,   masse, 1, finvmaold,     1 )
cym      CALL filtreg ( finvmaold ,jjp1, llm, -2,2, .TRUE., 1 )

cym  ne sert a rien
cym      call minmax(ijp1llm,q(:,:,3),zqmin,zqmax)

   2  CONTINUE ! Matsuno backward or leapfrog step begins here

c$OMP MASTER
      ItCount=ItCount+1
      if (MOD(ItCount,1)==1) then
        debug=.true.
      else
        debug=.false.
      endif
c$OMP END MASTER
c-----------------------------------------------------------------------

c   date:  (NB: only leapfrog step requires recomputing date)
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
        if (jH_cur > 1.0 ) then
          jD_cur = jD_cur +1.
          jH_cur = jH_cur -1.
        endif
      ENDIF


c   gestion des appels de la physique et des dissipations:
c   ------------------------------------------------------
c
c   ...    P.Le Van  ( 6/02/95 )  ....

      apphys = .FALSE.
      statcl = .FALSE.
      conser = .FALSE.
      apdiss = .FALSE.
c      idissip=1
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
         IF( MOD(itau+1,iphysiq).EQ.0.AND.physics) apphys=.TRUE.
      END IF

! Ehouarn: for Shallow Water case (ie: 1 vertical layer),
!          supress dissipation step
      if (llm.eq.1) then
        apdiss=.false.
      endif

#ifdef NODYN
      apdiss=.false.
#endif


cym    ---> Pour le moment      
cym      apphys = .FALSE.
      statcl = .FALSE.
      conser = .FALSE. ! ie: no output of control variables to stdout in //
      
      if (firstCaldyn) then
c$OMP MASTER
          call SetDistrib(jj_Nb_Caldyn)
c$OMP END MASTER
c$OMP BARRIER
          firstCaldyn=.FALSE.
cym          call InitTime
c$OMP MASTER
          call Init_timer
c$OMP END MASTER
      endif

c$OMP MASTER      
      IF (ok_start_timer) THEN
        CALL InitTime
        ok_start_timer=.FALSE.
      ENDIF      
c$OMP END MASTER      
     
      if (Adjust) then
c$OMP MASTER 
        AdjustCount=AdjustCount+1
        if (iapptrac==iapp_tracvl .and. (forward. OR . leapf)
     &         .and. itau/iphysiq>2 .and. Adjustcount>30) then
           AdjustCount=0
           call allgather_timer_average

        if (prt_level > 9) then
        
        print *,'*********************************'
        print *,'******    TIMER CALDYN     ******'
        do i=0,mpi_size-1
          print *,'proc',i,' :   Nb Bandes  :',jj_nb_caldyn(i),
     &            '  : temps moyen :',
     &             timer_average(jj_nb_caldyn(i),timer_caldyn,i),
     &            '+-',timer_delta(jj_nb_caldyn(i),timer_caldyn,i)
        enddo
      
        print *,'*********************************'
        print *,'******    TIMER VANLEER    ******'
        do i=0,mpi_size-1
          print *,'proc',i,' :   Nb Bandes  :',jj_nb_vanleer(i),
     &            '  : temps moyen :',
     &             timer_average(jj_nb_vanleer(i),timer_vanleer,i),
     &            '+-',timer_delta(jj_nb_vanleer(i),timer_vanleer,i)
        enddo
      
        print *,'*********************************'
        print *,'******    TIMER DISSIP    ******'
        do i=0,mpi_size-1
          print *,'proc',i,' :   Nb Bandes  :',jj_nb_dissip(i),
     &            '  : temps moyen :',
     &             timer_average(jj_nb_dissip(i),timer_dissip,i),
     &             '+-',timer_delta(jj_nb_dissip(i),timer_dissip,i)
        enddo
        
        if (mpi_rank==0) call WriteBands
        
       endif
       
         call AdjustBands_caldyn
         if (mpi_rank==0) call WriteBands
         
         call Register_SwapFieldHallo(ucov,ucov,ip1jmp1,llm,
     &                                jj_Nb_caldyn,0,0,TestRequest)
         call Register_SwapFieldHallo(ucovm1,ucovm1,ip1jmp1,llm,
     &                                jj_Nb_caldyn,0,0,TestRequest)
         call Register_SwapFieldHallo(vcov,vcov,ip1jm,llm,
     &                                jj_Nb_caldyn,0,0,TestRequest)
         call Register_SwapFieldHallo(vcovm1,vcovm1,ip1jm,llm,
     &                                jj_Nb_caldyn,0,0,TestRequest)
         call Register_SwapFieldHallo(teta,teta,ip1jmp1,llm,
     &                                jj_Nb_caldyn,0,0,TestRequest)
         call Register_SwapFieldHallo(tetam1,tetam1,ip1jmp1,llm,
     &                                jj_Nb_caldyn,0,0,TestRequest)
         call Register_SwapFieldHallo(masse,masse,ip1jmp1,llm,
     &                                jj_Nb_caldyn,0,0,TestRequest)
         call Register_SwapFieldHallo(massem1,massem1,ip1jmp1,llm,
     &                                jj_Nb_caldyn,0,0,TestRequest)
         call Register_SwapFieldHallo(ps,ps,ip1jmp1,1,
     &                                jj_Nb_caldyn,0,0,TestRequest)
         call Register_SwapFieldHallo(psm1,psm1,ip1jmp1,1,
     &                                jj_Nb_caldyn,0,0,TestRequest)
         call Register_SwapFieldHallo(pkf,pkf,ip1jmp1,llm,
     &                                jj_Nb_caldyn,0,0,TestRequest)
         call Register_SwapFieldHallo(pk,pk,ip1jmp1,llm,
     &                                jj_Nb_caldyn,0,0,TestRequest)
         call Register_SwapFieldHallo(pks,pks,ip1jmp1,1,
     &                                jj_Nb_caldyn,0,0,TestRequest)
         call Register_SwapFieldHallo(phis,phis,ip1jmp1,1,
     &                                jj_Nb_caldyn,0,0,TestRequest)
         call Register_SwapFieldHallo(phi,phi,ip1jmp1,llm,
     &                                jj_Nb_caldyn,0,0,TestRequest)
!         call Register_SwapFieldHallo(finvmaold,finvmaold,ip1jmp1,llm,
!     &                                jj_Nb_caldyn,0,0,TestRequest)

        do j=1,nqtot
         call Register_SwapFieldHallo(q(1,1,j),q(1,1,j),ip1jmp1,llm,
     &                                jj_nb_caldyn,0,0,TestRequest)
        enddo
! ADAPTATION GCM POUR CP(T)
         call Register_SwapFieldHallo(temp,temp,ip1jmp1,llm,
     &                                jj_Nb_caldyn,0,0,TestRequest)
         call Register_SwapFieldHallo(tsurpk,tsurpk,ip1jmp1,llm,
     &                                jj_Nb_caldyn,0,0,TestRequest)

         call SetDistrib(jj_nb_caldyn)
         call SendRequest(TestRequest)
         call WaitRequest(TestRequest)
         
        call AdjustBands_dissip
        call AdjustBands_physic

      endif
c$OMP END MASTER  
      endif ! of if (Adjust)
     
      
      
c-----------------------------------------------------------------------
c   calcul des tendances dynamiques:
c   --------------------------------
! ADAPTATION GCM POUR CP(T)
      call tpot2t_glo_p(teta,temp,pk)
      ijb=ij_begin
      ije=ij_end
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      do l=1,llm
        tsurpk(ijb:ije,l)=cpp*temp(ijb:ije,l)/pk(ijb:ije,l)
      enddo
!$OMP END DO

      if (debug) then
!$OMP BARRIER
!$OMP MASTER
        call WriteField_p('temp',reshape(temp,(/iip1,jmp1,llm/)))
        call WriteField_p('tsurpk',reshape(tsurpk,(/iip1,jmp1,llm/)))
!$OMP END MASTER        
!$OMP BARRIER     
      endif ! of if (debug)
      
c$OMP BARRIER
c$OMP MASTER
       call VTb(VThallo)
c$OMP END MASTER

       call Register_Hallo(ucov,ip1jmp1,llm,1,1,1,1,TestRequest)
       call Register_Hallo(vcov,ip1jm,llm,1,1,1,1,TestRequest)
       call Register_Hallo(teta,ip1jmp1,llm,1,1,1,1,TestRequest)
       call Register_Hallo(ps,ip1jmp1,1,1,2,2,1,TestRequest)
       call Register_Hallo(pkf,ip1jmp1,llm,1,1,1,1,TestRequest)
       call Register_Hallo(pk,ip1jmp1,llm,1,1,1,1,TestRequest)
       call Register_Hallo(pks,ip1jmp1,1,1,1,1,1,TestRequest)
       call Register_Hallo(p,ip1jmp1,llmp1,1,1,1,1,TestRequest)
! ADAPTATION GCM POUR CP(T)
       call Register_Hallo(temp,ip1jmp1,llm,1,1,1,1,TestRequest)
       call Register_Hallo(tsurpk,ip1jmp1,llm,1,1,1,1,TestRequest)
       
c       do j=1,nqtot
c         call Register_Hallo(q(1,1,j),ip1jmp1,llm,1,1,1,1,
c     *                       TestRequest)
c        enddo

       call SendRequest(TestRequest)
c$OMP BARRIER
       call WaitRequest(TestRequest)

c$OMP MASTER
       call VTe(VThallo)
c$OMP END MASTER
c$OMP BARRIER
      
      if (debug) then        
!$OMP BARRIER
!$OMP MASTER
        call WriteField_p('ucov',reshape(ucov,(/iip1,jmp1,llm/)))
        call WriteField_p('vcov',reshape(vcov,(/iip1,jjm,llm/)))
        call WriteField_p('teta',reshape(teta,(/iip1,jmp1,llm/)))
        call WriteField_p('ps',reshape(ps,(/iip1,jmp1/)))
        call WriteField_p('masse',reshape(masse,(/iip1,jmp1,llm/)))
        call WriteField_p('pk',reshape(pk,(/iip1,jmp1,llm/)))
        call WriteField_p('pks',reshape(pks,(/iip1,jmp1/)))
        call WriteField_p('pkf',reshape(pkf,(/iip1,jmp1,llm/)))
        call WriteField_p('phis',reshape(phis,(/iip1,jmp1/)))
        if (nqtot > 0) then
        do j=1,nqtot
          call WriteField_p('q'//trim(int2str(j)),
     .                reshape(q(:,:,j),(/iip1,jmp1,llm/)))
        enddo
        endif
!$OMP END MASTER        
c$OMP BARRIER
      endif

      
      True_itau=True_itau+1

c$OMP MASTER
      IF (prt_level>9) THEN
        WRITE(lunout,*)"leapfrog_p: Iteration No",True_itau
      ENDIF


      call start_timer(timer_caldyn)

      ! compute geopotential phi()
! ADAPTATION GCM POUR CP(T)
!      CALL geopot_p  ( ip1jmp1, teta  , pk , pks,  phis  , phi   )
      CALL geopot_p  ( ip1jmp1, tsurpk  , pk , pks,  phis  , phi   )
      
      call VTb(VTcaldyn)
c$OMP END MASTER
!      var_time=time+iday-day_ini

c$OMP BARRIER
!      CALL FTRACE_REGION_BEGIN("caldyn")
      time = jD_cur + jH_cur 
           rdaym_ini  = itau * dtvr / daysec

#ifdef NODYN
      !WRITE(lunout,*)"NO DYN !!!!!"
      dv(:,:) = 0.D+0
      du(:,:) = 0.D+0
      dteta(:,:) = 0.D+0
      dq(:,:,:) = 0.D+0
      dp(:) = 0.D+0
      if (planet_type.eq."generic") then
      if (ok_guide) then
       DO l=1,llm
       ucov(:,l) = ucov(:,l) + dt*attenua(l)*
     &   ((uforc(:,l)-ucov(:,l))/facwind)
       ENDDO
      endif
      endif
#else
! ADAPTATION GCM POUR CP(T)
!      CALL caldyn_p 
!     $  ( itau,ucov,vcov,teta,ps,masse,pk,pkf,phis ,
!     $    phi,conser,du,dv,dteta,dp,w, pbaru,pbarv, time )
      CALL caldyn_p 
     $  ( itau,ucov,vcov,teta,ps,masse,pk,pkf,tsurpk,phis,
     $    phi,conser,du,dv,dteta,dp,w, pbaru,pbarv, time )

!      CALL FTRACE_REGION_END("caldyn")

c$OMP MASTER
      call VTe(VTcaldyn)
c$OMP END MASTER      

cc$OMP BARRIER
cc$OMP MASTER
!      call WriteField_p('du',reshape(du,(/iip1,jmp1,llm/)))
!      call WriteField_p('dv',reshape(dv,(/iip1,jjm,llm/)))
!      call WriteField_p('dteta',reshape(dteta,(/iip1,jmp1,llm/)))
!      call WriteField_p('dp',reshape(dp,(/iip1,jmp1/)))
!      call WriteField_p('w',reshape(w,(/iip1,jmp1,llm/)))
!      call WriteField_p('pbaru',reshape(pbaru,(/iip1,jmp1,llm/)))
!      call WriteField_p('pbarv',reshape(pbarv,(/iip1,jjm,llm/)))
!      call WriteField_p('p',reshape(p,(/iip1,jmp1,llmp1/)))
!      call WriteField_p('masse',reshape(masse,(/iip1,jmp1,llm/)))
!      call WriteField_p('pk',reshape(pk,(/iip1,jmp1,llm/)))
cc$OMP END MASTER


      ! Simple zonal wind nudging for generic planetary model
      ! AS 09/2013
      ! ---------------------------------------------------
      if (planet_type.eq."generic") then
       if (ok_guide) then
         DO l=1,llm
          du(:,l)=du(:,l)+attenua(l)*((uforc(:,l)-ucov(:,l))/facwind) 
         ENDDO
       endif
      endif

c-----------------------------------------------------------------------
c   calcul des tendances advection des traceurs (dont l'humidite)
c   -------------------------------------------------------------

      IF( forward. OR . leapf )  THEN
! Ehouarn: NB: fields sent to advtrac are those at the beginning of the time step
         CALL advtrac_p( pbaru,pbarv, 
     *             p,  masse,q,iapptrac, teta,
     .             flxw, pk)

C        Stokage du flux de masse pour traceurs OFF-LINE
         IF (offline .AND. .NOT. adjust) THEN
            CALL fluxstokenc_p(pbaru,pbarv,masse,teta,phi,phis,
     .           dtvr, itau)
         ENDIF

      ENDIF ! of IF( forward. OR . leapf )

c-----------------------------------------------------------------------
c   integrations dynamique et traceurs:
c   ----------------------------------

c$OMP MASTER 
       call VTb(VTintegre)
c$OMP END MASTER
c      call WriteField_p('ucovm1',reshape(ucovm1,(/iip1,jmp1,llm/)))
c      call WriteField_p('vcovm1',reshape(vcovm1,(/iip1,jjm,llm/)))
c      call WriteField_p('tetam1',reshape(tetam1,(/iip1,jmp1,llm/)))
c      call WriteField_p('psm1',reshape(psm1,(/iip1,jmp1/)))
c      call WriteField_p('ucov',reshape(ucov,(/iip1,jmp1,llm/)))
c      call WriteField_p('vcov',reshape(vcov,(/iip1,jjm,llm/)))
c      call WriteField_p('teta',reshape(teta,(/iip1,jmp1,llm/)))
c      call WriteField_p('ps',reshape(ps,(/iip1,jmp1/)))
cc$OMP PARALLEL DEFAULT(SHARED)
c$OMP BARRIER
!       CALL FTRACE_REGION_BEGIN("integrd")

       CALL integrd_p (nqtot,vcovm1,ucovm1,tetam1,psm1,massem1 ,
     $         dv,du,dteta,dq,dp,vcov,ucov,teta,q,ps,masse,phis )
!     $              finvmaold                                    )

       IF ((planet_type.eq."titan").and.(tidal)) then
c-----------------------------------------------------------------------
c   Marées gravitationnelles causées par Saturne
c   B. Charnay (28/10/2010)
c   ----------------------------------------------------------
            CALL tidal_forces(rdaym_ini, dutidal, dvtidal)
            ucov=ucov+dutidal*dt
            vcov=vcov+dvtidal*dt
       ENDIF

!       CALL FTRACE_REGION_END("integrd")
c$OMP BARRIER
cc$OMP MASTER
c      call WriteField_p('ucovm1',reshape(ucovm1,(/iip1,jmp1,llm/)))
c      call WriteField_p('vcovm1',reshape(vcovm1,(/iip1,jjm,llm/)))
c      call WriteField_p('tetam1',reshape(tetam1,(/iip1,jmp1,llm/)))
c      call WriteField_p('psm1',reshape(psm1,(/iip1,jmp1/)))
c      call WriteField_p('ucov',reshape(ucov,(/iip1,jmp1,llm/)))
c      call WriteField_p('vcov',reshape(vcov,(/iip1,jjm,llm/)))
c      call WriteField_p('teta',reshape(teta,(/iip1,jmp1,llm/)))
c      call WriteField_p('dteta',reshape(dteta,(/iip1,jmp1,llm/)))
c
c      call WriteField_p('ps',reshape(ps,(/iip1,jmp1/)))
c      do j=1,nqtot
c        call WriteField_p('q'//trim(int2str(j)),
c     .                reshape(q(:,:,j),(/iip1,jmp1,llm/)))
c        call WriteField_p('dq'//trim(int2str(j)),
c     .                reshape(dq(:,:,j),(/iip1,jmp1,llm/)))
c      enddo
cc$OMP END MASTER

! NODYN precompiling flag
#endif

c$OMP MASTER 
       call VTe(VTintegre)
c$OMP END MASTER
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

cc$OMP END PARALLEL

c
c
       IF( apphys )  THEN
c
c     .......   Ajout   P.Le Van ( 17/04/96 )   ...........
c
cc$OMP PARALLEL DEFAULT(SHARED)
cc$OMP+         PRIVATE(rdaym_ini,rdayvrai,ijb,ije)

c$OMP MASTER
         call suspend_timer(timer_caldyn)

        if (prt_level >= 10) then
         write(lunout,*)
     &   'leapfrog_p: Entree dans la physique : Iteration No ',true_itau
        endif
c$OMP END MASTER

         CALL pression_p (  ip1jmp1, ap, bp, ps,  p      )

c$OMP BARRIER
         if (pressure_exner) then
           CALL exner_hyb_p(  ip1jmp1, ps, p,pks, pk, pkf )
         else
           CALL exner_milieu_p( ip1jmp1, ps, p, pks, pk, pkf )
         endif
c$OMP BARRIER
! Compute geopotential (physics might need it)
!====
! GEOP CORRECTION
! ADAPTATION GCM POUR CP(T)
         call tpot2t_glo_p(teta,temp,pk)
         ijb=ij_begin
         ije=ij_end
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
         do l=1,llm
           tsurpk(ijb:ije,l)=cpp*temp(ijb:ije,l)/pk(ijb:ije,l)
         enddo
!$OMP END DO
c$OMP MASTER
!         CALL geopot_p  ( ip1jmp1, teta  , pk , pks,  phis  , phi   )
         CALL geopot_p( ip1jmp1, tsurpk, pk, pks, phis, phi )
c$OMP END MASTER
c$OMP BARRIER
!====

           jD_cur = jD_ref + day_ini - day_ref
     $        + (itau+1)/day_step

           IF ((planet_type .eq."generic").or.
     &         (planet_type.eq."mars")) THEN
              ! AS: we make jD_cur to be pday
              jD_cur = int(day_ini + itau/day_step)
           ENDIF

           IF (planet_type .eq. "mars") THEN
             jH_cur = jH_ref + hour_ini +                               &
     &                mod(itau,day_step)/float(day_step)
           ELSE IF (planet_type .eq. "generic") THEN
             jH_cur = jH_ref + start_time +                             &
     &                mod(itau,day_step)/float(day_step)
           ELSE
             jH_cur = jH_ref + start_time +                             &
     &                mod(itau+1,day_step)/float(day_step)
           ENDIF
           if (jH_cur > 1.0 ) then
             jD_cur = jD_cur +1.
             jH_cur = jH_cur -1.
           endif

c rajout debug
c       lafin = .true.


c   Interface avec les routines de phylmd (phymars ... )
c   -----------------------------------------------------

c+jld

c  Diagnostique de conservation de l'energie : initialisation
      IF (ip_ebil_dyn.ge.1 ) THEN 
          ztit='bil dyn'
! Ehouarn: be careful, diagedyn is Earth-specific!
           IF (planet_type.eq."earth") THEN
            CALL diagedyn(ztit,2,1,1,dtphys
     &    , ucov    , vcov , ps, p ,pk , teta , q(:,:,1), q(:,:,2))
           ENDIF
      ENDIF 
c-jld
c$OMP BARRIER
c$OMP MASTER
        call VTb(VThallo)
c$OMP END MASTER

        call SetTag(Request_physic,800)
        
        call Register_SwapFieldHallo(ucov,ucov,ip1jmp1,llm,
     *                               jj_Nb_physic,2,2,Request_physic)
        
        call Register_SwapFieldHallo(vcov,vcov,ip1jm,llm,
     *                               jj_Nb_physic,2,2,Request_physic)
        
        call Register_SwapFieldHallo(teta,teta,ip1jmp1,llm,
     *                               jj_Nb_physic,2,2,Request_physic)
        
        call Register_SwapFieldHallo(masse,masse,ip1jmp1,llm,
     *                               jj_Nb_physic,1,2,Request_physic)

        call Register_SwapFieldHallo(ps,ps,ip1jmp1,1,
     *                               jj_Nb_physic,2,2,Request_physic)

        call Register_SwapFieldHallo(p,p,ip1jmp1,llmp1,
     *                               jj_Nb_physic,2,2,Request_physic)
        
        call Register_SwapFieldHallo(pk,pk,ip1jmp1,llm,
     *                               jj_Nb_physic,2,2,Request_physic)
        
        call Register_SwapFieldHallo(phis,phis,ip1jmp1,1,
     *                               jj_Nb_physic,2,2,Request_physic)
        
        call Register_SwapFieldHallo(phi,phi,ip1jmp1,llm,
     *                               jj_Nb_physic,2,2,Request_physic)
        
        call Register_SwapFieldHallo(w,w,ip1jmp1,llm,
     *                               jj_Nb_physic,2,2,Request_physic)
        
c        call SetDistrib(jj_nb_vanleer)
        do j=1,nqtot
 
          call Register_SwapFieldHallo(q(1,1,j),q(1,1,j),ip1jmp1,llm,
     *                               jj_Nb_physic,2,2,Request_physic)
        enddo

        call Register_SwapFieldHallo(flxw,flxw,ip1jmp1,llm,
     *                               jj_Nb_physic,2,2,Request_physic)
        
        call SendRequest(Request_Physic)
c$OMP BARRIER
        call WaitRequest(Request_Physic)       

c$OMP BARRIER
c$OMP MASTER
        call SetDistrib(jj_nb_Physic)
        call VTe(VThallo)
        
        call VTb(VTphysiq)
c$OMP END MASTER
c$OMP BARRIER

cc$OMP MASTER        
c      call WriteField_p('ucovfi',reshape(ucov,(/iip1,jmp1,llm/)))
c      call WriteField_p('vcovfi',reshape(vcov,(/iip1,jjm,llm/)))
c      call WriteField_p('tetafi',reshape(teta,(/iip1,jmp1,llm/)))
c      call WriteField_p('pfi',reshape(p,(/iip1,jmp1,llmp1/)))
c      call WriteField_p('pkfi',reshape(pk,(/iip1,jmp1,llm/)))
cc$OMP END MASTER
cc$OMP BARRIER
!        CALL FTRACE_REGION_BEGIN("calfis")
        CALL calfis_p(lafin ,jD_cur, jH_cur,
     $               ucov,vcov,teta,q,masse,ps,p,pk,phis,phi ,
     $               du,dv,dteta,dq,
     $               flxw,
     $               dufi,dvfi,dtetafi,dqfi,dpfi  )
!        CALL FTRACE_REGION_END("calfis")
        ijb=ij_begin
        ije=ij_end  
        if ( .not. pole_nord) then
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l=1,llm
          dufi_tmp(1:iip1,l)   = dufi(ijb:ijb+iim,l) 
          dvfi_tmp(1:iip1,l)   = dvfi(ijb:ijb+iim,l)  
          dtetafi_tmp(1:iip1,l)= dtetafi(ijb:ijb+iim,l)  
          dqfi_tmp(1:iip1,l,:) = dqfi(ijb:ijb+iim,l,:)  
          ENDDO
c$OMP END DO NOWAIT

c$OMP MASTER
          dpfi_tmp(1:iip1)     = dpfi(ijb:ijb+iim)  
c$OMP END MASTER
        endif ! of if ( .not. pole_nord)

c$OMP BARRIER
c$OMP MASTER
        call SetDistrib(jj_nb_Physic_bis)

        call VTb(VThallo)
c$OMP END MASTER
c$OMP BARRIER
 
        call Register_Hallo(dufi,ip1jmp1,llm,
     *                      1,0,0,1,Request_physic)
        
        call Register_Hallo(dvfi,ip1jm,llm,
     *                      1,0,0,1,Request_physic)
        
        call Register_Hallo(dtetafi,ip1jmp1,llm,
     *                      1,0,0,1,Request_physic)

        call Register_Hallo(dpfi,ip1jmp1,1,
     *                      1,0,0,1,Request_physic)

        if (nqtot > 0) then
        do j=1,nqtot
          call Register_Hallo(dqfi(1,1,j),ip1jmp1,llm,
     *                        1,0,0,1,Request_physic)
        enddo
        endif
        
        call SendRequest(Request_Physic)
c$OMP BARRIER
        call WaitRequest(Request_Physic)
             
c$OMP BARRIER
c$OMP MASTER
        call VTe(VThallo)
 
        call SetDistrib(jj_nb_Physic)
c$OMP END MASTER
c$OMP BARRIER        
                ijb=ij_begin
        if (.not. pole_nord) then
        
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l=1,llm
            dufi(ijb:ijb+iim,l) = dufi(ijb:ijb+iim,l)+dufi_tmp(1:iip1,l)
            dvfi(ijb:ijb+iim,l) = dvfi(ijb:ijb+iim,l)+dvfi_tmp(1:iip1,l) 
            dtetafi(ijb:ijb+iim,l) = dtetafi(ijb:ijb+iim,l)
     &                              +dtetafi_tmp(1:iip1,l)
            dqfi(ijb:ijb+iim,l,:) = dqfi(ijb:ijb+iim,l,:)
     &                              + dqfi_tmp(1:iip1,l,:)
          ENDDO
c$OMP END DO NOWAIT

c$OMP MASTER
          dpfi(ijb:ijb+iim)   = dpfi(ijb:ijb+iim)+ dpfi_tmp(1:iip1)
c$OMP END MASTER
          
        endif ! of if (.not. pole_nord)
c$OMP BARRIER
cc$OMP MASTER        
c      call WriteField_p('dufi',reshape(dufi,(/iip1,jmp1,llm/)))
c      call WriteField_p('dvfi',reshape(dvfi,(/iip1,jjm,llm/)))
c      call WriteField_p('dtetafi',reshape(dtetafi,(/iip1,jmp1,llm/)))
cc$OMP END MASTER

c      ajout des tendances physiques:
c      ------------------------------
          CALL addfi_p( dtphys, leapf, forward   ,
     $                  ucov, vcov, teta , q   ,ps ,
     $                 dufi, dvfi, dtetafi , dqfi ,dpfi  )
          ! since addfi updates ps(), also update p(), masse() and pk()
          CALL pression_p(ip1jmp1,ap,bp,ps,p)
c$OMP BARRIER
          CALL massdair_p(p,masse)
c$OMP BARRIER
          if (pressure_exner) then
            CALL exner_hyb_p(ip1jmp1,ps,p,pks,pk,pkf)
          else
            CALL exner_milieu_p(ip1jmp1,ps,p,pks,pk,pkf)
          endif
c$OMP BARRIER

c      Couche superieure :
c      -------------------
         IF (iflag_top_bound > 0) THEN
           CALL top_bound_p(vcov,ucov,teta,masse,dtphys,
     $                       dutop)
        ijb=ij_begin
        ije=ij_end
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)        
        DO l=1,llm
          dutop(ijb:ije,l)=dutop(ijb:ije,l)/dtphys ! convert to a tendency in (m/s)/s
        ENDDO
c$OMP END DO NOWAIT        
         ENDIF ! of IF (ok_strato)
       
c$OMP BARRIER
c$OMP MASTER
        call VTe(VTphysiq)

        call VTb(VThallo)
c$OMP END MASTER

        call SetTag(Request_physic,800)
        call Register_SwapField(ucov,ucov,ip1jmp1,llm,
     *                               jj_Nb_caldyn,Request_physic)
        
        call Register_SwapField(vcov,vcov,ip1jm,llm,
     *                               jj_Nb_caldyn,Request_physic)
        
        call Register_SwapField(teta,teta,ip1jmp1,llm,
     *                               jj_Nb_caldyn,Request_physic)
        
        call Register_SwapField(masse,masse,ip1jmp1,llm,
     *                               jj_Nb_caldyn,Request_physic)

        call Register_SwapField(ps,ps,ip1jmp1,1,
     *                               jj_Nb_caldyn,Request_physic)

        call Register_SwapField(p,p,ip1jmp1,llmp1,
     *                               jj_Nb_caldyn,Request_physic)
        
        call Register_SwapField(pk,pk,ip1jmp1,llm,
     *                               jj_Nb_caldyn,Request_physic)
        
        call Register_SwapField(phis,phis,ip1jmp1,1,
     *                               jj_Nb_caldyn,Request_physic)
        
        call Register_SwapField(phi,phi,ip1jmp1,llm,
     *                               jj_Nb_caldyn,Request_physic)
        
        call Register_SwapField(w,w,ip1jmp1,llm,
     *                               jj_Nb_caldyn,Request_physic)

        do j=1,nqtot
        
          call Register_SwapField(q(1,1,j),q(1,1,j),ip1jmp1,llm,
     *                               jj_Nb_caldyn,Request_physic)
        
        enddo

        call SendRequest(Request_Physic)
c$OMP BARRIER
        call WaitRequest(Request_Physic)     

c$OMP BARRIER
c$OMP MASTER
       call VTe(VThallo)
       call SetDistrib(jj_Nb_caldyn)
c$OMP END MASTER
c$OMP BARRIER
c
c  Diagnostique de conservation de l'energie : difference
      IF ((ip_ebil_dyn.ge.1 ) .and. (nqtot > 1)) THEN 
          ztit='bil phys'
          CALL diagedyn(ztit,2,1,1,dtphys
     e  , ucov    , vcov , ps, p ,pk , teta , q(:,:,1), q(:,:,2))
      ENDIF 

cc$OMP MASTER      
c      if (debug) then
c       call WriteField_p('ucovfi',reshape(ucov,(/iip1,jmp1,llm/)))
c       call WriteField_p('vcovfi',reshape(vcov,(/iip1,jjm,llm/)))
c       call WriteField_p('tetafi',reshape(teta,(/iip1,jmp1,llm/)))
c      endif
cc$OMP END MASTER


c-jld
c$OMP MASTER
         call resume_timer(timer_caldyn)
         if (FirstPhysic) then
           ok_start_timer=.TRUE.
           FirstPhysic=.false.
         endif
c$OMP END MASTER
       ENDIF ! of IF( apphys )

      IF(iflag_phys.EQ.2) THEN ! "Newtonian" case
!   Academic case : Simple friction and Newtonan relaxation 
!   -------------------------------------------------------
c$OMP MASTER
         if (FirstPhysic) then
           ok_start_timer=.TRUE.
           FirstPhysic=.false.
         endif
c$OMP END MASTER

       ijb=ij_begin
       ije=ij_end
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK) 
       do l=1,llm
        teta(ijb:ije,l)=teta(ijb:ije,l)-dtvr*
     &         (teta(ijb:ije,l)-tetarappel(ijb:ije,l))*
     &                  (knewt_g+knewt_t(l)*clat4(ijb:ije))
       enddo ! of do l=1,llm
!$OMP END DO

       if (planet_type.eq."giant") then
          ! Intrinsic heat flux
          ! Aymeric -- for giant planets
          if (ihf .gt. 1.e-6) then
          !print *, '**** INTRINSIC HEAT FLUX ****', ihf 
          teta(ijb:ije,1) = teta(ijb:ije,1)
     &        + dtvr * aire(ijb:ije) * ihf / cpp / masse(ijb:ije,1)
          !print *, '**** d teta ' 
          !print *, dtvr * aire(ijb:ije) * ihf / cpp / masse(ijb:ije,1)
          endif
       endif

       call Register_Hallo(ucov,ip1jmp1,llm,0,1,1,0,Request_Physic)
       call Register_Hallo(vcov,ip1jm,llm,1,1,1,1,Request_Physic)
       call SendRequest(Request_Physic)
c$OMP BARRIER
       call WaitRequest(Request_Physic)     
c$OMP BARRIER
       call friction_p(ucov,vcov,dtvr)
!$OMP BARRIER

        ! Sponge layer (if any)
        IF (ok_strato) THEN
          CALL top_bound_p(vcov,ucov,teta,masse,dtvr,
     $                     dutop)
          ijb=ij_begin
          ije=ij_end
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)        
          DO l=1,llm
            dutop(ijb:ije,l)=dutop(ijb:ije,l)/dtvr ! convert to a tendency in (m/s)/s
          ENDDO
c$OMP END DO NOWAIT        
!$OMP BARRIER
        ENDIF ! of IF (ok_strato) 
      ENDIF ! of IF(iflag_phys.EQ.2)


        CALL pression_p ( ip1jmp1, ap, bp, ps, p                  )
c$OMP BARRIER
        if (pressure_exner) then
          CALL exner_hyb_p( ip1jmp1, ps, p, pks, pk, pkf )
        else
          CALL exner_milieu_p( ip1jmp1, ps, p, pks, pk, pkf )
        endif
c$OMP BARRIER
        CALL massdair_p(p,masse)
c$OMP BARRIER

cc$OMP END PARALLEL

c-----------------------------------------------------------------------
c   dissipation horizontale et verticale  des petites echelles:
c   ----------------------------------------------------------

      IF(apdiss) THEN

c$OMP MASTER
        call suspend_timer(timer_caldyn)
        
c       print*,'Entree dans la dissipation : Iteration No ',true_itau
c   calcul de l'energie cinetique avant dissipation
c       print *,'Passage dans la dissipation'

        call VTb(VThallo)
c$OMP END MASTER

        ! sponge layer
        if (callsponge) then
          call Register_Hallo(ps,ip1jm,1,1,1,1,1,Request_Dissip)
          call SendRequest(Request_Dissip)
c$OMP BARRIER
          call WaitRequest(Request_Dissip)
c$OMP BARRIER
c$OMP MASTER
          call VTe(VThallo)
          call VTb(VThallo)
c$OMP END MASTER
c$OMP BARRIER
          CALL sponge_p(ucov,vcov,teta,ps,dtdiss,mode_sponge)
        endif


c$OMP BARRIER

        call Register_SwapFieldHallo(ucov,ucov,ip1jmp1,llm,
     *                          jj_Nb_dissip,1,1,Request_dissip)

        call Register_SwapFieldHallo(vcov,vcov,ip1jm,llm,
     *                          jj_Nb_dissip,1,1,Request_dissip)

        call Register_SwapField(teta,teta,ip1jmp1,llm,
     *                          jj_Nb_dissip,Request_dissip)

        call Register_SwapField(p,p,ip1jmp1,llmp1,
     *                          jj_Nb_dissip,Request_dissip)

        call Register_SwapField(pk,pk,ip1jmp1,llm,
     *                          jj_Nb_dissip,Request_dissip)

        call SendRequest(Request_dissip)       
c$OMP BARRIER
        call WaitRequest(Request_dissip)       

c$OMP BARRIER
c$OMP MASTER
        call SetDistrib(jj_Nb_dissip)
        call VTe(VThallo)
        call VTb(VTdissipation)
        call start_timer(timer_dissip)
c$OMP END MASTER
c$OMP BARRIER

        call covcont_p(llm,ucov,vcov,ucont,vcont)
        call enercin_p(vcov,ucov,vcont,ucont,ecin0)

c   dissipation
! ADAPTATION GCM POUR CP(T)
        call tpot2t_glo_p(teta,temp,pk)

!        CALL FTRACE_REGION_BEGIN("dissip")
        CALL dissip_p(vcov,ucov,teta,p,dvdis,dudis,dtetadis)
!        CALL FTRACE_REGION_END("dissip")
         
        ijb=ij_begin
        ije=ij_end
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)        
        DO l=1,llm
          ucov(ijb:ije,l)=ucov(ijb:ije,l)+dudis(ijb:ije,l)
          dudis(ijb:ije,l)=dudis(ijb:ije,l)/dtdiss   ! passage en (m/s)/s
        ENDDO
c$OMP END DO NOWAIT        
        if (pole_sud) ije=ije-iip1
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)        
        DO l=1,llm
          vcov(ijb:ije,l)=vcov(ijb:ije,l)+dvdis(ijb:ije,l)
          dvdis(ijb:ije,l)=dvdis(ijb:ije,l)/dtdiss   ! passage en (m/s)/s
        ENDDO
c$OMP END DO NOWAIT        


c------------------------------------------------------------------------
        if (dissip_conservative) then
C       On rajoute la tendance due a la transform. Ec -> E therm. cree
C       lors de la dissipation
c$OMP BARRIER
c$OMP MASTER
            call suspend_timer(timer_dissip)
            call VTb(VThallo)
c$OMP END MASTER
            call Register_Hallo(ucov,ip1jmp1,llm,1,1,1,1,Request_Dissip)
            call Register_Hallo(vcov,ip1jm,llm,1,1,1,1,Request_Dissip)
            call SendRequest(Request_Dissip)
c$OMP BARRIER
            call WaitRequest(Request_Dissip)
c$OMP MASTER
            call VTe(VThallo)
            call resume_timer(timer_dissip)
c$OMP END MASTER
c$OMP BARRIER            
            call covcont_p(llm,ucov,vcov,ucont,vcont)
            call enercin_p(vcov,ucov,vcont,ucont,ecin)
            
            ijb=ij_begin
            ije=ij_end
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)            
            do l=1,llm
              do ij=ijb,ije
! ADAPTATION GCM POUR CP(T)
!                dtetaecdt(ij,l)= (ecin0(ij,l)-ecin(ij,l))/ pk(ij,l)
!                dtetadis(ij,l)=dtetadis(ij,l)+dtetaecdt(ij,l)
                temp(ij,l)=temp(ij,l) +
     &                      (ecin0(ij,l)-ecin(ij,l))/cpdet(temp(ij,l))
              enddo
            enddo
c$OMP END DO 
!        call t2tpot_p(ije-ijb+1,llm,temp(ijb:ije,:),ztetaec(ijb:ije,:),
!     &                            pk(ijb:ije,:))
         call t2tpot_glo_p(temp,ztetaec,pk)
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)            
            do l=1,llm
              do ij=ijb,ije
                dtetaecdt(ij,l)=ztetaec(ij,l)-teta(ij,l)
                dtetadis(ij,l)=dtetadis(ij,l)+dtetaecdt(ij,l)
              enddo
            enddo
c$OMP END DO NOWAIT
       endif ! of if (dissip_conservative)

       ijb=ij_begin
       ije=ij_end
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)            
         do l=1,llm
           do ij=ijb,ije
              teta(ij,l)=teta(ij,l)+dtetadis(ij,l)
              dtetadis(ij,l)=dtetadis(ij,l)/dtdiss   ! passage en K/s
           enddo
         enddo
c$OMP END DO NOWAIT         
c------------------------------------------------------------------------


c    .......        P. Le Van (  ajout  le 17/04/96  )   ...........
c   ...      Calcul de la valeur moyenne, unique de h aux poles  .....
c

        ijb=ij_begin
        ije=ij_end
         
        if (pole_nord) then
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l  =  1, llm
            DO ij =  1,iim
             tppn(ij)  = aire(  ij    ) * teta(  ij    ,l)
            ENDDO
             tpn  = SSUM(iim,tppn,1)/apoln

            DO ij = 1, iip1
             teta(  ij    ,l) = tpn
            ENDDO
          ENDDO
c$OMP END DO NOWAIT

         if (1 == 0) then
!!! Ehouarn: lines here 1) kill 1+1=2 in the dynamics
!!!                     2) should probably not be here anyway
!!! but are kept for those who would want to revert to previous behaviour
c$OMP MASTER               
          DO ij =  1,iim
            tppn(ij)  = aire(  ij    ) * ps (  ij    )
          ENDDO
            tpn  = SSUM(iim,tppn,1)/apoln
  
          DO ij = 1, iip1
            ps(  ij    ) = tpn
          ENDDO
c$OMP END MASTER
         endif ! of if (1 == 0)
        endif ! of of (pole_nord)
        
        if (pole_sud) then
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
          DO l  =  1, llm
            DO ij =  1,iim
             tpps(ij)  = aire(ij+ip1jm) * teta(ij+ip1jm,l)
            ENDDO
             tps  = SSUM(iim,tpps,1)/apols

            DO ij = 1, iip1
             teta(ij+ip1jm,l) = tps
            ENDDO
          ENDDO
c$OMP END DO NOWAIT

         if (1 == 0) then
!!! Ehouarn: lines here 1) kill 1+1=2 in the dynamics
!!!                     2) should probably not be here anyway
!!! but are kept for those who would want to revert to previous behaviour
c$OMP MASTER               
          DO ij =  1,iim
            tpps(ij)  = aire(ij+ip1jm) * ps (ij+ip1jm)
          ENDDO
            tps  = SSUM(iim,tpps,1)/apols
  
          DO ij = 1, iip1
            ps(ij+ip1jm) = tps
          ENDDO
c$OMP END MASTER
         endif ! of if (1 == 0)
        endif ! of if (pole_sud)


c$OMP BARRIER
c$OMP MASTER
        call VTe(VTdissipation)

        call stop_timer(timer_dissip)
        
        call VTb(VThallo)
c$OMP END MASTER
        call Register_SwapField(ucov,ucov,ip1jmp1,llm,
     *                          jj_Nb_caldyn,Request_dissip)

        call Register_SwapField(vcov,vcov,ip1jm,llm,
     *                          jj_Nb_caldyn,Request_dissip)

        call Register_SwapField(teta,teta,ip1jmp1,llm,
     *                          jj_Nb_caldyn,Request_dissip)

        call Register_SwapField(p,p,ip1jmp1,llmp1,
     *                          jj_Nb_caldyn,Request_dissip)

        call Register_SwapField(pk,pk,ip1jmp1,llm,
     *                          jj_Nb_caldyn,Request_dissip)

        call SendRequest(Request_dissip)       
c$OMP BARRIER
        call WaitRequest(Request_dissip)       

c$OMP BARRIER
c$OMP MASTER
        call SetDistrib(jj_Nb_caldyn)
        call VTe(VThallo)
        call resume_timer(timer_caldyn)
c        print *,'fin dissipation'
c$OMP END MASTER
c$OMP BARRIER
      END IF ! of IF(apdiss)

cc$OMP END PARALLEL

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
cym      call WriteField('ucov',reshape(ucov,(/iip1,jmp1,llm/)))
cym      call WriteField('vcov',reshape(vcov,(/iip1,jjm,llm/)))
c$OMP MASTER      
      call stop_timer(timer_caldyn)
c$OMP END MASTER
      IF (itau==itaumax) then
c$OMP MASTER
            call allgather_timer_average

      if (mpi_rank==0) then
        
        print *,'*********************************'
        print *,'******    TIMER CALDYN     ******'
        do i=0,mpi_size-1
          print *,'proc',i,' :   Nb Bandes  :',jj_nb_caldyn(i),
     &            '  : temps moyen :',
     &             timer_average(jj_nb_caldyn(i),timer_caldyn,i)
        enddo
      
        print *,'*********************************'
        print *,'******    TIMER VANLEER    ******'
        do i=0,mpi_size-1
          print *,'proc',i,' :   Nb Bandes  :',jj_nb_vanleer(i),
     &            '  : temps moyen :',
     &             timer_average(jj_nb_vanleer(i),timer_vanleer,i)
        enddo
      
        print *,'*********************************'
        print *,'******    TIMER DISSIP    ******'
        do i=0,mpi_size-1
          print *,'proc',i,' :   Nb Bandes  :',jj_nb_dissip(i),
     &            '  : temps moyen :',
     &             timer_average(jj_nb_dissip(i),timer_dissip,i)
        enddo
        
        print *,'*********************************'
        print *,'******    TIMER PHYSIC    ******'
        do i=0,mpi_size-1
          print *,'proc',i,' :   Nb Bandes  :',jj_nb_physic(i),
     &            '  : temps moyen :',
     &             timer_average(jj_nb_physic(i),timer_physic,i)
        enddo
        
      endif  
      
      print *,'Taille du Buffer MPI (REAL*8)',MaxBufferSize
      print *,'Taille du Buffer MPI utilise (REAL*8)',MaxBufferSize_Used
      print *, 'Temps total ecoule sur la parallelisation :',DiffTime()
      print *, 'Temps CPU ecoule sur la parallelisation :',DiffCpuTime()
      CALL print_filtre_timer
      call fin_getparam
        call finalize_parallel
c$OMP END MASTER
c$OMP BARRIER
        RETURN
      ENDIF ! of IF (itau==itaumax)
      
      IF ( .NOT.purmats ) THEN
c       ........................................................
c       ..............  schema matsuno + leapfrog  ..............
c       ........................................................

            IF(forward. OR. leapf) THEN
              itau= itau + 1
!              iday= day_ini+itau/day_step
!              time= REAL(itau-(iday-day_ini)*day_step)/day_step+time_0
!                IF(time.GT.1.) THEN
!                  time = time-1.
!                  iday = iday+1
!                ENDIF
            ENDIF


            IF( itau. EQ. itaufinp1 ) then

              if (flag_verif) then
                write(79,*) 'ucov',ucov
                write(80,*) 'vcov',vcov
                write(81,*) 'teta',teta
                write(82,*) 'ps',ps
                write(83,*) 'q',q
                if (nqtot > 2) then
                 WRITE(85,*) 'q1 = ',q(:,:,1)
                 WRITE(86,*) 'q3 = ',q(:,:,3)
                endif
              endif
  

c$OMP MASTER
              call fin_getparam
c$OMP END MASTER
#ifdef INCA
                 call finalize_inca
#endif
c$OMP MASTER
               call finalize_parallel
c$OMP END MASTER
              abort_message = 'Simulation finished'
              call abort_gcm(modname,abort_message,0)
              RETURN
            ENDIF
c-----------------------------------------------------------------------
c   ecriture du fichier histoire moyenne:
c   -------------------------------------

            IF(MOD(itau,iperiod).EQ.0 .OR. itau.EQ.itaufin) THEN
c$OMP BARRIER
               IF(itau.EQ.itaufin) THEN
                  iav=1
               ELSE
                  iav=0
               ENDIF
#ifdef CPP_IOIPSL
             IF (ok_dynzon) THEN 
              CALL bilan_dyn_p(2,dtvr*iperiod,dtvr*day_step*periodav, 
     ,             ps,masse,pk,pbaru,pbarv,teta,phi,ucov,vcov,q) 
c les traceurs ne sont pas sortis, trop lourd. 
c Peut changer eventuellement si besoin.
!                 CALL bilan_dyn_p(dtvr*iperiod,dtvr*day_step*periodav,
!     &                 ps,masse,pk,pbaru,pbarv,teta,phi,ucov,vcov,
!     &                 du,dudis,dutop,dufi)
              ENDIF !ok_dynzon
#endif
               IF (ok_dyn_ave) THEN
!$OMP MASTER
#ifdef CPP_IOIPSL
! Ehouarn: Gather fields and make master send to output
                call Gather_Field(vcov,ip1jm,llm,0)
                call Gather_Field(ucov,ip1jmp1,llm,0)
                call Gather_Field(teta,ip1jmp1,llm,0)
                call Gather_Field(pk,ip1jmp1,llm,0)
		call Gather_Field(phi,ip1jmp1,llm,0)
                 do iq=1,nqtot
                  call Gather_Field(q(1,1,iq),ip1jmp1,llm,0)
                 enddo
		call Gather_Field(masse,ip1jmp1,llm,0)
                call Gather_Field(ps,ip1jmp1,1,0)
		call Gather_Field(phis,ip1jmp1,1,0)
                if (mpi_rank==0) then
                 CALL writedynav(itau,vcov,
     &                 ucov,teta,pk,phi,q,masse,ps,phis)
		endif
#endif
!$OMP END MASTER
               ENDIF ! of IF (ok_dyn_ave)
            ENDIF ! of IF((MOD(itau,iperiod).EQ.0).OR.(itau.EQ.itaufin))

c-----------------------------------------------------------------------
c   ecriture de la bande histoire:
c   ------------------------------

            IF( MOD(itau,iecri).EQ.0) THEN
             ! Ehouarn: output only during LF or Backward Matsuno
	     if (leapf.or.(.not.leapf.and.(.not.forward))) then
c$OMP BARRIER

! ADAPTATION GCM POUR CP(T)
      call tpot2t_glo_p(teta,temp,pk)
      ijb=ij_begin
      ije=ij_end
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      do l=1,llm
        tsurpk(ijb:ije,l)=cpp*temp(ijb:ije,l)/pk(ijb:ije,l)
      enddo
!$OMP END DO

!$OMP MASTER
!              CALL geopot_p(ip1jmp1,teta,pk,pks,phis,phi)
      CALL geopot_p  ( ip1jmp1, tsurpk  , pk , pks,  phis  , phi   )
       
cym        unat=0.
        
              ijb=ij_begin
              ije=ij_end
        
              if (pole_nord) then
                ijb=ij_begin+iip1
                unat(1:iip1,:)=0.
              endif
        
              if (pole_sud) then 
                ije=ij_end-iip1
                unat(ij_end-iip1+1:ij_end,:)=0.
              endif
            
              do l=1,llm
                unat(ijb:ije,l)=ucov(ijb:ije,l)/cu(ijb:ije)
              enddo

              ijb=ij_begin
              ije=ij_end
              if (pole_sud) ije=ij_end-iip1
        
              do l=1,llm
                vnat(ijb:ije,l)=vcov(ijb:ije,l)/cv(ijb:ije)
              enddo
        
#ifdef CPP_IOIPSL
              if (ok_dyn_ins) then
! Ehouarn: Gather fields and make master write to output
                call Gather_Field(vcov,ip1jm,llm,0)
                call Gather_Field(ucov,ip1jmp1,llm,0)
                call Gather_Field(teta,ip1jmp1,llm,0)
		call Gather_Field(phi,ip1jmp1,llm,0)
                 do iq=1,nqtot
                  call Gather_Field(q(1,1,iq),ip1jmp1,llm,0)
                 enddo
		call Gather_Field(masse,ip1jmp1,llm,0)
                call Gather_Field(ps,ip1jmp1,1,0)
		call Gather_Field(phis,ip1jmp1,1,0)
                if (mpi_rank==0) then
	         CALL writehist(itau,vcov,ucov,teta,phi,q,masse,ps,phis)
		endif
!              CALL writehist_p(histid,histvid, itau,vcov, 
!     &                         ucov,teta,phi,q,masse,ps,phis)
! or use writefield_p
!      call WriteField_p('ucov',reshape(ucov,(/iip1,jmp1,llm/)))
!      call WriteField_p('vcov',reshape(vcov,(/iip1,jjm,llm/)))
!      call WriteField_p('teta',reshape(teta,(/iip1,jmp1,llm/)))
!      call WriteField_p('ps',reshape(ps,(/iip1,jmp1/)))
              endif ! of if (ok_dyn_ins)
#endif
! For some Grads outputs of fields
              if (output_grads_dyn) then
! Ehouarn: hope this works the way I think it does:
                  call Gather_Field(unat,ip1jmp1,llm,0)
                  call Gather_Field(vnat,ip1jm,llm,0)
                  call Gather_Field(teta,ip1jmp1,llm,0)
                  call Gather_Field(ps,ip1jmp1,1,0)
                  do iq=1,nqtot
                    call Gather_Field(q(1,1,iq),ip1jmp1,llm,0)
                  enddo
                  if (mpi_rank==0) then
#include "write_grads_dyn.h"
                  endif
              endif ! of if (output_grads_dyn)
!$OMP END MASTER
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

c           Write to restart.nc if required
            IF (lrestart) THEN
c$OMP BARRIER
c$OMP MASTER
              if (planet_type=="mars") then
                CALL dynredem1_p("restart.nc",REAL(itau)/REAL(day_step),
     &                           vcov,ucov,teta,q,masse,ps)
              else
                CALL dynredem1_p("restart.nc",start_time,
     &                           vcov,ucov,teta,q,masse,ps)
              endif
!              CLOSE(99)
c$OMP END MASTER
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
                   END IF
            ELSE

c      ......   pas leapfrog  .....

                 leapf = .TRUE.
                 dt  = 2.*dtvr
                 GO TO 2
            END IF ! of IF (MOD(itau,iperiod).EQ.0)
                   !    ELSEIF (MOD(itau-1,iperiod).EQ.0)


      ELSE ! of IF (.not.purmats)

c       ........................................................
c       ..............       schema  matsuno        ...............
c       ........................................................
            IF( forward )  THEN

             itau =  itau + 1
!             iday = day_ini+itau/day_step
!             time = REAL(itau-(iday-day_ini)*day_step)/day_step+time_0
!
!                  IF(time.GT.1.) THEN
!                   time = time-1.
!                   iday = iday+1
!                  ENDIF

               forward =  .FALSE.
               IF( itau. EQ. itaufinp1 ) then  
c$OMP MASTER
                 call fin_getparam
                 call finalize_parallel
c$OMP END MASTER
                 abort_message = 'Simulation finished'
                 call abort_gcm(modname,abort_message,0)
                 RETURN
               ENDIF
               GO TO 2

            ELSE ! of IF(forward) i.e. backward step

              IF(MOD(itau,iperiod).EQ.0 .OR. itau.EQ.itaufin) THEN
               IF(itau.EQ.itaufin) THEN
                  iav=1
               ELSE
                  iav=0
               ENDIF
#ifdef CPP_IOIPSL
               IF (ok_dynzon) THEN
               CALL bilan_dyn_p(2,dtvr*iperiod,dtvr*day_step*periodav,
     ,           ps,masse,pk,pbaru,pbarv,teta,phi,ucov,vcov,q)
c les traceurs ne sont pas sortis, trop lourd. 
c Peut changer eventuellement si besoin.
!                 CALL bilan_dyn_p(dtvr*iperiod,dtvr*day_step*periodav,
!     &                 ps,masse,pk,pbaru,pbarv,teta,phi,ucov,vcov,
!     &                 du,dudis,dutop,dufi)
               END IF !ok_dynzon
#endif
               IF (ok_dyn_ave) THEN
!$OMP MASTER
#ifdef CPP_IOIPSL
! Ehouarn: Gather fields and make master send to output
                call Gather_Field(vcov,ip1jm,llm,0)
                call Gather_Field(ucov,ip1jmp1,llm,0)
                call Gather_Field(teta,ip1jmp1,llm,0)
                call Gather_Field(pk,ip1jmp1,llm,0)
		call Gather_Field(phi,ip1jmp1,llm,0)
                do iq=1,nqtot
                  call Gather_Field(q(1,1,iq),ip1jmp1,llm,0)
                enddo
		call Gather_Field(masse,ip1jmp1,llm,0)
                call Gather_Field(ps,ip1jmp1,1,0)
		call Gather_Field(phis,ip1jmp1,1,0)
                if (mpi_rank==0) then
                 CALL writedynav(itau,vcov,
     &                 ucov,teta,pk,phi,q,masse,ps,phis)
		endif
#endif
!$OMP END MASTER
               ENDIF ! of IF (ok_dyn_ave)

              ENDIF ! of IF(MOD(itau,iperiod).EQ.0 .OR. itau.EQ.itaufin)


               IF(MOD(itau,iecri         ).EQ.0) THEN
c              IF(MOD(itau,iecri*day_step).EQ.0) THEN
c$OMP BARRIER

! ADAPTATION GCM POUR CP(T)
                call tpot2t_glo_p(teta,temp,pk)
                ijb=ij_begin
                ije=ij_end
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
                do l=1,llm      
                  tsurpk(ijb:ije,l)=cpp*temp(ijb:ije,l)/
     &                                             pk(ijb:ije,l)
                enddo
!$OMP END DO

!$OMP MASTER
!                CALL geopot_p(ip1jmp1,teta,pk,pks,phis,phi)
                CALL geopot_p(ip1jmp1,tsurpk,pk,pks,phis,phi)

cym        unat=0.
                ijb=ij_begin
                ije=ij_end
        
                if (pole_nord) then
                  ijb=ij_begin+iip1
                  unat(1:iip1,:)=0.
                endif
        
                if (pole_sud) then 
                  ije=ij_end-iip1
                  unat(ij_end-iip1+1:ij_end,:)=0.
                endif
            
                do l=1,llm
                  unat(ijb:ije,l)=ucov(ijb:ije,l)/cu(ijb:ije)
                enddo

                ijb=ij_begin
                ije=ij_end
                if (pole_sud) ije=ij_end-iip1
        
                do l=1,llm
                  vnat(ijb:ije,l)=vcov(ijb:ije,l)/cv(ijb:ije)
                enddo

#ifdef CPP_IOIPSL
              if (ok_dyn_ins) then
! Ehouarn: Gather fields and make master send to output
                call Gather_Field(vcov,ip1jm,llm,0)
                call Gather_Field(ucov,ip1jmp1,llm,0)
                call Gather_Field(teta,ip1jmp1,llm,0)
		call Gather_Field(phi,ip1jmp1,llm,0)
                 do iq=1,nqtot
                  call Gather_Field(q(1,1,iq),ip1jmp1,llm,0)
                 enddo
		call Gather_Field(masse,ip1jmp1,llm,0)
                call Gather_Field(ps,ip1jmp1,1,0)
		call Gather_Field(phis,ip1jmp1,1,0)
                if (mpi_rank==0) then
                 CALL writehist(itau,vcov,ucov,teta,phi,q,masse,ps,phis)
		endif
!                CALL writehist_p(histid, histvid, itau,vcov , 
!     &                           ucov,teta,phi,q,masse,ps,phis)
              endif ! of if (ok_dyn_ins)
#endif
! For some Grads output (but does it work?)
                if (output_grads_dyn) then
                  call Gather_Field(unat,ip1jmp1,llm,0)
                  call Gather_Field(vnat,ip1jm,llm,0)
                  call Gather_Field(teta,ip1jmp1,llm,0)
                  call Gather_Field(ps,ip1jmp1,1,0)
                   do iq=1,nqtot
                    call Gather_Field(q(1,1,iq),ip1jmp1,llm,0)
                   enddo
c      
                  if (mpi_rank==0) then
#include "write_grads_dyn.h"
                  endif
                endif ! of if (output_grads_dyn)

!$OMP END MASTER
              ENDIF ! of IF(MOD(itau,iecri).EQ.0)

c             Determine whether to write to the restart.nc file
c             Decision can't be made in one IF statement as if
c             ecritstart==0 there will be a divide-by-zero error
              lrestart = .false.
              IF (itau.EQ.itaufin) THEN
                lrestart = .true.
              ELSE IF (ecritstart.GT.0) THEN
                IF (MOD(itau,ecritstart).EQ.0) lrestart  = .true.
              ENDIF

c             Write to restart.nc if required
              IF (lrestart) THEN
c$OMP MASTER
                if (planet_type=="mars") then
                  CALL dynredem1_p("restart.nc",
     &                              REAL(itau)/REAL(day_step),
     &                               vcov,ucov,teta,q,masse,ps)
                else
                  CALL dynredem1_p("restart.nc",start_time,
     &                               vcov,ucov,teta,q,masse,ps)
                endif
c$OMP END MASTER
              ENDIF ! of IF (lrestart)

              forward = .TRUE.
              GO TO  1

            ENDIF ! of IF (forward)

      END IF ! of IF(.not.purmats)
c$OMP MASTER
      call fin_getparam
      call finalize_parallel
c$OMP END MASTER
      RETURN
      END
