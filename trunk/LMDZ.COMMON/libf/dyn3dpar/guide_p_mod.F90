!
! $Id$
!
MODULE guide_p_mod

!=======================================================================
!   Auteur:  F.Hourdin
!            F. Codron 01/09
!=======================================================================

  USE getparam
  USE Write_Field_p
  USE netcdf, ONLY: nf90_nowrite, nf90_open, nf90_inq_varid, nf90_close
  USE pres2lev_mod

  IMPLICIT NONE

! ---------------------------------------------
! Declarations des cles logiques et parametres 
! ---------------------------------------------
  INTEGER, PRIVATE, SAVE  :: iguide_read,iguide_int,iguide_sav
  INTEGER, PRIVATE, SAVE  :: nlevnc, guide_plevs
  LOGICAL, PRIVATE, SAVE  :: guide_u,guide_v,guide_T,guide_Q,guide_P
  LOGICAL, PRIVATE, SAVE  :: guide_hr,guide_teta  
  LOGICAL, PRIVATE, SAVE  :: guide_BL,guide_reg,guide_add,gamma4,guide_zon 
  LOGICAL, PRIVATE, SAVE  :: invert_p,invert_y,ini_anal
  LOGICAL, PRIVATE, SAVE  :: guide_2D,guide_sav,guide_modele
  
  REAL, PRIVATE, SAVE     :: tau_min_u,tau_max_u
  REAL, PRIVATE, SAVE     :: tau_min_v,tau_max_v
  REAL, PRIVATE, SAVE     :: tau_min_T,tau_max_T
  REAL, PRIVATE, SAVE     :: tau_min_Q,tau_max_Q
  REAL, PRIVATE, SAVE     :: tau_min_P,tau_max_P

  REAL, PRIVATE, SAVE     :: lat_min_g,lat_max_g
  REAL, PRIVATE, SAVE     :: lon_min_g,lon_max_g
  REAL, PRIVATE, SAVE     :: tau_lon,tau_lat

  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE, SAVE     :: alpha_u,alpha_v 
  REAL, ALLOCATABLE, DIMENSION(:, :), PRIVATE, SAVE     :: alpha_T,alpha_Q 
  REAL, ALLOCATABLE, DIMENSION(:), PRIVATE, SAVE     :: alpha_P,alpha_pcor
  
! ---------------------------------------------
! Variables de guidage
! ---------------------------------------------
! Variables des fichiers de guidage
  REAL, ALLOCATABLE, DIMENSION(:,:,:), PRIVATE, SAVE   :: unat1,unat2
  REAL, ALLOCATABLE, DIMENSION(:,:,:), PRIVATE, SAVE   :: vnat1,vnat2
  REAL, ALLOCATABLE, DIMENSION(:,:,:), PRIVATE, SAVE   :: tnat1,tnat2
  REAL, ALLOCATABLE, DIMENSION(:,:,:), PRIVATE, SAVE   :: qnat1,qnat2
  REAL, ALLOCATABLE, DIMENSION(:,:,:), PRIVATE, SAVE   :: pnat1,pnat2
  REAL, ALLOCATABLE, DIMENSION(:,:),   PRIVATE, SAVE   :: psnat1,psnat2
  REAL, ALLOCATABLE, DIMENSION(:),     PRIVATE, SAVE   :: apnc,bpnc
! Variables aux dimensions du modele
  REAL, ALLOCATABLE, DIMENSION(:,:),   PRIVATE, SAVE   :: ugui1,ugui2
  REAL, ALLOCATABLE, DIMENSION(:,:),   PRIVATE, SAVE   :: vgui1,vgui2
  REAL, ALLOCATABLE, DIMENSION(:,:),   PRIVATE, SAVE   :: tgui1,tgui2
  REAL, ALLOCATABLE, DIMENSION(:,:),   PRIVATE, SAVE   :: qgui1,qgui2
  REAL, ALLOCATABLE, DIMENSION(:),   PRIVATE, SAVE   :: psgui1,psgui2
  
  INTEGER,SAVE,PRIVATE :: ijb_u,ijb_v,ije_u,ije_v,ijn_u,ijn_v
  INTEGER,SAVE,PRIVATE :: jjb_u,jjb_v,jje_u,jje_v,jjn_u,jjn_v


CONTAINS
!=======================================================================

  SUBROUTINE guide_init

    USE control_mod, ONLY: day_step
    USE serre_mod, ONLY: grossismx

    IMPLICIT NONE
  
    INCLUDE "dimensions.h"
    INCLUDE "paramet.h"
    INCLUDE "netcdf.inc"

    INTEGER                :: error,ncidpl,rid,rcod
    CHARACTER (len = 80)   :: abort_message
    CHARACTER (len = 20)   :: modname = 'guide_init'

! ---------------------------------------------
! Lecture des parametres:  
! ---------------------------------------------
! Variables guidees
    CALL getpar('guide_u',.true.,guide_u,'guidage de u')
    CALL getpar('guide_v',.true.,guide_v,'guidage de v')
    CALL getpar('guide_T',.true.,guide_T,'guidage de T')
    CALL getpar('guide_P',.true.,guide_P,'guidage de P')
    CALL getpar('guide_Q',.true.,guide_Q,'guidage de Q')
    CALL getpar('guide_hr',.true.,guide_hr,'guidage de Q par H.R')
    CALL getpar('guide_teta',.false.,guide_teta,'guidage de T par Teta')

    CALL getpar('guide_add',.false.,guide_add,'for�age constant?')
    CALL getpar('guide_zon',.false.,guide_zon,'guidage moy zonale')
    if (guide_zon .and. abs(grossismx - 1.) > 0.01) &
         call abort_gcm("guide_init", &
         "zonal nudging requires grid regular in longitude", 1)

!   Constantes de rappel. Unite : fraction de jour
    CALL getpar('tau_min_u',0.02,tau_min_u,'Cste de rappel min, u')
    CALL getpar('tau_max_u', 10.,tau_max_u,'Cste de rappel max, u')
    CALL getpar('tau_min_v',0.02,tau_min_v,'Cste de rappel min, v')
    CALL getpar('tau_max_v', 10.,tau_max_v,'Cste de rappel max, v')
    CALL getpar('tau_min_T',0.02,tau_min_T,'Cste de rappel min, T')
    CALL getpar('tau_max_T', 10.,tau_max_T,'Cste de rappel max, T')
    CALL getpar('tau_min_Q',0.02,tau_min_Q,'Cste de rappel min, Q')
    CALL getpar('tau_max_Q', 10.,tau_max_Q,'Cste de rappel max, Q')
    CALL getpar('tau_min_P',0.02,tau_min_P,'Cste de rappel min, P')
    CALL getpar('tau_max_P', 10.,tau_max_P,'Cste de rappel max, P')
    CALL getpar('gamma4',.false.,gamma4,'Zone sans rappel elargie')
    CALL getpar('guide_BL',.true.,guide_BL,'guidage dans C.Lim')
    
! Sauvegarde du for�age
    CALL getpar('guide_sav',.false.,guide_sav,'sauvegarde guidage')
    CALL getpar('iguide_sav',4,iguide_sav,'freq. sauvegarde guidage')
    ! frequences f>0: fx/jour; f<0: tous les f jours; f=0: 1 seule fois.
    IF (iguide_sav.GT.0) THEN
       iguide_sav=day_step/iguide_sav
    ELSE if (iguide_sav == 0) then
       iguide_sav = huge(0)
    ELSE
       iguide_sav=day_step*iguide_sav
    ENDIF

! Guidage regional seulement (sinon constant ou suivant le zoom)
    CALL getpar('guide_reg',.false.,guide_reg,'guidage regional')
    CALL getpar('lat_min_g',-90.,lat_min_g,'Latitude mini guidage ')
    CALL getpar('lat_max_g', 90.,lat_max_g,'Latitude maxi guidage ')
    CALL getpar('lon_min_g',-180.,lon_min_g,'longitude mini guidage ')
    CALL getpar('lon_max_g', 180.,lon_max_g,'longitude maxi guidage ')
    CALL getpar('tau_lat', 5.,tau_lat,'raideur lat guide regional ')
    CALL getpar('tau_lon', 5.,tau_lon,'raideur lon guide regional ')

! Parametres pour lecture des fichiers
    CALL getpar('iguide_read',4,iguide_read,'freq. lecture guidage')
    CALL getpar('iguide_int',4,iguide_int,'freq. interpolation vert')
    IF (iguide_int.EQ.0) THEN
        iguide_int=1
    ELSEIF (iguide_int.GT.0) THEN
        iguide_int=day_step/iguide_int
    ELSE
        iguide_int=day_step*iguide_int
    ENDIF
    CALL getpar('guide_plevs',0,guide_plevs,'niveaux pression fichiers guidage')
    ! Pour compatibilite avec ancienne version avec guide_modele
    CALL getpar('guide_modele',.false.,guide_modele,'niveaux pression ap+bp*psol')
    IF (guide_modele) THEN
        guide_plevs=1
    ENDIF
    ! Fin raccord
    CALL getpar('ini_anal',.false.,ini_anal,'Etat initial = analyse')
    CALL getpar('guide_invertp',.true.,invert_p,'niveaux p inverses')
    CALL getpar('guide_inverty',.true.,invert_y,'inversion N-S')
    CALL getpar('guide_2D',.false.,guide_2D,'fichier guidage lat-P')

! ---------------------------------------------
! Determination du nombre de niveaux verticaux
! des fichiers guidage
! ---------------------------------------------
    ncidpl=-99
    if (guide_plevs.EQ.1) then
       if (ncidpl.eq.-99) then
          rcod=nf90_open('apbp.nc',Nf90_NOWRITe, ncidpl)
          if (rcod.NE.NF_NOERR) THEN
             print *,'Guide: probleme -> pas de fichier apbp.nc'
             CALL abort_gcm(modname,abort_message,1)
          endif
       endif
    elseif (guide_plevs.EQ.2) then
       if (ncidpl.EQ.-99) then
          rcod=nf90_open('P.nc',Nf90_NOWRITe,ncidpl)
          if (rcod.NE.NF_NOERR) THEN
             print *,'Guide: probleme -> pas de fichier P.nc'
             CALL abort_gcm(modname,abort_message,1)
          endif
       endif
    elseif (guide_u) then
       if (ncidpl.eq.-99) then
          rcod=nf90_open('u.nc',Nf90_NOWRITe,ncidpl)
          if (rcod.NE.NF_NOERR) THEN
             print *,'Guide: probleme -> pas de fichier u.nc'
             CALL abort_gcm(modname,abort_message,1)
          endif
       endif
    elseif (guide_v) then
       if (ncidpl.eq.-99) then
          rcod=nf90_open('v.nc',nf90_nowrite,ncidpl)
          if (rcod.NE.NF_NOERR) THEN
             print *,'Guide: probleme -> pas de fichier v.nc'
             CALL abort_gcm(modname,abort_message,1)
          endif
       endif
    elseif (guide_T) then
       if (ncidpl.eq.-99) then
          rcod=nf90_open('T.nc',nf90_nowrite,ncidpl)
          if (rcod.NE.NF_NOERR) THEN
             print *,'Guide: probleme -> pas de fichier T.nc'
             CALL abort_gcm(modname,abort_message,1)
          endif
       endif
    elseif (guide_Q) then
       if (ncidpl.eq.-99) then
          rcod=nf90_open('hur.nc',nf90_nowrite, ncidpl)
          if (rcod.NE.NF_NOERR) THEN
             print *,'Guide: probleme -> pas de fichier hur.nc'
             CALL abort_gcm(modname,abort_message,1)
          endif
       endif
    endif 
    error=NF_INQ_DIMID(ncidpl,'LEVEL',rid)
    IF (error.NE.NF_NOERR) error=NF_INQ_DIMID(ncidpl,'PRESSURE',rid)
    IF (error.NE.NF_NOERR) THEN
        print *,'Guide: probleme lecture niveaux pression'
        CALL abort_gcm(modname,abort_message,1)
    ENDIF
    error=NF_INQ_DIMLEN(ncidpl,rid,nlevnc)
    print *,'Guide: nombre niveaux vert. nlevnc', nlevnc 
    rcod = nf90_close(ncidpl)

! ---------------------------------------------
! Allocation des variables
! ---------------------------------------------
    abort_message='pb in allocation guide'

    ALLOCATE(apnc(nlevnc), stat = error)
    IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
    ALLOCATE(bpnc(nlevnc), stat = error)
    IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
    apnc=0.;bpnc=0.

    ALLOCATE(alpha_pcor(llm), stat = error)
    IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
    ALLOCATE(alpha_u(ip1jmp1), stat = error)
    IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
    ALLOCATE(alpha_v(ip1jm), stat = error)
    IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
    ALLOCATE(alpha_T(iip1, jjp1), stat = error)
    IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
    ALLOCATE(alpha_Q(iip1, jjp1), stat = error)
    IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
    ALLOCATE(alpha_P(ip1jmp1), stat = error)
    IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
    alpha_u=0.;alpha_v=0;alpha_T=0;alpha_Q=0;alpha_P=0
    
    IF (guide_u) THEN
        ALLOCATE(unat1(iip1,jjp1,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(ugui1(ip1jmp1,llm), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(unat2(iip1,jjp1,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(ugui2(ip1jmp1,llm), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        unat1=0.;unat2=0.;ugui1=0.;ugui2=0.
    ENDIF

    IF (guide_T) THEN
        ALLOCATE(tnat1(iip1,jjp1,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(tgui1(ip1jmp1,llm), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(tnat2(iip1,jjp1,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(tgui2(ip1jmp1,llm), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        tnat1=0.;tnat2=0.;tgui1=0.;tgui2=0.
    ENDIF
     
    IF (guide_Q) THEN
        ALLOCATE(qnat1(iip1,jjp1,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(qgui1(ip1jmp1,llm), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(qnat2(iip1,jjp1,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(qgui2(ip1jmp1,llm), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        qnat1=0.;qnat2=0.;qgui1=0.;qgui2=0.
    ENDIF

    IF (guide_v) THEN
        ALLOCATE(vnat1(iip1,jjm,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(vgui1(ip1jm,llm), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(vnat2(iip1,jjm,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(vgui2(ip1jm,llm), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        vnat1=0.;vnat2=0.;vgui1=0.;vgui2=0.
    ENDIF

    IF (guide_plevs.EQ.2) THEN
        ALLOCATE(pnat1(iip1,jjp1,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(pnat2(iip1,jjp1,nlevnc), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        pnat1=0.;pnat2=0.;
    ENDIF

    IF (guide_P.OR.guide_plevs.EQ.1) THEN
        ALLOCATE(psnat1(iip1,jjp1), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(psnat2(iip1,jjp1), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        psnat1=0.;psnat2=0.;
    ENDIF
    IF (guide_P) THEN
        ALLOCATE(psgui2(ip1jmp1), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        ALLOCATE(psgui1(ip1jmp1), stat = error)
        IF (error /= 0) CALL abort_gcm(modname,abort_message,1)
        psgui1=0.;psgui2=0.
    ENDIF

! ---------------------------------------------
!   Lecture du premier etat de guidage.
! ---------------------------------------------
    IF (guide_2D) THEN
        CALL guide_read2D(1)
    ELSE
        CALL guide_read(1)
    ENDIF
    IF (guide_v) vnat1=vnat2
    IF (guide_u) unat1=unat2
    IF (guide_T) tnat1=tnat2
    IF (guide_Q) qnat1=qnat2
    IF (guide_plevs.EQ.2) pnat1=pnat2
    IF (guide_P.OR.guide_plevs.EQ.1) psnat1=psnat2

  END SUBROUTINE guide_init

!=======================================================================
  SUBROUTINE guide_main(itau,ucov,vcov,teta,q,masse,ps)
    use exner_hyb_p_m, only: exner_hyb_p
    use exner_milieu_p_m, only: exner_milieu_p
    USE parallel_lmdz
    USE control_mod
    USE comconst_mod, ONLY: daysec, dtvr, cpp, kappa
    USE comvert_mod, ONLY: ap, bp, preff, presnivs, pressure_exner
    
    IMPLICIT NONE
  
    INCLUDE "dimensions.h"
    INCLUDE "paramet.h"

    ! Variables entree
    INTEGER,                       INTENT(IN)    :: itau !pas de temps
    REAL, DIMENSION (ip1jmp1,llm), INTENT(INOUT) :: ucov,teta,q,masse
    REAL, DIMENSION (ip1jm,llm),   INTENT(INOUT) :: vcov
    REAL, DIMENSION (ip1jmp1),     INTENT(INOUT) :: ps

    ! Variables locales
    LOGICAL, SAVE :: first=.TRUE.
    LOGICAL       :: f_out ! sortie guidage
    REAL, DIMENSION (ip1jmp1,llm) :: f_add ! var aux: champ de guidage
    ! Variables pour fonction Exner (P milieu couche)
    REAL, DIMENSION (iip1,jjp1,llm)    :: pk
    REAL, DIMENSION (iip1,jjp1)        :: pks    
    REAL                               :: unskap
    REAL, DIMENSION (ip1jmp1,llmp1)    :: p ! besoin si guide_P
    ! Compteurs temps:
    INTEGER, SAVE :: step_rea,count_no_rea,itau_test ! lecture guidage
    REAL          :: ditau, dday_step
    REAL          :: tau,reste ! position entre 2 etats de guidage
    REAL, SAVE    :: factt ! pas de temps en fraction de jour
    
    INTEGER       :: i,j,l
    
    ijb_u=ij_begin ; ije_u=ij_end ; ijn_u=ije_u-ijb_u+1  
    jjb_u=jj_begin ; jje_u=jj_end ; jjn_u=jje_u-jjb_u+1 
    ijb_v=ij_begin ; ije_v=ij_end ; ijn_v=ije_v-ijb_v+1   
    jjb_v=jj_begin ; jje_v=jj_end ; jjn_v=jje_v-jjb_v+1 
    IF (pole_sud) THEN
      ije_v=ij_end-iip1
      jje_v=jj_end-1
      ijn_v=ije_v-ijb_v+1
      jjn_v=jje_v-jjb_v+1 
    ENDIF
      
!    CALL AllGather_Field(ucov,ip1jmp1,llm)
!    CALL AllGather_Field(vcov,ip1jm,llm)
!    CALL AllGather_Field(teta,ip1jmp1,llm)
!    CALL AllGather_Field(ps,ip1jmp1,1)
!    CALL AllGather_Field(q,ip1jmp1,llm)
    
!-----------------------------------------------------------------------
! Initialisations au premier passage
!-----------------------------------------------------------------------

    IF (first) THEN
        first=.FALSE.
        CALL guide_init 
        itau_test=1001
        step_rea=1
        count_no_rea=0
! Calcul des constantes de rappel
        factt=dtvr*iperiod/daysec 
        call tau2alpha(3,iip1,jjm ,factt,tau_min_v,tau_max_v,alpha_v)
        call tau2alpha(2,iip1,jjp1,factt,tau_min_u,tau_max_u,alpha_u)
        call tau2alpha(1,iip1,jjp1,factt,tau_min_T,tau_max_T,alpha_T)
        call tau2alpha(1,iip1,jjp1,factt,tau_min_P,tau_max_P,alpha_P)
        call tau2alpha(1,iip1,jjp1,factt,tau_min_Q,tau_max_Q,alpha_Q)
! correction de rappel dans couche limite
        if (guide_BL) then
             alpha_pcor(:)=1.
        else
            do l=1,llm
                alpha_pcor(l)=(1.+tanh((0.85-presnivs(l)/preff)/0.05))/2.
            enddo
        endif
! ini_anal: etat initial egal au guidage        
        IF (ini_anal) THEN
            CALL guide_interp(ps,teta)
            IF (guide_u) ucov(ijb_u:ije_u,:)=ugui2(ijb_u:ije_u,:)
            IF (guide_v) vcov(ijb_v:ije_v,:)=ugui2(ijb_v:ije_v,:)
            IF (guide_T) teta(ijb_u:ije_u,:)=tgui2(ijb_u:ije_u,:)
            IF (guide_Q) q(ijb_u:ije_u,:)=qgui2(ijb_u:ije_u,:)
            IF (guide_P) THEN
                ps(ijb_u:ije_u)=psgui2(ijb_u:ije_u)
                CALL pression_p(ip1jmp1,ap,bp,ps,p)
                CALL massdair_p(p,masse)
            ENDIF
            RETURN
        ENDIF
! Verification structure guidage
        IF (guide_u) THEN
            CALL writefield_p('unat',unat1)
            CALL writefield_p('ucov',RESHAPE(ucov,(/iip1,jjp1,llm/)))
        ENDIF
        IF (guide_T) THEN
            CALL writefield_p('tnat',tnat1)
            CALL writefield_p('teta',RESHAPE(teta,(/iip1,jjp1,llm/)))
        ENDIF

    ENDIF !first

!-----------------------------------------------------------------------
! Lecture des fichiers de guidage ?
!-----------------------------------------------------------------------
    IF (iguide_read.NE.0) THEN
      ditau=real(itau)
      dday_step=real(day_step)
      IF (iguide_read.LT.0) THEN
          tau=ditau/dday_step/REAL(iguide_read)
      ELSE
          tau=REAL(iguide_read)*ditau/dday_step
      ENDIF
      reste=tau-AINT(tau)
      IF (reste.EQ.0.) THEN
          IF (itau_test.EQ.itau) THEN
              write(*,*)'deuxieme passage de advreel a itau=',itau
              stop
          ELSE
              IF (guide_v) vnat1(:,jjb_v:jje_v,:)=vnat2(:,jjb_v:jje_v,:)
              IF (guide_u) unat1(:,jjb_u:jje_u,:)=unat2(:,jjb_u:jje_u,:)
              IF (guide_T) tnat1(:,jjb_u:jje_u,:)=tnat2(:,jjb_u:jje_u,:)
              IF (guide_Q) qnat1(:,jjb_u:jje_u,:)=qnat2(:,jjb_u:jje_u,:)
              IF (guide_plevs.EQ.2) pnat1(:,jjb_u:jje_u,:)=pnat2(:,jjb_u:jje_u,:)
              IF (guide_P.OR.guide_plevs.EQ.1) psnat1(:,jjb_u:jje_u)=psnat2(:,jjb_u:jje_u)
              step_rea=step_rea+1
              itau_test=itau
              print*,'Lecture fichiers guidage, pas ',step_rea, &
                    'apres ',count_no_rea,' non lectures'
              IF (guide_2D) THEN
                  CALL guide_read2D(step_rea)
              ELSE
                  CALL guide_read(step_rea)
              ENDIF
              count_no_rea=0
          ENDIF
      ELSE
        count_no_rea=count_no_rea+1

      ENDIF
    ENDIF !iguide_read=0

!-----------------------------------------------------------------------
! Interpolation et conversion des champs de guidage
!-----------------------------------------------------------------------
    IF (MOD(itau,iguide_int).EQ.0) THEN
        CALL guide_interp(ps,teta)
    ENDIF
! Repartition entre 2 etats de guidage
    IF (iguide_read.NE.0) THEN
        tau=reste
    ELSE
        tau=1.
    ENDIF

!-----------------------------------------------------------------------
!   Ajout des champs de guidage 
!-----------------------------------------------------------------------
! Sauvegarde du guidage?
    f_out=((MOD(itau,iguide_sav).EQ.0).AND.guide_sav)  
    IF (f_out) THEN
!       Calcul niveaux pression milieu de couches 
        CALL pression_p( ip1jmp1, ap, bp, ps, p )
        if (pressure_exner) then
          CALL exner_hyb_p(ip1jmp1,ps,p,pks,pk)
        else
          CALL exner_milieu_p(ip1jmp1,ps,p,pks,pk)
        endif
        unskap=1./kappa
        DO l = 1, llm
           DO j=jjb_u,jje_u
              DO i =1, iip1
                 p(i+(j-1)*iip1,l) = preff * ( pk(i,j,l)/cpp) ** unskap
              ENDDO
           ENDDO
        ENDDO
        CALL guide_out("SP",jjp1,llm,p(1:ip1jmp1,1:llm),1.)
    ENDIF
    
    if (guide_u) then
        if (guide_add) then
           f_add(ijb_u:ije_u,:)=(1.-tau)*ugui1(ijb_u:ije_u,:)+tau*ugui2(ijb_u:ije_u,:)
        else
           f_add(ijb_u:ije_u,:)=(1.-tau)*ugui1(ijb_u:ije_u,:)+tau*ugui2(ijb_u:ije_u,:)-ucov(ijb_u:ije_u,:)
        endif 

        if (guide_zon) CALL guide_zonave(1,jjp1,llm,f_add)
        CALL guide_addfield(ip1jmp1,llm,f_add,alpha_u)
        IF (f_out) CALL guide_out("u",jjp1,llm,ucov,factt)
        IF (f_out) CALL guide_out("ua",jjp1,llm,(1.-tau)*ugui1(:,:)+tau*ugui2(:,:),factt)
        IF (f_out) CALL guide_out("ucov",jjp1,llm,f_add(:,:)/factt,factt)
        ucov(ijb_u:ije_u,:)=ucov(ijb_u:ije_u,:)+f_add(ijb_u:ije_u,:)
    endif

    if (guide_T) then
        if (guide_add) then
           f_add(ijb_u:ije_u,:)=(1.-tau)*tgui1(ijb_u:ije_u,:)+tau*tgui2(ijb_u:ije_u,:)
        else
           f_add(ijb_u:ije_u,:)=(1.-tau)*tgui1(ijb_u:ije_u,:)+tau*tgui2(ijb_u:ije_u,:)-teta(ijb_u:ije_u,:)
        endif 
        if (guide_zon) CALL guide_zonave(2,jjp1,llm,f_add)
        CALL guide_addfield(ip1jmp1,llm,f_add,alpha_T)
        IF (f_out) CALL guide_out("teta",jjp1,llm,f_add(:,:)/factt,factt)
        teta(ijb_u:ije_u,:)=teta(ijb_u:ije_u,:)+f_add(ijb_u:ije_u,:)
    endif

    if (guide_P) then
        if (guide_add) then
           f_add(ijb_u:ije_u,1)=(1.-tau)*psgui1(ijb_u:ije_u)+tau*psgui2(ijb_u:ije_u)
        else
           f_add(ijb_u:ije_u,1)=(1.-tau)*psgui1(ijb_u:ije_u)+tau*psgui2(ijb_u:ije_u)-ps(ijb_u:ije_u)
        endif 
        if (guide_zon) CALL guide_zonave(2,jjp1,1,f_add(1:ip1jmp1,1))
        CALL guide_addfield(ip1jmp1,1,f_add(1:ip1jmp1,1),alpha_P)
        IF (f_out) CALL guide_out("ps",jjp1,1,f_add(1:ip1jmp1,1)/factt,factt)
        ps(ijb_u:ije_u)=ps(ijb_u:ije_u)+f_add(ijb_u:ije_u,1)
        CALL pression_p(ip1jmp1,ap,bp,ps,p)
        CALL massdair_p(p,masse)
    endif

    if (guide_Q) then
        if (guide_add) then
           f_add(ijb_u:ije_u,:)=(1.-tau)*qgui1(ijb_u:ije_u,:)+tau*qgui2(ijb_u:ije_u,:)
        else
           f_add(ijb_u:ije_u,:)=(1.-tau)*qgui1(ijb_u:ije_u,:)+tau*qgui2(ijb_u:ije_u,:)-q(ijb_u:ije_u,:)
        endif 
        if (guide_zon) CALL guide_zonave(2,jjp1,llm,f_add)
        CALL guide_addfield(ip1jmp1,llm,f_add,alpha_Q)
        IF (f_out) CALL guide_out("q",jjp1,llm,f_add(:,:)/factt,factt)
        q(ijb_u:ije_u,:)=q(ijb_u:ije_u,:)+f_add(ijb_u:ije_u,:)
    endif

    if (guide_v) then
        if (guide_add) then
           f_add(ijb_v:ije_v,:)=(1.-tau)*vgui1(ijb_v:ije_v,:)+tau*vgui2(ijb_v:ije_v,:)
        else
           f_add(ijb_v:ije_v,:)=(1.-tau)*vgui1(ijb_v:ije_v,:)+tau*vgui2(ijb_v:ije_v,:)-vcov(ijb_v:ije_v,:)
        endif 
        
        if (guide_zon) CALL guide_zonave(2,jjm,llm,f_add(1:ip1jm,:))
        CALL guide_addfield(ip1jm,llm,f_add(1:ip1jm,:),alpha_v)
        IF (f_out) CALL guide_out("v",jjm,llm,vcov(1:ip1jm,:),factt)
        IF (f_out) CALL guide_out("va",jjm,llm,(1.-tau)*vgui1(:,:)+tau*vgui2(:,:),factt)
        IF (f_out) CALL guide_out("vcov",jjm,llm,f_add(1:ip1jm,:)/factt,factt)
        vcov(ijb_v:ije_v,:)=vcov(ijb_v:ije_v,:)+f_add(ijb_v:ije_v,:)
    endif

  END SUBROUTINE guide_main

!=======================================================================
  SUBROUTINE guide_addfield(hsize,vsize,field,alpha)
! field1=a*field1+alpha*field2

    IMPLICIT NONE
    INCLUDE "dimensions.h"
    INCLUDE "paramet.h"

    ! input variables
    INTEGER,                      INTENT(IN)    :: hsize
    INTEGER,                      INTENT(IN)    :: vsize
    REAL, DIMENSION(hsize),       INTENT(IN)    :: alpha 
    REAL, DIMENSION(hsize,vsize), INTENT(INOUT) :: field

    ! Local variables
    INTEGER :: l

    IF (hsize==ip1jm) THEN
      do l=1,vsize
        field(ijb_v:ije_v,l)=alpha(ijb_v:ije_v)*field(ijb_v:ije_v,l)*alpha_pcor(l)
      enddo
    ELSE
      do l=1,vsize
        field(ijb_u:ije_u,l)=alpha(ijb_u:ije_u)*field(ijb_u:ije_u,l)*alpha_pcor(l)
      enddo
    ENDIF    

  END SUBROUTINE guide_addfield

!=======================================================================
  SUBROUTINE guide_zonave(typ,hsize,vsize,field)

    USE comconst_mod, ONLY: pi
    
    IMPLICIT NONE

    INCLUDE "dimensions.h"
    INCLUDE "paramet.h"
    INCLUDE "comgeom.h"
    
    ! input/output variables
    INTEGER,                           INTENT(IN)    :: typ
    INTEGER,                           INTENT(IN)    :: vsize
    INTEGER,                           INTENT(IN)    :: hsize
    REAL, DIMENSION(hsize*iip1,vsize), INTENT(INOUT) :: field

    ! Local variables
    LOGICAL, SAVE                :: first=.TRUE.
    INTEGER, DIMENSION (2), SAVE :: imin, imax ! averaging domain
    INTEGER                      :: i,j,l,ij
    REAL, DIMENSION (iip1)       :: lond       ! longitude in Deg.
    REAL, DIMENSION (hsize,vsize):: fieldm     ! zon-averaged field

    IF (first) THEN
        first=.FALSE.
!Compute domain for averaging
        lond=rlonu*180./pi
        imin(1)=1;imax(1)=iip1;
        imin(2)=1;imax(2)=iip1;
        IF (guide_reg) THEN
            DO i=1,iim
                IF (lond(i).LT.lon_min_g) imin(1)=i
                IF (lond(i).LE.lon_max_g) imax(1)=i
            ENDDO
            lond=rlonv*180./pi
            DO i=1,iim
                IF (lond(i).LT.lon_min_g) imin(2)=i
                IF (lond(i).LE.lon_max_g) imax(2)=i
            ENDDO
        ENDIF
    ENDIF

    fieldm=0.
    
    IF (hsize==jjm) THEN
      DO l=1,vsize
      ! Compute zonal average
          DO j=jjb_v,jje_v
              DO i=imin(typ),imax(typ)
                  ij=(j-1)*iip1+i
                  fieldm(j,l)=fieldm(j,l)+field(ij,l)
              ENDDO
          ENDDO 
          fieldm(:,l)=fieldm(:,l)/REAL(imax(typ)-imin(typ)+1)
    ! Compute forcing
          DO j=jjb_v,jje_v
              DO i=1,iip1
                  ij=(j-1)*iip1+i
                  field(ij,l)=fieldm(j,l)
              ENDDO
          ENDDO
      ENDDO
    ELSE
      DO l=1,vsize
      ! Compute zonal average
          DO j=jjb_v,jje_v
              DO i=imin(typ),imax(typ)
                  ij=(j-1)*iip1+i
                  fieldm(j,l)=fieldm(j,l)+field(ij,l)
              ENDDO
          ENDDO 
          fieldm(:,l)=fieldm(:,l)/REAL(imax(typ)-imin(typ)+1)
    ! Compute forcing
          DO j=jjb_u,jje_u
              DO i=1,iip1
                  ij=(j-1)*iip1+i
                  field(ij,l)=fieldm(j,l)
              ENDDO
          ENDDO
      ENDDO
    ENDIF    

  END SUBROUTINE guide_zonave

!=======================================================================
  SUBROUTINE guide_interp(psi,teta)
    use exner_hyb_p_m, only: exner_hyb_p
    use exner_milieu_p_m, only: exner_milieu_p
  USE parallel_lmdz
  USE mod_hallo
  USE Bands
  USE comconst_mod, ONLY: cpp, kappa
  USE comvert_mod, ONLY: preff, pressure_exner, bp, ap
  IMPLICIT NONE

  include "dimensions.h"
  include "paramet.h"
  include "comgeom2.h"

  REAL, DIMENSION (iip1,jjp1),     INTENT(IN) :: psi ! Psol gcm
  REAL, DIMENSION (iip1,jjp1,llm), INTENT(IN) :: teta ! Temp. Pot. gcm

  LOGICAL, SAVE                      :: first=.TRUE.
  ! Variables pour niveaux pression:
  REAL, DIMENSION (iip1,jjp1,nlevnc) :: plnc1,plnc2 !niveaux pression guidage
  REAL, DIMENSION (iip1,jjp1,llm)    :: plunc,plsnc !niveaux pression modele
  REAL, DIMENSION (iip1,jjm,llm)     :: plvnc       !niveaux pression modele
  REAL, DIMENSION (iip1,jjp1,llmp1)  :: p           ! pression intercouches 
  REAL, DIMENSION (iip1,jjp1,llm)    :: pls, pext   ! var intermediaire
  REAL, DIMENSION (iip1,jjp1,llm)    :: pbarx 
  REAL, DIMENSION (iip1,jjm,llm)     :: pbary 
  ! Variables pour fonction Exner (P milieu couche)
  REAL, DIMENSION (iip1,jjp1,llm)    :: pk
  REAL, DIMENSION (iip1,jjp1)        :: pks    
  REAL                               :: unskap
  ! Pression de vapeur saturante
  REAL, DIMENSION (ip1jmp1,llm)      :: qsat
  !Variables intermediaires interpolation
  REAL, DIMENSION (iip1,jjp1,llm)    :: zu1,zu2 
  REAL, DIMENSION (iip1,jjm,llm)     :: zv1,zv2
  
  INTEGER                            :: i,j,l,ij
  TYPE(Request) :: Req  

    print *,'Guide: conversion variables guidage'
! -----------------------------------------------------------------
! Calcul des niveaux de pression champs guidage (pour T et Q)
! -----------------------------------------------------------------
    IF (guide_plevs.EQ.0) THEN
        DO l=1,nlevnc
            DO j=jjb_u,jje_u
                DO i=1,iip1
                    plnc2(i,j,l)=apnc(l)
                    plnc1(i,j,l)=apnc(l)
               ENDDO
            ENDDO
        ENDDO
    ENDIF   

    if (first) then
        first=.FALSE.
        print*,'Guide: verification ordre niveaux verticaux'
        print*,'LMDZ :'
        do l=1,llm
            print*,'PL(',l,')=',(ap(l)+ap(l+1))/2. &
                  +psi(1,jje_u)*(bp(l)+bp(l+1))/2.
        enddo
        print*,'Fichiers guidage'
        SELECT CASE (guide_plevs)
        CASE (0) 
            do l=1,nlevnc
                 print*,'PL(',l,')=',plnc2(1,jjb_u,l)
            enddo
        CASE (1)
            DO l=1,nlevnc
                 print*,'PL(',l,')=',apnc(l)+bpnc(l)*psnat2(i,jjb_u)
             ENDDO
        CASE (2)
            do l=1,nlevnc
                 print*,'PL(',l,')=',pnat2(1,jjb_u,l)
            enddo
        END SELECT
        print *,'inversion de l''ordre: invert_p=',invert_p
        if (guide_u) then
            do l=1,nlevnc
                print*,'U(',l,')=',unat2(1,jjb_u,l)
            enddo
        endif
        if (guide_T) then
            do l=1,nlevnc
                print*,'T(',l,')=',tnat2(1,jjb_u,l)
            enddo
        endif
    endif
    
! -----------------------------------------------------------------
! Calcul niveaux pression modele 
! -----------------------------------------------------------------

!    ....  Calcul de pls , pression au milieu des couches ,en Pascals
    IF (guide_plevs.EQ.1) THEN
        DO l=1,llm
            DO j=jjb_u,jje_u
                DO i =1, iip1
                    pls(i,j,l)=(ap(l)+ap(l+1))/2.+psi(i,j)*(bp(l)+bp(l+1))/2.
                ENDDO
            ENDDO
        ENDDO
    ELSE
        CALL pression_p( ip1jmp1, ap, bp, psi, p )
        if (pressure_exner) then
          CALL exner_hyb_p(ip1jmp1,psi,p,pks,pk)
        else
          CALL exner_milieu_p(ip1jmp1,psi,p,pks,pk)
        endif
        unskap=1./kappa
        DO l = 1, llm
            DO j=jjb_u,jje_u
                DO i =1, iip1
                    pls(i,j,l) = preff * ( pk(i,j,l)/cpp) ** unskap
                ENDDO
            ENDDO
        ENDDO
    ENDIF

!   calcul des pressions pour les grilles u et v
    do l=1,llm
        do j=jjb_u,jje_u
            do i=1,iip1
                pext(i,j,l)=pls(i,j,l)*aire(i,j)
            enddo
        enddo
    enddo

     CALL Register_SwapFieldHallo(pext,pext,ip1jmp1,llm,jj_Nb_caldyn,1,2,Req)
     CALL SendRequest(Req)
     CALL WaitRequest(Req)

     call massbar_p(pext, pbarx, pbary )
    do l=1,llm
        do j=jjb_u,jje_u
            do i=1,iip1
                plunc(i,j,l)=pbarx(i,j,l)/aireu(i,j)
                plsnc(i,j,l)=pls(i,j,l)
            enddo
        enddo
    enddo
    do l=1,llm
        do j=jjb_v,jje_v
            do i=1,iip1
                plvnc(i,j,l)=pbary(i,j,l)/airev(i,j)
            enddo
        enddo
    enddo

! -----------------------------------------------------------------
! Interpolation verticale champs guidage sur niveaux modele
! Conversion en variables gcm (ucov, vcov...)
! -----------------------------------------------------------------
    if (guide_P) then
        do j=jjb_u,jje_u
            do i=1,iim
                ij=(j-1)*iip1+i
                psgui1(ij)=psnat1(i,j)
                psgui2(ij)=psnat2(i,j)
            enddo
            psgui1(iip1*j)=psnat1(1,j)
            psgui2(iip1*j)=psnat2(1,j)
        enddo
    endif

    IF (guide_T) THEN
        ! Calcul des nouvelles valeurs des niveaux de pression du guidage
        IF (guide_plevs.EQ.1) THEN
            DO l=1,nlevnc
                DO j=jjb_u,jje_u
                    DO i=1,iip1
                        plnc2(i,j,l)=apnc(l)+bpnc(l)*psnat2(i,j)
                        plnc1(i,j,l)=apnc(l)+bpnc(l)*psnat1(i,j)
                    ENDDO
                ENDDO
            ENDDO
        ELSE IF (guide_plevs.EQ.2) THEN
            DO l=1,nlevnc
                DO j=jjb_u,jje_u
                    DO i=1,iip1
                        plnc2(i,j,l)=pnat2(i,j,l)
                        plnc1(i,j,l)=pnat1(i,j,l)
                    ENDDO
                ENDDO
            ENDDO
        ENDIF

        ! Interpolation verticale
        CALL pres2lev(tnat1(:,jjb_u:jje_u,:),zu1(:,jjb_u:jje_u,:),nlevnc,llm,           &
                    plnc1(:,jjb_u:jje_u,:),plsnc(:,jjb_u:jje_u,:),iip1,jjn_u,invert_p)
        CALL pres2lev(tnat2(:,jjb_u:jje_u,:),zu2(:,jjb_u:jje_u,:),nlevnc,llm,           &
                    plnc2(:,jjb_u:jje_u,:),plsnc(:,jjb_u:jje_u,:),iip1,jjn_u,invert_p)

        ! Conversion en variables GCM
        do l=1,llm
            do j=jjb_u,jje_u
                IF (guide_teta) THEN
                    do i=1,iim
                        ij=(j-1)*iip1+i
                        tgui1(ij,l)=zu1(i,j,l)
                        tgui2(ij,l)=zu2(i,j,l)
                    enddo
                ELSE
                    do i=1,iim
                        ij=(j-1)*iip1+i
                        tgui1(ij,l)=zu1(i,j,l)*cpp/pk(i,j,l)
                        tgui2(ij,l)=zu2(i,j,l)*cpp/pk(i,j,l)
                    enddo
                ENDIF
                tgui1(j*iip1,l)=tgui1((j-1)*iip1+1,l)    
                tgui2(j*iip1,l)=tgui2((j-1)*iip1+1,l)    
            enddo
            do i=1,iip1
                tgui1(i,l)=tgui1(1,l)
                tgui1(ip1jm+i,l)=tgui1(ip1jm+1,l) 
                tgui2(i,l)=tgui2(1,l)
                tgui2(ip1jm+i,l)=tgui2(ip1jm+1,l) 
            enddo
        enddo
    ENDIF

    IF (guide_Q) THEN
        ! Calcul des nouvelles valeurs des niveaux de pression du guidage
        IF (guide_plevs.EQ.1) THEN
            DO l=1,nlevnc
                DO j=jjb_u,jje_u
                    DO i=1,iip1
                        plnc2(i,j,l)=apnc(l)+bpnc(l)*psnat2(i,j)
                        plnc1(i,j,l)=apnc(l)+bpnc(l)*psnat1(i,j)
                    ENDDO
                ENDDO
            ENDDO
        ELSE IF (guide_plevs.EQ.2) THEN
            DO l=1,nlevnc
                DO j=jjb_u,jje_u
                    DO i=1,iip1
                        plnc2(i,j,l)=pnat2(i,j,l)
                        plnc1(i,j,l)=pnat1(i,j,l)
                    ENDDO
                ENDDO
            ENDDO
        ENDIF

        ! Interpolation verticale
        CALL pres2lev(qnat1(:,jjb_u:jje_u,:),zu1(:,jjb_u:jje_u,:),nlevnc,llm,             &
                      plnc1(:,jjb_u:jje_u,:),plsnc(:,jjb_u:jje_u,:),iip1,jjn_u,invert_p)
        CALL pres2lev(qnat2(:,jjb_u:jje_u,:),zu2(:,jjb_u:jje_u,:),nlevnc,llm,             &
                      plnc2(:,jjb_u:jje_u,:),plsnc(:,jjb_u:jje_u,:),iip1,jjn_u,invert_p)

        ! Conversion en variables GCM
        ! On suppose qu'on a la bonne variable dans le fichier de guidage:
        ! Hum.Rel si guide_hr, Hum.Spec. sinon.
        do l=1,llm
            do j=jjb_u,jje_u
                do i=1,iim
                    ij=(j-1)*iip1+i
                    qgui1(ij,l)=zu1(i,j,l)
                    qgui2(ij,l)=zu2(i,j,l)
                enddo
                qgui1(j*iip1,l)=qgui1((j-1)*iip1+1,l)    
                qgui2(j*iip1,l)=qgui2((j-1)*iip1+1,l)    
            enddo
            do i=1,iip1
                qgui1(i,l)=qgui1(1,l)
                qgui1(ip1jm+i,l)=qgui1(ip1jm+1,l) 
                qgui2(i,l)=qgui2(1,l)
                qgui2(ip1jm+i,l)=qgui2(ip1jm+1,l) 
            enddo
        enddo
        IF (guide_hr) THEN
            CALL q_sat(iip1*jjn_u*llm,teta(:,jjb_u:jje_u,:)*pk(:,jjb_u:jje_u,:)/cpp,       &
                       plsnc(:,jjb_u:jje_u,:),qsat(ijb_u:ije_u,:))
            qgui1(ijb_u:ije_u,:)=qgui1(ijb_u:ije_u,:)*qsat(ijb_u:ije_u,:)*0.01 !hum. rel. en %
            qgui2(ijb_u:ije_u,:)=qgui2(ijb_u:ije_u,:)*qsat(ijb_u:ije_u,:)*0.01 
        ENDIF
    ENDIF

    IF (guide_u) THEN
        ! Calcul des nouvelles valeurs des niveaux de pression du guidage
        IF (guide_plevs.EQ.1) THEN
            DO l=1,nlevnc
                DO j=jjb_u,jje_u
                    DO i=1,iim
                        plnc2(i,j,l)=apnc(l)+bpnc(l)*(psnat2(i,j)*aire(i,j)*alpha1p2(i,j) &
                       &           +psnat2(i+1,j)*aire(i+1,j)*alpha3p4(i+1,j))/aireu(i,j)
                        plnc1(i,j,l)=apnc(l)+bpnc(l)*(psnat1(i,j)*aire(i,j)*alpha1p2(i,j) &
                       &           +psnat1(i+1,j)*aire(i+1,j)*alpha3p4(i+1,j))/aireu(i,j)
                    ENDDO
                    plnc2(iip1,j,l)=plnc2(1,j,l)
                    plnc1(iip1,j,l)=plnc1(1,j,l)
                ENDDO
            ENDDO
        ELSE IF (guide_plevs.EQ.2) THEN
            DO l=1,nlevnc
                DO j=jjb_u,jje_u
                    DO i=1,iim
                        plnc2(i,j,l)=(pnat2(i,j,l)*aire(i,j)*alpha1p2(i,j) &
                       & +pnat2(i+1,j,l)*aire(i,j)*alpha3p4(i+1,j))/aireu(i,j)
                        plnc1(i,j,l)=(pnat1(i,j,l)*aire(i,j)*alpha1p2(i,j) &
                       & +pnat1(i+1,j,l)*aire(i,j)*alpha3p4(i+1,j))/aireu(i,j)
                    ENDDO
                    plnc2(iip1,j,l)=plnc2(1,j,l)
                    plnc1(iip1,j,l)=plnc1(1,j,l)
                ENDDO
            ENDDO
        ENDIF
        
        ! Interpolation verticale
        CALL pres2lev(unat1(:,jjb_u:jje_u,:),zu1(:,jjb_u:jje_u,:),nlevnc,llm,            &
                      plnc1(:,jjb_u:jje_u,:),plunc(:,jjb_u:jje_u,:),iip1,jjn_u,invert_p)
        CALL pres2lev(unat2(:,jjb_u:jje_u,:),zu2(:,jjb_u:jje_u,:),nlevnc,llm,            &
                      plnc2(:,jjb_u:jje_u,:),plunc(:,jjb_u:jje_u,:),iip1,jjn_u,invert_p)

        ! Conversion en variables GCM
        do l=1,llm
            do j=jjb_u,jje_u
                do i=1,iim
                    ij=(j-1)*iip1+i
                    ugui1(ij,l)=zu1(i,j,l)*cu(i,j)
                    ugui2(ij,l)=zu2(i,j,l)*cu(i,j)
                enddo
                ugui1(j*iip1,l)=ugui1((j-1)*iip1+1,l)    
                ugui2(j*iip1,l)=ugui2((j-1)*iip1+1,l)    
            enddo
            do i=1,iip1
                ugui1(i,l)=0.
                ugui1(ip1jm+i,l)=0.
                ugui2(i,l)=0.
                ugui2(ip1jm+i,l)=0.
            enddo
        enddo
    ENDIF
    
    IF (guide_v) THEN
        ! Calcul des nouvelles valeurs des niveaux de pression du guidage
        IF (guide_plevs.EQ.1) THEN
         CALL Register_SwapFieldHallo(psnat1,psnat1,ip1jmp1,1,jj_Nb_caldyn,1,2,Req)
         CALL SendRequest(Req)
         CALL WaitRequest(Req)
         CALL Register_SwapFieldHallo(psnat2,psnat2,ip1jmp1,1,jj_Nb_caldyn,1,2,Req)
         CALL SendRequest(Req)
         CALL WaitRequest(Req)
            DO l=1,nlevnc
                DO j=jjb_v,jje_v
                    DO i=1,iip1
                        plnc2(i,j,l)=apnc(l)+bpnc(l)*(psnat2(i,j)*aire(i,j)*alpha2p3(i,j) &
                       &           +psnat2(i,j+1)*aire(i,j+1)*alpha1p4(i,j+1))/airev(i,j)
                        plnc1(i,j,l)=apnc(l)+bpnc(l)*(psnat1(i,j)*aire(i,j)*alpha2p3(i,j) &
                       &           +psnat1(i,j+1)*aire(i,j+1)*alpha1p4(i,j+1))/airev(i,j)
                    ENDDO
                ENDDO
            ENDDO
        ELSE IF (guide_plevs.EQ.2) THEN
         CALL Register_SwapFieldHallo(pnat1,pnat1,ip1jmp1,llm,jj_Nb_caldyn,1,2,Req)
         CALL SendRequest(Req)
         CALL WaitRequest(Req)
         CALL Register_SwapFieldHallo(pnat2,pnat2,ip1jmp1,llm,jj_Nb_caldyn,1,2,Req)
         CALL SendRequest(Req)
         CALL WaitRequest(Req)
            DO l=1,nlevnc
                DO j=jjb_v,jje_v
                    DO i=1,iip1
                        plnc2(i,j,l)=(pnat2(i,j,l)*aire(i,j)*alpha2p3(i,j) &
                       & +pnat2(i,j+1,l)*aire(i,j)*alpha1p4(i,j+1))/airev(i,j)
                        plnc1(i,j,l)=(pnat1(i,j,l)*aire(i,j)*alpha2p3(i,j) &
                       & +pnat1(i,j+1,l)*aire(i,j)*alpha1p4(i,j+1))/airev(i,j)
                    ENDDO
                ENDDO
            ENDDO
        ENDIF
        ! Interpolation verticale
        CALL pres2lev(vnat1(:,jjb_v:jje_v,:),zv1(:,jjb_v:jje_v,:),nlevnc,llm,             &
                      plnc1(:,jjb_v:jje_v,:),plvnc(:,jjb_v:jje_v,:),iip1,jjn_v,invert_p)
        CALL pres2lev(vnat2(:,jjb_v:jje_v,:),zv2(:,jjb_v:jje_v,:),nlevnc,llm,             &
                      plnc2(:,jjb_v:jje_v,:),plvnc(:,jjb_v:jje_v,:),iip1,jjn_v,invert_p)
        ! Conversion en variables GCM
        do l=1,llm
            do j=jjb_v,jje_v
                do i=1,iim
                    ij=(j-1)*iip1+i
                    vgui1(ij,l)=zv1(i,j,l)*cv(i,j)
                    vgui2(ij,l)=zv2(i,j,l)*cv(i,j)
                enddo
                vgui1(j*iip1,l)=vgui1((j-1)*iip1+1,l)    
                vgui2(j*iip1,l)=vgui2((j-1)*iip1+1,l)    
            enddo
        enddo
    ENDIF
    

  END SUBROUTINE guide_interp

!=======================================================================
  SUBROUTINE tau2alpha(typ,pim,pjm,factt,taumin,taumax,alpha)

! Calcul des constantes de rappel alpha (=1/tau)

    use comconst_mod, only: pi
    use serre_mod, only: clat, clon, grossismx, grossismy
    
    implicit none

    include "dimensions.h"
    include "paramet.h"
    include "comgeom2.h"

! input arguments :
    INTEGER, INTENT(IN) :: typ    ! u(2),v(3), ou scalaire(1)
    INTEGER, INTENT(IN) :: pim,pjm ! dimensions en lat, lon
    REAL, INTENT(IN)    :: factt   ! pas de temps en fraction de jour
    REAL, INTENT(IN)    :: taumin,taumax
! output arguments:
    REAL, DIMENSION(pim,pjm), INTENT(OUT) :: alpha 
  
!  local variables:
    LOGICAL, SAVE               :: first=.TRUE.
    REAL, SAVE                  :: gamma,dxdy_min,dxdy_max
    REAL, DIMENSION (iip1,jjp1) :: zdx,zdy
    REAL, DIMENSION (iip1,jjp1) :: dxdys,dxdyu
    REAL, DIMENSION (iip1,jjm)  :: dxdyv
    real dxdy_
    real zlat,zlon
    real alphamin,alphamax,xi
    integer i,j,ilon,ilat


    alphamin=factt/taumax
    alphamax=factt/taumin
    IF (guide_reg.OR.guide_add) THEN
        alpha=alphamax
!-----------------------------------------------------------------------
! guide_reg: alpha=alpha_min dans region, 0. sinon.
!-----------------------------------------------------------------------
        IF (guide_reg) THEN
            do j=1,pjm
                do i=1,pim
                    if (typ.eq.2) then
                       zlat=rlatu(j)*180./pi
                       zlon=rlonu(i)*180./pi
                    elseif (typ.eq.1) then
                       zlat=rlatu(j)*180./pi
                       zlon=rlonv(i)*180./pi
                    elseif (typ.eq.3) then
                       zlat=rlatv(j)*180./pi
                       zlon=rlonv(i)*180./pi
                    endif
                    alpha(i,j)=alphamax/16.* &
                              (1.+tanh((zlat-lat_min_g)/tau_lat))* &
                              (1.+tanh((lat_max_g-zlat)/tau_lat))* &
                              (1.+tanh((zlon-lon_min_g)/tau_lon))* &
                              (1.+tanh((lon_max_g-zlon)/tau_lon))
                enddo
            enddo
        ENDIF
    ELSE
!-----------------------------------------------------------------------
! Sinon, alpha varie entre alpha_min et alpha_max suivant le zoom.
!-----------------------------------------------------------------------
!Calcul de l'aire des mailles
        do j=2,jjm
            do i=2,iip1
               zdx(i,j)=0.5*(cu(i-1,j)+cu(i,j))/cos(rlatu(j))
            enddo
            zdx(1,j)=zdx(iip1,j)
        enddo
        do j=2,jjm
            do i=1,iip1
               zdy(i,j)=0.5*(cv(i,j-1)+cv(i,j))
            enddo
        enddo
        do i=1,iip1
            zdx(i,1)=zdx(i,2)
            zdx(i,jjp1)=zdx(i,jjm)
            zdy(i,1)=zdy(i,2)
            zdy(i,jjp1)=zdy(i,jjm)
        enddo
        do j=1,jjp1
            do i=1,iip1
               dxdys(i,j)=sqrt(zdx(i,j)*zdx(i,j)+zdy(i,j)*zdy(i,j))
            enddo
        enddo
        IF (typ.EQ.2) THEN
            do j=1,jjp1
                do i=1,iim
                   dxdyu(i,j)=0.5*(dxdys(i,j)+dxdys(i+1,j))
                enddo
                dxdyu(iip1,j)=dxdyu(1,j)
            enddo
        ENDIF
        IF (typ.EQ.3) THEN
            do j=1,jjm
                do i=1,iip1
                   dxdyv(i,j)=0.5*(dxdys(i,j)+dxdys(i,j+1))
                enddo
            enddo
        ENDIF
! Premier appel: calcul des aires min et max et de gamma.
        IF (first) THEN 
            first=.FALSE.
            ! coordonnees du centre du zoom
            CALL coordij(clon,clat,ilon,ilat) 
            ! aire de la maille au centre du zoom
            dxdy_min=dxdys(ilon,ilat)
            ! dxdy maximale de la maille
            dxdy_max=0.
            do j=1,jjp1
                do i=1,iip1
                     dxdy_max=max(dxdy_max,dxdys(i,j))
                enddo
            enddo
            ! Calcul de gamma
            if (abs(grossismx-1.).lt.0.1.or.abs(grossismy-1.).lt.0.1) then
                 print*,'ATTENTION modele peu zoome'
                 print*,'ATTENTION on prend une constante de guidage cste'
                 gamma=0.
            else
                gamma=(dxdy_max-2.*dxdy_min)/(dxdy_max-dxdy_min)
                print*,'gamma=',gamma
                if (gamma.lt.1.e-5) then
                  print*,'gamma =',gamma,'<1e-5'
                  stop
                endif
                gamma=log(0.5)/log(gamma)
                if (gamma4) then 
                  gamma=min(gamma,4.)
                endif
                print*,'gamma=',gamma
            endif
        ENDIF !first

        do j=1,pjm
            do i=1,pim
                if (typ.eq.1) then
                   dxdy_=dxdys(i,j)
                   zlat=rlatu(j)*180./pi
                elseif (typ.eq.2) then
                   dxdy_=dxdyu(i,j)
                   zlat=rlatu(j)*180./pi
                elseif (typ.eq.3) then
                   dxdy_=dxdyv(i,j)
                   zlat=rlatv(j)*180./pi
                endif
                if (abs(grossismx-1.).lt.0.1.or.abs(grossismy-1.).lt.0.1) then
                ! pour une grille reguliere, xi=xxx**0=1 -> alpha=alphamin
                    alpha(i,j)=alphamin
                else
                    xi=((dxdy_max-dxdy_)/(dxdy_max-dxdy_min))**gamma
                    xi=min(xi,1.)
                    if(lat_min_g.le.zlat .and. zlat.le.lat_max_g) then
                        alpha(i,j)=xi*alphamin+(1.-xi)*alphamax
                    else
                        alpha(i,j)=0.
                    endif
                endif
            enddo
        enddo
    ENDIF ! guide_reg

    if (.not. guide_add) alpha = 1. - exp(- alpha)

  END SUBROUTINE tau2alpha

!=======================================================================
  SUBROUTINE guide_read(timestep)

    IMPLICIT NONE

#include "netcdf.inc"
#include "dimensions.h"
#include "paramet.h"

    INTEGER, INTENT(IN)   :: timestep

    LOGICAL, SAVE         :: first=.TRUE.
! Identification fichiers et variables NetCDF:
    INTEGER, SAVE         :: ncidu,varidu,ncidv,varidv,ncidp,varidp
    INTEGER, SAVE         :: ncidQ,varidQ,ncidt,varidt,ncidps,varidps
    INTEGER               :: ncidpl,varidpl,varidap,varidbp
! Variables auxiliaires NetCDF:
    INTEGER, DIMENSION(4) :: start,count
    INTEGER               :: status,rcode

    CHARACTER (len = 80)   :: abort_message
    CHARACTER (len = 20)   :: modname = 'guide_read'
! -----------------------------------------------------------------
! Premier appel: initialisation de la lecture des fichiers
! -----------------------------------------------------------------
    if (first) then
         ncidpl=-99
         print*,'Guide: ouverture des fichiers guidage '
! Ap et Bp si Niveaux de pression hybrides
         if (guide_plevs.EQ.1) then
             print *,'Lecture du guidage sur niveaux modele'
             rcode = nf90_open('apbp.nc', nf90_nowrite, ncidpl)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de fichier apbp.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidpl, 'AP', varidap)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de variable AP, fichier apbp.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidpl, 'BP', varidbp)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de variable BP, fichier apbp.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             print*,'ncidpl,varidap',ncidpl,varidap
         endif
! Pression si guidage sur niveaux P variables
         if (guide_plevs.EQ.2) then
             rcode = nf90_open('P.nc', nf90_nowrite, ncidp)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de fichier P.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidp, 'PRES', varidp)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de variable PRES, fichier P.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             print*,'ncidp,varidp',ncidp,varidp
             if (ncidpl.eq.-99) ncidpl=ncidp
         endif
! Vent zonal
         if (guide_u) then
             rcode = nf90_open('u.nc', nf90_nowrite, ncidu)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de fichier u.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF 
             rcode = nf90_inq_varid(ncidu, 'UWND', varidu)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de variable UWND, fichier u.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             print*,'ncidu,varidu',ncidu,varidu
             if (ncidpl.eq.-99) ncidpl=ncidu
         endif
! Vent meridien
         if (guide_v) then
             rcode = nf90_open('v.nc', nf90_nowrite, ncidv)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de fichier v.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidv, 'VWND', varidv)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de variable VWND, fichier v.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             print*,'ncidv,varidv',ncidv,varidv
             if (ncidpl.eq.-99) ncidpl=ncidv
         endif
! Temperature
         if (guide_T) then
             rcode = nf90_open('T.nc', nf90_nowrite, ncidt)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de fichier T.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidt, 'AIR', varidt)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de variable AIR, fichier T.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             print*,'ncidT,varidT',ncidt,varidt
             if (ncidpl.eq.-99) ncidpl=ncidt
         endif
! Humidite
         if (guide_Q) then
             rcode = nf90_open('hur.nc', nf90_nowrite, ncidQ)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de fichier hur.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidQ, 'RH', varidQ)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de variable RH, fichier hur.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             print*,'ncidQ,varidQ',ncidQ,varidQ
             if (ncidpl.eq.-99) ncidpl=ncidQ
         endif
! Pression de surface
         if ((guide_P).OR.(guide_plevs.EQ.1)) then
             rcode = nf90_open('ps.nc', nf90_nowrite, ncidps)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de fichier ps.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidps, 'SP', varidps)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de variable SP, fichier ps.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             print*,'ncidps,varidps',ncidps,varidps
         endif
! Coordonnee verticale
         if (guide_plevs.EQ.0) then
              rcode = nf90_inq_varid(ncidpl, 'LEVEL', varidpl)
              IF (rcode.NE.0) rcode = nf90_inq_varid(ncidpl, 'PRESSURE', varidpl)
              print*,'ncidpl,varidpl',ncidpl,varidpl
         endif
! Coefs ap, bp pour calcul de la pression aux differents niveaux
         IF (guide_plevs.EQ.1) THEN
#ifdef NC_DOUBLE
             status=NF_GET_VARA_DOUBLE(ncidpl,varidap,1,nlevnc,apnc)
             status=NF_GET_VARA_DOUBLE(ncidpl,varidbp,1,nlevnc,bpnc)
#else
             status=NF_GET_VARA_REAL(ncidpl,varidap,1,nlevnc,apnc)
             status=NF_GET_VARA_REAL(ncidpl,varidbp,1,nlevnc,bpnc)
#endif
         ELSEIF (guide_plevs.EQ.0) THEN
#ifdef NC_DOUBLE
             status=NF_GET_VARA_DOUBLE(ncidpl,varidpl,1,nlevnc,apnc)
#else
             status=NF_GET_VARA_REAL(ncidpl,varidpl,1,nlevnc,apnc)
#endif
             apnc=apnc*100.! conversion en Pascals
             bpnc(:)=0.
         ENDIF
         first=.FALSE.
     ENDIF ! (first)

! -----------------------------------------------------------------
!   lecture des champs u, v, T, Q, ps
! -----------------------------------------------------------------

!  dimensions pour les champs scalaires et le vent zonal
     start(1)=1
     start(2)=1
     start(3)=1
     start(4)=timestep

     count(1)=iip1
     count(2)=jjp1
     count(3)=nlevnc
     count(4)=1

! Pression 
     if (guide_plevs.EQ.2) then
#ifdef NC_DOUBLE
         status=NF_GET_VARA_DOUBLE(ncidp,varidp,start,count,pnat2)
#else
         status=NF_GET_VARA_REAL(ncidp,varidp,start,count,pnat2)
#endif
         IF (invert_y) THEN
           CALL invert_lat(iip1,jjp1,nlevnc,pnat2)
         ENDIF
     endif

!  Vent zonal
     if (guide_u) then
#ifdef NC_DOUBLE
         status=NF_GET_VARA_DOUBLE(ncidu,varidu,start,count,unat2)
#else
         status=NF_GET_VARA_REAL(ncidu,varidu,start,count,unat2)
#endif
         IF (invert_y) THEN
           CALL invert_lat(iip1,jjp1,nlevnc,unat2)
         ENDIF

     endif

!  Temperature
     if (guide_T) then
#ifdef NC_DOUBLE
         status=NF_GET_VARA_DOUBLE(ncidt,varidt,start,count,tnat2)
#else
         status=NF_GET_VARA_REAL(ncidt,varidt,start,count,tnat2)
#endif
         IF (invert_y) THEN
           CALL invert_lat(iip1,jjp1,nlevnc,tnat2)
         ENDIF
     endif

!  Humidite
     if (guide_Q) then
#ifdef NC_DOUBLE
         status=NF_GET_VARA_DOUBLE(ncidQ,varidQ,start,count,qnat2)
#else
         status=NF_GET_VARA_REAL(ncidQ,varidQ,start,count,qnat2)
#endif
         IF (invert_y) THEN
           CALL invert_lat(iip1,jjp1,nlevnc,qnat2)
         ENDIF

     endif

!  Vent meridien
     if (guide_v) then
         count(2)=jjm
#ifdef NC_DOUBLE
         status=NF_GET_VARA_DOUBLE(ncidv,varidv,start,count,vnat2)
#else
         status=NF_GET_VARA_REAL(ncidv,varidv,start,count,vnat2)
#endif
         IF (invert_y) THEN
           CALL invert_lat(iip1,jjm,nlevnc,vnat2)
         ENDIF
     endif

!  Pression de surface
     if ((guide_P).OR.(guide_plevs.EQ.1))  then
         start(3)=timestep
         start(4)=0
         count(2)=jjp1
         count(3)=1
         count(4)=0
#ifdef NC_DOUBLE
         status=NF_GET_VARA_DOUBLE(ncidps,varidps,start,count,psnat2)
#else
         status=NF_GET_VARA_REAL(ncidps,varidps,start,count,psnat2)
#endif
         IF (invert_y) THEN
           CALL invert_lat(iip1,jjp1,1,psnat2)
         ENDIF
     endif

  END SUBROUTINE guide_read

!=======================================================================
  SUBROUTINE guide_read2D(timestep)

    IMPLICIT NONE

#include "netcdf.inc"
#include "dimensions.h"
#include "paramet.h"

    INTEGER, INTENT(IN)   :: timestep

    LOGICAL, SAVE         :: first=.TRUE.
! Identification fichiers et variables NetCDF:
    INTEGER, SAVE         :: ncidu,varidu,ncidv,varidv,ncidp,varidp
    INTEGER, SAVE         :: ncidQ,varidQ,ncidt,varidt,ncidps,varidps
    INTEGER               :: ncidpl,varidpl,varidap,varidbp
! Variables auxiliaires NetCDF:
    INTEGER, DIMENSION(4) :: start,count
    INTEGER               :: status,rcode
! Variables for 3D extension:
    REAL, DIMENSION (jjp1,llm) :: zu
    REAL, DIMENSION (jjm,llm)  :: zv
    INTEGER               :: i

    CHARACTER (len = 80)   :: abort_message
    CHARACTER (len = 20)   :: modname = 'guide_read2D'
! -----------------------------------------------------------------
! Premier appel: initialisation de la lecture des fichiers
! -----------------------------------------------------------------
    if (first) then
         ncidpl=-99
         print*,'Guide: ouverture des fichiers guidage '
! Ap et Bp si niveaux de pression hybrides
         if (guide_plevs.EQ.1) then
             print *,'Lecture du guidage sur niveaux mod�le'
             rcode = nf90_open('apbp.nc', nf90_nowrite, ncidpl)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de fichier apbp.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidpl, 'AP', varidap)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de variable AP, fichier apbp.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidpl, 'BP', varidbp)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de variable BP, fichier apbp.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             print*,'ncidpl,varidap',ncidpl,varidap
         endif
! Pression
         if (guide_plevs.EQ.2) then
             rcode = nf90_open('P.nc', nf90_nowrite, ncidp)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de fichier P.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidp, 'PRES', varidp)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de variable PRES, fichier P.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             print*,'ncidp,varidp',ncidp,varidp
             if (ncidpl.eq.-99) ncidpl=ncidp
         endif
! Vent zonal
         if (guide_u) then
             rcode = nf90_open('u.nc', nf90_nowrite, ncidu)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de fichier u.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidu, 'UWND', varidu)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de variable UWND, fichier u.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             print*,'ncidu,varidu',ncidu,varidu
             if (ncidpl.eq.-99) ncidpl=ncidu
         endif
! Vent meridien
         if (guide_v) then
             rcode = nf90_open('v.nc', nf90_nowrite, ncidv)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de fichier v.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF 
             rcode = nf90_inq_varid(ncidv, 'VWND', varidv)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de variable VWND, fichier v.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             print*,'ncidv,varidv',ncidv,varidv
             if (ncidpl.eq.-99) ncidpl=ncidv
         endif
! Temperature
         if (guide_T) then
             rcode = nf90_open('T.nc', nf90_nowrite, ncidt)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de fichier T.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidt, 'AIR', varidt)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de variable AIR, fichier T.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             print*,'ncidT,varidT',ncidt,varidt
             if (ncidpl.eq.-99) ncidpl=ncidt
         endif
! Humidite
         if (guide_Q) then
             rcode = nf90_open('hur.nc', nf90_nowrite, ncidQ)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de fichier hur.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidQ, 'RH', varidQ)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de variable RH, fichier hur.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             print*,'ncidQ,varidQ',ncidQ,varidQ
             if (ncidpl.eq.-99) ncidpl=ncidQ
         endif
! Pression de surface
         if ((guide_P).OR.(guide_plevs.EQ.1)) then
             rcode = nf90_open('ps.nc', nf90_nowrite, ncidps)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de fichier ps.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             rcode = nf90_inq_varid(ncidps, 'SP', varidps)
             IF (rcode.NE.NF_NOERR) THEN
              print *,'Guide: probleme -> pas de variable SP, fichier ps.nc'
              CALL abort_gcm(modname,abort_message,1)
             ENDIF
             print*,'ncidps,varidps',ncidps,varidps
         endif
! Coordonnee verticale
         if (guide_plevs.EQ.0) then
              rcode = nf90_inq_varid(ncidpl, 'LEVEL', varidpl)
              IF (rcode.NE.0) rcode = nf90_inq_varid(ncidpl, 'PRESSURE', varidpl)
              print*,'ncidpl,varidpl',ncidpl,varidpl
         endif
! Coefs ap, bp pour calcul de la pression aux differents niveaux
         if (guide_plevs.EQ.1) then
#ifdef NC_DOUBLE
             status=NF_GET_VARA_DOUBLE(ncidpl,varidap,1,nlevnc,apnc)
             status=NF_GET_VARA_DOUBLE(ncidpl,varidbp,1,nlevnc,bpnc)
#else
             status=NF_GET_VARA_REAL(ncidpl,varidap,1,nlevnc,apnc)
             status=NF_GET_VARA_REAL(ncidpl,varidbp,1,nlevnc,bpnc)
#endif
         elseif (guide_plevs.EQ.0) THEN
#ifdef NC_DOUBLE
             status=NF_GET_VARA_DOUBLE(ncidpl,varidpl,1,nlevnc,apnc)
#else
             status=NF_GET_VARA_REAL(ncidpl,varidpl,1,nlevnc,apnc)
#endif
             apnc=apnc*100.! conversion en Pascals
             bpnc(:)=0.
         endif
         first=.FALSE.
     endif ! (first)

! -----------------------------------------------------------------
!   lecture des champs u, v, T, Q, ps
! -----------------------------------------------------------------

!  dimensions pour les champs scalaires et le vent zonal
     start(1)=1
     start(2)=1
     start(3)=1
     start(4)=timestep

     count(1)=1
     count(2)=jjp1
     count(3)=nlevnc
     count(4)=1

!  Pression
     if (guide_plevs.EQ.2) then
#ifdef NC_DOUBLE
         status=NF_GET_VARA_DOUBLE(ncidp,varidp,start,count,zu)
#else
         status=NF_GET_VARA_REAL(ncidp,varidp,start,count,zu)
#endif
         DO i=1,iip1
             pnat2(i,:,:)=zu(:,:)
         ENDDO

         IF (invert_y) THEN
           CALL invert_lat(iip1,jjp1,nlevnc,pnat2)
         ENDIF
     endif
!  Vent zonal
     if (guide_u) then
#ifdef NC_DOUBLE
         status=NF_GET_VARA_DOUBLE(ncidu,varidu,start,count,zu)
#else
         status=NF_GET_VARA_REAL(ncidu,varidu,start,count,zu)
#endif
         DO i=1,iip1
             unat2(i,:,:)=zu(:,:)
         ENDDO

         IF (invert_y) THEN
           CALL invert_lat(iip1,jjp1,nlevnc,unat2)
         ENDIF
     endif

!  Temperature
     if (guide_T) then
#ifdef NC_DOUBLE
         status=NF_GET_VARA_DOUBLE(ncidt,varidt,start,count,zu)
#else
         status=NF_GET_VARA_REAL(ncidt,varidt,start,count,zu)
#endif
         DO i=1,iip1
             tnat2(i,:,:)=zu(:,:)
         ENDDO

         IF (invert_y) THEN
           CALL invert_lat(iip1,jjp1,nlevnc,tnat2)
         ENDIF
     endif

!  Humidite
     if (guide_Q) then
#ifdef NC_DOUBLE
         status=NF_GET_VARA_DOUBLE(ncidQ,varidQ,start,count,zu)
#else
         status=NF_GET_VARA_REAL(ncidQ,varidQ,start,count,zu)
#endif
         DO i=1,iip1
             qnat2(i,:,:)=zu(:,:)
         ENDDO
         
         IF (invert_y) THEN
           CALL invert_lat(iip1,jjp1,nlevnc,qnat2)
         ENDIF
     endif

!  Vent meridien
     if (guide_v) then
         count(2)=jjm
#ifdef NC_DOUBLE
         status=NF_GET_VARA_DOUBLE(ncidv,varidv,start,count,zv)
#else
         status=NF_GET_VARA_REAL(ncidv,varidv,start,count,zv)
#endif
         DO i=1,iip1
             vnat2(i,:,:)=zv(:,:)
         ENDDO

         IF (invert_y) THEN
           CALL invert_lat(iip1,jjm,nlevnc,vnat2)
         ENDIF
     endif

!  Pression de surface
     if ((guide_P).OR.(guide_plevs.EQ.1))  then
         start(3)=timestep
         start(4)=0
         count(2)=jjp1
         count(3)=1
         count(4)=0
#ifdef NC_DOUBLE
         status=NF_GET_VARA_DOUBLE(ncidps,varidps,start,count,zu(:,1))
#else
         status=NF_GET_VARA_REAL(ncidps,varidps,start,count,zu(:,1))
#endif
         DO i=1,iip1
             psnat2(i,:)=zu(:,1)
         ENDDO

         IF (invert_y) THEN
           CALL invert_lat(iip1,jjp1,1,psnat2)
         ENDIF
     endif

  END SUBROUTINE guide_read2D
  
!=======================================================================
  SUBROUTINE guide_out(varname,hsize,vsize,field,factt)
    USE parallel_lmdz
    USE comconst_mod, ONLY: pi
    USE comvert_mod, ONLY: presnivs
    use netcdf95, only: nf95_def_var, nf95_put_var
    use netcdf, only: nf90_float

    IMPLICIT NONE

    INCLUDE "dimensions.h"
    INCLUDE "paramet.h"
    INCLUDE "netcdf.inc"
    INCLUDE "comgeom2.h"
    
    ! Variables entree
    CHARACTER*(*), INTENT(IN)                          :: varname
    INTEGER,   INTENT (IN)                         :: hsize,vsize
    REAL, DIMENSION (iip1,hsize,vsize), INTENT(IN) :: field
    REAL, INTENT (IN)                              :: factt

    ! Variables locales
    INTEGER, SAVE :: timestep=0
    ! Identites fichier netcdf
    INTEGER       :: nid, id_lonu, id_lonv, id_latu, id_latv, id_tim, id_lev
    INTEGER       :: vid_lonu,vid_lonv,vid_latu,vid_latv,vid_cu,vid_cv,vid_lev
    INTEGER       :: vid_au,vid_av, varid_alpha_t, varid_alpha_q
    INTEGER       :: l
    INTEGER, DIMENSION (3) :: dim3
    INTEGER, DIMENSION (4) :: dim4,count,start
    INTEGER                :: ierr, varid
    REAL, DIMENSION (iip1,hsize,vsize) :: field2
    
    CALL gather_field(field,iip1*hsize,vsize,0)
    
    IF (mpi_rank /= 0) RETURN
    
    print *,'Guide: output timestep',timestep,'var ',varname
    IF (timestep.EQ.0) THEN 
! ----------------------------------------------
! initialisation fichier de sortie
! ----------------------------------------------
! Ouverture du fichier
        ierr=NF_CREATE("guide_ins.nc",NF_CLOBBER,nid)
! Definition des dimensions
        ierr=NF_DEF_DIM(nid,"LONU",iip1,id_lonu) 
        print*,'id_lonu 1 ',id_lonu
        ierr=NF_DEF_DIM(nid,"LONV",iip1,id_lonv) 
        ierr=NF_DEF_DIM(nid,"LATU",jjp1,id_latu) 
        ierr=NF_DEF_DIM(nid,"LATV",jjm,id_latv) 
        ierr=NF_DEF_DIM(nid,"LEVEL",llm,id_lev)
        ierr=NF_DEF_DIM(nid,"TIME",NF_UNLIMITED,id_tim)

! Creation des variables dimensions
        ierr=NF_DEF_VAR(nid,"LONU",NF_FLOAT,1,id_lonu,vid_lonu)
        print*,'id_lonu 2 ',id_lonu
        ierr=NF_DEF_VAR(nid,"LONV",NF_FLOAT,1,id_lonv,vid_lonv)
        ierr=NF_DEF_VAR(nid,"LATU",NF_FLOAT,1,id_latu,vid_latu)
        ierr=NF_DEF_VAR(nid,"LATV",NF_FLOAT,1,id_latv,vid_latv)
        ierr=NF_DEF_VAR(nid,"LEVEL",NF_FLOAT,1,id_lev,vid_lev)
        ierr=NF_DEF_VAR(nid,"cu",NF_FLOAT,2,(/id_lonu,id_latu/),vid_cu)
        ierr=NF_DEF_VAR(nid,"cv",NF_FLOAT,2,(/id_lonv,id_latv/),vid_cv)
        ierr=NF_DEF_VAR(nid,"au",NF_FLOAT,2,(/id_lonu,id_latu/),vid_au)
        ierr=NF_DEF_VAR(nid,"av",NF_FLOAT,2,(/id_lonv,id_latv/),vid_av)
        call nf95_def_var(nid, "alpha_T", nf90_float, (/id_lonv, id_latu/), &
             varid_alpha_t)
        call nf95_def_var(nid, "alpha_q", nf90_float, (/id_lonv, id_latu/), &
             varid_alpha_q)
        
        ierr=NF_ENDDEF(nid)

! Enregistrement des variables dimensions
#ifdef NC_DOUBLE
        print*,'id_lonu DOUBLE ',id_lonu,rlonu*180./pi
        ierr = NF_PUT_VAR_DOUBLE(nid,vid_lonu,rlonu*180./pi)
        ierr = NF_PUT_VAR_DOUBLE(nid,vid_lonv,rlonv*180./pi)
        ierr = NF_PUT_VAR_DOUBLE(nid,vid_latu,rlatu*180./pi)
        ierr = NF_PUT_VAR_DOUBLE(nid,vid_latv,rlatv*180./pi)
        ierr = NF_PUT_VAR_DOUBLE(nid,vid_lev,presnivs)
        ierr = NF_PUT_VAR_DOUBLE(nid,vid_cu,cu)
        ierr = NF_PUT_VAR_DOUBLE(nid,vid_cv,cv)
        ierr = NF_PUT_VAR_DOUBLE(nid,vid_au,alpha_u)
        ierr = NF_PUT_VAR_DOUBLE(nid,vid_av,alpha_v)
#else
        print*,'id_lonu 3 ',id_lonu,rlonu*180./pi
        ierr = NF_PUT_VAR_REAL(nid,vid_lonu,rlonu*180./pi)
        ierr = NF_PUT_VAR_REAL(nid,vid_lonv,rlonv*180./pi)
        ierr = NF_PUT_VAR_REAL(nid,vid_latu,rlatu*180./pi)
        ierr = NF_PUT_VAR_REAL(nid,vid_latv,rlatv*180./pi)
        ierr = NF_PUT_VAR_REAL(nid,vid_lev,presnivs)
        ierr = NF_PUT_VAR_REAL(nid,vid_cu,cu)
        ierr = NF_PUT_VAR_REAL(nid,vid_cv,cv)
        ierr = NF_PUT_VAR_REAL(nid,vid_au,alpha_u)
        ierr = NF_PUT_VAR_REAL(nid,vid_av,alpha_v)
#endif
        call nf95_put_var(nid, varid_alpha_t, alpha_t)
        call nf95_put_var(nid, varid_alpha_q, alpha_q)
! --------------------------------------------------------------------
! Cr�ation des variables sauvegard�es
! --------------------------------------------------------------------
        ierr = NF_REDEF(nid)
! Pressure (GCM)
        dim4=(/id_lonv,id_latu,id_lev,id_tim/)
        ierr = NF_DEF_VAR(nid,"SP",NF_FLOAT,4,dim4,varid)
! Surface pressure (guidage)
        IF (guide_P) THEN
            dim3=(/id_lonv,id_latu,id_tim/)
            ierr = NF_DEF_VAR(nid,"ps",NF_FLOAT,3,dim3,varid)
        ENDIF
! Zonal wind
        IF (guide_u) THEN
        print*,'id_lonu 4 ',id_lonu,varname
            dim4=(/id_lonu,id_latu,id_lev,id_tim/)
            ierr = NF_DEF_VAR(nid,"u",NF_FLOAT,4,dim4,varid)
            ierr = NF_DEF_VAR(nid,"ua",NF_FLOAT,4,dim4,varid)
            ierr = NF_DEF_VAR(nid,"ucov",NF_FLOAT,4,dim4,varid)
        ENDIF
! Merid. wind
        IF (guide_v) THEN
            dim4=(/id_lonv,id_latv,id_lev,id_tim/)
            ierr = NF_DEF_VAR(nid,"v",NF_FLOAT,4,dim4,varid)
            ierr = NF_DEF_VAR(nid,"va",NF_FLOAT,4,dim4,varid)
            ierr = NF_DEF_VAR(nid,"vcov",NF_FLOAT,4,dim4,varid)
        ENDIF
! Pot. Temperature
        IF (guide_T) THEN
            dim4=(/id_lonv,id_latu,id_lev,id_tim/)
            ierr = NF_DEF_VAR(nid,"teta",NF_FLOAT,4,dim4,varid)
        ENDIF
! Specific Humidity
        IF (guide_Q) THEN
            dim4=(/id_lonv,id_latu,id_lev,id_tim/)
            ierr = NF_DEF_VAR(nid,"q",NF_FLOAT,4,dim4,varid)
        ENDIF
        
        ierr = NF_ENDDEF(nid)
        ierr = NF_CLOSE(nid)
    ENDIF ! timestep=0

! --------------------------------------------------------------------
! Enregistrement du champ
! --------------------------------------------------------------------
 
    ierr=NF_OPEN("guide_ins.nc",NF_WRITE,nid)

    IF (varname=="SP") timestep=timestep+1

    IF (varname=="SP") THEN
      print*,'varname=SP=',varname
    ELSE
      print*,'varname diff SP=',varname
    ENDIF


    ierr = NF_INQ_VARID(nid,varname,varid)
    SELECT CASE (varname)
    CASE ("SP","ps")
        start=(/1,1,1,timestep/)
        count=(/iip1,jjp1,llm,1/)
    CASE ("v","va","vcov")
        start=(/1,1,1,timestep/)
        count=(/iip1,jjm,llm,1/)
    CASE DEFAULT
        start=(/1,1,1,timestep/)
        count=(/iip1,jjp1,llm,1/)
    END SELECT

    SELECT CASE (varname)
    CASE("u","ua")
        DO l=1,llm ; field2(:,2:jjm,l)=field(:,2:jjm,l)/cu(:,2:jjm) ; ENDDO
        field2(:,1,:)=0. ; field2(:,jjp1,:)=0.
    CASE("v","va")
        DO l=1,llm ; field2(:,:,l)=field(:,:,l)/cv(:,:) ; ENDDO
    CASE DEFAULT
        field2=field
    END SELECT

#ifdef NC_DOUBLE
    ierr = NF_PUT_VARA_DOUBLE(nid,varid,start,count,field2)
#else
    ierr = NF_PUT_VARA_REAL(nid,varid,start,count,field2)
#endif
 
    ierr = NF_CLOSE(nid)

  END SUBROUTINE guide_out
    
  
!===========================================================================
  subroutine correctbid(iim,nl,x)
    integer iim,nl
    real x(iim+1,nl)
    integer i,l
    real zz

    do l=1,nl
        do i=2,iim-1
            if(abs(x(i,l)).gt.1.e10) then
               zz=0.5*(x(i-1,l)+x(i+1,l))
              print*,'correction ',i,l,x(i,l),zz
               x(i,l)=zz
            endif
         enddo
     enddo
     return
  end subroutine correctbid

!===========================================================================
END MODULE guide_p_mod
