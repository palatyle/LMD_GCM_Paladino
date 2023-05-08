!
MODULE surf_land_orchidee_mod
#ifndef ORCHIDEE_NOOPENMP
!
! This module controles the interface towards the model ORCHIDEE
!
! Subroutines in this module : surf_land_orchidee
!                              Init_orchidee_index
!                              Get_orchidee_communicator
!                              Init_neighbours

  USE dimphy
#ifdef CPP_VEGET
  USE intersurf     ! module d'ORCHIDEE
#endif
  USE cpl_mod,      ONLY : cpl_send_land_fields
  USE surface_data, ONLY : type_ocean
  USE comgeomphy,   ONLY : cuphy, cvphy
  USE mod_grid_phy_lmdz
  USE mod_phys_lmdz_para, mpi_root_rank=>mpi_root

  IMPLICIT NONE

  PRIVATE
  PUBLIC  :: surf_land_orchidee

  LOGICAL, ALLOCATABLE, SAVE :: flag_omp(:)
CONTAINS
!
!****************************************************************************************
!  
  SUBROUTINE surf_land_orchidee(itime, dtime, date0, knon, &
       knindex, rlon, rlat, pctsrf, &
       debut, lafin, &
       plev,  u1_lay, v1_lay, temp_air, spechum, epot_air, ccanopy, & 
       tq_cdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
       precip_rain, precip_snow, lwdown, swnet, swdown, &
       ps, q2m, t2m, &
       evap, fluxsens, fluxlat, &              
       tsol_rad, tsurf_new, alb1_new, alb2_new, &
       emis_new, z0_new, qsurf)

    USE mod_surf_para
    USE mod_synchro_omp
    USE carbon_cycle_mod, ONLY : carbon_cycle_cpl

!    
! Cette routine sert d'interface entre le modele atmospherique et le 
! modele de sol continental. Appel a sechiba
!
! L. Fairhead 02/2000
!
! input:
!   itime        numero du pas de temps
!   dtime        pas de temps de la physique (en s)
!   nisurf       index de la surface a traiter (1 = sol continental)
!   knon         nombre de points de la surface a traiter
!   knindex      index des points de la surface a traiter
!   rlon         longitudes de la grille entiere
!   rlat         latitudes de la grille entiere
!   pctsrf       tableau des fractions de surface de chaque maille
!   debut        logical: 1er appel a la physique (lire les restart)
!   lafin        logical: dernier appel a la physique (ecrire les restart)
!                     (si false calcul simplifie des fluxs sur les continents)
!   plev         hauteur de la premiere couche (Pa)      
!   u1_lay       vitesse u 1ere couche
!   v1_lay       vitesse v 1ere couche
!   temp_air     temperature de l'air 1ere couche
!   spechum      humidite specifique 1ere couche
!   epot_air     temp pot de l'air
!   ccanopy      concentration CO2 canopee, correspond au co2_send de 
!                carbon_cycle_mod ou valeur constant co2_ppm
!   tq_cdrag     cdrag
!   petAcoef     coeff. A de la resolution de la CL pour t
!   peqAcoef     coeff. A de la resolution de la CL pour q
!   petBcoef     coeff. B de la resolution de la CL pour t
!   peqBcoef     coeff. B de la resolution de la CL pour q
!   precip_rain  precipitation liquide
!   precip_snow  precipitation solide
!   lwdown       flux IR descendant a la surface
!   swnet        flux solaire net
!   swdown       flux solaire entrant a la surface
!   ps           pression au sol
!   radsol       rayonnement net aus sol (LW + SW)
!   
!
! output:
!   evap         evaporation totale
!   fluxsens     flux de chaleur sensible
!   fluxlat      flux de chaleur latente
!   tsol_rad     
!   tsurf_new    temperature au sol
!   alb1_new     albedo in visible SW interval
!   alb2_new     albedo in near IR interval
!   emis_new     emissivite
!   z0_new       surface roughness
!   qsurf        air moisture at surface
!
    INCLUDE "indicesol.h"
    INCLUDE "temps.h"
    INCLUDE "YOMCST.h"
    INCLUDE "iniprint.h"
    INCLUDE "dimensions.h"
  
!
! Parametres d'entree
!****************************************************************************************
    INTEGER, INTENT(IN)                       :: itime
    REAL, INTENT(IN)                          :: dtime
    REAL, INTENT(IN)                          :: date0
    INTEGER, INTENT(IN)                       :: knon
    INTEGER, DIMENSION(klon), INTENT(IN)      :: knindex
    LOGICAL, INTENT(IN)                       :: debut, lafin
    REAL, DIMENSION(klon,nbsrf), INTENT(IN)   :: pctsrf
    REAL, DIMENSION(klon), INTENT(IN)         :: rlon, rlat
    REAL, DIMENSION(klon), INTENT(IN)         :: plev
    REAL, DIMENSION(klon), INTENT(IN)         :: u1_lay, v1_lay
    REAL, DIMENSION(klon), INTENT(IN)         :: temp_air, spechum
    REAL, DIMENSION(klon), INTENT(IN)         :: epot_air, ccanopy
    REAL, DIMENSION(klon), INTENT(IN)         :: tq_cdrag
    REAL, DIMENSION(klon), INTENT(IN)         :: petAcoef, peqAcoef
    REAL, DIMENSION(klon), INTENT(IN)         :: petBcoef, peqBcoef
    REAL, DIMENSION(klon), INTENT(IN)         :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)         :: lwdown, swnet, swdown, ps
    REAL, DIMENSION(klon), INTENT(IN)         :: q2m, t2m

! Parametres de sortie
!****************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)        :: evap, fluxsens, fluxlat, qsurf
    REAL, DIMENSION(klon), INTENT(OUT)        :: tsol_rad, tsurf_new
    REAL, DIMENSION(klon), INTENT(OUT)        :: alb1_new, alb2_new
    REAL, DIMENSION(klon), INTENT(OUT)        :: emis_new, z0_new

! Local
!****************************************************************************************
    INTEGER                                   :: ij, jj, igrid, ireal, index
    INTEGER                                   :: error
    REAL, DIMENSION(klon)                     :: swdown_vrai
    CHARACTER (len = 20)                      :: modname = 'surf_land_orchidee'
    CHARACTER (len = 80)                      :: abort_message
    LOGICAL,SAVE                              :: check = .FALSE.
    !$OMP THREADPRIVATE(check)

! type de couplage dans sechiba
!  character (len=10)   :: coupling = 'implicit' 
! drapeaux controlant les appels dans SECHIBA
!  type(control_type), save   :: control_in
! Preserved albedo
    REAL, ALLOCATABLE, DIMENSION(:), SAVE     :: albedo_keep, zlev
    !$OMP THREADPRIVATE(albedo_keep,zlev)
! coordonnees geographiques
    REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: lalo
    !$OMP THREADPRIVATE(lalo)
! pts voisins
    INTEGER,ALLOCATABLE, DIMENSION(:,:), SAVE :: neighbours
    !$OMP THREADPRIVATE(neighbours)
! fractions continents
    REAL,ALLOCATABLE, DIMENSION(:), SAVE      :: contfrac
    !$OMP THREADPRIVATE(contfrac)
! resolution de la grille
    REAL, ALLOCATABLE, DIMENSION (:,:), SAVE  :: resolution
    !$OMP THREADPRIVATE(resolution)

    REAL, ALLOCATABLE, DIMENSION (:,:), SAVE  :: lon_scat, lat_scat  
    !$OMP THREADPRIVATE(lon_scat,lat_scat)

    LOGICAL, SAVE                             :: lrestart_read = .TRUE.
    !$OMP THREADPRIVATE(lrestart_read)
    LOGICAL, SAVE                             :: lrestart_write = .FALSE.
    !$OMP THREADPRIVATE(lrestart_write)

    REAL, DIMENSION(knon,2)                   :: albedo_out

! Pb de nomenclature
    REAL, DIMENSION(klon)                     :: petA_orc, peqA_orc
    REAL, DIMENSION(klon)                     :: petB_orc, peqB_orc
! Pb de correspondances de grilles
    INTEGER, DIMENSION(:), SAVE, ALLOCATABLE  :: ig, jg
    !$OMP THREADPRIVATE(ig,jg)
    INTEGER :: indi, indj
    INTEGER, SAVE, ALLOCATABLE,DIMENSION(:)   :: ktindex
    !$OMP THREADPRIVATE(ktindex)

! Essai cdrag
    REAL, DIMENSION(klon)                     :: cdrag
    INTEGER,SAVE                              :: offset
    !$OMP THREADPRIVATE(offset)

    REAL, DIMENSION(klon_glo)                 :: rlon_g,rlat_g
    INTEGER, SAVE                             :: orch_comm
    !$OMP THREADPRIVATE(orch_comm)

    REAL, ALLOCATABLE, DIMENSION(:), SAVE     :: coastalflow
    !$OMP THREADPRIVATE(coastalflow)
    REAL, ALLOCATABLE, DIMENSION(:), SAVE     :: riverflow
    !$OMP THREADPRIVATE(riverflow)
    
    INTEGER :: orch_omp_rank
    INTEGER :: orch_omp_size
!
! Fin definition
!****************************************************************************************

    IF (check) WRITE(lunout,*)'Entree ', modname
  
! Initialisation
  
    IF (debut) THEN
! Test of coherence between variable ok_veget and cpp key CPP_VEGET
#ifndef CPP_VEGET
       abort_message='Pb de coherence: ok_veget = .true. mais CPP_VEGET = .false.'
       CALL abort_gcm(modname,abort_message,1)
#endif

       CALL Init_surf_para(knon)
       ALLOCATE(ktindex(knon))
       IF ( .NOT. ALLOCATED(albedo_keep)) THEN
!ym          ALLOCATE(albedo_keep(klon))
!ym bizarre que non alloué en knon precedement
          ALLOCATE(albedo_keep(knon))
          ALLOCATE(zlev(knon))
       ENDIF
! Pb de correspondances de grilles
       ALLOCATE(ig(klon))
       ALLOCATE(jg(klon))
       ig(1) = 1
       jg(1) = 1
       indi = 0
       indj = 2
       DO igrid = 2, klon - 1
          indi = indi + 1
          IF ( indi > iim) THEN
             indi = 1
             indj = indj + 1
          ENDIF
          ig(igrid) = indi
          jg(igrid) = indj
       ENDDO
       ig(klon) = 1
       jg(klon) = jjm + 1

       IF ((.NOT. ALLOCATED(lalo))) THEN
          ALLOCATE(lalo(knon,2), stat = error)
          IF (error /= 0) THEN
             abort_message='Pb allocation lalo'
             CALL abort_gcm(modname,abort_message,1)
          ENDIF
       ENDIF
       IF ((.NOT. ALLOCATED(lon_scat))) THEN
          ALLOCATE(lon_scat(iim,jjm+1), stat = error)
          IF (error /= 0) THEN
             abort_message='Pb allocation lon_scat'
             CALL abort_gcm(modname,abort_message,1)
          ENDIF
       ENDIF
       IF ((.NOT. ALLOCATED(lat_scat))) THEN
          ALLOCATE(lat_scat(iim,jjm+1), stat = error)
          IF (error /= 0) THEN
             abort_message='Pb allocation lat_scat'
             CALL abort_gcm(modname,abort_message,1)
          ENDIF
       ENDIF
       lon_scat = 0.
       lat_scat = 0.
       DO igrid = 1, knon
          index = knindex(igrid)
          lalo(igrid,2) = rlon(index)
          lalo(igrid,1) = rlat(index)
       ENDDO

       
       
       CALL Gather(rlon,rlon_g)
       CALL Gather(rlat,rlat_g)

       IF (is_mpi_root) THEN
          index = 1
          DO jj = 2, jjm
             DO ij = 1, iim
                index = index + 1
                lon_scat(ij,jj) = rlon_g(index)
                lat_scat(ij,jj) = rlat_g(index)
             ENDDO
          ENDDO
          lon_scat(:,1) = lon_scat(:,2)
          lat_scat(:,1) = rlat_g(1)
          lon_scat(:,jjm+1) = lon_scat(:,2)
          lat_scat(:,jjm+1) = rlat_g(klon_glo)
       ENDIF
   
       CALL bcast(lon_scat)
       CALL bcast(lat_scat)
!
! Allouer et initialiser le tableau des voisins et des fraction de continents
!
       IF ( (.NOT.ALLOCATED(neighbours))) THEN
          ALLOCATE(neighbours(knon,8), stat = error)
          IF (error /= 0) THEN
             abort_message='Pb allocation neighbours'
             CALL abort_gcm(modname,abort_message,1)
          ENDIF
       ENDIF
       neighbours = -1.
       IF (( .NOT. ALLOCATED(contfrac))) THEN
          ALLOCATE(contfrac(knon), stat = error)
          IF (error /= 0) THEN
             abort_message='Pb allocation contfrac'
             CALL abort_gcm(modname,abort_message,1)
          ENDIF
       ENDIF

       DO igrid = 1, knon
          ireal = knindex(igrid)
          contfrac(igrid) = pctsrf(ireal,is_ter)
       ENDDO


       CALL Init_neighbours(knon,neighbours,knindex,pctsrf(:,is_ter))

!
!  Allocation et calcul resolutions
       IF ( (.NOT.ALLOCATED(resolution))) THEN
          ALLOCATE(resolution(knon,2), stat = error)
          IF (error /= 0) THEN
             abort_message='Pb allocation resolution'
             CALL abort_gcm(modname,abort_message,1)
          ENDIF
       ENDIF
       DO igrid = 1, knon
          ij = knindex(igrid)
          resolution(igrid,1) = cuphy(ij)
          resolution(igrid,2) = cvphy(ij)
       ENDDO
     
       ALLOCATE(coastalflow(klon), stat = error)
       IF (error /= 0) THEN
          abort_message='Pb allocation coastalflow'
          CALL abort_gcm(modname,abort_message,1)
       ENDIF
       
       ALLOCATE(riverflow(klon), stat = error)
       IF (error /= 0) THEN
          abort_message='Pb allocation riverflow'
          CALL abort_gcm(modname,abort_message,1)
       ENDIF
!
! carbon_cycle_cpl not possible with this interface and version of ORHCHIDEE
!
       IF (carbon_cycle_cpl) THEN
          abort_message='carbon_cycle_cpl not yet possible with this interface of ORCHIDEE'
          CALL abort_gcm(modname,abort_message,1)
       END IF
       
    ENDIF                          ! (fin debut) 
 

! 
! Appel a la routine sols continentaux
!
    IF (lafin) lrestart_write = .TRUE.
    IF (check) WRITE(lunout,*)'lafin ',lafin,lrestart_write
     
    petA_orc(1:knon) = petBcoef(1:knon) * dtime
    petB_orc(1:knon) = petAcoef(1:knon)
    peqA_orc(1:knon) = peqBcoef(1:knon) * dtime
    peqB_orc(1:knon) = peqAcoef(1:knon)

    cdrag = 0.
    cdrag(1:knon) = tq_cdrag(1:knon)

! zlev(1:knon) = (100.*plev(1:knon))/((ps(1:knon)/287.05*temp_air(1:knon))*9.80665)
    zlev(1:knon) = (100.*plev(1:knon))/((ps(1:knon)/RD*temp_air(1:knon))*RG)


! PF et PASB
!   where(cdrag > 0.01) 
!     cdrag = 0.01
!   endwhere
!  write(*,*)'Cdrag = ',minval(cdrag),maxval(cdrag)

  
    IF (debut) THEN
       CALL Init_orchidee_index(knon,knindex,offset,ktindex)
       CALL Get_orchidee_communicator(orch_comm,orch_omp_size,orch_omp_rank)
       CALL Init_synchro_omp
       
       IF (knon > 0) THEN
#ifdef CPP_VEGET
         CALL Init_intersurf(nbp_lon,nbp_lat,knon,ktindex,offset,orch_omp_size,orch_omp_rank,orch_comm)
#endif
       ENDIF

       
       IF (knon > 0) THEN

#ifdef CPP_VEGET
          CALL intersurf_main (itime+itau_phy-1, iim, jjm+1, knon, ktindex, dtime, &
               lrestart_read, lrestart_write, lalo, &
               contfrac, neighbours, resolution, date0, &
               zlev,  u1_lay, v1_lay, spechum, temp_air, epot_air, ccanopy, &
               cdrag, petA_orc, peqA_orc, petB_orc, peqB_orc, &
               precip_rain, precip_snow, lwdown, swnet, swdown, ps, &
               evap, fluxsens, fluxlat, coastalflow, riverflow, &
               tsol_rad, tsurf_new, qsurf, albedo_out, emis_new, z0_new, &
               lon_scat, lat_scat, q2m, t2m)
#endif         
       ENDIF

       CALL Synchro_omp

       albedo_keep(1:knon) = (albedo_out(1:knon,1)+albedo_out(1:knon,2))/2.

    ENDIF

    
!  swdown_vrai(1:knon) = swnet(1:knon)/(1. - albedo_keep(1:knon))
    swdown_vrai(1:knon) = swdown(1:knon)

    IF (knon > 0) THEN
#ifdef CPP_VEGET    
       CALL intersurf_main (itime+itau_phy, iim, jjm+1, knon, ktindex, dtime,  &
            lrestart_read, lrestart_write, lalo, &
            contfrac, neighbours, resolution, date0, &
            zlev,  u1_lay(1:knon), v1_lay(1:knon), spechum(1:knon), temp_air(1:knon), epot_air(1:knon), ccanopy(1:knon), &
            cdrag(1:knon), petA_orc(1:knon), peqA_orc(1:knon), petB_orc(1:knon), peqB_orc(1:knon), &
            precip_rain(1:knon), precip_snow(1:knon), lwdown(1:knon), swnet(1:knon), swdown_vrai(1:knon), ps(1:knon), &
            evap(1:knon), fluxsens(1:knon), fluxlat(1:knon), coastalflow(1:knon), riverflow(1:knon), &
            tsol_rad(1:knon), tsurf_new(1:knon), qsurf(1:knon), albedo_out(1:knon,:), emis_new(1:knon), z0_new(1:knon), &
            lon_scat, lat_scat, q2m, t2m)
#endif       
    ENDIF

    CALL Synchro_omp
    
    albedo_keep(1:knon) = (albedo_out(1:knon,1)+albedo_out(1:knon,2))/2.

!* Send to coupler
!
    IF (type_ocean=='couple') THEN
       CALL cpl_send_land_fields(itime, knon, knindex, &
            riverflow, coastalflow)
    ENDIF

    alb1_new(1:knon) = albedo_out(1:knon,1) 
    alb2_new(1:knon) = albedo_out(1:knon,2)

! Convention orchidee: positif vers le haut
    fluxsens(1:knon) = -1. * fluxsens(1:knon)
    fluxlat(1:knon)  = -1. * fluxlat(1:knon)
    
!  evap     = -1. * evap

    IF (debut) lrestart_read = .FALSE.
    
    IF (debut) CALL Finalize_surf_para

    
  END SUBROUTINE surf_land_orchidee
!
!****************************************************************************************
!
  SUBROUTINE Init_orchidee_index(knon,knindex,offset,ktindex)
  USE mod_surf_para
  USE mod_grid_phy_lmdz
  
    INTEGER,INTENT(IN)    :: knon
    INTEGER,INTENT(IN)    :: knindex(klon)    
    INTEGER,INTENT(OUT)   :: offset
    INTEGER,INTENT(OUT)   :: ktindex(klon)
    
    INTEGER               :: ktindex_glo(knon_glo)
    INTEGER               :: offset_para(0:omp_size*mpi_size-1)
    INTEGER               :: LastPoint
    INTEGER               :: task
    
    ktindex(1:knon)=knindex(1:knon)+(klon_mpi_begin-1)+(klon_omp_begin-1)+nbp_lon-1
    
    CALL gather_surf(ktindex(1:knon),ktindex_glo) 
    
    IF (is_mpi_root .AND. is_omp_root) THEN
      LastPoint=0
      DO Task=0,mpi_size*omp_size-1
        IF (knon_glo_para(Task)>0) THEN
           offset_para(task)= LastPoint-MOD(LastPoint,nbp_lon)
           LastPoint=ktindex_glo(knon_glo_end_para(task))
        ENDIF
      ENDDO
    ENDIF
    
    CALL bcast(offset_para)
    
    offset=offset_para(omp_size*mpi_rank+omp_rank)
    
    ktindex(1:knon)=ktindex(1:knon)-offset

  END SUBROUTINE Init_orchidee_index

!
!************************* ***************************************************************
! 

  SUBROUTINE Get_orchidee_communicator(orch_comm,orch_omp_size,orch_omp_rank)
  USE  mod_surf_para
      
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
#endif    

    INTEGER,INTENT(OUT) :: orch_comm
    INTEGER,INTENT(OUT) :: orch_omp_size
    INTEGER,INTENT(OUT) :: orch_omp_rank
    INTEGER             :: color
    INTEGER             :: i,ierr
!
! End definition
!****************************************************************************************
    
    
    IF (is_omp_root) THEN          
      
      IF (knon_mpi==0) THEN 
         color = 0
      ELSE 
         color = 1
      ENDIF
    
#ifdef CPP_MPI    
      CALL MPI_COMM_SPLIT(COMM_LMDZ_PHY,color,mpi_rank,orch_comm,ierr)
#endif
    
    ENDIF
    CALL bcast_omp(orch_comm)
    
    IF (knon_mpi /= 0) THEN
      orch_omp_size=0
      DO i=0,omp_size-1
        IF (knon_omp_para(i) /=0) THEN
          orch_omp_size=orch_omp_size+1
          IF (i==omp_rank) orch_omp_rank=orch_omp_size-1
        ENDIF
      ENDDO
    ENDIF
   
    
  END SUBROUTINE Get_orchidee_communicator
!
!****************************************************************************************
!  

  SUBROUTINE Init_neighbours(knon,neighbours,knindex,pctsrf)
    USE mod_grid_phy_lmdz
    USE mod_surf_para    
    INCLUDE "indicesol.h"

#ifdef CPP_MPI
    INCLUDE 'mpif.h'
#endif    

! Input arguments
!****************************************************************************************
    INTEGER, INTENT(IN)                     :: knon
    INTEGER, DIMENSION(klon), INTENT(IN)    :: knindex
    REAL, DIMENSION(klon), INTENT(IN)       :: pctsrf
    
! Output arguments
!****************************************************************************************
    INTEGER, DIMENSION(knon,8), INTENT(OUT) :: neighbours

! Local variables
!****************************************************************************************
    INTEGER                              :: i, igrid, jj, ij, iglob
    INTEGER                              :: ierr, ireal, index
    INTEGER, DIMENSION(8,3)              :: off_ini
    INTEGER, DIMENSION(8)                :: offset  
    INTEGER, DIMENSION(nbp_lon,nbp_lat)  :: correspond
    INTEGER, DIMENSION(knon_glo)         :: ktindex_glo
    INTEGER, DIMENSION(knon_glo,8)       :: neighbours_glo
    REAL, DIMENSION(klon_glo)            :: pctsrf_glo
    INTEGER                              :: ktindex(klon)
!
! End definition
!****************************************************************************************

    ktindex(1:knon)=knindex(1:knon)+(klon_mpi_begin-1)+(klon_omp_begin-1)+nbp_lon-1
    
    CALL gather_surf(ktindex(1:knon),ktindex_glo)
    CALL gather(pctsrf,pctsrf_glo)
    
    IF (is_mpi_root .AND. is_omp_root) THEN
      neighbours_glo(:,:)=-1
!  Initialisation des offset    
!
! offset bord ouest
       off_ini(1,1) = - nbp_lon   ; off_ini(2,1) = - nbp_lon + 1     ; off_ini(3,1) = 1
       off_ini(4,1) = nbp_lon + 1 ; off_ini(5,1) = nbp_lon           ; off_ini(6,1) = 2 * nbp_lon - 1
       off_ini(7,1) = nbp_lon -1  ; off_ini(8,1) = - 1 
! offset point normal
       off_ini(1,2) = - nbp_lon   ; off_ini(2,2) = - nbp_lon + 1     ; off_ini(3,2) = 1
       off_ini(4,2) = nbp_lon + 1 ; off_ini(5,2) = nbp_lon           ; off_ini(6,2) = nbp_lon - 1
       off_ini(7,2) = -1          ; off_ini(8,2) = - nbp_lon - 1
! offset bord   est
       off_ini(1,3) = - nbp_lon   ; off_ini(2,3) = - 2 * nbp_lon + 1 ; off_ini(3,3) = - nbp_lon + 1
       off_ini(4,3) =  1          ; off_ini(5,3) = nbp_lon           ; off_ini(6,3) = nbp_lon - 1
       off_ini(7,3) = -1          ; off_ini(8,3) = - nbp_lon - 1
!
!
! Attention aux poles
!
       DO igrid = 1, knon_glo
          index = ktindex_glo(igrid)
          jj = INT((index - 1)/nbp_lon) + 1
          ij = index - (jj - 1) * nbp_lon
          correspond(ij,jj) = igrid
       ENDDO
       
       DO igrid = 1, knon_glo
          iglob = ktindex_glo(igrid)
          
          IF (MOD(iglob, nbp_lon) == 1) THEN
             offset = off_ini(:,1)
          ELSE IF(MOD(iglob, nbp_lon) == 0) THEN
             offset = off_ini(:,3)
          ELSE
             offset = off_ini(:,2)
          ENDIF
          
          DO i = 1, 8
             index = iglob + offset(i)
             ireal = (MIN(MAX(1, index - nbp_lon + 1), klon_glo))
             IF (pctsrf_glo(ireal) > EPSFRA) THEN
                jj = INT((index - 1)/nbp_lon) + 1
                ij = index - (jj - 1) * nbp_lon
                neighbours_glo(igrid, i) = correspond(ij, jj)
             ENDIF
          ENDDO
       ENDDO

    ENDIF
    
    DO i = 1, 8
      CALL scatter_surf(neighbours_glo(:,i),neighbours(1:knon,i))
    ENDDO
  END SUBROUTINE Init_neighbours

!
!****************************************************************************************
!
#endif
END MODULE surf_land_orchidee_mod
