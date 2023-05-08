       MODULE watercloud_mod

       IMPLICIT NONE

       REAL,SAVE,ALLOCATABLE :: zdqcloud(:,:,:) ! tendencies on pq due to condensation of H2O(kg/kg.s-1)
       REAL,SAVE,ALLOCATABLE :: zdqscloud(:,:) ! tendencies on qsurf (calculated only by calchim but declared here)

       CONTAINS

       SUBROUTINE watercloud(ngrid,nlay,ptimestep, 
     &                pplev,pplay,pdpsrf,pzlay,pt,pdt,
     &                pq,pdq,pdqcloud,pdtcloud,
     &                nq,tau,tauscaling,rdust,rice,nuice,
     &                rsedcloud,rhocloud,totcloudfrac)
      USE ioipsl_getincom, ONLY: getin
      USE updaterad, ONLY: updaterdust, updaterice_micro,
     &                     updaterice_typ
      USE improvedclouds_mod, ONLY: improvedclouds
      USE watersat_mod, ONLY: watersat
      use tracer_mod, only: nqmx, igcm_h2o_vap, igcm_h2o_ice,
     &                      igcm_dust_mass, igcm_dust_number,
     &                      igcm_ccn_mass, igcm_ccn_number,
     &                      rho_dust, nuice_sed, nuice_ref
      use dimradmars_mod, only: naerkind
      IMPLICIT NONE


c=======================================================================
c  Water-ice cloud formation
c  
c  Includes two different schemes:
c    - A simplified scheme (see simpleclouds.F)
c    - An improved microphysical scheme (see improvedclouds.F)
c
c  There is a time loop specific to cloud formation 
c  due to timescales smaller than the GCM integration timestep.
c
c  Authors: Franck Montmessin, Francois Forget, Ehouarn Millour, 
c           J.-B. Madeleine, Thomas Navarro
c
c  2004 - 2012
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

      include "callkeys.h"

c   Inputs/outputs:
c   ------

      INTEGER, INTENT(IN) :: ngrid,nlay
      INTEGER, INTENT(IN) ::  nq                 ! nombre de traceurs 
      REAL, INTENT(IN) ::  ptimestep             ! pas de temps physique (s)
      REAL, INTENT(IN) ::  pplev(ngrid,nlay+1)   ! pression aux inter-couches (Pa)
      REAL, INTENT(IN) ::  pplay(ngrid,nlay)     ! pression au milieu des couches (Pa)
      REAL, INTENT(IN) ::  pdpsrf(ngrid)         ! tendence surf pressure
      REAL, INTENT(IN) ::  pzlay(ngrid,nlay)     ! altitude at the middle of the layers
      REAL, INTENT(IN) ::  pt(ngrid,nlay)        ! temperature at the middle of the layers (K)
      REAL, INTENT(IN) ::  pdt(ngrid,nlay)       ! tendence temperature des autres param.

      REAL, INTENT(IN) ::  pq(ngrid,nlay,nq)     ! traceur (kg/kg)
      rEAL, INTENT(IN) ::  pdq(ngrid,nlay,nq)    ! tendence avant condensation  (kg/kg.s-1)

      REAL, INTENT(IN) ::  tau(ngrid,naerkind) ! Column dust optical depth at each point
      REAL, INTENT(IN) ::  tauscaling(ngrid)   ! Convertion factor for dust amount
      REAL, INTENT(INOUT) ::  rdust(ngrid,nlay)   ! Dust geometric mean radius (m)

      REAL, INTENT(OUT) ::  pdqcloud(ngrid,nlay,nq) ! tendence de la condensation H2O(kg/kg.s-1)
      REAL, INTENT(OUT) ::  pdtcloud(ngrid,nlay)    ! tendence temperature due
                                   ! a la chaleur latente
      REAL, INTENT(INOUT) ::  rice(ngrid,nlay)    ! Ice mass mean radius (m)
                               ! (r_c in montmessin_2004)
      REAL, INTENT(OUT) ::  nuice(ngrid,nlay)   ! Estimated effective variance
                               !   of the size distribution
      REAL, INTENT(OUT) ::  rsedcloud(ngrid,nlay) ! Cloud sedimentation radius
      REAL, INTENT(OUT) ::  rhocloud(ngrid,nlay)  ! Cloud density (kg.m-3)

      REAL, INTENT(INOUT):: totcloudfrac(ngrid) ! Cloud fraction (A. Pottier 2013)
      
c   Locals:
c   ------
  
      ! for ice radius computation
      REAL Mo,No
      REAl ccntyp
      
      ! for time loop
      INTEGER microstep  ! time subsampling step variable
      INTEGER,SAVE :: imicro ! time subsampling for coupled water microphysics & sedimentation
      REAL,SAVE :: microtimestep ! integration timestep for coupled water microphysics & sedimentation
      REAL,SAVE :: microtimestep_prev=-999
      
      ! tendency given by clouds (inside the micro loop)
      REAL subpdqcloud(ngrid,nlay,nq) ! cf. pdqcloud
      REAL subpdtcloud(ngrid,nlay)    ! cf. pdtcloud

      ! global tendency (clouds+physics)
      REAL sum_subpdq(ngrid,nlay,nq)      ! cf. pdqcloud
      REAL sum_subpdt(ngrid,nlay)         ! cf. pdtcloud

      ! no supersaturation when option supersat is false
      REAL zt(ngrid,nlay)       ! local value of temperature
      REAL zqsat(ngrid,nlay)    ! saturation

      INTEGER iq,ig,l
      LOGICAL,SAVE :: firstcall=.true.

! Representation of sub-grid water ice clouds A. Pottier 2013
      REAL :: ztclf(ngrid, nlay)
      REAL :: zqclf(ngrid, nlay,nq)
      REAL :: zdelt  
      REAL :: norm
      REAL :: ponder
      REAL :: tcond(ngrid,nlay)
      REAL :: zqvap(ngrid,nlay)
      REAL :: zqice(ngrid,nlay)
      REAL ::  spant ! delta T for the temperature distribution
!      REAL :: zqsat(ngrid,nlay) ! saturation
      REAL :: pteff(ngrid, nlay)! effective temperature in the cloud,neb
      REAL :: pqeff(ngrid, nlay, nq)! effective tracers quantities in the cloud
      REAL :: cloudfrac(ngrid,nlay) ! cloud fraction
      REAL :: mincloud ! min cloud frac
      LOGICAL, save :: flagcloud=.true. 
c    ** un petit test de coherence
c       --------------------------

      IF (firstcall) THEN
         
        if (nq.gt.nqmx) then
           write(*,*) 'stop in watercloud (nq.gt.nqmx)!'
           write(*,*) 'nq=',nq,' nqmx=',nqmx
           stop
        endif
         
        write(*,*) "watercloud: igcm_h2o_vap=",igcm_h2o_vap
        write(*,*) "            igcm_h2o_ice=",igcm_h2o_ice
                
        write(*,*) "time subsampling for microphysic ?"
#ifdef MESOSCALE
        imicro = 2
#else
        imicro = 30
#endif
        call getin("imicro",imicro)
        write(*,*)"watercloud: imicro = ",imicro
        
        firstcall=.false.
      ENDIF ! of IF (firstcall)

      !! AS: moved out of firstcall to allow nesting+evoluting timestep
      !!     TBD: consider possible diff imicro with domains?
      microtimestep = ptimestep/real(imicro)
      if (microtimestep/=microtimestep_prev) then
        ! only tell the world if microtimestep has changed
        write(*,*)"watercloud: Physical timestep is ",ptimestep
        write(*,*)"watercloud: Microphysics timestep is ",microtimestep
        microtimestep_prev=microtimestep
      endif
      
c-----Initialization
      sum_subpdq(1:ngrid,1:nlay,1:nq) = 0
      sum_subpdt(1:ngrid,1:nlay)      = 0
      
      ! default value if no ice
      rhocloud(1:ngrid,1:nlay) = rho_dust

c-------------------------------------------------------------------
c   0.  Representation of sub-grid water ice clouds
c------------------
c-----Initialization
      pteff(1:ngrid,1:nlay) = 0
      pqeff(1:ngrid,1:nlay,1:nq) = 0
      DO l=1,nlay
        DO ig=1,ngrid
             pteff(ig,l)=pt(ig,l)
        END DO
      END DO
      DO l=1,nlay
        DO ig=1,ngrid
          DO iq=1,nq
             pqeff(ig,l,iq)=pq(ig,l,iq)
          ENDDO
        ENDDO
      ENDDO
c-----Tendencies
      DO l=1,nlay
        DO ig=1,ngrid
          ztclf(ig,l)=pt(ig,l)+ pdt(ig,l)*ptimestep
        ENDDO
      ENDDO
      DO l=1,nlay
        DO ig=1,ngrid
          DO iq=1,nq
             zqclf(ig,l,iq)=pq(ig,l,iq)+pdq(ig,l,iq)*ptimestep 
          ENDDO
        ENDDO
      ENDDO
c-----Effective temperature calculation
      IF (CLFvarying) THEN
         spant=3.0 ! delta T for the temprature distribution
         mincloud=0.1 ! min cloudfrac when there is ice
         IF (flagcloud) THEN
             write(*,*) "Delta T", spant
             write(*,*) "mincloud", mincloud
             flagcloud=.false.
         END IF
         !CALL watersat(ngrid*nlay,ztclf,pplay,zqsat) !MV17: we dont need zqsat in the CLFvarying scheme
         zqvap=zqclf(:,:,igcm_h2o_vap)
         zqice=zqclf(:,:,igcm_h2o_ice)
         CALL tcondwater(ngrid*nlay,pplay,zqvap+zqice,tcond)
         DO l=1,nlay
           DO ig=1,ngrid
              zdelt=spant !MAX(spant*ztclf(ig,l),1.e-12), now totally in K. Fixed width
              IF (tcond(ig,l) .ge. (ztclf(ig,l)+zdelt)) THEN
                 pteff(ig,l)=ztclf(ig,l)
                 cloudfrac(ig,l)=1.
              ELSE IF (tcond(ig,l) .le. (ztclf(ig,l)-zdelt)) THEN
                 pteff(ig,l)=ztclf(ig,l)-zdelt
                 cloudfrac(ig,l)=mincloud
              ELSE
                 cloudfrac(ig,l)=(tcond(ig,l)-ztclf(ig,l)+zdelt)/
     &                           (2.0*zdelt)
                 pteff(ig,l)=(tcond(ig,l)+ztclf(ig,l)-zdelt)/2.
              END IF
              pteff(ig,l)=pteff(ig,l)-pdt(ig,l)*ptimestep
              IF (cloudfrac(ig,l).le.mincloud) THEN !MV17: replaced .le.0 by .le.mincloud
                 cloudfrac(ig,l)=mincloud
              ELSE IF (cloudfrac(ig,l).gt.1) THEN
                 cloudfrac(ig,l)=1.
              END IF
           ENDDO
         ENDDO
c-----Calculation of the total cloud coverage of the column
         DO ig=1,ngrid
            totcloudfrac(ig) = 0.
            norm=0.
            DO l=1,nlay
               ponder=zqice(ig,l)*(pplev(ig,l) - pplev(ig,l+1))
               totcloudfrac(ig) = totcloudfrac(ig) 
     &                   + cloudfrac(ig,l)*ponder
               norm=norm+ponder
            ENDDO
            totcloudfrac(ig)=MAX(totcloudfrac(ig)/norm,1.e-12) ! min value if NaNs
         ENDDO
c-----Effective tracers quantities in the cloud fraction
         IF (microphys) THEN
            pqeff(:,:,igcm_ccn_mass)=pq(:,:,igcm_ccn_mass)/
     &                              cloudfrac(:,:)
            pqeff(:,:,igcm_ccn_number)=pq(:,:,igcm_ccn_number)/
     &                              cloudfrac(:,:)
         END IF ! end if (microphys)
         pqeff(:,:,igcm_h2o_ice)=pq(:,:,igcm_h2o_ice)/
     &                           cloudfrac(:,:)
      !! CLFvarying outputs
      CALL WRITEDIAGFI(ngrid,'pqeffice','pqeffice',
     &             'kg/kg',3,pqeff(:,:,igcm_h2o_ice))
      CALL WRITEDIAGFI(ngrid,'pteff','pteff',
     &             'K',3,pteff(:,:))
      CALL WRITEDIAGFI(ngrid,'tcond','tcond',
     &             'K',3,tcond(:,:))
      CALL WRITEDIAGFI(ngrid,'cloudfrac','cloudfrac',
     &             'K',3,cloudfrac(:,:))
      END IF ! end if (CLFvarying)
c------------------------------------------------------------------
c Time subsampling for microphysics
c------------------------------------------------------------------ 
      rhocloud(1:ngrid,1:nlay) = rho_dust
      DO microstep=1,imicro 
      
c-------------------------------------------------------------------
c   1.  Tendencies: 
c------------------


c------ Temperature tendency subpdt
        ! Each microtimestep we give the cloud scheme a stepped entry subpdt instead of pdt
        ! If imicro=1 subpdt is the same as pdt
        DO l=1,nlay
          DO ig=1,ngrid
             sum_subpdt(ig,l) = sum_subpdt(ig,l)
     &        + pdt(ig,l) ! At each micro timestep we add pdt in order to have a stepped entry
          ENDDO
        ENDDO
c------ Tracers tendencies subpdq are additionned
c------ At each micro timestep we add pdq in order to have a stepped entry
        IF (microphys) THEN
          DO l=1,nlay
            DO ig=1,ngrid
              sum_subpdq(ig,l,igcm_dust_mass) = 
     &            sum_subpdq(ig,l,igcm_dust_mass)
     &          + pdq(ig,l,igcm_dust_mass)
              sum_subpdq(ig,l,igcm_dust_number) = 
     &            sum_subpdq(ig,l,igcm_dust_number)
     &          + pdq(ig,l,igcm_dust_number)
              sum_subpdq(ig,l,igcm_ccn_mass) = 
     &            sum_subpdq(ig,l,igcm_ccn_mass)
     &          + pdq(ig,l,igcm_ccn_mass)
              sum_subpdq(ig,l,igcm_ccn_number) = 
     &            sum_subpdq(ig,l,igcm_ccn_number)
     &          + pdq(ig,l,igcm_ccn_number)
            ENDDO
          ENDDO
        ENDIF
        DO l=1,nlay
          DO ig=1,ngrid
            sum_subpdq(ig,l,igcm_h2o_ice) = 
     &          sum_subpdq(ig,l,igcm_h2o_ice)
     &        + pdq(ig,l,igcm_h2o_ice)
            sum_subpdq(ig,l,igcm_h2o_vap) = 
     &          sum_subpdq(ig,l,igcm_h2o_vap)
     &        + pdq(ig,l,igcm_h2o_vap)
          ENDDO
        ENDDO      
        
c-------------------------------------------------------------------
c   2.  Main call to the different cloud schemes:
c------------------------------------------------
        IF (microphys) THEN
           CALL improvedclouds(ngrid,nlay,microtimestep,
     &             pplay,pteff,sum_subpdt, 
     &             pqeff,sum_subpdq,subpdqcloud,subpdtcloud,
     &             nq,tauscaling)

        ELSE
           CALL simpleclouds(ngrid,nlay,microtimestep,
     &             pplay,pzlay,pteff,sum_subpdt,
     &             pqeff,sum_subpdq,subpdqcloud,subpdtcloud,
     &             nq,tau,rice)
        ENDIF

c-------------------------------------------------------------------
c   3.  Updating tendencies after cloud scheme:
c-----------------------------------------------

        IF (microphys) THEN
          DO l=1,nlay
            DO ig=1,ngrid
              sum_subpdq(ig,l,igcm_dust_mass) =
     &            sum_subpdq(ig,l,igcm_dust_mass)
     &          + subpdqcloud(ig,l,igcm_dust_mass)
              sum_subpdq(ig,l,igcm_dust_number) =
     &            sum_subpdq(ig,l,igcm_dust_number)
     &          + subpdqcloud(ig,l,igcm_dust_number)
              sum_subpdq(ig,l,igcm_ccn_mass) =
     &            sum_subpdq(ig,l,igcm_ccn_mass)
     &          + subpdqcloud(ig,l,igcm_ccn_mass)
              sum_subpdq(ig,l,igcm_ccn_number) =
     &            sum_subpdq(ig,l,igcm_ccn_number)
     &          + subpdqcloud(ig,l,igcm_ccn_number)
            ENDDO
          ENDDO
        ENDIF
        DO l=1,nlay
          DO ig=1,ngrid
            sum_subpdq(ig,l,igcm_h2o_ice) =
     &          sum_subpdq(ig,l,igcm_h2o_ice)
     &        + subpdqcloud(ig,l,igcm_h2o_ice)
            sum_subpdq(ig,l,igcm_h2o_vap) =
     &          sum_subpdq(ig,l,igcm_h2o_vap)
     &        + subpdqcloud(ig,l,igcm_h2o_vap)
          ENDDO
        ENDDO
        
        IF (activice) THEN
          DO l=1,nlay
            DO ig=1,ngrid
              sum_subpdt(ig,l) =
     &            sum_subpdt(ig,l) + subpdtcloud(ig,l)
            ENDDO
          ENDDO
        ENDIF
     
 
      ENDDO ! of DO microstep=1,imicro
      
c-------------------------------------------------------------------
c   6.  Compute final tendencies after time loop:
c------------------------------------------------

c------ Temperature tendency pdtcloud
       DO l=1,nlay
         DO ig=1,ngrid
             pdtcloud(ig,l) =
     &         sum_subpdt(ig,l)/real(imicro)-pdt(ig,l)
          ENDDO
       ENDDO
       
c------ Tracers tendencies pdqcloud
       DO l=1,nlay
         DO ig=1,ngrid
            pdqcloud(ig,l,igcm_h2o_ice) = 
     &        sum_subpdq(ig,l,igcm_h2o_ice)/real(imicro) 
     &       - pdq(ig,l,igcm_h2o_ice)
            pdqcloud(ig,l,igcm_h2o_vap) = 
     &        sum_subpdq(ig,l,igcm_h2o_vap)/real(imicro) 
     &       - pdq(ig,l,igcm_h2o_vap)
         ENDDO
       ENDDO
       
       IF(microphys) THEN
        DO l=1,nlay
         DO ig=1,ngrid
            pdqcloud(ig,l,igcm_ccn_mass) = 
     &        sum_subpdq(ig,l,igcm_ccn_mass)/real(imicro) 
     &       - pdq(ig,l,igcm_ccn_mass)
            pdqcloud(ig,l,igcm_ccn_number) = 
     &        sum_subpdq(ig,l,igcm_ccn_number)/real(imicro) 
     &       - pdq(ig,l,igcm_ccn_number)
         ENDDO
        ENDDO
       ENDIF
       
       IF(scavenging) THEN
        DO l=1,nlay
         DO ig=1,ngrid
            pdqcloud(ig,l,igcm_dust_mass) = 
     &        sum_subpdq(ig,l,igcm_dust_mass)/real(imicro) 
     &       - pdq(ig,l,igcm_dust_mass)
            pdqcloud(ig,l,igcm_dust_number) = 
     &        sum_subpdq(ig,l,igcm_dust_number)/real(imicro)
     &       - pdq(ig,l,igcm_dust_number)
         ENDDO
        ENDDO
       ENDIF

c------- Due to stepped entry, other processes tendencies can add up to negative values
c------- Therefore, enforce positive values and conserve mass
       IF(microphys) THEN
        DO l=1,nlay
         DO ig=1,ngrid
          IF ((pq(ig,l,igcm_ccn_number) + 
     &      ptimestep* (pdq(ig,l,igcm_ccn_number) + 
     &        pdqcloud(ig,l,igcm_ccn_number)) .le. 1.)
     &   .or. (pq(ig,l,igcm_ccn_mass) + 
     &      ptimestep* (pdq(ig,l,igcm_ccn_mass) + 
     &        pdqcloud(ig,l,igcm_ccn_mass)) .le. 1.e-20)) THEN
         pdqcloud(ig,l,igcm_ccn_number) =
     &     - pq(ig,l,igcm_ccn_number)/ptimestep 
     &     - pdq(ig,l,igcm_ccn_number) + 1.
         pdqcloud(ig,l,igcm_dust_number) =  
     &     -pdqcloud(ig,l,igcm_ccn_number)
         pdqcloud(ig,l,igcm_ccn_mass) =
     &     - pq(ig,l,igcm_ccn_mass)/ptimestep
     &     - pdq(ig,l,igcm_ccn_mass) + 1.e-20
         pdqcloud(ig,l,igcm_dust_mass) = 
     &     -pdqcloud(ig,l,igcm_ccn_mass)
          ENDIF
         ENDDO
        ENDDO
       ENDIF

       IF(scavenging) THEN
        DO l=1,nlay
         DO ig=1,ngrid
          IF ((pq(ig,l,igcm_dust_number) + 
     &      ptimestep* (pdq(ig,l,igcm_dust_number) + 
     &        pdqcloud(ig,l,igcm_dust_number)) .le. 1.)
     &   .or. (pq(ig,l,igcm_dust_mass) + 
     &      ptimestep* (pdq(ig,l,igcm_dust_mass) + 
     &        pdqcloud(ig,l,igcm_dust_mass)) .le. 1.e-20)) THEN
         pdqcloud(ig,l,igcm_dust_number) =
     &     - pq(ig,l,igcm_dust_number)/ptimestep 
     &     - pdq(ig,l,igcm_dust_number) + 1.
         pdqcloud(ig,l,igcm_ccn_number) =  
     &     -pdqcloud(ig,l,igcm_dust_number)
         pdqcloud(ig,l,igcm_dust_mass) =
     &     - pq(ig,l,igcm_dust_mass)/ptimestep
     &     - pdq(ig,l,igcm_dust_mass) + 1.e-20
         pdqcloud(ig,l,igcm_ccn_mass) = 
     &     -pdqcloud(ig,l,igcm_dust_mass)
          ENDIF
         ENDDO
        ENDDO
       ENDIF

        DO l=1,nlay
         DO ig=1,ngrid
          IF (pq(ig,l,igcm_h2o_ice) + ptimestep*
     &       (pdq(ig,l,igcm_h2o_ice) + pdqcloud(ig,l,igcm_h2o_ice)) 
     &       .le. 1.e-8) THEN
           pdqcloud(ig,l,igcm_h2o_ice) = 
     &     - pq(ig,l,igcm_h2o_ice)/ptimestep - pdq(ig,l,igcm_h2o_ice)
           pdqcloud(ig,l,igcm_h2o_vap) = -pdqcloud(ig,l,igcm_h2o_ice)
          ENDIF 
          IF (pq(ig,l,igcm_h2o_vap) + ptimestep*
     &       (pdq(ig,l,igcm_h2o_vap) + pdqcloud(ig,l,igcm_h2o_vap)) 
     &       .le. 1.e-8) THEN
           pdqcloud(ig,l,igcm_h2o_vap) = 
     &     - pq(ig,l,igcm_h2o_vap)/ptimestep - pdq(ig,l,igcm_h2o_vap)
           pdqcloud(ig,l,igcm_h2o_ice) = -pdqcloud(ig,l,igcm_h2o_vap)
          ENDIF
         ENDDO
        ENDDO


c------Update the ice and dust particle size "rice" for output or photochemistry
c------Only rsedcloud is used for the water cycle

      IF(scavenging) THEN 
        DO l=1, nlay
         DO ig=1,ngrid

        call updaterdust(
     &    pq(ig,l,igcm_dust_mass) +                   ! dust mass
     &   (pdq(ig,l,igcm_dust_mass) +                  ! dust mass
     &    pdqcloud(ig,l,igcm_dust_mass))*ptimestep,   ! dust mass
     &    pq(ig,l,igcm_dust_number) +                 ! dust number
     &   (pdq(ig,l,igcm_dust_number) +                ! dust number
     &    pdqcloud(ig,l,igcm_dust_number))*ptimestep, ! dust number
     &    rdust(ig,l))

         ENDDO
        ENDDO
      ENDIF

      IF(microphys) THEN

       ! In case one does not want to allow supersatured water when using microphysics.
       ! Not done by default.
       IF(.not.supersat) THEN     
        zt  = pt + (pdt+pdtcloud)*ptimestep
        call watersat(ngrid*nlay,zt,pplay,zqsat)
        DO l=1, nlay
         DO ig=1,ngrid
          IF (pq(ig,l,igcm_h2o_vap)
     &      + (pdq(ig,l,igcm_h2o_vap) + pdqcloud(ig,l,igcm_h2o_vap))
     &      * ptimestep .ge. zqsat(ig,l)) THEN
             pdqcloud(ig,l,igcm_h2o_vap) = 
     &         (zqsat(ig,l) - pq(ig,l,igcm_h2o_vap))/ptimestep
     &        - pdq(ig,l,igcm_h2o_vap)
             pdqcloud(ig,l,igcm_h2o_ice) = 
     &         -pdqcloud(ig,l,igcm_h2o_vap)
             ! no need to correct ccn_number, updaterad can handle this properly.
          ENDIF
         ENDDO
        ENDDO        
       ENDIF
       
       DO l=1, nlay
         DO ig=1,ngrid

        call updaterice_micro(
     &    pq(ig,l,igcm_h2o_ice) +                    ! ice mass
     &   (pdq(ig,l,igcm_h2o_ice) +                   ! ice mass
     &    pdqcloud(ig,l,igcm_h2o_ice))*ptimestep,    ! ice mass
     &    pq(ig,l,igcm_ccn_mass) +                   ! ccn mass
     &   (pdq(ig,l,igcm_ccn_mass) +                  ! ccn mass
     &    pdqcloud(ig,l,igcm_ccn_mass))*ptimestep,   ! ccn mass
     &    pq(ig,l,igcm_ccn_number) +                 ! ccn number
     &   (pdq(ig,l,igcm_ccn_number) +                ! ccn number
     &    pdqcloud(ig,l,igcm_ccn_number))*ptimestep, ! ccn number
     &    tauscaling(ig),rice(ig,l),rhocloud(ig,l))
          
         ENDDO
       ENDDO
        
      ELSE ! no microphys
        
        DO l=1,nlay
          DO ig=1,ngrid
          
        call updaterice_typ(
     &    pq(ig,l,igcm_h2o_ice) +                    ! ice mass
     &   (pdq(ig,l,igcm_h2o_ice) +                   ! ice mass
     &    pdqcloud(ig,l,igcm_h2o_ice))*ptimestep,    ! ice mass
     &    tau(ig,1),pzlay(ig,l),rice(ig,l)) 

          ENDDO
         ENDDO
       
       ENDIF ! of IF(microphys)
      
      
      
c     A correction if a lot of subliming CO2 fills the 1st layer FF04/2005
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     Then that should not affect the ice particle radius
      do ig=1,ngrid
        if(pdpsrf(ig)*ptimestep.gt.0.9*(pplev(ig,1)-pplev(ig,2)))then
          if(pdpsrf(ig)*ptimestep.gt.0.9*(pplev(ig,1)-pplev(ig,3)))
     &    rice(ig,2)=rice(ig,3) 
          rice(ig,1)=rice(ig,2)
        end if
      end do
       
       
       DO l=1,nlay
         DO ig=1,ngrid
           rsedcloud(ig,l)=max(rice(ig,l)*
     &                 (1.+nuice_sed)*(1.+nuice_sed)*(1.+nuice_sed),
     &                    rdust(ig,l))
!          rsedcloud(ig,l)=min(rsedcloud(ig,l),1.e-4)
         ENDDO
       ENDDO
       
! used for rad. transfer calculations
! nuice is constant because a lognormal distribution is prescribed
      nuice(1:ngrid,1:nlay)=nuice_ref 

c------Update tendencies for sub-grid water ice clouds
      IF (CLFvarying) THEN
        DO ig=1,ngrid
          DO l=1,nlay
            pdqcloud(ig,l,igcm_dust_mass)=pdqcloud(ig,l,igcm_dust_mass)
     &          *cloudfrac(ig,l)
            pdqcloud(ig,l,igcm_ccn_mass)=pdqcloud(ig,l,igcm_ccn_mass)
     &          *cloudfrac(ig,l)
            pdqcloud(ig,l,igcm_dust_number)=pdqcloud(ig,l,
     &           igcm_dust_number) *cloudfrac(ig,l)
            pdqcloud(ig,l,igcm_ccn_number)=pdqcloud(ig,l,
     &           igcm_ccn_number) *cloudfrac(ig,l)
            pdqcloud(ig,l,igcm_h2o_vap)=pdqcloud(ig,l,
     &           igcm_h2o_vap) *cloudfrac(ig,l)
            pdqcloud(ig,l,igcm_h2o_ice)=pdqcloud(ig,l,
     &           igcm_h2o_ice) *cloudfrac(ig,l)
          ENDDO
        ENDDO   
        pdtcloud(:,:)=pdtcloud(:,:)*cloudfrac(:,:)
      ENDIF
#ifndef MESOSCALE
c=======================================================================
      call WRITEDIAGFI(ngrid,"pdqice2","pdqcloudice apres microphysique"
     &      ,"kg/kg.s-1",3,pdqcloud(1:ngrid,1:nlay,igcm_h2o_ice))
      call WRITEDIAGFI(ngrid,"pdqvap2","pdqcloudvap apres microphysique"
     &      ,"kg/kg.s-1",3,pdqcloud(1:ngrid,1:nlay,
     &      igcm_h2o_vap))
      call WRITEDIAGFI(ngrid,"pdqccn2","pdqcloudccn apres microphysique"
     &      ,"kg/kg.s-1",3,pdqcloud(1:ngrid,1:nlay,
     &      igcm_ccn_mass)) 
      call WRITEDIAGFI(ngrid,"pdqccnN2","pdqcloudccnN apres 
     &      microphysique","nb/kg.s-1",3,pdqcloud(1:ngrid,1:nlay,
     &      igcm_ccn_number))
      call WRITEDIAGFI(ngrid,"pdqdust2", "pdqclouddust apres
     &      microphysique","kg/kg.s-1",3,pdqcloud(1:ngrid,1:nlay,
     &      igcm_dust_mass))
      call WRITEDIAGFI(ngrid,"pdqdustN2", "pdqclouddustN apres
     &      microphysique","nb/kg.s-1",3,pdqcloud(1:ngrid,1:nlay,
     &      igcm_dust_number)) 
c=======================================================================
#endif

      END SUBROUTINE watercloud
      
      subroutine ini_watercloud_mod(ngrid,nlayer,nq)
        implicit none
  
        integer,intent(in) :: ngrid ! number of atmospheric columns
        integer,intent(in) :: nlayer ! number of atmospheric layers
        integer,intent(in) :: nq ! number of tracers

        allocate(zdqcloud(ngrid,nlayer,nq))
        zdqcloud(:,:,:)=0
        allocate(zdqscloud(ngrid,nq))
        zdqscloud(:,:)=0

       end subroutine ini_watercloud_mod


       subroutine end_watercloud_mod
         implicit none

         if (allocated(zdqcloud))      deallocate(zdqcloud)
         if (allocated(zdqscloud))      deallocate(zdqscloud)

       end subroutine end_watercloud_mod

      END MODULE watercloud_mod 
