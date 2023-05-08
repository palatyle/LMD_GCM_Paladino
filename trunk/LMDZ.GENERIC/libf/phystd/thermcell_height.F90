!
!
!
SUBROUTINE thermcell_height(ngrid,nlay,lmin,lmix,lmax,zlev,                   &
                            zmin,zmix,zmax,zw2,wmax,f_star)
      
      
!===============================================================================
!  Purpose: calcul des caracteristiques du thermique: zmax,wmax,zmix
!===============================================================================
      
      IMPLICIT NONE
      
      
!===============================================================================
! Declaration
!===============================================================================
      
!     Inputs:
!     -------
      
      INTEGER, INTENT(in) :: ngrid
      INTEGER, INTENT(in) :: nlay
      INTEGER, INTENT(in) :: lmin(ngrid)
      
      REAL, INTENT(in) :: zlev(ngrid,nlay+1)
      REAL, INTENT(in) :: f_star(ngrid,nlay+1)
      
!     Outputs:
!     --------
      
      INTEGER, INTENT(out) :: lmix(ngrid)
      INTEGER, INTENT(out) :: lmax(ngrid)
      
      REAL, INTENT(out) :: zmin(ngrid)                ! Altitude of the plume bottom
      REAL, INTENT(out) :: zmix(ngrid)                ! Altitude of maximal vertical speed
      REAL, INTENT(out) :: zmax(ngrid)                ! Altitude of the plume top
      REAL, INTENT(out) :: wmax(ngrid)                ! Maximal vertical speed
      
      REAL, INTENT(inout) :: zw2(ngrid,nlay+1)        ! Square vertical speed (becomes its square root)
      
!     Local:
!     ------
      
      INTEGER ig, k
      
      REAL linter(ngrid)                              ! Level (continuous) of maximal vertical speed
      
!===============================================================================
! Initialization
!===============================================================================
      
      linter(:) = lmin(:)
      lmix(:) = lmin(:)
      lmax(:) = lmin(:)
      
      wmax(:) = 0.
      
      zmin(:) = 0.
      zmix(:) = 0.
      zmax(:) = 0.
      
!===============================================================================
! Thermal plume height
!===============================================================================
      
!-------------------------------------------------------------------------------
! Calcul de zmin
!-------------------------------------------------------------------------------
      
      DO ig=1,ngrid
         zmin(ig) = zlev(ig,lmin(ig))
      ENDDO
      
!-------------------------------------------------------------------------------
! Calcul de lmax
!-------------------------------------------------------------------------------
      
      DO ig=1,ngrid
         DO k=nlay,lmin(ig)+1,-1
            IF ((zw2(ig,k) <= 0.).or.(f_star(ig,k) <= 0.)) THEN
               lmax(ig) = k - 1
            ENDIF
         ENDDO
      ENDDO
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! AB: Problematic case where thermal plume reaches the top of the model...
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DO ig=1,ngrid
         IF (zw2(ig,nlay) > 0.) THEN
            print *, 'WARNING: a thermal plume reaches the highest layer!'
            print *, 'ig,k', ig, nlay
            print *, 'zw2', zw2(ig,nlay)
            lmax(ig) = nlay
         ENDIF
      ENDDO
      
!-------------------------------------------------------------------------------
! Calcul de zmax continu via linter
!-------------------------------------------------------------------------------
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! AB: lmin=lmax means the plume is not active and then zw2=0 at each level.
!     Otherwise we have lmax>lmin.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DO ig=1,ngrid
         k = lmax(ig)
         IF (k == lmin(ig)) THEN
            linter(ig) = k
         ELSE
            linter(ig) = k - zw2(ig,k) / (zw2(ig,k+1) - zw2(ig,k))
         ENDIF
         zmax(ig) = zlev(ig,lmax(ig)) + (linter(ig) - lmax(ig))               &
         &        * (zlev(ig,lmax(ig)+1) - zlev(ig,lmax(ig)))
      ENDDO
      
!===============================================================================
! Thermal plume maximal speed and inversion layer
!===============================================================================
      
!-------------------------------------------------------------------------------
! Calcul de lmix et wmax
!-------------------------------------------------------------------------------
      
      DO k=1,nlay
         DO ig=1,ngrid
            IF((k <= lmax(ig)).and.(k > lmin(ig))) THEN
               IF (zw2(ig,k) < 0.) THEN
                  print *, 'ERROR: zw2 has negative value(s)!'
                  print *, 'ig,k', ig, k
                  print *, 'zw2', zw2(ig,k)
                  CALL abort
               ENDIF
               
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! AB: WARNING zw2 becomes its square root!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               zw2(ig,k) = sqrt(zw2(ig,k))
               
               IF (zw2(ig,k) > wmax(ig)) THEN
                  wmax(ig)  = zw2(ig,k)
                  lmix(ig)  = k
               ENDIF
            ELSE
               zw2(ig,k) = 0.
            ENDIF
         ENDDO
      ENDDO
      
!-------------------------------------------------------------------------------
! Calcul de zmix continu (profil parabolique des vitesses)
!-------------------------------------------------------------------------------
      
      DO ig=1,ngrid
         IF (lmix(ig) > lmin(ig)) THEN
            IF (((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))                        &
            &        *((zlev(ig,lmix(ig)))-(zlev(ig,lmix(ig)+1)))             &
            &        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))                   &
            &        *((zlev(ig,lmix(ig)-1))-(zlev(ig,lmix(ig))))) > 1e-10)   &
            &        THEN
               zmix(ig) = ((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))              &
               &        *((zlev(ig,lmix(ig)))**2-(zlev(ig,lmix(ig)+1))**2)    &
               &        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))                &
               &        *((zlev(ig,lmix(ig)-1))**2-(zlev(ig,lmix(ig)))**2))   &
               &        /(2.*((zw2(ig,lmix(ig)-1)-zw2(ig,lmix(ig)))           &
               &        *((zlev(ig,lmix(ig)))-(zlev(ig,lmix(ig)+1)))          &
               &        -(zw2(ig,lmix(ig))-zw2(ig,lmix(ig)+1))                &
               &        *((zlev(ig,lmix(ig)-1))-(zlev(ig,lmix(ig))))))
            ELSE
               zmix(ig) = zlev(ig,lmix(ig))
               print *, 'WARNING: problematic zmix value!'
            ENDIF
         ELSE
            zmix(ig) = 0.
         ENDIF
      ENDDO
      
!-------------------------------------------------------------------------------
! Check consistency between zmax and zmix
!-------------------------------------------------------------------------------
      
      DO ig=1,ngrid
         IF ((zmax(ig)-zmix(ig)) < 0.) THEN
            print *, 'WARNING: we have zmix > zmax!'
            print *, 'zmax,zmix_old,zmix_new', zmax(ig), zmix(ig), 0.9 * zmax(ig)
            zmix(ig) = 0.9 * zmax(ig)
            DO k=1,nlay
               IF ((zmix(ig) >= zlev(ig,k)).and.(zmix(ig) < zlev(ig,k+1))) THEN
                  lmix(ig) = k
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      
      
RETURN
END
