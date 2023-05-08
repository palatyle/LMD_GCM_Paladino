










      SUBROUTINE aeroptproperties(ngrid,nlayer,reffrad,nueffrad,   &
                                  QVISsQREF3d,omegaVIS3d,gVIS3d,   &
                                  QIRsQREF3d,omegaIR3d,gIR3d,      &
                                  QREFvis3d,QREFir3d)!,		   &
!                                  omegaREFvis3d,omegaREFir3d)

      use radinc_h,    only: L_NSPECTI,L_NSPECTV,nsizemax,naerkind
      use radcommon_h, only: QVISsQREF,omegavis,gvis,QIRsQREF,omegair,gir
      use radcommon_h, only: qrefvis,qrefir,omegarefir !,omegarefvis
      use radcommon_h, only: radiustab,nsize

      implicit none

!     =============================================================
!     Aerosol Optical Properties
!
!     Description:
!       Compute the scattering parameters in each grid
!       box, depending on aerosol grain sizes. Log-normal size
!       distribution and Gauss-Legendre integration are used.

!     Parameters:
!       Don't forget to set the value of varyingnueff below; If
!       the effective variance of the distribution for the given
!       aerosol is considered homogeneous in the atmosphere, please
!       set varyingnueff(iaer) to .false. Resulting computational
!       time will be much better.

!     Authors: J.-B. Madeleine, F. Forget, F. Montmessin
!     Slightly modified and converted to F90 by R. Wordsworth (2009)
!     Varying nueff section removed by R. Wordsworth for simplicity
!     ==============================================================

!     Local variables 
!     ---------------



!     =============================================================
      LOGICAL, PARAMETER :: varyingnueff(naerkind) = .false.
!     =============================================================

!     Min. and max radius of the interpolation grid (in METERS)
      REAL, PARAMETER :: refftabmin = 2e-8 !2e-8
!      REAL, PARAMETER :: refftabmax = 35e-6
      REAL, PARAMETER :: refftabmax = 1e-3
!     Log of the min and max variance of the interpolation grid
      REAL, PARAMETER :: nuefftabmin = -4.6
      REAL, PARAMETER :: nuefftabmax = 0.
!     Number of effective radius of the interpolation grid
      INTEGER, PARAMETER :: refftabsize = 200
!     Number of effective variances of the interpolation grid
!      INTEGER, PARAMETER :: nuefftabsize = 100
      INTEGER, PARAMETER :: nuefftabsize = 1
!     Interpolation grid indices (reff,nueff)
      INTEGER :: grid_i,grid_j
!     Intermediate variable
      REAL :: var_tmp,var3d_tmp(ngrid,nlayer)
!     Bilinear interpolation factors
      REAL :: kx,ky,k1,k2,k3,k4
!     Size distribution parameters
      REAL :: sizedistk1,sizedistk2
!     Pi!
      REAL,SAVE :: pi
!$OMP THREADPRIVATE(pi)
!     Variables used by the Gauss-Legendre integration:
      INTEGER radius_id,gausind
      REAL kint
      REAL drad
      INTEGER, PARAMETER :: ngau = 10
      REAL weightgaus(ngau),radgaus(ngau)
      SAVE weightgaus,radgaus
!      DATA weightgaus/.2955242247,.2692667193,.2190863625,.1494513491,.0666713443/
!      DATA radgaus/.1488743389,.4333953941,.6794095682,.8650633666,.9739065285/
      DATA    radgaus/0.07652652113350,0.22778585114165, &
                      0.37370608871528,0.51086700195146, &
                      0.63605368072468,0.74633190646476, &
                      0.83911697181213,0.91223442826796, &
                      0.96397192726078,0.99312859919241/

      DATA weightgaus/0.15275338723120,0.14917298659407, &
                      0.14209610937519,0.13168863843930, &
                      0.11819453196154,0.10193011980823, &
                      0.08327674160932,0.06267204829828, &
                      0.04060142982019,0.01761400714091/
!$OMP THREADPRIVATE(radgaus,weightgaus)
!     Indices
      INTEGER :: i,j,k,l,m,iaer,idomain
      INTEGER :: ig,lg,chg

!     Local saved variables
!     ---------------------

!     Radius axis of the interpolation grid
      REAL,SAVE :: refftab(refftabsize)
!     Variance axis of the interpolation grid
      REAL,SAVE :: nuefftab(nuefftabsize)
!     Volume ratio of the grid
      REAL,SAVE :: logvratgrid,vratgrid
!     Grid used to remember which calculation is done
      LOGICAL,SAVE :: checkgrid(refftabsize,nuefftabsize,naerkind,2) = .false.
!$OMP THREADPRIVATE(refftab,nuefftab,logvratgrid,vratgrid,checkgrid)
!     Optical properties of the grid (VISIBLE)
      REAL,SAVE :: qsqrefVISgrid(refftabsize,nuefftabsize,L_NSPECTV,naerkind)
      REAL,SAVE :: qextVISgrid(refftabsize,nuefftabsize,L_NSPECTV,naerkind)
      REAL,SAVE :: qscatVISgrid(refftabsize,nuefftabsize,L_NSPECTV,naerkind)
      REAL,SAVE :: omegVISgrid(refftabsize,nuefftabsize,L_NSPECTV,naerkind)
      REAL,SAVE :: gVISgrid(refftabsize,nuefftabsize,L_NSPECTV,naerkind)
!$OMP THREADPRIVATE(qsqrefVISgrid,qextVISgrid,qscatVISgrid,omegVISgrid,gVISgrid)
!     Optical properties of the grid (INFRARED)
      REAL,SAVE :: qsqrefIRgrid(refftabsize,nuefftabsize,L_NSPECTI,naerkind)
      REAL,SAVE :: qextIRgrid(refftabsize,nuefftabsize,L_NSPECTI,naerkind)
      REAL,SAVE :: qscatIRgrid(refftabsize,nuefftabsize,L_NSPECTI,naerkind)
      REAL,SAVE :: omegIRgrid(refftabsize,nuefftabsize,L_NSPECTI,naerkind)
      REAL,SAVE :: gIRgrid(refftabsize,nuefftabsize,L_NSPECTI,naerkind)
!$OMP THREADPRIVATE(qsqrefIRgrid,qextIRgrid,qscatIRgrid,omegIRgrid,gIRgrid)
!     Optical properties of the grid (REFERENCE WAVELENGTHS)
      REAL,SAVE :: qrefVISgrid(refftabsize,nuefftabsize,naerkind)
      REAL,SAVE :: qscatrefVISgrid(refftabsize,nuefftabsize,naerkind)
      REAL,SAVE :: qrefIRgrid(refftabsize,nuefftabsize,naerkind)
      REAL,SAVE :: qscatrefIRgrid(refftabsize,nuefftabsize,naerkind)
      REAL,SAVE :: omegrefVISgrid(refftabsize,nuefftabsize,naerkind)
      REAL,SAVE :: omegrefIRgrid(refftabsize,nuefftabsize,naerkind)
!$OMP THREADPRIVATE(qrefVISgrid,qscatrefVISgrid,qrefIRgrid,qscatrefIRgrid,omegrefVISgrid,&
	!$OMP omegrefIRgrid)
!     Firstcall
      LOGICAL,SAVE :: firstcall = .true.
!$OMP THREADPRIVATE(firstcall)
!     Variables used by the Gauss-Legendre integration:
      REAL,SAVE :: normd(refftabsize,nuefftabsize,naerkind,2)
      REAL,SAVE :: dista(refftabsize,nuefftabsize,naerkind,2,ngau)
      REAL,SAVE :: distb(refftabsize,nuefftabsize,naerkind,2,ngau)
!$OMP THREADPRIVATE(normd,dista,distb)

      REAL,SAVE :: radGAUSa(ngau,naerkind,2)
      REAL,SAVE :: radGAUSb(ngau,naerkind,2)
!$OMP THREADPRIVATE(radGAUSa,radGAUSb)

      REAL,SAVE :: qsqrefVISa(L_NSPECTV,ngau,naerkind)
      REAL,SAVE :: qrefVISa(ngau,naerkind)
      REAL,SAVE :: qsqrefVISb(L_NSPECTV,ngau,naerkind)
      REAL,SAVE :: qrefVISb(ngau,naerkind)
      REAL,SAVE :: omegVISa(L_NSPECTV,ngau,naerkind)
      REAL,SAVE :: omegrefVISa(ngau,naerkind)
      REAL,SAVE :: omegVISb(L_NSPECTV,ngau,naerkind)
      REAL,SAVE :: omegrefVISb(ngau,naerkind)
      REAL,SAVE :: gVISa(L_NSPECTV,ngau,naerkind)
      REAL,SAVE :: gVISb(L_NSPECTV,ngau,naerkind)
!$OMP THREADPRIVATE(qsqrefVISa,qrefVISa,qsqrefVISb,qrefVISb,omegVISa, &
	!$OMP omegrefVISa,omegVISb,omegrefVISb,gVISa,gVISb)

      REAL,SAVE :: qsqrefIRa(L_NSPECTI,ngau,naerkind)
      REAL,SAVE :: qrefIRa(ngau,naerkind)
      REAL,SAVE :: qsqrefIRb(L_NSPECTI,ngau,naerkind)
      REAL,SAVE :: qrefIRb(ngau,naerkind)
      REAL,SAVE :: omegIRa(L_NSPECTI,ngau,naerkind)
      REAL,SAVE :: omegrefIRa(ngau,naerkind)
      REAL,SAVE :: omegIRb(L_NSPECTI,ngau,naerkind)
      REAL,SAVE :: omegrefIRb(ngau,naerkind)
      REAL,SAVE :: gIRa(L_NSPECTI,ngau,naerkind)
      REAL,SAVE :: gIRb(L_NSPECTI,ngau,naerkind)
!$OMP THREADPRIVATE(qsqrefIRa,qrefIRa,qsqrefIRb,qrefIRb,omegIRa,omegrefIRa,&
	!$OMP omegIRb,omegrefIRb,gIRa,gIRb)

      REAL :: radiusm
      REAL :: radiusr

!     Inputs
!     ------

      INTEGER :: ngrid,nlayer
!     Aerosol effective radius used for radiative transfer (meter)
      REAL,INTENT(IN) :: reffrad(ngrid,nlayer,naerkind)
!     Aerosol effective variance used for radiative transfer (n.u.)
      REAL,INTENT(IN) :: nueffrad(ngrid,nlayer,naerkind)

!     Outputs
!     -------

      REAL,INTENT(OUT) :: QVISsQREF3d(ngrid,nlayer,L_NSPECTV,naerkind)
      REAL,INTENT(OUT) :: omegaVIS3d(ngrid,nlayer,L_NSPECTV,naerkind)
      REAL,INTENT(OUT) :: gVIS3d(ngrid,nlayer,L_NSPECTV,naerkind)

      REAL,INTENT(OUT) :: QIRsQREF3d(ngrid,nlayer,L_NSPECTI,naerkind)
      REAL,INTENT(OUT) :: omegaIR3d(ngrid,nlayer,L_NSPECTI,naerkind)
      REAL,INTENT(OUT) :: gIR3d(ngrid,nlayer,L_NSPECTI,naerkind)

      REAL,INTENT(OUT) :: QREFvis3d(ngrid,nlayer,naerkind)
      REAL,INTENT(OUT) :: QREFir3d(ngrid,nlayer,naerkind)

!      REAL :: omegaREFvis3d(ngrid,nlayer,naerkind)
!      REAL :: omegaREFir3d(ngrid,nlayer,naerkind)

      DO iaer = 1, naerkind ! Loop on aerosol kind
        IF ( (nsize(iaer,1).EQ.1).AND.(nsize(iaer,2).EQ.1) ) THEN
!==================================================================
!       If there is one single particle size, optical
!         properties of the considered aerosol are homogeneous
          DO lg = 1, nlayer
            DO ig = 1, ngrid
              DO chg = 1, L_NSPECTV
                QVISsQREF3d(ig,lg,chg,iaer)=QVISsQREF(chg,iaer,1)
                omegaVIS3d(ig,lg,chg,iaer)=omegaVIS(chg,iaer,1)
                gVIS3d(ig,lg,chg,iaer)=gVIS(chg,iaer,1)
              ENDDO
              DO chg = 1, L_NSPECTI
                QIRsQREF3d(ig,lg,chg,iaer)=QIRsQREF(chg,iaer,1)
                omegaIR3d(ig,lg,chg,iaer)=omegaIR(chg,iaer,1)
                gIR3d(ig,lg,chg,iaer)=gIR(chg,iaer,1)
              ENDDO
              QREFvis3d(ig,lg,iaer)=QREFvis(iaer,1)
              QREFir3d(ig,lg,iaer)=QREFir(iaer,1)
!              omegaREFvis3d(ig,lg,iaer)=omegaREFvis(iaer,1)
!              omegaREFir3d(ig,lg,iaer)=omegaREFir(iaer,1)
            ENDDO
          ENDDO


          if (firstcall) then
             print*,'Optical prop. of the aerosol are homogenous for:'
             print*,'iaer = ',iaer
          endif

        ELSE ! Varying effective radius and variance
      DO idomain = 1, 2 ! Loop on visible or infrared channel
!==================================================================
!     1. Creating the effective radius and variance grid
!     --------------------------------------------------
      IF (firstcall) THEN

!       1.1 Pi!
        pi = 2. * asin(1.e0)

!       1.2 Effective radius
        refftab(1)    = refftabmin
        refftab(refftabsize) = refftabmax

        logvratgrid = log(refftabmax/refftabmin) / float(refftabsize-1)*3.
        vratgrid = exp(logvratgrid)

        do i = 2, refftabsize-1
          refftab(i) = refftab(i-1)*vratgrid**(1./3.)
        enddo

!       1.3 Effective variance
        if(nuefftabsize.eq.1)then ! addded by RDW
           print*,'Warning: no variance range in aeroptproperties'
           nuefftab(1)=0.2
        else
           do i = 0, nuefftabsize-1
              nuefftab(i+1) = exp( nuefftabmin + i*(nuefftabmax-nuefftabmin)/(nuefftabsize-1) )
           enddo
        endif

        firstcall = .false.
      ENDIF

!       1.4 Radius middle point and range for Gauss integration
        radiusm=0.5*(radiustab(iaer,idomain,nsize(iaer,idomain)) + radiustab(iaer,idomain,1))
        radiusr=0.5*(radiustab(iaer,idomain,nsize(iaer,idomain)) - radiustab(iaer,idomain,1))

!       1.5 Interpolating data at the Gauss quadrature points:
        DO gausind=1,ngau
          drad=radiusr*radgaus(gausind)
          radGAUSa(gausind,iaer,idomain)=radiusm-drad

          radius_id=minloc(abs(radiustab(iaer,idomain,:) - (radiusm-drad)),DIM=1)
          IF ((radiustab(iaer,idomain,radius_id) - (radiusm-drad)).GT.0) THEN
            radius_id=radius_id-1
          ENDIF
          IF (radius_id.GE.nsize(iaer,idomain)) THEN
            radius_id=nsize(iaer,idomain)-1
            kint = 1.
          ELSEIF (radius_id.LT.1) THEN
            radius_id=1
            kint = 0.
          ELSE
          kint = ( (radiusm-drad) -				&
                   radiustab(iaer,idomain,radius_id) ) /	&
                 ( radiustab(iaer,idomain,radius_id+1) -	&
                   radiustab(iaer,idomain,radius_id) )
          ENDIF
          IF (idomain.EQ.1) THEN ! VISIBLE DOMAIN -----------------
            DO m=1,L_NSPECTV
               qsqrefVISa(m,gausind,iaer)=                      &
                    (1-kint)*QVISsQREF(m,iaer,radius_id) +      &
                    kint*QVISsQREF(m,iaer,radius_id+1)
            omegVISa(m,gausind,iaer)=                           &
                    (1-kint)*omegaVIS(m,iaer,radius_id) +       &
                    kint*omegaVIS(m,iaer,radius_id+1)
            gVISa(m,gausind,iaer)=                              &
                    (1-kint)*gVIS(m,iaer,radius_id) +           &
                    kint*gVIS(m,iaer,radius_id+1)
            ENDDO
            qrefVISa(gausind,iaer)=                             &
                    (1-kint)*QREFvis(iaer,radius_id) +          &
                    kint*QREFvis(iaer,radius_id+1)
            omegrefVISa(gausind,iaer)= 0
!            omegrefVISa(gausind,iaer)=                          &
!                    (1-kint)*omegaREFvis(iaer,radius_id) +      &
!                    kint*omegaREFvis(iaer,radius_id+1)
          ELSE ! INFRARED DOMAIN ----------------------------------
            DO m=1,L_NSPECTI
            qsqrefIRa(m,gausind,iaer)=                          &
                    (1-kint)*QIRsQREF(m,iaer,radius_id) +       &
                    kint*QIRsQREF(m,iaer,radius_id+1)
            omegIRa(m,gausind,iaer)=                            &
                    (1-kint)*omegaIR(m,iaer,radius_id) +        &
                    kint*omegaIR(m,iaer,radius_id+1)
            gIRa(m,gausind,iaer)=                               &
                    (1-kint)*gIR(m,iaer,radius_id) +            &
                    kint*gIR(m,iaer,radius_id+1)
            ENDDO
            qrefIRa(gausind,iaer)=                              &
                    (1-kint)*QREFir(iaer,radius_id) +           &
                    kint*QREFir(iaer,radius_id+1)
            omegrefIRa(gausind,iaer)=                           &
                    (1-kint)*omegaREFir(iaer,radius_id) +       &
                    kint*omegaREFir(iaer,radius_id+1)
          ENDIF
        ENDDO

        DO gausind=1,ngau
          drad=radiusr*radgaus(gausind)
          radGAUSb(gausind,iaer,idomain)=radiusm+drad

          radius_id=minloc(abs(radiustab(iaer,idomain,:) -      &
                               (radiusm+drad)),DIM=1)
          IF ((radiustab(iaer,idomain,radius_id) -              &
               (radiusm+drad)).GT.0) THEN
            radius_id=radius_id-1
          ENDIF
          IF (radius_id.GE.nsize(iaer,idomain)) THEN
            radius_id=nsize(iaer,idomain)-1
            kint = 1.
          ELSEIF (radius_id.LT.1) THEN
            radius_id=1
            kint = 0.
          ELSE
            kint = ( (radiusm+drad) -                           &
                     radiustab(iaer,idomain,radius_id) ) /      &
                   ( radiustab(iaer,idomain,radius_id+1) -      &
                     radiustab(iaer,idomain,radius_id) )
          ENDIF
          IF (idomain.EQ.1) THEN ! VISIBLE DOMAIN -----------------
            DO m=1,L_NSPECTV
            qsqrefVISb(m,gausind,iaer)=                         &
                    (1-kint)*QVISsQREF(m,iaer,radius_id) +      &
                    kint*QVISsQREF(m,iaer,radius_id+1)    
            omegVISb(m,gausind,iaer)=                           &
                    (1-kint)*omegaVIS(m,iaer,radius_id) +       &
                    kint*omegaVIS(m,iaer,radius_id+1)
            gVISb(m,gausind,iaer)=                              &
                    (1-kint)*gVIS(m,iaer,radius_id) +           &
                    kint*gVIS(m,iaer,radius_id+1)
            ENDDO
            qrefVISb(gausind,iaer)=                             &
                    (1-kint)*QREFvis(iaer,radius_id) +          &
                    kint*QREFvis(iaer,radius_id+1)
            omegrefVISb(gausind,iaer)= 0
!            omegrefVISb(gausind,iaer)=                          &
!                    (1-kint)*omegaREFvis(iaer,radius_id) +      &
!                    kint*omegaREFvis(iaer,radius_id+1)
          ELSE ! INFRARED DOMAIN ----------------------------------
            DO m=1,L_NSPECTI
            qsqrefIRb(m,gausind,iaer)=                          &
                    (1-kint)*QIRsQREF(m,iaer,radius_id) +       &
                    kint*QIRsQREF(m,iaer,radius_id+1)
            omegIRb(m,gausind,iaer)=                            &
                    (1-kint)*omegaIR(m,iaer,radius_id) +        &
                    kint*omegaIR(m,iaer,radius_id+1)
            gIRb(m,gausind,iaer)=                               &
                    (1-kint)*gIR(m,iaer,radius_id) +            &
                    kint*gIR(m,iaer,radius_id+1)
            ENDDO
            qrefIRb(gausind,iaer)=                              &
                    (1-kint)*QREFir(iaer,radius_id) +           &
                    kint*QREFir(iaer,radius_id+1)
            omegrefIRb(gausind,iaer)=                           &
                    (1-kint)*omegaREFir(iaer,radius_id) +       &
                    kint*omegaREFir(iaer,radius_id+1)
          ENDIF
        ENDDO

!==================================================================
! CONSTANT NUEFF FROM HERE
!==================================================================

!     2. Compute the scattering parameters using linear
!       interpolation over grain sizes and constant nueff
!     ---------------------------------------------------

      DO lg = 1,nlayer
        DO ig = 1, ngrid
!         2.1 Effective radius index and kx calculation
          var_tmp=reffrad(ig,lg,iaer)/refftabmin
          var_tmp=log(var_tmp)*3.
          var_tmp=var_tmp/logvratgrid+1.
          grid_i=floor(var_tmp)
          IF (grid_i.GE.refftabsize) THEN
!           WRITE(*,*) 'Warning: particle size in grid box #'
!           WRITE(*,*) ig,' is too large to be used by the '
!           WRITE(*,*) 'radiative transfer; please extend the '
!           WRITE(*,*) 'interpolation grid to larger grain sizes.'
            grid_i=refftabsize-1
            kx = 1.
          ELSEIF (grid_i.LT.1) THEN
!           WRITE(*,*) 'Warning: particle size in grid box #'
!           WRITE(*,*) ig,' is too small to be used by the '
!           WRITE(*,*) 'radiative transfer; please extend the '
!           WRITE(*,*) 'interpolation grid to smaller grain sizes.'
            grid_i=1
            kx = 0.
          ELSE
            kx = ( reffrad(ig,lg,iaer)-refftab(grid_i) ) /            &
                 ( refftab(grid_i+1)-refftab(grid_i) )
          ENDIF
!         2.3 Integration
          DO j=grid_i,grid_i+1
!             2.3.1 Check if the calculation has been done
              IF (.NOT.checkgrid(j,1,iaer,idomain)) THEN
!               2.3.2 Log-normal dist., r_g and sigma_g are defined
!                 in [hansen_1974], "Light scattering in planetary
!                 atmospheres", Space Science Reviews 16 527-610.
!                 Here, sizedistk1=r_g and sizedistk2=sigma_g^2
                sizedistk2 = log(1.+nueffrad(1,1,iaer))
                sizedistk1 = exp(2.5*sizedistk2)
                sizedistk1 = refftab(j) / sizedistk1

                normd(j,1,iaer,idomain) = 1e-30
                DO gausind=1,ngau
                  drad=radiusr*radgaus(gausind)
                  dista(j,1,iaer,idomain,gausind) = LOG((radiusm-drad)/sizedistk1)
                  dista(j,1,iaer,idomain,gausind) =                   &
                    EXP(-dista(j,1,iaer,idomain,gausind) *            &
                    dista(j,1,iaer,idomain,gausind) *                 &
                    0.5e0/sizedistk2)/(radiusm-drad)                  
                  dista(j,1,iaer,idomain,gausind) =                   &
                    dista(j,1,iaer,idomain,gausind) /                 &
                    (sqrt(2e0*pi*sizedistk2))

                  distb(j,1,iaer,idomain,gausind) = LOG((radiusm+drad)/sizedistk1)
                  distb(j,1,iaer,idomain,gausind) =                   &
                    EXP(-distb(j,1,iaer,idomain,gausind) *            &
                    distb(j,1,iaer,idomain,gausind) *                 &
                    0.5e0/sizedistk2)/(radiusm+drad)
                  distb(j,1,iaer,idomain,gausind) =                   &
                    distb(j,1,iaer,idomain,gausind) /                 &
                    (sqrt(2e0*pi*sizedistk2))

                  normd(j,1,iaer,idomain)=normd(j,1,iaer,idomain) +   &
                    weightgaus(gausind) *                             &
                    (                                                 &
                    distb(j,1,iaer,idomain,gausind) * pi *            &
                    radGAUSb(gausind,iaer,idomain) *                  &
                    radGAUSb(gausind,iaer,idomain) +                  &
                    dista(j,1,iaer,idomain,gausind) * pi *            &
                    radGAUSa(gausind,iaer,idomain) *                  &
                    radGAUSa(gausind,iaer,idomain)                    &
                    )
                ENDDO
                IF (idomain.EQ.1) THEN ! VISIBLE DOMAIN -----------
!                 2.3.3.vis Initialization
                  qsqrefVISgrid(j,1,:,iaer)=0.
                  qextVISgrid(j,1,:,iaer)=0.
                  qscatVISgrid(j,1,:,iaer)=0.
                  omegVISgrid(j,1,:,iaer)=0.
                  gVISgrid(j,1,:,iaer)=0.
                  qrefVISgrid(j,1,iaer)=0.
                  qscatrefVISgrid(j,1,iaer)=0.
                  omegrefVISgrid(j,1,iaer)=0.

                  DO gausind=1,ngau
                    DO m=1,L_NSPECTV
!                     Convolution:
                      qextVISgrid(j,1,m,iaer) =              &
                        qextVISgrid(j,1,m,iaer) +            & 
                        weightgaus(gausind) *                &
                        (                                    &
                        qsqrefVISb(m,gausind,iaer) *         &
                        qrefVISb(gausind,iaer) *             &
                        pi*radGAUSb(gausind,iaer,idomain) *  &
                        radGAUSb(gausind,iaer,idomain) *     &
                        distb(j,1,iaer,idomain,gausind) +    &
                        qsqrefVISa(m,gausind,iaer) *         &
                        qrefVISa(gausind,iaer) *             &
                        pi*radGAUSa(gausind,iaer,idomain) *  &
                        radGAUSa(gausind,iaer,idomain) *     &
                        dista(j,1,iaer,idomain,gausind)      &
                        )
                      qscatVISgrid(j,1,m,iaer) =             &
                        qscatVISgrid(j,1,m,iaer) +           &
                        weightgaus(gausind) *                &
                        (                                    &
                        omegVISb(m,gausind,iaer) *           &
                        qsqrefVISb(m,gausind,iaer) *         &
                        qrefVISb(gausind,iaer) *             &
                        pi*radGAUSb(gausind,iaer,idomain) *  &
                        radGAUSb(gausind,iaer,idomain) *     &
                        distb(j,1,iaer,idomain,gausind) +    &
                        omegVISa(m,gausind,iaer) *           &
                        qsqrefVISa(m,gausind,iaer) *         &
                        qrefVISa(gausind,iaer) *             &
                        pi*radGAUSa(gausind,iaer,idomain) *  &
                        radGAUSa(gausind,iaer,idomain) *     &
                        dista(j,1,iaer,idomain,gausind)      &
                        )
                      gVISgrid(j,1,m,iaer) =                 &
                        gVISgrid(j,1,m,iaer) +               &
                        weightgaus(gausind) *                &
                        (                                    &
                        omegVISb(m,gausind,iaer) *           &
                        qsqrefVISb(m,gausind,iaer) *         &
                        qrefVISb(gausind,iaer) *             &
                        gVISb(m,gausind,iaer) *              &
                        pi*radGAUSb(gausind,iaer,idomain) *  &
                        radGAUSb(gausind,iaer,idomain) *     &
                        distb(j,1,iaer,idomain,gausind) +    &
                        omegVISa(m,gausind,iaer) *           &
                        qsqrefVISa(m,gausind,iaer) *         &
                        qrefVISa(gausind,iaer) *             &
                        gVISa(m,gausind,iaer) *              &
                        pi*radGAUSa(gausind,iaer,idomain) *  &
                        radGAUSa(gausind,iaer,idomain) *     &
                        dista(j,1,iaer,idomain,gausind)      &
                        )
                    ENDDO
                    qrefVISgrid(j,1,iaer) =                  &
                      qrefVISgrid(j,1,iaer) +                &
                      weightgaus(gausind) *                  &
                      (                                      &
                      qrefVISb(gausind,iaer) *               &
                      pi*radGAUSb(gausind,iaer,idomain) *    &
                      radGAUSb(gausind,iaer,idomain) *       &
                      distb(j,1,iaer,idomain,gausind) +      &
                      qrefVISa(gausind,iaer) *               &
                      pi*radGAUSa(gausind,iaer,idomain) *    &
                      radGAUSa(gausind,iaer,idomain) *       &
                      dista(j,1,iaer,idomain,gausind)        &
                      )
                    qscatrefVISgrid(j,1,iaer) =              &
                      qscatrefVISgrid(j,1,iaer) +            &
                      weightgaus(gausind) *                  &
                      (                                      &
                      omegrefVISb(gausind,iaer) *            &
                      qrefVISb(gausind,iaer) *               & 
                      pi*radGAUSb(gausind,iaer,idomain) *    &
                      radGAUSb(gausind,iaer,idomain) *       &
                      distb(j,1,iaer,idomain,gausind) +      &
                      omegrefVISa(gausind,iaer) *            &
                      qrefVISa(gausind,iaer) *               &
                      pi*radGAUSa(gausind,iaer,idomain) *    &
                      radGAUSa(gausind,iaer,idomain) *       &
                      dista(j,1,iaer,idomain,gausind)        &
                      )
                  ENDDO

                  qrefVISgrid(j,1,iaer)=qrefVISgrid(j,1,iaer) /          &
                                normd(j,1,iaer,idomain)       
                  qscatrefVISgrid(j,1,iaer)=qscatrefVISgrid(j,1,iaer) /  &
                                normd(j,1,iaer,idomain)
                  omegrefVISgrid(j,1,iaer)=qscatrefVISgrid(j,1,iaer) /   &
                               qrefVISgrid(j,1,iaer)
                  DO m=1,L_NSPECTV
                    qextVISgrid(j,1,m,iaer)=qextVISgrid(j,1,m,iaer) /    &
                                normd(j,1,iaer,idomain)
                    qscatVISgrid(j,1,m,iaer)=qscatVISgrid(j,1,m,iaer) /  &
                                normd(j,1,iaer,idomain)
                    gVISgrid(j,1,m,iaer)=gVISgrid(j,1,m,iaer) /          &
                                qscatVISgrid(j,1,m,iaer) /               &
                                normd(j,1,iaer,idomain)

                    qsqrefVISgrid(j,1,m,iaer)=qextVISgrid(j,1,m,iaer) /  &
                                qrefVISgrid(j,1,iaer)
                    omegVISgrid(j,1,m,iaer)=qscatVISgrid(j,1,m,iaer) /   &
                                qextVISgrid(j,1,m,iaer)
                  ENDDO
                ELSE                   ! INFRARED DOMAIN ----------
!                 2.3.3.ir Initialization
                  qsqrefIRgrid(j,1,:,iaer)=0.
                  qextIRgrid(j,1,:,iaer)=0.
                  qscatIRgrid(j,1,:,iaer)=0.
                  omegIRgrid(j,1,:,iaer)=0.
                  gIRgrid(j,1,:,iaer)=0.
                  qrefIRgrid(j,1,iaer)=0.
                  qscatrefIRgrid(j,1,iaer)=0.
                  omegrefIRgrid(j,1,iaer)=0.

                  DO gausind=1,ngau
                    DO m=1,L_NSPECTI
!                     Convolution:
                      qextIRgrid(j,1,m,iaer) =                  &
                        qextIRgrid(j,1,m,iaer) +                &
                        weightgaus(gausind) *                   &
                        (                                       &
                        qsqrefIRb(m,gausind,iaer) *             &
                        qrefVISb(gausind,iaer) *                &
                        pi*radGAUSb(gausind,iaer,idomain) *     &
                        radGAUSb(gausind,iaer,idomain) *        &
                        distb(j,1,iaer,idomain,gausind) +       &
                        qsqrefIRa(m,gausind,iaer) *             &
                        qrefVISa(gausind,iaer) *                &
                        pi*radGAUSa(gausind,iaer,idomain) *     &
                        radGAUSa(gausind,iaer,idomain) *        &
                        dista(j,1,iaer,idomain,gausind)         &
                        )
                      qscatIRgrid(j,1,m,iaer) =                 &
                        qscatIRgrid(j,1,m,iaer) +               &
                        weightgaus(gausind) *                   &
                        (                                       &
                        omegIRb(m,gausind,iaer) *               &
                        qsqrefIRb(m,gausind,iaer) *             &
                        qrefVISb(gausind,iaer) *                &
                        pi*radGAUSb(gausind,iaer,idomain) *     &
                        radGAUSb(gausind,iaer,idomain) *        &
                        distb(j,1,iaer,idomain,gausind) +       &
                        omegIRa(m,gausind,iaer) *               &
                        qsqrefIRa(m,gausind,iaer) *             &
                        qrefVISa(gausind,iaer) *                &
                        pi*radGAUSa(gausind,iaer,idomain) *     &
                        radGAUSa(gausind,iaer,idomain) *        &
                        dista(j,1,iaer,idomain,gausind)         &
                        )
                      gIRgrid(j,1,m,iaer) =                     &
                        gIRgrid(j,1,m,iaer) +                   &
                        weightgaus(gausind) *                   &
                        (                                       &
                        omegIRb(m,gausind,iaer) *               &
                        qsqrefIRb(m,gausind,iaer) *             &
                        qrefVISb(gausind,iaer) *                &
                        gIRb(m,gausind,iaer) *                  &
                        pi*radGAUSb(gausind,iaer,idomain) *     &
                        radGAUSb(gausind,iaer,idomain) *        &
                        distb(j,1,iaer,idomain,gausind) +       &
                        omegIRa(m,gausind,iaer) *               &
                        qsqrefIRa(m,gausind,iaer) *             &
                        qrefVISa(gausind,iaer) *                &
                        gIRa(m,gausind,iaer) *                  &
                        pi*radGAUSa(gausind,iaer,idomain) *     &
                        radGAUSa(gausind,iaer,idomain) *        &
                        dista(j,1,iaer,idomain,gausind)         &
                        )
                    ENDDO
                    qrefIRgrid(j,1,iaer) =                      &
                      qrefIRgrid(j,1,iaer) +                    &
                      weightgaus(gausind) *                     &
                      (                                         &
                      qrefIRb(gausind,iaer) *                   &
                      pi*radGAUSb(gausind,iaer,idomain) *       &
                      radGAUSb(gausind,iaer,idomain) *          &
                      distb(j,1,iaer,idomain,gausind) +         &
                      qrefIRa(gausind,iaer) *                   &
                      pi*radGAUSa(gausind,iaer,idomain) *       &
                      radGAUSa(gausind,iaer,idomain) *          &
                      dista(j,1,iaer,idomain,gausind)           &
                      )
                    qscatrefIRgrid(j,1,iaer) =                  &
                      qscatrefIRgrid(j,1,iaer) +                &
                      weightgaus(gausind) *                     &
                      (                                         &
                      omegrefIRb(gausind,iaer) *                &
                      qrefIRb(gausind,iaer) *                   &
                      pi*radGAUSb(gausind,iaer,idomain) *       &
                      radGAUSb(gausind,iaer,idomain) *          &
                      distb(j,1,iaer,idomain,gausind) +         &
                      omegrefIRa(gausind,iaer) *                &
                      qrefIRa(gausind,iaer) *                   &
                      pi*radGAUSa(gausind,iaer,idomain) *       &
                      radGAUSa(gausind,iaer,idomain) *          &
                      dista(j,1,iaer,idomain,gausind)           &
                      )
                  ENDDO
 
                  qrefIRgrid(j,1,iaer)=qrefIRgrid(j,1,iaer) /          &
                                normd(j,1,iaer,idomain)
                  qscatrefIRgrid(j,1,iaer)=qscatrefIRgrid(j,1,iaer) /  &
                                normd(j,1,iaer,idomain)
                  omegrefIRgrid(j,1,iaer)=qscatrefIRgrid(j,1,iaer) /   &
                               qrefIRgrid(j,1,iaer)
                  DO m=1,L_NSPECTI
                    qextIRgrid(j,1,m,iaer)=qextIRgrid(j,1,m,iaer) /    &
                                normd(j,1,iaer,idomain)
                    qscatIRgrid(j,1,m,iaer)=qscatIRgrid(j,1,m,iaer) /  &
                                normd(j,1,iaer,idomain)
                    gIRgrid(j,1,m,iaer)=gIRgrid(j,1,m,iaer) /          &
                                qscatIRgrid(j,1,m,iaer) /              &
                                normd(j,1,iaer,idomain)

                    qsqrefIRgrid(j,1,m,iaer)=qextIRgrid(j,1,m,iaer) /  &
                                qrefVISgrid(j,1,iaer)
                    omegIRgrid(j,1,m,iaer)=qscatIRgrid(j,1,m,iaer) /   &
                                qextIRgrid(j,1,m,iaer)
                  ENDDO
                ENDIF                  ! --------------------------
                checkgrid(j,1,iaer,idomain) = .true.
              ENDIF !checkgrid
          ENDDO !grid_i
!         2.4 Linear interpolation
          k1 = (1-kx)
          k2 = kx
          IF (idomain.EQ.1) THEN ! VISIBLE ------------------------
          DO m=1,L_NSPECTV
             QVISsQREF3d(ig,lg,m,iaer) =                           &
                        k1*qsqrefVISgrid(grid_i,1,m,iaer) +        &
                        k2*qsqrefVISgrid(grid_i+1,1,m,iaer)
            omegaVIS3d(ig,lg,m,iaer) =                             &
                        k1*omegVISgrid(grid_i,1,m,iaer) +          &
                        k2*omegVISgrid(grid_i+1,1,m,iaer)
            gVIS3d(ig,lg,m,iaer) =                                 &
                        k1*gVISgrid(grid_i,1,m,iaer) +             &
                        k2*gVISgrid(grid_i+1,1,m,iaer)
          ENDDO !L_NSPECTV
          QREFvis3d(ig,lg,iaer) =                                  &
                        k1*qrefVISgrid(grid_i,1,iaer) +            &
                        k2*qrefVISgrid(grid_i+1,1,iaer)
!          omegaREFvis3d(ig,lg,iaer) =                              &
!                        k1*omegrefVISgrid(grid_i,1,iaer) +         &
!                        k2*omegrefVISgrid(grid_i+1,1,iaer)
          ELSE                   ! INFRARED -----------------------
          DO m=1,L_NSPECTI
            QIRsQREF3d(ig,lg,m,iaer) =                             &
                        k1*qsqrefIRgrid(grid_i,1,m,iaer) +         &
                        k2*qsqrefIRgrid(grid_i+1,1,m,iaer)
            omegaIR3d(ig,lg,m,iaer) =                              &
                        k1*omegIRgrid(grid_i,1,m,iaer) +           &
                        k2*omegIRgrid(grid_i+1,1,m,iaer) 
            gIR3d(ig,lg,m,iaer) =                                  & 
                        k1*gIRgrid(grid_i,1,m,iaer) +              &
                        k2*gIRgrid(grid_i+1,1,m,iaer)
          ENDDO !L_NSPECTI
          QREFir3d(ig,lg,iaer) =                                   &
                        k1*qrefIRgrid(grid_i,1,iaer) +             &
                        k2*qrefIRgrid(grid_i+1,1,iaer)
!          omegaREFir3d(ig,lg,iaer) =                               &
!                        k1*omegrefIRgrid(grid_i,1,iaer) +          &
!                        k2*omegrefIRgrid(grid_i+1,1,iaer)
          ENDIF                  ! --------------------------------
        ENDDO !nlayer
      ENDDO !ngrid

!==================================================================



      ENDDO ! idomain

      ENDIF ! nsize = 1

      ENDDO ! iaer (loop on aerosol kind)

      RETURN
    END SUBROUTINE aeroptproperties



      
