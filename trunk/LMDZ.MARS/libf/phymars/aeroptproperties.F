      SUBROUTINE aeroptproperties(ngrid,nlayer,reffrad,nueffrad,
     &                            QVISsQREF3d,omegaVIS3d,gVIS3d,
     &                            QIRsQREF3d,omegaIR3d,gIR3d,
     &                            QREFvis3d,QREFir3d,
     &                            omegaREFvis3d,omegaREFir3d)
      use dimradmars_mod, only: nir, nsun, naerkind,
     &                radiustab, nsize, QVISsQREF, omegavis, gvis,
     &                QIRsQREF, omegaIR, gIR, QREFvis, QREFir,
     &                omegaREFvis, omegaREFir
      IMPLICIT NONE
c     =============================================================
c     Aerosol Optical Properties
c
c     Description:
c       Compute the scattering parameters in each grid
c       box, depending on aerosol grain sizes. Log-normal size
c       distribution and Gauss-Legendre integration are used.

c     Parameters:
c       Don't forget to set the value of varyingnueff below; If
c       the effective variance of the distribution for the given
c       aerosol is considered homogeneous in the atmosphere, please
c       set varyingnueff(iaer) to .false. Resulting computational
c       time will be much better.

c     Authors: J.-B. Madeleine, F. Forget, F. Montmessin
c     =============================================================

#include "callkeys.h"

c     Local variables 
c     ---------------

c     =============================================================
      LOGICAL,SAVE,ALLOCATABLE :: varyingnueff(:)
c     =============================================================

c     Min. and max radius of the interpolation grid (in METERS)
      REAL, SAVE, ALLOCATABLE :: refftabmin(:,:)
      REAL, SAVE, ALLOCATABLE :: refftabmax(:,:)

c     Log of the min and max variance of the interpolation grid
      REAL, PARAMETER :: nuefftabmin = -4.6
      REAL, PARAMETER :: nuefftabmax = 0.
c     Number of effective radius of the interpolation grid
      INTEGER, PARAMETER :: refftabsize = 100
c     Number of effective variances of the interpolation grid
      INTEGER, PARAMETER :: nuefftabsize = 100
c     Interpolation grid indices (reff,nueff)
      INTEGER :: grid_i,grid_j
c     Intermediate variable
      REAL :: var_tmp,var3d_tmp(ngrid,nlayer)
c     Bilinear interpolation factors
      REAL :: kx,ky,k1,k2,k3,k4
c     Size distribution parameters
      REAL :: sizedistk1,sizedistk2
c     Pi!
      REAL,SAVE :: pi
c     Variables used by the Gauss-Legendre integration:
      INTEGER radius_id,gausind
      REAL kint
      REAL drad
      INTEGER, PARAMETER :: ngau = 10
      REAL weightgaus(ngau),radgaus(ngau)
      SAVE weightgaus,radgaus
c     DATA weightgaus/.2955242247,.2692667193,.2190863625,
c    &                .1494513491,.0666713443/
c     DATA radgaus/.1488743389,.4333953941,.6794095682,
c    &             .8650633666,.9739065285/
      DATA    radgaus/0.07652652113350,0.22778585114165,
     &                0.37370608871528,0.51086700195146,
     &                0.63605368072468,0.74633190646476,
     &                0.83911697181213,0.91223442826796,
     &                0.96397192726078,0.99312859919241/

      DATA weightgaus/0.15275338723120,0.14917298659407,
     &                0.14209610937519,0.13168863843930,
     &                0.11819453196154,0.10193011980823,
     &                0.08327674160932,0.06267204829828,
     &                0.04060142982019,0.01761400714091/

c     Indices
      INTEGER :: i,j,k,l,m,iaer,idomain
      INTEGER :: ig,lg,chg

c     Local saved variables
c     ---------------------

c     Radius axis of the interpolation grid
      REAL,SAVE,ALLOCATABLE :: refftab(:,:,:)
c     Variance axis of the interpolation grid
      REAL,SAVE,ALLOCATABLE :: nuefftab(:,:,:)
c     Volume ratio of the grid
      REAL,SAVE,ALLOCATABLE :: logvratgrid(:,:)
c     Grid used to remember which calculation is done
      LOGICAL,SAVE,ALLOCATABLE :: checkgrid(:,:,:,:)
c     Optical properties of the grid (VISIBLE)
      REAL,SAVE,ALLOCATABLE :: qsqrefVISgrid(:,:,:,:)
      REAL,SAVE,ALLOCATABLE :: qextVISgrid(:,:,:,:)
      REAL,SAVE,ALLOCATABLE :: qscatVISgrid(:,:,:,:)
      REAL,SAVE,ALLOCATABLE :: omegVISgrid(:,:,:,:)
      REAL,SAVE,ALLOCATABLE :: gVISgrid(:,:,:,:)
c     Optical properties of the grid (INFRARED)
      REAL,SAVE,ALLOCATABLE :: qsqrefIRgrid(:,:,:,:)
      REAL,SAVE,ALLOCATABLE :: qextIRgrid(:,:,:,:)
      REAL,SAVE,ALLOCATABLE :: qscatIRgrid(:,:,:,:)
      REAL,SAVE,ALLOCATABLE :: omegIRgrid(:,:,:,:)
      REAL,SAVE,ALLOCATABLE :: gIRgrid(:,:,:,:)
c     Optical properties of the grid (REFERENCE WAVELENGTHS)
      REAL,SAVE,ALLOCATABLE :: qrefVISgrid(:,:,:)
      REAL,SAVE,ALLOCATABLE :: qscatrefVISgrid(:,:,:)
      REAL,SAVE,ALLOCATABLE :: qrefIRgrid(:,:,:)
      REAL,SAVE,ALLOCATABLE :: qscatrefIRgrid(:,:,:)
      REAL,SAVE,ALLOCATABLE :: omegrefVISgrid(:,:,:)
      REAL,SAVE,ALLOCATABLE :: omegrefIRgrid(:,:,:)
c     Firstcall
      LOGICAL,SAVE :: firstcall = .true.
c     Variables used by the Gauss-Legendre integration:
      REAL,SAVE,ALLOCATABLE :: normd(:,:,:,:)
      REAL,SAVE,ALLOCATABLE :: dista(:,:,:,:,:)
      REAL,SAVE,ALLOCATABLE :: distb(:,:,:,:,:)
      REAL,SAVE,ALLOCATABLE :: radGAUSa(:,:,:)
      REAL,SAVE,ALLOCATABLE :: radGAUSb(:,:,:)

      REAL,SAVE,ALLOCATABLE :: qsqrefVISa(:,:,:)
      REAL,SAVE,ALLOCATABLE :: qrefVISa(:,:)
      REAL,SAVE,ALLOCATABLE :: qsqrefVISb(:,:,:)
      REAL,SAVE,ALLOCATABLE :: qrefVISb(:,:)
      REAL,SAVE,ALLOCATABLE :: omegVISa(:,:,:)
      REAL,SAVE,ALLOCATABLE :: omegrefVISa(:,:)
      REAL,SAVE,ALLOCATABLE :: omegVISb(:,:,:)
      REAL,SAVE,ALLOCATABLE :: omegrefVISb(:,:)
      REAL,SAVE,ALLOCATABLE :: gVISa(:,:,:)
      REAL,SAVE,ALLOCATABLE :: gVISb(:,:,:)

      REAL,SAVE,ALLOCATABLE :: qsqrefIRa(:,:,:)
      REAL,SAVE,ALLOCATABLE :: qrefIRa(:,:)
      REAL,SAVE,ALLOCATABLE :: qsqrefIRb(:,:,:)
      REAL,SAVE,ALLOCATABLE :: qrefIRb(:,:)
      REAL,SAVE,ALLOCATABLE :: omegIRa(:,:,:)
      REAL,SAVE,ALLOCATABLE :: omegrefIRa(:,:)
      REAL,SAVE,ALLOCATABLE :: omegIRb(:,:,:)
      REAL,SAVE,ALLOCATABLE :: omegrefIRb(:,:)
      REAL,SAVE,ALLOCATABLE :: gIRa(:,:,:)
      REAL,SAVE,ALLOCATABLE :: gIRb(:,:,:)

      REAL :: radiusm
      REAL :: radiusr

c     Inputs
c     ------

      INTEGER,INTENT(IN) :: ngrid,nlayer
c     Aerosol effective radius used for radiative transfer (meter)
      REAL,INTENT(IN) :: reffrad(ngrid,nlayer,naerkind)
c     Aerosol effective variance used for radiative transfer (n.u.)
      REAL,INTENT(IN) :: nueffrad(ngrid,nlayer,naerkind)

c     Outputs
c     -------

      REAL,INTENT(OUT) :: QVISsQREF3d(ngrid,nlayer,nsun,naerkind)
      REAL,INTENT(OUT) :: omegaVIS3d(ngrid,nlayer,nsun,naerkind)
      REAL,INTENT(OUT) :: gVIS3d(ngrid,nlayer,nsun,naerkind)

      REAL,INTENT(OUT) :: QIRsQREF3d(ngrid,nlayer,nir,naerkind)
      REAL,INTENT(OUT) :: omegaIR3d(ngrid,nlayer,nir,naerkind)
      REAL,INTENT(OUT) :: gIR3d(ngrid,nlayer,nir,naerkind)

      REAL,INTENT(OUT) :: QREFvis3d(ngrid,nlayer,naerkind)
      REAL,INTENT(OUT) :: QREFir3d(ngrid,nlayer,naerkind)

      REAL,INTENT(OUT) :: omegaREFvis3d(ngrid,nlayer,naerkind)
      REAL,INTENT(OUT) :: omegaREFir3d(ngrid,nlayer,naerkind)

c     Tests
c     -----

      LOGICAL,SAVE :: out_qwg = .false.
      INTEGER, PARAMETER :: out_iaer = 2
      INTEGER :: out_ndim
      REAL :: out_qext(ngrid,nlayer)
      REAL :: out_omeg(ngrid,nlayer)
      REAL :: out_g(ngrid,nlayer)
      INTEGER :: out_nchannel
      CHARACTER*1 :: out_str

c     Creating the effective radius and variance grid
c     -----------------------------------------------
!     AS: OK firstcall absolute
      IF (firstcall) THEN
 
c       0.0 Allocate all local saved arrays:
        allocate(refftab(refftabsize,naerkind,2))
        allocate(nuefftab(nuefftabsize,naerkind,2))
        ! Optical properties of the grid (VISIBLE)
        allocate(qsqrefVISgrid(refftabsize,nuefftabsize,nsun,naerkind))
        allocate(qextVISgrid(refftabsize,nuefftabsize,nsun,naerkind))
        allocate(qscatVISgrid(refftabsize,nuefftabsize,nsun,naerkind))
        allocate(omegVISgrid(refftabsize,nuefftabsize,nsun,naerkind))
        allocate(gVISgrid(refftabsize,nuefftabsize,nsun,naerkind))
        ! Optical properties of the grid (INFRARED)
        allocate(qsqrefIRgrid(refftabsize,nuefftabsize,nir,naerkind))
        allocate(qextIRgrid(refftabsize,nuefftabsize,nir,naerkind))
        allocate(qscatIRgrid(refftabsize,nuefftabsize,nir,naerkind))
        allocate(omegIRgrid(refftabsize,nuefftabsize,nir,naerkind))
        allocate(gIRgrid(refftabsize,nuefftabsize,nir,naerkind))
        
        allocate(qsqrefVISa(nsun,ngau,naerkind))
        allocate(qrefVISa(ngau,naerkind))
        allocate(qsqrefVISb(nsun,ngau,naerkind))
        allocate(qrefVISb(ngau,naerkind))
        allocate(omegVISa(nsun,ngau,naerkind))
        allocate(omegrefVISa(ngau,naerkind))
        allocate(omegVISb(nsun,ngau,naerkind))
        allocate(omegrefVISb(ngau,naerkind))
        allocate(gVISa(nsun,ngau,naerkind))
        allocate(gVISb(nsun,ngau,naerkind))
        
        allocate(qsqrefIRa(nir,ngau,naerkind))
        allocate(qrefIRa(ngau,naerkind))
        allocate(qsqrefIRb(nir,ngau,naerkind))
        allocate(qrefIRb(ngau,naerkind))
        allocate(omegIRa(nir,ngau,naerkind))
        allocate(omegrefIRa(ngau,naerkind))
        allocate(omegIRb(nir,ngau,naerkind))
        allocate(omegrefIRb(ngau,naerkind))
        allocate(gIRa(nir,ngau,naerkind))
        allocate(gIRb(nir,ngau,naerkind))
       
        allocate(qrefVISgrid(refftabsize,nuefftabsize,naerkind))
        allocate(qscatrefVISgrid(refftabsize,nuefftabsize,naerkind))
        allocate(qrefIRgrid(refftabsize,nuefftabsize,naerkind))
        allocate(qscatrefIRgrid(refftabsize,nuefftabsize,naerkind))
        allocate(omegrefVISgrid(refftabsize,nuefftabsize,naerkind))
        allocate(omegrefIRgrid(refftabsize,nuefftabsize,naerkind))

        allocate(normd(refftabsize,nuefftabsize,naerkind,2))
        allocate(dista(refftabsize,nuefftabsize,naerkind,2,ngau))
        allocate(distb(refftabsize,nuefftabsize,naerkind,2,ngau))
        allocate(radGAUSa(ngau,naerkind,2))
        allocate(radGAUSb(ngau,naerkind,2))

        allocate(checkgrid(refftabsize,nuefftabsize,naerkind,2))

        allocate(logvratgrid(naerkind,2))

        allocate(refftabmin(naerkind,2))
        allocate(refftabmax(naerkind,2))

        allocate(varyingnueff(naerkind))

        checkgrid(1:refftabsize,1:nuefftabsize,1:naerkind,1:2) = .false.
        varyingnueff(1:naerkind) = .false.

c       0.1 Pi!
        pi = 2. * asin(1.e0)

        WRITE(*,*) "aeroptproperties: interpolation grid"
        DO iaer = 1, naerkind ! Loop on aerosol kind
          DO idomain = 1, 2 ! Loop on visible or infrared channel

c           0.2 Effective radius
            radiusm=
     &             0.5*(radiustab(iaer,idomain,nsize(iaer,idomain))+
     &                  radiustab(iaer,idomain,1))
            radiusr=
     &             0.5*(radiustab(iaer,idomain,nsize(iaer,idomain))-
     &                  radiustab(iaer,idomain,1))
            refftabmin(iaer,idomain) = 
     &        radiusm-radiusr*radgaus(ngau)
            refftabmax(iaer,idomain) = 
     &        radiustab(iaer,idomain,nsize(iaer,idomain))

            WRITE(*,*) "Scatterer: ",iaer
            WRITE(*,*) "Domain: ",idomain
            WRITE(*,*) "Min radius (m): ", refftabmin(iaer,idomain)
            WRITE(*,*) "Max radius (m): ", refftabmax(iaer,idomain)

            refftab(1,iaer,idomain) = 
     &        refftabmin(iaer,idomain)
            refftab(refftabsize,iaer,idomain) = 
     &        refftabmax(iaer,idomain)

            logvratgrid(iaer,idomain) = 
     &                    log(refftabmax(iaer,idomain)/
     &                      refftabmin(iaer,idomain)) /
     &                    float(refftabsize-1)*3.
            do i = 2, refftabsize-1
              refftab(i,iaer,idomain) = refftab(i-1,iaer,idomain)*
     &                             exp(1./3.*logvratgrid(iaer,idomain))
            enddo

c           0.3 Effective variance
            do i = 0, nuefftabsize-1
              nuefftab(i+1,iaer,idomain) = exp( nuefftabmin +
     &                 i*(nuefftabmax-nuefftabmin)/(nuefftabsize-1) )
            enddo

          ENDDO
        ENDDO
        firstcall = .false.
      ENDIF

      DO iaer = 1, naerkind ! Loop on aerosol kind
        IF ( (nsize(iaer,1).EQ.1).AND.(nsize(iaer,2).EQ.1) ) THEN
c==================================================================
c       If there is one single particle size, optical
c         properties of the considered aerosol are homogeneous
          DO lg = 1, nlayer
            DO ig = 1, ngrid
              DO chg = 1, nsun
                QVISsQREF3d(ig,lg,chg,iaer)=QVISsQREF(chg,iaer,1)
                omegaVIS3d(ig,lg,chg,iaer)=omegaVIS(chg,iaer,1)
                gVIS3d(ig,lg,chg,iaer)=gVIS(chg,iaer,1)
              ENDDO
              DO chg = 1, nir
                QIRsQREF3d(ig,lg,chg,iaer)=QIRsQREF(chg,iaer,1)
                omegaIR3d(ig,lg,chg,iaer)=omegaIR(chg,iaer,1)
                gIR3d(ig,lg,chg,iaer)=gIR(chg,iaer,1)
              ENDDO
              QREFvis3d(ig,lg,iaer)=QREFvis(iaer,1)
              QREFir3d(ig,lg,iaer)=QREFir(iaer,1)
              omegaREFvis3d(ig,lg,iaer)=omegaREFvis(iaer,1)
              omegaREFir3d(ig,lg,iaer)=omegaREFir(iaer,1)
            ENDDO
          ENDDO
        ELSE ! Varying effective radius and variance
      DO idomain = 1, 2 ! Loop on visible or infrared channel
c==================================================================

c       1.1 Radius middle point and range for Gauss integration
        radiusm=
     &         0.5*(radiustab(iaer,idomain,nsize(iaer,idomain)) +
     &              radiustab(iaer,idomain,1))
        radiusr=
     &         0.5*(radiustab(iaer,idomain,nsize(iaer,idomain)) -
     &              radiustab(iaer,idomain,1))

c       1.2 Interpolating data at the Gauss quadrature points:
        DO gausind=1,ngau
          drad=radiusr*radgaus(gausind)
          radGAUSa(gausind,iaer,idomain)=radiusm-drad

          radius_id=minloc(abs(radiustab(iaer,idomain,:) -
     &                         (radiusm-drad)),DIM=1)
          IF ((radiustab(iaer,idomain,radius_id) - 
     &         (radiusm-drad)).GT.0) THEN
            radius_id=radius_id-1
          ENDIF
          IF (radius_id.GE.nsize(iaer,idomain)) THEN
            radius_id=nsize(iaer,idomain)-1
            kint = 1.
          ELSEIF (radius_id.LT.1) THEN
            radius_id=1
            kint = 0.
          ELSE
          kint = ( (radiusm-drad) - 
     &             radiustab(iaer,idomain,radius_id) ) /
     &           ( radiustab(iaer,idomain,radius_id+1) - 
     &             radiustab(iaer,idomain,radius_id) )
          ENDIF
          IF (idomain.EQ.1) THEN ! VISIBLE DOMAIN -----------------
            DO m=1,nsun
            qsqrefVISa(m,gausind,iaer)=
     &              (1-kint)*QVISsQREF(m,iaer,radius_id) + 
     &              kint*QVISsQREF(m,iaer,radius_id+1)
            omegVISa(m,gausind,iaer)=
     &              (1-kint)*omegaVIS(m,iaer,radius_id) + 
     &              kint*omegaVIS(m,iaer,radius_id+1)
            gVISa(m,gausind,iaer)=
     &              (1-kint)*gVIS(m,iaer,radius_id) + 
     &              kint*gVIS(m,iaer,radius_id+1)
            ENDDO
            qrefVISa(gausind,iaer)=
     &              (1-kint)*QREFvis(iaer,radius_id) + 
     &              kint*QREFvis(iaer,radius_id+1)
            omegrefVISa(gausind,iaer)=
     &              (1-kint)*omegaREFvis(iaer,radius_id) + 
     &              kint*omegaREFvis(iaer,radius_id+1)
          ELSE ! INFRARED DOMAIN ----------------------------------
            DO m=1,nir
            qsqrefIRa(m,gausind,iaer)=
     &              (1-kint)*QIRsQREF(m,iaer,radius_id) + 
     &              kint*QIRsQREF(m,iaer,radius_id+1)
            omegIRa(m,gausind,iaer)=
     &              (1-kint)*omegaIR(m,iaer,radius_id) + 
     &              kint*omegaIR(m,iaer,radius_id+1)
            gIRa(m,gausind,iaer)=
     &              (1-kint)*gIR(m,iaer,radius_id) + 
     &              kint*gIR(m,iaer,radius_id+1)
            ENDDO
            qrefIRa(gausind,iaer)=
     &              (1-kint)*QREFir(iaer,radius_id) + 
     &              kint*QREFir(iaer,radius_id+1)
            omegrefIRa(gausind,iaer)=
     &              (1-kint)*omegaREFir(iaer,radius_id) + 
     &              kint*omegaREFir(iaer,radius_id+1)
          ENDIF
        ENDDO

        DO gausind=1,ngau
          drad=radiusr*radgaus(gausind)
          radGAUSb(gausind,iaer,idomain)=radiusm+drad

          radius_id=minloc(abs(radiustab(iaer,idomain,:) - 
     &                         (radiusm+drad)),DIM=1)
          IF ((radiustab(iaer,idomain,radius_id) - 
     &         (radiusm+drad)).GT.0) THEN
            radius_id=radius_id-1
          ENDIF
          IF (radius_id.GE.nsize(iaer,idomain)) THEN
            radius_id=nsize(iaer,idomain)-1
            kint = 1.
          ELSEIF (radius_id.LT.1) THEN
            radius_id=1
            kint = 0.
          ELSE
            kint = ( (radiusm+drad) -
     &               radiustab(iaer,idomain,radius_id) ) /
     &             ( radiustab(iaer,idomain,radius_id+1) - 
     &               radiustab(iaer,idomain,radius_id) )
          ENDIF
          IF (idomain.EQ.1) THEN ! VISIBLE DOMAIN -----------------
            DO m=1,nsun
            qsqrefVISb(m,gausind,iaer)=
     &              (1-kint)*QVISsQREF(m,iaer,radius_id) + 
     &              kint*QVISsQREF(m,iaer,radius_id+1)
            omegVISb(m,gausind,iaer)=
     &              (1-kint)*omegaVIS(m,iaer,radius_id) + 
     &              kint*omegaVIS(m,iaer,radius_id+1)
            gVISb(m,gausind,iaer)=
     &              (1-kint)*gVIS(m,iaer,radius_id) + 
     &              kint*gVIS(m,iaer,radius_id+1)
            ENDDO
            qrefVISb(gausind,iaer)=
     &              (1-kint)*QREFvis(iaer,radius_id) + 
     &              kint*QREFvis(iaer,radius_id+1)
            omegrefVISb(gausind,iaer)=
     &              (1-kint)*omegaREFvis(iaer,radius_id) + 
     &              kint*omegaREFvis(iaer,radius_id+1)
          ELSE ! INFRARED DOMAIN ----------------------------------
            DO m=1,nir
            qsqrefIRb(m,gausind,iaer)=
     &              (1-kint)*QIRsQREF(m,iaer,radius_id) + 
     &              kint*QIRsQREF(m,iaer,radius_id+1)
            omegIRb(m,gausind,iaer)=
     &              (1-kint)*omegaIR(m,iaer,radius_id) + 
     &              kint*omegaIR(m,iaer,radius_id+1)
            gIRb(m,gausind,iaer)=
     &              (1-kint)*gIR(m,iaer,radius_id) + 
     &              kint*gIR(m,iaer,radius_id+1)
            ENDDO
            qrefIRb(gausind,iaer)=
     &              (1-kint)*QREFir(iaer,radius_id) + 
     &              kint*QREFir(iaer,radius_id+1)
            omegrefIRb(gausind,iaer)=
     &              (1-kint)*omegaREFir(iaer,radius_id) + 
     &              kint*omegaREFir(iaer,radius_id+1)
          ENDIF
        ENDDO
c==================================================================
      IF ( .NOT.varyingnueff(iaer) ) THEN          ! CONSTANT NUEFF
c==================================================================
c     2. Compute the scattering parameters using linear
c       interpolation over grain sizes and constant nueff
c     ---------------------------------------------------

      DO lg = 1,nlayer
        DO ig = 1,ngrid
c         2.1 Effective radius index and kx calculation
          var_tmp=reffrad(ig,lg,iaer)/refftabmin(iaer,idomain)
          var_tmp=log(var_tmp)*3.
          var_tmp=var_tmp/logvratgrid(iaer,idomain)+1.
          grid_i=floor(var_tmp)
          IF (grid_i.GE.refftabsize) THEN
c           WRITE(*,*) 'Warning: particle size in grid box #'
c           WRITE(*,*) ig,' is too large to be used by the '
c           WRITE(*,*) 'radiative transfer; please extend the '
c           WRITE(*,*) 'interpolation grid to larger grain sizes.'
            grid_i=refftabsize-1
            kx = 1.
          ELSEIF (grid_i.LT.1) THEN
c           WRITE(*,*) 'Warning: particle size in grid box #'
c           WRITE(*,*) ig,' is too small to be used by the '
c           WRITE(*,*) 'radiative transfer; please extend the '
c           WRITE(*,*) 'interpolation grid to smaller grain sizes.'
            grid_i=1
            kx = 0.
          ELSE
            kx = ( reffrad(ig,lg,iaer)-
     &             refftab(grid_i,iaer,idomain) ) /
     &           ( refftab(grid_i+1,iaer,idomain)-
     &             refftab(grid_i,iaer,idomain) )
          ENDIF
c         2.3 Integration
          DO j=grid_i,grid_i+1
c             2.3.1 Check if the calculation has been done
              IF (.NOT.checkgrid(j,1,iaer,idomain)) THEN
c               2.3.2 Log-normal dist., r_g and sigma_g are defined
c                 in [hansen_1974], "Light scattering in planetary
c                 atmospheres", Space Science Reviews 16 527-610.
c                 Here, sizedistk1=r_g and sizedistk2=sigma_g^2
                sizedistk2 = log(1.+nueffrad(1,1,iaer))
                sizedistk1 = exp(2.5*sizedistk2)
                sizedistk1 = refftab(j,iaer,idomain) / sizedistk1

                normd(j,1,iaer,idomain) = 1e-30
                DO gausind=1,ngau
                  drad=radiusr*radgaus(gausind)
                  dista(j,1,iaer,idomain,gausind) =
     &              LOG((radiusm-drad)/sizedistk1)
                  dista(j,1,iaer,idomain,gausind) =
     &              EXP(-dista(j,1,iaer,idomain,gausind) * 
     &              dista(j,1,iaer,idomain,gausind) *
     &              0.5e0/sizedistk2)/(radiusm-drad)
                  dista(j,1,iaer,idomain,gausind) =
     &              dista(j,1,iaer,idomain,gausind) / 
     &              (sqrt(2e0*pi*sizedistk2))

                  distb(j,1,iaer,idomain,gausind) = 
     &              LOG((radiusm+drad)/sizedistk1)
                  distb(j,1,iaer,idomain,gausind) = 
     &              EXP(-distb(j,1,iaer,idomain,gausind) * 
     &              distb(j,1,iaer,idomain,gausind) *
     &              0.5e0/sizedistk2)/(radiusm+drad)
                  distb(j,1,iaer,idomain,gausind) = 
     &              distb(j,1,iaer,idomain,gausind) / 
     &              (sqrt(2e0*pi*sizedistk2))

                  normd(j,1,iaer,idomain)=normd(j,1,iaer,idomain) +
     &              weightgaus(gausind) *
     &              (
     &              distb(j,1,iaer,idomain,gausind) * pi * 
     &              radGAUSb(gausind,iaer,idomain) * 
     &              radGAUSb(gausind,iaer,idomain) +
     &              dista(j,1,iaer,idomain,gausind) * pi * 
     &              radGAUSa(gausind,iaer,idomain) * 
     &              radGAUSa(gausind,iaer,idomain)
     &              )
                ENDDO
                IF (normd(j,1,iaer,idomain).EQ.1e-30) THEN
                  WRITE(*,*)"normd:", normd(j,1,iaer,idomain)
                  WRITE(*,*)"Risk of division by 0 (aeroptproperties.F)"
                  WRITE(*,*)"Check the size of the interpolation grid."
                  CALL ABORT
                ENDIF
                IF (idomain.EQ.1) THEN ! VISIBLE DOMAIN -----------
c                 2.3.3.vis Initialization
                  qsqrefVISgrid(j,1,:,iaer)=0.
                  qextVISgrid(j,1,:,iaer)=0.
                  qscatVISgrid(j,1,:,iaer)=0.
                  omegVISgrid(j,1,:,iaer)=0.
                  gVISgrid(j,1,:,iaer)=0.
                  qrefVISgrid(j,1,iaer)=0.
                  qscatrefVISgrid(j,1,iaer)=0.
                  omegrefVISgrid(j,1,iaer)=0.

                  DO gausind=1,ngau
                    DO m=1,nsun
c                     Convolution:
                      qextVISgrid(j,1,m,iaer) = 
     &                  qextVISgrid(j,1,m,iaer) + 
     &                  weightgaus(gausind) *
     &                  (
     &                  qsqrefVISb(m,gausind,iaer) * 
     &                  qrefVISb(gausind,iaer) *
     &                  pi*radGAUSb(gausind,iaer,idomain) * 
     &                  radGAUSb(gausind,iaer,idomain) * 
     &                  distb(j,1,iaer,idomain,gausind) +
     &                  qsqrefVISa(m,gausind,iaer) * 
     &                  qrefVISa(gausind,iaer) *
     &                  pi*radGAUSa(gausind,iaer,idomain) * 
     &                  radGAUSa(gausind,iaer,idomain) * 
     &                  dista(j,1,iaer,idomain,gausind)
     &                  )
                      qscatVISgrid(j,1,m,iaer) = 
     &                  qscatVISgrid(j,1,m,iaer) + 
     &                  weightgaus(gausind) *
     &                  (
     &                  omegVISb(m,gausind,iaer) *
     &                  qsqrefVISb(m,gausind,iaer) *
     &                  qrefVISb(gausind,iaer) *
     &                  pi*radGAUSb(gausind,iaer,idomain) *
     &                  radGAUSb(gausind,iaer,idomain) *
     &                  distb(j,1,iaer,idomain,gausind) +
     &                  omegVISa(m,gausind,iaer) *
     &                  qsqrefVISa(m,gausind,iaer) *
     &                  qrefVISa(gausind,iaer) *
     &                  pi*radGAUSa(gausind,iaer,idomain) *
     &                  radGAUSa(gausind,iaer,idomain) *
     &                  dista(j,1,iaer,idomain,gausind)
     &                  )
                      gVISgrid(j,1,m,iaer) =
     &                  gVISgrid(j,1,m,iaer) +
     &                  weightgaus(gausind) *
     &                  (
     &                  omegVISb(m,gausind,iaer) *
     &                  qsqrefVISb(m,gausind,iaer) *
     &                  qrefVISb(gausind,iaer) *
     &                  gVISb(m,gausind,iaer) *
     &                  pi*radGAUSb(gausind,iaer,idomain) *
     &                  radGAUSb(gausind,iaer,idomain) *
     &                  distb(j,1,iaer,idomain,gausind) +
     &                  omegVISa(m,gausind,iaer) *
     &                  qsqrefVISa(m,gausind,iaer) *
     &                  qrefVISa(gausind,iaer) *
     &                  gVISa(m,gausind,iaer) *
     &                  pi*radGAUSa(gausind,iaer,idomain) *
     &                  radGAUSa(gausind,iaer,idomain) *
     &                  dista(j,1,iaer,idomain,gausind)
     &                  )
                    ENDDO
                    qrefVISgrid(j,1,iaer) = 
     &                qrefVISgrid(j,1,iaer) + 
     &                weightgaus(gausind) *
     &                (
     &                qrefVISb(gausind,iaer) *
     &                pi*radGAUSb(gausind,iaer,idomain) * 
     &                radGAUSb(gausind,iaer,idomain) * 
     &                distb(j,1,iaer,idomain,gausind) +
     &                qrefVISa(gausind,iaer) *
     &                pi*radGAUSa(gausind,iaer,idomain) * 
     &                radGAUSa(gausind,iaer,idomain) * 
     &                dista(j,1,iaer,idomain,gausind)
     &                )
                    qscatrefVISgrid(j,1,iaer) = 
     &                qscatrefVISgrid(j,1,iaer) + 
     &                weightgaus(gausind) *
     &                (
     &                omegrefVISb(gausind,iaer) *
     &                qrefVISb(gausind,iaer) *
     &                pi*radGAUSb(gausind,iaer,idomain) * 
     &                radGAUSb(gausind,iaer,idomain) * 
     &                distb(j,1,iaer,idomain,gausind) +
     &                omegrefVISa(gausind,iaer) *
     &                qrefVISa(gausind,iaer) *
     &                pi*radGAUSa(gausind,iaer,idomain) * 
     &                radGAUSa(gausind,iaer,idomain) * 
     &                dista(j,1,iaer,idomain,gausind)
     &                )
                  ENDDO

                  qrefVISgrid(j,1,iaer)=qrefVISgrid(j,1,iaer) /
     &                          normd(j,1,iaer,idomain)
                  qscatrefVISgrid(j,1,iaer)=qscatrefVISgrid(j,1,iaer) /
     &                          normd(j,1,iaer,idomain)
                  omegrefVISgrid(j,1,iaer)=qscatrefVISgrid(j,1,iaer) /
     &                         qrefVISgrid(j,1,iaer)
                  DO m=1,nsun
                    qextVISgrid(j,1,m,iaer)=qextVISgrid(j,1,m,iaer) /
     &                          normd(j,1,iaer,idomain)
                    qscatVISgrid(j,1,m,iaer)=qscatVISgrid(j,1,m,iaer) /
     &                          normd(j,1,iaer,idomain)
                    gVISgrid(j,1,m,iaer)=gVISgrid(j,1,m,iaer) / 
     &                          qscatVISgrid(j,1,m,iaer) /
     &                          normd(j,1,iaer,idomain)

                    qsqrefVISgrid(j,1,m,iaer)=qextVISgrid(j,1,m,iaer) /
     &                          qrefVISgrid(j,1,iaer)
                    omegVISgrid(j,1,m,iaer)=qscatVISgrid(j,1,m,iaer) / 
     &                          qextVISgrid(j,1,m,iaer)
                  ENDDO
                ELSE                   ! INFRARED DOMAIN ----------
c                 2.3.3.ir Initialization
                  qsqrefIRgrid(j,1,:,iaer)=0.
                  qextIRgrid(j,1,:,iaer)=0.
                  qscatIRgrid(j,1,:,iaer)=0.
                  omegIRgrid(j,1,:,iaer)=0.
                  gIRgrid(j,1,:,iaer)=0.
                  qrefIRgrid(j,1,iaer)=0.
                  qscatrefIRgrid(j,1,iaer)=0.
                  omegrefIRgrid(j,1,iaer)=0.

                  DO gausind=1,ngau
                    DO m=1,nir
c                     Convolution:
                      qextIRgrid(j,1,m,iaer) = 
     &                  qextIRgrid(j,1,m,iaer) + 
     &                  weightgaus(gausind) *
     &                  (
     &                  qsqrefIRb(m,gausind,iaer) * 
     &                  qrefVISb(gausind,iaer) *
     &                  pi*radGAUSb(gausind,iaer,idomain) * 
     &                  radGAUSb(gausind,iaer,idomain) * 
     &                  distb(j,1,iaer,idomain,gausind) +
     &                  qsqrefIRa(m,gausind,iaer) * 
     &                  qrefVISa(gausind,iaer) *
     &                  pi*radGAUSa(gausind,iaer,idomain) * 
     &                  radGAUSa(gausind,iaer,idomain) * 
     &                  dista(j,1,iaer,idomain,gausind)
     &                  )
                      qscatIRgrid(j,1,m,iaer) = 
     &                  qscatIRgrid(j,1,m,iaer) + 
     &                  weightgaus(gausind) *
     &                  (
     &                  omegIRb(m,gausind,iaer) *
     &                  qsqrefIRb(m,gausind,iaer) *
     &                  qrefVISb(gausind,iaer) *
     &                  pi*radGAUSb(gausind,iaer,idomain) *
     &                  radGAUSb(gausind,iaer,idomain) *
     &                  distb(j,1,iaer,idomain,gausind) +
     &                  omegIRa(m,gausind,iaer) *
     &                  qsqrefIRa(m,gausind,iaer) *
     &                  qrefVISa(gausind,iaer) *
     &                  pi*radGAUSa(gausind,iaer,idomain) *
     &                  radGAUSa(gausind,iaer,idomain) *
     &                  dista(j,1,iaer,idomain,gausind)
     &                  )
                      gIRgrid(j,1,m,iaer) =
     &                  gIRgrid(j,1,m,iaer) +
     &                  weightgaus(gausind) *
     &                  (
     &                  omegIRb(m,gausind,iaer) *
     &                  qsqrefIRb(m,gausind,iaer) *
     &                  qrefVISb(gausind,iaer) *
     &                  gIRb(m,gausind,iaer) *
     &                  pi*radGAUSb(gausind,iaer,idomain) *
     &                  radGAUSb(gausind,iaer,idomain) *
     &                  distb(j,1,iaer,idomain,gausind) +
     &                  omegIRa(m,gausind,iaer) *
     &                  qsqrefIRa(m,gausind,iaer) *
     &                  qrefVISa(gausind,iaer) *
     &                  gIRa(m,gausind,iaer) *
     &                  pi*radGAUSa(gausind,iaer,idomain) *
     &                  radGAUSa(gausind,iaer,idomain) *
     &                  dista(j,1,iaer,idomain,gausind)
     &                  )
                    ENDDO
                    qrefIRgrid(j,1,iaer) = 
     &                qrefIRgrid(j,1,iaer) + 
     &                weightgaus(gausind) *
     &                (
     &                qrefIRb(gausind,iaer) *
     &                pi*radGAUSb(gausind,iaer,idomain) * 
     &                radGAUSb(gausind,iaer,idomain) * 
     &                distb(j,1,iaer,idomain,gausind) +
     &                qrefIRa(gausind,iaer) *
     &                pi*radGAUSa(gausind,iaer,idomain) * 
     &                radGAUSa(gausind,iaer,idomain) * 
     &                dista(j,1,iaer,idomain,gausind)
     &                )
                    qscatrefIRgrid(j,1,iaer) = 
     &                qscatrefIRgrid(j,1,iaer) + 
     &                weightgaus(gausind) *
     &                (
     &                omegrefIRb(gausind,iaer) *
     &                qrefIRb(gausind,iaer) *
     &                pi*radGAUSb(gausind,iaer,idomain) * 
     &                radGAUSb(gausind,iaer,idomain) * 
     &                distb(j,1,iaer,idomain,gausind) +
     &                omegrefIRa(gausind,iaer) *
     &                qrefIRa(gausind,iaer) *
     &                pi*radGAUSa(gausind,iaer,idomain) * 
     &                radGAUSa(gausind,iaer,idomain) * 
     &                dista(j,1,iaer,idomain,gausind)
     &                )
                  ENDDO

                  qrefIRgrid(j,1,iaer)=qrefIRgrid(j,1,iaer) /
     &                          normd(j,1,iaer,idomain)
                  qscatrefIRgrid(j,1,iaer)=qscatrefIRgrid(j,1,iaer) /
     &                          normd(j,1,iaer,idomain)
                  omegrefIRgrid(j,1,iaer)=qscatrefIRgrid(j,1,iaer) /
     &                         qrefIRgrid(j,1,iaer)
                  DO m=1,nir
                    qextIRgrid(j,1,m,iaer)=qextIRgrid(j,1,m,iaer) /
     &                          normd(j,1,iaer,idomain)
                    qscatIRgrid(j,1,m,iaer)=qscatIRgrid(j,1,m,iaer) /
     &                          normd(j,1,iaer,idomain)
                    gIRgrid(j,1,m,iaer)=gIRgrid(j,1,m,iaer) / 
     &                          qscatIRgrid(j,1,m,iaer) /
     &                          normd(j,1,iaer,idomain)

                    qsqrefIRgrid(j,1,m,iaer)=qextIRgrid(j,1,m,iaer) /
     &                          qrefVISgrid(j,1,iaer)
                    omegIRgrid(j,1,m,iaer)=qscatIRgrid(j,1,m,iaer) / 
     &                          qextIRgrid(j,1,m,iaer)
                  ENDDO
                ENDIF                  ! --------------------------
                checkgrid(j,1,iaer,idomain) = .true.
              ENDIF !checkgrid
          ENDDO !grid_i
c         2.4 Linear interpolation
          k1 = (1-kx)
          k2 = kx
          IF (idomain.EQ.1) THEN ! VISIBLE ------------------------
          DO m=1,nsun
            QVISsQREF3d(ig,lg,m,iaer) =
     &                  k1*qsqrefVISgrid(grid_i,1,m,iaer) +
     &                  k2*qsqrefVISgrid(grid_i+1,1,m,iaer)
            omegaVIS3d(ig,lg,m,iaer) =
     &                  k1*omegVISgrid(grid_i,1,m,iaer) +
     &                  k2*omegVISgrid(grid_i+1,1,m,iaer)
            gVIS3d(ig,lg,m,iaer) =
     &                  k1*gVISgrid(grid_i,1,m,iaer) +
     &                  k2*gVISgrid(grid_i+1,1,m,iaer)
          ENDDO !nsun
          QREFvis3d(ig,lg,iaer) =
     &                  k1*qrefVISgrid(grid_i,1,iaer) +
     &                  k2*qrefVISgrid(grid_i+1,1,iaer)
          omegaREFvis3d(ig,lg,iaer) =
     &                  k1*omegrefVISgrid(grid_i,1,iaer) +
     &                  k2*omegrefVISgrid(grid_i+1,1,iaer)
          ELSE                   ! INFRARED -----------------------
          DO m=1,nir
            QIRsQREF3d(ig,lg,m,iaer) =
     &                  k1*qsqrefIRgrid(grid_i,1,m,iaer) +
     &                  k2*qsqrefIRgrid(grid_i+1,1,m,iaer)
            omegaIR3d(ig,lg,m,iaer) =
     &                  k1*omegIRgrid(grid_i,1,m,iaer) +
     &                  k2*omegIRgrid(grid_i+1,1,m,iaer)
            gIR3d(ig,lg,m,iaer) =
     &                  k1*gIRgrid(grid_i,1,m,iaer) +
     &                  k2*gIRgrid(grid_i+1,1,m,iaer)
          ENDDO !nir
          QREFir3d(ig,lg,iaer) =
     &                  k1*qrefIRgrid(grid_i,1,iaer) +
     &                  k2*qrefIRgrid(grid_i+1,1,iaer)
          omegaREFir3d(ig,lg,iaer) =
     &                  k1*omegrefIRgrid(grid_i,1,iaer) +
     &                  k2*omegrefIRgrid(grid_i+1,1,iaer)
          ENDIF                  ! --------------------------------
        ENDDO !nlayer
      ENDDO !ngrid
c==================================================================
      ELSE                                          ! VARYING NUEFF
c==================================================================
c     3. Compute the scattering parameters in each grid box
c       using bilinear interpolation over the grain sizes
c       and the effective variances;
c     -----------------------------------------------------

      DO lg = 1,nlayer
        DO ig = 1,ngrid
c         3.1 Effective variance index and ky calculation
          var_tmp=log(nueffrad(ig,lg,iaer))
          grid_j=floor( (nuefftabsize-1)/(nuefftabmax-nuefftabmin)*
     &       (var_tmp-nuefftabmin)+1. )
          IF (grid_j.GE.nuefftabsize) THEN
c           WRITE(*,*) 'Warning: effective variance '
c           WRITE(*,*) 'is too large to be used by the '
c           WRITE(*,*) 'radiative transfer; please extend the '
c           WRITE(*,*) 'interpolation grid to larger values.'
            grid_j=nuefftabsize-1
            ky = 1.
          ELSEIF (grid_j.LT.1) THEN
c           WRITE(*,*) 'Warning: effective variance '
c           WRITE(*,*) 'is too small to be used by the '
c           WRITE(*,*) 'radiative transfer; please extend the '
c           WRITE(*,*) 'interpolation grid to smaller values.'
            grid_j=1
            ky = 0.
          ELSE
            ky = ( nueffrad(ig,lg,iaer)-
     &             nuefftab(grid_j,iaer,idomain) ) /
     &           ( nuefftab(grid_j+1,iaer,idomain)-
     &             nuefftab(grid_j,iaer,idomain) )
          ENDIF
c         3.2 Effective radius index and kx calculation
          var_tmp=reffrad(ig,lg,iaer)/refftabmin(iaer,idomain)
          var_tmp=log(var_tmp)*3.
          var_tmp=var_tmp/logvratgrid(iaer,idomain)+1.
          grid_i=floor(var_tmp)
          IF (grid_i.GE.refftabsize) THEN
c           WRITE(*,*) 'Warning: particle size in grid box #'
c           WRITE(*,*) ig,' is too large to be used by the '
c           WRITE(*,*) 'radiative transfer; please extend the '
c           WRITE(*,*) 'interpolation grid to larger grain sizes.'
            grid_i=refftabsize-1
            kx = 1.
          ELSEIF (grid_i.LT.1) THEN
c           WRITE(*,*) 'Warning: particle size in grid box #'
c           WRITE(*,*) ig,' is too small to be used by the '
c           WRITE(*,*) 'radiative transfer; please extend the '
c           WRITE(*,*) 'interpolation grid to smaller grain sizes.'
            grid_i=1
            kx = 0.
          ELSE
            kx = ( reffrad(ig,lg,iaer)-
     &             refftab(grid_i,iaer,idomain) ) /
     &           ( refftab(grid_i+1,iaer,idomain)-
     &             refftab(grid_i,iaer,idomain) )
          ENDIF
c         3.3 Integration
          DO j=grid_i,grid_i+1
            DO k=grid_j,grid_j+1
c             3.3.1 Check if the calculation has been done
              IF (.NOT.checkgrid(j,k,iaer,idomain)) THEN

c               3.3.2 Log-normal dist., r_g and sigma_g are defined
c                 in [hansen_1974], "Light scattering in planetary
c                 atmospheres", Space Science Reviews 16 527-610.
c                 Here, sizedistk1=r_g and sizedistk2=sigma_g^2
                sizedistk2 = log(1.+nuefftab(k,iaer,idomain))
                sizedistk1 = exp(2.5*sizedistk2)
                sizedistk1 = refftab(j,iaer,idomain) / sizedistk1

                normd(j,k,iaer,idomain) = 1e-30
                DO gausind=1,ngau
                  drad=radiusr*radgaus(gausind)

                  dista(j,k,iaer,idomain,gausind) =
     &              LOG((radiusm-drad)/sizedistk1)
                  dista(j,k,iaer,idomain,gausind) =
     &              EXP(-dista(j,k,iaer,idomain,gausind) *
     &              dista(j,k,iaer,idomain,gausind) *
     &              0.5e0/sizedistk2)/(radiusm-drad)
                  dista(j,k,iaer,idomain,gausind) =
     &              dista(j,k,iaer,idomain,gausind) /
     &              (sqrt(2e0*pi*sizedistk2))

                  distb(j,k,iaer,idomain,gausind) =
     &              LOG((radiusm+drad)/sizedistk1)
                  distb(j,k,iaer,idomain,gausind) =
     &              EXP(-distb(j,k,iaer,idomain,gausind) *
     &              distb(j,k,iaer,idomain,gausind) *
     &              0.5e0/sizedistk2)/(radiusm+drad)
                  distb(j,k,iaer,idomain,gausind) =
     &              distb(j,k,iaer,idomain,gausind) /
     &              (sqrt(2e0*pi*sizedistk2))

                  normd(j,k,iaer,idomain)=normd(j,k,iaer,idomain) +
     &              weightgaus(gausind) *
     &              (
     &              distb(j,k,iaer,idomain,gausind) * pi *
     &              radGAUSb(gausind,iaer,idomain) *
     &              radGAUSb(gausind,iaer,idomain) +
     &              dista(j,k,iaer,idomain,gausind) * pi *
     &              radGAUSa(gausind,iaer,idomain) *
     &              radGAUSa(gausind,iaer,idomain)
     &              )
                ENDDO
                IF (normd(j,k,iaer,idomain).EQ.1e-30) THEN
                  WRITE(*,*)"normd:", normd(j,k,iaer,idomain)
                  WRITE(*,*)"Risk of division by 0 (aeroptproperties.F)"
                  WRITE(*,*)"Check the size of the interpolation grid."
                  CALL ABORT
                ENDIF
                IF (idomain.EQ.1) THEN ! VISIBLE DOMAIN -----------
c                 2.3.3.vis Initialization
                  qsqrefVISgrid(j,k,:,iaer)=0.
                  qextVISgrid(j,k,:,iaer)=0.
                  qscatVISgrid(j,k,:,iaer)=0.
                  omegVISgrid(j,k,:,iaer)=0.
                  gVISgrid(j,k,:,iaer)=0.
                  qrefVISgrid(j,k,iaer)=0.
                  qscatrefVISgrid(j,k,iaer)=0.
                  omegrefVISgrid(j,k,iaer)=0.

                  DO gausind=1,ngau
                    DO m=1,nsun
c                     Convolution:
                      qextVISgrid(j,k,m,iaer) =
     &                  qextVISgrid(j,k,m,iaer) +
     &                  weightgaus(gausind) *
     &                  (
     &                  qsqrefVISb(m,gausind,iaer) *
     &                  qrefVISb(gausind,iaer) *
     &                  pi*radGAUSb(gausind,iaer,idomain) *
     &                  radGAUSb(gausind,iaer,idomain) *
     &                  distb(j,k,iaer,idomain,gausind) +
     &                  qsqrefVISa(m,gausind,iaer) *
     &                  qrefVISa(gausind,iaer) *
     &                  pi*radGAUSa(gausind,iaer,idomain) *
     &                  radGAUSa(gausind,iaer,idomain) *
     &                  dista(j,k,iaer,idomain,gausind)
     &                  )
                      qscatVISgrid(j,k,m,iaer) =
     &                  qscatVISgrid(j,k,m,iaer) +
     &                  weightgaus(gausind) *
     &                  (
     &                  omegVISb(m,gausind,iaer) *
     &                  qsqrefVISb(m,gausind,iaer) *
     &                  qrefVISb(gausind,iaer) *
     &                  pi*radGAUSb(gausind,iaer,idomain) *
     &                  radGAUSb(gausind,iaer,idomain) *
     &                  distb(j,k,iaer,idomain,gausind) +
     &                  omegVISa(m,gausind,iaer) *
     &                  qsqrefVISa(m,gausind,iaer) *
     &                  qrefVISa(gausind,iaer) *
     &                  pi*radGAUSa(gausind,iaer,idomain) *
     &                  radGAUSa(gausind,iaer,idomain) *
     &                  dista(j,k,iaer,idomain,gausind)
     &                  )
                      gVISgrid(j,k,m,iaer) =
     &                  gVISgrid(j,k,m,iaer) +
     &                  weightgaus(gausind) *
     &                  (
     &                  omegVISb(m,gausind,iaer) *
     &                  qsqrefVISb(m,gausind,iaer) *
     &                  qrefVISb(gausind,iaer) *
     &                  gVISb(m,gausind,iaer) *
     &                  pi*radGAUSb(gausind,iaer,idomain) *
     &                  radGAUSb(gausind,iaer,idomain) *
     &                  distb(j,k,iaer,idomain,gausind) +
     &                  omegVISa(m,gausind,iaer) *
     &                  qsqrefVISa(m,gausind,iaer) *
     &                  qrefVISa(gausind,iaer) *
     &                  gVISa(m,gausind,iaer) *
     &                  pi*radGAUSa(gausind,iaer,idomain) *
     &                  radGAUSa(gausind,iaer,idomain) *
     &                  dista(j,k,iaer,idomain,gausind)
     &                  )
                    ENDDO
                    qrefVISgrid(j,k,iaer) =
     &                qrefVISgrid(j,k,iaer) +
     &                weightgaus(gausind) *
     &                (
     &                qrefVISb(gausind,iaer) *
     &                pi*radGAUSb(gausind,iaer,idomain) *
     &                radGAUSb(gausind,iaer,idomain) *
     &                distb(j,k,iaer,idomain,gausind) +
     &                qrefVISa(gausind,iaer) *
     &                pi*radGAUSa(gausind,iaer,idomain) *
     &                radGAUSa(gausind,iaer,idomain) *
     &                dista(j,k,iaer,idomain,gausind)
     &                )
                    qscatrefVISgrid(j,k,iaer) =
     &                qscatrefVISgrid(j,k,iaer) +
     &                weightgaus(gausind) *
     &                (
     &                omegrefVISb(gausind,iaer) *
     &                qrefVISb(gausind,iaer) *
     &                pi*radGAUSb(gausind,iaer,idomain) *
     &                radGAUSb(gausind,iaer,idomain) *
     &                distb(j,k,iaer,idomain,gausind) +
     &                omegrefVISa(gausind,iaer) *
     &                qrefVISa(gausind,iaer) *
     &                pi*radGAUSa(gausind,iaer,idomain) *
     &                radGAUSa(gausind,iaer,idomain) *
     &                dista(j,k,iaer,idomain,gausind)
     &                )
                  ENDDO
                  qrefVISgrid(j,k,iaer)=qrefVISgrid(j,k,iaer) /
     &                          normd(j,k,iaer,idomain)
                  qscatrefVISgrid(j,k,iaer)=qscatrefVISgrid(j,k,iaer) /
     &                          normd(j,k,iaer,idomain)
                  omegrefVISgrid(j,k,iaer)=qscatrefVISgrid(j,k,iaer) /
     &                         qrefVISgrid(j,k,iaer)
                  DO m=1,nsun
                    qextVISgrid(j,k,m,iaer)=qextVISgrid(j,k,m,iaer) /
     &                          normd(j,k,iaer,idomain)
                    qscatVISgrid(j,k,m,iaer)=qscatVISgrid(j,k,m,iaer) /
     &                          normd(j,k,iaer,idomain)
                    gVISgrid(j,k,m,iaer)=gVISgrid(j,k,m,iaer) /
     &                          qscatVISgrid(j,k,m,iaer) /
     &                          normd(j,k,iaer,idomain)

                    qsqrefVISgrid(j,k,m,iaer)=qextVISgrid(j,k,m,iaer) /
     &                          qrefVISgrid(j,k,iaer)
                    omegVISgrid(j,k,m,iaer)=qscatVISgrid(j,k,m,iaer) /
     &                          qextVISgrid(j,k,m,iaer)
                  ENDDO
                ELSE                   ! INFRARED DOMAIN ----------
c                 2.3.3.ir Initialization
                  qsqrefIRgrid(j,k,:,iaer)=0.
                  qextIRgrid(j,k,:,iaer)=0.
                  qscatIRgrid(j,k,:,iaer)=0.
                  omegIRgrid(j,k,:,iaer)=0.
                  gIRgrid(j,k,:,iaer)=0.
                  qrefIRgrid(j,k,iaer)=0.
                  qscatrefIRgrid(j,k,iaer)=0.
                  omegrefIRgrid(j,k,iaer)=0.

                  DO gausind=1,ngau
                    DO m=1,nir
c                     Convolution:
                      qextIRgrid(j,k,m,iaer) =
     &                  qextIRgrid(j,k,m,iaer) +
     &                  weightgaus(gausind) *
     &                  (
     &                  qsqrefIRb(m,gausind,iaer) *
     &                  qrefVISb(gausind,iaer) *
     &                  pi*radGAUSb(gausind,iaer,idomain) *
     &                  radGAUSb(gausind,iaer,idomain) *
     &                  distb(j,k,iaer,idomain,gausind) +
     &                  qsqrefIRa(m,gausind,iaer) *
     &                  qrefVISa(gausind,iaer) *
     &                  pi*radGAUSa(gausind,iaer,idomain) *
     &                  radGAUSa(gausind,iaer,idomain) *
     &                  dista(j,k,iaer,idomain,gausind)
     &                  )
                      qscatIRgrid(j,k,m,iaer) =
     &                  qscatIRgrid(j,k,m,iaer) +
     &                  weightgaus(gausind) *
     &                  (
     &                  omegIRb(m,gausind,iaer) *
     &                  qsqrefIRb(m,gausind,iaer) *
     &                  qrefVISb(gausind,iaer) *
     &                  pi*radGAUSb(gausind,iaer,idomain) *
     &                  radGAUSb(gausind,iaer,idomain) *
     &                  distb(j,k,iaer,idomain,gausind) +
     &                  omegIRa(m,gausind,iaer) *
     &                  qsqrefIRa(m,gausind,iaer) *
     &                  qrefVISa(gausind,iaer) *
     &                  pi*radGAUSa(gausind,iaer,idomain) *
     &                  radGAUSa(gausind,iaer,idomain) *
     &                  dista(j,k,iaer,idomain,gausind)
     &                  )
                      gIRgrid(j,k,m,iaer) =
     &                  gIRgrid(j,k,m,iaer) +
     &                  weightgaus(gausind) *
     &                  (
     &                  omegIRb(m,gausind,iaer) *
     &                  qsqrefIRb(m,gausind,iaer) *
     &                  qrefVISb(gausind,iaer) *
     &                  gIRb(m,gausind,iaer) *
     &                  pi*radGAUSb(gausind,iaer,idomain) *
     &                  radGAUSb(gausind,iaer,idomain) *
     &                  distb(j,k,iaer,idomain,gausind) +
     &                  omegIRa(m,gausind,iaer) *
     &                  qsqrefIRa(m,gausind,iaer) *
     &                  qrefVISa(gausind,iaer) *
     &                  gIRa(m,gausind,iaer) *
     &                  pi*radGAUSa(gausind,iaer,idomain) *
     &                  radGAUSa(gausind,iaer,idomain) *
     &                  dista(j,k,iaer,idomain,gausind)
     &                  )
                    ENDDO
                    qrefIRgrid(j,k,iaer) =
     &                qrefIRgrid(j,k,iaer) +
     &                weightgaus(gausind) *
     &                (
     &                qrefIRb(gausind,iaer) *
     &                pi*radGAUSb(gausind,iaer,idomain) *
     &                radGAUSb(gausind,iaer,idomain) *
     &                distb(j,k,iaer,idomain,gausind) +
     &                qrefIRa(gausind,iaer) *
     &                pi*radGAUSa(gausind,iaer,idomain) *
     &                radGAUSa(gausind,iaer,idomain) *
     &                dista(j,k,iaer,idomain,gausind)
     &                )
                    qscatrefIRgrid(j,k,iaer) =
     &                qscatrefIRgrid(j,k,iaer) +
     &                weightgaus(gausind) *
     &                (
     &                omegrefIRb(gausind,iaer) *
     &                qrefIRb(gausind,iaer) *
     &                pi*radGAUSb(gausind,iaer,idomain) *
     &                radGAUSb(gausind,iaer,idomain) *
     &                distb(j,k,iaer,idomain,gausind) +
     &                omegrefIRa(gausind,iaer) *
     &                qrefIRa(gausind,iaer) *
     &                pi*radGAUSa(gausind,iaer,idomain) *
     &                radGAUSa(gausind,iaer,idomain) *
     &                dista(j,k,iaer,idomain,gausind)
     &                )
                  ENDDO
                  qrefIRgrid(j,k,iaer)=qrefIRgrid(j,k,iaer) /
     &                          normd(j,k,iaer,idomain)
                  qscatrefIRgrid(j,k,iaer)=qscatrefIRgrid(j,k,iaer) /
     &                          normd(j,k,iaer,idomain)
                  omegrefIRgrid(j,k,iaer)=qscatrefIRgrid(j,k,iaer) /
     &                         qrefIRgrid(j,k,iaer)
                  DO m=1,nir
                    qextIRgrid(j,k,m,iaer)=qextIRgrid(j,k,m,iaer) /
     &                          normd(j,k,iaer,idomain)
                    qscatIRgrid(j,k,m,iaer)=qscatIRgrid(j,k,m,iaer) /
     &                          normd(j,k,iaer,idomain)
                    gIRgrid(j,k,m,iaer)=gIRgrid(j,k,m,iaer) /
     &                          qscatIRgrid(j,k,m,iaer) /
     &                          normd(j,k,iaer,idomain)

                    qsqrefIRgrid(j,k,m,iaer)=qextIRgrid(j,k,m,iaer) /
     &                          qrefVISgrid(j,k,iaer)
                    omegIRgrid(j,k,m,iaer)=qscatIRgrid(j,k,m,iaer) /
     &                          qextIRgrid(j,k,m,iaer)
                  ENDDO
                ENDIF                  ! --------------------------
                checkgrid(j,k,iaer,idomain) = .true.
              ENDIF !checkgrid
            ENDDO !grid_j
          ENDDO !grid_i
c         3.4 Bilinear interpolation
          k1 = (1-kx)*(1-ky)
          k2 = kx*(1-ky)
          k3 = kx*ky
          k4 = (1-kx)*ky
          IF (idomain.EQ.1) THEN ! VISIBLE ------------------------
          DO m=1,nsun
            QVISsQREF3d(ig,lg,m,iaer) = 
     &                  k1*qsqrefVISgrid(grid_i,grid_j,m,iaer) +
     &                  k2*qsqrefVISgrid(grid_i+1,grid_j,m,iaer) +
     &                  k3*qsqrefVISgrid(grid_i+1,grid_j+1,m,iaer) +
     &                  k4*qsqrefVISgrid(grid_i,grid_j+1,m,iaer)
            omegaVIS3d(ig,lg,m,iaer) = 
     &                  k1*omegVISgrid(grid_i,grid_j,m,iaer) +
     &                  k2*omegVISgrid(grid_i+1,grid_j,m,iaer) +
     &                  k3*omegVISgrid(grid_i+1,grid_j+1,m,iaer) +
     &                  k4*omegVISgrid(grid_i,grid_j+1,m,iaer)
            gVIS3d(ig,lg,m,iaer) = 
     &                  k1*gVISgrid(grid_i,grid_j,m,iaer) +
     &                  k2*gVISgrid(grid_i+1,grid_j,m,iaer) +
     &                  k3*gVISgrid(grid_i+1,grid_j+1,m,iaer) +
     &                  k4*gVISgrid(grid_i,grid_j+1,m,iaer)
          ENDDO !nsun
          QREFvis3d(ig,lg,iaer) = 
     &                  k1*qrefVISgrid(grid_i,grid_j,iaer) +
     &                  k2*qrefVISgrid(grid_i+1,grid_j,iaer) +
     &                  k3*qrefVISgrid(grid_i+1,grid_j+1,iaer) +
     &                  k4*qrefVISgrid(grid_i,grid_j+1,iaer)
          omegaREFvis3d(ig,lg,iaer) = 
     &                  k1*omegrefVISgrid(grid_i,grid_j,iaer) +
     &                  k2*omegrefVISgrid(grid_i+1,grid_j,iaer) +
     &                  k3*omegrefVISgrid(grid_i+1,grid_j+1,iaer) +
     &                  k4*omegrefVISgrid(grid_i,grid_j+1,iaer)
          ELSE                   ! INFRARED -----------------------
          DO m=1,nir
            QIRsQREF3d(ig,lg,m,iaer) = 
     &                  k1*qsqrefIRgrid(grid_i,grid_j,m,iaer) + 
     &                  k2*qsqrefIRgrid(grid_i+1,grid_j,m,iaer) +
     &                  k3*qsqrefIRgrid(grid_i+1,grid_j+1,m,iaer) +
     &                  k4*qsqrefIRgrid(grid_i,grid_j+1,m,iaer)
            omegaIR3d(ig,lg,m,iaer) = 
     &                  k1*omegIRgrid(grid_i,grid_j,m,iaer) + 
     &                  k2*omegIRgrid(grid_i+1,grid_j,m,iaer) +
     &                  k3*omegIRgrid(grid_i+1,grid_j+1,m,iaer) +
     &                  k4*omegIRgrid(grid_i,grid_j+1,m,iaer)
            gIR3d(ig,lg,m,iaer) =
     &                  k1*gIRgrid(grid_i,grid_j,m,iaer) + 
     &                  k2*gIRgrid(grid_i+1,grid_j,m,iaer) +
     &                  k3*gIRgrid(grid_i+1,grid_j+1,m,iaer) +
     &                  k4*gIRgrid(grid_i,grid_j+1,m,iaer)
          ENDDO !nir
          QREFir3d(ig,lg,iaer) =
     &                  k1*qrefIRgrid(grid_i,grid_j,iaer) +
     &                  k2*qrefIRgrid(grid_i+1,grid_j,iaer) +
     &                  k3*qrefIRgrid(grid_i+1,grid_j+1,iaer) +
     &                  k4*qrefIRgrid(grid_i,grid_j+1,iaer)
          omegaREFir3d(ig,lg,iaer) =
     &                  k1*omegrefIRgrid(grid_i,grid_j,iaer) +
     &                  k2*omegrefIRgrid(grid_i+1,grid_j,iaer) +
     &                  k3*omegrefIRgrid(grid_i+1,grid_j+1,iaer) +
     &                  k4*omegrefIRgrid(grid_i,grid_j+1,iaer)
          ENDIF                  ! --------------------------------
        ENDDO !nlayer
      ENDDO !ngrid

      ENDIF ! varyingnueff
c==================================================================
      ENDDO ! idomain

      ENDIF ! nsize = 1

      ENDDO ! iaer (loop on aerosol kind)

c=====Radiative properties - TESTS=================================
      IF (out_qwg) THEN
c     -------------------------------------------------------------
        IF (ngrid.NE.1) THEN
          out_ndim = 3
        ELSE
          out_ndim = 1
        ENDIF
c     -------------------------------------------------------------
        DO out_nchannel = 1, 2
c     -------------------------------------------------------------
          DO lg = 1, nlayer
            DO ig = 1, ngrid
              out_qext(ig,lg) = 
     &           QVISsQREF3d(ig,lg,out_nchannel,out_iaer)*
     &           QREFvis3d(ig,lg,out_iaer)
              out_omeg(ig,lg) =
     &           omegaVIS3d(ig,lg,out_nchannel,out_iaer)
              out_g(ig,lg)    = gVIS3d(ig,lg,out_nchannel,out_iaer)
            ENDDO ! ig
          ENDDO ! lg
#ifndef MESOSCALE
          write(out_str(1:1),'(i1.1)') out_nchannel
          call WRITEDIAGFI(ngrid,'qextvis'//out_str,"Ext.efficiency","",
     &                     out_ndim,out_qext)
          call WRITEDIAGFI(ngrid,'omegvis'//out_str,"Sing.Scat.Alb.","",
     &                     out_ndim,out_omeg)
          call WRITEDIAGFI(ngrid,'gvis'//out_str,"Asym.Factor","",
     &                     out_ndim,out_g)
#endif
c     -------------------------------------------------------------
        ENDDO ! out_nchannel 
        DO out_nchannel = 2, 4
c     -------------------------------------------------------------
          DO lg = 1, nlayer
            DO ig = 1, ngrid
              out_qext(ig,lg) =
     &          QIRsQREF3d(ig,lg,out_nchannel,out_iaer)*
     &          QREFir3d(ig,lg,out_iaer)
              out_omeg(ig,lg) =
     &          omegaIR3d(ig,lg,out_nchannel,out_iaer)
              out_g(ig,lg)    = gIR3d(ig,lg,out_nchannel,out_iaer)
            ENDDO ! ig
          ENDDO ! lg
#ifndef MESOSCALE
          write(out_str(1:1),'(i1.1)') out_nchannel
       call WRITEDIAGFI(ngrid,'qextir'//out_str,"Ext.efficiency","",
     &                     out_ndim,out_qext)
       call WRITEDIAGFI(ngrid,'omegir'//out_str,"Sing.Scat.Alb.","",
     &                     out_ndim,out_omeg)
       call WRITEDIAGFI(ngrid,'gir'//out_str,"Asym.Factor","",
     &                     out_ndim,out_g)
#endif
c     -------------------------------------------------------------
        ENDDO ! out_nchannel 
#ifndef MESOSCALE
        call WRITEDIAGFI(ngrid,"omegvisref","Sing.Scat.Alb.","",
     &                   out_ndim,omegaREFvis3d(1,1,out_iaer))
        call WRITEDIAGFI(ngrid,"omegirref","Sing.Scat.Alb.","",
     &                   out_ndim,omegaREFir3d(1,1,out_iaer))
#endif
      ENDIF ! out_qwg
c==================================================================

      RETURN
      END
