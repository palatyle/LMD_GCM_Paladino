      subroutine new_cloud_sedim(n_lon, n_lev, ptimestep,
     $                           pmidlay, pbndlay, pt, pq,
     $                           d_tr_chem, pdqsed,
     $                           nq, F_sed)

      USE ioipsl
      USE dimphy
      USE chemparam_mod
      IMPLICIT NONE

c-----------------------------------------------------------------------
c   declarations:
c   -------------
#include "YOMCST.h"     
c#include "dimphys.h"
c#include "comcstfi.h"
c#include "tracer.h"
c#include "callkeys.h"
c
c   arguments:
c   ----------

      INTEGER n_lon                 ! number of horizontal grid points
      INTEGER n_lev                 ! number of atmospheric layers
      REAL ptimestep                ! physics time step (s)
      REAL pmidlay(n_lon,n_lev)     ! pressure at middle layers (Pa)
      REAL pt(n_lon,n_lev)          ! temperature at mid-layer (l)
      REAL pbndlay(n_lon,n_lev+1)   ! pressure at layer boundaries

c    Traceurs :
      integer nq                    ! number of tracers
      real pq(n_lon,n_lev,nq)       ! tracers (kg/kg)
      real pdqsed(n_lon,n_lev,2)    ! tendency due to sedimentation (kg/kg)
      real d_tr_chem(n_lon,n_lev,nq)! tendency due to chemistry and clouds (kg/kg)
      
c   local:
c   ------
      integer imode
      integer ig
      integer iq
      integer l

      real zlev(n_lon,n_lev+1)      ! altitude at layer boundaries
      real zlay(n_lon,n_lev)        ! altitude at the midlle layer
      real zqi_wv(n_lon,n_lev)      ! to locally store H2O tracer
      real zqi_sa(n_lon,n_lev)      ! to locally store H2SO4 tracer
      real m_lay (n_lon,n_lev)      ! Layer Pressure over gravity (Dp/g == kg.m-2)
      real wq(n_lon,n_lev+1)        ! displaced tracer mass (kg.m-2)

c    Physical constant
c    ~~~~~~~~~~~~~~~~~
c     Gas molecular viscosity (N.s.m-2)
c      real,parameter :: visc=1.e-5       ! CO2
      REAL :: VISCOSITY_CO2
c     Effective gas molecular radius (m)
      real,parameter :: molrad=2.2e-10   ! CO2
      
c     Ratio radius shell model du mode 3
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	Ce ratio correspond aux mesures effectuées par J. Cimino (1982), Icarus
c     Fixer ce parametre a 0 revient a une gouttelette pure en liquide acide sulfurique
c     ATTENTION ! DOIT ETRE COHERENT AVEC new_cloud_venus !
      REAL, PARAMETER :: qrad = 0.97
      REAL :: qmass
c	masse volumique du coeur (kg.m-3)
c     ATTENTION ! DOIT ETRE COHERENT AVEC new_cloud_venus !
      REAL, PARAMETER :: rho_core = 2500.0

      REAL, DIMENSION(n_lon,n_lev+1) :: wgt_SA   ! Fraction of H2SO4 in droplet local

c     Stokes speed and sedimentation flux variable
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      REAL :: A1,A2,A3,A4,            ! coeff du DL du Flux de sedimentation
     + D_stokes,                      ! coeff de la vitesse de Stokes
     + Rp_DL,                         ! "Point" du DL
     + l_mean,                        ! libre parcours moyen (m)
     + a,b_exp,c                      ! coeff du calcul du Flux de sedimentation
      REAL, DIMENSION(n_lon,n_lev+1) ::
     + F_sed                          ! Flux de sedimentation (kg.m-2.s-1 puis en output kg.m-2)
      
      REAL :: R_mode0                 ! Rayon mode 0 (m), rayon le plus frequent

!      PRINT*,'RHO_DROPLET new_cloud_sedim.F'
!      PRINT*,'rho_droplet',rho_droplet(16,21)
!      PRINT*,'T',pt(16,21),'WSA',WH2SO4(16,21)

c-----------------------------------------------------------------------
c    1. Initialization
c    -----------------
    
!     update water vapour and sulfuric acid mixing ratios

      zqi_wv(:,:) = pq(:,:,i_h2oliq) + d_tr_chem(:,:,i_h2oliq)*ptimestep
      zqi_sa(:,:) = pq(:,:,i_h2so4liq) 
     $            + d_tr_chem(:,:,i_h2so4liq)*ptimestep

      wgt_SA(:,:) = wh2so4(:,:)

c     Init F_sed
      F_sed(:,:) = 0.0E+0

c     Au niveau top+1 , tout égal a 0      
      wgt_SA(:,n_lev+1) = 0.0E+0   

c    Computing the different layer properties
c    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c    m_lay (kg.m-2)
c    Ici g=8.87, conflit pour g entre #include "YOMCST.h" 
c	et #include "comcstfi.h"

      do  l=1,n_lev
         do ig=1, n_lon
         m_lay(ig,l)=(pbndlay(ig,l) - pbndlay(ig,l+1)) /8.87E+0 
            IF (m_lay(ig,l).LE.0.0) THEN
            PRINT*,'!!!! STOP PROBLEME SEDIMENTATION!!!!'
            PRINT*,'!!!!          m_lay <= 0        !!!!'
            PRINT*,'!!!! STOP PROBLEME SEDIMENTATION!!!!'
            ENDIF
         end do
      end do
	
c         Computing sedimentation for droplet "tracer"
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	    pbndlay(:,51)=0 (en parallèle c'est sûr), ne pas l'utiliser pour Fse
	
c     Sedimentation pour une gouttelette mode 3 de type J. Cimino, 1982, Icarus
c     c.a.d 97% radius due a solide 3% radius acide sulfurique	
	DO imode=1, nbr_mode - 1
         DO l = cloudmin, cloudmax
            DO ig=1,n_lon

c     RD=1000.*RNAVO*RKBOL/RMD avec RMD=43.44 Masse molaire atm venus en g.mol-1     
           D_stokes=((rho_droplet(ig,l)-pmidlay(ig,l)/(RD*pt(ig,l))))
     &	*(2./9.)*(RG/VISCOSITY_CO2(pt(ig,l)))
      
           l_mean=(pt(ig,l)/pmidlay(ig,l))*
     &	(0.707*R/(4.*RPI* molrad*molrad * RNAVO))
      
           R_mode0=R_MEDIAN(ig,l,imode)*
     &     EXP(-LOG(STDDEV(ig,l,imode))**2.)
              IF ((l_mean/(R_mode0)).GT.10.) THEN
              Rp_DL=R_MEDIAN(ig,l,imode)*
     &        EXP(3.*LOG(STDDEV(ig,l,imode))**2.)
              ELSE
              Rp_DL=R_MEDIAN(ig,l,imode)*
     &        EXP(4.*LOG(STDDEV(ig,l,imode))**2.)
              ENDIF
		
           a=1.246*l_mean
	
           c=0.87/l_mean
	
           b_exp=0.42*l_mean*EXP(-c*Rp_DL)
	
           A1=a+b_exp*(1.+c*Rp_DL
     &	+0.5*(Rp_DL*c)**2
     &	+1./6.*(Rp_DL*c)**3)
     
           A2=1.-b_exp*(c
     &	+Rp_DL*c**2
     &	+0.5*(Rp_DL**2)*(c**3))
	
           A3=0.5*b_exp*(c**2+Rp_DL*c**3)
	
           A4=-b_exp*1./6.*c**3

c     Addition des Flux de tous les modes presents	
       F_sed(ig,l)=F_sed(ig,l)+(rho_droplet(ig,l)*4./3.*RPI*
     &  NBRTOT(ig,l,imode)*1.0E6*D_stokes*(
     &  A1*R_MEDIAN(ig,l,imode)**4
     &  *EXP(8.0*LOG(STDDEV(ig,l,imode))**2.)
     &  +A2*R_MEDIAN(ig,l,imode)**5
     &  *EXP(12.5*LOG(STDDEV(ig,l,imode))**2.)
     &  +A3*R_MEDIAN(ig,l,imode)**6
     &  *EXP(18.0*LOG(STDDEV(ig,l,imode))**2.)
     &  +A4*R_MEDIAN(ig,l,imode)**7
     &  *EXP(24.5*LOG(STDDEV(ig,l,imode))**2.)))
      
c      PRINT*,' APRES dTime: F_sed=',F_sed(ig,l), ig, l
      
      	IF (F_sed(ig,l).GT.m_lay(ig,l)) THEN
      	PRINT*,'==============================================='
      	PRINT*,'WARNING On a epuise la couche', ig, l
      	PRINT*,'On epuise pas une couche avec une espèce
     & 	 minoritaire, c est pas bien maaaaaal'
            PRINT*,'Water',zqi_wv(ig,l),'Sulfuric Acid',zqi_sa(ig,l)
      	PRINT*,'F_sed:',F_sed(ig,l),'m_lay:',m_lay(ig,l)
      	PRINT*,'F_sed/dtphy',F_sed(ig,l)/ptimestep
      	PRINT*,'Pbnd top',pbndlay(ig,l+1),'Temp',pt(ig,l),'Rho',
     & 	rho_droplet(ig,l)
     		PRINT*,'Ntot',NBRTOT(ig,l,:)
     		PRINT*,'StdDev',STDDEV(ig,l,:),'Rmed',R_MEDIAN(ig,l,:)
     		PRINT*,'K_MASS',K_MASS(ig,l,:)
     		PRINT*,'WSA',WH2SO4(ig,l),'RHO',rho_droplet(ig,l)
     		
c     		ELSE
c     		
c     		PRINT*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c      	PRINT*,'WARNING On a PAS epuise la couche', ig, l
c      	PRINT*,'F_sed:',F_sed(ig,l),'m_lay:',m_lay(ig,l)
c      	PRINT*,'F_sed/dtphy',F_sed(ig,l)/ptimestep
c      	PRINT*,'Pbnd top',pbndlay(ig,l+1),'Temp',pt(ig,l),'Rho',
c     & 	rho_droplet(ig,l)(ig,l)
c     		PRINT*,'Ntot',NBRTOT(ig,l),'Ntot m3',NBRTOT(ig,l)*1.0e6
c     		PRINT*,'StdDev',STDDEV(ig,l),'Rmed',R_MEDIAN(ig,l) 
            STOP    		
     	      ENDIF
	
           IF (F_sed(ig,l).LT.0.0d0) THEN
	      PRINT*,"F_sed est négatif !!!"
	      PRINT*,'F_sed:',F_sed(ig,l),'m_lay:',m_lay(ig,l)
      	PRINT*,'F_sed/dtphy',F_sed(ig,l)/ptimestep
      	PRINT*,'Pbnd top',pbndlay(ig,l+1),'Pmid',pmidlay(ig,l)
      	PRINT*,'Temp',pt(ig,l),'Rho',
     & 	rho_droplet(ig,l)
     		PRINT*,'Ntot',NBRTOT(ig,l,imode),'Ntot m3',
     &	NBRTOT(ig,l,imode)*1.0e6
     		PRINT*,'StdDev',STDDEV(ig,l,imode),'Rmed',
     &		R_MEDIAN(ig,l,imode)
     		PRINT*,'A1',A1,'A2',A2
     		PRINT*,'A3',A1,'A4',A2
     		PRINT*,'D_stokes',D_stokes
     		STOP
           ENDIF
           
	      ENDDO
	      
c           ELSE           
c           F_sed(:,l)=0.0d0           
c           ENDIF
           
         ENDDO
      ENDDO

c****************************************************************
c        On calcule le F_sed du mode 3 + coeff*(Fsed1 + Fsed2)
c****************************************************************
         DO l = cloudmin, cloudmax
            DO ig=1,n_lon

c	calcul de qmass
	    qmass=(rho_core*qrad**3)/
     &    (rho_core*qrad**3+rho_droplet(ig,l)*(1.-qrad**3))

c     RD=1000.*RNAVO*RKBOL/RMD avec RMD=43.44 Masse molaire atm venus en g.mol-1    
           D_stokes=(((qmass*rho_core+(1.-qmass)*rho_droplet(ig,l))
     &      -pmidlay(ig,l)/(RD*pt(ig,l))))
     &	*(2./9.)*(RG/VISCOSITY_CO2(pt(ig,l)))
      
           l_mean=(pt(ig,l)/pmidlay(ig,l))*
     &	(0.707*R/(4.*RPI* molrad*molrad * RNAVO))
      
           R_mode0=R_MEDIAN(ig,l,3)*
     &     EXP(-LOG(STDDEV(ig,l,3))**2.)
              IF ((l_mean/(R_mode0)).GT.10.) THEN
              Rp_DL=R_MEDIAN(ig,l,3)*
     &        EXP(3.*LOG(STDDEV(ig,l,3))**2.)
              ELSE
              Rp_DL=R_MEDIAN(ig,l,3)*
     &        EXP(4.*LOG(STDDEV(ig,l,3))**2.)
              ENDIF
		
           a=1.246*l_mean
	
           c=0.87/l_mean
	
           b_exp=0.42*l_mean*EXP(-c*Rp_DL)
	
           A1=a+b_exp*(1.+c*Rp_DL
     &	+0.5*(Rp_DL*c)**2
     &	+1./6.*(Rp_DL*c)**3)
     
           A2=1.-b_exp*(c
     &	+Rp_DL*c**2
     &	+0.5*(Rp_DL**2)*(c**3))
	
           A3=0.5*b_exp*(c**2+Rp_DL*c**3)
	
           A4=-b_exp*1./6.*c**3

c     Addition des Flux de tous les modes presents	
       F_sed(ig,l)=F_sed(ig,l)
     &  +((1.-qmass)/(1.-qmass*K_MASS(ig,l,3)))*(
     &  (qmass*rho_core+(1.-qmass)*rho_droplet(ig,l))*4./3.*RPI*
     &  NBRTOT(ig,l,3)*1.0E6*D_stokes*(
     &  A1*R_MEDIAN(ig,l,3)**4
     &  *EXP(8.0*LOG(STDDEV(ig,l,3))**2.)
     &  +A2*R_MEDIAN(ig,l,3)**5
     &  *EXP(12.5*LOG(STDDEV(ig,l,3))**2.)
     &  +A3*R_MEDIAN(ig,l,3)**6
     &  *EXP(18.0*LOG(STDDEV(ig,l,3))**2.)
     &  +A4*R_MEDIAN(ig,l,3)**7
     &  *EXP(24.5*LOG(STDDEV(ig,l,3))**2.)))
      
c      PRINT*,' APRES dTime: F_sed=',F_sed(ig,l), ig, l
      
      	IF (F_sed(ig,l).GT.m_lay(ig,l)) THEN
      	PRINT*,'==============================================='
      	PRINT*,'WARNING On a epuise la couche', ig, l
      	PRINT*,'On epuise pas une couche avec une espèce
     & 	 minoritaire, c est pas bien maaaaaal'
            PRINT*,'Water',zqi_wv(ig,l),'Sulfuric Acid',zqi_sa(ig,l)
      	PRINT*,'F_sed:',F_sed(ig,l),'m_lay:',m_lay(ig,l)
      	PRINT*,'F_sed/dtphy',F_sed(ig,l)/ptimestep
      	PRINT*,'Pbnd top',pbndlay(ig,l+1),'Temp',pt(ig,l),'Rho',
     & 	rho_droplet(ig,l)
     		PRINT*,'Ntot',NBRTOT(ig,l,:)
     		PRINT*,'StdDev',STDDEV(ig,l,:),'Rmed',R_MEDIAN(ig,l,:)
     		PRINT*,'K_MASS',K_MASS(ig,l,:)
     		PRINT*,'WSA',WH2SO4(ig,l),'RHO',rho_droplet(ig,l)
     		
c     		ELSE
c     		
c     		PRINT*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
c      	PRINT*,'WARNING On a PAS epuise la couche', ig, l
c      	PRINT*,'F_sed:',F_sed(ig,l),'m_lay:',m_lay(ig,l)
c      	PRINT*,'F_sed/dtphy',F_sed(ig,l)/ptimestep
c      	PRINT*,'Pbnd top',pbndlay(ig,l+1),'Temp',pt(ig,l),'Rho',
c     & 	rho_droplet(ig,l)(ig,l)
c     		PRINT*,'Ntot',NBRTOT(ig,l),'Ntot m3',NBRTOT(ig,l)*1.0e6
c     		PRINT*,'StdDev',STDDEV(ig,l),'Rmed',R_MEDIAN(ig,l) 
            STOP    		
     	      ENDIF
	
           IF (F_sed(ig,l).LT.0.0d0) THEN
	      PRINT*,"F_sed est négatif !!!"
	      PRINT*,'F_sed:',F_sed(ig,l),'m_lay:',m_lay(ig,l)
      	PRINT*,'F_sed/dtphy',F_sed(ig,l)/ptimestep
      	PRINT*,'Pbnd top',pbndlay(ig,l+1),'Pmid',pmidlay(ig,l)
      	PRINT*,'Temp',pt(ig,l),'Rho',
     & 	rho_droplet(ig,l)
     		PRINT*,'Ntot',NBRTOT(ig,l,imode),'Ntot m3',
     &	NBRTOT(ig,l,imode)*1.0e6
     		PRINT*,'StdDev',STDDEV(ig,l,imode),'Rmed',
     &		R_MEDIAN(ig,l,imode)
     		PRINT*,'A1',A1,'A2',A2
     		PRINT*,'A3',A1,'A4',A2
     		PRINT*,'D_stokes',D_stokes
     		STOP
           ENDIF
           
	      ENDDO
	      
c           ELSE           
c           F_sed(:,l)=0.0d0           
c           ENDIF
           
         ENDDO

c     Passage du Flux au Flux pour un pas de temps (== kg.m-2)     

      F_sed(:,:) = F_sed(:,:)*ptimestep

!=========================================================
!     compute tendency due to sedimentation
!=========================================================

!     h2so4

      do l = 1,n_lev
         do ig = 1,n_lon
            zqi_sa(ig,l) = zqi_sa(ig,l) 
     $                   + (F_sed(ig,l+1)*wgt_SA(ig,l+1)
     $                    - F_sed(ig,l)*wgt_SA(ig,l))/m_lay(ig,l)
!           if (zqi_sa(ig,l) < 0.) THEN
!              print*,'STOP sedim on epuise tout le H2SO4l present'
!              print*,'point ',ig,'level ',l
!              print*,'zqi_sa = ', zqi_sa(ig,l)
!              STOP
!              zqi_sa(ig,l) = 0.
!           end if
            zqi_sa(ig,l) = max(zqi_sa(ig,l), 0.)
            pdqsed(ig,l,1) = zqi_sa(ig,l) - pq(ig,l,i_h2so4liq)
         end do
      end do

!     h2o
                     
      do l = 1, n_lev
         do ig=1,n_lon
            zqi_wv(ig,l) = zqi_wv(ig,l) 
     $                   + (F_sed(ig,l+1)*(1. - wgt_SA(ig,l+1))
     &                    - F_sed(ig,l)*(1. - wgt_SA(ig,l)))
     &                    /m_lay(ig,l)
!           if (zqi_wv(ig,l) < 0.) THEN
!              print*,'STOP sedim on epuise tout le H2Ol present'
!              print*,'point ',ig,'level ',l
!              print*,'zqi_wv = ', zqi_wv(ig,l)
!              STOP
!              zqi_wv(ig,l) = 0.
!           end if
            zqi_wv(ig,l) = max(zqi_wv(ig,l), 0.)
            pdqsed(ig,l,2) = zqi_wv(ig,l) - pq(ig,l,i_h2oliq) 
         end do
      end do

c		Save output file in 1D model
c		============================ 
c      IF (n_lon .EQ. 1) THEN
c      PRINT*,'Save output sedim' 	
c      DO l = 1, n_lev
c      	DO ig=1,n_lon
c      	WRITE(77,"(i4,','11(e15.8,','))") l,pdqsed(ig,l),zqi(ig,l),
c     & 	(WH2SO4(ig,l)*pq(ig,l,i_h2so4liq)+
c     & 	(1.-WH2SO4(ig,l))*pq(ig,l,i_h2oliq)),
c     & 	pq(ig,l,i_h2so4liq),pq(ig,l,i_h2oliq)
c      ENDDO
c      	ENDDO
c      ENDIF   

      RETURN
      END

*******************************************************************************
      REAL FUNCTION VISCOSITY_CO2(temp)
c	Aurélien Stolzenbach 2015
c	Calcul de la viscosité dynamique du CO2 80°K -> 300°K
c	Viscosité dynamique en Pa.s
c	Source: Johnston & Grilly (1942)

c	température en °K
	REAL, INTENT(IN) :: temp
	
	REAL :: denom, numer
	
c	Calcul de la viscosité dynamique grâce à la formule de Jones (Lennard-Jones (1924))
	
	numer = 200.**(2.27/4.27)-0.435
	denom = temp**(2.27/4.27)-0.435
	
	VISCOSITY_CO2 = (numer/denom)*1015.*(temp/200.)**(3./2.)

c	convertion de Poises*1e7 -> Pa.s	
	VISCOSITY_CO2 = VISCOSITY_CO2*1.e-8	

      END FUNCTION VISCOSITY_CO2
*******************************************************************************


