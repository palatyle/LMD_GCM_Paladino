










!
! $Id: calfis_p.F 1407 2010-07-07 10:31:52Z fairhead $
!
C
C
      SUBROUTINE calfis_p(lafin,
     $                  jD_cur, jH_cur,
     $                  pucov,
     $                  pvcov,
     $                  pteta,
     $                  pq,
     $                  pmasse,
     $                  pps,
     $                  pp,
     $                  ppk,
     $                  pphis,
     $                  pphi,
     $                  pducov,
     $                  pdvcov,
     $                  pdteta,
     $                  pdq,
     $                  flxw,
     $                  pdufi,
     $                  pdvfi,
     $                  pdhfi,
     $                  pdqfi,
     $                  pdpsfi)
! Ehouarn: if using (parallelized) physics
      USE dimphy
      USE mod_phys_lmdz_mpi_data, mpi_root_xx=>mpi_master
      USE mod_phys_lmdz_omp_data, ONLY: klon_omp, klon_omp_begin
      USE mod_const_mpi, ONLY: COMM_LMDZ
      USE mod_interface_dyn_phys
!      USE IOPHY
      USE infotrac, ONLY: nqtot, niadv, tname
      USE control_mod, ONLY: planet_type, nsplit_phys
      USE comvert_mod, ONLY: preff,presnivs
      USE comconst_mod, ONLY: daysec,dtvr,dtphys,kappa,cpp,g,rad,pi
      USE logic_mod, ONLY: moyzon_ch,moyzon_mu
      USE callphysiq_mod, ONLY: call_physiq

      IMPLICIT NONE
c=======================================================================
c
c   1. rearrangement des tableaux et transformation
c      variables dynamiques  >  variables physiques
c   2. calcul des termes physiques
c   3. retransformation des tendances physiques en tendances dynamiques
c
c   remarques:
c   ----------
c
c    - les vents sont donnes dans la physique par leurs composantes 
c      naturelles.
c    - la variable thermodynamique de la physique est une variable
c      intensive :   T 
c      pour la dynamique on prend    T * ( preff / p(l) ) **kappa
c    - les deux seules variables dependant de la geometrie necessaires
c      pour la physique sont la latitude pour le rayonnement et 
c      l'aire de la maille quand on veut integrer une grandeur 
c      horizontalement.
c    - les points de la physique sont les points scalaires de la 
c      la dynamique; numerotation:
c          1 pour le pole nord
c          (jjm-1)*iim pour l'interieur du domaine
c          ngridmx pour le pole sud
c      ---> ngridmx=2+(jjm-1)*iim
c
c     Input :
c     -------
c       ecritphy        frequence d'ecriture (en jours)de histphy
c       pucov           covariant zonal velocity
c       pvcov           covariant meridional velocity 
c       pteta           potential temperature
c       pps             surface pressure
c       pmasse          masse d'air dans chaque maille
c       pts             surface temperature  (K)
c       callrad         clef d'appel au rayonnement
c
c    Output :
c    --------
c        pdufi          tendency for the natural zonal velocity (ms-1)
c        pdvfi          tendency for the natural meridional velocity 
c        pdhfi          tendency for the potential temperature (K/s)
c        pdtsfi         tendency for the surface temperature
c
c        pdtrad         radiative tendencies  \  both input
c        pfluxrad       radiative fluxes      /  and output
c
c=======================================================================
c
c-----------------------------------------------------------------------
c
c    0.  Declarations :
c    ------------------

      include "dimensions.h"
      include "paramet.h"

      INTEGER ngridmx
      PARAMETER( ngridmx = 2+(jjm-1)*iim - 1/jjm   )

      include "comgeom2.h"
      include "iniprint.h"
c    Arguments :
c    -----------
      LOGICAL,INTENT(IN) ::  lafin ! .true. for the very last call to physics
      REAL,INTENT(IN) :: jD_cur, jH_cur
      REAL,INTENT(IN) :: pvcov(iip1,jjm,llm) ! covariant meridional velocity
      REAL,INTENT(IN) :: pucov(iip1,jjp1,llm) ! covariant zonal velocity
      REAL,INTENT(IN) :: pteta(iip1,jjp1,llm) ! potential temperature
      REAL,INTENT(IN) :: pmasse(iip1,jjp1,llm) ! mass in each cell ! not used
      REAL,INTENT(IN) :: pq(iip1,jjp1,llm,nqtot) ! tracers
      REAL,INTENT(IN) :: pphis(iip1,jjp1) ! surface geopotential
      REAL,INTENT(IN) :: pphi(iip1,jjp1,llm) ! geopotential

      REAL,INTENT(IN) :: pdvcov(iip1,jjm,llm) ! dynamical tendency on vcov
      REAL,INTENT(IN) :: pducov(iip1,jjp1,llm) ! dynamical tendency on ucov
      REAL,INTENT(IN) :: pdteta(iip1,jjp1,llm) ! dynamical tendency on teta
! commentaire SL: pdq ne sert que pour le calcul de pcvgq,
! qui lui meme ne sert a rien dans la routine telle qu'elle est
! ecrite, et que j'ai donc commente....
      REAL,INTENT(IN) :: pdq(iip1,jjp1,llm,nqtot) ! dynamical tendency on tracers
      ! NB: pdq is only used to compute pcvgq which is in fact not used...

      REAL,INTENT(IN) :: pps(iip1,jjp1) ! surface pressure (Pa)
      REAL,INTENT(IN) :: pp(iip1,jjp1,llmp1) ! pressure at mesh interfaces (Pa)
      REAL,INTENT(IN) :: ppk(iip1,jjp1,llm) ! Exner at mid-layer
      REAL,INTENT(IN) :: flxw(iip1,jjp1,llm)  ! Vertical mass flux on lower mesh interfaces (kg/s) (on llm because flxw(:,:,llm+1)=0)

      ! tendencies (in */s) from the physics
      REAL,INTENT(OUT) :: pdvfi(iip1,jjm,llm) ! tendency on covariant meridional wind
      REAL,INTENT(OUT) :: pdufi(iip1,jjp1,llm) ! tendency on covariant zonal wind
      REAL,INTENT(OUT) :: pdhfi(iip1,jjp1,llm) ! tendency on potential temperature (K/s)
      REAL,INTENT(OUT) :: pdqfi(iip1,jjp1,llm,nqtot) ! tendency on tracers
      REAL,INTENT(OUT) :: pdpsfi(iip1,jjp1) ! tendency on surface pressure (Pa/s)

! of #ifdef CPP_PARA
      END
