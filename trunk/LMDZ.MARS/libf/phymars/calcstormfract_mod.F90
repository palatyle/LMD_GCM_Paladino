      MODULE calcstormfract_mod

      IMPLICIT NONE

      CONTAINS
!=======================================================================
! ROCKET DUST STORM - Defining dust storm fraction of the mesh
!=======================================================================
! before calling radiative transfer
! SUBGRID PARAMETRIZATION 
! -----------------------------------------------------------------------
! Author : T. Bertrand 2014-05 
! Institution : Laboratoire de Meteorologie Dynamique (LMD) Paris, France
! -----------------------------------------------------------------------
! defining the fraction of the mesh containing the rocket dust storm.
! file using the fraction: aeropacity.F, before calling radiative transfer
! the radiative transfer is performed in the mesh without the storm and in the mesh with the storm
! the fraction is used to calculate the true opacity of dust within the storm

      SUBROUTINE calcstormfract(ngrid,nlayer,nq,pq,totstormfract)
                                 
      use tracer_mod, only: igcm_stormdust_mass

      implicit none

!--------------------------------------------------------
! Input Variables
!--------------------------------------------------------

      INTEGER, INTENT(IN) :: ngrid ! number of horizontal grid points
      INTEGER, INTENT(IN) :: nlayer ! number of vertical grid points
      INTEGER, INTENT(IN) :: nq ! number of tracer species
      REAL, INTENT(IN) ::  pq(ngrid,nlayer,nq) ! advected field nq

!--------------------------------------------------------
! Output Variables
!--------------------------------------------------------

      REAL, INTENT(OUT) :: totstormfract(ngrid) ! fraction of the mesh containing the dust storm

!--------------------------------------------------------
! Local variables 
!--------------------------------------------------------

      REAL, PARAMETER :: mmr_ref=5.e-4 ! mass reference mixing ratio (corresponding
                                       ! to the tau=10 reference OMEGA local dust storm)
      REAL, PARAMETER :: fracmax=0.6    ! maximal fraction of the storm
      REAL, PARAMETER :: fracmin=1.e-2 ! minimal fraction of the storm
      INTEGER ig     
 
! **********************************************************************
! Definition of stormfraction
! **********************************************************************

      DO ig=1, ngrid
       totstormfract(ig)=abs(maxval(pq(ig,2:nlayer,  & !why from 2nd layer?
     &                     igcm_stormdust_mass)))/mmr_ref
       totstormfract(ig)=max(totstormfract(ig),fracmin)  
       totstormfract(ig)=min(totstormfract(ig),fracmax) 
      ENDDO

      END SUBROUTINE calcstormfract

      END MODULE calcstormfract_mod
