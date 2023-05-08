!=======================================================================
! THERMCELL_DQUP
!=======================================================================
!
!   Compute the thermals contribution of explicit thermals 
!   to vertical transport in the PBL.
!   dq is computed once upward, entrainment and detrainment mass fluxes
!   are known.
!
!   Version with sub-timestep for Martian thin layers
!
!=======================================================================
! Author : A. Colaitis 2011-01-05 (with updates 2011-2013)
! Institution : Laboratoire de Meteorologie Dynamique (LMD) Paris, France
! -----------------------------------------------------------------------
! Corresponding author : A. Spiga aymeric.spiga_AT_upmc.fr
! -----------------------------------------------------------------------
! Reference paper:
! A. Colaïtis, A. Spiga, F. Hourdin, C. Rio, F. Forget, and E. Millour. 
! A thermal plume model for the Martian convective boundary layer. 
! Journal of Geophysical Research (Planets), 118:1468-1487, July 2013.
! http://dx.doi.org/10.1002/jgre.20104
! http://arxiv.org/abs/1306.6215
! -----------------------------------------------------------------------

      subroutine thermcell_dqup(ngrid,nlayer,ptimestep,fm,entr,detr,  &
     &    masse0,q_therm,dq_therm,ndt,limz)
      implicit none

! ============================ INPUTS ============================

      INTEGER, INTENT(IN) :: ngrid,nlayer ! number of grid points and number of levels
      REAL, INTENT(IN) :: ptimestep ! timestep (s)
      REAL, INTENT(IN) :: fm(ngrid,nlayer+1) ! upward mass flux
      REAL, INTENT(IN) :: entr(ngrid,nlayer) ! entrainment mass flux
      REAL, INTENT(IN) :: detr(ngrid,nlayer) ! detrainment mass flux
      REAL, INTENT(IN) :: q_therm(ngrid,nlayer) ! initial profil of q 
      REAL, INTENT(IN) :: masse0(ngrid,nlayer) ! mass of cells
      INTEGER, INTENT(IN) :: ndt ! number of subtimesteps
      INTEGER, INTENT(IN) :: limz ! index of maximum layer

! ============================ OUTPUTS ===========================

      REAL, INTENT(OUT) :: dq_therm(ngrid,nlayer)  ! dq/dt -> derivative

! ============================ LOCAL =============================

      REAL q(ngrid,nlayer)
      REAL qa(ngrid,nlayer)
      INTEGER ig,k,i
      REAL invflux0(ngrid,nlayer)
      REAL ztimestep

! =========== Init ==============================================

      qa(:,:)=q_therm(:,:) !q profile in the updraft
      q(:,:)=q_therm(:,:) !mean q profile

! ====== Computing q ============================================
! Based on equation 14 in appendix 4.2

      dq_therm(:,:)=0.
      ztimestep=ptimestep/real(ndt)
      invflux0(:,:)=ztimestep/masse0(:,:)      

      do i=1,ndt !subtimestep loop

        do ig=1,ngrid
           qa(ig,1)=q(ig,1)
       enddo

        do k=2,limz
           do ig=1,ngrid
              if ((fm(ig,k+1)+detr(ig,k))*ptimestep.gt.  &
     &        1.e-5*masse0(ig,k)) then
                 qa(ig,k)=(fm(ig,k)*qa(ig,k-1)+entr(ig,k)*q(ig,k))  &
     &           /(fm(ig,k+1)+detr(ig,k))
              else
                 qa(ig,k)=q(ig,k)
              endif
           enddo
        enddo

        do k=1,limz
          q(:,k)=q(:,k)+         &
     &    (detr(:,k)*qa(:,k)-entr(:,k)*q(:,k) &
     &    -fm(:,k)*q(:,k)+fm(:,k+1)*q(:,k+1))  &
     &    *invflux0(:,k)
        enddo

      enddo !of do i=1,ndt

! ====== Derivative ==============================================

         do k=1,limz
          dq_therm(:,k)=(q(:,k)-q_therm(:,k))/ptimestep
         enddo

! ==============

      return
      end

