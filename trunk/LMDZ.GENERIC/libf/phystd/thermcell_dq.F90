!
!
!
SUBROUTINE thermcell_dq(ngrid,nlay,ptimestep,fm,entr,detr,masse,              &
                        q,dq,qa)
      
      
!===============================================================================
!  Purpose: Calcul du transport verticale dans la couche limite en presence de
!           "thermiques" explicitement representes
!           Calcul du dq/dt une fois qu'on connait les ascendances
!  
!  Modif 2013/01/04 (FH hourdin@lmd.jussieu.fr)
!  Introduction of an implicit computation of vertical advection in the environ-
!     ment of thermal plumes in thermcell_dq
!  
!  Modif 2019/04 (AB alexandre.boissinot@lmd.jussieu.fr)
!     dqimpl = true  : implicit scheme
!     dqimpl = false : explicit scheme
!  
!===============================================================================
      
      USE print_control_mod, ONLY: prt_level
      USE thermcell_mod, ONLY: dqimpl
      
      IMPLICIT NONE
      
      
!===============================================================================
! Declaration
!===============================================================================
      
!     Inputs:
!     -------
      
      INTEGER, INTENT(in) :: ngrid
      INTEGER, INTENT(in) :: nlay
      
      REAL, INTENT(in) :: ptimestep
      REAL, INTENT(in) :: masse(ngrid,nlay)
      REAL, INTENT(in) :: fm(ngrid,nlay+1)
      REAL, INTENT(in) :: entr(ngrid,nlay)
      REAL, INTENT(in) :: detr(ngrid,nlay)
      
!     Outputs:
!     --------
      
      REAL, INTENT(inout) :: q(ngrid,nlay)
      REAL, INTENT(out) :: dq(ngrid,nlay)
      REAL, INTENT(out) :: qa(ngrid,nlay)
      
!     Local:
!     ------
      
      INTEGER ig, l, k
      INTEGER niter, iter
      
      REAL cfl
      REAL qold(ngrid,nlay)
      REAL fqa(ngrid,nlay+1)
      REAL zzm
      
!===============================================================================
! Initialization
!===============================================================================
      
      qold(:,:) = q(:,:)
      
!===============================================================================
! Tracer variation computation
!===============================================================================
      
!-------------------------------------------------------------------------------
! CFL criterion computation for advection in downdraft
!-------------------------------------------------------------------------------
      
      cfl = 0.
      
      DO l=1,nlay
         DO ig=1,ngrid
            zzm = masse(ig,l) / ptimestep
            cfl = max(cfl, fm(ig,l) / zzm)
            
            IF (entr(ig,l) > zzm) THEN
               print *, 'ERROR: entrainment is greater than the layer mass!'
               print *, 'ig,l,entr', ig, l, entr(ig,l)
               print *, '-------------------------------'
               print *, 'entr*dt,mass', entr(ig,l)*ptimestep, masse(ig,l)
               print *, '-------------------------------'
               DO k=nlay,1,-1
                  print *, 'fm ', fm(ig,k+1)
                  print *, 'entr,detr', entr(ig,k), detr(ig,k)
               ENDDO
               print *, 'fm ', fm(ig,1)
               print *, '-------------------------------'
               CALL abort
            ENDIF
         ENDDO
      ENDDO
      
!-------------------------------------------------------------------------------
! Computation of tracer concentrations in the ascending plume
!-------------------------------------------------------------------------------
      
      DO ig=1,ngrid
         DO l=1,nlay
            IF ((fm(ig,l+1)+detr(ig,l))*ptimestep > 1.e-6*masse(ig,l)) THEN
               qa(ig,l) = (fm(ig,l) * qa(ig,l-1) + entr(ig,l) * q(ig,l))      &
               &        / (fm(ig,l+1) + detr(ig,l))
            ELSE
               qa(ig,l) = q(ig,l)
            ENDIF
         ENDDO
      ENDDO
      
!-------------------------------------------------------------------------------
! Plume vertical flux of tracer
!-------------------------------------------------------------------------------
      
      DO l=2,nlay-1
         fqa(:,l) = fm(:,l) * qa(:,l-1)
      ENDDO
      
      fqa(:,1) = 0.
      fqa(:,nlay) = 0.
      
!-------------------------------------------------------------------------------
! Trace species evolution
!-------------------------------------------------------------------------------
      
      IF (dqimpl) THEN
         DO l=nlay-1,1,-1
            q(:,l) = ( q(:,l) + ptimestep / masse(:,l)                        &
            &      * ( fqa(:,l) - fqa(:,l+1) + fm(:,l+1) * q(:,l+1) ) )       &
            &      / ( 1. + fm(:,l) * ptimestep / masse(:,l) )
         ENDDO
      ELSE
         DO l=1,nlay-1
            q(:,l) = q(:,l) + (fqa(:,l) - fqa(:,l+1) - fm(:,l) * q(:,l)       &
            &      + fm(:,l+1) * q(:,l+1)) * ptimestep / masse(:,l)
         ENDDO
      ENDIF
      
!===============================================================================
! Tendencies
!===============================================================================
      
      DO l=1,nlay
         DO ig=1,ngrid
            dq(ig,l) = (q(ig,l) - qold(ig,l)) / ptimestep
            q(ig,l) = qold(ig,l)
         ENDDO
      ENDDO
      
      
RETURN
END
