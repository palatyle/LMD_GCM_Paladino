!
!
!
SUBROUTINE thermcell_dv2(ngrid,nlay,ptimestep,fm,entr,detr,masse,fraca,       &
                         zmax,zmin,u,v,du,dv,ua,va)
     
      
!===============================================================================
!
!   Calcul du transport vertical dans la couche limite en presence
!   de "thermiques" explicitement representes
!   Calcul du dq/dt une fois qu'on connait les ascendances
!
! Vectorisation, FH : 2010/03/08
!
!===============================================================================
      
      IMPLICIT none
      
      
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
      REAL, INTENT(in) :: fraca(ngrid,nlay+1)
      REAL, INTENT(in) :: zmax(ngrid)
      REAL, INTENT(in) :: zmin(ngrid)
      REAL, INTENT(in) :: u(ngrid,nlay)
      REAL, INTENT(in) :: v(ngrid,nlay)
      
!     Outputs:
!     --------
      
      REAL, INTENT(out) :: ua(ngrid,nlay)       ! u in the plume
      REAL, INTENT(out) :: va(ngrid,nlay)       ! v in the plume
      REAL, INTENT(out) :: du(ngrid,nlay)       ! large scale u variation
      REAL, INTENT(out) :: dv(ngrid,nlay)       ! large scale v variation
      
!     Local:
!     ------
      
      INTEGER ig, l
      INTEGER iter
      
      REAL qa(ngrid,nlay)
      REAL zf
      REAL zf2
      REAL wvd(ngrid,nlay+1)
      REAL wud(ngrid,nlay+1)
      REAL gamma0(ngrid,nlay+1)
      REAL gamma(ngrid,nlay+1)
      REAL ue(ngrid,nlay)
      REAL ve(ngrid,nlay)
      REAL dua(ngrid,nlay)
      REAL dva(ngrid,nlay)
      REAL plume_height(ngrid)
      
      LOGICAL ltherm(ngrid,nlay)
      
!===============================================================================
! Initialization
!===============================================================================
      
      DO ig=1,ngrid
         plume_height(ig) = zmax(ig) - zmin(ig)
      ENDDO
      
      DO ig=1,ngrid
         ua(ig,1)=u(ig,1)
         va(ig,1)=v(ig,1)
         ue(ig,1)=u(ig,1)
         ve(ig,1)=v(ig,1)
      ENDDO
      
      gamma0(1:ngrid,1)=0.
      
      DO l=2,nlay
         DO ig=1,ngrid
            ltherm(ig,l) = (fm(ig,l+1) + detr(ig,l)) * ptimestep > 1.e-5 * masse(ig,l)
!            IF (ltherm(ig,l).and.(zmax(ig) > 0.)) THEN
            IF (ltherm(ig,l).and.(plume_height(ig) > 0.)) THEN
!               gamma0(ig,l) = masse(ig,l) * 0.5 / zmax(ig)                    &
               gamma0(ig,l) = masse(ig,l) * 0.5 / plume_height(ig)            &
               &            * sqrt(0.5 * (fraca(ig,l+1) + fraca(ig,l)))
            ELSE
               gamma0(ig,l) = 0.
            ENDIF
         ENDDO
      ENDDO
      
      gamma(:,:) = 0.
      
!===============================================================================
! 
!===============================================================================
      
      DO l=2,nlay
         
         DO ig=1,ngrid
            IF (ltherm(ig,l)) THEN
               dua(ig,l) = ua(ig,l-1) - u(ig,l-1)
               dva(ig,l) = va(ig,l-1) - v(ig,l-1)
            ELSE
               ua(ig,l) = u(ig,l)
               va(ig,l) = v(ig,l)
               ue(ig,l) = u(ig,l)
               ve(ig,l) = v(ig,l)
            ENDIF
         ENDDO
         
         DO iter=1,5
            DO ig=1,ngrid
! FH: Calcul prenant en compte la fraction reelle
               zf = 0.5 * (fraca(ig,l) + fraca(ig,l+1))
               zf2 = 1./(1.-zf)
! FH: Calcul avec fraction infiniement petite
!               zf = 0.
!               zf2 = 1.
               
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! FH: La première fois on multiplie le coefficient de freinage par le module du
!     vent dans la couche en dessous. Mais pourquoi donc ?
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               IF (ltherm(ig,l)) THEN
! FH: On choisit une relaxation lineaire.
!                  gamma(ig,l) = gamma0(ig,l)
! FH: On choisit une relaxation quadratique.
                  gamma(ig,l) = gamma0(ig,l) * sqrt(dua(ig,l)**2 + dva(ig,l)**2)
                  
                  ua(ig,l) = (fm(ig,l) * ua(ig,l-1)                           &
                  &        + (zf2 * entr(ig,l) + gamma(ig,l)) * u(ig,l))      &
                  &        / (fm(ig,l+1) + detr(ig,l) + entr(ig,l) * zf * zf2 &
                  &        + gamma(ig,l))
                  va(ig,l) = (fm(ig,l) * va(ig,l-1)                           &
                  &        + (zf2 * entr(ig,l) + gamma(ig,l)) * v(ig,l))      &
                  &        / (fm(ig,l+1) + detr(ig,l) + entr(ig,l) * zf * zf2 &
                  &        + gamma(ig,l))
                  dua(ig,l) = ua(ig,l) - u(ig,l)
                  dva(ig,l) = va(ig,l) - v(ig,l)
                  ue(ig,l) = (u(ig,l) - zf * ua(ig,l)) * zf2
                  ve(ig,l) = (v(ig,l) - zf * va(ig,l)) * zf2
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      
!-------------------------------------------------------------------------------
! Calcul du flux vertical de moment dans l'environnement
!-------------------------------------------------------------------------------
      
      DO l=2,nlay
         DO ig=1,ngrid
            wud(ig,l) = fm(ig,l) * ue(ig,l)
            wvd(ig,l) = fm(ig,l) * ve(ig,l)
         ENDDO
      ENDDO
      
      DO ig=1,ngrid
         wud(ig,1) = 0.
         wvd(ig,1) = 0.
         wud(ig,nlay+1) = 0.
         wvd(ig,nlay+1) = 0.
      ENDDO
      
!===============================================================================
! Tendencies
!===============================================================================
      
      DO l=1,nlay
         DO ig=1,ngrid
            du(ig,l) = ((detr(ig,l) + gamma(ig,l)) * ua(ig,l)                 &
            &        - (entr(ig,l) + gamma(ig,l)) * ue(ig,l)                  &
            &        - wud(ig,l) + wud(ig,l+1))                               &
            &        / masse(ig,l)
            dv(ig,l) = ((detr(ig,l) + gamma(ig,l)) * va(ig,l)                 &
            &        - (entr(ig,l) + gamma(ig,l)) * ve(ig,l)                  &
            &        - wvd(ig,l) + wvd(ig,l+1))                               &
            &        / masse(ig,l)
         ENDDO
      ENDDO
      
      
RETURN
END
