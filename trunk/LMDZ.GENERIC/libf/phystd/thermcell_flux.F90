!
!
!
SUBROUTINE thermcell_flux(ngrid,nlay,ptimestep,masse,                         &
                          lmin,lmax,entr_star,detr_star,                      &
                          f,rhobarz,zlev,zw2,fm,entr,detr)
      
      
!===============================================================================
!  Purpose: deduction des flux
!  
!  Modif 2019/04 (AB alexandre.boissinot@lmd.jussieu.fr)
!  
!===============================================================================
      
      USE print_control_mod, ONLY: prt_level
      USE thermcell_mod, ONLY: fomass_max, alpha_max
      
      IMPLICIT NONE
      
      
!===============================================================================
! Declaration
!===============================================================================
      
!     Inputs:
!     ------- 
      
      INTEGER, INTENT(in) :: ngrid
      INTEGER, INTENT(in) :: nlay
      INTEGER, INTENT(in) :: lmin(ngrid)
      
      REAL, INTENT(in) :: entr_star(ngrid,nlay)
      REAL, INTENT(in) :: detr_star(ngrid,nlay)
      REAL, INTENT(in) :: zw2(ngrid,nlay+1)
      REAL, INTENT(in) :: zlev(ngrid,nlay+1)
      REAL, INTENT(in) :: masse(ngrid,nlay)
      REAL, INTENT(in) :: ptimestep
      REAL, INTENT(in) :: rhobarz(ngrid,nlay)
      REAL, INTENT(in) :: f(ngrid)
      
!     Outputs:
!     --------
      
      INTEGER, INTENT(inout) :: lmax(ngrid)
      
      REAL, INTENT(out) :: entr(ngrid,nlay)
      REAL, INTENT(out) :: detr(ngrid,nlay)
      REAL, INTENT(out) :: fm(ngrid,nlay+1)
      
!     Local:
!     ------
      
      INTEGER ig, l, k
      INTEGER igout, lout                 ! Error grid point and level
      
      REAL fmax                           ! Maximal authorized mass flux (alpha < alpha_max)
      REAL fff0                           ! Save initial value of mass flux
      REAL emax                           ! Maximal authorized entrainment (entr*dt < mass_max)
      REAL eee0                           ! Save initial value of entrainment
      REAL ddd0                           ! Save initial value of detrainment
      REAL eee                            ! eee0 - layer mass * maximal authorized mass fraction / dt
      REAL ddd                            ! ddd0 - eee
      REAL fact
      REAL test
      
      INTEGER ncorecentr
      INTEGER ncorecdetr
      INTEGER nerrorequa
      INTEGER ncorecfact
      INTEGER ncorecalpha
      
      LOGICAL labort_physic
      
!===============================================================================
! Initialization
!===============================================================================
      
      nerrorequa = 0
      ncorecentr = 0
      ncorecdetr = 0
      ncorecfact = 0
      ncorecalpha = 0
      
      entr(:,:) = 0.
      detr(:,:) = 0.
      fm(:,:)   = 0.
      
      labort_physic = .false.
      
      fact = 0.
      
!===============================================================================
! Calcul de l'entrainement, du detrainement et du flux de masse
!===============================================================================
      
!-------------------------------------------------------------------------------
! Multiplication par la norme issue de la relation de fermeture
!-------------------------------------------------------------------------------
      
      DO l=1,nlay
         entr(:,l) = f(:) * entr_star(:,l)
         detr(:,l) = f(:) * detr_star(:,l)
      ENDDO
      
!-------------------------------------------------------------------------------
! Mass flux
!-------------------------------------------------------------------------------
      
      DO l=1,nlay
         DO ig=1,ngrid
            IF (l < lmax(ig) .and. l >= lmin(ig)) THEN
               fm(ig,l+1) = fm(ig,l) + entr(ig,l) - detr(ig,l)
            ELSEIF (l == lmax(ig)) THEN
               fm(ig,l+1) = 0.
               entr(ig,l) = 0.
               detr(ig,l) = fm(ig,l) + entr(ig,l)
            ELSE
               fm(ig,l+1) = 0.
               entr(ig,l) = 0.
               detr(ig,l) = 0.
            ENDIF
         ENDDO
      ENDDO
      
!===============================================================================
! Checking
!===============================================================================
      
      DO l=1,nlay
         
!-------------------------------------------------------------------------------
! Is incoming mass flux positive ?
!-------------------------------------------------------------------------------
         
         DO ig=1,ngrid
            IF (fm(ig,l) < 0.) THEN
               labort_physic = .true.
               igout = ig
               lout = l
            ENDIF
         ENDDO
         
!-------------------------------------------------------------------------------
! Is entrainment positive ?
!-------------------------------------------------------------------------------
         
         DO ig=1,ngrid
            IF (entr(ig,l) < 0.) THEN
               labort_physic = .true.
               igout = ig
               lout = l
            ENDIF
         ENDDO
         
!-------------------------------------------------------------------------------
! Is detrainment positive ?
!-------------------------------------------------------------------------------
         
         DO ig=1,ngrid
            IF (detr(ig,l) < 0.) THEN
               labort_physic = .true.
               igout = ig
               lout = l
            ENDIF
         ENDDO
         
!-------------------------------------------------------------------------------
! Abort
!-------------------------------------------------------------------------------
         
         IF (labort_physic) THEN
            print *, '---------------------------------------------------------'
            print *, 'ERROR: mass flux has negative value(s)!'
            print *, 'ig,l,norm', igout, lout, f(igout)
            print *, 'lmin,lmax', lmin(igout), lmax(igout)
            print *, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
            DO k=nlay,1,-1
               print *, 'k, fm ', k+1, fm(igout,k+1)
               print *, 'entr,detr', entr(igout,k), detr(igout,k)
            ENDDO
            print *, 'k, fm ', 1, fm(igout,1)
            print *, '---------------------------------------------------------'
            CALL abort
         ENDIF
         
!-------------------------------------------------------------------------------
! Is entrained mass lesser than fomass_max ?
!-------------------------------------------------------------------------------
         
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! AB : Entrainment is bigger than the maximal authorized value.
!      If we consider that the excess entrainement is in fact plume air which
!      is not detrained then we compensate it by decreasing detr.
!      If it's not enough, we can increase entr in the layer above and decrease
!      the outgoing mass flux in the current layer.
!      If it's still insufficient, code will abort (now commented).
!      Else we reset entr to its intial value and we renormalize entrainment,
!      detrainment and mass flux profiles such as we do not exceed the maximal
!      authorized entrained mass.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         
         DO ig=1,ngrid
            eee0 = entr(ig,l)
            ddd0 = detr(ig,l)
            emax = masse(ig,l) * fomass_max / ptimestep
            IF (emax < 0.) THEN
               print *, 'ERROR: layer mass is negative!'
               print *, 'mass,emax', masse(ig,l), emax
               print *, 'ig,l', ig, l
            ENDIF
            IF (eee0 > emax) THEN
               entr(ig,l) = emax
               ddd = ddd0 + emax - eee0
               ncorecentr  = ncorecentr + 1
               IF (ddd > 0.) THEN
                  detr(ig,l) = ddd
               ELSEIF (l == lmax(ig)) THEN
                  detr(ig,l) = fm(ig,l) + entr(ig,l)
               ELSEIF (entr(ig,l+1) > ABS(ddd)) THEN
                  detr(ig,l) = 0.
                  fm(ig,l+1) = fm(ig,l) + entr(ig,l)
                  entr(ig,l+1) = entr(ig,l+1) + ddd
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! AB: Simulation abort (try to reduce the physical time step).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               ELSE
                  entr(ig,l) = entr(ig,l) + eee
                  igout = ig
                  lout = l
                  labort_physic = .true.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! AB: We can renormalize the plume mass fluxes. I think it does not work.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!               ELSE
!                  fact = max(fact, eee0 / emax)
!                  entr(ig,l) = eee0
!                  ncorecfact = ncorecfact + 1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! AB: The renormalization can be just applied in the local plume.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                  DO k=lmin(ig),lmax(ig)
!                     entr(ig,k) = entr(ig,k) * emax / eee0
!                     detr(ig,k) = detr(ig,k) * emax / eee0
!                     fm(ig,k) = fm(ig,k) * emax / eee0
!                  ENDDO
               ENDIF
            ENDIF
         ENDDO
         
         IF (labort_physic) THEN
            print *, '---------------------------------------------------------'
            print *, 'ERROR: Entrainment is greater than maximal authorized value!'
            print *, '       Nor detrainment neither entrainment can compensate it!'
            print *, 'ig,l,entr', igout, lout, entr(igout,lout)
            print *, 'lmin,lmax', lmin(igout), lmax(igout)
            print *, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
            print *, 'e_max :', masse(igout,lout)*fomass_max/ptimestep
            print *, '   fomass_max :', fomass_max
            print *, '   masse      :', masse(igout,lout)
            print *, '   ptimestep  :', ptimestep
            print *, 'norm  :', f(igout)
            print *, 'entr* :', entr_star(igout,lout)
            print *, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
            DO k=nlay,1,-1
               print *, 'fm ', fm(igout,k+1)
               print *, 'entr,detr', entr(igout,k), detr(igout,k)
            ENDDO
            print *, 'fm ', fm(igout,1)
            print *, '---------------------------------------------------------'
            CALL abort
         ENDIF
         
!-------------------------------------------------------------------------------
! Is updraft fraction lesser than alpha_max ?
!-------------------------------------------------------------------------------
         
         DO ig=1,ngrid
            fff0 = fm(ig,l+1)
            fmax = rhobarz(ig,l+1) * zw2(ig,l+1) * alpha_max
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! AB: The plume mass flux can be reduced.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!            IF (fff0 > fmax) THEN
!               fm(ig,l+1) = fmax
!               detr(ig,l) = detr(ig,l) + fff0 - fmax
!               ncorecalpha = ncorecalpha + 1
!               entr(ig,l+1) = entr(ig,l+1) + fff0 - fmax
!            ENDIF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! AB: The plume can be stopped here.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            IF (fff0 > fmax) THEN
               ncorecalpha = ncorecalpha + 1
               DO k=l+1,lmax(ig)
                  entr(ig,k) = 0.
                  detr(ig,k) = 0.
                  fm(ig,k) = 0.
               ENDDO
               lmax(ig) = l
               entr(ig,l) = 0.
               detr(ig,l) = fm(ig,l)
            ENDIF
         ENDDO
         
!-------------------------------------------------------------------------------
! Is detrainment lesser than incoming mass flux ?
!-------------------------------------------------------------------------------
         
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! AB : Even if fm has no negative value, it can be lesser than detr.
!      That is not suitable because when we will mix the plume with the
!      environment, it will detrain more mass than it is physically able to do.
!      When it occures, that imply that entr + fm is greater than detr,
!      otherwise we get a negative outgoing mass flux (cf. above).
!      That is why we decrease entrainment and detrainment as follows.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         
         DO ig=1,ngrid
            IF (detr(ig,l) > fm(ig,l)) THEN
               detr(ig,l) = fm(ig,l)
               entr(ig,l) = fm(ig,l+1)
               ncorecdetr  = ncorecdetr + 1
            ENDIF
         ENDDO
         
      ENDDO
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! AB: The renormalization can be applied everywhere.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      IF (fact > 0.) THEN
!         entr(:,:) = entr(:,:) / fact
!         detr(:,:) = detr(:,:) / fact
!         fm(:,:) = fm(:,:) / fact
!      ENDIF
      
!-------------------------------------------------------------------------------
! Is equation df/dz = e - d still verified ?
!-------------------------------------------------------------------------------
      
!      DO l=1,nlay
!         DO ig=1,ngrid
!            test = abs(fm(ig,l) + entr(ig,l) - detr(ig,l) - fm(ig,l+1))
!            IF (test > 1.e-10) THEN
!               nerrorequa = nerrorequa + 1
!            ENDIF
!         ENDDO
!      ENDDO
      
!-------------------------------------------------------------------------------
! Reset top boundary conditions
!-------------------------------------------------------------------------------
      
      DO ig=1,ngrid
         IF (lmax(ig) > 0) THEN
            detr(ig,lmax(ig)) = fm(ig,lmax(ig))
            fm(ig,lmax(ig)+1) = 0.
            entr(ig,lmax(ig)) = 0.
         ENDIF
      ENDDO
      
!===============================================================================
! Outputs
!===============================================================================
      
      IF (prt_level > 0) THEN
         
         IF (ncorecdetr > 0) THEN
            print *, 'WARNING: Detrainment is greater than mass flux!'
            print *, 'In', ncorecdetr, 'grid point(s) over', nlay, 'x', ngrid
         ENDIF
         
         IF (ncorecentr > 0) THEN
            print *, 'WARNING: Entrained mass is greater than maximal authorized value!'
            print *, 'In', ncorecentr, 'grid point(s) over', nlay, 'x', ngrid
         ENDIF
         
         IF (ncorecfact > 0) THEN
            print *, 'WARNING: Entrained mass needs renormalization to be ok!'
            print *, 'In', ncorecfact, 'grid point(s) over', nlay, 'x', ngrid
!            print *, 'WARNING: Entr fact:', fact
         ENDIF
         
!         IF (nerrorequa > 0) THEN
!            print *, 'WARNING: !'
!            print *, 'in', nerrorequa, 'grid point(s) over', nlay, 'x', ngrid
!         ENDIF
         
         IF (ncorecalpha > 0) THEN
            print *, 'WARNING: Updraft fraction is greater than maximal authorized value!'
            print *, 'In', ncorecalpha, 'grid point(s) over', nlay, 'x', ngrid
         ENDIF
         
      ENDIF
      
      
RETURN
END
