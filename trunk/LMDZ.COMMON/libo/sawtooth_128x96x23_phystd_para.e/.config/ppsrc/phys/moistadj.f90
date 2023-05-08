










subroutine moistadj(ngrid, nlayer, nq, pt, pq, pdq, pplev, pplay, pdtmana, pdqmana, ptimestep, rneb)

  use watercommon_h, only: T_h2O_ice_liq, RLVTT, RCPD, RCPV, Psat_water, Lcpdqsat_water
  USE tracer_h, only: igcm_h2o_vap, igcm_h2o_ice
  use comcstfi_mod, only: r

  implicit none


!=====================================================================
!     
!     Purpose
!     -------
!     Calculates moist convective adjustment by the method of Manabe.
!     
!     Authors
!     -------
!     Adapted from the LMDTERRE code by R. Wordsworth (2010)
!     Original author Z. X. Li (1993)
!     
!=====================================================================

      INTEGER,INTENT(IN) :: ngrid, nlayer, nq

      REAL,INTENT(IN) :: pt(ngrid,nlayer) ! temperature (K)
      REAL,INTENT(IN) :: pq(ngrid,nlayer,nq) ! tracer (kg/kg)
      REAL,INTENT(IN) :: pdq(ngrid,nlayer,nq)
      REAL,INTENT(IN) :: pplev(ngrid,nlayer+1) ! inter-layer pressure (Pa)
      REAL,INTENT(IN) :: pplay(ngrid,nlayer)   ! mid-layer pressure (Pa)
      REAL,INTENT(IN) :: ptimestep ! physics timestep (s)
      REAL,INTENT(OUT) :: pdqmana(ngrid,nlayer,nq)  ! tracer tendencies (kg/kg.s-1)
      REAL,INTENT(OUT) :: pdtmana(ngrid,nlayer) ! temperature increment(K/s)
      REAL,INTENT(OUT) :: rneb(ngrid,nlayer) ! cloud fraction 

!     local variables
      REAL zt(ngrid,nlayer)      ! temperature (K)
      REAL zq(ngrid,nlayer)      ! humidite specifique (kg/kg)

      REAL d_t(ngrid,nlayer)     ! temperature increment
      REAL d_q(ngrid,nlayer)     ! incrementation pour vapeur d'eau
      REAL d_ql(ngrid,nlayer)    ! incrementation pour l'eau liquide

!      REAL t_coup
!      PARAMETER (t_coup=234.0)
      REAL seuil_vap
      PARAMETER (seuil_vap=1.0E-10)

!     Local variables
      INTEGER i, k, iq, nn
      INTEGER, PARAMETER :: niter = 1
      INTEGER k1, k1p, k2, k2p
      LOGICAL itest(ngrid)
      REAL delta_q(ngrid, nlayer)
      DOUBLE PRECISION :: cp_new_t(nlayer), v_cptt(ngrid,nlayer)
      REAL cp_delta_t(nlayer)
      DOUBLE PRECISION :: v_cptj(nlayer), v_cptjk1, v_ssig
      REAL v_p, v_t, v_zqs,v_cptt2,v_pratio,v_dlnpsat
      REAL zqs(ngrid,nlayer), zdqs(ngrid,nlayer),zpsat(ngrid,nlayer),zdlnpsat(ngrid,nlayer)
      REAL zq1(ngrid), zq2(ngrid)
      DOUBLE PRECISION :: gamcpdz(ngrid,2:nlayer)
      DOUBLE PRECISION :: zdp, zdpm

      REAL zsat ! super-saturation
      REAL zflo ! flotabilite

      DOUBLE PRECISION :: local_q(ngrid,nlayer),local_t(ngrid,nlayer)

      REAL zdelta, zcor, zcvm5

      REAL dEtot, dqtot, masse ! conservation diagnostics
      real dL1tot, dL2tot

!     Indices of water vapour and water ice tracers
      INTEGER,SAVE :: i_h2o=0  ! water vapour
      INTEGER,SAVE :: i_ice=0  ! water ice
!$OMP THREADPRIVATE(i_h2o,i_ice)

      LOGICAL,SAVE :: firstcall=.TRUE.
!$OMP THREADPRIVATE(firstcall)

      IF (firstcall) THEN

         i_h2o=igcm_h2o_vap
         i_ice=igcm_h2o_ice
        
         write(*,*) "rain: i_ice=",i_ice
         write(*,*) "      i_h2o=",i_h2o

         firstcall = .FALSE.
      ENDIF

!     GCM -----> subroutine variables
      zq(1:ngrid,1:nlayer)    = pq(1:ngrid,1:nlayer,i_h2o)+ pdq(1:ngrid,1:nlayer,i_h2o)*ptimestep
      zt(1:ngrid,1:nlayer)    = pt(1:ngrid,1:nlayer)
      pdqmana(1:ngrid,1:nlayer,1:nq)=0.0

      DO k = 1, nlayer
       DO i = 1, ngrid
         if(zq(i,k).lt.0.)then
            zq(i,k)=0.0
         endif
       ENDDO
      ENDDO
      
      local_q(1:ngrid,1:nlayer) = zq(1:ngrid,1:nlayer)
      local_t(1:ngrid,1:nlayer) = zt(1:ngrid,1:nlayer)
      rneb(1:ngrid,1:nlayer) = 0.0
      d_ql(1:ngrid,1:nlayer) = 0.0
      d_t(1:ngrid,1:nlayer)  = 0.0
      d_q(1:ngrid,1:nlayer)  = 0.0

!     Calculate v_cptt
      DO k = 1, nlayer
         DO i = 1, ngrid
            v_cptt(i,k) = RCPD * local_t(i,k)
            v_t = MAX(local_t(i,k),15.)
            v_p = pplay(i,k)

            call Psat_water(v_t,v_p,zpsat(i,k),zqs(i,k))
	    call Lcpdqsat_water(v_t,v_p,zpsat(i,k),zqs(i,k),zdqs(i,k),zdlnpsat(i,k))
         ENDDO
      ENDDO

!     Calculate Gamma * Cp * dz: (gamma is the critical gradient)
      DO k = 2, nlayer
         DO i = 1, ngrid
            zdp = pplev(i,k)-pplev(i,k+1)
            zdpm = pplev(i,k-1)-pplev(i,k)
!            gamcpdz(i,k) = ( ( R/RCPD /(zdpm+zdp) * (v_cptt(i,k-1)*zdpm + v_cptt(i,k)*zdp)          &
!                +  RLVTT /(zdpm+zdp) * (zqs(i,k-1)*zdpm + zqs(i,k)*zdp) )                           &
!	        * (pplay(i,k-1)-pplay(i,k)) / pplev(i,k) )                                          &
!                / (1.0+ (zdqs(i,k-1)*zdpm + zdqs(i,k)*zdp)/(zdpm+zdp) )                    
! general case where water is not a trace gas (JL13)
            v_zqs     = (zqs(i,k-1)*zdpm + zqs(i,k)*zdp)/(zdpm+zdp)
	    v_cptt2   = (v_cptt(i,k-1)*zdpm + v_cptt(i,k)*zdp)/(zdpm+zdp)
	    v_pratio  = ((1.-zpsat(i,k-1)/pplay(i,k-1))*zdpm + (1.-zpsat(i,k)/pplay(i,k))*zdp)/(zdpm+zdp)
	    v_dlnpsat = (zdlnpsat(i,k-1)*zdpm + zdlnpsat(i,k)*zdp)/(zdpm+zdp)
            gamcpdz(i,k) = ( (R/RCPD*v_cptt2*(1.- v_zqs) + RLVTT*v_zqs) * (pplay(i,k-1)-pplay(i,k))/pplev(i,k) )  &
                / (((1.- v_zqs) + v_zqs * RCPV/RCPD)*v_pratio + v_zqs  * v_dlnpsat)                
         ENDDO
      ENDDO

!------------------------------------ modification of unstable profile
      DO 9999 i = 1, ngrid

      itest(i) = .FALSE.

!        print*,'we in the loop'
!        stop    

      k1 = 0
      k2 = 1

  810 CONTINUE ! look for k1, the base of the column
      k2 = k2 + 1
      IF (k2 .GT. nlayer) GOTO 9999
      zflo = v_cptt(i,k2-1) - v_cptt(i,k2) - gamcpdz(i,k2)
      zsat=(local_q(i,k2-1)-zqs(i,k2-1))*(pplev(i,k2-1)-pplev(i,k2))   &
         +(local_q(i,k2)-zqs(i,k2))*(pplev(i,k2)-pplev(i,k2+1))

      IF ( zflo.LE.0.0 .OR. zsat.LE.0.0 ) GOTO 810
      k1 = k2 - 1
      itest(i) = .TRUE.

  820 CONTINUE !! look for k2, the top of the column
      IF (k2 .EQ. nlayer) GOTO 821
      k2p = k2 + 1
      zsat=zsat+(pplev(i,k2p)-pplev(i,k2p+1))*(local_q(i,k2p)-zqs(i,k2p))
      zflo = v_cptt(i,k2p-1) - v_cptt(i,k2p) - gamcpdz(i,k2p)

      IF (zflo.LE.0.0 .OR. zsat.LE.0.0) GOTO 821
      k2 = k2p
      GOTO 820
  821 CONTINUE

!------------------------------------------------------ local adjustment
  830 CONTINUE ! actual adjustment
    Do nn=1,niter
      v_cptj(k1) = 0.0
      zdp = pplev(i,k1)-pplev(i,k1+1)
      v_cptjk1 = ( (1.0+zdqs(i,k1))*(v_cptt(i,k1)+v_cptj(k1)) + RLVTT*(local_q(i,k1)-zqs(i,k1)) ) * zdp
      v_ssig = zdp * (1.0+zdqs(i,k1))

      k1p = k1 + 1
      DO k = k1p, k2
         zdp = pplev(i,k)-pplev(i,k+1)
         v_cptj(k) = v_cptj(k-1) + gamcpdz(i,k)
         v_cptjk1 = v_cptjk1 + zdp * ( (1.0+zdqs(i, k))*(v_cptt(i,k)+v_cptj(k)) + RLVTT*(local_q(i,k)-zqs(i,k)) )        
         v_ssig = v_ssig + zdp *(1.0+zdqs(i,k))
      ENDDO


      ! this right here is where the adjustment is done???
      DO k = k1, k2
         cp_new_t(k) = v_cptjk1/v_ssig - v_cptj(k)
         cp_delta_t(k) = cp_new_t(k) - v_cptt(i,k)
         v_cptt(i,k)=cp_new_t(k)
         local_q(i,k) = zqs(i,k) + zdqs(i,k)*cp_delta_t(k)/RLVTT
         local_t(i,k) = cp_new_t(k) / RCPD

         v_t = MAX(local_t(i,k),15.)
         v_p = pplay(i,k)
         
         call Psat_water(v_t,v_p,zpsat(i,k),zqs(i,k))
         call Lcpdqsat_water(v_t,v_p,zpsat(i,k),zqs(i,k),zdqs(i,k),zdlnpsat(i,k))

      ENDDO
    Enddo ! nn=1,niter


!--------------------------------------------------- sounding downwards
!              -- we refine the prognostic variables in
!              -- the layer about to be adjusted

!      DO k = k1, k2
!         v_cptt(i,k) = RCPD * local_t(i,k)
!         v_t = local_t(i,k)
!         v_p = pplay(i,k)
!       
!         call Psat_water(v_t,v_p,zpsat,zqs(i,k))
!         call Lcpdqsat_water(v_t,v_p,zpsat,zqs(i,k),zdqs(i,k))
!      ENDDO

      DO k = 2, nlayer
         zdpm = pplev(i,k-1) - pplev(i,k)
         zdp = pplev(i,k) - pplev(i,k+1)
!         gamcpdz(i,k) = ( ( R/RCPD /(zdpm+zdp) * (v_cptt(i,k-1)*zdpm + v_cptt(i,k)*zdp)          &
!             +  RLVTT /(zdpm+zdp) * (zqs(i,k-1)*zdpm + zqs(i,k)*zdp) )                           &
!	      * (pplay(i,k-1)-pplay(i,k)) / pplev(i,k) )                                          &
!             / (1.0+ (zdqs(i,k-1)*zdpm + zdqs(i,k)*zdp)/(zdpm+zdp) )                    
! general case where water is not a trace gas (JL13)
         v_zqs     = (zqs(i,k-1)*zdpm + zqs(i,k)*zdp)/(zdpm+zdp)
	 v_cptt2   = (v_cptt(i,k-1)*zdpm + v_cptt(i,k)*zdp)/(zdpm+zdp)
	 v_pratio  = ((1.-zpsat(i,k-1)/pplay(i,k-1))*zdpm + (1.-zpsat(i,k)/pplay(i,k))*zdp)/(zdpm+zdp)
	 v_dlnpsat = (zdlnpsat(i,k-1)*zdpm + zdlnpsat(i,k)*zdp)/(zdpm+zdp)
         gamcpdz(i,k) = ( (R/RCPD*v_cptt2*(1.- v_zqs) + RLVTT*v_zqs) * (pplay(i,k-1)-pplay(i,k))/pplev(i,k) )  &
                / (((1.- v_zqs) + v_zqs * RCPV/RCPD)*v_pratio + v_zqs  * v_dlnpsat)                
      ENDDO

!     Test to see if we've reached the bottom

      IF (k1 .EQ. 1) GOTO 841 ! yes we have!
      zflo = v_cptt(i,k1-1) - v_cptt(i,k1) - gamcpdz(i,k1)
      zsat=(local_q(i,k1-1)-zqs(i,k1-1))*(pplev(i,k1-1)-pplev(i,k1))   &
        + (local_q(i,k1)-zqs(i,k1))*(pplev(i,k1)-pplev(i,k1+1))
      IF (zflo.LE.0.0 .OR. zsat.LE.0.0) GOTO 841 ! yes we have!

  840 CONTINUE
      k1 = k1 - 1
      IF (k1 .EQ. 1) GOTO 830 ! GOTO 820 (a tester, Z.X.Li, mars 1995)
      zsat = zsat + (local_q(i,k1-1)-zqs(i,k1-1))               &
                  *(pplev(i,k1-1)-pplev(i,k1))
      zflo = v_cptt(i,k1-1) - v_cptt(i,k1) - gamcpdz(i,k1)
      IF (zflo.GT.0.0 .AND. zsat.GT.0.0) THEN
         GOTO 840
      ELSE
         GOTO 830 ! GOTO 820 (a tester, Z.X.Li, mars 1995)
      ENDIF
  841 CONTINUE

      GOTO 810 ! look for other layers higher up

 9999 CONTINUE ! loop over all the points

!-----------------------------------------------------------------------
! Determine the cloud fraction (hypothese: la nebulosite a lieu
! a l'endroit ou la vapeur d'eau est diminuee par l'ajustement):

      DO k = 1, nlayer
      DO i = 1, ngrid
         IF (itest(i)) THEN
         delta_q(i,k) = local_q(i,k) - zq(i,k)
         IF (delta_q(i,k).LT.0.) rneb(i,k)  = 1.0
         ENDIF
      ENDDO
      ENDDO

! Distribuer l'eau condensee en eau liquide nuageuse (hypothese:
! l'eau liquide est distribuee aux endroits ou la vapeur d'eau
! diminue et d'une maniere proportionnelle a cet diminution):

      DO i = 1, ngrid
         IF (itest(i)) THEN
         zq1(i) = 0.0
         zq2(i) = 0.0
         ENDIF
      ENDDO
      DO k = 1, nlayer
      DO i = 1, ngrid
         IF (itest(i)) THEN
         zdp = pplev(i,k)-pplev(i,k+1)
         zq1(i) = zq1(i) - delta_q(i,k) * zdp
         zq2(i) = zq2(i) - MIN(0.0, delta_q(i,k)) * zdp
         ENDIF
      ENDDO
      ENDDO
      DO k = 1, nlayer
      DO i = 1, ngrid
         IF (itest(i)) THEN
           IF (zq2(i).NE.0.0) d_ql(i,k) = - MIN(0.0,delta_q(i,k))*zq1(i)/zq2(i)
         ENDIF
      ENDDO
      ENDDO

      DO k = 1, nlayer
      DO i = 1, ngrid
          local_q(i, k) = MAX(local_q(i, k), seuil_vap)
      ENDDO
      ENDDO

      DO k = 1, nlayer
      DO i = 1, ngrid
         d_t(i,k) = local_t(i,k) - zt(i,k)
         d_q(i,k) = local_q(i,k) - zq(i,k)
      ENDDO
      ENDDO

!     now subroutine -----> GCM variables
      DO k = 1, nlayer
         DO i = 1, ngrid
            
            pdtmana(i,k)       = d_t(i,k)/ptimestep
            pdqmana(i,k,i_h2o) = d_q(i,k)/ptimestep
            pdqmana(i,k,i_ice) = d_ql(i,k)/ptimestep
         
         ENDDO
      ENDDO


   END
