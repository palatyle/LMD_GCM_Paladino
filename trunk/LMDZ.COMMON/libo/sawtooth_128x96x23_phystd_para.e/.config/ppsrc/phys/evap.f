










      subroutine evap(ngrid,nlayer,nq,dtime,pt, pq, pdq, pdt,  
     $     dqevap,dtevap, qevap, tevap) 
        
      use watercommon_h 
      USE tracer_h

      implicit none

!==================================================================
!     
!     Purpose
!     -------
!     Evaporation of all water in the atmopshere.
!     
!     Authors
!     -------
!     Adapted from the LMDTERRE code by B. Charnay (2010)
!     Original author Z. X. Li (1993)
!     
!==================================================================

      INTEGER ngrid,nlayer,nq

! Arguments:
      REAL pt(ngrid,nlayer) 
      REAL pq(ngrid,nlayer,nq) 
      REAL pdt(ngrid,nlayer) 
      REAL pdq(ngrid,nlayer,nq) 
      REAL dqevap(ngrid,nlayer)
      REAL dtevap(ngrid,nlayer)
      REAL qevap(ngrid,nlayer,nq)
      REAL dtime
 
! Local:
      REAL tevap(ngrid,nlayer) 
      REAL zlvdcp
      REAL zlsdcp
      REAL zdelta
      INTEGER l,ig

!
! Re-evaporer l'eau liquide nuageuse
!

      DO l=1,nlayer
        DO ig=1,ngrid
         qevap(ig,l,igcm_h2o_vap)=pq(ig,l,igcm_h2o_vap)	
     s          	          +pdq(ig,l,igcm_h2o_vap)*dtime
         qevap(ig,l,igcm_h2o_ice)=pq(ig,l,igcm_h2o_ice)	
     s                            +pdq(ig,l,igcm_h2o_ice)*dtime
         tevap(ig,l)=pt(ig,l)+pdt(ig,l)*dtime
	ENDDO
      ENDDO

      DO l = 1, nlayer
      	DO ig = 1, ngrid
         zlvdcp=RLVTT/RCPD!/(1.0+RVTMP2*qevap(ig,l,igcm_h2o_vap))
         zlsdcp=RLSTT/RCPD!/(1.0+RVTMP2*qevap(ig,l,igcm_h2o_vap))
         ! ignoring qevap term creates huge difference when qevap large!!!

         zdelta = MAX(0.,SIGN(1.,T_h2O_ice_liq-tevap(ig,l))) ! what is this?
                                                  ! for division between water / liquid
	 dqevap(ig,l) = MAX(0.0,qevap(ig,l,igcm_h2o_ice))/dtime
         dtevap(ig,l) = - dqevap(ig,l)*RLVTT/RCPD ! exactly as in largescale.F
!         dtevap(ig,l) = - dqevap(ig,l)	
!     s                       * (zlvdcp*(1.-zdelta)+zlsdcp*zdelta)
         qevap(ig,l,igcm_h2o_vap) = qevap(ig,l,igcm_h2o_vap)	
     s     	                    +dqevap(ig,l)*dtime
	 qevap(ig,l,igcm_h2o_ice) = 0.0
         tevap(ig,l) = tevap(ig,l)+dtevap(ig,l)*dtime

      	ENDDO
      ENDDO

      END
