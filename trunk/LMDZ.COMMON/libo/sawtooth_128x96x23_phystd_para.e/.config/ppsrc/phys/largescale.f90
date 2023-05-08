










      subroutine largescale(ngrid,nlayer,nq,ptimestep, pplev, pplay,    &
                    pt, pq, pdt, pdq, pdtlsc, pdqvaplsc, pdqliqlsc, rneb)


      use ioipsl_getin_p_mod, only: getin_p
      use watercommon_h, only : RLVTT, RCPD, RVTMP2,  &
          T_h2O_ice_clouds,T_h2O_ice_liq,Psat_water,Lcpdqsat_water
      USE tracer_h
      IMPLICIT none

!==================================================================
!     
!     Purpose
!     -------
!     Calculates large-scale (stratiform) H2O condensation.
!     
!     Authors
!     -------
!     Adapted from the LMDTERRE code by R. Wordsworth (2009)
!     Original author Z. X. Li (1993)
!     
!==================================================================

      INTEGER ngrid,nlayer,nq

!     Arguments
      REAL ptimestep                 ! intervalle du temps (s)
      REAL pplev(ngrid,nlayer+1) ! pression a inter-couche
      REAL pplay(ngrid,nlayer)   ! pression au milieu de couche
      REAL pt(ngrid,nlayer)      ! temperature (K)
      REAL pq(ngrid,nlayer,nq) ! tracer mixing ratio (kg/kg)
      REAL pdt(ngrid,nlayer)     ! physical temperature tenedency (K/s)
      REAL pdq(ngrid,nlayer,nq)! physical tracer tenedency (K/s)
      REAL pdtlsc(ngrid,nlayer)  ! incrementation de la temperature (K)
      REAL pdqvaplsc(ngrid,nlayer) ! incrementation de la vapeur d'eau
      REAL pdqliqlsc(ngrid,nlayer) ! incrementation de l'eau liquide
      REAL rneb(ngrid,nlayer)    ! fraction nuageuse


!     Options du programme
      REAL, SAVE :: ratqs   ! determine largeur de la distribution de vapeur
!$OMP THREADPRIVATE(ratqs)

!     Variables locales
      REAL CBRT
      EXTERNAL CBRT
      INTEGER i, k , nn
      INTEGER,PARAMETER :: nitermax=5000
      DOUBLE PRECISION,PARAMETER :: alpha=.1,qthreshold=1.d-8
      ! JL13: if "careful, T<Tmin in psat water" appears often, you may want to stabilise the model by
      !                   decreasing alpha and increasing nitermax accordingly
      DOUBLE PRECISION zq(ngrid)
      DOUBLE PRECISION zcond(ngrid),zcond_iter
      DOUBLE PRECISION zdelq(ngrid)
      DOUBLE PRECISION zqs(ngrid)
      real zt(ngrid),local_p,psat_tmp,dlnpsat_tmp,Lcp,zqs_temp,zdqs
      
! evaporation calculations
      REAL dqevap(ngrid,nlayer),dtevap(ngrid,nlayer)     
      REAL qevap(ngrid,nlayer,nq)
      REAL tevap(ngrid,nlayer)

      DOUBLE PRECISION zx_q(ngrid)
      LOGICAL,SAVE :: firstcall=.true.
!$OMP THREADPRIVATE(firstcall)


      IF (firstcall) THEN

	 write(*,*) "value for ratqs? "
         ratqs=0.2 ! default value
         call getin_p("ratqs",ratqs)
         write(*,*) " ratqs = ",ratqs

         firstcall = .false.
      ENDIF

!     GCM -----> subroutine variables, initialisation of outputs

      pdtlsc(1:ngrid,1:nlayer)  = 0.0
      pdqvaplsc(1:ngrid,1:nlayer)  = 0.0
      pdqliqlsc(1:ngrid,1:nlayer) = 0.0
      rneb(1:ngrid,1:nlayer) = 0.0
      Lcp=RLVTT/RCPD


      ! Evaporate cloud water/ice
      call evap(ngrid,nlayer,nq,ptimestep,pt,pq,pdq,pdt,dqevap,dtevap,qevap,tevap)
      ! note: we use qevap but not tevap in largescale/moistadj
            ! otherwise is a big mess


!  Boucle verticale (du haut vers le bas)
   DO k = nlayer, 1, -1

      zt(1:ngrid)=pt(1:ngrid,k)+(pdt(1:ngrid,k)+dtevap(1:ngrid,k))*ptimestep
      zq(1:ngrid)=qevap(1:ngrid,k,igcm_h2o_vap) !liquid water is included in qevap

!     Calculer la vapeur d'eau saturante et 
!     determiner la condensation partielle
      DO i = 1, ngrid

         local_p=pplay(i,k)
         if(zt(i).le.15.) then
	    print*,'in lsc',i,k,zt(i)
!	    zt(i)=15.   ! check too low temperatures
         endif
         call Psat_water(zt(i),local_p,psat_tmp,zqs_temp)
	 zqs(i)=zqs_temp
 
         zdelq(i) = MAX(MIN(ratqs * zq(i),1.-zq(i)),1.d-12)
	 rneb(i,k) = (zq(i)+zdelq(i)-zqs(i)) / (2.0*zdelq(i))
	 if (rneb(i,k).lt.0.) then  !no clouds

	    rneb(i,k)=0.
	    zcond(i)=0.

	 else if ((rneb(i,k).gt.0.99).or.(ratqs.lt.1.e-6)) then    !complete cloud cover, we start without evaporating
	    rneb(i,k)=1.
            zt(i)=pt(i,k)+pdt(i,k)*ptimestep
	    zx_q(i) = pq(i,k,igcm_h2o_vap)+pdq(i,k,igcm_h2o_vap)*ptimestep
	    dqevap(i,k)=0.
!           iterative process to stabilize the scheme when large water amounts JL12
            zcond(i) = 0.0d0
            Do nn=1,nitermax  
               call Psat_water(zt(i),local_p,psat_tmp,zqs_temp)
	       zqs(i)=zqs_temp
	       call Lcpdqsat_water(zt(i),local_p,psat_tmp,zqs_temp,zdqs,dlnpsat_tmp)
               zcond_iter = alpha*(zx_q(i)-zqs(i))/(1.d0+zdqs)	   
                  !zcond can be negative here
               zx_q(i) = zx_q(i) - zcond_iter
	       zcond(i) = zcond(i) + zcond_iter
	       zt(i) = zt(i) + zcond_iter*Lcp
	       if (ABS(zcond_iter/alpha/zqs(i)).lt.qthreshold) exit
!	       if (ABS(zcond_iter/alpha).lt.qthreshold) exit
	       if (nn.eq.nitermax) print*,'itermax in largescale'
	    End do ! niter
	    zcond(i)=MAX(zcond(i),-(pq(i,k,igcm_h2o_ice)+pdq(i,k,igcm_h2o_ice)*ptimestep))

	 else   !standard case	    
	    zx_q(i) = (zq(i)+zdelq(i)+zqs(i))/2.0d0 !water vapor in cloudy sky
!           iterative process to stabilize the scheme when large water amounts JL12
            zcond(i) = 0.0d0
            Do nn=1,nitermax 
               ! use zqs_temp and not zqs(i) to force type conversion
               ! -- might not be a good solution, actually
               ! but this is compliant with "complete cloud cover" case above 
	       call Lcpdqsat_water(zt(i),local_p,psat_tmp,zqs_temp,zdqs,dlnpsat_tmp)
               zcond_iter = MAX(0.0d0,alpha*(zx_q(i)-zqs(i))/(1.d0+zdqs))	   
                  !zcond always postive! cannot evaporate clouds!
                  !this is why we must reevaporate before largescale
               zx_q(i) = zx_q(i) - zcond_iter
	       zcond(i) = zcond(i) + zcond_iter
	       if (ABS(zcond_iter/alpha/zqs(i)).lt.qthreshold) exit
!	       if (ABS(zcond_iter/alpha).lt.qthreshold) exit
	       zt(i) = zt(i) + zcond_iter*Lcp*rneb(i,k)
               call Psat_water(zt(i),local_p,psat_tmp,zqs_temp)
	       zqs(i)=zqs_temp
	       if (nn.eq.nitermax) print*,'itermax in largescale'
	    End do ! niter

	 Endif

         zcond(i) = zcond(i)*rneb(i,k)/ptimestep ! JL12

      ENDDO

!     Tendances de t et q
         pdqvaplsc(1:ngrid,k)  = dqevap(1:ngrid,k) - zcond(1:ngrid)
         pdqliqlsc(1:ngrid,k) = - pdqvaplsc(1:ngrid,k)
         pdtlsc(1:ngrid,k)  = pdqliqlsc(1:ngrid,k)*Lcp

   Enddo ! k= nlayer, 1, -1
 

      end
