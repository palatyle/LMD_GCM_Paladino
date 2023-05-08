SUBROUTINE OPTCV(PQMO,NLAY,PLEV,TMID,PMID,  &
     DTAUV,TAUV,TAUCUMV,WBARV,COSBV,TAURAY,TAUGSURF,SEASHAZEFACT)

  use radinc_h
  use radcommon_h, only: gasv,gasv_recomb,tlimit,Cmk,gzlat_ig, &
                         tgasref,pfgasref,wnov,scalep,indv
  use gases_h
  use datafile_mod, only: haze_opt_file
  use comcstfi_mod, only: r
  use callkeys_mod, only: continuum,graybody,callgasvis,corrk_recombin,     &
                          callclouds,callmufi,seashaze,uncoupl_optic_haze
  use tracer_h, only: nmicro,nice
  use MMP_OPTICS

  implicit none

  !==================================================================
  !     
  !     Purpose
  !     -------
  !     Calculates shortwave optical constants at each level.
  !     
  !     Authors
  !     -------
  !     Adapted from the NASA Ames code by R. Wordsworth (2009)
  !     Clean and adaptation to Titan by J. Vatant d'Ollone (2016-17)
  !     
  !==================================================================
  !     
  !     THIS SUBROUTINE SETS THE OPTICAL CONSTANTS IN THE VISUAL  
  !     IT CALCULATES FOR EACH LAYER, FOR EACH SPECTRAL INTERVAL IN THE VISUAL
  !     LAYER: WBAR, DTAU, COSBAR
  !     LEVEL: TAU
  !     
  !     TAUV(L,NW,NG) is the cumulative optical depth at the top of radiation code
  !     layer L. NW is spectral wavelength interval, ng the Gauss point index.
  !     
  !     TLEV(L) - Temperature at the layer boundary
  !     PLEV(L) - Pressure at the layer boundary (i.e. level)
  !     GASV(NT,NPS,NW,NG) - Visible k-coefficients 
  !     
  !-------------------------------------------------------------------


  !==========================================================
  ! Input/Output
  !==========================================================
  REAL*8, INTENT(IN)  :: PQMO(nlay,nmicro)  ! Tracers for microphysics optics (X/m2).
  INTEGER, INTENT(IN) :: NLAY               ! Number of pressure layers (for pqmo)
  REAL*8, INTENT(IN)  :: PLEV(L_LEVELS)
  REAL*8, INTENT(IN)  :: TMID(L_LEVELS), PMID(L_LEVELS)
  REAL*8, INTENT(IN)  :: TAURAY(L_NSPECTV)
  REAL*8, INTENT(IN)  :: SEASHAZEFACT(L_LEVELS)
  
  REAL*8, INTENT(OUT) :: DTAUV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
  REAL*8, INTENT(OUT) :: TAUV(L_NLEVRAD,L_NSPECTV,L_NGAUSS)
  REAL*8, INTENT(OUT) :: TAUCUMV(L_LEVELS,L_NSPECTV,L_NGAUSS)
  REAL*8, INTENT(OUT) :: COSBV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
  REAL*8, INTENT(OUT) :: WBARV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
  REAL*8, INTENT(OUT) :: TAUGSURF(L_NSPECTV,L_NGAUSS-1)
  ! ==========================================================
  
  real*8 DTAUKV(L_LEVELS,L_NSPECTV,L_NGAUSS)

  ! Titan customisation
  ! J. Vatant d'Ollone (2016)
  real*8 DHAZE_T(L_LEVELS,L_NSPECTI)
  real*8 DHAZES_T(L_LEVELS,L_NSPECTI)
  real*8 SSA_T(L_LEVELS,L_NSPECTI)
  real*8 ASF_T(L_LEVELS,L_NSPECTI)
  ! ==========================

  integer L, NW, NG, K, LK, IAER
  integer MT(L_LEVELS), MP(L_LEVELS), NP(L_LEVELS)
  real*8  ANS, TAUGAS
  real*8  TRAY(L_LEVELS,L_NSPECTV)
  real*8  DPR(L_LEVELS), U(L_LEVELS)
  real*8  LCOEF(4), LKCOEF(L_LEVELS,4)

  real*8 DCONT
  real*8 DRAYAER
  double precision wn_cont, p_cont, p_air, T_cont, dtemp, dtempc
  double precision p_cross

  real*8  KCOEF(4)
  
  ! temporary variable to reduce memory access time to gasv
  real*8 tmpk(2,2)

  ! temporary variables for multiple aerosol calculation
  real*8 atemp(L_NLAYRAD,L_NSPECTV)
  real*8 btemp(L_NLAYRAD,L_NSPECTV)
  real*8 ctemp(L_NLAYRAD,L_NSPECTV)

  ! variables for k in units m^-1
  real*8 dz(L_LEVELS)

  integer igas, jgas, ilay

  integer interm

  ! Variables for haze optics
  character(len=200) file_path
  logical file_ok
  integer dumch
  real*8  dumwvl

  real*8 m3as,m3af
  real*8 dtauaer_s,dtauaer_f
  real*8,save :: rhoaer_s(L_NSPECTV),ssa_s(L_NSPECTV),asf_s(L_NSPECTV)
  real*8,save :: rhoaer_f(L_NSPECTV),ssa_f(L_NSPECTV),asf_f(L_NSPECTV)
!$OMP THREADPRIVATE(rhoaer_s,rhoaer_f,ssa_s,ssa_f,asf_s,asf_f)
  
  logical,save :: firstcall=.true.
!$OMP THREADPRIVATE(firstcall)


  !! AS: to save time in computing continuum (see bilinearbig)
  IF (.not.ALLOCATED(indv)) THEN
      ALLOCATE(indv(L_NSPECTV,ngasmx,ngasmx))
      indv = -9999 ! this initial value means "to be calculated"
  ENDIF
  
  ! Some initialisation beacause there's a pb with disr_haze at the limits (nw=1)
  ! I should check this - For now we set vars to zero : better than nans - JVO 2017
  DHAZE_T(:,:) = 0.0
  SSA_T(:,:)   = 0.0
  ASF_T(:,:)   = 0.0
  
  ! Load tabulated haze optical properties if needed.
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  IF (firstcall .AND. callmufi .AND. (.NOT. uncoupl_optic_haze)) THEN
     OPEN(12,file=TRIM(haze_opt_file),form='formatted') ! The file has been inquired in physiq_mod firstcall
     READ(12,*) ! dummy header
     DO NW=1,L_NSPECTI
       READ(12,*) ! there's IR 1st
     ENDDO
     DO NW=1,L_NSPECTV
       READ(12,*) dumch, dumwvl, rhoaer_f(nw), ssa_f(nw), asf_f(nw), rhoaer_s(nw), ssa_s(nw), asf_s(nw)
     ENDDO
     CLOSE(12)
  ENDIF

  !=======================================================================
  !     Determine the total gas opacity throughout the column, for each
  !     spectral interval, NW, and each Gauss point, NG.
  !     Calculate the continuum opacities, i.e., those that do not depend on
  !     NG, the Gauss index.

  taugsurf(:,:) = 0.0
  dpr(:)        = 0.0
  lkcoef(:,:)   = 0.0

  do K=2,L_LEVELS
  
     ilay = L_NLAYRAD+1 - k/2 ! int. arithmetic => gives the gcm layer index (reversed)
  
     DPR(k) = PLEV(K)-PLEV(K-1)

     ! if we have continuum opacities, we need dz

      dz(k) = dpr(k)*R*TMID(K)/(gzlat_ig(ilay)*PMID(K))
      U(k)  = Cmk(ilay)*DPR(k)     ! only Cmk line in optcv.F     

     call tpindex(PMID(K),TMID(K),pfgasref,tgasref,LCOEF,MT(K),MP(K))

     do LK=1,4
        LKCOEF(K,LK) = LCOEF(LK)
     end do
  end do                    ! levels

  ! Rayleigh scattering
  do NW=1,L_NSPECTV
     TRAY(1:4,NW)   = 1.d-30
     do K=5,L_LEVELS
        TRAY(K,NW)   = TAURAY(NW) * DPR(K)
     end do                    ! levels
  end do
  
  !     we ignore K=1...
  do K=2,L_LEVELS
  
     ilay = L_NLAYRAD+1 - k/2 ! int. arithmetic => gives the gcm layer index (reversed)

     do NW=1,L_NSPECTV
     
        IF (callmufi .AND. (.NOT. uncoupl_optic_haze)) THEN
          m3as = pqmo(ilay,2) / 2.0
          m3af = pqmo(ilay,4) / 2.0
          
          IF ( ilay .lt. 18 ) THEN
            m3as = pqmo(18,2) / 2.0 * dz(k) / dz(18)
            m3af = pqmo(18,4) / 2.0 * dz(k) / dz(18)
          ENDIF

          dtauaer_s     = m3as*rhoaer_s(nw)
          dtauaer_f     = m3af*rhoaer_f(nw)
          DHAZE_T(k,nw) = dtauaer_s + dtauaer_f

          IF ( dtauaer_s + dtauaer_f .GT. 1.D-30 ) THEN
            SSA_T(k,nw)   = ( dtauaer_s*ssa_s(nw) + dtauaer_f*ssa_f(nw) ) / ( dtauaer_s+dtauaer_f )
            ASF_T(k,nw)   = ( dtauaer_s*ssa_s(nw)*asf_s(nw) + dtauaer_f*ssa_f(nw)*asf_f(nw) )  &
                            / ( ssa_s(nw)*dtauaer_s + ssa_f(nw)*dtauaer_f )
          ELSE
             DHAZE_T(k,nw) = 0.D0
             SSA_T(k,nw)   = 1.0
             ASF_T(k,nw)   = 1.0
          ENDIF
          
          IF (callclouds.and.firstcall) & 
            WRITE(*,*) 'WARNING: In optcv, optical properties &
                       &calculations are not implemented yet'
        ELSE
          ! Call fixed vertical haze profile of extinction - same for all columns
          call disr_haze(dz(k),plev(k),wnov(nw),DHAZE_T(k,nw),SSA_T(k,nw),ASF_T(k,nw))
          if (seashaze) DHAZE_T(k,nw) = DHAZE_T(k,nw)*seashazefact(k)
        ENDIF
        
        !JL18 It seems to be good to have aerosols in the first "radiative layer" of the gcm in the IR
        !   but visible does not handle very well diffusion in first layer.
        !   The tauaero and tauray are thus set to 0 (a small value for rayleigh because the code crashes otherwise)
        !   in the 4 first semilayers in optcv, but not optci.
        !   This solves random variations of the sw heating at the model top. 
        if (k<5)  DHAZE_T(K,:) = 0.0
         
        DRAYAER = TRAY(K,NW)
        !     DRAYAER is Tau RAYleigh scattering, plus AERosol opacity
        DRAYAER = DRAYAER + DHAZE_T(K,NW) ! Titan's aerosol

        DCONT = 0.0 ! continuum absorption

        if(continuum.and.(.not.graybody).and.callgasvis)then
           ! include continua if necessary
           wn_cont = dble(wnov(nw))
           T_cont  = dble(TMID(k))
           do igas=1,ngasmx

              p_cont  = dble(PMID(k)*scalep*gfrac(igas,ilay))

              dtemp=0.0
              if(igas.eq.igas_N2)then

                 interm = indv(nw,igas,igas)
!                 call interpolateN2N2(wn_cont,T_cont,p_cont,dtemp,.false.,interm)
                 indv(nw,igas,igas) = interm
                 ! only goes to 500 cm^-1, so unless we're around a cold brown dwarf, this is irrelevant in the visible

              elseif(igas.eq.igas_H2)then

                 ! first do self-induced absorption
                 interm = indv(nw,igas,igas)
                 call interpolateH2H2(wn_cont,T_cont,p_cont,dtemp,.false.,interm)
                 indv(nw,igas,igas) = interm

                 ! then cross-interactions with other gases
                 do jgas=1,ngasmx
                    p_cross = dble(PMID(k)*scalep*gfrac(jgas,ilay))
                    dtempc  = 0.0
                    if(jgas.eq.igas_N2)then 
                       interm = indv(nw,igas,jgas)
                       call interpolateN2H2(wn_cont,T_cont,p_cross,p_cont,dtempc,.false.,interm)
                       indv(nw,igas,jgas) = interm
                       ! should be irrelevant in the visible
                    endif
                    dtemp = dtemp + dtempc
                 enddo

               elseif(igas.eq.igas_CH4)then

                 ! first do self-induced absorption
                 interm = indv(nw,igas,igas)
                 call interpolateCH4CH4(wn_cont,T_cont,p_cont,dtemp,.false.,interm)
                 indv(nw,igas,igas) = interm

                 ! then cross-interactions with other gases
                 do jgas=1,ngasmx
                    p_cross = dble(PMID(k)*scalep*gfrac(jgas,ilay))
                    dtempc  = 0.0
                    if(jgas.eq.igas_N2)then 
                       interm = indv(nw,igas,jgas)
                       call interpolateN2CH4(wn_cont,T_cont,p_cross,p_cont,dtempc,.false.,interm)
                       indv(nw,igas,jgas) = interm
                    endif
                    dtemp = dtemp + dtempc
                 enddo

              endif

              DCONT = DCONT + dtemp

           enddo

           DCONT = DCONT*dz(k)

        endif

        do ng=1,L_NGAUSS-1

           ! Now compute TAUGAS

           ! JVO 2017 : added tmpk because the repeated calls to gasi/v increased dramatically
           ! the execution time of optci/v -> ~ factor 2 on the whole radiative
           ! transfer on the tested simulations !

           if (corrk_recombin) then 
             tmpk = GASV_RECOMB(MT(K):MT(K)+1,MP(K):MP(K)+1,NW,NG)
           else
             tmpk = GASV(MT(K):MT(K)+1,MP(K):MP(K)+1,1,NW,NG)
           endif
              
           KCOEF(1) = tmpk(1,1) ! KCOEF(1) = GASV(MT(K),MP(K),1,NW,NG)
           KCOEF(2) = tmpk(1,2) ! KCOEF(2) = GASV(MT(K),MP(K)+1,1,NW,NG)
           KCOEF(3) = tmpk(2,2) ! KCOEF(3) = GASV(MT(K)+1,MP(K)+1,1,NW,NG)
           KCOEF(4) = tmpk(2,1) ! KCOEF(4) = GASV(MT(K)+1,MP(K),1,NW,NG)

           ! Interpolate the gaseous k-coefficients to the requested T,P values

           ANS = LKCOEF(K,1)*KCOEF(1) + LKCOEF(K,2)*KCOEF(2) +            &
                LKCOEF(K,3)*KCOEF(3) + LKCOEF(K,4)*KCOEF(4)


           TAUGAS  = U(k)*ANS

           TAUGSURF(NW,NG) = TAUGSURF(NW,NG) + TAUGAS + DCONT
           DTAUKV(K,nw,ng) = TAUGAS & 
                             + DRAYAER & ! DRAYAER includes all scattering contributions
                             + DCONT ! For parameterized continuum aborption

        end do

        ! Now fill in the "clear" part of the spectrum (NG = L_NGAUSS),
        ! which holds continuum opacity only

        NG              = L_NGAUSS
        DTAUKV(K,nw,ng) = DRAYAER + DCONT ! Scattering + parameterized continuum absorption, including Titan's haze

     end do
  end do


  !=======================================================================
  !     Now the full treatment for the layers, where besides the opacity
  !     we need to calculate the scattering albedo and asymmetry factors
  ! ======================================================================

  ! Haze scattering
            !JL18 It seems to be good to have aerosols in the first "radiative layer" of the gcm in the IR
            !   but not in the visible
            !   The dhaze_s is thus set to 0 in the 4 first semilayers in optcv, but not optci.
            !   This solves random variations of the sw heating at the model top. 
  DO NW=1,L_NSPECTV
    DHAZES_T(1:4,NW) = 0.d0
    DO K=5,L_LEVELS
      DHAZES_T(K,NW) = DHAZE_T(K,NW) * SSA_T(K,NW) ! effect of scattering albedo on haze
    ENDDO
  ENDDO


  DO NW=1,L_NSPECTV
     DO L=1,L_NLAYRAD-1
        K              = 2*L+1
	atemp(L,NW) = ASF_T(K,NW)*DHAZES_T(K,NW) + ASF_T(K+1,NW)*DHAZES_T(K+1,NW)
        btemp(L,NW) = DHAZES_T(K,NW) + DHAZES_T(K+1,NW)
	ctemp(L,NW) = btemp(L,NW) + 0.9999*(TRAY(K,NW) + TRAY(K+1,NW)) ! JVO 2017 : does this 0.999 is really meaningful ?
	btemp(L,NW) = btemp(L,NW) + TRAY(K,NW) + TRAY(K+1,NW)
	COSBV(L,NW,1:L_NGAUSS) = atemp(L,NW)/btemp(L,NW)
     END DO ! L vertical loop
     
     ! Last level
     L           = L_NLAYRAD
     K           = 2*L+1
     atemp(L,NW) = ASF_T(K,NW)*DHAZES_T(K,NW)
     btemp(L,NW) = DHAZES_T(K,NW)
     ctemp(L,NW) = btemp(L,NW) + 0.9999*TRAY(K,NW) ! JVO 2017 : does this 0.999 is really meaningful ?
     btemp(L,NW) = btemp(L,NW) + TRAY(K,NW)
     COSBV(L,NW,1:L_NGAUSS) = atemp(L,NW)/btemp(L,NW)
     
     
  END DO                    ! NW spectral loop

  DO NG=1,L_NGAUSS
    DO NW=1,L_NSPECTV
     DO L=1,L_NLAYRAD-1

        K              = 2*L+1
        DTAUV(L,nw,ng) = DTAUKV(K,NW,NG) + DTAUKV(K+1,NW,NG)
        WBARV(L,nw,ng) = ctemp(L,NW) / DTAUV(L,nw,ng)

      END DO ! L vertical loop

        ! Last level

        L              = L_NLAYRAD
        K              = 2*L+1
	DTAUV(L,nw,ng) = DTAUKV(K,NW,NG)

        WBARV(L,NW,NG) = ctemp(L,NW) / DTAUV(L,NW,NG)

     END DO                 ! NW spectral loop
  END DO                    ! NG Gauss loop

  ! Total extinction optical depths

  DO NG=1,L_NGAUSS       ! full gauss loop
     DO NW=1,L_NSPECTV       
        TAUCUMV(1,NW,NG)=0.0D0
        DO K=2,L_LEVELS
           TAUCUMV(K,NW,NG)=TAUCUMV(K-1,NW,NG)+DTAUKV(K,NW,NG)
        END DO

        DO L=1,L_NLAYRAD
           TAUV(L,NW,NG)=TAUCUMV(2*L,NW,NG)
        END DO
        TAUV(L,NW,NG)=TAUCUMV(2*L_NLAYRAD+1,NW,NG)
     END DO            
  END DO                 ! end full gauss loop

  if(firstcall) firstcall = .false.

  return


end subroutine optcv
