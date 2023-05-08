subroutine optci(PQMO,NLAY,PLEV,TLEV,TMID,PMID,      &
     DTAUI,TAUCUMI,COSBI,WBARI,TAUGSURF,SEASHAZEFACT)

  use radinc_h
  use radcommon_h, only: gasi,gasi_recomb,tlimit,Cmk,gzlat_ig, &
                         tgasref,pfgasref,wnoi,scalep,indi
  use gases_h
  use datafile_mod, only: haze_opt_file
  use comcstfi_mod, only: r
  use callkeys_mod, only: continuum,graybody,corrk_recombin,               &
                          callclouds,callmufi,seashaze,uncoupl_optic_haze
  use tracer_h, only : nmicro,nice

  implicit none

  !==================================================================
  !     
  !     Purpose
  !     -------
  !     Calculates longwave optical constants at each level. For each
  !     layer and spectral interval in the IR it calculates WBAR, DTAU
  !     and COSBAR. For each level it calculates TAU.
  !     
  !     TAUCUMI(L,LW) is the cumulative optical depth at level L (or alternatively
  !     at the *bottom* of layer L), LW is the spectral wavelength interval.
  !     
  !     TLEV(L) - Temperature at the layer boundary (i.e., level)
  !     PLEV(L) - Pressure at the layer boundary (i.e., level)
  !
  !     Authors
  !     -------
  !     Adapted from the NASA Ames code by R. Wordsworth (2009)
  !     Clean and adaptation to Titan by J. Vatant d'Ollone (2016-17)
  !     
  !==================================================================


  !==========================================================
  ! Input/Output
  !==========================================================
  REAL*8, INTENT(IN)  :: PQMO(nlay,nmicro)  ! Tracers for microphysics optics (X/m2).
  INTEGER, INTENT(IN) :: NLAY               ! Number of pressure layers (for pqmo)
  REAL*8, INTENT(IN)  :: PLEV(L_LEVELS), TLEV(L_LEVELS)
  REAL*8, INTENT(IN)  :: TMID(L_LEVELS), PMID(L_LEVELS)
  REAL*8, INTENT(IN)  :: SEASHAZEFACT(L_LEVELS)
  
  REAL*8, INTENT(OUT) :: DTAUI(L_NLAYRAD,L_NSPECTI,L_NGAUSS) 
  REAL*8, INTENT(OUT) :: TAUCUMI(L_LEVELS,L_NSPECTI,L_NGAUSS)
  REAL*8, INTENT(OUT) :: COSBI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
  REAL*8, INTENT(OUT) :: WBARI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
  REAL*8, INTENT(OUT) :: TAUGSURF(L_NSPECTI,L_NGAUSS-1)
  ! ==========================================================
  
  real*8 DTAUKI(L_LEVELS,L_NSPECTI,L_NGAUSS)

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
  real*8  DPR(L_LEVELS), U(L_LEVELS)
  real*8  LCOEF(4), LKCOEF(L_LEVELS,4)

  real*8 DCONT
  double precision wn_cont, p_cont, p_air, T_cont, dtemp, dtempc
  double precision p_cross

  real*8  KCOEF(4)
   
  ! temporary variable to reduce memory access time to gasi
  real*8 tmpk(2,2)
  
  ! temporary variables for multiple aerosol calculation
  real*8 atemp
  real*8 btemp(L_NLAYRAD,L_NSPECTI)

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
  real*8,save :: rhoaer_s(L_NSPECTI),ssa_s(L_NSPECTI),asf_s(L_NSPECTI)
  real*8,save :: rhoaer_f(L_NSPECTI),ssa_f(L_NSPECTI),asf_f(L_NSPECTI)
!$OMP THREADPRIVATE(rhoaer_s,rhoaer_f,ssa_s,ssa_f,asf_s,asf_f)

  logical,save :: firstcall=.true.
!$OMP THREADPRIVATE(firstcall)


  !! AS: to save time in computing continuum (see bilinearbig)
  IF (.not.ALLOCATED(indi)) THEN
      ALLOCATE(indi(L_NSPECTI,ngasmx,ngasmx))
      indi = -9999  ! this initial value means "to be calculated"
  ENDIF

  ! Some initialisation because there's a pb with disr_haze at the limits (nw=1)
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
       READ(12,*) dumch, dumwvl, rhoaer_f(nw), ssa_f(nw), asf_f(nw), rhoaer_s(nw), ssa_s(nw), asf_s(nw)
     ENDDO
     CLOSE(12)
  ENDIF

  !=======================================================================
  !     Determine the total gas opacity throughout the column, for each
  !     spectral interval, NW, and each Gauss point, NG.

  taugsurf(:,:) = 0.0
  dpr(:)        = 0.0
  lkcoef(:,:)   = 0.0

  do K=2,L_LEVELS
  
     ilay = L_NLAYRAD+1 - k/2 ! int. arithmetic => gives the gcm layer index (reversed)
     
     DPR(k) = PLEV(K)-PLEV(K-1)

     ! if we have continuum opacities, we need dz
      dz(k) = dpr(k)*R*TMID(K)/(gzlat_ig(ilay)*PMID(K))
      U(k)  = Cmk(ilay)*DPR(k)     ! only Cmk line in optci.F
      
     call tpindex(PMID(K),TMID(K),pfgasref,tgasref,LCOEF,MT(K),MP(K))

     do LK=1,4
        LKCOEF(K,LK) = LCOEF(LK)
     end do
  end do                    ! levels

  do NW=1,L_NSPECTI

     do K=2,L_LEVELS
     
        ilay = L_NLAYRAD+1 - k/2 ! int. arithmetic => gives the gcm layer index (reversed)
        
        IF (callmufi .AND. (.NOT. uncoupl_optic_haze)) THEN
          m3as = pqmo(ilay,2) / 2.0
          m3af = pqmo(ilay,4) / 2.0
          
          IF ( ilay .lt. 18 ) THEN
            m3as = pqmo(18,2) / 2.0 *dz(k)/dz(18)
            m3af = pqmo(18,4) / 2.0 *dz(k)/dz(18)
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
            WRITE(*,*) 'WARNING: In optci, optical properties &
                       &calculations are not implemented yet'
        ELSE
          ! Call fixed vertical haze profile of extinction - same for all columns
          call disr_haze(dz(k),plev(k),wnoi(nw),DHAZE_T(k,nw),SSA_T(k,nw),ASF_T(k,nw))
          if (seashaze) DHAZE_T(k,nw) = DHAZE_T(k,nw)*seashazefact(k)
        ENDIF

        DCONT = 0.0d0 ! continuum absorption

        if(continuum.and.(.not.graybody))then
           ! include continua if necessary
           wn_cont = dble(wnoi(nw))
           T_cont  = dble(TMID(k))
           do igas=1,ngasmx

              p_cont  = dble(PMID(k)*scalep*gfrac(igas,ilay))

              dtemp=0.0d0
              if(igas.eq.igas_N2)then

                 interm = indi(nw,igas,igas)
                 call interpolateN2N2(wn_cont,T_cont,p_cont,dtemp,.false.,interm)
                 indi(nw,igas,igas) = interm

              elseif(igas.eq.igas_H2)then

                 ! first do self-induced absorption
                 interm = indi(nw,igas,igas)
                 call interpolateH2H2(wn_cont,T_cont,p_cont,dtemp,.false.,interm)
                 indi(nw,igas,igas) = interm

                 ! then cross-interactions with other gases
                 do jgas=1,ngasmx
                    p_cross = dble(PMID(k)*scalep*gfrac(jgas,ilay))
                    dtempc  = 0.0d0
                    if(jgas.eq.igas_N2)then 
                       interm = indi(nw,igas,jgas)
                       call interpolateN2H2(wn_cont,T_cont,p_cross,p_cont,dtempc,.false.,interm)
                       indi(nw,igas,jgas) = interm
                    endif
                    dtemp = dtemp + dtempc
                 enddo

               elseif(igas.eq.igas_CH4)then

                 ! first do self-induced absorption
                 interm = indi(nw,igas,igas)
                 call interpolateCH4CH4(wn_cont,T_cont,p_cont,dtemp,.false.,interm)
                 indi(nw,igas,igas) = interm

                 ! then cross-interactions with other gases
                 do jgas=1,ngasmx
                    p_cross = dble(PMID(k)*scalep*gfrac(jgas,ilay))
                    dtempc  = 0.0d0
                    if(jgas.eq.igas_N2)then 
                       interm = indi(nw,igas,jgas)
                       call interpolateN2CH4(wn_cont,T_cont,p_cross,p_cont,dtempc,.false.,interm)
                       indi(nw,igas,jgas) = interm
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
             tmpk = GASI_RECOMB(MT(K):MT(K)+1,MP(K):MP(K)+1,NW,NG)
           else
             tmpk = GASI(MT(K):MT(K)+1,MP(K):MP(K)+1,1,NW,NG)
           endif

           KCOEF(1) = tmpk(1,1) ! KCOEF(1) = GASI(MT(K),MP(K),1,NW,NG)
           KCOEF(2) = tmpk(1,2) ! KCOEF(2) = GASI(MT(K),MP(K)+1,1,NW,NG)
           KCOEF(3) = tmpk(2,2) ! KCOEF(3) = GASI(MT(K)+1,MP(K)+1,1,NW,NG)
           KCOEF(4) = tmpk(2,1) ! KCOEF(4) = GASI(MT(K)+1,MP(K),1,NW,NG)


           ! Interpolate the gaseous k-coefficients to the requested T,P values

           ANS = LKCOEF(K,1)*KCOEF(1) + LKCOEF(K,2)*KCOEF(2) +            &
                LKCOEF(K,3)*KCOEF(3) + LKCOEF(K,4)*KCOEF(4)

           TAUGAS  = U(k)*ANS

           TAUGSURF(NW,NG) = TAUGSURF(NW,NG) + TAUGAS + DCONT
           DTAUKI(K,nw,ng) = TAUGAS    & 
                             + DCONT   & ! For parameterized continuum absorption
			     + DHAZE_T(K,NW)  ! For Titan haze

        end do

        ! Now fill in the "clear" part of the spectrum (NG = L_NGAUSS),
        ! which holds continuum opacity only

        NG              = L_NGAUSS
        DTAUKI(K,nw,ng) = 0.d0      & 
                          + DCONT   & ! For parameterized continuum absorption
	                  + DHAZE_T(K,NW)     ! For Titan Haze

     end do
  end do

  !=======================================================================
  !     Now the full treatment for the layers, where besides the opacity
  !     we need to calculate the scattering albedo and asymmetry factors
  ! ======================================================================

  ! Haze scattering
  DO NW=1,L_NSPECTI
    DO K=2,L_LEVELS
      DHAZES_T(K,NW) = DHAZE_T(K,NW) * SSA_T(K,NW)
    ENDDO
  ENDDO

  DO NW=1,L_NSPECTI
     DO L=1,L_NLAYRAD-1
        K              = 2*L+1
        btemp(L,NW) = DHAZES_T(K,NW) + DHAZES_T(K+1,NW)
     END DO ! L vertical loop
     
     ! Last level
     L           = L_NLAYRAD
     K           = 2*L+1
     btemp(L,NW) = DHAZES_T(K,NW)
     
  END DO                    ! NW spectral loop
  

  DO NW=1,L_NSPECTI
     NG = L_NGAUSS
     DO L=1,L_NLAYRAD-1

        K              = 2*L+1
        DTAUI(L,nw,ng) = DTAUKI(K,NW,NG) + DTAUKI(K+1,NW,NG)! + 1.e-50

        atemp = 0.
        if(DTAUI(L,NW,NG) .GT. 1.0D-9) then
           atemp = atemp +                   &
                ASF_T(K,NW)*DHAZES_T(K,NW) + &
                ASF_T(K+1,NW)*DHAZES_T(K+1,NW)

           WBARI(L,nw,ng) = btemp(L,nw)  / DTAUI(L,NW,NG)
        else
           WBARI(L,nw,ng) = 0.0D0
           DTAUI(L,NW,NG) = 1.0D-9
        endif

        if(btemp(L,nw) .GT. 0.0d0) then
           cosbi(L,NW,NG) = atemp/btemp(L,nw)
        else
           cosbi(L,NW,NG) = 0.0D0
        end if

     END DO ! L vertical loop
     
     ! Last level
     
     L              = L_NLAYRAD
     K              = 2*L+1
     DTAUI(L,nw,ng) = DTAUKI(K,NW,NG) ! + 1.e-50

     atemp = 0.
     if(DTAUI(L,NW,NG) .GT. 1.0D-9) then
        atemp = atemp + ASF_T(K,NW)*DHAZES_T(K,NW)
        WBARI(L,nw,ng) = btemp(L,nw)  / DTAUI(L,NW,NG)
     else
        WBARI(L,nw,ng) = 0.0D0
        DTAUI(L,NW,NG) = 1.0D-9
     endif

     if(btemp(L,nw) .GT. 0.0d0) then
        cosbi(L,NW,NG) = atemp/btemp(L,nw)
     else
        cosbi(L,NW,NG) = 0.0D0
     end if


     ! Now the other Gauss points, if needed.

     DO NG=1,L_NGAUSS-1
        IF(TAUGSURF(NW,NG) .gt. TLIMIT) THEN

           DO L=1,L_NLAYRAD-1
              K              = 2*L+1
              DTAUI(L,nw,ng) = DTAUKI(K,NW,NG)+DTAUKI(K+1,NW,NG)! + 1.e-50

              if(DTAUI(L,NW,NG) .GT. 1.0D-9) then

                 WBARI(L,nw,ng) = btemp(L,nw)  / DTAUI(L,NW,NG)

              else
                 WBARI(L,nw,ng) = 0.0D0
                 DTAUI(L,NW,NG) = 1.0D-9
              endif

              cosbi(L,NW,NG) = cosbi(L,NW,L_NGAUSS)
           END DO ! L vertical loop
           
           ! Last level 
           L              = L_NLAYRAD
           K              = 2*L+1
           DTAUI(L,nw,ng) = DTAUKI(K,NW,NG)! + 1.e-50

           if(DTAUI(L,NW,NG) .GT. 1.0D-9) then

              WBARI(L,nw,ng) = btemp(L,nw)  / DTAUI(L,NW,NG)

           else
              WBARI(L,nw,ng) = 0.0D0
              DTAUI(L,NW,NG) = 1.0D-9
           endif

           cosbi(L,NW,NG) = cosbi(L,NW,L_NGAUSS)
           
        END IF

     END DO                 ! NG Gauss loop
  END DO                    ! NW spectral loop

  ! Total extinction optical depths

  DO NG=1,L_NGAUSS       ! full gauss loop
     DO NW=1,L_NSPECTI       
        TAUCUMI(1,NW,NG)=0.0D0
        DO K=2,L_LEVELS
           TAUCUMI(K,NW,NG)=TAUCUMI(K-1,NW,NG)+DTAUKI(K,NW,NG)
        END DO
     END DO                 ! end full gauss loop
  END DO

  ! be aware when comparing with textbook results 
  ! (e.g. Pierrehumbert p. 218) that 
  ! taucumi does not take the <cos theta>=0.5 factor into
  ! account. It is the optical depth for a vertically 
  ! ascending ray with angle theta = 0.

  if(firstcall) firstcall = .false.

  return


end subroutine optci



