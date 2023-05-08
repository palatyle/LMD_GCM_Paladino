MODULE optcv_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE OPTCV(DTAUV,TAUV,TAUCUMV,PLEV,  &
     QXVAER,QSVAER,GVAER,WBARV,COSBV,       &
     TAURAY,TAUAERO,TMID,PMID,TAUGSURF,QVAR,MUVAR)

  use radinc_h, only: L_NLAYRAD, L_NLEVRAD, L_LEVELS, L_NSPECTV, L_NGAUSS, L_REFVAR, NAERKIND
  use radcommon_h, only: gasv, tlimit, wrefVAR, Cmk, tgasref, pfgasref,wnov,scalep,indv,glat_ig
  use gases_h, only: gfrac, ngasmx, igas_H2, igas_H2O, igas_He, igas_N2
  use comcstfi_mod, only: g, r, mugaz
  use callkeys_mod, only: kastprof,continuum,graybody,H2Ocont_simple,callgasvis

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


  real*8 DTAUV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
  real*8 DTAUKV(L_LEVELS,L_NSPECTV,L_NGAUSS)
  real*8 TAUV(L_NLEVRAD,L_NSPECTV,L_NGAUSS)
  real*8 TAUCUMV(L_LEVELS,L_NSPECTV,L_NGAUSS)
  real*8 PLEV(L_LEVELS)
  real*8 TMID(L_LEVELS), PMID(L_LEVELS)
  real*8 COSBV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
  real*8 WBARV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)

  ! for aerosols
  real*8  QXVAER(L_LEVELS,L_NSPECTV,NAERKIND)
  real*8  QSVAER(L_LEVELS,L_NSPECTV,NAERKIND)
  real*8  GVAER(L_LEVELS,L_NSPECTV,NAERKIND)
  real*8  TAUAERO(L_LEVELS,NAERKIND)
  real*8  TAUAEROLK(L_LEVELS,L_NSPECTV,NAERKIND)
  real*8  TAEROS(L_LEVELS,L_NSPECTV,NAERKIND)

  integer L, NW, NG, K, LK, IAER
  integer MT(L_LEVELS), MP(L_LEVELS), NP(L_LEVELS)
  real*8  ANS, TAUGAS
  real*8  TAURAY(L_NSPECTV)
  real*8  TRAY(L_LEVELS,L_NSPECTV)
  real*8  DPR(L_LEVELS), U(L_LEVELS)
  real*8  LCOEF(4), LKCOEF(L_LEVELS,4)

  real*8 taugsurf(L_NSPECTV,L_NGAUSS-1)
  real*8 DCONT,DAERO
  real*8 DRAYAER
  double precision wn_cont, p_cont, p_air, T_cont, dtemp, dtempc
  double precision p_cross

  ! variable species mixing ratio variables
  real*8  QVAR(L_LEVELS), WRATIO(L_LEVELS), MUVAR(L_LEVELS)
  real*8  KCOEF(4)
  integer NVAR(L_LEVELS)
  
  ! temporary variables to reduce memory access time to gasv
  real*8 tmpk(2,2)
  real*8 tmpkvar(2,2,2)

  ! temporary variables for multiple aerosol calculation
  real*8 atemp(L_NLAYRAD,L_NSPECTV)
  real*8 btemp(L_NLAYRAD,L_NSPECTV)
  real*8 ctemp(L_NLAYRAD,L_NSPECTV)

  ! variables for k in units m^-1
  real*8 dz(L_LEVELS)


  integer igas, jgas

  integer interm

  !! AS: to save time in computing continuum (see bilinearbig)
  IF (.not.ALLOCATED(indv)) THEN
      ALLOCATE(indv(L_NSPECTV,ngasmx,ngasmx))
      indv = -9999 ! this initial value means "to be calculated"
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
     DPR(k) = PLEV(K)-PLEV(K-1)

     ! if we have continuum opacities, we need dz
     if(kastprof)then
        dz(k) = dpr(k)*(1000.0d0*8.3145d0/muvar(k))*TMID(K)/(g*PMID(K))
        U(k)  = Cmk*DPR(k)*mugaz/muvar(k) 
     else
        dz(k) = dpr(k)*R*TMID(K)/(glat_ig*PMID(K))*mugaz/muvar(k)
        U(k)  = Cmk*DPR(k)*mugaz/muvar(k)     ! only Cmk line in optci.F  
	    !JL13 the mugaz/muvar factor takes into account water meanmolecular weight if water is present
     endif

     call tpindex(PMID(K),TMID(K),QVAR(K),pfgasref,tgasref,WREFVAR, &
          LCOEF,MT(K),MP(K),NVAR(K),WRATIO(K))

     do LK=1,4
        LKCOEF(K,LK) = LCOEF(LK)
     end do
  end do                    ! levels

  ! Spectral dependance of aerosol absorption
            !JL18 It seems to be good to have aerosols in the first "radiative layer" of the gcm in the IR
	    !   but visible does not handle very well diffusion in first layer.
	    !   The tauaero and tauray are thus set to 0 (a small value for rayleigh because the code crashes otherwise)
	    !   in the 4 first semilayers in optcv, but not optci.
	    !   This solves random variations of the sw heating at the model top. 
  do iaer=1,naerkind
     do NW=1,L_NSPECTV
        TAEROS(1:4,NW,IAER)=0.d0
        do K=5,L_LEVELS
           TAEROS(K,NW,IAER) = TAUAERO(K,IAER) * QXVAER(K,NW,IAER)
        end do                    ! levels
     end do
  end do
  
  ! Rayleigh scattering
  do NW=1,L_NSPECTV
     TRAY(1:4,NW)   = 1d-30
     do K=5,L_LEVELS
        TRAY(K,NW)   = TAURAY(NW) * DPR(K)
     end do                    ! levels
  end do
  
  !     we ignore K=1...
  do K=2,L_LEVELS

     do NW=1,L_NSPECTV

        DRAYAER = TRAY(K,NW)
        !     DRAYAER is Tau RAYleigh scattering, plus AERosol opacity
        do iaer=1,naerkind
           DRAYAER = DRAYAER + TAEROS(K,NW,IAER)
        end do

        DCONT = 0.0 ! continuum absorption

        if(continuum.and.(.not.graybody).and.callgasvis)then
           ! include continua if necessary
           wn_cont = dble(wnov(nw))
           T_cont  = dble(TMID(k))
           do igas=1,ngasmx

              if(gfrac(igas).eq.-1)then ! variable
                 p_cont  = dble(PMID(k)*scalep*QVAR(k)) ! qvar = mol/mol
              else
                 p_cont  = dble(PMID(k)*scalep*gfrac(igas)*(1.-QVAR(k)))
              endif

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
                    p_cross = dble(PMID(k)*scalep*gfrac(jgas)*(1.-QVAR(k)))
                    dtempc  = 0.0
                    if(jgas.eq.igas_N2)then 
                       interm = indv(nw,igas,jgas)
                       call interpolateN2H2(wn_cont,T_cont,p_cross,p_cont,dtempc,.false.,interm)
                       indv(nw,igas,jgas) = interm
                       ! should be irrelevant in the visible
                    elseif(jgas.eq.igas_He)then 
                       interm = indv(nw,igas,jgas)
                       call interpolateH2He(wn_cont,T_cont,p_cross,p_cont,dtempc,.false.,interm)
                       indv(nw,igas,jgas) = interm
                    endif
                    dtemp = dtemp + dtempc
                 enddo

              elseif(igas.eq.igas_H2O.and.T_cont.gt.200.0)then

                 p_air = dble(PMID(k)*scalep) - p_cont ! note assumes background is air!
                 if(H2Ocont_simple)then
                    call interpolateH2Ocont_PPC(wn_cont,T_cont,p_cont,p_air,dtemp,.false.)
                 else
                    interm = indv(nw,igas,igas)
                    call interpolateH2Ocont_CKD(wn_cont,T_cont,p_cont,p_air,dtemp,.false.,interm)
                    indv(nw,igas,igas) = interm
                 endif

              endif

              DCONT = DCONT + dtemp

           enddo

           DCONT = DCONT*dz(k)

        endif

        do ng=1,L_NGAUSS-1

           ! Now compute TAUGAS

           ! Interpolate between water mixing ratios
           ! WRATIO = 0.0 if the requested water amount is equal to, or outside the
           ! the water data range

           if(L_REFVAR.eq.1)then ! added by RW for special no variable case
           
              ! JVO 2017 : added tmpk because the repeated calls to gasi/v increased dramatically
              ! the execution time of optci/v -> ~ factor 2 on the whole radiative
              ! transfer on the tested simulations !

              tmpk = GASV(MT(K):MT(K)+1,MP(K):MP(K)+1,1,NW,NG)
              
              KCOEF(1) = tmpk(1,1) ! KCOEF(1) = GASV(MT(K),MP(K),1,NW,NG)
              KCOEF(2) = tmpk(1,2) ! KCOEF(2) = GASV(MT(K),MP(K)+1,1,NW,NG)
              KCOEF(3) = tmpk(2,2) ! KCOEF(3) = GASV(MT(K)+1,MP(K)+1,1,NW,NG)
              KCOEF(4) = tmpk(2,1) ! KCOEF(4) = GASV(MT(K)+1,MP(K),1,NW,NG)

           else

              tmpkvar = GASV(MT(K):MT(K)+1,MP(K):MP(K)+1,NVAR(K):NVAR(K)+1,NW,NG)

              KCOEF(1) = tmpkvar(1,1,1) + WRATIO(K) *  &
                        ( tmpkvar(1,1,2)-tmpkvar(1,1,1) )

              KCOEF(2) = tmpkvar(1,2,1) + WRATIO(K) *  &
                        ( tmpkvar(1,2,2)-tmpkvar(1,2,1) )

              KCOEF(3) = tmpkvar(2,2,1) + WRATIO(K) *  &
                        ( tmpkvar(2,2,2)-tmpkvar(2,2,1) )
              
              KCOEF(4) = tmpkvar(2,1,1) + WRATIO(K) *  &
                        ( tmpkvar(2,1,2)-tmpkvar(2,1,1) )


           endif

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
        DTAUKV(K,nw,ng) = DRAYAER + DCONT ! Scattering + parameterized continuum absorption

     end do
  end do


  !=======================================================================
  !     Now the full treatment for the layers, where besides the opacity
  !     we need to calculate the scattering albedo and asymmetry factors

            !JL18 It seems to be good to have aerosols in the first "radiative layer" of the gcm in the IR
	    !   but not in the visible
	    !   The tauaero is thus set to 0 in the 4 first semilayers in optcv, but not optci.
	    !   This solves random variations of the sw heating at the model top. 
  do iaer=1,naerkind
    DO NW=1,L_NSPECTV
      TAUAEROLK(1:4,NW,IAER)=0.d0
      DO K=5,L_LEVELS
           TAUAEROLK(K,NW,IAER) = TAUAERO(K,IAER) * QSVAER(K,NW,IAER) ! effect of scattering albedo
      ENDDO
    ENDDO
  end do

  DO NW=1,L_NSPECTV
     DO L=1,L_NLAYRAD-1
        K              = 2*L+1
	atemp(L,NW) = SUM(GVAER(K,NW,1:naerkind) * TAUAEROLK(K,NW,1:naerkind))+SUM(GVAER(K+1,NW,1:naerkind) * TAUAEROLK(K+1,NW,1:naerkind))
        btemp(L,NW) = SUM(TAUAEROLK(K,NW,1:naerkind)) + SUM(TAUAEROLK(K+1,NW,1:naerkind))
	ctemp(L,NW) = btemp(L,NW) + 0.9999*(TRAY(K,NW) + TRAY(K+1,NW))  ! JVO 2017 : does this 0.999 is really meaningful ?
	btemp(L,NW) = btemp(L,NW) + TRAY(K,NW) + TRAY(K+1,NW)
	COSBV(L,NW,1:L_NGAUSS) = atemp(L,NW)/btemp(L,NW)
     END DO ! L vertical loop
     
     ! Last level
     L           = L_NLAYRAD
     K           = 2*L+1
     atemp(L,NW) = SUM(GVAER(K,NW,1:naerkind) * TAUAEROLK(K,NW,1:naerkind))
     btemp(L,NW) = SUM(TAUAEROLK(K,NW,1:naerkind))
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




end subroutine optcv

END MODULE optcv_mod
