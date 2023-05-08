MODULE optci_mod

IMPLICIT NONE

CONTAINS

subroutine optci(PLEV,TLEV,DTAUI,TAUCUMI,      &
     QXIAER,QSIAER,GIAER,COSBI,WBARI,TAUAERO,  &
     TMID,PMID,TAUGSURF,QVAR,MUVAR)

  use radinc_h, only: L_LEVELS, L_NLAYRAD, L_NSPECTI, L_NGAUSS, &
                      L_NLEVRAD, L_REFVAR, naerkind
  use radcommon_h, only: gasi,tlimit,wrefVAR,Cmk,tgasref,pfgasref,wnoi,scalep,indi,glat_ig
  use gases_h, only: gfrac, ngasmx, igas_N2, igas_He, igas_H2O, igas_H2
  use comcstfi_mod, only: g, r, mugaz
  use callkeys_mod, only: kastprof,continuum,graybody,H2Ocont_simple
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
  !     
  !==================================================================


  real*8 DTAUI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
  real*8 DTAUKI(L_LEVELS,L_NSPECTI,L_NGAUSS)
  real*8 TAUI(L_NLEVRAD,L_NSPECTI,L_NGAUSS)
  real*8 TAUCUMI(L_LEVELS,L_NSPECTI,L_NGAUSS)
  real*8 PLEV(L_LEVELS)
  real*8 TLEV(L_LEVELS)
  real*8 TMID(L_LEVELS), PMID(L_LEVELS)
  real*8 COSBI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
  real*8 WBARI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)

  ! for aerosols
  real*8  QXIAER(L_LEVELS,L_NSPECTI,NAERKIND)
  real*8  QSIAER(L_LEVELS,L_NSPECTI,NAERKIND)
  real*8  GIAER(L_LEVELS,L_NSPECTI,NAERKIND)
  real*8  TAUAERO(L_LEVELS,NAERKIND)
  real*8  TAUAEROLK(L_LEVELS,L_NSPECTI,NAERKIND)
  real*8  TAEROS(L_LEVELS,L_NSPECTI,NAERKIND)

  integer L, NW, NG, K, LK, IAER
  integer MT(L_LEVELS), MP(L_LEVELS), NP(L_LEVELS)
  real*8  ANS, TAUGAS
  real*8  DPR(L_LEVELS), U(L_LEVELS)
  real*8  LCOEF(4), LKCOEF(L_LEVELS,4)

  real*8 taugsurf(L_NSPECTI,L_NGAUSS-1)
  real*8 DCONT,DAERO
  double precision wn_cont, p_cont, p_air, T_cont, dtemp, dtempc
  double precision p_cross

  ! variable species mixing ratio variables
  real*8  QVAR(L_LEVELS), WRATIO(L_LEVELS), MUVAR(L_LEVELS)
  real*8  KCOEF(4)
  integer NVAR(L_LEVELS)
  
  ! temporary variables to reduce memory access time to gasi
  real*8 tmpk(2,2)
  real*8 tmpkvar(2,2,2)

  ! temporary variables for multiple aerosol calculation
  real*8 atemp
  real*8 btemp(L_NLAYRAD,L_NSPECTI)

  ! variables for k in units m^-1
  real*8 dz(L_LEVELS)
  !real*8 rho !! see test below

  integer igas, jgas

  integer interm

  !--- Kasting's CIA ----------------------------------------
  !real*8, parameter :: Ci(L_NSPECTI)=[                         &
  !     3.8E-5, 1.2E-5, 2.8E-6, 7.6E-7, 4.5E-7, 2.3E-7,    &
  !     5.4E-7, 1.6E-6, 0.0,                               &
  !     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,            & 
  !     0.0, 4.0E-7, 4.0E-6, 1.4E-5,    &
  !     1.0E-5, 1.2E-6, 2.0E-7, 5.0E-8, 3.0E-8, 0.0 ] 
  !real*8, parameter :: Ti(L_NSPECTI)=[ -2.2, -1.9,             &
  !     -1.7, -1.7, -1.7, -1.7, -1.7, -1.7,                &
  !     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
  !     -1.7,-1.7,-1.7,-1.7,-1.7,-1.7,-1.7, -1.7,0.0 ]
  !----------------------------------------------------------

  !! AS: to save time in computing continuum (see bilinearbig)
  IF (.not.ALLOCATED(indi)) THEN
      ALLOCATE(indi(L_NSPECTI,ngasmx,ngasmx))
      indi = -9999  ! this initial value means "to be calculated"
  ENDIF

  !=======================================================================
  !     Determine the total gas opacity throughout the column, for each
  !     spectral interval, NW, and each Gauss point, NG.

  taugsurf(:,:) = 0.0
  dpr(:)        = 0.0
  lkcoef(:,:)   = 0.0

  do K=2,L_LEVELS
     DPR(k) = PLEV(K)-PLEV(K-1)

     !--- Kasting's CIA ----------------------------------------
     !dz(k)=dpr(k)*189.02*TMID(K)/(0.03720*PMID(K))
     ! this is CO2 path length (in cm) as written by Francois
     ! delta_z = delta_p * R_specific * T / (g * P)
     ! But Kasting states that W is in units of _atmosphere_ cm
     ! So we do
     !dz(k)=dz(k)*(PMID(K)/1013.25)
     !dz(k)=dz(k)/100.0 ! in m for SI calc
     !----------------------------------------------------------

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
  do iaer=1,naerkind
     DO NW=1,L_NSPECTI
        do K=2,L_LEVELS
           TAEROS(K,NW,IAER) = TAUAERO(K,IAER) * QXIAER(K,NW,IAER)
        end do                    ! levels
     END DO
  end do

  do NW=1,L_NSPECTI

     do K=2,L_LEVELS
     
     	DAERO=SUM(TAEROS(K,NW,1:naerkind)) ! aerosol absorption

        DCONT = 0.0d0 ! continuum absorption

        if(continuum.and.(.not.graybody))then
           ! include continua if necessary
           wn_cont = dble(wnoi(nw))
           T_cont  = dble(TMID(k))
           do igas=1,ngasmx

              if(gfrac(igas).eq.-1)then ! variable
                 p_cont  = dble(PMID(k)*scalep*QVAR(k)) ! qvar = mol/mol
              else
                 p_cont  = dble(PMID(k)*scalep*gfrac(igas)*(1.-QVAR(k)))
              endif

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
                    p_cross = dble(PMID(k)*scalep*gfrac(jgas)*(1.-QVAR(k)))
                    dtempc  = 0.0d0
                    if(jgas.eq.igas_N2)then 
                       interm = indi(nw,igas,jgas)
                       call interpolateN2H2(wn_cont,T_cont,p_cross,p_cont,dtempc,.false.,interm)
                       indi(nw,igas,jgas) = interm
                    elseif(jgas.eq.igas_He)then 
                       interm = indi(nw,igas,jgas)
                       call interpolateH2He(wn_cont,T_cont,p_cross,p_cont,dtempc,.false.,interm)
                       indi(nw,igas,jgas) = interm
                    endif
                    dtemp = dtemp + dtempc
                 enddo

              elseif(igas.eq.igas_H2O.and.T_cont.gt.200.0)then

                 p_air = dble(PMID(k)*scalep) - p_cont ! note assumes background is air!
                 if(H2Ocont_simple)then
                    call interpolateH2Ocont_PPC(wn_cont,T_cont,p_cont,p_air,dtemp,.false.)
                 else
                    interm = indi(nw,igas,igas)
                    call interpolateH2Ocont_CKD(wn_cont,T_cont,p_cont,p_air,dtemp,.false.,interm)
                    indi(nw,igas,igas) = interm
                 endif

              endif

              DCONT = DCONT + dtemp

           enddo

           ! Oobleck test
           !rho = PMID(k)*scalep / (TMID(k)*286.99)
           !if(WNOI(nw).gt.300.0 .and. WNOI(nw).lt.500.0)then
           !   DCONT = rho * 0.125 * 4.6e-4
           !elseif(WNOI(nw).gt.500.0 .and. WNOI(nw).lt.700.0)then
           !   DCONT = 1000*dpr(k) * 1.0 * 4.6e-4 / g
           !   DCONT = rho * 1.0 * 4.6e-4
           !elseif(WNOI(nw).gt.700.0 .and. WNOI(nw).lt.900.0)then
           !   DCONT = rho * 0.125 * 4.6e-4
           !endif

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

              tmpk = GASI(MT(K):MT(K)+1,MP(K):MP(K)+1,1,NW,NG)

              KCOEF(1) = tmpk(1,1) ! KCOEF(1) = GASI(MT(K),MP(K),1,NW,NG)
              KCOEF(2) = tmpk(1,2) ! KCOEF(2) = GASI(MT(K),MP(K)+1,1,NW,NG)
              KCOEF(3) = tmpk(2,2) ! KCOEF(3) = GASI(MT(K)+1,MP(K)+1,1,NW,NG)
              KCOEF(4) = tmpk(2,1) ! KCOEF(4) = GASI(MT(K)+1,MP(K),1,NW,NG)

           else

              tmpkvar = GASI(MT(K):MT(K)+1,MP(K):MP(K)+1,NVAR(K):NVAR(K)+1,NW,NG)

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
           DTAUKI(K,nw,ng) = TAUGAS    & 
                             + DCONT   & ! For parameterized continuum absorption
			     + DAERO     ! For aerosol absorption

        end do

        ! Now fill in the "clear" part of the spectrum (NG = L_NGAUSS),
        ! which holds continuum opacity only

        NG              = L_NGAUSS
        DTAUKI(K,nw,ng) = 0.d0      & 
                          + DCONT   & ! For parameterized continuum absorption
	                  + DAERO     ! For aerosol absorption

     end do
  end do

  !=======================================================================
  !     Now the full treatment for the layers, where besides the opacity
  !     we need to calculate the scattering albedo and asymmetry factors

  do iaer=1,naerkind
    DO NW=1,L_NSPECTI
     DO K=2,L_LEVELS
           TAUAEROLK(K,NW,IAER) = TAUAERO(K,IAER)*QSIAER(K,NW,IAER) ! effect of scattering albedo
     ENDDO
    ENDDO
  end do
  
  DO NW=1,L_NSPECTI
     DO L=1,L_NLAYRAD-1
        K              = 2*L+1
        btemp(L,NW) = SUM(TAUAEROLK(K,NW,1:naerkind)) + SUM(TAUAEROLK(K+1,NW,1:naerkind))
     END DO ! L vertical loop
     
     ! Last level
     L           = L_NLAYRAD
     K           = 2*L+1    
     btemp(L,NW) = SUM(TAUAEROLK(K,NW,1:naerkind))
     
  END DO                    ! NW spectral loop
  

  DO NW=1,L_NSPECTI
     NG = L_NGAUSS
     DO L=1,L_NLAYRAD-1

        K              = 2*L+1
        DTAUI(L,nw,ng) = DTAUKI(K,NW,NG) + DTAUKI(K+1,NW,NG)! + 1.e-50

        atemp = 0.
        if(DTAUI(L,NW,NG) .GT. 1.0D-9) then
           do iaer=1,naerkind
              atemp = atemp +                                     &
                   GIAER(K,NW,IAER)   * TAUAEROLK(K,NW,IAER) +    &
                   GIAER(K+1,NW,IAER) * TAUAEROLK(K+1,NW,IAER)
           end do
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
        do iaer=1,naerkind
           atemp = atemp + GIAER(K,NW,IAER)   * TAUAEROLK(K,NW,IAER)
        end do
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

  !open(127,file='taucum.out')
  !do nw=1,L_NSPECTI
  !   write(127,*) taucumi(L_LEVELS,nw,L_NGAUSS)
  !enddo
  !close(127)
  
!  print*,'WBARI'
!  print*,WBARI
!  print*,'DTAUI'
!  print*,DTAUI
!  call abort

end subroutine optci

END MODULE optci_mod

