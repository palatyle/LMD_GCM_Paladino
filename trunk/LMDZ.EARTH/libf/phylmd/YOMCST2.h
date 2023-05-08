
      INTEGER choice, iflag_mix
      REAL  gammas, alphas, betas, Fmax, qqa1, qqa2, qqa3, scut
      REAL  Qcoef1max,Qcoef2max,Supcrit1,Supcrit2
!
      COMMON/YOMCST2/gammas,    alphas, betas, Fmax, scut,              &
     &               qqa1, qqa2, qqa3,                                  &
     &               Qcoef1max,Qcoef2max,                               &
     &               Supcrit1, Supcrit2,                                &
     &               choice,iflag_mix
!$OMP THREADPRIVATE(/YOMCST2/)
!    --------------------------------------------------------------------

