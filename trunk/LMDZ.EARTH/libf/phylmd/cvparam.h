!
! $Header$
!
c------------------------------------------------------------
c Parameters for convectL:
c (includes - microphysical parameters, 
c			- parameters that control the rate of approach 
c               to quasi-equilibrium)
c			- noff & minorig (previously in input of convect1)
c------------------------------------------------------------

      integer noff, minorig, nl, nlp, nlm
      real elcrit, tlcrit
      real entp
      real sigs, sigd
      real omtrain, omtsnow, coeffr, coeffs
      real dtmax
      real cu
      real betad
      real alpha, damp
      real delta

      COMMON /cvparam/ noff, minorig, nl, nlp, nlm
     :                ,elcrit, tlcrit
     :                ,entp, sigs, sigd
     :                ,omtrain, omtsnow, coeffr, coeffs
     :                ,dtmax, cu, betad, alpha, damp, delta

c$OMP THREADPRIVATE(/cvparam/)
