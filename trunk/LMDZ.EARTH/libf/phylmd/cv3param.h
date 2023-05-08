c------------------------------------------------------------
c Parameters for convectL, iflag_con=3:
c (includes - microphysical parameters,
c			- parameters that control the rate of approach
c               to quasi-equilibrium)
c			- noff & minorig (previously in input of convect1)
c------------------------------------------------------------

      integer noff, minorig, nl, nlp, nlm
      real sigdz, spfac
      real pbcrit, ptcrit
      real omtrain
      real dtovsh, dpbase, dttrig
      real dtcrit, tau, beta, alpha, alpha1
      real delta
      real betad

      COMMON /cv3param/  noff, minorig, nl, nlp, nlm
     :                ,  sigdz, spfac
     :                ,pbcrit, ptcrit
     :                ,omtrain
     :                ,dtovsh, dpbase, dttrig
     :                ,dtcrit, tau, beta, alpha, alpha1, delta, betad
!$OMP THREADPRIVATE(/cv3param/)

