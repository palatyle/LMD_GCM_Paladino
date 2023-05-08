!
! $Header$
!
c------------------------------------------------------------
c Parameters for convectL, iflag_con=30:
c (includes - microphysical parameters, 
c			- parameters that control the rate of approach 
c               to quasi-equilibrium)
c			- noff & minorig (previously in input of convect1)
c------------------------------------------------------------

      integer noff, minorig, nl, nlp, nlm
      real sigd, spfac
cIM cf. FH : pour compatibilite avec conema3 TEMPORAIRE   real pbcrit, ptcrit, epmax
      real pbcrit, ptcrit
      real omtrain
      real dtovsh, dpbase, dttrig
      real dtcrit, tau, beta, alpha
      real delta
      real betad

      COMMON /cv30param/  noff, minorig, nl, nlp, nlm
     :                ,  sigd, spfac
cIM cf. FH : pour compatibilite avec conema3 TEMPORAIRE  :                ,pbcrit, ptcrit, epmax
     :                ,pbcrit, ptcrit
     :                ,omtrain
     :                ,dtovsh, dpbase, dttrig
     :                ,dtcrit, tau, beta, alpha, delta, betad

c$OMP THREADPRIVATE(/cv30param/)
