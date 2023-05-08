!
! $Header$
!
C
C================================================================================
C
      SUBROUTINE CLOUDS_GNO(klon,ND,R,RS,QSUB,PTCONV,RATQSC,CLDF)
      IMPLICIT NONE
C     
C--------------------------------------------------------------------------------
C
C Inputs:
C
C  ND----------: Number of vertical levels
C  R--------ND-: Domain-averaged mixing ratio of total water 
C  RS-------ND-: Mean saturation humidity mixing ratio within the gridbox
C  QSUB-----ND-: Mixing ratio of condensed water within clouds associated
C                with SUBGRID-SCALE condensation processes (here, it is
C                predicted by the convection scheme)
C Outputs:
C
C  PTCONV-----ND-: Point convectif = TRUE
C  RATQSC-----ND-: Largeur normalisee de la distribution
C  CLDF-----ND-: Fraction nuageuse
C
C--------------------------------------------------------------------------------


      INTEGER klon,ND
      REAL R(klon,ND),  RS(klon,ND), QSUB(klon,ND)
      LOGICAL PTCONV(klon,ND)
      REAL RATQSC(klon,ND)
      REAL CLDF(klon,ND)

c -- parameters controlling the iteration:
c --    nmax    : maximum nb of iterations (hopefully never reached)
c --    epsilon : accuracy of the numerical resolution 
c --    vmax    : v-value above which we use an asymptotic expression for ERF(v)

      INTEGER nmax
      PARAMETER ( nmax = 10) 
      REAL epsilon, vmax0, vmax(klon)
      PARAMETER ( epsilon = 0.02, vmax0 = 2.0 ) 

      REAL min_mu, min_Q
      PARAMETER ( min_mu =  1.e-12, min_Q=1.e-12 )
     
      INTEGER i,K, n, m
      REAL mu(klon), qsat, delta(klon), beta(klon) 
      real zu2,zv2
      REAL xx(klon), aux(klon), coeff, block
      REAL  dist, fprime, det
      REAL pi, u, v, erfcu, erfcv
      REAL  xx1, xx2
      real erf,hsqrtlog_2,v2
      real sqrtpi,sqrt2,zx1,zx2,exdel
c lconv = true si le calcul a converge (entre autre si qsub < min_q)
       LOGICAL lconv(klon)

!cdir arraycomb
      cldf  (1:klon,1:ND)=0.0        ! cym
      ratqsc(1:klon,1:ND)=0.0
      ptconv(1:klon,1:ND)=.false.
!cdir end arraycomb
      
      pi = ACOS(-1.)
      sqrtpi=sqrt(pi)
      sqrt2=sqrt(2.)
      hsqrtlog_2=0.5*SQRT(log(2.))

      DO 500 K = 1, ND

                                    do i=1,klon ! vector
      mu(i) = R(i,K)
      mu(i) = MAX(mu(i),min_mu)
      qsat = RS(i,K) 
      qsat = MAX(qsat,min_mu)
      delta(i) = log(mu(i)/qsat)
c                                   enddo ! vector

C
C ***          There is no subgrid-scale condensation;        ***
C ***   the scheme becomes equivalent to an "all-or-nothing"  *** 
C ***             large-scale condensation scheme.            ***
C

C
C ***     Some condensation is produced at the subgrid-scale       ***
C ***                                                              ***
C ***       PDF = generalized log-normal distribution (GNO)        ***
C ***   (k<0 because a lower bound is considered for the PDF)      ***
C ***                                                              ***
C ***  -> Determine x (the parameter k of the GNO PDF) such        ***
C ***  that the contribution of subgrid-scale processes to         ***
C ***  the in-cloud water content is equal to QSUB(K)              ***
C ***  (equations (13), (14), (15) + Appendix B of the paper)      ***
C ***                                                              ***
C ***    Here, an iterative method is used for this purpose        ***
C ***    (other numerical methods might be more efficient)         ***
C ***                                                              ***
C ***          NB: the "error function" is called ERF              ***
C ***                 (ERF in double precision)                   ***
C

c  On commence par eliminer les cas pour lesquels on n'a pas
c  suffisamment d'eau nuageuse.

c                                   do i=1,klon ! vector

      IF ( QSUB(i,K) .lt. min_Q ) THEN
        ptconv(i,k)=.false.
        ratqsc(i,k)=0.
        lconv(i)  = .true.

c   Rien on a deja initialise

      ELSE 

        lconv(i)  = .FALSE. 
        vmax(i) = vmax0

        beta(i) = QSUB(i,K)/mu(i) + EXP( -MIN(0.0,delta(i)) )

c --  roots of equation v > vmax:

        det = delta(i) + vmax(i)*vmax(i)
        if (det.LE.0.0) vmax(i) = vmax0 + 1.0
        det = delta(i) + vmax(i)*vmax(i)

        if (det.LE.0.) then
          xx(i) = -0.0001
        else 
         zx1=-sqrt2*vmax(i)
         zx2=SQRT(1.0+delta(i)/(vmax(i)*vmax(i)))
         xx1=zx1*(1.0-zx2)
         xx2=zx1*(1.0+zx2)
         xx(i) = 1.01 * xx1
         if ( xx1 .GE. 0.0 ) xx(i) = 0.5*xx2
        endif
        if (delta(i).LT.0.) xx(i) = -hsqrtlog_2

      ENDIF

                                    enddo       ! vector

c----------------------------------------------------------------------
c   Debut des nmax iterations pour trouver la solution.
c----------------------------------------------------------------------

      DO n = 1, nmax 

                                    do i=1,klon ! vector
        if (.not.lconv(i)) then

          u = delta(i)/(xx(i)*sqrt2) + xx(i)/(2.*sqrt2)
          v = delta(i)/(xx(i)*sqrt2) - xx(i)/(2.*sqrt2)
          v2 = v*v

          IF ( v .GT. vmax(i) ) THEN 

            IF (     ABS(u)  .GT. vmax(i) 
     :          .AND.  delta(i) .LT. 0. ) THEN

c -- use asymptotic expression of erf for u and v large:
c ( -> analytic solution for xx )
             exdel=beta(i)*EXP(delta(i))
             aux(i) = 2.0*delta(i)*(1.-exdel)
     :                       /(1.+exdel)
             if (aux(i).lt.0.) then
c                print*,'AUX(',i,',',k,')<0',aux(i),delta(i),beta(i)
                aux(i)=0.
             endif
             xx(i) = -SQRT(aux(i))
             block = EXP(-v*v) / v / sqrtpi
             dist = 0.0
             fprime = 1.0

            ELSE

c -- erfv -> 1.0, use an asymptotic expression of erfv for v large:

             erfcu = 1.0-ERF(u)
c  !!! ATTENTION : rajout d'un seuil pour l'exponentiel
             aux(i) = sqrtpi*erfcu*EXP(min(v2,100.))
             coeff = 1.0 - 0.5/(v2) + 0.75/(v2*v2)
             block = coeff * EXP(-v2) / v / sqrtpi
             dist = v * aux(i) / coeff - beta(i)
             fprime = 2.0 / xx(i) * (v2)
     :           * ( EXP(-delta(i)) - u * aux(i) / coeff )
     :           / coeff
            
            ENDIF ! ABS(u)

          ELSE

c -- general case:

           erfcu = 1.0-ERF(u)
           erfcv = 1.0-ERF(v)
           block = erfcv
           dist = erfcu / erfcv - beta(i)
           zu2=u*u
           zv2=v2
           if(zu2.gt.20..or. zv2.gt.20.) then
c              print*,'ATTENTION !!! xx(',i,') =', xx(i)
c           print*,'ATTENTION !!! klon,ND,R,RS,QSUB,PTCONV,RATQSC,CLDF',
c     .klon,ND,R(i,k),RS(i,k),QSUB(i,k),PTCONV(i,k),RATQSC(i,k),
c     .CLDF(i,k)
c              print*,'ATTENTION !!! zu2 zv2 =',zu2(i),zv2(i)
              zu2=20.
              zv2=20.
             fprime = 0.
           else
             fprime = 2. /sqrtpi /xx(i) /(erfcv*erfcv)
     :           * (   erfcv*v*EXP(-zu2)
     :               - erfcu*u*EXP(-zv2) )
           endif
          ENDIF ! x

c -- test numerical convergence:

!          if (beta(i).lt.1.e-10) then
!              print*,'avant test ',i,k,lconv(i),u(i),v(i),beta(i)
!              stop
!          endif
          if (abs(fprime).lt.1.e-11) then
!              print*,'avant test fprime<.e-11 '
!     s        ,i,k,lconv(i),u(i),v(i),beta(i),fprime(i)
!              print*,'klon,ND,R,RS,QSUB',
!     s        klon,ND,R(i,k),rs(i,k),qsub(i,k)
              fprime=sign(1.e-11,fprime)
          endif


          if ( ABS(dist/beta(i)) .LT. epsilon ) then 
c           print*,'v-u **2',(v(i)-u(i))**2
c           print*,'exp v-u **2',exp((v(i)-u(i))**2)
            ptconv(i,K) = .TRUE. 
            lconv(i)=.true.
c  borne pour l'exponentielle
            ratqsc(i,k)=min(2.*(v-u)*(v-u),20.)
            ratqsc(i,k)=sqrt(exp(ratqsc(i,k))-1.)
            CLDF(i,K) = 0.5 * block
          else
            xx(i) = xx(i) - dist/fprime
          endif
c         print*,'apres test ',i,k,lconv(i)

        endif ! lconv
                                    enddo       ! vector

c----------------------------------------------------------------------
c   Fin des nmax iterations pour trouver la solution.
        ENDDO ! n
c----------------------------------------------------------------------

500   CONTINUE  ! K

       RETURN
       END

 

