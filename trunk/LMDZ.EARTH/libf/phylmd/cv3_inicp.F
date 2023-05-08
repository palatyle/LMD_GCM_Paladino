        SUBROUTINE cv3_inicp()
*
***************************************************************
*                                                             *
* CV3_INIP Lecture des choix de lois de probabilité de mélange*
*          et calcul de leurs coefficients normalisés.        *
*                                                             *
* written by   : Jean-Yves Grandpeix, 06/06/2006, 19.39.27    *
* modified by :                                               *
***************************************************************
*
#include "YOMCST2.h"
c
      INTEGER iflag_clos
      CHARACTER (LEN=20) :: modname='cv3_inicp'
      CHARACTER (LEN=80) :: abort_message
c
c --   Mixing probability distribution functions
c
      real Qcoef1,Qcoef2,QFF,QFFF,Qmix,Rmix,Qmix1,Rmix1,Qmix2,Rmix2,F
      Qcoef1(F) = tanh(F/gammas)
      Qcoef2(F) = ( tanh(F/gammas) + gammas *
     $            log(cosh((1.- F)/gammas)/cosh(F/gammas)))
      QFF(F) = Max(Min(F,1.),0.)
      QFFF(F) = Min(QFF(F),scut)
      Qmix1(F) = ( tanh((QFF(F) - Fmax)/gammas)+Qcoef1max )/
     $           Qcoef2max
      Rmix1(F) = ( gammas*log(cosh((QFF(F)-Fmax)/gammas))
     1             +QFF(F)*Qcoef1max ) / Qcoef2max
      Qmix2(F) = -Log(1.-QFFF(F))/scut
      Rmix2(F) = (QFFF(F)+(1.-QFF(F))*Log(1.-QFFF(F)))/scut
      Qmix(F) = qqa1*Qmix1(F) + qqa2*Qmix2(F)
      Rmix(F) = qqa1*Rmix1(F) + qqa2*Rmix2(F)
C
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C
C===========================================================================
C    	READ IN PARAMETERS FOR THE MIXING DISTRIBUTION
C	AND PASS THESE THROUGH A COMMON BLOCK TO SUBROUTINE CONVECT etc.
C       (Written by V.T.J. Phillips, 20-30/Jan/99)
C===========================================================================
C
C   line 1:  a flag (0 or 1) to decide whether P(F) = 1 or the general P(F) is to be
C         used, followed by SCUT, which is the cut-off value of F in CONVECT
C   line 2:  blank
C   line 3:  the coefficients for the linear combination of P(F)s to
C                 make the general P(F)
C   line 4:  blank
C   line 5:  gammas, Fmax for the cosh^2 component of P(F)
C   line 6:  blank
C   line 7:  alphas for the 1st irrational P(F)
C   line 8:  blank
C   line 9:  betas  for the 2nd irrational P(F)
C

c        open(57,file='parameter_mix.data')

c        read(57,*) iflag_clos
c        read(57,*) iflag_mix, scut
c        read(57,*)
c        if(iflag_mix .gt. 0) then
c	      read(57,*) qqa1, qqa2
c              read(57,*)
c              read(57,*) gammas, Fmax
c              read(57,*)
c              read(57,*) alphas
c         endif
c	 close(57)
c
      if(iflag_mix .gt. 0) then
c
c--      Normalize Pdf weights
c
        sumcoef=qqa1+qqa2
        qqa1=qqa1/sumcoef
        qqa2=qqa2/sumcoef
c
        Qcoef1max = Qcoef1(Fmax)
        Qcoef2max = Qcoef2(Fmax)
c
        sigma = 0.
        aire=0.0
        pdf=0.0
        mu=0.0
        df = 0.0001
c
c        do ff = 0.0 + df, 1.0 - 2.*df, df
         ff=df
         dowhile ( ff .le. 1.0 - 2.*df )
              pdf = (Qmix(ff+df) -  Qmix(ff)) * (1.-ff) / df
              aire=aire+(Qmix(ff+df) - Qmix(ff)) * (1.-ff)
              mu = mu + pdf * ff * df
cc              write(*,*) pdf,  Qmix(ff), aire, ff
         ff=ff+df
         enddo
c
c         do ff=0.0+df,1.0 - 2.*df,df
         ff=df
         dowhile ( ff .le. 1.0 - 2.*df )
              pdf = (Qmix(ff+df)- Qmix(ff)) * (1.-ff) / df
              sigma = sigma+pdf*(ff - mu)*(ff - mu)*df
         ff=ff+df
         enddo
         sigma = sqrt(sigma)
c
        if (abs(aire-1.0) .gt. 0.02) then
            print *,'WARNING:: AREA OF MIXING PDF IS::', aire
            abort_message = ''
            CALL abort_gcm (modname,abort_message,1)
        else
            print *,'Area, mean & std deviation are ::', aire,mu,sigma
        endif
      endif     !  (iflag_mix .gt. 0)

      RETURN
      END
