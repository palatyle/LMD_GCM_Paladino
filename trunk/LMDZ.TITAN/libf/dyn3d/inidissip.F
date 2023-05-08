      SUBROUTINE inidissip ( lstardis,nitergdiv,nitergrot,niterh  ,
     *                       tetagdiv,tetagrot,tetatemp             )
c=======================================================================
c   initialisation de la dissipation horizontale
c=======================================================================
c-----------------------------------------------------------------------
c   declarations:
c   -------------

      use control_mod, only: dissip_period, iperiod
      USE comvert_mod, ONLY: aps,bps,pseudoalt,preff
      USE comconst_mod, ONLY: dtdiss,dtvr
      IMPLICIT NONE
#include "dimensions.h"
#include "paramet.h"
#include "comdissipn.h"
!#include "control.h"

      LOGICAL lstardis
      INTEGER nitergdiv,nitergrot,niterh
      REAL    tetagdiv,tetagrot,tetatemp
      REAL zvert(llm),zz
      REAL zh(ip1jmp1),zu(ip1jmp1),zv(ip1jm),deltap(ip1jmp1,llm)
      REAL ullm,vllm,umin,vmin,zhmin,zhmax
      REAL zllm,z1llm

      INTEGER l,ij,idum,ii
      REAL tetamin

      EXTERNAL ran1
      REAL ran1
      real sig_s(llm)
      save sig_s
      logical firstcall
      data firstcall/.true./
      save firstcall 

      REAL fac_mid
      REAL fac_up
      REAL delta
      REAL middle,startalt
      SAVE fac_mid, fac_up, delta, startalt, middle

c ------------------------------------------------------
      if (firstcall) then
         firstcall=.false.
         do l=1,llm
            sig_s(l)=aps(l)/preff + bps(l)
         enddo

c        COMPUTING THE VARIATION OF DISSIPATION AS A FUNCTION OF MODEL TOP :
c         FF 2004       
         if (pseudoalt(llm).lt.160.) then 
!     currently disabled for the universal model!!!
             fac_mid=2  ! coeff  pour dissipation aux basses/moyennes altitudes
             fac_up=10 ! coeff multiplicateur pour dissipation hautes altitudes
             startalt=90. ! altitude en Km de la transition mid / up
             delta=20.! Intervalle (km) pour le changement mid / up

!             fac_mid=1  ! coeff  pour dissipation aux basses/moyennes altitudes
!             fac_up=1 ! coeff multiplicateur pour dissipation hautes altitudes
!             startalt=100. ! altitude en Km de la transition mid / up
!             delta=20.! Intervalle (km) pour le changement mid / up

         else ! thermosphere model
             fac_mid=2 ! coeff pour dissipation aux basses/moyennes altitudes
             fac_up=25 ! coeff multiplicateur pour dissipation hautes altitudes
c            startalt: 95 OK for MY24
             startalt=95. ! altitude en Km de la transition mid / up
             delta=30.! Intervalle (km) pour le changement mid /up
         end if
         middle=startalt+delta/2
         write(*,*) 'Dissipation : '
         write(*,*) 'Multiplication de la dissipation en altitude :',
     &          ' fac_mid, fac_up =', fac_mid, fac_up 
         write(*,*) 'Transition mid /up : startalt, delta =',
     &             startalt, delta , '(km)'
      endif

c-----------------------------------------------------------------------
c
c   calcul des valeurs propres des operateurs par methode iterrative:
c   -----------------------------------------------------------------

      crot     = -1.
      cdivu    = -1.
      cdivh    = -1.

c   calcul de la valeur propre de divgrad:
c   --------------------------------------
      idum = 0
      DO l = 1, llm
       DO ij = 1, ip1jmp1
        deltap(ij,l) = 1.
       ENDDO
      ENDDO

      idum  = -1
      zh(1) = RAN1(idum)-.5
      idum  = 0
      DO ij = 2, ip1jmp1
        zh(ij) = RAN1(idum) -.5
      ENDDO

      CALL filtreg (zh,jjp1,1,2,1,.TRUE.,1)

      CALL minmax(iip1*jjp1,zh,zhmin,zhmax )

      IF ( zhmin .GE. zhmax  )     THEN
         PRINT*,'  Inidissip  zh min max  ',zhmin,zhmax
         STOP'probleme generateur alleatoire dans inidissip'
      ENDIF

      zllm = ABS( zhmax )
      DO l = 1,50
         IF(lstardis) THEN
            CALL divgrad2(1,zh,deltap,niterh,zh)
         ELSE
            CALL divgrad (1,zh,niterh,zh)
         ENDIF

        CALL minmax(iip1*jjp1,zh,zhmin,zhmax )

         zllm  = ABS( zhmax )
         z1llm = 1./zllm
         DO ij = 1,ip1jmp1
            zh(ij) = zh(ij)* z1llm
         ENDDO
      ENDDO

      IF(lstardis) THEN
         cdivh = 1./ zllm
      ELSE
         cdivh = zllm ** ( -1./niterh )
      ENDIF

c   calcul des valeurs propres de gradiv (ii =1) et  nxgrarot(ii=2)
c   -----------------------------------------------------------------
      print*,'calcul des valeurs propres'

      DO  20  ii = 1, 2
c
         DO ij = 1, ip1jmp1
           zu(ij)  = RAN1(idum) -.5
         ENDDO
         CALL filtreg (zu,jjp1,1,2,1,.TRUE.,1)
         DO ij = 1, ip1jm
            zv(ij) = RAN1(idum) -.5
         ENDDO
         CALL filtreg (zv,jjm,1,2,1,.FALSE.,1)

         CALL minmax(iip1*jjp1,zu,umin,ullm )
         CALL minmax(iip1*jjm, zv,vmin,vllm )

         ullm = ABS ( ullm )
         vllm = ABS ( vllm )

         DO  5  l = 1, 50
            IF(ii.EQ.1) THEN
               IF(lstardis) THEN
                  CALL gradiv2( 1,zu,zv,nitergdiv,zu,zv )
               ELSE
                  CALL gradiv ( 1,zu,zv,nitergdiv,zu,zv )
               ENDIF
            ELSE
               IF(lstardis) THEN
                  CALL nxgraro2( 1,zu,zv,nitergrot,zu,zv )
               ELSE
                  CALL nxgrarot( 1,zu,zv,nitergrot,zu,zv )
               ENDIF
            ENDIF

            CALL minmax(iip1*jjp1,zu,umin,ullm )
            CALL minmax(iip1*jjm, zv,vmin,vllm )

            ullm = ABS  ( ullm )
            vllm = ABS  ( vllm )

            zllm  = MAX( ullm,vllm )
            z1llm = 1./ zllm
            DO ij = 1, ip1jmp1
              zu(ij) = zu(ij)* z1llm
            ENDDO
            DO ij = 1, ip1jm
               zv(ij) = zv(ij)* z1llm
            ENDDO
 5       CONTINUE

         IF ( ii.EQ.1 ) THEN
            IF(lstardis) THEN
               cdivu  = 1./zllm
            ELSE
               cdivu  = zllm **( -1./nitergdiv )
            ENDIF
         ELSE
            IF(lstardis) THEN
               crot   = 1./ zllm
            ELSE
               crot   = zllm **( -1./nitergrot )
            ENDIF
         ENDIF

 20   CONTINUE

c   petit test pour les operateurs non star:
c   ----------------------------------------

c     IF(.NOT.lstardis) THEN
c        fac_mid    = rad*24./float(jjm)
c        fac_mid    = fac_mid*fac_mid
c        PRINT*,'coef u ', fac_mid/cdivu, 1./cdivu
c        PRINT*,'coef r ', fac_mid/crot , 1./crot
c        PRINT*,'coef h ', fac_mid/cdivh, 1./cdivh
c     ENDIF

c-----------------------------------------------------------------------
c   variation verticale du coefficient de dissipation:
c   --------------------------------------------------

      DO l=1,llm
         zvert(l)=1.
      ENDDO

c
      DO l = 1, llm
         zz   = 1. -1./sig_s(l)
         zvert(l)= fac_mid -( fac_mid-1.)/( 1.+zz*zz )

c ---------------------------------------------------------------
c Utilisation de la fonction tangente hyperbolique pour augmenter
c arbitrairement la dissipation et donc la stabilite du modele a 
c partir d'une certaine altitude.

c   Le facteur multiplicatif de basse atmosphere etant deja pris 
c   en compte, il faut diviser le facteur multiplicatif de haute 
c   atmosphere par celui-ci.
c   ============================================================

         zvert(l)= zvert(l)*(1.0+((fac_up/fac_mid-1)
     &   *(1-(0.5*(1+tanh(-6./delta*    
     &   (10.*(-log(sig_s(l)))-middle)))))
     &   ))
      ENDDO
 
c -----------------------------------------------------------------------------

c

      PRINT*,'Constantes de temps de la diffusion horizontale'

      tetamin =  1.e+6


      DO l=1,llm
        tetaudiv(l)   = zvert(l)/tetagdiv
        tetaurot(l)   = zvert(l)/tetagrot
        tetah(l)      = zvert(l)/tetatemp

        IF( tetamin.GT. (1./tetaudiv(l)) ) tetamin = 1./ tetaudiv(l)
        IF( tetamin.GT. (1./tetaurot(l)) ) tetamin = 1./ tetaurot(l)
        IF( tetamin.GT. (1./   tetah(l)) ) tetamin = 1./    tetah(l)
      ENDDO

c     Calcul automatique de dissip_period
c     -----------------------------
c    :::::::::::::::::::::
c     A Commenter pour garder la valeur de run.def :
c     dissip_period = INT( tetamin/( 2.*dtvr*iperiod) ) * iperiod
c     dissip_period = MAX(iperiod,idissip)
c    :::::::::::::::::::::
      dtdiss  = dissip_period * dtvr
 
      PRINT *,' INIDI tetamin dtvr ',tetamin,dtvr,iperiod
      PRINT *,' INIDI tetamin dissip_period ',tetamin,dissip_period
      PRINT *,' INIDI dissip_period dtdiss ',dissip_period,dtdiss

      PRINT*,'pseudoZ(km)  zvert  dt(tetagdiv) dt(tetagrot) dt(divgrad)'
      DO l = 1,llm
         PRINT*,pseudoalt(l),zvert(l),dtdiss*tetaudiv(l),
     *        dtdiss*tetaurot(l),dtdiss*tetah(l)
         if ( (dtdiss*tetaudiv(l).gt.1.9).or.
     &        (dtdiss*tetaurot(l).gt.1.9).or.
     &        (dtdiss*tetah(l).gt.1.9)) then
          PRINT *,'STOP : your dissipation is too intense for the '      
          PRINT *,'dissipation timestep : unstable model !'
          PRINT *,'in run.def, you must increase tetagdiv,' 
          PRINT *,'(or tetagrot and tetatemp if they are smaller than'
          PRINT *,'tetagdiv) OR decrease dissip_period'
          PRINT *,'OR increase day_step)' 
          stop
         end if
  
      ENDDO
c
      RETURN
      END
