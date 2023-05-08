! $Id: $
      SUBROUTINE disvert_noterre

c    Auteur :  F. Forget Y. Wanherdrick, P. Levan
c    Nouvelle version 100% Mars !!
c    On l'utilise aussi pour Venus et Titan, legerment modifiee.

#ifdef CPP_IOIPSL
      use IOIPSL
#else
! if not using IOIPSL, we still need to use (a local version of) getin
      use ioipsl_getincom
#endif
      USE comvert_mod, ONLY: ap,bp,aps,bps,presnivs,pseudoalt,
     .			nivsig,nivsigs,pa,preff,scaleheight
      USE comconst_mod, ONLY: kappa
      USE logic_mod, ONLY: hybrid

      IMPLICIT NONE

#include "dimensions.h"
#include "paramet.h"
#include "iniprint.h"
c
c=======================================================================
c    Discretisation verticale en coordonnée hybride (ou sigma)
c
c=======================================================================
c
c   declarations:
c   -------------
c
c
      INTEGER l,ll
      REAL snorm
      REAL alpha,beta,gama,delta,deltaz
      real quoi,quand
      REAL zsig(llm),sig(llm+1)
      INTEGER np,ierr
      integer :: ierr1,ierr2,ierr3,ierr4
      REAL x

      REAL SSUM
      EXTERNAL SSUM
      real newsig 
      REAL dz0,dz1,nhaut,sig1,esig,csig,zz
      real tt,rr,gg, prevz
      real s(llm),dsig(llm) 

      integer iz 
      real z, ps,p
      character(len=*),parameter :: modname="disvert_noterre"

c
c-----------------------------------------------------------------------
c
! Initializations:
!      pi=2.*ASIN(1.) ! already done in iniconst
      
      hybrid=.true. ! default value for hybrid (ie: use hybrid coordinates)
      CALL getin('hybrid',hybrid)
      write(lunout,*) trim(modname),': hybrid=',hybrid

! Ouverture possible de fichiers typiquement E.T.

         open(99,file="esasig.def",status='old',form='formatted',
     s   iostat=ierr2)
         if(ierr2.ne.0) then
              close(99)
              open(99,file="z2sig.def",status='old',form='formatted',
     s        iostat=ierr4)
         endif

c-----------------------------------------------------------------------
c   cas 1 on lit les options dans esasig.def:
c   ----------------------------------------

      IF(ierr2.eq.0) then

c        Lecture de esasig.def :
c        Systeme peu souple, mais qui respecte en theorie
c        La conservation de l'energie (conversion Energie potentielle
c        <-> energie cinetique, d'apres la note de Frederic Hourdin...

         write(lunout,*)'*****************************'
         write(lunout,*)'WARNING reading esasig.def'
         write(lunout,*)'*****************************'
         READ(99,*) scaleheight
         READ(99,*) dz0
         READ(99,*) dz1
         READ(99,*) nhaut
         CLOSE(99)

         dz0=dz0/scaleheight
         dz1=dz1/scaleheight

         sig1=(1.-dz1)/tanh(.5*(llm-1)/nhaut)

         esig=1.

         do l=1,20
            esig=-log((1./sig1-1.)*exp(-dz0)/esig)/(llm-1.)
         enddo
         csig=(1./sig1-1.)/(exp(esig)-1.)

         DO L = 2, llm
            zz=csig*(exp(esig*(l-1.))-1.)
            sig(l) =1./(1.+zz)
     &      * tanh(.5*(llm+1-l)/nhaut)
         ENDDO
         sig(1)=1.
         sig(llm+1)=0.
         quoi      = 1. + 2.* kappa
         s( llm )  = 1.
         s(llm-1) = quoi
         IF( llm.gt.2 )  THEN
            DO  ll = 2, llm-1
               l         = llm+1 - ll
               quand     = sig(l+1)/ sig(l)
               s(l-1)    = quoi * (1.-quand) * s(l)  + quand * s(l+1)
            ENDDO
         END IF
c
         snorm=(1.-.5*sig(2)+kappa*(1.-sig(2)))*s(1)+.5*sig(2)*s(2)
         DO l = 1, llm
            s(l)    = s(l)/ snorm
         ENDDO

c-----------------------------------------------------------------------
c   cas 2 on lit les options dans z2sig.def:
c   ----------------------------------------

      ELSE IF(ierr4.eq.0) then
         write(lunout,*)'****************************'
         write(lunout,*)'Reading z2sig.def'
         write(lunout,*)'****************************'

         READ(99,*) scaleheight
         do l=1,llm
            read(99,*) zsig(l)
         end do
         CLOSE(99)

         sig(1) =1
         do l=2,llm
           sig(l) = 0.5 * ( exp(-zsig(l)/scaleheight) + 
     &                      exp(-zsig(l-1)/scaleheight) )
         end do
         sig(llm+1) =0

c-----------------------------------------------------------------------
      ELSE
         write(lunout,*) 'didn t you forget something ??? '
         write(lunout,*) 'We need file  z2sig.def ! (OR esasig.def)'
         stop
      ENDIF
c-----------------------------------------------------------------------

      DO l=1,llm
        nivsigs(l) = REAL(l)
      ENDDO

      DO l=1,llmp1
        nivsig(l)= REAL(l)
      ENDDO

 
c-----------------------------------------------------------------------
c    ....  Calculs  de ap(l) et de bp(l)  ....
c    .........................................
c
c   .....  pa et preff sont lus  sur les fichiers start par dynetat0 .....
c-----------------------------------------------------------------------
c

      if (hybrid) then  ! use hybrid coordinates
         write(lunout,*) "*********************************"
         write(lunout,*) "Using hybrid vertical coordinates"
         write(lunout,*) 
c        Coordonnees hybrides avec mod
         DO l = 1, llm

         call sig_hybrid(sig(l),pa,preff,newsig)
            bp(l) = EXP( 1. - 1./(newsig**2)  )
            ap(l) = pa * (newsig - bp(l) )
         enddo
         bp(llmp1) = 0.
         ap(llmp1) = 0.
      else ! use sigma coordinates
         write(lunout,*) "********************************"
         write(lunout,*) "Using sigma vertical coordinates"
         write(lunout,*) 
c        Pour ne pas passer en coordonnees hybrides
         DO l = 1, llm
            ap(l) = 0.
            bp(l) = sig(l)
         ENDDO
         ap(llmp1) = 0.
      endif

      bp(llmp1) =   0.

      write(lunout,*) trim(modname),': BP '
      write(lunout,*)  bp
      write(lunout,*) trim(modname),': AP '
      write(lunout,*)  ap

c     Calcul au milieu des couches :
c     WARNING : le choix de placer le milieu des couches au niveau de
c     pression intermédiaire est arbitraire et pourrait etre modifié.
c     Le calcul du niveau pour la derniere couche 
c     (on met la meme distance (en log pression)  entre P(llm)
c     et P(llm -1) qu'entre P(llm-1) et P(llm-2) ) est
c     Specifique.  Ce choix est spécifié ici ET dans exner_milieu.F

      DO l = 1, llm-1
       aps(l) =  0.5 *( ap(l) +ap(l+1)) 
       bps(l) =  0.5 *( bp(l) +bp(l+1)) 
      ENDDO
     
      if (hybrid) then
         aps(llm) = aps(llm-1)**2 / aps(llm-2) 
         bps(llm) = 0.5*(bp(llm) + bp(llm+1))
      else
         bps(llm) = bps(llm-1)**2 / bps(llm-2) 
         aps(llm) = 0.
      end if

      write(lunout,*) trim(modname),': BPs '
      write(lunout,*)  bps
      write(lunout,*) trim(modname),': APs'
      write(lunout,*)  aps

      DO l = 1, llm
       presnivs(l) = aps(l)+bps(l)*preff
       pseudoalt(l) = -scaleheight*log(presnivs(l)/preff)
      ENDDO

      write(lunout,*)trim(modname),' : PRESNIVS' 
      write(lunout,*)presnivs 
      write(lunout,*)'Pseudo altitude of Presnivs : (for a scale ',
     &                'height of ',scaleheight,' km)' 
      write(lunout,*)pseudoalt

c     --------------------------------------------------
c     This can be used to plot the vertical discretization
c     (> xmgrace -nxy testhybrid.tab )
c     --------------------------------------------------
c     open (53,file='testhybrid.tab')
c     scaleheight=15.5
c     do iz=0,34
c       z = -5 + min(iz,34-iz)
c     approximation of scale height for Venus
c       scaleheight = 15.5 - z/55.*10.
c       ps = preff*exp(-z/scaleheight)
c       zsig(1)= -scaleheight*log((aps(1) + bps(1)*ps)/preff)
c       do l=2,llm
c     approximation of scale height for Venus
c          if (zsig(l-1).le.55.) then
c             scaleheight = 15.5 - zsig(l-1)/55.*10.
c          else
c             scaleheight = 5.5 - (zsig(l-1)-55.)/35.*2.
c          endif
c          zsig(l)= zsig(l-1)-scaleheight*
c    .    log((aps(l) + bps(l)*ps)/(aps(l-1) + bps(l-1)*ps))
c       end do
c       write(53,'(I3,50F10.5)') iz, zsig
c      end do
c      close(53)
c     --------------------------------------------------


      RETURN
      END

c ************************************************************
      subroutine sig_hybrid(sig,pa,preff,newsig)
c     ----------------------------------------------
c     Subroutine utilisee pour calculer des valeurs de sigma modifie
c     pour conserver les coordonnees verticales decrites dans
c     esasig.def/z2sig.def lors du passage en coordonnees hybrides
c     F. Forget 2002
c     Connaissant sig (niveaux "sigma" ou on veut mettre les couches)
c     L'objectif est de calculer newsig telle que
c       (1 -pa/preff)*exp(1-1./newsig**2)+(pa/preff)*newsig = sig
c     Cela ne se résoud pas analytiquement: 
c     => on résoud par iterration bourrine 
c     ----------------------------------------------
c     Information  : where exp(1-1./x**2) become << x
c           x      exp(1-1./x**2) /x
c           1           1
c           0.68       0.5
c           0.5        1.E-1
c           0.391      1.E-2
c           0.333      1.E-3
c           0.295      1.E-4
c           0.269      1.E-5
c           0.248      1.E-6
c        => on peut utiliser newsig = sig*preff/pa si sig*preff/pa < 0.25


      implicit none
      real x1, x2, sig,pa,preff, newsig, F
      integer j

      newsig = sig
      x1=0
      x2=1
      if (sig.ge.1) then
            newsig= sig
      else if (sig*preff/pa.ge.0.25) then
        DO J=1,9999  ! nombre d''iteration max
          F=((1 -pa/preff)*exp(1-1./newsig**2)+(pa/preff)*newsig)/sig
c         write(0,*) J, ' newsig =', newsig, ' F= ', F
          if (F.gt.1) then
              X2 = newsig
              newsig=(X1+newsig)*0.5
          else
              X1 = newsig
              newsig=(X2+newsig)*0.5
          end if
c         Test : on arete lorsque on approxime sig à moins de 0.01 m près 
c         (en pseudo altitude) :
          IF(abs(10.*log(F)).LT.1.E-5) goto 999
        END DO
       else   !    if (sig*preff/pa.le.0.25) then
             newsig= sig*preff/pa
       end if
 999   continue
       Return
      END
