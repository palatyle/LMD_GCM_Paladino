      SUBROUTINE calfis(nq, lafin, rdayvrai,rday_ecri, heure,
     $            pucov,pvcov,pteta,pq,pmasse,pps,pp,ppk,pphis,pphi,
     $            pducov,pdvcov,pdteta,pdq,pw,
     $            pdufi,pdvfi,pdhfi,pdqfi,pdpsfi )
c
c    Auteur :  P. Le Van, F. Hourdin 
c   .........
      USE infotrac, ONLY: tname, nqtot
      USE comvert_mod, ONLY: preff
      USE comconst_mod, ONLY: dtphys,kappa,cpp,pi
      USE physiq_mod, ONLY: physiq

      IMPLICIT NONE
c=======================================================================
c
c   1. rearrangement des tableaux et transformation
c      variables dynamiques  >  variables physiques
c   2. calcul des termes physiques
c   3. retransformation des tendances physiques en tendances dynamiques
c
c   remarques:
c   ----------
c
c    - les vents sont donnes dans la physique par leurs composantes 
c      naturelles.
c    - la variable thermodynamique de la physique est une variable
c      intensive :   T 
c      pour la dynamique on prend    T * ( preff / p(l) ) **kappa
c    - les deux seules variables dependant de la geometrie necessaires
c      pour la physique sont la latitude pour le rayonnement et 
c      l'aire de la maille quand on veut integrer une grandeur 
c      horizontalement.
c    - les points de la physique sont les points scalaires de la 
c      la dynamique; numerotation:
c          1 pour le pole nord
c          (jjm-1)*iim pour l'interieur du domaine
c          ngridmx pour le pole sud
c      ---> ngridmx=2+(jjm-1)*iim
c
c     Input :
c     -------
c       ecritphy        frequence d'ecriture (en jours)de histphy
c       pucov           covariant zonal velocity
c       pvcov           covariant meridional velocity 
c       pteta           potential temperature
c       pps             surface pressure
c       pmasse          masse d'air dans chaque maille
c       pts             surface temperature  (K)
c       pw              flux vertical (kg/s)
c
c    Output :
c    --------
c        pdufi          tendency for the natural zonal velocity (ms-1)
c        pdvfi          tendency for the natural meridional velocity 
c        pdhfi          tendency for the potential temperature
c        pdtsfi         tendency for the surface temperature
c
c=======================================================================
c
c-----------------------------------------------------------------------
c
c    0.  Declarations :
c    ------------------

#include "dimensions.h"
#include "paramet.h"

      INTEGER ngridmx,nq
      PARAMETER( ngridmx = 2+(jjm-1)*iim - 1/jjm   )

#include "comgeom2.h"
!#include "control.h"

!#include "advtrac.h"
!! this is to get tnom (tracers name)

c    Arguments :
c    -----------
      LOGICAL  lafin
      REAL heure

      REAL pvcov(iip1,jjm,llm)
      REAL pucov(iip1,jjp1,llm)
      REAL pteta(iip1,jjp1,llm)
      REAL pmasse(iip1,jjp1,llm)
      REAL pq(iip1,jjp1,llm,nqtot)
      REAL pphis(iip1,jjp1)
      REAL pphi(iip1,jjp1,llm)
c
      REAL pdvcov(iip1,jjm,llm)
      REAL pducov(iip1,jjp1,llm)
      REAL pdteta(iip1,jjp1,llm)
      REAL pdq(iip1,jjp1,llm,nqtot)
c
      REAL pw(iip1,jjp1,llm)
c
      REAL pps(iip1,jjp1)
      REAL pp(iip1,jjp1,llmp1)
      REAL ppk(iip1,jjp1,llm)
c
      REAL pdvfi(iip1,jjm,llm)
      REAL pdufi(iip1,jjp1,llm)
      REAL pdhfi(iip1,jjp1,llm)
      REAL pdqfi(iip1,jjp1,llm,nqtot)
      REAL pdpsfi(iip1,jjp1)

c    Local variables :
c    -----------------

      INTEGER i,j,l,ig0,ig,iq
      REAL zpsrf(ngridmx)
      REAL zplev(ngridmx,llm+1),zplay(ngridmx,llm)
      REAL zphi(ngridmx,llm),zphis(ngridmx)
c
      REAL zufi(ngridmx,llm), zvfi(ngridmx,llm)
      REAL ztfi(ngridmx,llm),zqfi(ngridmx,llm,nqtot)
c
!      REAL zvervel(ngridmx,llm)
      REAL flxwfi(ngridmx,llm) ! vertical mass flux (kg/s) on physics grid
c
      REAL zdufi(ngridmx,llm),zdvfi(ngridmx,llm)
      REAL zdtfi(ngridmx,llm),zdqfi(ngridmx,llm,nqtot)
      REAL zdpsrf(ngridmx)
c
      REAL zsin(iim),zcos(iim),z1(iim)
      REAL zsinbis(iim),zcosbis(iim),z1bis(iim)
      REAL unskap, pksurcp
c
      
      EXTERNAL gr_dyn_fi,gr_fi_dyn
      REAL SSUM
      EXTERNAL SSUM

      REAL latfi(ngridmx),lonfi(ngridmx)
      REAL airefi(ngridmx)
      SAVE latfi, lonfi, airefi

      LOGICAL firstcal, debut
      DATA firstcal/.true./
      SAVE firstcal,debut
      REAL rdayvrai,rday_ecri
c
c-----------------------------------------------------------------------
c
c    1. Initialisations :
c    --------------------
c


      IF (ngridmx.NE.2+(jjm-1)*iim) THEN
         PRINT*,'STOP dans calfis'
         PRINT*,'La dimension ngridmx doit etre egale a 2 + (jjm-1)*iim'
         PRINT*,'  ngridmx  jjm   iim   '
         PRINT*,ngridmx,jjm,iim
         STOP
      ENDIF

c-----------------------------------------------------------------------
c   latitude, longitude et aires des mailles pour la physique:
c   ----------------------------------------------------------

c
      IF ( firstcal )  THEN
          debut = .TRUE.
      ELSE
          debut = .FALSE.
      ENDIF

c
!      IF (firstcal) THEN
!         latfi(1)=rlatu(1)
!         lonfi(1)=0.
!         DO j=2,jjm
!            DO i=1,iim
!               latfi((j-2)*iim+1+i)= rlatu(j)
!               lonfi((j-2)*iim+1+i)= rlonv(i)
!            ENDDO
!         ENDDO
!         latfi(ngridmx)= rlatu(jjp1)
!         lonfi(ngridmx)= 0.
!
!         ! build airefi(), mesh area on physics grid
!         CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,aire,airefi)
!         ! Poles are single points on physics grid
!         airefi(1)=airefi(1)*iim
!         airefi(ngridmx)=airefi(ngridmx)*iim
!
!         CALL inifis(ngridmx,llm,day_ini,daysec,dtphys,
!     .                latfi,lonfi,airefi,rad,g,r,cpp)
!      ENDIF

c
c-----------------------------------------------------------------------
c   40. transformation des variables dynamiques en variables physiques:
c   ---------------------------------------------------------------

c   41. pressions au sol (en Pascals)
c   ----------------------------------

       
      zpsrf(1) = pps(1,1)

      ig0  = 2
      DO j = 2,jjm
         CALL SCOPY( iim,pps(1,j),1,zpsrf(ig0), 1 )
         ig0 = ig0+iim
      ENDDO

      zpsrf(ngridmx) = pps(1,jjp1)


c   42. pression intercouches :
c
c   -----------------------------------------------------------------
c     .... zplev  definis aux (llm +1) interfaces des couches  ....
c     .... zplay  definis aux (  llm )    milieux des couches  .... 
c   -----------------------------------------------------------------

c    ...    Exner = cp * ( p(l) / preff ) ** kappa     ....
c
       unskap   = 1./ kappa
c
      DO l = 1, llmp1
        zplev( 1,l ) = pp(1,1,l)
        ig0 = 2
          DO j = 2, jjm
             DO i =1, iim
              zplev( ig0,l ) = pp(i,j,l)
              ig0 = ig0 +1
             ENDDO
          ENDDO
        zplev( ngridmx,l ) = pp(1,jjp1,l)
      ENDDO
c
c

c   43. temperature naturelle (en K) et pressions milieux couches .
c   ---------------------------------------------------------------

      DO l=1,llm

         pksurcp     =  ppk(1,1,l) / cpp
         zplay(1,l)  =  preff * pksurcp ** unskap
         ztfi(1,l)   =  pteta(1,1,l) *  pksurcp
         ig0         =  2

         DO j = 2, jjm
            DO i = 1, iim
              pksurcp        = ppk(i,j,l) / cpp
              zplay(ig0,l)   = preff * pksurcp ** unskap
              ztfi(ig0,l)    = pteta(i,j,l)  * pksurcp
              ig0            = ig0 + 1
            ENDDO
         ENDDO

         pksurcp       = ppk(1,jjp1,l) / cpp
         zplay(ig0,l)  = preff * pksurcp ** unskap
         ztfi (ig0,l)  = pteta(1,jjp1,l)  * pksurcp

      ENDDO

      DO l=1, llm
        DO ig=1,ngridmx
             if (ztfi(ig,l).lt.15) then
                  write(*,*) 'New Temperature below 15 K !!! '
	              write(*,*) 'Stop in calfis.F '
	              write(*,*) 'ig=', ig, ' l=', l
                      write(*,*) 'ztfi(ig,l)=',ztfi(ig,l)
                  stop
             end if
        ENDDO
      ENDDO



c   43.bis Taceurs (en kg/kg)
c   --------------------------
      DO iq=1,nqtot
         DO l=1,llm
            zqfi(1,l,iq) = pq(1,1,l,iq)
            ig0          = 2
            DO j=2,jjm
               DO i = 1, iim
                  zqfi(ig0,l,iq)  = pq(i,j,l,iq)
                  ig0             = ig0 + 1
               ENDDO
            ENDDO
            zqfi(ig0,l,iq) = pq(1,jjp1,l,iq)
         ENDDO
      ENDDO

c   Geopotentiel calcule par rapport a la surface locale:
c   -----------------------------------------------------

      CALL gr_dyn_fi(llm,iip1,jjp1,ngridmx,pphi,zphi)
      CALL gr_dyn_fi(1,iip1,jjp1,ngridmx,pphis,zphis)
      DO l=1,llm
         DO ig=1,ngridmx
            zphi(ig,l)=zphi(ig,l)-zphis(ig)
         ENDDO
      ENDDO

c   Calcul de la vitesse  verticale (m/s) pour diagnostique 
c   -------------------------------
c     pw est en kg/s
c On interpole "lineairement" la temperature entre les couches(FF,10/95)

!      DO ig=1,ngridmx
!         zvervel(ig,1)=0.
!      END DO
!      DO l=2,llm
!        zvervel(1,l)=(pw(1,1,l)/apoln)
!     &  * r *0.5*(ztfi(1,l)+ztfi(1,l-1)) /zplev(1,l)              
!        ig0=2
!       DO j=2,jjm
!           DO i = 1, iim
!              zvervel(ig0,l) = pw(i,j,l) * unsaire(i,j)
!     &        * r *0.5*(ztfi(ig0,l)+ztfi(ig0,l-1)) /zplev(ig0,l)              
!              ig0 = ig0 + 1
!           ENDDO
!       ENDDO
!        zvervel(ig0,l)=(pw(1,jjp1,l)/apols)
!     &  * r *0.5*(ztfi(ig0,l)+ztfi(ig0,l-1)) /zplev(ig0,l)              
!      ENDDO

c    .........  Reindexation : calcul de zvervel au MILIEU des couches
!       DO l=1,llm-1
!	      DO ig=1,ngridmx
!		     zvervel(ig,l) = 0.5*(zvervel(ig,l)+zvervel(ig,l+1))
!          END DO 
!       END DO 
c      (dans la couche llm, on garde la valeur à la limite inférieure llm)

! vertical mass flux
      ! tranfer values from dynamics grid to physics grid:
      CALL gr_dyn_fi(llm,iip1,jjp1,ngridmx,pw,flxwfi)
      ! but mass flux is an extensive variable, so take the sum at the poles
      DO l=1,llm
        flxwfi(1,l)=sum(pw(1:iim,1,l))
        flxwfi(ngridmx,l)=sum(pw(1:iim,jjp1,l))
      ENDDO

c   45. champ u:
c   ------------

      DO l=1,llm
         DO j=2,jjm
            ig0 = 1+(j-2)*iim
            zufi(ig0+1,l)= 0.5 * 
     $      ( pucov(iim,j,l)/cu(iim,j) + pucov(1,j,l)/cu(1,j) )
            DO i=2,iim
               zufi(ig0+i,l)= 0.5 *
     $         ( pucov(i-1,j,l)/cu(i-1,j) + pucov(i,j,l)/cu(i,j) )
            ENDDO
        ENDDO
      ENDDO


c   46.champ v:
c   -----------

      DO l=1,llm
         DO j=2,jjm
            ig0=1+(j-2)*iim
            DO i=1,iim
               zvfi(ig0+i,l)= 0.5 *
     $         ( pvcov(i,j-1,l)/cv(i,j-1) + pvcov(i,j,l)/cv(i,j) )
            ENDDO
         ENDDO
      ENDDO


c   47. champs de vents aux pole nord   
c   ------------------------------
c        U = 1 / pi  *  integrale [ v * cos(long) * d long ]
c        V = 1 / pi  *  integrale [ v * sin(long) * d long ]

      DO l=1,llm

         z1(1)   =(rlonu(1)-rlonu(iim)+2.*pi)*pvcov(1,1,l)/cv(1,1)
         z1bis(1)=(rlonu(1)-rlonu(iim)+2.*pi)*pdvcov(1,1,l)/cv(1,1)
         DO i=2,iim
            z1(i)   =(rlonu(i)-rlonu(i-1))*pvcov(i,1,l)/cv(i,1)
            z1bis(i)=(rlonu(i)-rlonu(i-1))*pdvcov(i,1,l)/cv(i,1)
         ENDDO

         DO i=1,iim
            zcos(i)   = COS(rlonv(i))*z1(i)
            zcosbis(i)= COS(rlonv(i))*z1bis(i)
            zsin(i)   = SIN(rlonv(i))*z1(i)
            zsinbis(i)= SIN(rlonv(i))*z1bis(i)
         ENDDO

         zufi(1,l)  = SSUM(iim,zcos,1)/pi
         zvfi(1,l)  = SSUM(iim,zsin,1)/pi

      ENDDO


c   48. champs de vents aux pole sud:
c   ---------------------------------
c        U = 1 / pi  *  integrale [ v * cos(long) * d long ]
c        V = 1 / pi  *  integrale [ v * sin(long) * d long ]

      DO l=1,llm

         z1(1)   =(rlonu(1)-rlonu(iim)+2.*pi)*pvcov(1,jjm,l)/cv(1,jjm)
         z1bis(1)=(rlonu(1)-rlonu(iim)+2.*pi)*pdvcov(1,jjm,l)/cv(1,jjm)
         DO i=2,iim
            z1(i)   =(rlonu(i)-rlonu(i-1))*pvcov(i,jjm,l)/cv(i,jjm)
            z1bis(i)=(rlonu(i)-rlonu(i-1))*pdvcov(i,jjm,l)/cv(i,jjm)
      ENDDO

         DO i=1,iim
            zcos(i)    = COS(rlonv(i))*z1(i)
            zcosbis(i) = COS(rlonv(i))*z1bis(i)
            zsin(i)    = SIN(rlonv(i))*z1(i)
            zsinbis(i) = SIN(rlonv(i))*z1bis(i)
      ENDDO

         zufi(ngridmx,l)  = SSUM(iim,zcos,1)/pi
         zvfi(ngridmx,l)  = SSUM(iim,zsin,1)/pi

      ENDDO

c-----------------------------------------------------------------------
c   Appel de la physique:
c   ---------------------


      CALL physiq (ngridmx,llm,nq,
     .     tname,
     ,     debut,lafin,
     ,     rday_ecri,heure,dtphys,
     ,     zplev,zplay,zphi,
     ,     zufi, zvfi,ztfi, zqfi,  
!     ,     zvervel,
     ,     flxwfi,
C - sorties
     s     zdufi, zdvfi, zdtfi, zdqfi,zdpsrf)


c-----------------------------------------------------------------------
c   transformation des tendances physiques en tendances dynamiques:
c   ---------------------------------------------------------------

c  tendance sur la pression :
c  -----------------------------------

      CALL gr_fi_dyn(1,ngridmx,iip1,jjp1,zdpsrf,pdpsfi)
c
ccc     CALL multipl(ip1jmp1,aire,pdpsfi,pdpsfi)

c   62. enthalpie potentielle
c   ---------------------

      DO l=1,llm

         DO i=1,iip1
          pdhfi(i,1,l)    = cpp *  zdtfi(1,l)      / ppk(i, 1  ,l)
          pdhfi(i,jjp1,l) = cpp *  zdtfi(ngridmx,l)/ ppk(i,jjp1,l)
         ENDDO

         DO j=2,jjm
            ig0=1+(j-2)*iim
            DO i=1,iim
               pdhfi(i,j,l) = cpp * zdtfi(ig0+i,l) / ppk(i,j,l)
            ENDDO
               pdhfi(iip1,j,l) =  pdhfi(1,j,l)
         ENDDO

      ENDDO


c   62. traceurs
c   ---------------------

      DO iq=1,nqtot
         DO l=1,llm
            DO i=1,iip1
               pdqfi(i,1,l,iq)    = zdqfi(1,l,iq)
               pdqfi(i,jjp1,l,iq) = zdqfi(ngridmx,l,iq)
            ENDDO
            DO j=2,jjm
               ig0=1+(j-2)*iim
               DO i=1,iim
                  pdqfi(i,j,l,iq) = zdqfi(ig0+i,l,iq)
               ENDDO
               pdqfi(iip1,j,l,iq) = pdqfi(1,j,l,iq)
            ENDDO
         ENDDO
      ENDDO

c   65. champ u:
c   ------------

      DO l=1,llm

         DO i=1,iip1
            pdufi(i,1,l)    = 0.
            pdufi(i,jjp1,l) = 0.
         ENDDO

         DO j=2,jjm
            ig0=1+(j-2)*iim
            DO i=1,iim-1
               pdufi(i,j,l)=
     $         0.5*(zdufi(ig0+i,l)+zdufi(ig0+i+1,l))*cu(i,j)
            ENDDO
            pdufi(iim,j,l)=
     $      0.5*(zdufi(ig0+1,l)+zdufi(ig0+iim,l))*cu(iim,j)
            pdufi(iip1,j,l)=pdufi(1,j,l)
         ENDDO

      ENDDO


c   67. champ v:
c   ------------

      DO l=1,llm

         DO j=2,jjm-1
            ig0=1+(j-2)*iim
            DO i=1,iim
               pdvfi(i,j,l)=
     $         0.5*(zdvfi(ig0+i,l)+zdvfi(ig0+i+iim,l))*cv(i,j)
            ENDDO
            pdvfi(iip1,j,l) = pdvfi(1,j,l)
         ENDDO
      ENDDO


c   68. champ v pres des poles:
c   ---------------------------
c      v = U * cos(long) + V * SIN(long)

      DO l=1,llm

         DO i=1,iim
            pdvfi(i,1,l)=
     $      zdufi(1,l)*COS(rlonv(i))+zdvfi(1,l)*SIN(rlonv(i))
            pdvfi(i,jjm,l)=zdufi(ngridmx,l)*COS(rlonv(i))
     $      +zdvfi(ngridmx,l)*SIN(rlonv(i))
            pdvfi(i,1,l)=
     $      0.5*(pdvfi(i,1,l)+zdvfi(i+1,l))*cv(i,1)
            pdvfi(i,jjm,l)=
     $      0.5*(pdvfi(i,jjm,l)+zdvfi(ngridmx-iip1+i,l))*cv(i,jjm)
          ENDDO

         pdvfi(iip1,1,l)  = pdvfi(1,1,l)
         pdvfi(iip1,jjm,l)= pdvfi(1,jjm,l)

      ENDDO

c-----------------------------------------------------------------------

700   CONTINUE

      firstcal = .FALSE.

      RETURN
      END
