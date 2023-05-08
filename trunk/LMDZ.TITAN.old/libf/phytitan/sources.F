         SUBROUTINE sources(ngrid,nlay,
     $                      ptimestep,pz0,pu,pv,
     $                      pplev,pzlay,pzlev,
     $                      gaz1,gaz2,gaz3,
     $                      ptsrf,evapch4,reserv)

c=======================================================================
c     Calcul des flux aux interfaces pour les sources
c     CH4 a la surface
c     Production de C2H6 en haut du modele.
c     Production de C2H2 en haut du modele.
c     
c     NOTE :
c     Les gaz ont la tete en haut. 
c     ils ne suivent pas la meme convention que muphys :
c     (1 -> sol  / klev = haut du modele)
c=======================================================================

c-----------------------------------------------------------------------
c   declarations:
c   -------------

         use dimphy
         IMPLICIT NONE
#include "YOMCST.h"
c
c   arguments:
c   ----------

         INTEGER ngrid,nlay,nq,ihor
         REAL ptimestep
         REAL pplev(ngrid,nlay+1)
         REAL pzlay(ngrid,nlay),pzlev(ngrid,nlay+1)
         REAL pu(ngrid),pv(ngrid)
         REAL gaz1(ngrid,nlay),gaz2(ngrid,nlay),gaz3(ngrid,nlay)
         REAL ptsrf(ngrid)
         REAL evapch4(ngrid)
c
c   local:
c   ------
 
         INTEGER ilev,ig,ilay,nlev,k,inch4,inc2h6

         REAL zgz1(klon,klev),zgz2(klon,klev),zgz3(klon,klev)
         REAL zcdv(klon),zu2,pz0
         REAL xmair,gg,zrho,ws,ch,qch4,flux
         REAL effg               ! effg est une fonction(z), z en m.
         REAL xmuair
         REAL zmem,zmem2,zmem3
         REAL prodc2h6,prodc2h2
         real reserv(ngrid),restemp,drestemp
         REAL zevapch4

         real umin
         data umin/1.e-12/
         save umin
c
c
c-----------------------------------------------------------------------
c   initialisations:
c   -----------------

         nlev=nlay+1

         if(nlay.ne.klev) THEN
           PRINT*,'STOP dans sources.F'
           PRINT*,'probleme de dimensions :'
           PRINT*,'nlay  =',nlay
           PRINT*,'klev  =',klev
           STOP
         endif

         IF(ngrid.NE.klon) THEN
           PRINT*,'STOP dans sources.F'
           PRINT*,'probleme de dimensions :'
           PRINT*,'ngrid  =',ngrid
           PRINT*,'klon  =',klon
           STOP
         ENDIF

         zgz1 = gaz1
         zgz2 = gaz2
         zgz3 = gaz3

         evapch4 = 0.

c-----------------------------------------------------------------------
c     2. calcul de  cd :
c     ----------------
c
         DO ig=1,ngrid
           zu2=pu(ig)*pu(ig)+pv(ig)*pv(ig)+umin
           zcdv(ig)=pz0*(sqrt(zu2))
c           write(99,'(I4,3(ES24.17,1X))') ig,
c     &     pz0,zu2,(sqrt(zu2))
         ENDDO
c          write(99,*) ""

c-----------------------------------------------------------------------
c    4. Conditions aux limites pour CH4 et C2H6
c    -------------------------------------------
c

*   Conditions CH4
         DO ig=1,ngrid
           zevapch4=0.
           restemp=0.
           gg=effg(pzlay(ig,1))
           zrho=(pplev(ig,1)-pplev(ig,2))/gg
           zrho=zrho/(pzlev(ig,2)-pzlev(ig,1))
           ws=sqrt(pu(ig)**2.+pv(ig)**2.)*(10./pzlay(ig,1))**0.2
           ch=1.5*sqrt(zcdv(ig))
           call ch4sat(ptsrf(ig),pplev(ig,1),qch4) ! qch4=kg/kg
           qch4=qch4*0.50 ! ici on impose 50% d'humidité au sol 
        
           if(reserv(ig).le. 1.e-10 ) then  
             flux=0.
             reserv(ig)=1.e-10
           else
             flux=zrho*ch*ws
             flux=flux*0.1 ! fraction occupée par les lacs 
           endif

           zmem=zgz1(ig,1)
           zgz1(ig,1)=(zgz1(ig,1)+flux*ptimestep*qch4*28./16.)
     &                /(1.+flux*ptimestep)

           gg=effg(pzlay(ig,1))
           xmair=(pplev(ig,1)-pplev(ig,1+1))/gg
           xmair=xmair/(pzlev(ig,1+1)-pzlev(ig,1))
           xmuair=28.!*(1.-zmem)+zmem*16.

           drestemp = - (zgz1(ig,1)-zmem)*xmair ! en m3/m2=m
     &     *(pzlev(ig,2)-pzlev(ig,1))*16./xmuair/425.

c           ici on peut fixer un seuil sur drestemp 
c           (ie limiter l'echange atm/surface)
      
           restemp=reserv(ig) +drestemp

           IF (restemp.ge.0.) THEN 
             reserv(ig) = reserv(ig) + drestemp
             zevapch4   = zevapch4   + drestemp
           ELSE
*          Il n'y a pas suffisamment de méthane; on re-évalue le flux d'évaporation
*          Quelle nouvelle concentration zgz1(ig,1) atteint-on en évaporant tout ?
             zgz1(ig,1)= reserv(ig)/(xmair*(pzlev(ig,1+1)-pzlev(ig,1))
     &                 *16./xmuair/425.)+zmem
             zevapch4  = zevapch4-reserv(ig)

             if(reserv(ig).eq.0.) 
     &       print*,'assechement du sol en ig=', ig,reserv(ig),flux

             reserv(ig)=0.  ! on a tout évaporé...
           ENDIF
c         
           evapch4(ig) = zevapch4  ! < 0 si volume évaporé (m3/m2)

         ENDDO

*   Conditions C2H6

         prodc2h6=6.e-12/5. ! kg/m2/s
 
         IF (1.EQ.1) THEN
           DO ig=1,ngrid
             DO ilev=nlay,nlay-4,-1
*            calcule de zrho (kg/m3) pour la couche...
               gg=effg(pzlay(ig,ilev))
               zrho=(pplev(ig,ilev)-pplev(ig,ilev+1))/gg
               zrho=zrho/(pzlev(ig,ilev+1)-pzlev(ig,ilev))

*              passage taux de production --> Dx_c2h6 a rajouter au niveau ilev
               zmem2=zgz2(ig,ilev)
               zgz2(ig,ilev)=zgz2(ig,ilev)+
     &         prodc2h6*ptimestep/
     &         abs(pzlev(ig,ilev+1)-pzlev(ig,ilev))    !  kg/m3/s
     &         /zrho*28./30.

             ENDDO
           ENDDO

         ELSE

           DO ig=1,ngrid
             DO ilev=nlay,nlay-8,-1
               zgz2(ig,ilev)=1.2e-5
             ENDDO
           ENDDO

         ENDIF  ! (fin 1.eq.0)

*-------------------------------------
*   Conditions C2H2

         prodc2h2=1.6e-12/5. ! kg/m2/s

         IF(1  .EQ.  1) THEN 

           DO ig=1,ngrid
             DO ilev=nlay,nlay-4,-1
*            calcule de zrho (kg/m3) pour la couche...
                gg=effg(pzlay(ig,ilev))
                zrho=(pplev(ig,ilev)-pplev(ig,ilev+1))/gg
                zrho=zrho/(pzlev(ig,ilev+1)-pzlev(ig,ilev))

*            passage taux de production --> Dx_c2h2 a rajouter au niveau ilev
                zmem3=zgz3(ig,ilev)
                zgz3(ig,ilev)=zgz3(ig,ilev)+
     &          (prodc2h2)*ptimestep/
     &          abs(pzlev(ig,ilev+1)-pzlev(ig,ilev))  !  kg/m3/s
     &          /zrho*28./26.

             ENDDO

           ENDDO

         ENDIF

c-----------------------------------------------------------------------
         DO ig=1,ngrid
           DO ilev=1,nlay
             gaz1(ig,ilev)=zgz1(ig,ilev)
             gaz2(ig,ilev)=zgz2(ig,ilev)
             gaz3(ig,ilev)=zgz3(ig,ilev)
           ENDDO
         ENDDO

        RETURN
        END
