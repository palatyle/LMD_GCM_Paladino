      SUBROUTINE grid_noro1(imdep, jmdep, xdata, ydata, entree,
     .                 imar, jmar, x, y, zmea,zstd,zsig,zgam,zthe)
c=======================================================================
c (F. Lott) (voir aussi z.x. Li, A. Harzallah et L. Fairhead)
c
c      Calcul des parametres de l'orographie sous-maille necessaires
c      au nouveau shema de representation des montagnes meso-echelles
c      dans le modele.  Les points sont mis sur une grille rectangulaire
c      pseudo-physique.  Typiquement, il y a iim+1 latitudes incluant
c      le pole nord et le pole sud.  Il y a jjm+1 longitudes, y compris
c      aux poles.  Aux poles les champs peuvent ont une valeurs repetee
c      jjm+1 fois.....  La valeur du champs en jjm+1 (jmar) est celle
c      en j=1.  
c      Les parametres a,b,c,d representent les limites de la region
c      de point de grille correspondant a un point decrit precedemment.
c      Les moyennes sur ces regions des valeurs calculees a partir de
c      l'USN, sont ponderees par un poids, fonction de la surface
c      occuppe par ces donnees a l'interieure de la grille du modele.
c      Dans la plupart des cas ce poid est le rapport entre la surface
c      de la region de point de grille USN et la surface de la region
c      de point de grille du modele.
c       
c
c           (c)
c        ----d-----
c        | . . . .|
c        |        |
c     (b)a . * . .b(a)
c        |        |
c        | . . . .|
c        ----c-----
c           (d)
C=======================================================================
c INPUT:
c        imdep, jmdep: dimensions X et Y pour depart
c        xdata, ydata: coordonnees X et Y pour depart
c        entree: champ d'entree a transformer
c        dans ce programme, on assume que les donnees sont les altitudes
c        de l'USNavy: imdep=iusn=2160, jmdep=jusn=1080.
c OUTPUT:
c        imar, jmar: dimensions X et Y d'arrivee
c        x, y: coordonnees X et Y d'arrivee
c        les champs de sorties sont sur une grille physique:
c             zmea:  orographie moyenne
c             zstd:  deviation standard de l'orographie sous-maille
c             zsig:  pente de l'orographie sous-maille 
c             zgam:  anisotropy de l'orographie sous maille
c             zthe:  orientation de l'axe oriente dans la direction
c                    de plus grande pente de l'orographie sous maille
C=======================================================================
c     IMPLICIT INTEGER (I,J)
c     IMPLICIT REAL(X,Z) 

       USE comconst_mod, ONLY: rad

       implicit none
       integer iusn,jusn,iext
       parameter(iusn=360,jusn=180,iext=40)
c!-*-      include 'param1'
c!-*-      include 'comcstfi.h'
#include "dimensions.h"
c!-*-
c!-*-      parameter(iim=cols,jjm=rows)
      REAL xusn(iusn+2*iext),yusn(jusn+2)	
      REAL zusn(iusn+2*iext,jusn+2),zusnfi(iusn+2*iext,jusn+2)

c   modif declarations pour implicit none
      real zmeanor,zmeasud,zstdnor,zstdsud,zsignor
      real zsigsud,zweinor,zweisud
      real xk,xl,xm,xw,xp,xq
      real zmaxmea,zmaxstd,zmaxsig,zmaxgam,zmaxthe,zminthe
      real zbordnor,zbordsud,zbordest,zbordoue,xpi
      real zdeltax,zdeltay,zlenx,zleny,weighx,weighy,xincr
      integer i,j,ii,jj,ideltax,ihalph

      INTEGER imdep, jmdep
      REAL xdata(imdep),ydata(jmdep) 
      REAL entree(imdep,jmdep)
c
      INTEGER imar, jmar
  
      REAL ztz(iim+1,jjm+1),zxtzx(iim+1,jjm+1)
      REAL zytzy(iim+1,jjm+1),zxtzy(iim+1,jjm+1)
      REAL zxtzxusn(iusn+2*iext,jusn+2),zytzyusn(iusn+2*iext,jusn+2)
      REAL zxtzyusn(iusn+2*iext,jusn+2)
      REAL weight(iim+1,jjm+1)
      REAL x(imar+1),y(jmar)
      REAL zmea(imar+1,jmar),zstd(imar+1,jmar)
      REAL zsig(imar+1,jmar),zgam(imar+1,jmar),zthe(imar+1,jmar)
c
      REAL a(2200),b(2200),c(1100),d(1100)
c
c  quelques constantes:
c
      print *,' parametres de l orographie a l echelle sous maille' 
      print*,'rad =',rad
      print*,'Long et lat entree'
      print*,(x(i),i=1,imar+1)
      print*,(y(j),j=1,jmar)
       print*,'Long et lat donnees'
      print*,(xdata(i),i=1,imdep)
      print*,(ydata(j),j=1,jmdep)

      xpi=acos(-1.)
      zdeltay=2.*xpi/real(jusn)*rad
c
c  quelques tests de dimensions:
c    
      IF (imar.GT.2200 .OR. jmar.GT.1100) THEN
         PRINT*, 'imar ou jmar trop grand', imar, jmar
         CALL ABORT
      ENDIF

      IF(imdep.ne.iusn.or.jmdep.ne.jusn)then
         print *,' imdep ou jmdep mal dimensionnes:',imdep,jmdep
         call abort
      ENDIF

      IF(imar+1.gt.iim+1.or.jmar.gt.jjm+1)THEN
        print *,' imar ou jmar mal dimensionnes:',imar,jmar
        call abort
      ENDIF
c
C  Extension de la base de donnee de l'USN pour faciliter
C  les calculs ulterieurs:
c
      DO j=1,jusn
        yusn(j+1)=ydata(j)
      DO i=1,iusn
        zusn(i+iext,j+1)=entree(i,j)
        xusn(i+iext)=xdata(i)
      ENDDO
      DO i=1,iext
        zusn(i,j+1)=entree(iusn-iext+i,j)
        xusn(i)=xdata(iusn-iext+i)-2.*xpi
        zusn(iusn+iext+i,j+1)=entree(i,j)
        xusn(iusn+iext+i)=xdata(i)+2.*xpi
      ENDDO
      ENDDO

        yusn(1)=ydata(1)+(ydata(1)-ydata(2))
        yusn(jusn+2)=ydata(jusn)+(ydata(jusn)-ydata(jusn-1))
       DO i=1,iusn/2+iext
        zusn(i,1)=zusn(i+iusn/2,2)
        zusn(i+iusn/2+iext,1)=zusn(i,2)
        zusn(i,jusn+2)=zusn(i+iusn/2,jusn+1)
        zusn(i+iusn/2+iext,jusn+2)=zusn(i,jusn+1)
       ENDDO
c
c  Calcul d'une orographie filtree aux hautes latitudes
c  pour permettre des calculs plus isotropiques sur la pente
c  des montagnes
c
       DO i=1,IUSN+2*iext
       DO J=1,JUSN+2
          zusnfi(i,j)=0.0
       ENDDO
       ENDDO

      DO j=1,jusn
            ideltax=1./cos(yusn(j+1))
            ideltax=min(iusn/2-1,ideltax)
            IF(MOD(IDELTAX,2).EQ.0)THEN
              IDELTAX=IDELTAX+1
            ENDIF
            IHALPH=(IDELTAX-1)/2 
c           print *,' ideltax=',ideltax
         IF(ideltax.eq.1)THEN
            DO i=1,iusn
               zusnfi(i+iext,j+1)=entree(i,j)
            ENDDO   
         ELSE
            DO i=1,ihalph
               DO ii=1,i+ihalph
               zusnfi(i+iext,j+1)=zusnfi(i+iext,j+1)+entree(ii,j)
               ENDDO
               DO ii=ihalph-i,0,-1
               zusnfi(i+iext,j+1)=zusnfi(i+iext,j+1)+entree(iusn-ii,j)
               ENDDO  
               zusnfi(i+iext,j+1)=zusnfi(i+iext,j+1)/real(ideltax)
            ENDDO   
            DO i=iusn-ihalph+1,iusn
               DO ii = i-ihalph,iusn
               zusnfi(i+iext,j+1)=zusnfi(i+iext,j+1)+entree(ii,j)
               ENDDO 
               DO ii = 1,ihalph+i-iusn
               zusnfi(i+iext,j+1)=zusnfi(i+iext,j+1)+entree(ii,j)
               ENDDO
               zusnfi(i+iext,j+1)=zusnfi(i+iext,j+1)/real(ideltax)
            ENDDO
            DO i=ihalph+1,iusn-ihalph
               DO ii=-ihalph,ihalph
               zusnfi(i+iext,j+1)=zusnfi(i+iext,j+1)+entree(i+ii,j)
               ENDDO
               zusnfi(i+iext,j+1)=zusnfi(i+iext,j+1)/real(ideltax)
            ENDDO
         ENDIF
            DO i=1,iext
               zusnfi(i,j+1)=zusnfi(iusn-iext+i,j+1)
               zusnfi(i+iusn+iext,j+1)=zusnfi(i,j+1)
            ENDDO
      ENDDO
c  
c Calculer les limites des zones des nouveaux points
c
      a(1) = x(1) - (x(2)-x(1))/2.0
      b(1) = (x(1)+x(2))/2.0
      DO i = 2, imar-1
         a(i) = b(i-1)
         b(i) = (x(i)+x(i+1))/2.0
      ENDDO
      a(imar) = b(imar-1)
      b(imar) = x(imar) + (x(imar)-x(imar-1))/2.0

      c(1) = y(1) - (y(2)-y(1))/2.0
      d(1) = (y(1)+y(2))/2.0
      DO j = 2, jmar-1
         c(j) = d(j-1)
         d(j) = (y(j)+y(j+1))/2.0
      ENDDO
      c(jmar) = d(jmar-1)
      d(jmar) = y(jmar) + (y(jmar)-y(jmar-1))/2.0
c
c      quelques initialisations:
      print*,'OKM1'
c
      DO i = 1, imar
      DO j = 1, jmar
         weight(i,j) = 0.0
         zxtzx(i,j) = 0.0
         zytzy(i,j) = 0.0
         zxtzy(i,j) = 0.0
         ztz(i,j) = 0.0
         zmea(i,j) = 0.0
         zstd(i,j)=0.0
      ENDDO
      ENDDO
c
c  calculs des correlations de pentes sur la grille de l'USN.
c
         DO j = 2,jusn+1 
         DO i = 1, iusn+2*iext
            zytzyusn(i,j)=0.0
            zxtzxusn(i,j)=0.0
            zxtzyusn(i,j)=0.0
         ENDDO
         ENDDO


         DO j = 2,jusn+1 
            zdeltax=zdeltay*cos(yusn(j))
         DO i = 2, iusn+2*iext-1
            zytzyusn(i,j)=(zusn(i,j+1)-zusn(i,j-1))**2/zdeltay**2
            zxtzxusn(i,j)=(zusnfi(i+1,j)-zusnfi(i-1,j))**2/zdeltax**2
            zxtzyusn(i,j)=(zusn(i,j+1)-zusn(i,j-1))/zdeltay
     *                   *(zusnfi(i+1,j)-zusnfi(i-1,j))/zdeltax
         ENDDO

         ENDDO

 

      print*,'OK0'
c
c  sommations des differentes quantites definies precedemment
c  sur une grille du modele.
c 
      zleny=xpi/real(jusn)*rad
      xincr=xpi/2./real(jusn)
       DO ii = 1, imar
       DO jj = 1, jmar
c        PRINT *,' iteration ii jj:',ii,jj
         DO j = 2,jusn+1 
c         DO j = 3,jusn 
            zlenx=zleny*cos(yusn(j))
            zdeltax=zdeltay*cos(yusn(j))
            zbordnor=(c(jj)-yusn(j)+xincr)*rad
            zbordsud=(yusn(j)-d(jj)+xincr)*rad
            weighy=amax1(0.,
     *             amin1(zbordnor,zbordsud,zleny))
         IF(weighy.ne.0)THEN
         DO i = 2, iusn+2*iext-1
            zbordest=(xusn(i)-a(ii)+xincr)*rad*cos(yusn(j))
            zbordoue=(b(ii)+xincr-xusn(i))*rad*cos(yusn(j))
            weighx=amax1(0.,
     *             amin1(zbordest,zbordoue,zlenx))
            IF(weighx.ne.0)THEN
            weight(ii,jj)=weight(ii,jj)+weighx*weighy
            zxtzx(ii,jj)=zxtzx(ii,jj)+zxtzxusn(i,j)*weighx*weighy
            zytzy(ii,jj)=zytzy(ii,jj)+zytzyusn(i,j)*weighx*weighy
            zxtzy(ii,jj)=zxtzy(ii,jj)+zxtzyusn(i,j)*weighx*weighy
            ztz(ii,jj)  =ztz(ii,jj)  +zusn(i,j)*zusn(i,j)*weighx*weighy
            zmea(ii,jj) =zmea(ii,jj)+zusn(i,j)*weighx*weighy
            ENDIF
         ENDDO
         ENDIF
         ENDDO
       ENDDO
       ENDDO
c
c  calculs des differents parametres necessaires au programme
c  de parametrisation de l'orographie a l'echelle moyenne:
c
      zmaxmea=0.
      zmaxstd=0.
      zmaxsig=0.
      zmaxgam=0.
      zmaxthe=0.
      zminthe=0.
c     print 100,' '
c100  format(1X,A1,'II JJ',4X,'H',8X,'SD',8X,'SI',3X,'GA',3X,'TH') 
       print*,'OK1'
       DO ii = 1, imar
       DO jj = 1, jmar
c       print*,'ok0'
         IF (weight(ii,jj) .NE. 0.0) THEN
c  Orography moyenne:
c         print*,'ok1'
           zmea (ii,jj)=zmea (ii,jj)/weight(ii,jj)
           zxtzx(ii,jj)=zxtzx(ii,jj)/weight(ii,jj)
           zytzy(ii,jj)=zytzy(ii,jj)/weight(ii,jj)
           zxtzy(ii,jj)=zxtzy(ii,jj)/weight(ii,jj)
           ztz(ii,jj)  =ztz(ii,jj)/weight(ii,jj)
c         print*,'ok2'
c  Deviation standard:
           zstd(ii,jj)=sqrt(amax1(0.,ztz(ii,jj)-zmea(ii,jj)**2))
c  Coefficients K, L et M:
           xk=(zxtzx(ii,jj)+zytzy(ii,jj))/2.
           xl=(zxtzx(ii,jj)-zytzy(ii,jj))/2.
           xm=zxtzy(ii,jj)
           xp=xk-sqrt(xl**2+xm**2)
           xq=xk+sqrt(xl**2+xm**2)
           xw=1.e-8
           if(xp.le.xw) xp=0.
           if(xq.le.xw) xq=xw
           if(abs(xm).le.xw) xm=xw*sign(1.,xm)
c          print*,'ok3'
c pente: 
           zsig(ii,jj)=sqrt(xq)
c           zsig(ii,jj)=sqrt(2.*xk)
c isotropy:
           zgam(ii,jj)=xp/xq
c angle theta:
           zthe(ii,jj)=57.29577951*atan2(xm,xl)/2.

c          print 101,ii,jj,
c    *           zmea(ii,jj),zstd(ii,jj),zsig(ii,jj),zgam(ii,jj),
c    *           zthe(ii,jj)
c101  format(1x,2(1x,i2),2(1x,f7.1),1x,f7.4,2x,f4.2,1x,f5.1)     
c          print*,'ok4'
         ELSE
c           PRINT*, 'probleme,ii,jj=', ii,jj
c          print*,'ok1b'
         ENDIF
      zmaxmea=amax1(zmea(ii,jj),zmaxmea)
c         print*,'oka'
      zmaxstd=amax1(zstd(ii,jj),zmaxstd)
c         print*,'okb'
      zmaxsig=amax1(zsig(ii,jj),zmaxsig)
c         print*,'okc'
      zmaxgam=amax1(zgam(ii,jj),zmaxgam)
c         print*,'okd'
      zmaxthe=amax1(zthe(ii,jj),zmaxthe)
c         print*,'oke'
      zminthe=amin1(zthe(ii,jj),zminthe)
c      print*,'ok5'
       ENDDO
       ENDDO

      print *,'  MEAN ORO:',zmaxmea
	  print *,'  ST. DEV.:',zmaxstd
      print *,'  PENTE:',zmaxsig
      print *,' ANISOTROP:',zmaxgam
      print *,'  ANGLE:',zminthe,zmaxthe	
      
C
c  On passe ce donnees sur la grille dite physique....(?)
c  On met gamma et theta a 1. et 0. aux poles ou ces quantites
c  n'ont pas vraiment de sens
c
      DO jj=1,jmar
      zmea(imar+1,jj)=zmea(1,jj)
      zstd(imar+1,jj)=zstd(1,jj)
      zsig(imar+1,jj)=zsig(1,jj)
      zgam(imar+1,jj)=zgam(1,jj)
      zthe(imar+1,jj)=zthe(1,jj)
      ENDDO


      zmeanor=0.0
      zmeasud=0.0
      zstdnor=0.0
      zstdsud=0.0
      zsignor=0.0
      zsigsud=0.0
      zweinor=0.0
      zweisud=0.0

      DO ii=1,imar
      zweinor=zweinor+              weight(ii,   1)
      zweisud=zweisud+              weight(ii,jmar)
      zmeanor=zmeanor+zmea(ii,   1)*weight(ii,   1)
      zmeasud=zmeasud+zmea(ii,jmar)*weight(ii,jmar)
      zstdnor=zstdnor+zstd(ii,   1)*weight(ii,   1)
      zstdsud=zstdsud+zstd(ii,jmar)*weight(ii,jmar)
      zsignor=zsignor+zsig(ii,   1)*weight(ii,   1)
      zsigsud=zsigsud+zsig(ii,jmar)*weight(ii,jmar)
      ENDDO

      DO ii=1,imar+1
      zmea(ii,   1)=zmeanor/zweinor
      zmea(ii,jmar)=zmeasud/zweisud
      zstd(ii,   1)=zstdnor/zweinor
      zstd(ii,jmar)=zstdsud/zweisud
      zsig(ii,   1)=zsignor/zweinor
      zsig(ii,jmar)=zsigsud/zweisud
      zgam(ii,   1)=1.
      zgam(ii,jmar)=1.
      zthe(ii,   1)=0.
      zthe(ii,jmar)=0.
      ENDDO


      RETURN
      END
