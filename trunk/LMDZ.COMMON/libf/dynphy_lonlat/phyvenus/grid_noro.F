!
! $Id: grid_noro.F 1442 2010-10-18 08:31:31Z jghattas $
!
c
c
      SUBROUTINE grid_noro(imdep, jmdep, xdata, ydata, zdata,
     .             imar, jmar, x, y,
     .             zphi,zmea,zstd,zsig,zgam,zthe,
     .             zpic,zval)
c=======================================================================
c (F. Lott) (voir aussi z.x. Li, A. Harzallah et L. Fairhead)
c
c      Compute the Parameters of the SSO scheme as described in
c      LOTT & MILLER (1997) and LOTT(1999).
c      Target points are on a rectangular grid:
c      iim+1 latitudes including North and South Poles;
c      jjm+1 longitudes, with periodicity jjm+1=1.
c      aux poles.  At the poles the fields value is repeated
c      jjm+1 time.
c      The parameters a,b,c,d represent the limite of the target
c      gridpoint region. The means over this region are calculated
c      from USN data, ponderated by a weight proportional to the 
c      surface occupated by the data inside the model gridpoint area.
c      In most circumstances, this weight is the ratio between the
c      surface of the USN gridpoint area and the surface of the
c      model gridpoint area. 
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
c        imdep, jmdep: dimensions X and Y input field
c        xdata, ydata: coordinates X and Y input field
c        zdata: Input field
c OUTPUT:
c        imar, jmar: dimensions X and Y Output field
c        x, y: ccordinates  X and Y Output field.
c             zmea:  Mean orographie   
c             zstd:  Standard deviation
c             zsig:  Slope
c             zgam:  Anisotropy
c             zthe:  Orientation of the small axis
c             zpic:  Maximum altitude
c             zval:  Minimum altitude
C=======================================================================

      IMPLICIT none 
      
#include "YOMCST.h"

      INTEGER imdep, jmdep
      REAL xdata(imdep),ydata(jmdep) 
      REAL zdata(imdep,jmdep)
c
      INTEGER imar, jmar
c parametres lies au fichier d entree... A documenter...
      integer iext
      parameter(iext=216)
      REAL xusn(imdep+2*iext),yusn(jmdep+2)
      REAL zusn(imdep+2*iext,jmdep+2)
  
c local var
      real zdeltax,zdeltay,zlenx,zleny,xincr
      real zbordnor,zbordsud,zbordest,zbordoue,weighx,weighy
      real zllmmea,zllmstd,zllmsig,zllmgam,zllmpic,zllmval,zllmthe
      real zminthe,xk,xl,xm,xp,xq,xw
      real zmeanor,zmeasud,zstdnor,zstdsud,zsignor,zsigsud
      real zweinor,zweisud,zpicnor,zpicsud,zvalnor,zvalsud
      integer i,j,ii,jj

C INTERMEDIATE FIELDS  (CORRELATIONS OF OROGRAPHY GRADIENT)

      REAL ztz(imar+1,jmar),zxtzx(imar+1,jmar)
      REAL zytzy(imar+1,jmar),zxtzy(imar+1,jmar)
      REAL weight(imar+1,jmar)

C CORRELATIONS OF USN OROGRAPHY GRADIENTS

      REAL zxtzxusn(imdep+2*iext,jmdep+2)
      REAL zytzyusn(imdep+2*iext,jmdep+2)
      REAL zxtzyusn(imdep+2*iext,jmdep+2)
      REAL x(imar+1),y(jmar),zphi(imar+1,jmar)
      REAL zmea(imar+1,jmar),zstd(imar+1,jmar)
      REAL zsig(imar+1,jmar),zgam(imar+1,jmar),zthe(imar+1,jmar)
      REAL zpic(imar+1,jmar),zval(imar+1,jmar)
      real num_tot(2200,1100),num_lan(2200,1100)

      REAL a(2200),b(2200),c(1100),d(1100)

c pas defini puisque pas de physique dans newstart...
      RPI=2.*ASIN(1.)
      RA=6051300.

      print *,' parametres de l orographie a l echelle sous maille' 

      zdeltay=2.*RPI/REAL(jmdep)*RA
c
c  quelques tests de dimensions:
c    
c
      IF (imar.GT.2200 .OR. jmar.GT.1100) THEN
         PRINT*, 'imar or jmar too big', imar, jmar
         CALL ABORT
      ENDIF

c      print *,'xdata:',xdata
c      print *,'ydata:',ydata
c      print *,'x:',x
c      print *,'y:',y
c
C  EXTENSION OF THE USN DATABASE TO POCEED COMPUTATIONS AT
C  BOUNDARIES:
c
      DO j=1,jmdep
        yusn(j+1)=ydata(j)
      DO i=1,imdep
        zusn(i+iext,j+1)=zdata(i,j)
        xusn(i+iext)=xdata(i)
      ENDDO
      DO i=1,iext
        zusn(i,j+1)=zdata(imdep-iext+i,j)
        xusn(i)=xdata(imdep-iext+i)-2.*RPI
        zusn(imdep+iext+i,j+1)=zdata(i,j)
        xusn(imdep+iext+i)=xdata(i)+2.*RPI
      ENDDO
      ENDDO

        yusn(1)=ydata(1)+(ydata(1)-ydata(2))
        yusn(jmdep+2)=ydata(jmdep)+(ydata(jmdep)-ydata(jmdep-1))
       DO i=1,imdep/2+iext
        zusn(i,1)=zusn(i+imdep/2,2)
        zusn(i+imdep/2+iext,1)=zusn(i,2)
        zusn(i,jmdep+2)=zusn(i+imdep/2,jmdep+1)
        zusn(i+imdep/2+iext,jmdep+2)=zusn(i,jmdep+1)
       ENDDO
c  
c COMPUTE LIMITS OF MODEL GRIDPOINT AREA
C     ( REGULAR GRID)
c
      a(1) = x(1) - (x(2)-x(1))/2.0
      b(1) = (x(1)+x(2))/2.0
      DO i = 2, imar
         a(i) = b(i-1)
         b(i) = (x(i)+x(i+1))/2.0
      ENDDO
      a(imar+1) = b(imar)
      b(imar+1) = x(imar+1) + (x(imar+1)-x(imar))/2.0

      c(1) = y(1) - (y(2)-y(1))/2.0
      d(1) = (y(1)+y(2))/2.0
      DO j = 2, jmar-1
         c(j) = d(j-1)
         d(j) = (y(j)+y(j+1))/2.0
      ENDDO
      c(jmar) = d(jmar-1)
      d(jmar) = y(jmar) + (y(jmar)-y(jmar-1))/2.0
c
c  initialisations:
c
      DO i = 1, imar+1
      DO j = 1, jmar
         weight(i,j) = 0.0
         zxtzx(i,j)  = 0.0
         zytzy(i,j)  = 0.0
         zxtzy(i,j)  = 0.0
         ztz(i,j)    = 0.0
         zmea(i,j)   = 0.0
         zpic(i,j)  =-1.E+10
         zval(i,j)  = 1.E+10
      ENDDO
      ENDDO
c
c  COMPUTE SLOPES CORRELATIONS ON USN GRID
c
         DO j = 1,jmdep+2 
         DO i = 1, imdep+2*iext
            zytzyusn(i,j)=0.0
            zxtzxusn(i,j)=0.0
            zxtzyusn(i,j)=0.0
         ENDDO
         ENDDO


         DO j = 2,jmdep+1 
            zdeltax=zdeltay*cos(yusn(j))
         DO i = 2, imdep+2*iext-1
            zytzyusn(i,j)=(zusn(i,j+1)-zusn(i,j-1))**2/zdeltay**2
            zxtzxusn(i,j)=(zusn(i+1,j)-zusn(i-1,j))**2/zdeltax**2
            zxtzyusn(i,j)=(zusn(i,j+1)-zusn(i,j-1))/zdeltay
     *                   *(zusn(i+1,j)-zusn(i-1,j))/zdeltax
         ENDDO
         ENDDO
c
c  SUMMATION OVER GRIDPOINT AREA
c 
      zleny=RPI/REAL(jmdep)*RA
      xincr=RPI/2./REAL(jmdep)
       DO ii = 1, imar+1
       DO jj = 1, jmar
       num_tot(ii,jj)=0.
       num_lan(ii,jj)=0.
c        PRINT *,' iteration ii jj:',ii,jj
         DO j = 2,jmdep+1 
c         DO j = 3,jmdep 
            zlenx=zleny*cos(yusn(j))
            zdeltax=zdeltay*cos(yusn(j))
            zbordnor=(c(jj)-yusn(j)+xincr)*RA
            zbordsud=(yusn(j)-d(jj)+xincr)*RA
            weighy=AMAX1(0.,
     *             amin1(zbordnor,zbordsud,zleny))
         IF(weighy.ne.0)THEN
         DO i = 2, imdep+2*iext-1
            zbordest=(xusn(i)-a(ii)+xincr)*RA*cos(yusn(j))
            zbordoue=(b(ii)+xincr-xusn(i))*RA*cos(yusn(j))
            weighx=AMAX1(0.,
     *             amin1(zbordest,zbordoue,zlenx))
            IF(weighx.ne.0)THEN
            num_tot(ii,jj)=num_tot(ii,jj)+1.0
            if(zusn(i,j).ge.1.)num_lan(ii,jj)=num_lan(ii,jj)+1.0
            weight(ii,jj)=weight(ii,jj)+weighx*weighy
            zxtzx(ii,jj)=zxtzx(ii,jj)+zxtzxusn(i,j)*weighx*weighy
            zytzy(ii,jj)=zytzy(ii,jj)+zytzyusn(i,j)*weighx*weighy
            zxtzy(ii,jj)=zxtzy(ii,jj)+zxtzyusn(i,j)*weighx*weighy
            ztz(ii,jj)  =ztz(ii,jj)  +zusn(i,j)*zusn(i,j)*weighx*weighy
c mean
            zmea(ii,jj) =zmea(ii,jj)+zusn(i,j)*weighx*weighy
c peacks
            zpic(ii,jj)=amax1(zpic(ii,jj),zusn(i,j))
c valleys
            zval(ii,jj)=amin1(zval(ii,jj),zusn(i,j))
            ENDIF
         ENDDO
         ENDIF
         ENDDO
       ENDDO
       ENDDO
c
c  COMPUTE PARAMETERS NEEDED BY THE LOTT & MILLER (1997) AND
C  LOTT (1999) SSO SCHEME.
c
      zllmmea=0.
      zllmstd=0.
      zllmsig=0.
      zllmgam=0.
      zllmpic=0.
      zllmval=0.
      zllmthe=0.
      zminthe=0.
c     print 100,' '
c100  format(1X,A1,'II JJ',4X,'H',8X,'SD',8X,'SI',3X,'GA',3X,'TH') 
       DO ii = 1, imar+1
       DO jj = 1, jmar
         IF (weight(ii,jj) .NE. 0.0) THEN
c  Mean Orography:
           zmea (ii,jj)=zmea (ii,jj)/weight(ii,jj)
           zxtzx(ii,jj)=zxtzx(ii,jj)/weight(ii,jj)
           zytzy(ii,jj)=zytzy(ii,jj)/weight(ii,jj)
           zxtzy(ii,jj)=zxtzy(ii,jj)/weight(ii,jj)
           ztz(ii,jj)  =ztz(ii,jj)/weight(ii,jj)
c  Standard deviation:
           zstd(ii,jj)=sqrt(AMAX1(0.,ztz(ii,jj)-zmea(ii,jj)**2))
         ELSE
            PRINT*, 'probleme,ii,jj=', ii,jj
         ENDIF
       ENDDO
       ENDDO

C CORRECT VALUES OF HORIZONTAL SLOPE NEAR THE POLES:

       DO ii = 1, imar+1
         zxtzx(ii,1)=zxtzx(ii,2)
         zxtzx(ii,jmar)=zxtzx(ii,jmar-1)
         zxtzy(ii,1)=zxtzy(ii,2)
         zxtzy(ii,jmar)=zxtzy(ii,jmar-1)
         zytzy(ii,1)=zytzy(ii,2)
         zytzy(ii,jmar)=zytzy(ii,jmar-1)
       ENDDO

C  FILTERS TO SMOOTH OUT FIELDS FOR INPUT INTO SSO SCHEME.

C  FIRST FILTER, MOVING AVERAGE OVER 9 POINTS.

       CALL MVA9(zmea,imar+1,jmar)
       CALL MVA9(zstd,imar+1,jmar)
       CALL MVA9(zpic,imar+1,jmar)
       CALL MVA9(zval,imar+1,jmar)
       CALL MVA9(zxtzx,imar+1,jmar)
       CALL MVA9(zxtzy,imar+1,jmar) 
       CALL MVA9(zytzy,imar+1,jmar)

       DO ii = 1, imar
       DO jj = 1, jmar
         IF (weight(ii,jj) .NE. 0.0) THEN
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
c slope: 
           zsig(ii,jj)=sqrt(xq)
c isotropy:
           zgam(ii,jj)=xp/xq
c angle theta:
           zthe(ii,jj)=57.29577951*atan2(xm,xl)/2.
           zphi(ii,jj)=zmea(ii,jj)
           !
           zmea(ii,jj)=zmea(ii,jj)
           zpic(ii,jj)=zpic(ii,jj)
           zval(ii,jj)=zval(ii,jj)
           zstd(ii,jj)=zstd(ii,jj)
c          print 101,ii,jj,
c    *           zmea(ii,jj),zstd(ii,jj),zsig(ii,jj),zgam(ii,jj),
c    *           zthe(ii,jj)
c101  format(1x,2(1x,i2),2(1x,f7.1),1x,f7.4,2x,f4.2,1x,f5.1)     
         ELSE
c           PRINT*, 'probleme,ii,jj=', ii,jj
         ENDIF
      zllmmea=AMAX1(zmea(ii,jj),zllmmea)
      zllmstd=AMAX1(zstd(ii,jj),zllmstd)
      zllmsig=AMAX1(zsig(ii,jj),zllmsig)
      zllmgam=AMAX1(zgam(ii,jj),zllmgam)
      zllmthe=AMAX1(zthe(ii,jj),zllmthe)
      zminthe=amin1(zthe(ii,jj),zminthe)
      zllmpic=AMAX1(zpic(ii,jj),zllmpic)
      zllmval=AMAX1(zval(ii,jj),zllmval)
       ENDDO
       ENDDO
      print *,'  MEAN ORO:',zllmmea
      print *,'  ST. DEV.:',zllmstd
      print *,'  PENTE:',zllmsig
      print *,' ANISOTROP:',zllmgam
      print *,'  ANGLE:',zminthe,zllmthe
      print *,'  pic:',zllmpic
      print *,'  val:',zllmval
      
C
c gamma and theta a 1. and 0. at poles
c
      DO jj=1,jmar
      zmea(imar+1,jj)=zmea(1,jj)
      zphi(imar+1,jj)=zphi(1,jj)
      zpic(imar+1,jj)=zpic(1,jj)
      zval(imar+1,jj)=zval(1,jj)
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
      zpicnor=0.0
      zpicsud=0.0                                   
      zvalnor=0.0
      zvalsud=0.0 

      DO ii=1,imar
      zweinor=zweinor+              weight(ii,   1)
      zweisud=zweisud+              weight(ii,jmar)
      zmeanor=zmeanor+zmea(ii,   1)*weight(ii,   1)
      zmeasud=zmeasud+zmea(ii,jmar)*weight(ii,jmar)
      zstdnor=zstdnor+zstd(ii,   1)*weight(ii,   1)
      zstdsud=zstdsud+zstd(ii,jmar)*weight(ii,jmar)
      zsignor=zsignor+zsig(ii,   1)*weight(ii,   1)
      zsigsud=zsigsud+zsig(ii,jmar)*weight(ii,jmar)
      zpicnor=zpicnor+zpic(ii,   1)*weight(ii,   1)
      zpicsud=zpicsud+zpic(ii,jmar)*weight(ii,jmar)
      zvalnor=zvalnor+zval(ii,   1)*weight(ii,   1)
      zvalsud=zvalsud+zval(ii,jmar)*weight(ii,jmar)
      ENDDO

      DO ii=1,imar+1
      zmea(ii,   1)=zmeanor/zweinor
      zmea(ii,jmar)=zmeasud/zweisud
      zphi(ii,   1)=zmeanor/zweinor
      zphi(ii,jmar)=zmeasud/zweisud
      zpic(ii,   1)=zpicnor/zweinor
      zpic(ii,jmar)=zpicsud/zweisud
      zval(ii,   1)=zvalnor/zweinor
      zval(ii,jmar)=zvalsud/zweisud
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

      SUBROUTINE MVA9(X,IMAR,JMAR)

C MAKE A MOVING AVERAGE OVER 9 GRIDPOINTS OF THE X FIELDS

      REAL X(IMAR,JMAR),XF(IMAR,JMAR)
      real WEIGHTpb(-1:1,-1:1)


      SUM=0.
      DO IS=-1,1
        DO JS=-1,1
          WEIGHTpb(IS,JS)=1./REAL((1+IS**2)*(1+JS**2))
          SUM=SUM+WEIGHTpb(IS,JS)
        ENDDO
      ENDDO
      
c     WRITE(*,*) 'MVA9 ', IMAR, JMAR
c     WRITE(*,*) 'MVA9 ', WEIGHTpb
c     WRITE(*,*) 'MVA9 SUM ', SUM
      DO IS=-1,1
        DO JS=-1,1
          WEIGHTpb(IS,JS)=WEIGHTpb(IS,JS)/SUM
        ENDDO
      ENDDO

      DO J=2,JMAR-1
        DO I=2,IMAR-1
          XF(I,J)=0.
          DO IS=-1,1
            DO JS=-1,1
              XF(I,J)=XF(I,J)+X(I+IS,J+JS)*WEIGHTpb(IS,JS)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO J=2,JMAR-1
        XF(1,J)=0.
        IS=IMAR-1
        DO JS=-1,1 
          XF(1,J)=XF(1,J)+X(IS,J+JS)*WEIGHTpb(-1,JS)
        ENDDO
        DO IS=0,1 
          DO JS=-1,1 
            XF(1,J)=XF(1,J)+X(1+IS,J+JS)*WEIGHTpb(IS,JS)
          ENDDO
        ENDDO
        XF(IMAR,J)=XF(1,J)
      ENDDO

      DO I=1,IMAR
        XF(I,1)=XF(I,2)
        XF(I,JMAR)=XF(I,JMAR-1)
      ENDDO

      DO I=1,IMAR
        DO J=1,JMAR
          X(I,J)=XF(I,J)
        ENDDO
      ENDDO

      RETURN
      END



