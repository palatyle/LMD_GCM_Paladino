      SUBROUTINE to36(x,nbre,xmin,xmax,iunit,icase)
Ccray IMPLICIT NONE
      IMPLICIT NONE
c=======================================================================
c-- x    : le champ a compacter
c-- nbre : la taille du champ x
c-- iunit: numero logique du fichier sortie
c-- icase: nombre de cases a utiliser (qui determine la precision)
      INTEGER iunit, nbre, icase
      REAL x(nbre) 
c=======================================================================
      INTEGER ibase, maxnbr, ifact1, ifact2, ifact3
      PARAMETER (ibase=36, maxnbr=100000)
      PARAMETER (ifact1 = ibase**1-1)
      PARAMETER (ifact2 = ibase**2-1)
      PARAMETER (ifact3 = ibase**3-1)
c=======================================================================
      CHARACTER*1 word1(0:ifact1)
      COMMON /elbat1/ word1
      CHARACTER*1 chax1(maxnbr)
c======================================================================
      CHARACTER*2 word2(0:ifact2)
      COMMON /elbat2/ word2
      CHARACTER*2 chax2(maxnbr)
c======================================================================
      CHARACTER*3 word3(0:ifact3)
      COMMON /elbat3/ word3
      CHARACTER*3 chax3(maxnbr)
c=====================================================================
      INTEGER i, inte
      REAL xmin, xmax
      EXTERNAL ini36
c=====================================================================
c
c--- initialisation pour le codage
c la routine ini36 est appelee une seule fois
c     IF (word1(36) .NE. 'Z' ) THEN
        CALL ini36
c     ENDIF

      IF (nbre .GT. maxnbr) THEN
        PRINT*, 'Nombre elements trop grand', nbre, maxnbr
        STOP
      ENDIF
c
      WRITE(iunit,*) nbre, ibase, icase, '  (nbre,ibase,icase)'
      WRITE(iunit,*) xmin
      WRITE(iunit,*) xmax
c
      IF(xmin .EQ. xmax)THEN
      WRITE(iunit,'(a)')'Ce champ est constant'
      GO TO 999
      ENDIF 
c
      IF (icase .EQ. 1) THEN
        DO 10 i = 1, nbre
          inte = NINT( FLOAT(ifact1)/(xmax-xmin) * (x(i)-xmin) )
          chax1(i) = word1(inte)
   10   CONTINUE
        WRITE(iunit,91) (chax1(i),i=1,nbre)
      ELSE IF (icase .EQ. 2) THEN
        DO 20 i = 1, nbre
          inte = NINT( FLOAT(ifact2)/(xmax-xmin) * (x(i)-xmin) )
          chax2(i) = word2(inte)
   20   CONTINUE
        WRITE(iunit,92) (chax2(i),i=1,nbre)
      ELSE IF (icase.eq.3) THEN
        DO 30 i = 1, nbre
          inte = NINT( FLOAT(ifact3)/(xmax-xmin) * (x(i)-xmin) )
          chax3(i) = word3(inte)
   30   CONTINUE
        WRITE(iunit,93) (chax3(i),i=1,nbre)
      ELSE
        PRINT*, 'Verfier icase (1,2,3) dans les arguments de to36'
        STOP
      ENDIF
c
   91 FORMAT(120a1)
   92 FORMAT(64a2)
   93 FORMAT(40a3)
c
  999 RETURN
      END 

