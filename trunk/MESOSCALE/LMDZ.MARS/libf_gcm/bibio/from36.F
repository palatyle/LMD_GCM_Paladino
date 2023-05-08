      SUBROUTINE from36(x,nbre,iunit)
      IMPLICIT NONE
c====================================================================
      INTEGER iunit, nbre
      REAL x(nbre)
c====================================================================
      INTEGER i, l
      INTEGER nbreok, ibaseok, icase
      REAL xmax, xmin, errcod
c====================================================================
      INTEGER lbord
      PARAMETER (lbord=255)
      REAL conv(0:lbord)
c====================================================================     
      INTEGER ibase, maxcas, maxnbr, ifact1, ifact2, ifact3
      PARAMETER (ibase=36, maxcas=3, maxnbr=100000)
      PARAMETER(ifact1=ibase**1-1)
      PARAMETER(ifact2=ibase**2-1)
      PARAMETER(ifact3=ibase**3-1)
      CHARACTER*1 base(ibase), fac(maxcas,maxnbr)
      INTEGER inte(maxcas)
      DATA (base(l),l=1,ibase)
     .     /'0','1','2','3','4','5','6','7','8','9',
     .      'A','B','C','D','E','F','G','H','I','J',
     .      'K','L','M','N','O','P','Q','R','S','T',
     .      'U','V','W','X','Y','Z'/
c====================================================================
      IF (nbre .GT. maxnbr) THEN
      PRINT*, 'Verifier maxnbr dans from36 (trop petite)'
      STOP
      ENDIF
      DO 10 i = 0, lbord
   10 conv(i) = 1.0e+33
c    34  c  10 conv(i) = 1.0e+49
      DO 20 i = 1, ibase
      l = ICHAR(base(i))
   20 conv(l) = FLOAT(i-1)
c====================================================================
      READ (iunit,*) nbreok, ibaseok, icase
        IF (nbre .NE. nbreok) THEN
          PRINT*, 'Incoherence (Nombre elements)', nbreok, nbre
          STOP
        ENDIF
        IF (ibase .NE. ibaseok) THEN
          PRINT*, 'Incoherence (Nombre base)', ibaseok, ibase
          STOP
        ENDIF
        IF (icase.NE.1 .AND. icase.NE.2 .AND. icase.NE.3) THEN
          PRINT*, 'Verifier icase (1,2,3) dans le fichier a lire', 
     .             icase
          STOP
        ENDIF
c====================================================================
      READ (iunit,*) xmin
      READ (iunit,*) xmax
        PRINT*,'Xmin et Xmax=', xmin, xmax
c====================================================================
      IF (xmin .EQ. xmax) THEN
        READ(iunit,'(a)')
        PRINT*, 'Ce champ est constant'
        DO 30 i = 1, nbre
   30   x(i) = xmin
        GO TO 123
      END IF
c====================================================================
      IF (icase .EQ. 1) THEN
        errcod = (xmax-xmin) / FLOAT(ifact1)
        PRINT*,'La precision du codage est entre + / - ', errcod
        READ(iunit,91)((fac(l,i),l=1,icase),i=1,nbre)
        DO 190 i = 1, nbre
        x(i) = 0.0
        DO 191 l = 1, icase
        inte(l) = ICHAR( fac(l,i) )
        x(i) = x(i) + conv ( inte(l) ) * FLOAT( ibase**(icase-l) )
  191   CONTINUE
        x(i) = xmin + x(i)/FLOAT(ifact1)*(xmax-xmin)
  190   CONTINUE
      ELSE IF (icase .EQ. 2) THEN
        errcod = (xmax-xmin) / FLOAT(ifact2)
        PRINT*,'La precision du codage est entre + / - ', errcod
        READ(iunit,92)((fac(l,i),l=1,icase),i=1,nbre)
        DO 290 i = 1, nbre
        x(i) = 0.0
        DO 291 l = 1, icase
        inte(l) = ICHAR( fac(l,i) )
        x(i) = x(i) + conv ( inte(l) ) * FLOAT( ibase**(icase-l) )
  291   CONTINUE
        x(i) = xmin + x(i)/float(ifact2)*(xmax-xmin)
  290   CONTINUE
      ELSE IF (icase .EQ. 3) THEN
        errcod = (xmax-xmin) / FLOAT(ifact3)
        PRINT*,'La precision du codage est entre + / - ', errcod
        READ(iunit,93)((fac(l,i),l=1,icase),i=1,nbre)
        DO 390 i = 1, nbre
        x(i) = 0.0
        DO 391 l = 1, icase
        inte(l) = ICHAR( fac(l,i) )
        x(i) = x(i) + conv ( inte(l) ) * FLOAT( ibase**(icase-l) )
  391   CONTINUE
        x(i) = xmin + x(i)/FLOAT(ifact3)*(xmax-xmin)
  390   CONTINUE
      ELSE
        PRINT*, 'Verfier icase (1, 2 ou 3)'
        STOP
      ENDIF
   91 FORMAT(120a1)
   92 FORMAT(64(2a1))
   93 FORMAT(40(3a1))
c====================================================================
  123 RETURN
      END
