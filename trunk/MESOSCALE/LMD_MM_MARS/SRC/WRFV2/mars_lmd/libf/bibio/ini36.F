      SUBROUTINE ini36
Ccray IMPLICIT NONE
      IMPLICIT NONE
c====================================================================
c-- cette routine fait l'initialisation pour le codage sur base 36
c-- de la routine to36. on peut utiliser soit 1 soit 2 soit 3 cases
      INTEGER i,j,k,l, ijk, jk
c====================================================================
      INTEGER ibase, ifact1, ifact2, ifact3
      PARAMETER (ibase=36)
      PARAMETER (ifact1 = ibase-1)
      PARAMETER (ifact2 = ibase*ibase-1)
      PARAMETER (ifact3 = ibase*ibase*ibase-1)
c====================================================================
      CHARACTER*1 word1(0:ifact1)
      COMMON /elbat1/ word1
      CHARACTER*2 word2(0:ifact2)
      COMMON /elbat2/ word2
      CHARACTER*3 word3(0:ifact3)
      COMMON /elbat3/ word3
c====================================================================
      CHARACTER*1 base(ibase)
      DATA (base(l),l=1,ibase)
     .     /'0','1','2','3','4','5','6','7','8','9',
     .      'A','B','C','D','E','F','G','H','I','J',
     .      'K','L','M','N','O','P','Q','R','S','T',
     .      'U','V','W','X','Y','Z'/
c===================================================================
      DO 30 i = 1, ibase
      DO 30 j = 1, ibase
      DO 30 k = 1, ibase
      ijk = ( (i-1)*ibase + j-1 )*ibase + k-1
      word3(ijk) = base(i)//base(j)//base(k)
   30 CONTINUE
      DO 20 j = 1, ibase
      DO 20 k = 1, ibase
      jk = (j-1)*ibase + k-1
      word2(jk) = base(j)//base(k)
   20 CONTINUE
      DO 10 k = 1, ibase
      word1(k-1) = base(k)
   10 CONTINUE
c===================================================================
      RETURN
      END

