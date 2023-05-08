c-----------------------------------------------------------------------
c INCLUDE comdissnew.h


c****WRF: attention ligne trop longue

      COMMON/comdissnew/ lstardis,nitergdiv,nitergrot,niterh,
     1                   tetagdiv,tetagrot,tetatemp,coefdis 

      LOGICAL lstardis
      INTEGER nitergdiv, nitergrot, niterh
      REAL     tetagdiv, tetagrot,  tetatemp, coefdis

c
c ... Les parametres de ce common comdissnew sont  lues par defrun_new 
c              sur le fichier  run.def    ....
c
c-----------------------------------------------------------------------
