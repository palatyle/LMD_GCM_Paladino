c-----------------------------------------------------------------------
c INCLUDE comdissipn.h

      REAL  tetaudiv, tetaurot, tetah, cdivu, crot, cdivh
c
      COMMON/comdissipn/ tetaudiv(llm),tetaurot(llm),tetah(llm)   ,
     1                        cdivu,      crot,         cdivh

c
c    Les parametres de ce common proviennent des calculs effectues dans 
c             Inidissip  .
c
c-----------------------------------------------------------------------
