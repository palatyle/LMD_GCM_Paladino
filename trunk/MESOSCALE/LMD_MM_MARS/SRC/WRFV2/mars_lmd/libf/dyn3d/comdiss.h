c-----------------------------------------------------------------------
c INCLUDE 'comdiss.h'

      INTEGER niter
      REAL tetau(llm),tetah(llm),cdivu,cdivh,crot
      LOGICAL dissipst

      COMMON/comdissi/niter
      COMMON/comdissr/tetau,tetah,cdivu,cdivh,crot
      COMMON/comdissl/dissipst

c-----------------------------------------------------------------------
