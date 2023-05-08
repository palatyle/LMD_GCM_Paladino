c-----------------------------------------------------------------------
c   INCLUDE 'comvert.h'

      COMMON/comvert/ap(llm+1),bp(llm+1),presnivs(llm),dpres(llm) ,
     .               pa,preff,nivsigs(llm),nivsig(llm+1),
     .               sig(llm+1),ds(llm),aps(llm),bps(llm),pseudoalt(llm)

      REAL ap,bp,presnivs,dpres,pa,preff,nivsigs,nivsig
c Mars Ce qui suit vient de gcm
      REAL sig,ds,aps,bps,pseudoalt

c-----------------------------------------------------------------------
