!
! $Header$
!
C     -----------------------------------------------------------------
C*    *COMMON* *YOEGWD* - PARAMETERS FOR GRAVITY WAVE DRAG CALCULATIONS
C     -----------------------------------------------------------------
C
      integer NKTOPG,NSTRA
      real GFRCRIT,GKWAKE,GRCRIT,GVCRIT,GKDRAG,GKLIFT
      real GHMAX,GRAHILO,GSIGCR,GSSEC,GTSEC,GVSEC
      COMMON/YOEGWD/ GFRCRIT,GKWAKE,GRCRIT,GVCRIT,GKDRAG,GKLIFT
     *        ,GHMAX,GRAHILO,GSIGCR,NKTOPG,NSTRA,GSSEC,GTSEC,GVSEC
c$OMP THREADPRIVATE(/YOEGWD/)
C


