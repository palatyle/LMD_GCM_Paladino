!     -----------------------------------------------------------------
!*    *COMMON* *YOEGWD* - PARAMETERS FOR GRAVITY WAVE DRAG CALCULATIONS
!     -----------------------------------------------------------------
!
      integer :: NKTOPG,NTOP
      real    :: GFRCRIT,GKWAKE,GRCRIT,GVCRIT,GKDRAG,GKLIFT
      real    :: GHMAX,GRAHILO,GSIGCR,GSSEC,GTSEC,GVSEC
      real    :: TAUBS
      integer :: LEVBS
      COMMON/YOEGWD/ GFRCRIT,GKWAKE,GRCRIT,GVCRIT,GKDRAG,GKLIFT         &
     &   ,GHMAX,GRAHILO,GSIGCR,NKTOPG,NTOP,GSSEC,GTSEC,GVSEC            &
     &   ,TAUBS,LEVBS

