REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: &
       TSURF,&
       TAU_DUST,&
       SWDOWN,&
       LWDOWN,&
       SWUP,&
       LWUP,&
   PSFC,TSK
REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT)  :: &
   RTHBLTEN,RUBLTEN,RVBLTEN
