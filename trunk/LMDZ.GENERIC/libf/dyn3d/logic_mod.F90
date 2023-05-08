MODULE logic_mod

IMPLICIT NONE  

LOGICAL purmats,forward,leapf,apphys,statcl,conser, &
     & apdiss,apdelq,saison,ecripar,fxyhypb,ysinus,hybrid,autozlevs

INTEGER iflag_phys ! ==1 if calling a physics package

END MODULE logic_mod
