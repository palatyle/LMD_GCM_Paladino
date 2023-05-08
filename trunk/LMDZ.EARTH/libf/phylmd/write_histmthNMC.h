!
! $Id: write_histmthNMC.h 1403 2010-07-01 09:02:53Z fairhead $
!
      IF (ok_histNMC(1)) THEN
c
       ndex3d = 0
       itau_w = itau_phy + itap
ccc
c  Champs interpolles sur des niveaux de pression du NMC
c
c     PARAMETER(nout=3) 
c nout=1 : in=pdtphys,    out=mth
c nout=2 : in=pdtphys,    out=day
c nout=3 : in=pdtphys,    out=hf
ccc
       CALL histwrite_phy(nid_mthnmc,"tnondef",itau_w,tnondef(:,:,1))
c
       CALL histwrite_phy(nid_mthnmc,"ta",itau_w,twriteSTD(:,:,1))
c
       CALL histwrite_phy(nid_mthnmc,"zg",itau_w,phiwriteSTD(:,:,1))
c
       CALL histwrite_phy(nid_mthnmc,"hus",itau_w,qwriteSTD(:,:,1))
c
       CALL histwrite_phy(nid_mthnmc,"hur",itau_w,rhwriteSTD(:,:,1))
c
       CALL histwrite_phy(nid_mthnmc,"ua",itau_w,uwriteSTD(:,:,1))
c
       CALL histwrite_phy(nid_mthnmc,"va",itau_w,vwriteSTD(:,:,1))
c
       CALL histwrite_phy(nid_mthnmc,"wap",itau_w,wwriteSTD(:,:,1))
c
       DO k=1, nlevSTD
        DO i=1, klon
         IF(tnondef(i,k,1).NE.missing_val) THEN
          zx_tmp_fiNC(i,k) = (100.*tnondef(i,k,1))/freq_moyNMC(1)
         ELSE
          zx_tmp_fiNC(i,k) = missing_val
         ENDIF
        ENDDO
       ENDDO !k=1, nlevSTD
c
       CALL histwrite_phy(nid_mthnmc,"psbg",itau_w,zx_tmp_fiNC)
c
       CALL histwrite_phy(nid_mthnmc,"uv",itau_w,uvsumSTD(:,:,1))
c
       CALL histwrite_phy(nid_mthnmc,"vq",itau_w,vqsumSTD(:,:,1))
c
       CALL histwrite_phy(nid_mthnmc,"vT",itau_w,vTsumSTD(:,:,1))
c
       CALL histwrite_phy(nid_mthnmc,"wq",itau_w,wqsumSTD(:,:,1))
c
       CALL histwrite_phy(nid_mthnmc,"vphi",itau_w,vphisumSTD(:,:,1))
c
       CALL histwrite_phy(nid_mthnmc,"wT",itau_w,wTsumSTD(:,:,1))
c
       CALL histwrite_phy(nid_mthnmc,"uxu",itau_w,u2sumSTD(:,:,1))
c
       CALL histwrite_phy(nid_mthnmc,"vxv",itau_w,v2sumSTD(:,:,1))
c
       CALL histwrite_phy(nid_mthnmc,"TxT",itau_w,T2sumSTD(:,:,1))
c
       DO k=1, nlevSTD
        DO i=1, klon
         IF(O3sumSTD(i,k,1).NE.missing_val) THEN
          zx_tmp_fiNC(i,k) = O3sumSTD(i,k,1) * 1.e+9
         ELSE
          zx_tmp_fiNC(i,k) = missing_val
         ENDIF
        ENDDO
       ENDDO !k=1, nlevSTD
       CALL histwrite_phy(nid_mthnmc,"tro3",itau_w,
     $ zx_tmp_fiNC)
c
       if (read_climoz == 2) THEN
       DO k=1, nlevSTD
        DO i=1, klon
         IF(O3daysumSTD(i,k,1).NE.missing_val) THEN
          zx_tmp_fiNC(i,k) = O3daysumSTD(i,k,1) * 1.e+9
         ELSE
          zx_tmp_fiNC(i,k) = missing_val
         ENDIF
        ENDDO
       ENDDO !k=1, nlevSTD
c
        CALL histwrite_phy(nid_mthnmc,"tro3_daylight",itau_w,
     $  zx_tmp_fiNC)
       endif 
c
      if (ok_sync) then
c$OMP MASTER
        call histsync(nid_mthnmc)
c$OMP END MASTER
      endif
c
      ENDIF !(ok_histNMC(1)) THEN
