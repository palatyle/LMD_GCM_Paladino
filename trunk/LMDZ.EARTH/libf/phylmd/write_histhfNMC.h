!
! $Header$
!
      IF (ok_histNMC(3)) THEN
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
       CALL histwrite_phy(nid_hfnmc,"tnondef",itau_w,tnondef(:,:,3))
c
       CALL histwrite_phy(nid_hfnmc,"ta",itau_w,twriteSTD3)
c
       CALL histwrite_phy(nid_hfnmc,"zg",itau_w,phiwriteSTD3)
c
       CALL histwrite_phy(nid_hfnmc,"hus",itau_w,qwriteSTD3)
c
       CALL histwrite_phy(nid_hfnmc,"hur",itau_w,rhwriteSTD3)
c
       CALL histwrite_phy(nid_hfnmc,"ua",itau_w,uwriteSTD3)
c
       CALL histwrite_phy(nid_hfnmc,"va",itau_w,vwriteSTD3)
c
       CALL histwrite_phy(nid_hfnmc,"wap",itau_w,wwriteSTD3)
c
       DO k=1, nlevSTD
        DO i=1, klon
         IF(tnondef(i,k,3).NE.missing_val) THEN
          zx_tmp_fiNC(i,k) = (100.*tnondef(i,k,3))/freq_moyNMC(3)
         ELSE
          zx_tmp_fiNC(i,k) = missing_val
         ENDIF
        ENDDO
       ENDDO !k=1, nlevSTD
c
       CALL histwrite_phy(nid_hfnmc,"psbg",itau_w,zx_tmp_fiNC)
c
       CALL histwrite_phy(nid_hfnmc,"uv",itau_w,uvsumSTD(:,:,3))
c
       CALL histwrite_phy(nid_hfnmc,"vq",itau_w,vqsumSTD(:,:,3))
c
       CALL histwrite_phy(nid_hfnmc,"vT",itau_w,vTsumSTD(:,:,3))
c
       CALL histwrite_phy(nid_hfnmc,"wq",itau_w,wqsumSTD(:,:,3))
c
       CALL histwrite_phy(nid_hfnmc,"vphi",itau_w,vphisumSTD(:,:,3))
c
       CALL histwrite_phy(nid_hfnmc,"wT",itau_w,wTsumSTD(:,:,3))
c
       CALL histwrite_phy(nid_hfnmc,"uxu",itau_w,u2sumSTD(:,:,3))
c
       CALL histwrite_phy(nid_hfnmc,"vxv",itau_w,v2sumSTD(:,:,3))
c
       CALL histwrite_phy(nid_hfnmc,"TxT",itau_w,T2sumSTD(:,:,3))
c
c     ENDIF !type_run
c
      if (ok_sync) then
c$OMP MASTER
        call histsync(nid_hfnmc)
c$OMP END MASTER
      endif
c
      ENDIF !      (ok_histNMC(3)) THEN
