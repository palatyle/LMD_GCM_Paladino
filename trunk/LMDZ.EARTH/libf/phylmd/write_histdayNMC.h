!
! $Header$
!
      IF (ok_histNMC(2)) THEN
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
      IF(lev_histdayNMC.EQ.nlevSTD) THEN
       CALL histwrite_phy(nid_daynmc,"tnondef",itau_w,tnondef(:,:,2))
       CALL histwrite_phy(nid_daynmc,"ta",itau_w,twriteSTD(:,:,2))
       CALL histwrite_phy(nid_daynmc,"zg",itau_w,phiwriteSTD(:,:,2))
       CALL histwrite_phy(nid_daynmc,"hus",itau_w,qwriteSTD(:,:,2))
       CALL histwrite_phy(nid_daynmc,"hur",itau_w,rhwriteSTD(:,:,2))
       CALL histwrite_phy(nid_daynmc,"ua",itau_w,uwriteSTD(:,:,2))
       CALL histwrite_phy(nid_daynmc,"va",itau_w,vwriteSTD(:,:,2))
       CALL histwrite_phy(nid_daynmc,"wap",itau_w,wwriteSTD(:,:,2))
      ELSE IF(lev_histdayNMC.EQ.nlevSTD8) THEN
       CALL histwrite_phy(nid_daynmc,"tnondef",itau_w,tnondefSTD8)
       CALL histwrite_phy(nid_daynmc,"ta",itau_w,twriteSTD8)
       CALL histwrite_phy(nid_daynmc,"zg",itau_w,phiwriteSTD8)
       CALL histwrite_phy(nid_daynmc,"hus",itau_w,qwriteSTD8)
       CALL histwrite_phy(nid_daynmc,"hur",itau_w,rhwriteSTD8)
       CALL histwrite_phy(nid_daynmc,"ua",itau_w,uwriteSTD8)
       CALL histwrite_phy(nid_daynmc,"va",itau_w,vwriteSTD8)
       CALL histwrite_phy(nid_daynmc,"wap",itau_w,wwriteSTD8)
      ENDIF
c
      if (ok_sync) then
c$OMP MASTER
        call histsync(nid_daynmc)
c$OMP END MASTER
      endif
c
      ENDIF ! (ok_histNMC(2)) THEN
