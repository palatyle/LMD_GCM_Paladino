
c
c $Header$
c
c
      ndex2d = 0
      ndex3d = 0
c
      itau_w = itau_phy + itap
c
c Champs 3D:
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, t_seri, zx_tmp_3d)
      CALL histwrite_phy(nid_hf3d,"temp",itau_w,t_seri)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, qx(1,1,ivap), zx_tmp_3d)
      CALL histwrite_phy(nid_hf3d,"ovap",itau_w,qx(:,:,ivap))
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, u_seri, zx_tmp_3d)
      CALL histwrite_phy(nid_hf3d,"vitu",itau_w,u_seri)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, v_seri, zx_tmp_3d)
      CALL histwrite_phy(nid_hf3d,"vitv",itau_w,v_seri)
      if (ok_sync) then
c$OMP MASTER
        call histsync(nid_hf3d)
c$OMP END MASTER      
      endif
