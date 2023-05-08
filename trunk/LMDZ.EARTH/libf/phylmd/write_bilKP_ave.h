c
c $Header$
c
      IF (ok_journe) THEN
c
      ndex2d = 0
      ndex3d = 0
c
c Champs 2D:
c
      itau_w = itau_phy + itap
c
cym      CALL gr_fi_ecrit(klev, klon,iim,jjmp1, ue_lay,zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"ue",itau_w,ue_lay)
c
cym      CALL gr_fi_ecrit(klev, klon,iim,jjmp1, ve_lay,zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"ve",itau_w,ve_lay)
c
cym      CALL gr_fi_ecrit(klev, klon,iim,jjmp1, uq_lay,zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"uq",itau_w,uq_lay)
c
cym      CALL gr_fi_ecrit(klev, klon,iim,jjmp1, vq_lay,zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"vq",itau_w,vq_lay)
c
c Champs 3D:
C
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, t_seri, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"temp",itau_w,t_seri)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, qx(1,1,ivap), zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"ovap",itau_w,qx(:,:,ivap))
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, zphi, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"geop",itau_w,zphi)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, u_seri, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"vitu",itau_w,u_seri)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, v_seri, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"vitv",itau_w,v_seri)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, omega, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"vitw",itau_w,omega)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, pplay, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"pres",itau_w,pplay)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, paprs, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"play",itau_w,paprs)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, cldliq, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"oliq",itau_w,cldliq)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_t_dyn, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dtdyn",itau_w,d_t_dyn)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_q_dyn, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dqdyn",itau_w,d_q_dyn)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_t_con, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dtcon",itau_w,d_t_con)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_u_con, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"ducon",itau_w,d_u_con)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_v_con, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dvcon",itau_w,d_v_con)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_q_con, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dqcon",itau_w,d_q_con)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_t_lsc, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dtlsc",itau_w,d_t_lsc)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_q_lsc, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dqlsc",itau_w,d_q_lsc)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_t_vdf, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dtvdf",itau_w,d_t_vdf)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_q_vdf, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dqvdf",itau_w,d_q_vdf)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_t_ajs, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dtajs",itau_w,d_t_ajs)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_q_ajs, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dqajs",itau_w,d_q_ajs)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_t_eva, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dteva",itau_w,d_t_eva)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_q_eva, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dqeva",itau_w,d_q_eva)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, heat, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dtswr",itau_w,heat)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, heat0, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dtsw0",itau_w,heat0)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, cool, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dtlwr",itau_w,cool)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, cool0, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dtlw0",itau_w,cool0)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_u_vdf, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"duvdf",itau_w,d_u_vdf)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_v_vdf, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dvvdf",itau_w,d_v_vdf)
c
      IF (ok_orodr) THEN
      IF (ok_orolf) THEN
c
      DO k = 1, klev
      DO i = 1, klon
        d_u_oli(i,k) = d_u_oro(i,k) + d_u_lif(i,k)
        d_v_oli(i,k) = d_v_oro(i,k) + d_v_lif(i,k)
      ENDDO
      ENDDO
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_u_oli, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"duoli",d_u_oli)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_v_oli, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dvoli",itau_w,d_v_oli)
c
      ENDIF
      ENDIF
C
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_u, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"duphy",itau_w,d_u)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_v, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dvphy",itau_w,d_v)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_t, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dtphy",itau_w,d_t)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_qx(:,:,1), 
cymf     .zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dqphy",itau_w,d_qx(:,:,1))
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_qx(:,:,2), 
cym     .zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPave,"dqlphy",itau_w,d_qx(:,:,2))
c
C
      if (ok_sync) then
        call histsync(nid_bilKPave)
      endif
       ENDIF

