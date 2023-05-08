 c
c $Header$
c
      IF (ok_journe) THEN
c
      ndex2d = 0
      ndex3d = 0
c
      itau_w = itau_phy + itap
c
c Champs 3D:
c
cym      CALL gr_fi_ecrit(klev, klon,iim,jjmp1, ue_lay,zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"ue",itau_w,ue_lay)
c
cym      CALL gr_fi_ecrit(klev, klon,iim,jjmp1, ve_lay,zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"ve",itau_w,ve_lay)
c
cym      CALL gr_fi_ecrit(klev, klon,iim,jjmp1, uq_lay,zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"uq",itau_w,uq_lay)
c
cym      CALL gr_fi_ecrit(klev, klon,iim,jjmp1, vq_lay,zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"vq",itau_w,vq_lay)
c
c Champs 3D:
C
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, t_seri, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"temp",itau_w,t_seri)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, qx(1,1,ivap), zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"ovap",itau_w,qx(:,:,ivap))
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, zphi, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"geop",itau_w,zphi)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, u_seri, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"vitu",itau_w,u_seri)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, v_seri, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"vitv",itau_w,v_seri)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, omega, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"vitw",itau_w,omega)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, pplay, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"pres",itau_w,pplay)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, paprs, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"play",itau_w,paprs)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, cldliq, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"oliq",itau_w,cldliq)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_t_dyn, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dtdyn",itau_w,d_t_dyn)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_q_dyn, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dqdyn",itau_w,d_q_dyn)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_t_con, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dtcon",itau_w,d_t_con)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_u_con, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"ducon",itau_w,d_u_con)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_v_con, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dvcon",itau_w,d_v_con)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_q_con, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dqcon",itau_w,d_q_con)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_t_lsc, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dtlsc",itau_w,d_t_lsc)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_q_lsc, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dqlsc",itau_w,d_q_lsc)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_t_vdf, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dtvdf",itau_w,d_t_vdf)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_q_vdf, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dqvdf",itau_w,d_q_vdf)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_t_ajs, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dtajs",itau_w,d_t_ajs)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_q_ajs, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dqajs",itau_w,d_q_ajs)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_t_eva, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dteva",itau_w,d_t_eva)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_q_eva, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dqeva",itau_w,d_q_eva)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, heat, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dtswr",itau_w,heat)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, heat0, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dtsw0",itau_w,heat0)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, cool, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dtlwr",itau_w,cool)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, cool0, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dtlw0",itau_w,cool0)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_u_vdf, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"duvdf",itau_w,d_u_vdf)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_v_vdf, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dvvdf",itau_w,d_v_vdf)
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
      CALL histwrite_phy(nid_bilKPins,"duoli",itau_w,d_u_oli)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_v_oli, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dvoli",itau_w,d_v_oli)
c
      ENDIF
      ENDIF
C
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_u, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"duphy",itau_w,d_u)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_v, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dvphy",itau_w,d_v)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_t, zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dtphy",itau_w,d_t)
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_qx(:,:,1), 
cym     .zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dqphy",itau_w,d_qx(:,:,1))
c
cym      CALL gr_fi_ecrit(klev,klon,iim,jjmp1, d_qx(:,:,2), 
cym     .zx_tmp_3d)
      CALL histwrite_phy(nid_bilKPins,"dqlphy",itau_w,d_qx(:,:,2))
c
cIM 280405 BEG
c
c Champs 2D:
c
c   Ecriture de champs dynamiques sur des niveaux de pression
c     DO k=1, nlevSTD
      DO k=1, 12
c
       IF(k.GE.2.AND.k.LE.12) bb2=clevSTD(k)
       IF(k.GE.13.AND.k.LE.17) bb3=clevSTD(k)
c
       IF(bb2.EQ."850") THEN
c
cym        CALL gr_fi_ecrit(1, klon,iim,jjmp1,usumSTD(:,k,1),zx_tmp_2d)
        CALL histwrite_phy(nid_bilKPins,"u"//bb2,itau_w,usumSTD(:,k,1))
c
cym        CALL gr_fi_ecrit(1, klon,iim,jjmp1,vsumSTD(:,k,1),zx_tmp_2d)
        CALL histwrite_phy(nid_bilKPins,"v"//bb2,itau_w,vsumSTD(:,k,1))
c
       ENDIF !(bb2.EQ."850")
c
       ENDDO !k=1, 12
c
cIM 280405 END
C
      if (ok_sync) then
        call histsync(nid_bilKPins)
      endif
       ENDIF

