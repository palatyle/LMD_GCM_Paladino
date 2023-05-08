!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/write_histmth.h,v 1.2 2005/03/09 12:30:16 fairhead Exp $
!
      IF (ok_mensuel) THEN

         itau_w = itau_phy + itap

c-------------------------------------------------------
      IF(lev_histmth.GE.1) THEN

ccccccccccccc 2D fields, basics

      call histwrite_phy(nid_mth,.false.,"phis",itau_w,pphis)
c      call histwrite_phy(nid_mth,.false.,"aire",itau_w,cell_area)
      cell_area_out(:)=cell_area(:)
      if (is_north_pole_phy) cell_area_out(1)=cell_area(1)/nbp_lon
      if (is_south_pole_phy) cell_area_out(klon)=cell_area(klon)/nbp_lon
      call histwrite_phy(nid_mth,.false.,"aire",itau_w,cell_area_out)

      call histwrite_phy(nid_mth,.false.,"tsol",itau_w,ftsol)
      call histwrite_phy(nid_mth,.false.,"psol",itau_w,paprs(:,1))
c     call histwrite_phy(nid_mth,.false.,"ue",itau_w,ue)
c VENUS: regardee a l envers!!!!!!!!!!!!!!!
c     call histwrite_phy(nid_mth,.false.,"ve",itau_w,-1.*ve)
c     call histwrite_phy(nid_mth,.false.,"cdragh",itau_w,cdragh)
c     call histwrite_phy(nid_mth,.false.,"cdragm",itau_w,cdragm)

      ENDIF !lev_histmth.GE.1

c-------------------------------------------------------
      IF(lev_histmth.GE.2) THEN

ccccccccccccc 3D fields, basics

      call histwrite_phy(nid_mth,.false.,"temp",itau_w,t_seri)
      call histwrite_phy(nid_mth,.false.,"pres",itau_w,pplay)
      call histwrite_phy(nid_mth,.false.,"geop",itau_w,zphi)
      call histwrite_phy(nid_mth,.false.,"vitu",itau_w,u_seri)
c VENUS: regardee a l envers!!!!!!!!!!!!!!!
      call histwrite_phy(nid_mth,.false.,"vitv",itau_w,-1.*v_seri)
      call histwrite_phy(nid_mth,.false.,"vitw",itau_w,omega)
c en (m/s)/s      
      call histwrite_phy(nid_mth,.false.,"dudyn",itau_w,d_u_dyn)
c en (m/s)/s      
      call histwrite_phy(nid_mth,.false.,"duvdf",itau_w,d_u_vdf)
c     call histwrite_phy(nid_mth,.false.,"mang",itau_w,mang)
c     call histwrite_phy(nid_mth,.false.,"Kz",itau_w,ycoefh)
      call histwrite_phy(nid_mth,.false.,"mmean",itau_w,mmean)
c     call histwrite_phy(nid_mth,.false.,"rho",itau_w,rho)

c plusieurs traceurs  !!!outputs in [vmr]
       IF (iflag_trac.eq.1) THEN
         DO iq=1,nqmax
       call histwrite_phy(nid_mth,.false.,tname(iq),itau_w,qx(:,:,iq)
     &			 *mmean(:,:)/M_tr(iq))
         ENDDO
       ENDIF

       IF (callthermos .and. ok_chem) THEN
       call histwrite_phy(nid_mth,.false.,"d_qmoldifCO2",itau_w,
     .                 d_q_moldif(:,:,i_co2))
       call histwrite_phy(nid_mth,.false.,"d_qmoldifO3p",itau_w,
     .                  d_q_moldif(:,:,i_o))
       call histwrite_phy(nid_mth,.false.,"d_qmoldifN2",itau_w,
     .                  d_q_moldif(:,:,i_n2))
       ENDIF

      call histwrite_phy(nid_mth,.false.,"tops",itau_w,topsw)

      ENDIF !lev_histmth.GE.2

c-------------------------------------------------------
      IF(lev_histmth.GE.3) THEN

cccccccccccccccccc  Radiative transfer

c 2D

      call histwrite_phy(nid_mth,.false.,"topl",itau_w,toplw)
      call histwrite_phy(nid_mth,.false.,"sols",itau_w,solsw)
      call histwrite_phy(nid_mth,.false.,"soll",itau_w,sollw)

c 3D

      call histwrite_phy(nid_mth,.false.,"SWnet",itau_w,swnet)
      call histwrite_phy(nid_mth,.false.,"LWnet",itau_w,lwnet)
c     call histwrite_phy(nid_mth,.false.,"fluxvdf",itau_w,fluxt)
c     call histwrite_phy(nid_mth,.false.,"fluxdyn",itau_w,flux_dyn)
c     call histwrite_phy(nid_mth,.false.,"fluxajs",itau_w,flux_ajs)
c     call histwrite_phy(nid_mth,.false.,"fluxec",itau_w,flux_ec)

      ENDIF !lev_histmth.GE.3

c-------------------------------------------------------
      IF(lev_histmth.GE.4) THEN

c en K/s      
      call histwrite_phy(nid_mth,.false.,"dtdyn",itau_w,d_t_dyn)
c en K/s      
c     call histwrite_phy(nid_mth,.false.,"dtphy",itau_w,d_t)
c en K/s      
      call histwrite_phy(nid_mth,.false.,"dtvdf",itau_w,d_t_vdf)
c en K/s      
      call histwrite_phy(nid_mth,.false.,"dtajs",itau_w,d_t_ajs)
c en K/s      
      call histwrite_phy(nid_mth,.false.,"dtswr",itau_w,dtsw)
      call histwrite_phy(nid_mth,.false.,"dtswrNLTE",itau_w,d_t_nirco2)
      call histwrite_phy(nid_mth,.false.,"dtswrDCrisp",itau_w,heat)
c en K/s      
      call histwrite_phy(nid_mth,.false.,"dtlwr",itau_w,dtlw)
c en K/s      
c     call histwrite_phy(nid_mth,.false.,"dtlwrNLTE",itau_w,d_t_nlte)
c     call histwrite_phy(nid_mth,.false.,"dtlwrLTE",itau_w,-1.*cool)
c     call histwrite_phy(nid_mth,.false.,"dteuv",itau_w,d_t_euv)
c     call histwrite_phy(nid_mth,.false.,"dtcond",itau_w,d_t_conduc)
c     call histwrite_phy(nid_mth,.false.,"dumolvis",itau_w,d_u_molvis)
c     call histwrite_phy(nid_mth,.false.,"dvmolvis",itau_w,-1.*d_v_molvis)

c     call histwrite_phy(nid_mth,.false.,"dtec",itau_w,d_t_ec)
c en (m/s)/s      
      call histwrite_phy(nid_mth,.false.,"duajs",itau_w,d_u_ajs)
c en (m/s)/s      
      call histwrite_phy(nid_mth,.false.,"dugwo",itau_w,d_u_oro)
c en (m/s)/s      
      call histwrite_phy(nid_mth,.false.,"dugwno",itau_w,d_u_hin)
c en (m/s)/s      

    
c     VENUS: regardee a l envers!!!!!!!!!!!!!!!
c     call histwrite_phy(nid_mth,.false.,"dvvdf",itau_w,-1.*d_v_vdf)

      ENDIF !lev_histmth.GE.4

c-------------------------------------------------------
      IF(lev_histmth.GE.5) THEN

c     call histwrite_phy(nid_mth,.false.,"taux_",itau_w,fluxu)
c     call histwrite_phy(nid_mth,.false.,"tauy_",itau_w,fluxv)
c     call histwrite_phy(nid_mth,.false.,"cdrm",itau_w,cdragm)
c     call histwrite_phy(nid_mth,.false.,"cdrh",itau_w,cdragh)

      ENDIF !lev_histmth.GE.5
c-------------------------------------------------------

      if (ok_sync) then
        call histsync(nid_mth)
      endif

      ENDIF

