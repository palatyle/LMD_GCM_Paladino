!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/write_histday.h,v 1.2 2004/06/01 09:27:10 lmdzadmin Exp $
!
      IF (ok_journe) THEN

         itau_w = itau_phy + itap

c-------------------------------------------------------
      IF(lev_histday.GE.1) THEN

ccccccccccccc 2D fields, basics

      call histwrite_phy(nid_day,.false.,"phis",itau_w,pphis)
c      call histwrite_phy(nid_day,.false.,"aire",itau_w,cell_area)
      cell_area_out(:)=cell_area(:)
      if (is_north_pole_phy) cell_area_out(1)=cell_area(1)/nbp_lon
      if (is_south_pole_phy) cell_area_out(klon)=cell_area(klon)/nbp_lon
      call histwrite_phy(nid_day,.false.,"aire",itau_w,cell_area_out)

      call histwrite_phy(nid_day,.false.,"tsol",itau_w,ftsol)
      call histwrite_phy(nid_day,.false.,"psol",itau_w,paprs(:,1))
c     call histwrite_phy(nid_day,.false.,"ue",itau_w,ue)
c VENUS: regardee a l'envers!!!!!!!!!!!!!!!
c     call histwrite_phy(nid_day,.false.,"ve",itau_w,-1.*ve)
c     call histwrite_phy(nid_day,.false.,"cdragh",itau_w,cdragh)
c     call histwrite_phy(nid_day,.false.,"cdragm",itau_w,cdragm)

      ENDIF !lev_histday.GE.1

c-------------------------------------------------------
      IF(lev_histday.GE.2) THEN

ccccccccccccc 3D fields, basics

      call histwrite_phy(nid_day,.false.,"temp",itau_w,t_seri)
      call histwrite_phy(nid_day,.false.,"pres",itau_w,pplay)
      call histwrite_phy(nid_day,.false.,"geop",itau_w,zphi)
      call histwrite_phy(nid_day,.false.,"vitu",itau_w,u_seri)
c VENUS: regardee a l'envers!!!!!!!!!!!!!!!
      call histwrite_phy(nid_day,.false.,"vitv",itau_w,-1.*v_seri)
      call histwrite_phy(nid_day,.false.,"vitw",itau_w,omega)
c en (m/s)/s      
      call histwrite_phy(nid_day,.false.,"dudyn",itau_w,d_u_dyn)
c en (m/s)/s      
      call histwrite_phy(nid_day,.false.,"duvdf",itau_w,d_u_vdf)
c     call histwrite_phy(nid_day,.false.,"mang",itau_w,mang)
c     call histwrite_phy(nid_day,.false.,"Kz",itau_w,ycoefh)

c plusieurs traceurs
       IF (iflag_trac.eq.1) THEN
         DO iq=1,nqmax
       call histwrite_phy(nid_day,.false.,tname(iq),itau_w,qx(:,:,iq))
         ENDDO
       ENDIF

      call histwrite_phy(nid_day,.false.,"tops",itau_w,topsw)

      ENDIF !lev_histday.GE.2

c-------------------------------------------------------
      IF(lev_histday.GE.3) THEN

cccccccccccccccccc  Radiative transfer

c 2D

      call histwrite_phy(nid_day,.false.,"topl",itau_w,toplw)
      call histwrite_phy(nid_day,.false.,"sols",itau_w,solsw)
      call histwrite_phy(nid_day,.false.,"soll",itau_w,sollw)

c 3D

      call histwrite_phy(nid_day,.false.,"SWnet",itau_w,swnet)
      call histwrite_phy(nid_day,.false.,"LWnet",itau_w,lwnet)
c     call histwrite_phy(nid_day,.false.,"fluxvdf",itau_w,fluxt)
c     call histwrite_phy(nid_day,.false.,"fluxdyn",itau_w,flux_dyn)
c     call histwrite_phy(nid_day,.false.,"fluxajs",itau_w,flux_ajs)
c     call histwrite_phy(nid_day,.false.,"fluxec",itau_w,flux_ec)

      ENDIF !lev_histday.GE.3

c-------------------------------------------------------
      IF(lev_histday.GE.4) THEN

c en K/s      
      call histwrite_phy(nid_day,.false.,"dtdyn",itau_w,d_t_dyn)
c en K/s      
c     call histwrite_phy(nid_day,.false.,"dtphy",itau_w,d_t)
c en K/s      
      call histwrite_phy(nid_day,.false.,"dtvdf",itau_w,d_t_vdf)
c en K/s      
      call histwrite_phy(nid_day,.false.,"dtajs",itau_w,d_t_ajs)
c en K/s      
      call histwrite_phy(nid_day,.false.,"dtswr",itau_w,dtsw)
c en K/s      
      call histwrite_phy(nid_day,.false.,"dtlwr",itau_w,dtlw)
c en K/s      
c     call histwrite_phy(nid_day,.false.,"dtec",itau_w,d_t_ec)
c en (m/s)/s      
      call histwrite_phy(nid_day,.false.,"duajs",itau_w,d_u_ajs)
c en (m/s)/s      
      call histwrite_phy(nid_day,.false.,"dugwo",itau_w,d_u_oro)
c en (m/s)/s      
      call histwrite_phy(nid_day,.false.,"dugwno",itau_w,d_u_hin)
c en (m/s)/s      
c VENUS: regardee a l'envers!!!!!!!!!!!!!!!
c     call histwrite_phy(nid_day,.false.,"dvvdf",itau_w,-1.*d_v_vdf)

      ENDIF !lev_histday.GE.4

c-------------------------------------------------------
      IF(lev_histday.GE.5) THEN

c     call histwrite_phy(nid_day,.false.,"taux_",itau_w,fluxu)
c     call histwrite_phy(nid_day,.false.,"tauy_",itau_w,fluxv)
c     call histwrite_phy(nid_day,.false.,"cdrm",itau_w,cdragm)
c     call histwrite_phy(nid_day,.false.,"cdrh",itau_w,cdragh)

      ENDIF !lev_histday.GE.5
c-------------------------------------------------------

      if (ok_sync) then
        call histsync(nid_day)
      endif

      ENDIF
