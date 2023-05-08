!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/write_histins.h,v 1.1.1.1 2004/05/19 12:53:09 lmdzadmin Exp $
!
      IF (ok_instan) THEN

         itau_w = itau_phy + itap

c-------------------------------------------------------
      IF(lev_histins.GE.1) THEN

ccccccccccccc 2D fields, basics

      call histwrite_phy(nid_ins,.false.,"phis",itau_w,pphis)
c      call histwrite_phy(nid_ins,.false.,"aire",itau_w,cell_area)
      cell_area_out(:)=cell_area(:)
      if (is_north_pole_phy) cell_area_out(1)=cell_area(1)/nbp_lon
      if (is_south_pole_phy) cell_area_out(klon)=cell_area(klon)/nbp_lon
      call histwrite_phy(nid_ins,.false.,"aire",itau_w,cell_area_out)

      call histwrite_phy(nid_ins,.false.,"tsol",itau_w,ftsol)
      call histwrite_phy(nid_ins,.false.,"psol",itau_w,paprs(:,1))
c     call histwrite_phy(nid_ins,.false.,"ue",itau_w,ue)
c VENUS: regardee a l' envers!!!!!!!!!!!!!!!
c     call histwrite_phy(nid_ins,.false.,"ve",itau_w,-1.*ve)
c     call histwrite_phy(nid_ins,.false.,"cdragh",itau_w,cdragh)
c     call histwrite_phy(nid_ins,.false.,"cdragm",itau_w,cdragm)

      ENDIF !lev_histins.GE.1

c-------------------------------------------------------
      IF(lev_histins.GE.2) THEN

ccccccccccccc 3D fields, basics

      call histwrite_phy(nid_ins,.false.,"temp",itau_w,t_seri)
      call histwrite_phy(nid_ins,.false.,"pres",itau_w,pplay)
      call histwrite_phy(nid_ins,.false.,"geop",itau_w,zphi)
      call histwrite_phy(nid_ins,.false.,"vitu",itau_w,u_seri)
c VENUS: regardee a l' envers !!!!!!!!!!!!!!!
      call histwrite_phy(nid_ins,.false.,"vitv",itau_w,-1.*v_seri)
      call histwrite_phy(nid_ins,.false.,"vitw",itau_w,omega)
c en (m/s)/s      
      call histwrite_phy(nid_ins,.false.,"dudyn",itau_w,d_u_dyn)
c en (m/s)/s      
      call histwrite_phy(nid_ins,.false.,"duvdf",itau_w,d_u_vdf)
c     call histwrite_phy(nid_ins,.false.,"mang",itau_w,mang)
c     call histwrite_phy(nid_ins,.false.,"Kz",itau_w,ycoefh)
      call histwrite_phy(nid_ins,.false.,"mmean",itau_w,mmean)
      call histwrite_phy(nid_ins,.false.,"rho",itau_w,rho)

c plusieurs traceurs  !!!outputs in [vmr]
       IF (iflag_trac.eq.1) THEN
         DO iq=1,nqmax
       call histwrite_phy(nid_ins,.false.,tname(iq),itau_w,qx(:,:,iq)*
     &			  mmean(:,:)/M_tr(iq))
         ENDDO
       ENDIF

       IF (callthermos .and. ok_chem) THEN
       call histwrite_phy(nid_ins,.false.,"d_qmoldifCO2",itau_w,
     .                 d_q_moldif(:,:,i_co2))
       call histwrite_phy(nid_ins,.false.,"d_qmoldifO3p",itau_w,
     .                  d_q_moldif(:,:,i_o))
       call histwrite_phy(nid_ins,.false.,"d_qmoldifN2",itau_w,
     .                  d_q_moldif(:,:,i_n2))
       ENDIF
      
       call histwrite_phy(nid_ins,.false.,"tops",itau_w,topsw)
      
       IF (ok_cloud.and.(cl_scheme.eq.1)) THEN
       
       IF (nb_mode.GE.1) THEN
      call histwrite_phy(nid_ins,.false.,"NBRTOTm1",itau_w,
     & NBRTOT(:,:,1))
      
c      call histwrite_phy(nid_ins,.false.,"R_MEDIANm1",itau_w,
c     & R_MEDIAN(:,:,1))
      
c      call histwrite_phy(nid_ins,.false.,"STDDEVm1",itau_w,
c     & STDDEV(:,:,1))
       
       IF (nb_mode.GE.2) THEN
      call histwrite_phy(nid_ins,.false.,"NBRTOTm2",itau_w,
     & NBRTOT(:,:,2))
      
c      call histwrite_phy(nid_ins,.false.,"R_MEDIANm2",itau_w,
c     & R_MEDIAN(:,:,2)) 
      
c      call histwrite_phy(nid_ins,.false.,"STDDEVm2",itau_w,
c     & STDDEV(:,:,2))
          
       IF (nb_mode.GE.3) THEN
      call histwrite_phy(nid_ins,.false.,"NBRTOTm3",itau_w,
     & NBRTOT(:,:,3))
      
c      call histwrite_phy(nid_ins,.false.,"R_MEDIANm3",itau_w,
c     & R_MEDIAN(:,:,3))
     
c      call histwrite_phy(nid_ins,.false.,"STDDEVm3",itau_w,
c     & STDDEV(:,:,3))

       ENDIF
       ENDIF
       ENDIF
             
      call histwrite_phy(nid_ins,.false.,"WH2SO4",itau_w,WH2SO4)
      
      call histwrite_phy(nid_ins,.false.,"rho_droplet",itau_w,
     & rho_droplet)
		ENDIF

       IF (ok_sedim.and.(cl_scheme.eq.1)) THEN 
      
      call histwrite_phy(nid_ins,.false.,"d_tr_sed_H2SO4",itau_w,
     & d_tr_sed(:,:,1))      
      call histwrite_phy(nid_ins,.false.,"d_tr_sed_H2O",itau_w,
     & d_tr_sed(:,:,2))
      
      call histwrite_phy(nid_ins,.false.,"F_sedim",itau_w,Fsedim)
		ENDIF             

      ENDIF !lev_histins.GE.2

c-------------------------------------------------------
      IF(lev_histins.GE.3) THEN

cccccccccccccccccc  Radiative transfer

c 2D

      call histwrite_phy(nid_ins,.false.,"topl",itau_w,toplw)
      call histwrite_phy(nid_ins,.false.,"sols",itau_w,solsw)
      call histwrite_phy(nid_ins,.false.,"soll",itau_w,sollw)

c 3D

      call histwrite_phy(nid_ins,.false.,"SWnet",itau_w,swnet)
      call histwrite_phy(nid_ins,.false.,"LWnet",itau_w,lwnet)
c     call histwrite_phy(nid_ins,.false.,"fluxvdf",itau_w,fluxt)
c     call histwrite_phy(nid_ins,.false.,"fluxdyn",itau_w,flux_dyn)
c     call histwrite_phy(nid_ins,.false.,"fluxajs",itau_w,flux_ajs)
c     call histwrite_phy(nid_ins,.false.,"fluxec",itau_w,flux_ec)

      ENDIF !lev_histins.GE.3

c-------------------------------------------------------
      IF(lev_histins.GE.4) THEN

c en K/s      
      call histwrite_phy(nid_ins,.false.,"dtdyn",itau_w,d_t_dyn)
c en K/s      
c     call histwrite_phy(nid_ins,.false.,"dtphy",itau_w,d_t)
c en K/s      
      call histwrite_phy(nid_ins,.false.,"dtvdf",itau_w,d_t_vdf)
c en K/s      
      call histwrite_phy(nid_ins,.false.,"dtajs",itau_w,d_t_ajs)
c K/day ==> K/s
      call histwrite_phy(nid_ins,.false.,"dtswr",itau_w,heat/RDAY)
c K/day ==> K/s      
      call histwrite_phy(nid_ins,.false.,"dtlwr",itau_w,-1.*cool/RDAY)
c en K/s      
c     call histwrite_phy(nid_ins,.false.,"dtec",itau_w,d_t_ec)
c en (m/s)/s      
      call histwrite_phy(nid_ins,.false.,"duajs",itau_w,d_u_ajs)
c en (m/s)/s      
      call histwrite_phy(nid_ins,.false.,"dugwo",itau_w,d_u_oro)
c en (m/s)/s      
      call histwrite_phy(nid_ins,.false.,"dugwno",itau_w,d_u_hin)
c en (m/s)/s      
c VENUS: regardee a l envers!!!!!!!!!!!!!!!
c     call histwrite_phy(nid_ins,.false.,"dvvdf",itau_w,-1.*d_v_vdf)

      ENDIF !lev_histins.GE.4

c-------------------------------------------------------
      IF(lev_histins.GE.5) THEN

c     call histwrite_phy(nid_ins,.false.,"taux_",itau_w,fluxu)
c     call histwrite_phy(nid_ins,.false.,"tauy_",itau_w,fluxv)
c     call histwrite_phy(nid_ins,.false.,"cdrm",itau_w,cdragm)
c     call histwrite_phy(nid_ins,.false.,"cdrh",itau_w,cdragh)

      ENDIF !lev_histins.GE.5
c-------------------------------------------------------

      if (ok_sync) then
        call histsync(nid_ins)
      endif

      ENDIF
