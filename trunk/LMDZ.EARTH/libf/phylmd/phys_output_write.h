      itau_w = itau_phy + itap

      DO iff=1,nfiles

       IF (clef_files(iff)) THEN
             ndex2d = 0
             ndex3d = 0

!!! Champs 1D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF (o_phis%flag(iff)<=lev_files(iff)) THEN 
         CALL histwrite_phy(nid_files(iff),
     $                      o_phis%name,itau_w,pphis)
       ENDIF

       IF (o_aire%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),o_aire%name,itau_w,airephy)
       ENDIF

       IF (o_contfracATM%flag(iff)<=lev_files(iff)) THEN
      DO i=1, klon
       zx_tmp_fi2d(i)=pctsrf(i,is_ter)+pctsrf(i,is_lic)
      ENDDO
      CALL histwrite_phy(nid_files(iff),
     $             o_contfracATM%name,itau_w,zx_tmp_fi2d)
       ENDIF

       IF (o_contfracOR%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_contfracOR%name,itau_w,
     $                   pctsrf(:,is_ter))
       ENDIF

       IF (o_aireTER%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     $                  o_aireTER%name,itau_w,paire_ter)
       ENDIF

!!! Champs 2D !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       IF (o_flat%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_flat%name,itau_w,zxfluxlat)
       ENDIF

       IF (o_slp%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_slp%name,itau_w,slp)
       ENDIF

       IF (o_tsol%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_tsol%name,itau_w,zxtsol)
       ENDIF

       IF (o_t2m%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_t2m%name,itau_w,zt2m)
       ENDIF

       IF (o_t2m_min%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_t2m_min%name,itau_w,zt2m)
       ENDIF

       IF (o_t2m_max%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_t2m_max%name,itau_w,zt2m)
       ENDIF

       IF (o_wind10m%flag(iff)<=lev_files(iff)) THEN
      DO i=1, klon
       zx_tmp_fi2d(i)=SQRT(zu10m(i)*zu10m(i)+zv10m(i)*zv10m(i))
      ENDDO
      CALL histwrite_phy(nid_files(iff),
     s                  o_wind10m%name,itau_w,zx_tmp_fi2d)
       ENDIF

       IF (o_wind10max%flag(iff)<=lev_files(iff)) THEN
      DO i=1, klon
       zx_tmp_fi2d(i)=SQRT(zu10m(i)*zu10m(i)+zv10m(i)*zv10m(i))
      ENDDO
      CALL histwrite_phy(nid_files(iff),o_wind10max%name, 
     $                   itau_w,zx_tmp_fi2d)
       ENDIF

       IF (o_sicf%flag(iff)<=lev_files(iff)) THEN
      DO i = 1, klon
         zx_tmp_fi2d(i) = pctsrf(i,is_sic)
      ENDDO
      CALL histwrite_phy(nid_files(iff),
     $                   o_sicf%name,itau_w,zx_tmp_fi2d)
       ENDIF

       IF (o_q2m%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_q2m%name,itau_w,zq2m)
       ENDIF

       IF (o_u10m%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_u10m%name,itau_w,zu10m)
       ENDIF

       IF (o_v10m%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_v10m%name,itau_w,zv10m)
       ENDIF

       IF (o_psol%flag(iff)<=lev_files(iff)) THEN
      DO i = 1, klon
         zx_tmp_fi2d(i) = paprs(i,1)
      ENDDO
      CALL histwrite_phy(nid_files(iff),
     s                   o_psol%name,itau_w,zx_tmp_fi2d)
       ENDIF

       IF (o_mass%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_mass%name,itau_w,zmasse)
        ENDIF


       IF (o_qsurf%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_qsurf%name,itau_w,zxqsurf)
       ENDIF

       if (.not. ok_veget) then
         IF (o_qsol%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_qsol%name,itau_w,qsol)
         ENDIF
       endif

      IF (o_precip%flag(iff)<=lev_files(iff)) THEN
       DO i = 1, klon
         zx_tmp_fi2d(i) = rain_fall(i) + snow_fall(i)
       ENDDO
      CALL histwrite_phy(nid_files(iff),o_precip%name,
     s                   itau_w,zx_tmp_fi2d)
      ENDIF

       IF (o_ndayrain%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_ndayrain%name,
     s                   itau_w,nday_rain)
       ENDIF

      IF (o_plul%flag(iff)<=lev_files(iff)) THEN
       DO i = 1, klon
         zx_tmp_fi2d(i) = rain_lsc(i) + snow_lsc(i)
       ENDDO
      CALL histwrite_phy(nid_files(iff),o_plul%name,itau_w,zx_tmp_fi2d)
      ENDIF

      IF (o_pluc%flag(iff)<=lev_files(iff)) THEN
      DO i = 1, klon
         zx_tmp_fi2d(i) = rain_con(i) + snow_con(i)
      ENDDO
      CALL histwrite_phy(nid_files(iff),o_pluc%name,itau_w,zx_tmp_fi2d)
      ENDIF

       IF (o_snow%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_snow%name,itau_w,snow_fall)
       ENDIF

       IF (o_msnow%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_msnow%name,itau_w,snow_o)
       ENDIF

       IF (o_fsnow%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_fsnow%name,itau_w,zfra_o)
       ENDIF

       IF (o_evap%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_evap%name,itau_w,evap)
       ENDIF

       IF (o_tops%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_tops%name,itau_w,topsw)
       ENDIF

       IF (o_tops0%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_tops0%name,itau_w,topsw0)
       ENDIF

       IF (o_topl%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_topl%name,itau_w,toplw)
       ENDIF

       IF (o_topl0%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_topl0%name,itau_w,toplw0)
       ENDIF

       IF (o_SWupTOA%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(1 : klon) = swup ( 1 : klon, klevp1 )
      CALL histwrite_phy(nid_files(iff),o_SWupTOA%name,
     s                     itau_w,zx_tmp_fi2d)
       ENDIF

       IF (o_SWupTOAclr%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(1 : klon) = swup0 ( 1 : klon, klevp1 )
      CALL histwrite_phy(nid_files(iff), 
     $                  o_SWupTOAclr%name,itau_w,zx_tmp_fi2d)
       ENDIF

       IF (o_SWdnTOA%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(1 : klon) = swdn ( 1 : klon, klevp1 )
      CALL histwrite_phy(nid_files(iff),
     s                  o_SWdnTOA%name,itau_w,zx_tmp_fi2d)
       ENDIF

       IF (o_SWdnTOAclr%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(1 : klon) = swdn0 ( 1 : klon, klevp1 )
      CALL histwrite_phy(nid_files(iff), 
     $                  o_SWdnTOAclr%name,itau_w,zx_tmp_fi2d)
       ENDIF

       IF (o_nettop%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(:) = topsw(:)-toplw(:)
      CALL histwrite_phy(nid_files(iff),
     $                  o_nettop%name,itau_w,zx_tmp_fi2d)
       ENDIF

       IF (o_SWup200%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_SWup200%name,itau_w,SWup200)
       ENDIF

       IF (o_SWup200clr%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     s                   o_SWup200clr%name,itau_w,SWup200clr)
       ENDIF

       IF (o_SWdn200%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_SWdn200%name,itau_w,SWdn200)
       ENDIF

       IF (o_SWdn200clr%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     s                o_SWdn200clr%name,itau_w,SWdn200clr)
       ENDIF

       IF (o_LWup200%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_LWup200%name,itau_w,LWup200)
       ENDIF

       IF (o_LWup200clr%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     s                   o_LWup200clr%name,itau_w,LWup200clr)
       ENDIF

       IF (o_LWdn200%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     s                   o_LWdn200%name,itau_w,LWdn200)
       ENDIF

       IF (o_LWdn200clr%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     s                  o_LWdn200clr%name,itau_w,LWdn200clr)
       ENDIF

       IF (o_sols%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_sols%name,itau_w,solsw)
       ENDIF

       IF (o_sols0%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_sols0%name,itau_w,solsw0)
       ENDIF

       IF (o_soll%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_soll%name,itau_w,sollw)
       ENDIF

       IF (o_radsol%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_radsol%name,itau_w,radsol)
       ENDIF

       IF (o_soll0%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_soll0%name,itau_w,sollw0)
       ENDIF

       IF (o_SWupSFC%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(1 : klon) = swup ( 1 : klon, 1 )
      CALL histwrite_phy(nid_files(iff),
     s               o_SWupSFC%name,itau_w,zx_tmp_fi2d)
       ENDIF

       IF (o_SWupSFCclr%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(1 : klon) = swup0 ( 1 : klon, 1 )
      CALL histwrite_phy(nid_files(iff), 
     $                   o_SWupSFCclr%name,itau_w,zx_tmp_fi2d)
       ENDIF

       IF (o_SWdnSFC%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(1 : klon) = swdn ( 1 : klon, 1 )
      CALL histwrite_phy(nid_files(iff), 
     $                   o_SWdnSFC%name,itau_w,zx_tmp_fi2d)
       ENDIF

       IF (o_SWdnSFCclr%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(1 : klon) = swdn0 ( 1 : klon, 1 )
      CALL histwrite_phy(nid_files(iff), 
     $                  o_SWdnSFCclr%name,itau_w,zx_tmp_fi2d)
       ENDIF

       IF (o_LWupSFC%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(1:klon)=sollwdown(1:klon)-sollw(1:klon)
      CALL histwrite_phy(nid_files(iff),
     $                    o_LWupSFC%name,itau_w,zx_tmp_fi2d)
       ENDIF

       IF (o_LWdnSFC%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     $                   o_LWdnSFC%name,itau_w,sollwdown)
       ENDIF

       sollwdownclr(1:klon) = -1.*lwdn0(1:klon,1)
       IF (o_LWupSFCclr%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(1:klon)=sollwdownclr(1:klon)-sollw0(1:klon)
      CALL histwrite_phy(nid_files(iff),
     $                   o_LWupSFCclr%name,itau_w,zx_tmp_fi2d)
       ENDIF

       IF (o_LWdnSFCclr%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     $                   o_LWdnSFCclr%name,itau_w,sollwdownclr)
       ENDIF

       IF (o_bils%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_bils%name,itau_w,bils)
       ENDIF

       IF (o_sens%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(1:klon)=-1*sens(1:klon)
      CALL histwrite_phy(nid_files(iff),o_sens%name,itau_w,zx_tmp_fi2d)
       ENDIF

       IF (o_fder%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_fder%name,itau_w,fder)
       ENDIF

       IF (o_ffonte%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),o_ffonte%name,itau_w,zxffonte)
       ENDIF

       IF (o_fqcalving%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),
     $                    o_fqcalving%name,itau_w,zxfqcalving)
       ENDIF

       IF (o_fqfonte%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),
     $                   o_fqfonte%name,itau_w,zxfqfonte)
       ENDIF

       IF (o_taux%flag(iff)<=lev_files(iff)) THEN
         zx_tmp_fi2d=0.
         do nsrf=1,nbsrf
          zx_tmp_fi2d(:)=zx_tmp_fi2d(:)+pctsrf(:,nsrf)*fluxu(:,1,nsrf)
         enddo
         CALL histwrite_phy(nid_files(iff),
     $                   o_taux%name,itau_w,zx_tmp_fi2d)
       ENDIF

       IF (o_tauy%flag(iff)<=lev_files(iff)) THEN
         zx_tmp_fi2d=0.
         do nsrf=1,nbsrf
          zx_tmp_fi2d(:)=zx_tmp_fi2d(:)+pctsrf(:,nsrf)*fluxv(:,1,nsrf)
         enddo
         CALL histwrite_phy(nid_files(iff),
     $                   o_tauy%name,itau_w,zx_tmp_fi2d)
       ENDIF

         DO nsrf = 1, nbsrf
!           IF(nsrf.GE.2) THEN
            IF (o_pourc_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN
            zx_tmp_fi2d(1 : klon) = pctsrf( 1 : klon, nsrf)*100.
            CALL histwrite_phy(nid_files(iff),
     $                     o_pourc_srf(nsrf)%name,itau_w,
     $                     zx_tmp_fi2d)
            ENDIF

          IF (o_fract_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN
          zx_tmp_fi2d(1 : klon) = pctsrf( 1 : klon, nsrf)
          CALL histwrite_phy(nid_files(iff),
     $                  o_fract_srf(nsrf)%name,itau_w,
     $                  zx_tmp_fi2d)
          ENDIF
!         ENDIF !nsrf.GT.2

        IF (o_taux_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN
        zx_tmp_fi2d(1 : klon) = fluxu( 1 : klon, 1, nsrf)
        CALL histwrite_phy(nid_files(iff),
     $                     o_taux_srf(nsrf)%name,itau_w,
     $                     zx_tmp_fi2d)
        ENDIF

        IF (o_tauy_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN           
        zx_tmp_fi2d(1 : klon) = fluxv( 1 : klon, 1, nsrf)
        CALL histwrite_phy(nid_files(iff),
     $                    o_tauy_srf(nsrf)%name,itau_w,
     $                    zx_tmp_fi2d)
        ENDIF

        IF (o_tsol_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN
        zx_tmp_fi2d(1 : klon) = ftsol( 1 : klon, nsrf)
        CALL histwrite_phy(nid_files(iff),
     $                   o_tsol_srf(nsrf)%name,itau_w,
     $      zx_tmp_fi2d)
        ENDIF

      IF (o_u10m_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(1 : klon) = u10m(1 : klon, nsrf)
      CALL histwrite_phy(nid_files(iff),o_u10m_srf(nsrf)%name,
     $                 itau_w,zx_tmp_fi2d)
      ENDIF

      IF (o_v10m_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(1 : klon) = v10m(1 : klon, nsrf)
      CALL histwrite_phy(nid_files(iff),o_v10m_srf(nsrf)%name,
     $              itau_w,zx_tmp_fi2d)
      ENDIF
 
      IF (o_t2m_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(1 : klon) = t2m(1 : klon, nsrf)
      CALL histwrite_phy(nid_files(iff),o_t2m_srf(nsrf)%name,
     $           itau_w,zx_tmp_fi2d)
      ENDIF

      IF (o_evap_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(1 : klon) = fevap(1 : klon, nsrf)
      CALL histwrite_phy(nid_files(iff),o_evap_srf(nsrf)%name,
     $           itau_w,zx_tmp_fi2d)
      ENDIF

       IF (o_sens_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN
       zx_tmp_fi2d(1 : klon) = fluxt( 1 : klon, 1, nsrf)
       CALL histwrite_phy(nid_files(iff),
     $                    o_sens_srf(nsrf)%name,itau_w,
     $      zx_tmp_fi2d)
       ENDIF

        IF (o_lat_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN
        zx_tmp_fi2d(1 : klon) = fluxlat( 1 : klon, nsrf)
        CALL histwrite_phy(nid_files(iff),
     $                 o_lat_srf(nsrf)%name,itau_w,
     $                                   zx_tmp_fi2d)
          ENDIF

        IF (o_flw_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN
        zx_tmp_fi2d(1 : klon) = fsollw( 1 : klon, nsrf)
        CALL histwrite_phy(nid_files(iff),
     $                     o_flw_srf(nsrf)%name,itau_w,
     $      zx_tmp_fi2d)
        ENDIF

        IF (o_fsw_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN
        zx_tmp_fi2d(1 : klon) = fsolsw( 1 : klon, nsrf)
        CALL histwrite_phy(nid_files(iff),
     $                   o_fsw_srf(nsrf)%name,itau_w,
     $      zx_tmp_fi2d)
        ENDIF

        IF (o_wbils_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN
        zx_tmp_fi2d(1 : klon) = wfbils( 1 : klon, nsrf)
        CALL histwrite_phy(nid_files(iff),
     $                   o_wbils_srf(nsrf)%name,itau_w,
     $      zx_tmp_fi2d)
        ENDIF

        IF (o_wbilo_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN
        zx_tmp_fi2d(1 : klon) = wfbilo( 1 : klon, nsrf)
        CALL histwrite_phy(nid_files(iff),
     $                    o_wbilo_srf(nsrf)%name,itau_w,
     $      zx_tmp_fi2d)
        ENDIF

       if (iflag_pbl>1 .and. lev_histday.gt.10 ) then
        IF (o_tke_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),
     $                   o_tke_srf(nsrf)%name,itau_w,
     $                    pbl_tke(:,1:klev,nsrf))
       ENDIF

        IF (o_tke_max_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),
     $                    o_tke_max_srf(nsrf)%name,itau_w,
     $      pbl_tke(:,1:klev,nsrf))
        ENDIF
       endif
      ENDDO

        IF (o_cdrm%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_cdrm%name,itau_w,cdragm)
        ENDIF

        IF (o_cdrh%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_cdrh%name,itau_w,cdragh)
        ENDIF

        IF (o_cldl%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_cldl%name,itau_w,cldl)
        ENDIF

        IF (o_cldm%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_cldm%name,itau_w,cldm)
        ENDIF

        IF (o_cldh%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_cldh%name,itau_w,cldh)
        ENDIF

        IF (o_cldt%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_cldt%name, 
     &                   itau_w,cldt)
        ENDIF

        IF (o_cldq%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_cldq%name,itau_w,cldq)
        ENDIF

        IF (o_lwp%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(1:klon) = flwp(1:klon)
      CALL histwrite_phy(nid_files(iff),
     s                   o_lwp%name,itau_w,zx_tmp_fi2d)
        ENDIF

        IF (o_iwp%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(1:klon) = fiwp(1:klon)
      CALL histwrite_phy(nid_files(iff),
     s                    o_iwp%name,itau_w,zx_tmp_fi2d)
        ENDIF

        IF (o_ue%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_ue%name,itau_w,ue)
        ENDIF

        IF (o_ve%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_ve%name,itau_w,ve)
        ENDIF

        IF (o_uq%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_uq%name,itau_w,uq)
        ENDIF

        IF (o_vq%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_vq%name,itau_w,vq)
        ENDIF

      IF(iflag_con.GE.3) THEN ! sb
        IF (o_cape%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_cape%name,itau_w,cape)
        ENDIF

        IF (o_pbase%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_pbase%name,itau_w,ema_pcb)
        ENDIF

        IF (o_ptop%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_ptop%name,itau_w,ema_pct)
        ENDIF

        IF (o_fbase%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_fbase%name,itau_w,ema_cbmf)
        ENDIF

        IF (o_prw%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_prw%name,itau_w,prw)
        ENDIF

      IF (o_cape_max%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_cape_max%name,itau_w,cape)
      ENDIF

       IF (o_upwd%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_upwd%name,itau_w,upwd)
       ENDIF

       IF (o_Ma%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_Ma%name,itau_w,Ma)
       ENDIF

       IF (o_dnwd%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_dnwd%name,itau_w,dnwd)
       ENDIF

       IF (o_dnwd0%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_dnwd0%name,itau_w,dnwd0)
       ENDIF

       IF (o_ftime_con%flag(iff)<=lev_files(iff)) THEN
        zx_tmp_fi2d=float(itau_con)/float(itap)
      CALL histwrite_phy(nid_files(iff),o_ftime_con%name,
     s                   itau_w,zx_tmp_fi2d)
       ENDIF

       IF (o_mc%flag(iff)<=lev_files(iff)) THEN
        if(iflag_thermals.gt.1)then
         zx_tmp_fi3d=dnwd+dnwd0+upwd+fm_therm
        else
         zx_tmp_fi3d=dnwd+dnwd0+upwd
        endif 
      CALL histwrite_phy(nid_files(iff),o_mc%name,itau_w,zx_tmp_fi3d)
       ENDIF
      
      ENDIF !iflag_con .GE. 3

        IF (o_s_pblh%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_s_pblh%name,itau_w,s_pblh)
        ENDIF

        IF (o_s_pblt%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_s_pblt%name,itau_w,s_pblt)
        ENDIF

        IF (o_s_lcl%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_s_lcl%name,itau_w,s_lcl)
        ENDIF

        IF (o_s_therm%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_s_therm%name,itau_w,s_therm)
        ENDIF

!IM : Les champs suivants (s_capCL, s_oliqCL, s_cteiCL, s_trmb1, s_trmb2, s_trmb3) ne sont pas definis dans HBTM.F
!       IF (o_s_capCL%flag(iff)<=lev_files(iff)) THEN
!     CALL histwrite_phy(nid_files(iff),o_s_capCL%name,itau_w,s_capCL)
!       ENDIF

!       IF (o_s_oliqCL%flag(iff)<=lev_files(iff)) THEN
!     CALL histwrite_phy(nid_files(iff),o_s_oliqCL%name,itau_w,s_oliqCL)
!       ENDIF

!       IF (o_s_cteiCL%flag(iff)<=lev_files(iff)) THEN
!     CALL histwrite_phy(nid_files(iff),o_s_cteiCL%name,itau_w,s_cteiCL)
!       ENDIF

!       IF (o_s_trmb1%flag(iff)<=lev_files(iff)) THEN
!     CALL histwrite_phy(nid_files(iff),o_s_trmb1%name,itau_w,s_trmb1)
!       ENDIF

!       IF (o_s_trmb2%flag(iff)<=lev_files(iff)) THEN
!     CALL histwrite_phy(nid_files(iff),o_s_trmb2%name,itau_w,s_trmb2)
!       ENDIF

!       IF (o_s_trmb3%flag(iff)<=lev_files(iff)) THEN
!     CALL histwrite_phy(nid_files(iff),o_s_trmb3%name,itau_w,s_trmb3)
!       ENDIF

! Champs interpolles sur des niveaux de pression

        ll=0
        DO k=1, nlevSTD
!         IF(k.GE.2.AND.k.LE.12) bb2=clevSTD(k)
!         IF(k.GE.13.AND.k.LE.17) bb3=clevSTD(k)
         bb2=clevSTD(k) 
         IF(bb2.EQ."850".OR.bb2.EQ."700".OR.
     $      bb2.EQ."500".OR.bb2.EQ."200".OR.
     $      bb2.EQ."50".OR.bb2.EQ."10") THEN

! a refaire correctement !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ll=ll+1
       IF (o_uSTDlevs(ll)%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),o_uSTDlevs(ll)%name,
     &                    itau_w,uwriteSTD(:,k,iff))
       ENDIF

       IF (o_vSTDlevs(ll)%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_vSTDlevs(ll)%name,  
     &                   itau_w,vwriteSTD(:,k,iff))
       ENDIF

       IF (o_wSTDlevs(ll)%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_wSTDlevs(ll)%name,
     &                    itau_w,wwriteSTD(:,k,iff))
       ENDIF

       IF (o_zSTDlevs(ll)%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_zSTDlevs(ll)%name,
     &               itau_w,phiwriteSTD(:,k,iff))
       ENDIF

       IF (o_qSTDlevs(ll)%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_qSTDlevs(ll)%name,
     &                   itau_w, qwriteSTD(:,k,iff))
       ENDIF

       IF (o_tSTDlevs(ll)%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_tSTDlevs(ll)%name,
     &                   itau_w, twriteSTD(:,k,iff))
       ENDIF

       ENDIF !(bb2.EQ."850".OR.bb2.EQ."700".OR.
       ENDDO

      IF (o_t_oce_sic%flag(iff)<=lev_files(iff)) THEN
      DO i=1, klon
       IF (pctsrf(i,is_oce).GT.epsfra.OR.
     $     pctsrf(i,is_sic).GT.epsfra) THEN
        zx_tmp_fi2d(i) = (ftsol(i, is_oce) * pctsrf(i,is_oce)+
     $                   ftsol(i, is_sic) * pctsrf(i,is_sic))/
     $                   (pctsrf(i,is_oce)+pctsrf(i,is_sic))
       ELSE
        zx_tmp_fi2d(i) = 273.15
       ENDIF
      ENDDO
      CALL histwrite_phy(nid_files(iff),
     s                   o_t_oce_sic%name,itau_w,zx_tmp_fi2d)
      ENDIF

! Couplage convection-couche limite
      IF (iflag_con.GE.3) THEN
      IF (iflag_coupl.EQ.1) THEN
       IF (o_ale_bl%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),o_ale_bl%name,itau_w,ale_bl)
       ENDIF
       IF (o_alp_bl%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),o_alp_bl%name,itau_w,alp_bl)
       ENDIF
      ENDIF !iflag_coupl.EQ.1
      ENDIF !(iflag_con.GE.3)

! Wakes
      IF (iflag_con.EQ.3) THEN
      IF (iflag_wake.EQ.1) THEN
       IF (o_ale_wk%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),o_ale_wk%name,itau_w,ale_wake)
       ENDIF
       IF (o_alp_wk%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),o_alp_wk%name,itau_w,alp_wake)
       ENDIF

       IF (o_ale%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),o_ale%name,itau_w,ale)
       ENDIF
       IF (o_alp%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),o_alp%name,itau_w,alp)
       ENDIF
       IF (o_cin%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),o_cin%name,itau_w,cin)
       ENDIF
       IF (o_wape%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),o_WAPE%name,itau_w,wake_pe)
       ENDIF
       IF (o_wake_h%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_wake_h%name,itau_w,wake_h)
       ENDIF

       IF (o_wake_s%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_wake_s%name,itau_w,wake_s)
       ENDIF

        IF (o_wake_deltat%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),o_wake_deltat%name,
     $                     itau_w,wake_deltat)
        ENDIF

        IF (o_wake_deltaq%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),o_wake_deltaq%name,
     $                    itau_w,wake_deltaq)
        ENDIF

        IF (o_wake_omg%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),
     s                    o_wake_omg%name,itau_w,wake_omg)
        ENDIF

         IF (o_dtwak%flag(iff)<=lev_files(iff)) THEN
           zx_tmp_fi3d(1:klon,1:klev)=d_t_wake(1:klon,1:klev)
     &                                        /pdtphys
           CALL histwrite_phy(nid_files(iff),
     &                       o_dtwak%name,itau_w,zx_tmp_fi3d)
         ENDIF

        IF (o_dqwak%flag(iff)<=lev_files(iff)) THEN
        zx_tmp_fi3d(1:klon,1:klev)=d_q_wake(1:klon,1:klev)/pdtphys
        CALL histwrite_phy(nid_files(iff),
     &                     o_dqwak%name,itau_w,zx_tmp_fi3d)
        ENDIF
      ENDIF ! iflag_wake.EQ.1

        IF (o_Vprecip%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),o_Vprecip%name,itau_w,Vprecip)
        ENDIF

        IF (o_ftd%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),o_ftd%name,itau_w,ftd)
        ENDIF

        IF (o_fqd%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),o_fqd%name,itau_w,fqd)
        ENDIF
      ENDIF !(iflag_con.EQ.3) 
 
      IF (type_ocean=='slab ') THEN
      IF ( o_slab_bils%flag(iff)<=lev_files(iff)) 
     $     CALL histwrite_phy(
     $     nid_files(iff),o_slab_bils%name,itau_w,slab_wfbils)

      ENDIF !type_ocean == force/slab

      IF (o_weakinv%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     s                  o_weakinv%name,itau_w,weak_inversion)
      ENDIF

      IF (o_dthmin%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_dthmin%name,itau_w,dthmin)
      ENDIF

       IF (o_cldtau%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),o_cldtau%name,itau_w,cldtau)
       ENDIF

       IF (o_cldemi%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),o_cldemi%name,itau_w,cldemi)
       ENDIF

      IF (o_pr_con_l%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     s         o_pr_con_l%name,itau_w,pmflxr(:,1:klev))
      ENDIF

      IF (o_pr_con_i%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     s         o_pr_con_i%name,itau_w,pmflxs(:,1:klev))
      ENDIF

      IF (o_pr_lsc_l%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     s        o_pr_lsc_l%name,itau_w,prfl(:,1:klev))
      ENDIF

      IF (o_pr_lsc_i%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     s        o_pr_lsc_i%name,itau_w,psfl(:,1:klev))
      ENDIF

      IF (o_re%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_re%name,itau_w,re)
      ENDIF

      IF (o_fl%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_fl%name,itau_w,fl)
      ENDIF



      IF (o_rh2m%flag(iff)<=lev_files(iff)) THEN
      DO i=1, klon
       zx_tmp_fi2d(i)=MIN(100.,rh2m(i)*100.)
      ENDDO
      CALL histwrite_phy(nid_files(iff),o_rh2m%name,itau_w,zx_tmp_fi2d)
      ENDIF

      IF (o_rh2m_min%flag(iff)<=lev_files(iff)) THEN
      DO i=1, klon
       zx_tmp_fi2d(i)=MIN(100.,rh2m(i)*100.)
      ENDDO
      CALL histwrite_phy(nid_files(iff),o_rh2m_min%name,
     s               itau_w,zx_tmp_fi2d)
      ENDIF

      IF (o_rh2m_max%flag(iff)<=lev_files(iff)) THEN
      DO i=1, klon
       zx_tmp_fi2d(i)=MIN(100.,rh2m(i)*100.)
      ENDDO
      CALL histwrite_phy(nid_files(iff),o_rh2m_max%name,
     s              itau_w,zx_tmp_fi2d)
      ENDIF


      IF (o_qsat2m%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_qsat2m%name,itau_w,qsat2m)
      ENDIF

      IF (o_tpot%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_tpot%name,itau_w,tpot)
      ENDIF

       IF (o_tpote%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_tpote%name,itau_w,tpote)
       ENDIF

      IF (o_SWnetOR%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(1 : klon) = fsolsw( 1 : klon, is_ter)
      CALL histwrite_phy(nid_files(iff),
     s                   o_SWnetOR%name,itau_w, zx_tmp_fi2d)
      ENDIF

      IF (o_SWdownOR%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi2d(1:klon) = solsw(1:klon)/(1.-albsol1(1:klon))
      CALL histwrite_phy(nid_files(iff),
     s                   o_SWdownOR%name,itau_w, zx_tmp_fi2d)
      ENDIF

      IF (o_LWdownOR%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     s                  o_LWdownOR%name,itau_w,sollwdown)
      ENDIF

      IF (o_snowl%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_snowl%name,itau_w,snow_lsc)
      ENDIF

      IF (o_solldown%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     s                   o_solldown%name,itau_w,sollwdown)
      ENDIF

      IF (o_dtsvdfo%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     s                 o_dtsvdfo%name,itau_w,d_ts(:,is_oce))
      ENDIF

      IF (o_dtsvdft%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     s                   o_dtsvdft%name,itau_w,d_ts(:,is_ter))
      ENDIF

       IF (o_dtsvdfg%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),
     $                   o_dtsvdfg%name,itau_w, d_ts(:,is_lic))
       ENDIF

       IF (o_dtsvdfi%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     s                   o_dtsvdfi%name,itau_w,d_ts(:,is_sic))
       ENDIF

       IF (o_rugs%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_rugs%name,itau_w,zxrugs)
       ENDIF

! OD550 per species
      IF (new_aod .and. (.not. aerosol_couple)) THEN
          IF (ok_ade.OR.ok_aie) THEN

          IF (o_od550aer%flag(iff)<=lev_files(iff)) THEN
             CALL histwrite_phy(nid_files(iff),o_od550aer%name,itau_w,
     $            od550aer)
          ENDIF
          IF (o_od865aer%flag(iff)<=lev_files(iff)) THEN
             CALL histwrite_phy(nid_files(iff),o_od865aer%name,itau_w,
     $            od865aer)
          ENDIF
          IF (o_absvisaer%flag(iff)<=lev_files(iff)) THEN
             CALL histwrite_phy(nid_files(iff),o_absvisaer%name,itau_w,
     $            absvisaer)
          ENDIF
          IF (o_od550lt1aer%flag(iff)<=lev_files(iff)) THEN
            CALL histwrite_phy(nid_files(iff),o_od550lt1aer%name,itau_w,
     $            od550lt1aer)
          ENDIF

          IF (o_sconcso4%flag(iff)<=lev_files(iff)) THEN
              CALL histwrite_phy(nid_files(iff),o_sconcso4%name,itau_w,
     $            sconcso4)
          ENDIF
          IF (o_sconcoa%flag(iff)<=lev_files(iff)) THEN
              CALL histwrite_phy(nid_files(iff),o_sconcoa%name,itau_w,
     $            sconcoa)
          ENDIF
          IF (o_sconcbc%flag(iff)<=lev_files(iff)) THEN
              CALL histwrite_phy(nid_files(iff),o_sconcbc%name,itau_w,
     $            sconcbc)
          ENDIF
          IF (o_sconcss%flag(iff)<=lev_files(iff)) THEN
              CALL histwrite_phy(nid_files(iff),o_sconcss%name,itau_w,
     $            sconcss)
          ENDIF
          IF (o_sconcdust%flag(iff)<=lev_files(iff)) THEN
              CALL histwrite_phy(nid_files(iff),o_sconcdust%name,itau_w,
     $            sconcdust)
          ENDIF
          
          IF (o_concso4%flag(iff)<=lev_files(iff)) THEN
              CALL histwrite_phy(nid_files(iff),o_concso4%name,itau_w,
     $            concso4)
          ENDIF
          IF (o_concoa%flag(iff)<=lev_files(iff)) THEN
              CALL histwrite_phy(nid_files(iff),o_concoa%name,itau_w,
     $            concoa)
          ENDIF
          IF (o_concbc%flag(iff)<=lev_files(iff)) THEN
              CALL histwrite_phy(nid_files(iff),o_concbc%name,itau_w,
     $            concbc)
          ENDIF
          IF (o_concss%flag(iff)<=lev_files(iff)) THEN
              CALL histwrite_phy(nid_files(iff),o_concss%name,itau_w,
     $            concss)
          ENDIF
          IF (o_concdust%flag(iff)<=lev_files(iff)) THEN
              CALL histwrite_phy(nid_files(iff),o_concdust%name,itau_w,
     $            concdust)
          ENDIF
          
          IF (o_loadso4%flag(iff)<=lev_files(iff)) THEN
              CALL histwrite_phy(nid_files(iff),o_loadso4%name,itau_w,
     $            loadso4)
          ENDIF
          IF (o_loadoa%flag(iff)<=lev_files(iff)) THEN
              CALL histwrite_phy(nid_files(iff),o_loadoa%name,itau_w,
     $            loadoa)
          ENDIF
          IF (o_loadbc%flag(iff)<=lev_files(iff)) THEN
              CALL histwrite_phy(nid_files(iff),o_loadbc%name,itau_w,
     $            loadbc)
          ENDIF
          IF (o_loadss%flag(iff)<=lev_files(iff)) THEN
              CALL histwrite_phy(nid_files(iff),o_loadss%name,itau_w,
     $            loadss)
          ENDIF
          IF (o_loaddust%flag(iff)<=lev_files(iff)) THEN
              CALL histwrite_phy(nid_files(iff),o_loaddust%name,itau_w,
     $            loaddust)
          ENDIF
          
          DO naero = 1, naero_spc
            IF (o_tausumaero(naero)%flag(iff)<=lev_files(iff)) THEN
                CALL histwrite_phy(nid_files(iff),
     $              o_tausumaero(naero)%name,itau_w,
     $              tausum_aero(:,2,naero) )
            ENDIF
          END DO
          endif
      ENDIF
      
       IF (ok_ade) THEN
          IF (o_topswad%flag(iff)<=lev_files(iff)) THEN
             CALL histwrite_phy(nid_files(iff),o_topswad%name,itau_w,
     $            topswad_aero)
          ENDIF
          IF (o_solswad%flag(iff)<=lev_files(iff)) THEN
             CALL histwrite_phy(nid_files(iff),o_solswad%name,itau_w,
     $            solswad_aero)
          ENDIF

!====MS forcing diagnostics
        if (new_aod) then	       
        IF (o_swtoaas_nat%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_swtoaas_nat%name,itau_w,
     $      topsw_aero(:,1))
        ENDIF

        IF (o_swsrfas_nat%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_swsrfas_nat%name,itau_w,
     $      solsw_aero(:,1))
        ENDIF

        IF (o_swtoacs_nat%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_swtoacs_nat%name,itau_w,
     $      topsw0_aero(:,1))
        ENDIF

        IF (o_swsrfcs_nat%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_swsrfcs_nat%name,itau_w,
     $      solsw0_aero(:,1))
        ENDIF
  
!ant
        IF (o_swtoaas_ant%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_swtoaas_ant%name,itau_w,
     $      topsw_aero(:,2))
        ENDIF

        IF (o_swsrfas_ant%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_swsrfas_ant%name,itau_w,
     $      solsw_aero(:,2))
        ENDIF

        IF (o_swtoacs_ant%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_swtoacs_ant%name,itau_w,
     $      topsw0_aero(:,2))
        ENDIF

        IF (o_swsrfcs_ant%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_swsrfcs_ant%name,itau_w,
     $      solsw0_aero(:,2))
        ENDIF

!cf

        if (.not. aerosol_couple) then
        IF (o_swtoacf_nat%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_swtoacf_nat%name,itau_w,
     $      topswcf_aero(:,1))
        ENDIF

        IF (o_swsrfcf_nat%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_swsrfcf_nat%name,itau_w,
     $      solswcf_aero(:,1))
        ENDIF

        IF (o_swtoacf_ant%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_swtoacf_ant%name,itau_w,
     $      topswcf_aero(:,2))
        ENDIF

        IF (o_swsrfcf_ant%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_swsrfcf_ant%name,itau_w,
     $      solswcf_aero(:,2))
        ENDIF

        IF (o_swtoacf_zero%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_swtoacf_zero%name,itau_w,
     $      topswcf_aero(:,3))
        ENDIF

        IF (o_swsrfcf_zero%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_swsrfcf_zero%name,itau_w,
     $      solswcf_aero(:,3))
        ENDIF
        endif

	endif ! new_aod
!====MS forcing diagnostics

       ENDIF

       IF (ok_aie) THEN
          IF (o_topswai%flag(iff)<=lev_files(iff)) THEN
             CALL histwrite_phy(nid_files(iff),o_topswai%name,itau_w,
     $            topswai_aero)
          ENDIF
          IF (o_solswai%flag(iff)<=lev_files(iff)) THEN
             CALL histwrite_phy(nid_files(iff),o_solswai%name,itau_w,
     $            solswai_aero)
          ENDIF
          IF (o_scdnc%flag(iff)<=lev_files(iff)) THEN
             CALL histwrite_phy(nid_files(iff),o_scdnc%name,itau_w,
     $            scdnc)
          ENDIF
          IF (o_cldncl%flag(iff)<=lev_files(iff)) THEN
             CALL histwrite_phy(nid_files(iff),o_cldncl%name,itau_w,
     $            cldncl)
          ENDIF
          IF (o_reffclws%flag(iff)<=lev_files(iff)) THEN
             CALL histwrite_phy(nid_files(iff),o_reffclws%name,itau_w,
     $            reffclws)
          ENDIF
          IF (o_reffclwc%flag(iff)<=lev_files(iff)) THEN
             CALL histwrite_phy(nid_files(iff),o_reffclwc%name,itau_w,
     $            reffclwc)
          ENDIF
          IF (o_cldnvi%flag(iff)<=lev_files(iff)) THEN
             CALL histwrite_phy(nid_files(iff),o_cldnvi%name,itau_w,
     $            cldnvi)
          ENDIF
          IF (o_lcc%flag(iff)<=lev_files(iff)) THEN
             CALL histwrite_phy(nid_files(iff),o_lcc%name,itau_w,
     $            lcc)
          ENDIF
          IF (o_lcc3d%flag(iff)<=lev_files(iff)) THEN
             CALL histwrite_phy(nid_files(iff),o_lcc3d%name,itau_w,
     $            lcc3d)
          ENDIF
          IF (o_lcc3dcon%flag(iff)<=lev_files(iff)) THEN
             CALL histwrite_phy(nid_files(iff),o_lcc3dcon%name,itau_w,
     $            lcc3dcon)
          ENDIF
          IF (o_lcc3dstra%flag(iff)<=lev_files(iff)) THEN
             CALL histwrite_phy(nid_files(iff),o_lcc3dstra%name,itau_w,
     $            lcc3dstra)
          ENDIF
          IF (o_reffclwtop%flag(iff)<=lev_files(iff)) THEN
             CALL histwrite_phy(nid_files(iff),o_reffclwtop%name,itau_w,
     $            reffclwtop)
          ENDIF
       ENDIF

! Champs 3D:
       IF (ok_ade .OR. ok_aie) then
          IF (o_ec550aer%flag(iff)<=lev_files(iff)) THEN
             CALL histwrite_phy(nid_files(iff),o_ec550aer%name,itau_w,
     &            ec550aer)
          ENDIF
       ENDIF

       IF (o_lwcon%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_lwcon%name,itau_w,flwc)
       ENDIF

       IF (o_iwcon%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_iwcon%name,itau_w,fiwc)
       ENDIF

       IF (o_temp%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_temp%name,itau_w,t_seri)
       ENDIF

       IF (o_theta%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_theta%name,itau_w,theta)
       ENDIF

       IF (o_ovapinit%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_ovapinit%name,itau_w,
     $ qx(:,:,ivap))
       ENDIF

       IF (o_ovap%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     $                   o_ovap%name,itau_w,q_seri)
       ENDIF

       IF (o_geop%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_geop%name,itau_w,zphi)
       ENDIF

       IF (o_vitu%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_vitu%name,itau_w,u_seri)
       ENDIF

       IF (o_vitv%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_vitv%name,itau_w,v_seri)
       ENDIF

       IF (o_vitw%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_vitw%name,itau_w,omega)
       ENDIF

        IF (o_pres%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_pres%name,itau_w,pplay)
        ENDIF

        IF (o_paprs%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_paprs%name,
     s                    itau_w,paprs(:,1:klev))
        ENDIF

       IF (o_rneb%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_rneb%name,itau_w,cldfra)
       ENDIF

       IF (o_rnebcon%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_rnebcon%name,itau_w,rnebcon)
       ENDIF

       IF (o_rhum%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_rhum%name,itau_w,zx_rh)
       ENDIF

      IF (o_ozone%flag(iff)<=lev_files(iff)) THEN
         CALL histwrite_phy(nid_files(iff), o_ozone%name, itau_w,
     $        wo(:, :, 1) * dobson_u * 1e3 / zmasse / rmo3 * rmd)
      ENDIF

      IF (o_ozone_light%flag(iff)<=lev_files(iff) .and.
     $     read_climoz == 2) THEN
         CALL histwrite_phy(nid_files(iff), o_ozone_light%name, itau_w,
     $        wo(:, :, 2) * dobson_u * 1e3 / zmasse / rmo3 * rmd)
      ENDIF

       IF (o_dtphy%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_dtphy%name,itau_w,d_t)
       ENDIF

       IF (o_dqphy%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     s                  o_dqphy%name,itau_w, d_qx(:,:,ivap))
       ENDIF

        DO nsrf=1, nbsrf
        IF (o_albe_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN 
        zx_tmp_fi2d(1 : klon) = falb1( 1 : klon, nsrf)
        CALL histwrite_phy(nid_files(iff),
     s                    o_albe_srf(nsrf)%name,itau_w,
     $                     zx_tmp_fi2d)
        ENDIF

        IF (o_rugs_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN  
        zx_tmp_fi2d(1 : klon) = frugs( 1 : klon, nsrf)
        CALL histwrite_phy(nid_files(iff),
     s                     o_rugs_srf(nsrf)%name,itau_w,
     $      zx_tmp_fi2d)
        ENDIF

        IF (o_ages_srf(nsrf)%flag(iff)<=lev_files(iff)) THEN
        zx_tmp_fi2d(1 : klon) = agesno( 1 : klon, nsrf)
        CALL histwrite_phy(nid_files(iff),
     s                     o_ages_srf(nsrf)%name,itau_w
     $    ,zx_tmp_fi2d)
        ENDIF
        ENDDO !nsrf=1, nbsrf

       IF (o_alb1%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_alb1%name,itau_w,albsol1)
       ENDIF

       IF (o_alb2%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_alb2%name,itau_w,albsol2)
       ENDIF

!FH Sorties pour la couche limite
      if (iflag_pbl>1) then
      zx_tmp_fi3d=0.
      do nsrf=1,nbsrf
         do k=1,klev
          zx_tmp_fi3d(:,k)=zx_tmp_fi3d(:,k)
     $    +pctsrf(:,nsrf)*pbl_tke(:,k,nsrf)
         enddo
      enddo
       IF (o_tke%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_tke%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (o_tke_max%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),
     s                   o_tke_max%name,itau_w,zx_tmp_fi3d)
       ENDIF
      endif

       IF (o_kz%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_kz%name,itau_w,coefh)
       ENDIF

       IF (o_kz_max%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_kz_max%name,itau_w,coefh)
       ENDIF

       IF (o_clwcon%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_clwcon%name,itau_w,clwcon0)
       ENDIF

       IF (o_dtdyn%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_dtdyn%name,itau_w,d_t_dyn)
       ENDIF

       IF (o_dqdyn%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_dqdyn%name,itau_w,d_q_dyn)
       ENDIF

       IF (o_dudyn%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_dudyn%name,itau_w,d_u_dyn)
       ENDIF                                                    

       IF (o_dvdyn%flag(iff)<=lev_files(iff)) THEN                 
      CALL histwrite_phy(nid_files(iff),o_dvdyn%name,itau_w,d_v_dyn)  
       ENDIF                                                     

       IF (o_dtcon%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=d_t_con(1:klon,1:klev)/pdtphys
      CALL histwrite_phy(nid_files(iff),o_dtcon%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (o_ducon%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=d_u_con(1:klon,1:klev)/pdtphys
      CALL histwrite_phy(nid_files(iff),o_ducon%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (o_dqcon%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=d_q_con(1:klon,1:klev)/pdtphys
      CALL histwrite_phy(nid_files(iff),o_dqcon%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (o_dtlsc%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=d_t_lsc(1:klon,1:klev)/pdtphys
      CALL histwrite_phy(nid_files(iff),o_dtlsc%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (o_dtlschr%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon, 1:klev)=(d_t_lsc(1:klon,1:klev)+
     $                           d_t_eva(1:klon,1:klev))/pdtphys
      CALL histwrite_phy(nid_files(iff),
     s                   o_dtlschr%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (o_dqlsc%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=d_q_lsc(1:klon,1:klev)/pdtphys
      CALL histwrite_phy(nid_files(iff),o_dqlsc%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (o_dtvdf%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=d_t_vdf(1:klon,1:klev)/pdtphys
      CALL histwrite_phy(nid_files(iff),o_dtvdf%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (o_dqvdf%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=d_q_vdf(1:klon,1:klev)/pdtphys
      CALL histwrite_phy(nid_files(iff),o_dqvdf%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (o_dteva%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=d_t_eva(1:klon,1:klev)/pdtphys
      CALL histwrite_phy(nid_files(iff),o_dteva%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (o_dqeva%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=d_q_eva(1:klon,1:klev)/pdtphys
      CALL histwrite_phy(nid_files(iff),o_dqeva%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (o_ptconv%flag(iff)<=lev_files(iff)) THEN
      zpt_conv = 0.
      where (ptconv) zpt_conv = 1.
      CALL histwrite_phy(nid_files(iff),o_ptconv%name,itau_w,zpt_conv)
       ENDIF

       IF (o_ratqs%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_ratqs%name,itau_w,ratqs)
       ENDIF

       IF (o_dtthe%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=d_t_ajs(1:klon,1:klev)/pdtphys
      CALL histwrite_phy(nid_files(iff),o_dtthe%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (iflag_thermals.gt.1) THEN
        IF (o_ftime_th%flag(iff)<=lev_files(iff)) THEN
! Pour l instant 0 a y reflichir pour les thermiques
         zx_tmp_fi2d=0. 
        CALL histwrite_phy(nid_files(iff),o_ftime_th%name,
     s                     itau_w,zx_tmp_fi2d)
        ENDIF

        IF (o_f_th%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_f_th%name,itau_w,fm_therm)
        ENDIF

        IF (o_e_th%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_e_th%name,itau_w,entr_therm)
        ENDIF

        IF (o_w_th%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_w_th%name,itau_w,zw2)
        ENDIF

        IF (o_q_th%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_q_th%name,itau_w,zqasc)
        ENDIF

        IF (o_lambda_th%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),
     s                     o_lambda_th%name,itau_w,lambda_th)
        ENDIF

        IF (o_a_th%flag(iff)<=lev_files(iff)) THEN
        CALL histwrite_phy(nid_files(iff),o_a_th%name,itau_w,fraca)
        ENDIF

       IF (o_d_th%flag(iff)<=lev_files(iff)) THEN
       CALL histwrite_phy(nid_files(iff),o_d_th%name,itau_w,detr_therm)
       ENDIF


       IF (o_f0_th%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_f0_th%name,itau_w,f0)
       ENDIF

       IF (o_f0_th%flag(iff)<=lev_files(iff)) THEN
      CALL histwrite_phy(nid_files(iff),o_zmax_th%name,itau_w,zmax0)
       ENDIF

       IF (o_dqthe%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=d_q_ajs(1:klon,1:klev)/pdtphys
      CALL histwrite_phy(nid_files(iff),o_dqthe%name,itau_w,zx_tmp_fi3d)
       ENDIF

      ENDIF !iflag_thermals

       IF (o_dtajs%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=d_t_ajsb(1:klon,1:klev)/pdtphys
      CALL histwrite_phy(nid_files(iff),o_dtajs%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (o_dqajs%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=d_q_ajsb(1:klon,1:klev)/pdtphys
      CALL histwrite_phy(nid_files(iff),o_dqajs%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (o_dtswr%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=heat(1:klon,1:klev)/RDAY
      CALL histwrite_phy(nid_files(iff),o_dtswr%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (o_dtsw0%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=heat0(1:klon,1:klev)/RDAY
      CALL histwrite_phy(nid_files(iff),o_dtsw0%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (o_dtlwr%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=-1.*cool(1:klon,1:klev)/RDAY
      CALL histwrite_phy(nid_files(iff),o_dtlwr%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (o_dtlw0%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=-1.*cool0(1:klon,1:klev)/RDAY
      CALL histwrite_phy(nid_files(iff),o_dtlw0%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (o_dtec%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=d_t_ec(1:klon,1:klev)
      CALL histwrite_phy(nid_files(iff),o_dtec%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (o_duvdf%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=d_u_vdf(1:klon,1:klev)/pdtphys
      CALL histwrite_phy(nid_files(iff),o_duvdf%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (o_dvvdf%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=d_v_vdf(1:klon,1:klev)/pdtphys
      CALL histwrite_phy(nid_files(iff),o_dvvdf%name,itau_w,zx_tmp_fi3d)
       ENDIF

       IF (ok_orodr) THEN
      IF (o_duoro%flag(iff)<=lev_files(iff)) THEN 
      zx_tmp_fi3d(1:klon,1:klev)=d_u_oro(1:klon,1:klev)/pdtphys
      CALL histwrite_phy(nid_files(iff),o_duoro%name,itau_w,zx_tmp_fi3d)
       ENDIF

      IF (o_dvoro%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=d_v_oro(1:klon,1:klev)/pdtphys
      CALL histwrite_phy(nid_files(iff),o_dvoro%name,itau_w,zx_tmp_fi3d)
      ENDIF
       ENDIF

        IF (ok_orolf) THEN
       IF (o_dulif%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=d_u_lif(1:klon,1:klev)/pdtphys
      CALL histwrite_phy(nid_files(iff),o_dulif%name,itau_w,zx_tmp_fi3d)
       ENDIF

        IF (o_dvlif%flag(iff)<=lev_files(iff)) THEN
      zx_tmp_fi3d(1:klon,1:klev)=d_v_lif(1:klon,1:klev)/pdtphys
      CALL histwrite_phy(nid_files(iff),o_dvlif%name,itau_w,zx_tmp_fi3d)
       ENDIF
        ENDIF

        if (nqtot.GE.3) THEN
         DO iq=3,nqtot
       IF (o_trac(iq-2)%flag(iff)<=lev_files(iff)) THEN
         CALL histwrite_phy(nid_files(iff),
     s                  o_trac(iq-2)%name,itau_w,qx(:,:,iq))
       ENDIF
         ENDDO
        endif

      if (ok_sync) then
c$OMP MASTER
        call histsync(nid_files(iff))
c$OMP END MASTER
      endif

       ENDIF ! clef_files

      ENDDO ! iff=1,nfiles
