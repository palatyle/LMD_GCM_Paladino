! Ecriture des fichiers de sorties COSP
! Sorties journalierres
! Abderrahmane Idelkadi Septembre 2009

      IF (MOD(itap,NINT(freq_COSP/dtime)).EQ.0) THEN

       itau_wcosp = itau_phy + itap

! Sorties LIDAR
       if (cfg%Llidar_sim) then
         if (cfg%Lcllcalipso) then
          CALL histwrite_phy(nid_day_cosp,"cllcalipso",itau_wcosp,stlidar%cldlayer(:,1))
         endif
         if (cfg%Lclhcalipso) then
          CALL histwrite_phy(nid_day_cosp,"clhcalipso",itau_wcosp,stlidar%cldlayer(:,3))
         endif
         if (cfg%Lclmcalipso) then
          CALL histwrite_phy(nid_day_cosp,"clmcalipso",itau_wcosp,stlidar%cldlayer(:,2))
         endif
         if (cfg%Lcltcalipso) then
          CALL histwrite_phy(nid_day_cosp,"cltcalipso",itau_wcosp,stlidar%cldlayer(:,4)) 
         endif
         if (cfg%Lclcalipso) then
          CALL histwrite_phy(nid_day_cosp,"clcalipso",itau_wcosp,stlidar%lidarcld)
         endif
         if (cfg%Lcfad_lidarsr532) then
           do ii=1,SR_BINS
            CALL histwrite_phy(nid_day_cosp,"cfad_lidarsr532_"//chcol(ii),itau_wcosp,stlidar%cfad_sr(:,ii,:))
           enddo
         endif
         if (cfg%Lparasol_refl) then
           CALL histwrite_phy(nid_day_cosp,"parasol_refl",itau_wcosp,stlidar%parasolrefl)
         endif
         if (cfg%Latb532) then
           do ii=1,Ncolumns
            CALL histwrite_phy(nid_day_cosp,"atb532_"//chcol(ii),itau_wcosp,sglidar%beta_tot(:,ii,:))
           enddo
         endif
         if (cfg%Lbeta_mol532) then
           CALL histwrite_phy(nid_day_cosp,"beta_mol532",itau_wcosp,sglidar%beta_mol)
         endif
        endif ! Lidar

! Sorties RADAR
!Attention A FAIRE
!        if (cfg%Lradar_sim) then
!         print*,'Ecriture sorties Radar'
!          if (cfg%Lcfad_dbze94) then
!              print*,'Ecriture de cfad_dbze94.nc '
!              A revoir l axe vertical Nlvgrid
!               do ii=1,DBZE_BINS
!                   dbze_ax(ii) = CFAD_ZE_MIN + CFAD_ZE_WIDTH*(ii - 0.5)
!               enddo
!               call write_netcdf4d('cfad_dbze94.nc',use_vgrid,nlon,nlat,Nlevout,DBZE_BINS, &
!                                   x,y,out_levs,dbze_ax,i,ndays,time,stradar%cfad_ze)
!          endif
!          if (cfg%Lclcalipso2) then
!               call write_netcdf3d('clcalipso2.nc',use_vgrid,'clcalipso2', &
!                              nlon,nlat,Nlevout,x,y,out_levs,i,ndays,time,stradar%lidar_only_freq_cloud)
!          endif
!          if (cfg%Ldbze94) then
!             do ii=1,Ncolumns
!                xcol(ii)=float(i)
!             enddo
!             call write_netcdf4d('dbze94.nc',use_vgrid,nlon,nlat,Nlevout,Ncolumns, &
!                                 x,y,out_levs,xcol,i,ndays,time,sgradar%Ze_tot)
!          endif
!          if (cfg%Lcltlidarradar) then
!             call write_netcdf2d('cltlidarradar.nc','cltlidarradar', &
!                                 nlon,nlat,x,y,i,ndays,time,stradar%radar_lidar_tcc)
!          endif
!        endif  ! Radar

! Sorties MISR
!Attention A FAIRE
!        if (cfg%Lmisr_sim) then
!         print*,'Ecriture sorties Misr'
!            call write_netcdf4d('clMISR.nc',use_vgrid,nlon,nlat,MISR_N_CTH,7, &
!                                x,y,MISR_CTH,ISCCP_TAU,i,ndays,time,misr%fq_MISR)
!        endif

! Sorties ISCCP
        if (cfg%Lisccp_sim) then
          if (cfg%Lclisccp2) then
            do ii=1,7
              CALL histwrite_phy(nid_day_cosp,"clisccp2_"//chcol(ii),itau_wcosp,isccp%fq_isccp(:,ii,:))
            enddo
          endif
          if (cfg%Lboxtauisccp) then
             CALL histwrite_phy(nid_day_cosp,"boxtauisccp",itau_wcosp,isccp%boxtau)
          endif
          if (cfg%Lboxptopisccp) then
             CALL histwrite_phy(nid_day_cosp,"boxptopisccp",itau_wcosp,isccp%boxptop)
          endif
          if (cfg%Ltclisccp) then
             CALL histwrite_phy(nid_day_cosp,"tclisccp",itau_wcosp,isccp%totalcldarea)
          endif
          if (cfg%Lctpisccp) then
             CALL histwrite_phy(nid_day_cosp,"ctpisccp",itau_wcosp,isccp%meanptop)
          endif
          if (cfg%Ltauisccp) then
             CALL histwrite_phy(nid_day_cosp,"tauisccp",itau_wcosp,isccp%meantaucld)
          endif
          if (cfg%Lalbisccp) then
             CALL histwrite_phy(nid_day_cosp,"albisccp",itau_wcosp,isccp%meanalbedocld)
          endif
          if (cfg%Lmeantbisccp) then
             CALL histwrite_phy(nid_day_cosp,"meantbisccp",itau_wcosp,isccp%meantb)
          endif
          if (cfg%Lmeantbclrisccp) then
             CALL histwrite_phy(nid_day_cosp,"meantbclrisccp",itau_wcosp,isccp%meantbclr)
          endif
        endif ! Isccp

!       if (ok_sync) then
!$OMP MASTER
        call histsync(nid_day_cosp)
!$OMP END MASTER      
!       endif

      ENDIF ! if freq_COSP
