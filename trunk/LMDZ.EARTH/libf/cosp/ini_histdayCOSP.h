! Abderrahmane Idelkadi Septebmre 2009
! Sorties journalieres de COSP 
! Pour l'instant sorties Lidar et ISCCP
!
! sorties par jour
!
!$OMP MASTER
        zstoday = ecrit_day
        zout = freq_COSP
!
!       PRINT*, 'La frequence de sortie hf3d est de ', ecrit_hf
!

        idayref = day_ref
        CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)

        CALL histbeg_phy("histdayCOSP",itau_phy,zjulian,dtime,nhori,nid_day_cosp) 

! Definition de l'axe vertical
        CALL histvert(nid_day_cosp,"height","height","m",Nlevout,vgrid%z,nvert)
        print*,'Ok height Nlevout, height =',Nlevout,vgrid%z
        CALL histvert(nid_day_cosp,"height_mlev","height_mlev","m",Nlevlmdz,vgrid%mz,nvertm)
        print*,'Ok height_mlev Nlevout, height_mlev =',Nlevout,vgrid%mz
!        CALL histvert(nid_day_cosp,"presnivs","Vertical levels","mb",Nlevout,presnivs,nvert)

        CALL histvert(nid_day_cosp,"sza","solar_zenith_angle","degrees",PARASOL_NREFL,PARASOL_SZA,nvertp)

        CALL histvert(nid_day_cosp,"pressure2","pressure","mb",7,ISCCP_PC,nvertisccp)

        CALL histvert(nid_day_cosp,"column","column","count",Ncolumns,column_ax,nvertcol)

! Sorties LIDAR
       if (cfg%Llidar_sim) then
         if (cfg%Lcllcalipso) then
         CALL histdef(nid_day_cosp, "cllcalipso", &
                     "Lidar Low-level Cloud Fraction", "1", &
                     iim, jj_nb,nhori,1,1,1,-99,32, &
                     "ave(X)", zout,zstoday)
         endif
         if (cfg%Lclhcalipso) then
         CALL histdef(nid_day_cosp, "clhcalipso", &
                     "Lidar High-level Cloud Fraction", "1", &
                     iim, jj_nb,nhori,1,1,1,-99,32, &
                     "ave(X)", zout,zstoday)
         endif
         if (cfg%Lclmcalipso) then
         CALL histdef(nid_day_cosp, "clmcalipso", &
                     "Lidar Mid-level Cloud Fraction", "1", &
                     iim, jj_nb,nhori,1,1,1,-99,32, &
                     "ave(X)", zout,zstoday)
         endif
         if (cfg%Lcltcalipso) then
         CALL histdef(nid_day_cosp, "cltcalipso", &
                     "Lidar Total Cloud Fraction", "1", &
                     iim, jj_nb,nhori,1,1,1,-99,32, &
                     "ave(X)", zout,zstoday)
         endif
         if (cfg%Lclcalipso) then
         CALL histdef(nid_day_cosp, "clcalipso", &
                     "Lidar Cloud Fraction (532 nm)", "1", &
                     iim,jj_nb,nhori, Nlevout,1,Nlevout,nvert, 32, &
                     "ave(X)", zout,zstoday)
         endif
           if (cfg%Lcfad_lidarsr532) then
              do ii=1,SR_BINS
               CALL histdef(nid_day_cosp, "cfad_lidarsr532_"//chcol(ii), &
                           "Lidar Scattering Ratio CFAD (532 nm)","1", &
                           iim,jj_nb,nhori, Nlevout,1,Nlevout,nvert, 32, &
                           "ave(X)", zout,zstoday)   
              enddo
           endif
           if (cfg%Lparasol_refl) then
            CALL histdef(nid_day_cosp, "parasol_refl", &
                        "PARASOL-like mono-directional reflectance","1", &
                        iim,jj_nb,nhori, PARASOL_NREFL,1, PARASOL_NREFL, nvertp,32, &
                        "ave(X)", zout,zstoday)   
           endif
           if (cfg%Latb532) then
            do ii=1,Ncolumns
             CALL histdef(nid_day_cosp, "atb532_"//chcol(ii), &
                         "Lidar Attenuated Total Backscatter (532 nm)","1", &
                         iim,jj_nb,nhori, Nlevlmdz,1,Nlevlmdz,nvertm, 32, &
                         "ave(X)", zout,zstoday)
            enddo
           endif
           if (cfg%Lbeta_mol532) then
            CALL histdef(nid_day_cosp, "beta_mol532", &
                        "Lidar Molecular Backscatter (532 nm)","m-1 sr-1", &
                        iim,jj_nb,nhori, Nlevlmdz,1,Nlevlmdz,nvertm, 32, &
                         "ave(X)", zout,zstoday)
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
             CALL histdef(nid_day_cosp, "clisccp2_"//chcol(ii), &
                         "Cloud Fraction as Calculated by the ISCCP Simulator","1", &
                         iim,jj_nb,nhori,7,1,7,nvertisccp, 32, &
                         "ave(X)", zout,zstoday)
            enddo
          endif
          if (cfg%Lboxtauisccp) then
            CALL histdef(nid_day_cosp, "boxtauisccp", &
                         "Optical Depth in Each Column as Calculated by the ISCCP Simulator","1", &
                         iim,jj_nb,nhori,Ncolumns,1,Ncolumns,nvertcol, 32, &
                         "ave(X)", zout,zstoday)
          endif
          if (cfg%Lboxptopisccp) then
            CALL histdef(nid_day_cosp, "boxptopisccp", &
                         "Cloud Top Pressure in Each Column as Calculated by the ISCCP Simulator","Pa", &
                         iim,jj_nb,nhori,Ncolumns,1,Ncolumns,nvertcol, 32, &
                         "ave(X)", zout,zstoday)
          endif
          if (cfg%Ltclisccp) then
           CALL histdef(nid_day_cosp, "tclisccp", &
                     "Total Cloud Fraction as Calculated by the ISCCP Simulator", "1", &
                     iim, jj_nb,nhori,1,1,1,-99,32, &
                     "ave(X)", zout,zstoday)
          endif
          if (cfg%Lctpisccp) then
            CALL histdef(nid_day_cosp, "ctpisccp", &
                     "Mean Cloud Top Pressure as Calculated by the ISCCP Simulator", "Pa", &
                     iim, jj_nb,nhori,1,1,1,-99,32, &
                     "ave(X)", zout,zstoday)
          endif
          if (cfg%Ltauisccp) then
           CALL histdef(nid_day_cosp, "tauisccp", &
                     "Optical Depth as Calculated by the ISCCP Simulator", "1", &
                     iim, jj_nb,nhori,1,1,1,-99,32, &
                     "ave(X)", zout,zstoday)
          endif
          if (cfg%Lalbisccp) then
           CALL histdef(nid_day_cosp, "albisccp", &
                     "Mean Cloud Albedo as Calculated by the ISCCP Simulator", "1", &
                     iim, jj_nb,nhori,1,1,1,-99,32, &
                     "ave(X)", zout,zstoday) 
          endif
          if (cfg%Lmeantbisccp) then
            CALL histdef(nid_day_cosp, "meantbisccp", &
             " Mean all-sky 10.5 micron brightness temperature as calculated by the ISCCP Simulator","K", &
             iim, jj_nb,nhori,1,1,1,-99,32, &
             "ave(X)", zout,zstoday)
          endif
          if (cfg%Lmeantbclrisccp) then
           CALL histdef(nid_day_cosp, "meantbclrisccp", &
            "Mean clear-sky 10.5 micron brightness temperature as calculated by the ISCCP Simulator","K", &
             iim, jj_nb,nhori,1,1,1,-99,32, &
             "ave(X)", zout,zstoday) 
          endif
        endif ! Isccp


        CALL histend(nid_day_cosp)
!$OMP END MASTER
!$OMP BARRIER
