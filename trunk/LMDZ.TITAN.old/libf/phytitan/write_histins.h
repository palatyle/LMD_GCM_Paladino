!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/write_histins.h,v 1.1.1.1 2004/05/19 12:53:09 lmdzadmin Exp $
!
      IF (ok_instan) THEN

         zsto = dtime * REAL(ecrit_ins)
         zout = dtime * REAL(ecrit_ins)
         itau_w = itau_phy + itap

c-------------------------------------------------------
      IF(lev_histday.GE.1) THEN

ccccccccccccc 2D fields, invariables

      call histwrite_phy(nid_ins,.false.,"phis",itau_w,pphis)
c      call histwrite_phy(nid_ins,.false.,"aire",itau_w,cell_area)
      cell_area_out(:)=cell_area(:)
      if (is_north_pole_phy) cell_area_out(1)=cell_area(1)/nbp_lon
      if (is_south_pole_phy) cell_area_out(klon)=cell_area(klon)/nbp_lon
      call histwrite_phy(nid_ins,.false.,"aire",itau_w,cell_area_out)

ccccccc axe Ls ... Faudrait le reduire a axe temporel seulement...
      do i=1,klon
        tmpout(i,1) = zls*180./RPI
      enddo
      call histwrite_phy(nid_ins,.false.,"ls",itau_w,tmpout(:,1))

ccccccccccccc 2D fields, variables

      call histwrite_phy(nid_ins,.false.,"tsol",itau_w,ftsol)
      call histwrite_phy(nid_ins,.false.,"psol",itau_w,paprs(:,1))

c     call histwrite_phy(nid_ins,.false.,"ue",itau_w,ue)
c     call histwrite_phy(nid_ins,.false.,"ve",itau_w,ve)

      ENDIF !lev_histday.GE.1

c-------------------------------------------------------
      IF(lev_histday.GE.2) THEN

ccccccccccccc 3D fields, basics

      call histwrite_phy(nid_ins,.false.,"temp",itau_w,t_seri)
      call histwrite_phy(nid_ins,.false.,"pres",itau_w,pplay)
      call histwrite_phy(nid_ins,.false.,"geop",itau_w,zphi)
      call histwrite_phy(nid_ins,.false.,"vitu",itau_w,u_seri)
      call histwrite_phy(nid_ins,.false.,"vitv",itau_w,v_seri)
      call histwrite_phy(nid_ins,.false.,"vitw",itau_w,omega)
      call histwrite_phy(nid_ins,.false.,"tops",itau_w,topsw)
c      call histwrite_phy(nid_ins,.false.,"duvdf",itau_w,d_u_vdf)
c      call histwrite_phy(nid_ins,.false.,"dudyn",itau_w,d_u_dyn)

      ENDIF !lev_histday.GE.2

c-------------------------------------------------------
      IF(lev_histday.GE.3) THEN

cccccccccccccccccc  Tracers

         if (iflag_trac.eq.1) THEN
          if (microfi.eq.1) then
           DO iq=1,nmicro
       call histwrite_phy(nid_ins,.false.,tname(iq),
     .                    itau_w,qaer(1:klon,1:klev,iq))
           ENDDO
	  endif
	  if (nmicro.lt.nqmax) then
           DO iq=nmicro+1,nqmax
       call histwrite_phy(nid_ins,.false.,tname(iq),
     .                    itau_w,tr_seri(1:klon,1:klev,iq))
           ENDDO
	  endif
         endif

cccccccccccccccccc  Radiative transfer

c 2D

      call histwrite_phy(nid_ins,.false.,"topl",itau_w,toplw)
      call histwrite_phy(nid_ins,.false.,"sols",itau_w,solsw)
      call histwrite_phy(nid_ins,.false.,"soll",itau_w,sollw)

c 3D

      call histwrite_phy(nid_ins,.false.,"SWnet",
     .          itau_w,swnet(1:klon,1:klev))
      call histwrite_phy(nid_ins,.false.,"LWnet",
     .          itau_w,lwnet(1:klon,1:klev))

c --------------
c ----- OPACITE BRUME
       do k=7,NSPECV,10
         do i=1,klon
         do l=1,klev
           t_tauhvd(i,l)=TAUHVD(i,klev-l+1,k)
         enddo
         enddo
         write(str2,'(i2.2)') k
       call histwrite_phy(nid_ins,.false.,"thv"//str2,itau_w,t_tauhvd)
       enddo      ! fin boucle NSPECV 

       do k=8,NSPECI,10
         do i=1,klon
         do l=1,klev
           t_tauhvd(i,l)=TAUHID(i,klev-l+1,k)
         enddo
         enddo
         write(str2,'(i2.2)') k
       call histwrite_phy(nid_ins,.false.,"thi"//str2,itau_w,t_tauhvd)
       enddo      ! fin boucle NSPECI 
c --------------
c ----- EXTINCTION BRUME
       do k=7,NSPECV,10
         do i=1,klon
         do l=1,klev
          if(l.ne.klev)
     s     t_khvd(i,l)=TAUHVD(i,klev-l+1,k)
     s                -TAUHVD(i,klev-l+1-1,k)
          if(l.eq.klev)
     s     t_khvd(i,l)=TAUHVD(i,klev-l+1,k)

         t_khvd(i,l)=t_khvd(i,l)/(zzlev(i,l+1)-zzlev(i,l))
         enddo
         enddo
         write(str2,'(i2.2)') k
       call histwrite_phy(nid_ins,.false.,"khv"//str2,itau_w,t_khvd)
       enddo      ! fin boucle NSPECV 

       do k=8,NSPECI,10
         do i=1,klon
         do l=1,klev
          if(l.ne.klev)
     s     t_khvd(i,l)=TAUHID(i,klev-l+1,k)
     s                -TAUHID(i,klev-l+1-1,k)
          if(l.eq.klev)
     s     t_khvd(i,l)=TAUHID(i,klev-l+1,k)

         t_khvd(i,l)=t_khvd(i,l)/(zzlev(i,l+1)-zzlev(i,l))
         enddo
         enddo
         write(str2,'(i2.2)') k
       call histwrite_phy(nid_ins,.false.,"khi"//str2,itau_w,t_khvd)
       enddo      ! fin boucle NSPECI 
c --------------
c ----- OPACITE GAZ
       do k=7,NSPECV,10
         do i=1,klon
         do l=1,klev
           t_tauhvd(i,l)=TAUGVD(i,klev-l+1,k)
         enddo
         enddo
         write(str2,'(i2.2)') k
       call histwrite_phy(nid_ins,.false.,"tgv"//str2,itau_w,t_tauhvd)
       enddo      ! fin boucle NSPECV 

       do k=8,NSPECI,10
         do i=1,klon
         do l=1,klev
           t_tauhvd(i,l)=TAUGID(i,klev-l+1,k)
         enddo
         enddo
         write(str2,'(i2.2)') k
       call histwrite_phy(nid_ins,.false.,"tgi"//str2,itau_w,t_tauhvd)
       enddo      ! fin boucle NSPECI 
c --------------
c ----- EXTINCTION GAZ
       do k=7,NSPECV,10
         do i=1,klon
         do l=1,klev
          if(l.ne.klev)
     s     t_khvd(i,l)=TAUGVD(i,klev-l+1,k)
     s                -TAUGVD(i,klev-l+1-1,k)
          if(l.eq.klev)
     s     t_khvd(i,l)=TAUGVD(i,klev-l+1,k)

         t_khvd(i,l)=t_khvd(i,l)/(zzlev(i,l+1)-zzlev(i,l))
         enddo
         enddo
         write(str2,'(i2.2)') k
       call histwrite_phy(nid_ins,.false.,"kgv"//str2,itau_w,t_khvd)
       enddo      ! fin boucle NSPECV 

       do k=8,NSPECI,10
         do i=1,klon
         do l=1,klev
          if(l.ne.klev)
     s     t_khvd(i,l)=TAUGID(i,klev-l+1,k)
     s                -TAUGID(i,klev-l+1-1,k)

          if(l.eq.klev)
     s     t_khvd(i,l)=TAUGID(i,klev-l+1,k)

         t_khvd(i,l)=t_khvd(i,l)/(zzlev(i,l+1)-zzlev(i,l))
         enddo
         enddo
         write(str2,'(i2.2)') k
       call histwrite_phy(nid_ins,.false.,"kgi"//str2,itau_w,t_khvd)
       enddo      ! fin boucle NSPECI 

      ENDIF !lev_histday.GE.3

c-------------------------------------------------------
      IF(lev_histday.GE.4) THEN

      call histwrite_phy(nid_ins,.false.,"dtdyn",itau_w,d_t_dyn)
      call histwrite_phy(nid_ins,.false.,"dtphy",itau_w,d_t)
c K/s
      call histwrite_phy(nid_ins,.false.,"dtvdf",itau_w,d_t_vdf)
c K/s
      call histwrite_phy(nid_ins,.false.,"dtajs",itau_w,d_t_ajs)
c K/s
      call histwrite_phy(nid_ins,.false.,"dtswr",itau_w,heat)
c K/s
      call histwrite_phy(nid_ins,.false.,"dtlwr",itau_w,-1.*cool)
c K/s      
c      call histwrite_phy(nid_ins,.false.,"dtec",itau_w,d_t_ec)
c      call histwrite_phy(nid_ins,.false.,"dvvdf",itau_w,d_v_vdf)

      ENDIF !lev_histday.GE.4

c-------------------------------------------------------
      IF(lev_histday.GE.5) THEN

c      call histwrite_phy(nid_ins,.false.,"taux",itau_w,fluxu)
c      call histwrite_phy(nid_ins,.false.,"tauy",itau_w,fluxv)
c      call histwrite_phy(nid_ins,.false.,"cdrm",itau_w,cdragm)
c      call histwrite_phy(nid_ins,.false.,"cdrh",itau_w,cdragh)

      ENDIF !lev_histday.GE.5
c-------------------------------------------------------

      if (ok_sync) then
        call histsync(nid_ins)
      endif
      ENDIF
