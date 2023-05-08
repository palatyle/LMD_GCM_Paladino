!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/write_histday.h,v 1.2 2004/06/01 09:27:10 lmdzadmin Exp $
!
      IF (ok_journe) THEN

         zsto = dtime
         zout = dtime * REAL(ecrit_day)
         itau_w = itau_phy + itap

c-------------------------------------------------------
      IF(lev_histday.GE.1) THEN

ccccccccccccc 2D fields, invariables

      call histwrite_phy(nid_day,.false.,"phis",itau_w,pphis)
c      call histwrite_phy(nid_day,.false.,"aire",itau_w,cell_area)
      cell_area_out(:)=cell_area(:)
      if (is_north_pole_phy) cell_area_out(1)=cell_area(1)/nbp_lon
      if (is_south_pole_phy) cell_area_out(klon)=cell_area(klon)/nbp_lon
      call histwrite_phy(nid_day,.false.,"aire",itau_w,cell_area_out)

ccccccc axe Ls ... Faudrait le reduire a axe temporel seulement...
c Correction passage de 360 a 0... Sinon probleme avec moyenne
      if (zls.lt.zlsm1) then
        do i=1,klon
          tmpout(i,1) = zls*180./RPI+360.
        enddo
        zlsm1 = 2.*RPI
      else
        do i=1,klon
          tmpout(i,1) = zls*180./RPI
        enddo
        zlsm1 = zls
      endif
      call histwrite_phy(nid_day,.false.,"ls",itau_w,tmpout(:,1))

ccccccccccccc 2D fields, variables

      call histwrite_phy(nid_day,.false.,"tsol",itau_w,ftsol)
      call histwrite_phy(nid_day,.false.,"psol",itau_w,paprs(:,1))

c     call histwrite_phy(nid_day,.false.,"ue",itau_w,ue)
c     call histwrite_phy(nid_day,.false.,"ve",itau_w,ve)

      ENDIF !lev_histday.GE.1

c-------------------------------------------------------
      IF(lev_histday.GE.2) THEN

ccccccccccccc 3D fields, basics

      call histwrite_phy(nid_day,.false.,"temp",itau_w,t_seri)
      call histwrite_phy(nid_day,.false.,"pres",itau_w,pplay)
      call histwrite_phy(nid_day,.false.,"geop",itau_w,zphi)
      call histwrite_phy(nid_day,.false.,"vitu",itau_w,u_seri)
      call histwrite_phy(nid_day,.false.,"vitv",itau_w,v_seri)
      call histwrite_phy(nid_day,.false.,"vitw",itau_w,omega)
      call histwrite_phy(nid_day,.false.,"tops",itau_w,topsw)
      call histwrite_phy(nid_day,.false.,"duvdf",itau_w,d_u_vdf)
      call histwrite_phy(nid_day,.false.,"dudyn",itau_w,d_u_dyn)

cccccccccccccccccc  Tracers

         if (iflag_trac.eq.1) THEN
          if (microfi.ge.1) then
c          DO iq=1,nmicro
c      call histwrite_phy(nid_day,.false.,tname(iq),
c    .                    itau_w,qaer(1:klon,1:klev,iq))
c          ENDDO
c    -------   NB AER TOT
               do i=1,klon
                 do j=1,klev
                   tmpout(i,j)= SUM(qaer(i,j,1:nrad))
                 enddo
               enddo
       call histwrite_phy(nid_day,.false.,"qaer",itau_w,tmpout)

             if (clouds.eq.1) then
c    -------   NB NOY TOT
               do i=1,klon
                 do j=1,klev
                   tmpout(i,j)= SUM(qaer(i,j,nrad+1:2*nrad))
                 enddo
               enddo
       call histwrite_phy(nid_day,.false.,"qnoy",itau_w,tmpout)
c    -------   V GLA1 TOT
               do i=1,klon
                 do j=1,klev
                   tmpout(i,j)= SUM(qaer(i,j,2*nrad+1:3*nrad))
                 enddo
               enddo
       call histwrite_phy(nid_day,.false.,"qgl1",itau_w,tmpout)
c    -------   V GLA2 TOT
               do i=1,klon
                 do j=1,klev
                   tmpout(i,j)= SUM(qaer(i,j,3*nrad+1:4*nrad))
                 enddo
               enddo
       call histwrite_phy(nid_day,.false.,"qgl2",itau_w,tmpout)
c    -------   V GLA3 TOT
               do i=1,klon
                 do j=1,klev
                   tmpout(i,j)= SUM(qaer(i,j,4*nrad+1:5*nrad))
                 enddo
               enddo
       call histwrite_phy(nid_day,.false.,"qgl3",itau_w,tmpout)
c --------------
c ----- SATURATION ESP NUAGES
       call histwrite_phy(nid_day,.false.,"ch4sat",itau_w,satch4)
       call histwrite_phy(nid_day,.false.,"c2h6sat",itau_w,satc2h6)
       call histwrite_phy(nid_day,.false.,"c2h2sat",itau_w,satc2h2)
c --------------
c ----- RESERVOIR DE SURFACE
       call histwrite_phy(nid_day,.false.,"reserv",itau_w,reservoir)
c --------------
c ----- ECHANGE GAZ SURF/ATM (evaporation)
       call histwrite_phy(nid_day,.false.,"evapch4",itau_w,evapch4)
c --------------
c ----- PRECIPITATIONS
c       -----  CH4
       call histwrite_phy(nid_day,.false.,"prech4",
     .            itau_w,precip(1:klon,1))
c       -----  C2H6
       call histwrite_phy(nid_day,.false.,"prec2h6",
     .            itau_w,precip(1:klon,2))
c       -----  C2H2
       call histwrite_phy(nid_day,.false.,"prec2h2",
     .            itau_w,precip(1:klon,3))
c       -----  NOY
       call histwrite_phy(nid_day,.false.,"prenoy",
     .            itau_w,precip(1:klon,4))
c       -----  AER
       call histwrite_phy(nid_day,.false.,"preaer",
     .            itau_w,precip(1:klon,5))
c --------------
c ----- FLUX GLACE
c       -----  CH4
       call histwrite_phy(nid_day,.false.,"flxgl1",
     .            itau_w,flxesp_i(1:klon,1:klev,1))
c       -----  C2H6
       call histwrite_phy(nid_day,.false.,"flxgl2",
     .            itau_w,flxesp_i(1:klon,1:klev,2))
c       -----  C2H2
       call histwrite_phy(nid_day,.false.,"flxgl3",
     .            itau_w,flxesp_i(1:klon,1:klev,3))
c --------------
c ----- Source/puits GLACE
c       -----  CH4
       call histwrite_phy(nid_day,.false.,"solch4",
     .            itau_w,solesp(1:klon,1:klev,1))
c       -----  C2H6
       call histwrite_phy(nid_day,.false.,"solc2h6",
     .            itau_w,solesp(1:klon,1:klev,2))
c       -----  C2H2
       call histwrite_phy(nid_day,.false.,"solc2h2",
     .            itau_w,solesp(1:klon,1:klev,3))
c --------------
c ----- RAYON MOYEN GOUTTE
       call histwrite_phy(nid_day,.false.,"rcldbar",itau_w,rmcloud)

             endif
	  endif

c --------------
c ----- TRACEURS CHIMIQUES
	  if (nmicro.lt.nqmax) then
           DO iq=nmicro+1,nqmax
       call histwrite_phy(nid_day,.false.,tname(iq),
     .                    itau_w,tr_seri(1:klon,1:klev,iq))
           ENDDO
	  endif
         endif

      ENDIF !lev_histday.GE.2

c-------------------------------------------------------
      IF(lev_histday.GE.3) THEN

cccccccccccccccccc  Radiative transfer

c 2D

      call histwrite_phy(nid_day,.false.,"topl",itau_w,toplw)
      call histwrite_phy(nid_day,.false.,"sols",itau_w,solsw)
      call histwrite_phy(nid_day,.false.,"soll",itau_w,sollw)

c 3D

      call histwrite_phy(nid_day,.false.,"SWnet",
     .          itau_w,swnet(1:klon,1:klev))
      call histwrite_phy(nid_day,.false.,"LWnet",
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
       call histwrite_phy(nid_day,.false.,"thv"//str2,itau_w,t_tauhvd)
       enddo      ! fin boucle NSPECV 

       do k=8,NSPECI,10
         do i=1,klon
         do l=1,klev
           t_tauhvd(i,l)=TAUHID(i,klev-l+1,k)
         enddo
         enddo
         write(str2,'(i2.2)') k
       call histwrite_phy(nid_day,.false.,"thi"//str2,itau_w,t_tauhvd)
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
       call histwrite_phy(nid_day,.false.,"khv"//str2,itau_w,t_khvd)
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
       call histwrite_phy(nid_day,.false.,"khi"//str2,itau_w,t_khvd)
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
       call histwrite_phy(nid_day,.false.,"tgv"//str2,itau_w,t_tauhvd)
       enddo      ! fin boucle NSPECV 

       do k=8,NSPECI,10
         do i=1,klon
         do l=1,klev
           t_tauhvd(i,l)=TAUGID(i,klev-l+1,k)
         enddo
         enddo
         write(str2,'(i2.2)') k
       call histwrite_phy(nid_day,.false.,"tgi"//str2,itau_w,t_tauhvd)
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
       call histwrite_phy(nid_day,.false.,"kgv"//str2,itau_w,t_khvd)
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
       call histwrite_phy(nid_day,.false.,"kgi"//str2,itau_w,t_khvd)
       enddo      ! fin boucle NSPECI 

c --------------
         if (clouds.eq.1) then
c --------------
c ----- OPACITE NUAGES (ATTENTION PROXY)
         call histwrite_phy(nid_day,.false.,"tcld",itau_w,occcld)
c --------------
c ----- EXTINCTION NUAGES (ATTENTION PROXY)
           do i=1,klon
             t_kcld(i,klev)=occcld(i,klev)
     .       /(zzlev(i,klev+1)-zzlev(i,klev))
             do j=klev-1,1,-1
               t_kcld(i,j)=(occcld(i,j)-occcld(i,j+1))
     .         /(zzlev(i,j+1)-zzlev(i,j))
             enddo
           enddo
         call histwrite_phy(nid_day,.false.,"kcld",itau_w,t_kcld)
c --------------
        endif  
c --------------

      ENDIF !lev_histday.GE.3

c-------------------------------------------------------
      IF(lev_histday.GE.4) THEN

      call histwrite_phy(nid_day,.false.,"dtdyn",itau_w,d_t_dyn)
      call histwrite_phy(nid_day,.false.,"dtphy",itau_w,d_t)
c K/s
      call histwrite_phy(nid_day,.false.,"dtvdf",itau_w,d_t_vdf)
c K/s
      call histwrite_phy(nid_day,.false.,"dtajs",itau_w,d_t_ajs)
c K/s
      call histwrite_phy(nid_day,.false.,"dtswr",itau_w,heat)
c K/s
      call histwrite_phy(nid_day,.false.,"dtlwr",itau_w,-1.*cool)
c K/s      
c      call histwrite_phy(nid_day,.false.,"dtec",itau_w,d_t_ec)

      ENDIF !lev_histday.GE.4

c-------------------------------------------------------
      IF(lev_histday.GE.5) THEN

c      call histwrite_phy(nid_day,.false.,"taux",itau_w,fluxu)
c      call histwrite_phy(nid_day,.false.,"tauy",itau_w,fluxv)
c      call histwrite_phy(nid_day,.false.,"cdrm",itau_w,cdragm)
c      call histwrite_phy(nid_day,.false.,"cdrh",itau_w,cdragh)

      ENDIF !lev_histday.GE.5
c-------------------------------------------------------

      if (ok_sync) then
        call histsync(nid_day)
      endif

      ENDIF
