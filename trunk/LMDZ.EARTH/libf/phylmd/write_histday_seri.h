c
c $Header$
c
      IF (is_sequential) THEN
      
      IF (type_run.EQ."AMIP") THEN
c
      ndex2d = 0
      itau_w = itau_phy + itap
c
c Champs 2D:
c
      pi = ACOS(-1.)
      pir = 4.0*ATAN(1.0) / 180.0
c
      DO i=1, klon
       zx_tmp_fi2d(i)=(topsw(i)-toplw(i))
      ENDDO
c
      ok_msk=.FALSE.
      msk(1:klon)=pctsrf(1:klon,is_ter)
      CALL moyglo_pondaire(klon, zx_tmp_fi2d, airephy, 
     .     ok_msk, msk, moyglo)
      zx_tmp_fi2d(1:klon)=moyglo
c
      CALL gr_fi_ecrit(1, klon,iim,jjmp1, zx_tmp_fi2d,zx_tmp_2d)
      CALL histwrite(nid_day_seri,"bilTOA",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      ok_msk=.FALSE.
      CALL moyglo_pondaire(klon, bils, airephy, 
     .     ok_msk, msk, moyglo)
      zx_tmp_fi2d(1:klon)=moyglo
c
      CALL gr_fi_ecrit(1, klon,iim,jjmp1, zx_tmp_fi2d,zx_tmp_2d)
      CALL histwrite(nid_day_seri,"bils",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      DO k=1, klev
      DO i=1, klon
cIM 080904    zx_tmp_fi3d(i,k)=u(i,k)**2+v(i,k)**2
       zx_tmp_fi3d(i,k)=(u(i,k)**2+v(i,k)**2)/2.
      ENDDO
      ENDDO
c
      CALL moyglo_pondaima(klon, klev, zx_tmp_fi3d, 
     .     airephy, paprs, moyglo)
      zx_tmp_fi2d(1:klon)=moyglo
c
      CALL gr_fi_ecrit(1, klon,iim,jjmp1, zx_tmp_fi2d,zx_tmp_2d)
      CALL histwrite(nid_day_seri,"ecin",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d) 
c
cIM 151004 BEG
      IF(1.EQ.0) THEN
c
      DO k=1, klev
      DO i=1, klon
       zx_tmp_fi3d(i,k)=u_seri(i,k)*RA*cos(pir* rlat(i))
      ENDDO
      ENDDO
c
      CALL moyglo_pondaima(klon, klev, zx_tmp_fi3d, 
     .     airephy, paprs, moyglo)
      zx_tmp_fi2d(1:klon)=moyglo
c
      CALL gr_fi_ecrit(1, klon,iim,jjmp1, zx_tmp_fi2d,zx_tmp_2d)
      CALL histwrite(nid_day_seri,"momang",itau_w,zx_tmp_2d,
     .               iim*jjmp1,ndex2d)
c
c friction torque
c
      DO i=1, klon
       zx_tmp_fi2d(i)=zxfluxu(i,1)*RA* cos(pir* rlat(i))
      ENDDO
c
      ok_msk=.FALSE.
      CALL moyglo_pondaire(klon, zx_tmp_fi2d, airephy, 
     .     ok_msk, msk, moyglo)
      zx_tmp_fi2d(1:klon)=moyglo
c
      CALL gr_fi_ecrit(1, klon,iim,jjmp1, zx_tmp_fi2d,zx_tmp_2d)
      CALL histwrite(nid_day_seri,"frictor",itau_w,zx_tmp_2d,
     .               iim*jjmp1,ndex2d)
c
c mountain torque
c
cIM 190504 BEG
      CALL gr_fi_dyn(1,klon,iim+1,jjm+1,airephy,airedyn)
      CALL gr_fi_dyn(klev+1,klon,iim+1,jjm+1,paprs,padyn)
      CALL gr_fi_dyn(1,klon,iim+1,jjm+1,rlat,rlatdyn)
      mountor=0.
      airetot=0.
      DO j = 1, jjmp1
       DO i = 1, iim+1
        ij=i+(iim+1)*(j-1)
        zx_tmp(ij)=0.
        DO k = 1, klev
         zx_tmp(ij)=zx_tmp(ij)+dudyn(i,j,k)*airedyn(i,j)*
     $              (padyn(i,j,k+1)-padyn(i,j,k))/RG
         airetot=airetot+airedyn(i,j)
        ENDDO
cIM 190504 mountor=mountor+zx_tmp(ij)*airedyn(i,j)*RA*
        mountor=mountor+zx_tmp(ij)*RA*
     $           cos(pir* rlatdyn(i,j))
       ENDDO
      ENDDO
cIM 151004 BEG
      IF(itap.EQ.1) PRINT*,'airetot=',airetot,airetot/klev
cIM 151004 END
cIM 190504      mountor=mountor/(airetot*airetot)
      mountor=mountor/airetot
c
cIM 190504 END
      zx_tmp_2d(1:iim,1:jjmp1)=mountor
      CALL histwrite(nid_day_seri,"mountor",itau_w,zx_tmp_2d,
     .               iim*jjmp1,ndex2d)
c
      ENDIF !(1.EQ.0) THEN
c
c
      CALL gr_fi_dyn(1,klon,iim+1,jjm+1,airephy,airedyn)
      CALL gr_fi_ecrit(1,klon,iim,jjmp1,airephy,zx_tmp_2d)
      airetot=0.
c     DO j = 1, jjmp1
c      DO i = 1, iim+1
c       ij=i+(iim+1)*(j-1)
c       DO k = 1, klev
c        airetot=airetot+airedyn(i,j)
c        airetot=airetot+airedyn(i,j)
c       ENDDO !k
c      ENDDO !i
c     ENDDO !j
c
      DO i=1, klon
       airetot=airetot+airephy(i)
      ENDDO
c     IF(itap.EQ.1) PRINT*,'airetotphy=',airetot
c
      airetot=0.
      DO j=1, jjmp1
       DO i=1, iim
        airetot=airetot+zx_tmp_2d(i,j)
       ENDDO
      ENDDO
c
c     IF(itap.EQ.1) PRINT*,'airetotij=',airetot,
c    $ '4piR2',4.*pi*RA*RA
c
      zx_tmp_fi2d(1:klon)=aam/airetot
      CALL gr_fi_ecrit(1,klon,iim,jjmp1,zx_tmp_fi2d,zx_tmp_2d)
      CALL histwrite(nid_day_seri,"momang",itau_w,zx_tmp_2d,
     .               iim*jjmp1,ndex2d)
c
      zx_tmp_fi2d(1:klon)=torsfc/airetot
      CALL gr_fi_ecrit(1,klon,iim,jjmp1,zx_tmp_fi2d,zx_tmp_2d)
      CALL histwrite(nid_day_seri,"torsfc",itau_w,zx_tmp_2d,
     .               iim*jjmp1,ndex2d)
c
cIM 151004 END
c
      CALL moyglo_pondmass(klon, klev, t_seri,
     .     airephy, paprs, moyglo)
      zx_tmp_fi2d(1:klon)=moyglo
c
      CALL gr_fi_ecrit(1,klon,iim,jjmp1,zx_tmp_fi2d,zx_tmp_2d)
      CALL histwrite(nid_day_seri,"tamv",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      ok_msk=.FALSE.
      CALL moyglo_pondaire(klon, paprs(:,1), airephy, 
     .     ok_msk, msk, moyglo)
      zx_tmp_fi2d(1:klon)=moyglo
c
      CALL gr_fi_ecrit(1, klon,iim,jjmp1, zx_tmp_fi2d,zx_tmp_2d)
      CALL histwrite(nid_day_seri,"psol",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      ok_msk=.FALSE.
      CALL moyglo_pondaire(klon, evap, airephy, 
     .     ok_msk, msk, moyglo)
      zx_tmp_fi2d(1:klon)=moyglo
c
      CALL gr_fi_ecrit(1, klon,iim,jjmp1, zx_tmp_fi2d,zx_tmp_2d)
      CALL histwrite(nid_day_seri,"evap",itau_w,
     .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
c     DO i=1, klon
c      zx_tmp_fi2d(i)=SnowFrac(i,is_ter)
c     ENDDO
c
c     ok_msk=.TRUE.
c     msk(1:klon)=pctsrf(1:klon,is_ter)
c     CALL moyglo_pondaire(klon, zx_tmp_fi2d, airephy, 
c    .                     ok_msk, msk, moyglo)
c     zx_tmp_fi2d(1:klon)=moyglo
c
c     CALL gr_fi_ecrit(1, klon,iim,jjmp1,zx_tmp_fi2d,zx_tmp_2d)
c     CALL histwrite(nid_day_seri,"SnowFrac",
c    .               itau_w,zx_tmp_2d,iim*jjmp1,ndex2d) 
c
c     DO i=1, klon
cIM 080904    zx_tmp_fi2d(i)=zsnow_mass(i)/330.*rowl
c      zx_tmp_fi2d(i)=zsnow_mass(i)
c     ENDDO
c
cIM 140904   ok_msk=.FALSE.
c     ok_msk=.TRUE.
c     msk(1:klon)=pctsrf(1:klon,is_ter)
c     CALL moyglo_pondaire(klon, zx_tmp_fi2d, airephy, 
c    .     ok_msk, msk, moyglo)
c     zx_tmp_fi2d(1:klon)=moyglo
c
c     CALL gr_fi_ecrit(1, klon,iim,jjmp1,zx_tmp_fi2d,zx_tmp_2d)
c     CALL histwrite(nid_day_seri,"snow_depth",itau_w,
c    .               zx_tmp_2d,iim*jjmp1,ndex2d)
c
      DO i=1, klon
       zx_tmp_fi2d(i)=ftsol(i,is_oce)
      ENDDO
c
      ok_msk=.TRUE.
      msk(1:klon)=pctsrf(1:klon,is_oce)
      CALL moyglo_pondaire(klon, zx_tmp_fi2d, airephy, 
     .     ok_msk, msk, moyglo)
      zx_tmp_fi2d(1:klon)=moyglo
c
      CALL gr_fi_ecrit(1, klon,iim,jjmp1, zx_tmp_fi2d, zx_tmp_2d)
      CALL histwrite(nid_day_seri,"tsol_"//clnsurf(is_oce),
     $               itau_w,zx_tmp_2d,iim*jjmp1,ndex2d) 
c
c=================================================================
c=================================================================
c=================================================================
c
      if (ok_sync) then
        call histsync(nid_day_seri)
      endif
c
      ENDIF !fin test sur type_run.EQ."AMIP"
      
      ENDIF  ! mono_cpu
