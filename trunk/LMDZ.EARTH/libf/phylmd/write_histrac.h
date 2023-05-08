!$Id $
!***************************************
!  ECRITURE DU FICHIER :  histrac.nc
!***************************************
  IF (ecrit_tra > 0. .AND. config_inca == 'none') THEN
     
     itau_w = itau_phy + nstep
     
     CALL histwrite_phy(nid_tra,"phis",itau_w,pphis)
     CALL histwrite_phy(nid_tra,"aire",itau_w,airephy)
     CALL histwrite_phy(nid_tra,"zmasse",itau_w,zmasse)

!TRACEURS
!----------------
     DO it=1,nbtr
        iiq=niadv(it+2)

! CONCENTRATIONS
        CALL histwrite_phy(nid_tra,tname(iiq),itau_w,tr_seri(:,:,it))

! TD LESSIVAGE       
        IF (lessivage .AND. aerosol(it)) THEN
           CALL histwrite_phy(nid_tra,"fl"//tname(iiq),itau_w,flestottr(:,:,it))
        ENDIF

! TD THERMIQUES
        IF (iflag_thermals.gt.0) THEN
           CALL histwrite_phy(nid_tra,"d_tr_th_"//tname(iiq),itau_w,d_tr_th(:,:,it))
        ENDIF

! TD CONVECTION
        IF (iflag_con.GE.2) THEN
           CALL histwrite_phy(nid_tra,"d_tr_cv_"//tname(iiq),itau_w,d_tr_cv(:,:,it))
        ENDIF

! TD COUCHE-LIMITE
        CALL histwrite_phy(nid_tra,"d_tr_cl_"//tname(iiq),itau_w,d_tr_cl(:,:,it))
     ENDDO
!---------------
!
!
! VENT (niveau 1)   
     CALL histwrite_phy(nid_tra,"pyu1",itau_w,yu1)
     CALL histwrite_phy(nid_tra,"pyv1",itau_w,yv1)
!
! TEMPERATURE DU SOL
     zx_tmp_fi2d(:)=ftsol(:,1)         
     CALL histwrite_phy(nid_tra,"ftsol1",itau_w,zx_tmp_fi2d)
     zx_tmp_fi2d(:)=ftsol(:,2)
     CALL histwrite_phy(nid_tra,"ftsol2",itau_w,zx_tmp_fi2d)
     zx_tmp_fi2d(:)=ftsol(:,3)
     CALL histwrite_phy(nid_tra,"ftsol3",itau_w,zx_tmp_fi2d)
     zx_tmp_fi2d(:)=ftsol(:,4)
     CALL histwrite_phy(nid_tra,"ftsol4",itau_w,zx_tmp_fi2d)
!      
! NATURE DU SOL
     zx_tmp_fi2d(:)=pctsrf(:,1)
     CALL histwrite_phy(nid_tra,"psrf1",itau_w,zx_tmp_fi2d)
     zx_tmp_fi2d(:)=pctsrf(:,2)
     CALL histwrite_phy(nid_tra,"psrf2",itau_w,zx_tmp_fi2d)
     zx_tmp_fi2d(:)=pctsrf(:,3)
     CALL histwrite_phy(nid_tra,"psrf3",itau_w,zx_tmp_fi2d)
     zx_tmp_fi2d(:)=pctsrf(:,4)
     CALL histwrite_phy(nid_tra,"psrf4",itau_w,zx_tmp_fi2d)
 
! DIVERS    
     CALL histwrite_phy(nid_tra,"pplay",itau_w,pplay)     
     CALL histwrite_phy(nid_tra,"T",itau_w,t_seri)     
     CALL histwrite_phy(nid_tra,"mfu",itau_w,pmfu)
     CALL histwrite_phy(nid_tra,"mfd",itau_w,pmfd)
     CALL histwrite_phy(nid_tra,"en_u",itau_w,pen_u)
     CALL histwrite_phy(nid_tra,"en_d",itau_w,pen_d)
     CALL histwrite_phy(nid_tra,"de_d",itau_w,pde_d)
     CALL histwrite_phy(nid_tra,"de_u",itau_w,pde_u)
     CALL histwrite_phy(nid_tra,"coefh",itau_w,coefh)

     IF (ok_sync) THEN
!$OMP MASTER
        CALL histsync(nid_tra)
!$OMP END MASTER
     ENDIF

  ENDIF !ecrit_tra>0. .AND. config_inca == 'none'

