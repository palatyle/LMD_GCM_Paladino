!
! $Id: write_histISCCP.h 1403 2010-07-01 09:02:53Z fairhead $
!
      IF (ok_isccp) THEN
c
       IF (MOD(itap,NINT(freq_ISCCP/dtime)).EQ.0) THEN
c
       ndex2d = 0
       ndex3d = 0
c
       itau_w = itau_phy + itap
c
       IF(type_run.EQ."ENSP".OR.type_run.EQ."CLIM") THEN
c
        DO n=1, napisccp
c
        DO k=1,kmaxm1
         zx_tmp_fi3d(1:klon, 1:lmaxm1)=fq_isccp(1:klon,k,1:lmaxm1,n)*100.
cym         CALL gr_fi_ecrit(lmaxm1,klon,iim,jjmp1,zx_tmp_fi3d,
cym     .                    zx_tmp_3d)
c
cIM: champ 3d : (lon,lat,pres) pour un tau fixe
c
      CALL histwrite_phy(nid_isccp,"cldISCCP_"//taulev(k)//verticaxe(n),
     .                  itau_w,zx_tmp_fi3d)
        ENDDO !k
c
cym        CALL gr_fi_ecrit(1, klon,iim,jjmp1,nbsunlit(1,:,n),zx_tmp_2d)
        CALL histwrite_phy(nid_isccp,"nsunlit"//verticaxe(n),itau_w,
     .                 nbsunlit(1,:,n))
c
       CALL histwrite_phy(nid_isccp,"meantaucld"//verticaxe(n),itau_w,
     .                 meantaucld(:,n))
c
        ENDDO ! n=1, napisccp
        ELSE IF(type_run.EQ."AMIP".OR.type_run.EQ."CFMI") THEN
c
        DO n=1, napisccp
c        print*,'n=',n,' write_ISCCP avant fq_isccp'
         DO k=1, kmaxm1
          DO l=1, lmaxm1
c
         IF(top_height.LE.2) THEN
          DO i=1, klon
c281008 beg
c          print*,'write_ISCCP i n nbsunlit',i,n,nbsunlit(1,i,n)
c281008 end
c
           IF(nbsunlit(1,i,n).NE.0.) THEN
            fq_is_true(i,k,l,n)=
     $      fq_isccp(i,k,l,n)*100./nbsunlit(1,i,n)
           ELSE
            fq_is_true(i,k,l,n)=0
           ENDIF
          ENDDO 
         ELSE IF(top_height.EQ.3) THEN 
          DO i=1, klon
           fq_is_true(i,k,l,n) = fq_isccp(i,k,l,n)*100.
          ENDDO
         ENDIF
cym         CALL gr_fi_ecrit(1,klon,iim,jjmp1,fq_is_true,
cym     .                    zx_tmp_2d)
c
cIM: champ 2d : (lon,lat) pour un tau et une pc fixes
c
         CALL histwrite_phy(nid_isccp,pclev(l)//taulev(k)//verticaxe(n),
     .                  itau_w,fq_is_true(:,k,l,n))
         ENDDO !l
        ENDDO !k
c
c       print*,'n=',n,' write_ISCCP avant nbsunlit'
cym        CALL gr_fi_ecrit(1, klon,iim,jjmp1,nbsunlit(1,:,n),zx_tmp_2d)
        CALL histwrite_phy(nid_isccp,"nsunlit"//verticaxe(n),
     .                 itau_w,nbsunlit(1,:,n))
c
       CALL histwrite_phy(nid_isccp,"meantaucld"//verticaxe(n),itau_w,
     .                 meantaucld(:,n))
c
        zx_tmp_fi2d(1:klon)=REAL(seed(1:klon,n))
c
c       print*,'n=',n,' write_ISCCP avant seed'
cym        CALL gr_fi_ecrit(1, klon,iim,jjmp1,zx_tmp_fi2d,zx_tmp_2d)
        CALL histwrite_phy(nid_isccp,"seed"//verticaxe(n),
     .                 itau_w,zx_tmp_fi2d)
c
c 9types de nuages ISCCP-D2
c fq_isccp(1:klon,k,l,n)*100. <=> pc_tau(k)_pclev(l)
        DO i=1, klon
         zx_tmp_fi2d(i)=
     $ (fq_is_true(i,1,1,n)+ fq_is_true(i,2,1,n)+ fq_is_true(i,3,1,n) +
     $  fq_is_true(i,1,2,n)+ fq_is_true(i,2,2,n)+ fq_is_true(i,3,2,n) +
     $  fq_is_true(i,1,3,n)+ fq_is_true(i,2,3,n)+ fq_is_true(i,3,3,n) )
        ENDDO
cym       CALL gr_fi_ecrit(1, klon,iim,jjmp1,zx_tmp_fi2d,zx_tmp_2d)
        CALL histwrite_phy(nid_isccp,"cirr",itau_w,zx_tmp_fi2d)
c
        DO i=1, klon
         zx_tmp_fi2d(i)=
     $  (fq_is_true(i,4,1,n)+ fq_is_true(i,5,1,n) +
     $   fq_is_true(i,4,2,n)+ fq_is_true(i,5,2,n) +
     $   fq_is_true(i,4,3,n)+ fq_is_true(i,5,3,n) )
        ENDDO
cym	CALL gr_fi_ecrit(1, klon,iim,jjmp1,zx_tmp_fi2d,zx_tmp_2d)
        CALL histwrite_phy(nid_isccp,"cist",itau_w,zx_tmp_fi2d)
c
        DO i=1, klon
         zx_tmp_fi2d(i)=
     $  (fq_is_true(i,6,1,n)+ fq_is_true(i,7,1,n) +
     $   fq_is_true(i,6,2,n)+ fq_is_true(i,7,2,n) +
     $   fq_is_true(i,6,3,n)+ fq_is_true(i,7,3,n) )
        ENDDO
cym	CALL gr_fi_ecrit(1, klon,iim,jjmp1,zx_tmp_fi2d,zx_tmp_2d)
        CALL histwrite_phy(nid_isccp,"deep",itau_w,zx_tmp_fi2d)
c
        DO i=1, klon
         zx_tmp_fi2d(i)=
     $ (fq_is_true(i,1,4,n)+ fq_is_true(i,2,4,n)+ fq_is_true(i,3,4,n) +
     $  fq_is_true(i,1,5,n)+ fq_is_true(i,2,5,n)+ fq_is_true(i,3,5,n) )
        ENDDO
cym	CALL gr_fi_ecrit(1, klon,iim,jjmp1,zx_tmp_fi2d,zx_tmp_2d)
        CALL histwrite_phy(nid_isccp,"alcu",itau_w,zx_tmp_fi2d)
c
        DO i=1, klon
         zx_tmp_fi2d(i)=
     $  (fq_is_true(i,4,4,n)+ fq_is_true(i,5,4,n) +
     $   fq_is_true(i,4,5,n)+ fq_is_true(i,5,5,n) )
        ENDDO
cym	CALL gr_fi_ecrit(1, klon,iim,jjmp1,zx_tmp_fi2d,zx_tmp_2d)
        CALL histwrite_phy(nid_isccp,"alst",itau_w,zx_tmp_fi2d)
c
        DO i=1, klon
         zx_tmp_fi2d(i)=
     $  (fq_is_true(i,6,4,n)+ fq_is_true(i,7,4,n) +
     $   fq_is_true(i,6,5,n)+ fq_is_true(i,7,5,n) )
        ENDDO
cym	CALL gr_fi_ecrit(1, klon,iim,jjmp1,zx_tmp_fi2d,zx_tmp_2d)
        CALL histwrite_phy(nid_isccp,"nist",itau_w,zx_tmp_fi2d)
c
        DO i=1, klon
         zx_tmp_fi2d(i)=
     $ (fq_is_true(i,1,6,n)+ fq_is_true(i,2,6,n)+ fq_is_true(i,3,6,n) +
     $  fq_is_true(i,1,7,n)+ fq_is_true(i,2,7,n)+ fq_is_true(i,3,7,n) )
        ENDDO
cym	CALL gr_fi_ecrit(1, klon,iim,jjmp1,zx_tmp_fi2d,zx_tmp_2d)
        CALL histwrite_phy(nid_isccp,"cumu",itau_w,zx_tmp_fi2d)
c
        DO i=1, klon
         zx_tmp_fi2d(i)=
     $  (fq_is_true(i,4,6,n)+ fq_is_true(i,5,6,n) +
     $   fq_is_true(i,4,7,n)+ fq_is_true(i,5,7,n) )
        ENDDO
cym	CALL gr_fi_ecrit(1, klon,iim,jjmp1,zx_tmp_fi2d,zx_tmp_2d)
        CALL histwrite_phy(nid_isccp,"stcu",itau_w,zx_tmp_fi2d)
c
        DO i=1, klon
         zx_tmp_fi2d(i)=
     $  (fq_is_true(i,6,6,n)+ fq_is_true(i,7,6,n) +
     $   fq_is_true(i,6,7,n)+ fq_is_true(i,7,7,n) )
        ENDDO
cym	CALL gr_fi_ecrit(1, klon,iim,jjmp1,zx_tmp_fi2d,zx_tmp_2d)
        CALL histwrite_phy(nid_isccp,"stra",itau_w,zx_tmp_fi2d)
c
c 3_tau_nuages x 3_levels
c fq_is_true(1:klon,k,l,n)*100. <=> pc_tau(k)_pclev(l)
        DO i=1, klon
         cld_fi3d(i,1)= 
     $ (fq_is_true(i,1,1,n)+ fq_is_true(i,2,1,n)+ fq_is_true(i,3,1,n) +
     $  fq_is_true(i,1,2,n)+ fq_is_true(i,2,2,n)+ fq_is_true(i,3,2,n) +
     $  fq_is_true(i,1,3,n)+ fq_is_true(i,2,3,n)+ fq_is_true(i,3,3,n) )
	 cld_fi3d(i,2)=
     $ (fq_is_true(i,1,4,n)+ fq_is_true(i,2,4,n)+ fq_is_true(i,3,4,n) +
     $  fq_is_true(i,1,5,n)+ fq_is_true(i,2,5,n)+ fq_is_true(i,3,5,n) )
         cld_fi3d(i,3)=
     $ (fq_is_true(i,1,6,n)+ fq_is_true(i,2,6,n)+ fq_is_true(i,3,6,n) +
     $  fq_is_true(i,1,7,n)+ fq_is_true(i,2,7,n)+ fq_is_true(i,3,7,n) )
        ENDDO   
cym        CALL gr_fi_ecrit(lmax3,klon,iim,jjmp1,cld_fi3d,cld_3d)
        CALL histwrite_phy(nid_isccp,"thin",itau_w,cld_fi3d)
c
        DO i=1, klon
	 cld_fi3d(i,1)=
     $   (fq_is_true(i,4,1,n)+ fq_is_true(i,5,1,n) +
     $    fq_is_true(i,4,2,n)+ fq_is_true(i,5,2,n) +
     $    fq_is_true(i,4,3,n)+ fq_is_true(i,5,3,n) )
	 cld_fi3d(i,2)=
     $   (fq_is_true(i,4,4,n)+ fq_is_true(i,5,4,n) +
     $    fq_is_true(i,4,5,n)+ fq_is_true(i,5,5,n) )
	 cld_fi3d(i,3)=
     $   (fq_is_true(i,4,6,n)+ fq_is_true(i,5,6,n) +
     $    fq_is_true(i,4,7,n)+ fq_is_true(i,5,7,n) )
	ENDDO   
cym       CALL gr_fi_ecrit(lmax3, klon,iim,jjmp1,cld_fi3d,cld_3d)
        CALL histwrite_phy(nid_isccp,"mid",itau_w,cld_fi3d)
c
        DO i=1, klon
	 cld_fi3d(i,1)=
     $   (fq_is_true(i,6,1,n)+ fq_is_true(i,7,1,n) +
     $    fq_is_true(i,6,2,n)+ fq_is_true(i,7,2,n) +
     $    fq_is_true(i,6,3,n)+ fq_is_true(i,7,3,n) )
         cld_fi3d(i,2)=
     $   (fq_is_true(i,6,4,n)+ fq_is_true(i,7,4,n) +
     $    fq_is_true(i,6,5,n)+ fq_is_true(i,7,5,n) )
	 cld_fi3d(i,3)=
     $   (fq_is_true(i,6,6,n)+ fq_is_true(i,7,6,n) +
     $    fq_is_true(i,6,7,n)+ fq_is_true(i,7,7,n) )
        ENDDO   
cym       CALL gr_fi_ecrit(lmax3, klon,iim,jjmp1,cld_fi3d,cld_3d)
        CALL histwrite_phy(nid_isccp,"thick",itau_w,cld_fi3d)
c
        ENDDO ! n=1, napisccp
c
       ENDIF
c
       if (ok_sync) then
c$OMP MASTER
        call histsync(nid_isccp)
c$OMP END MASTER       
       endif

       ENDIF !(MOD(itap,NINT(freq_ISCCP/dtime)).EQ.0) THEN

      ENDIF !ok_isccp
