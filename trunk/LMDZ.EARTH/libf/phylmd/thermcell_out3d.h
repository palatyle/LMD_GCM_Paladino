!       if (sorties) then

!      print*,'16 OK convect8'
         call wrgradsfi(1,nlay,pt(igout,1:klev),'pt        ','pt        ')
         call wrgradsfi(1,nlay,fraca(igout,1:klev),'fraca     ','fraca     ')
         call wrgradsfi(1,nlay,zh(igout,1:klev),'zh        ','zh        ')
         call wrgradsfi(1,nlay,zha(igout,1:klev),'zha        ','zha        ')
         call wrgradsfi(1,nlay,zua(igout,1:klev),'zua        ','zua        ')
         call wrgradsfi(1,nlay,zva(igout,1:klev),'zva        ','zva        ')
         call wrgradsfi(1,nlay,zu(igout,1:klev),'zu        ','zu        ')
         call wrgradsfi(1,nlay,zv(igout,1:klev),'zv        ','zv        ')
         call wrgradsfi(1,nlay,zo(igout,1:klev),'zo        ','zo        ')
         call wrgradsfi(1,1,zmax(igout),'zmax      ','zmax      ')
!         call wrgradsfi(1,nlay,zdhadj(igout,1:klev),'zdhadj    ','zdhadj    ')
         call wrgradsfi(1,nlay,pduadj(igout,1:klev),'pduadj    ','pduadj    ')
         call wrgradsfi(1,nlay,pdvadj(igout,1:klev),'pdvadj    ','pdvadj    ')
         call wrgradsfi(1,nlay,pdoadj(igout,1:klev),'pdoadj    ','pdoadj    ')
         call wrgradsfi(1,nlay,entr(igout,1:klev),'entr      ','entr      ')
         call wrgradsfi(1,nlay,detr(igout,1:klev),'detr      ','detr      ')
         call wrgradsfi(1,nlay,fm(igout,1:klev),'fm        ','fm        ')
         call wrgradsfi(1,nlay,zw2(igout,1:klev),'zw2       ','zw2       ')
         call wrgradsfi(1,nlay,zw_est(igout,1:klev),'w_est      ','w_est      ')
!on sort les moments
         call wrgradsfi(1,nlay,thetath2(igout,1:klev),'zh2       ','zh2       ')
         call wrgradsfi(1,nlay,wth2(igout,1:klev),'w2       ','w2       ')
         call wrgradsfi(1,nlay,wth3(igout,1:klev),'w3       ','w3       ')
         call wrgradsfi(1,nlay,q2(igout,1:klev),'q2       ','q2       ')
!
!
         call wrgradsfi(1,nlay,wthl(igout,1:klev),'wthl       ','wthl       ')
         call wrgradsfi(1,nlay,wthv(igout,1:klev),'wthv       ','wthv       ')
         call wrgradsfi(1,nlay,wq(igout,1:klev),'wq       ','wq       ')
         
         call wrgradsfi(1,nlay,ztva(igout,1:klev),'ztva      ','ztva      ')
         call wrgradsfi(1,nlay,ztv(igout,1:klev),'ztv       ','ztv       ')

         call wrgradsfi(1,nlay,zo(igout,1:klev),'zo        ','zo        ')
         call wrgradsfi(1,nlay,zoa(igout,1:klev),'zoa        ','zoa        ')

!nouveaux diagnostiques
         call wrgradsfi(1,nlay,zthl(igout,1:klev),'zthl        ','zthl        ')
         call wrgradsfi(1,nlay,zta(igout,1:klev),'zta        ','zta        ')
         call wrgradsfi(1,nlay,zl(igout,1:klev),'zl        ','zl        ')
         call wrgradsfi(1,nlay,zdthladj(igout,1:klev),'zdthladj    ',  &
     &        'zdthladj    ')
         call wrgradsfi(1,nlay,ztla(igout,1:klev),'ztla      ','ztla      ')
         call wrgradsfi(1,nlay,zqta(igout,1:klev),'zqta      ','zqta      ')
         call wrgradsfi(1,nlay,zqla(igout,1:klev),'zqla      ','zqla      ')
         call wrgradsfi(1,nlay,deltaz(igout,1:klev),'deltaz      ','deltaz      ')
!nouveaux diagnostiques
      call wrgradsfi(1,nlay,entr_star  (igout,1:klev),'entr_star   ','entr_star   ')
      call wrgradsfi(1,nlay,detr_star  (igout,1:klev),'detr_star   ','detr_star   ')     
      call wrgradsfi(1,nlay,f_star    (igout,1:klev),'f_star   ','f_star   ')
      call wrgradsfi(1,nlay,zqsat    (igout,1:klev),'zqsat   ','zqsat   ')
      call wrgradsfi(1,nlay,zqsatth    (igout,1:klev),'qsath   ','qsath   ')
      call wrgradsfi(1,nlay,alim_star    (igout,1:klev),'alim_star   ','alim_star   ')
!      call wrgradsfi(1,nlay,alim    (igout,1:klev),'alim   ','alim   ')
      call wrgradsfi(1,1,f(igout),'f      ','f      ')
      call wrgradsfi(1,1,alim_star_tot(igout),'a_s_t      ','a_s_t      ')
      call wrgradsfi(1,1,zmax(igout),'zmax      ','zmax      ')
      call wrgradsfi(1,1,zmax_sec(igout),'z_sec      ','z_sec      ')
      call wrgradsfi(1,1,zmix(igout),'zmix      ','zmix      ') 
!      call wrgradsfi(1,1,nivcon(igout),'nivcon      ','nivcon      ')
      call wrgradsfi(1,1,zcon(igout),'zcon      ','zcon      ')
      call wrgradsfi(1,1,zcon2(igout),'zcon2      ','zcon2      ')
      zsortie1d(:)=lmax(:)
      call wrgradsfi(1,1,zsortie1d(igout),'lmax      ','lmax      ')
      call wrgradsfi(1,1,wmax(igout),'wmax      ','wmax      ')
      call wrgradsfi(1,1,wmax_sec(igout),'w_sec      ','w_sec      ')
!      zsortie1d(:)=lmix(:)
!      call wrgradsfi(1,1,zsortie1d(igout),'lmix      ','lmix      ')
!      zsortie1d(:)=lentr(:)
!      call wrgradsfi(1,1,zsortie1d(igout),'lentr      ','lentr     ')

      print*,'Fin des wrgradsfi'
