!***********************************************************************
! varmuphy.h : 
! fichier contenant toutes les variables communes a la microphysique.
! Necessite microtab.h (pour la definition de nrad).
! NOTE : certaines variables restent spécifiques a seulement quelques 
!        une des routines de la microphysique.  
!***********************************************************************
!***********************************************************************
!             variables communes a tous les routines
!***********************************************************************
       real p(llm),t(llm),rho(llm),pb(llm+1),tb(llm+1),rhob(llm+1)
       real ach4(llm),aar(llm),an2(llm) !var communes a brume.F et snuages.F
       real z(llm),dz(llm),zb(llm+1),dzb(llm+1)
       real v_e(nrad),r_e(nrad),vrat_e,dr_e(nrad),dv_e(nrad)

       common/grille/z,zb,dz,dzb
       common/donnees/p,pb,t,tb,rho,rhob,ach4,aar,an2
       common/aergrid/v_e,r_e,vrat_e,dr_e,dv_e
!***********************************************************************
!             constantes communes a la microphysique.
!***********************************************************************
       real pi,                                                         &
     &      nav,rgp,kbz,                                                &
     &      rtit,g0,                                                    &
     &      mch4,mc2h6,mc2h2,mar,mn2,mair,                              &
     &      rhol,rhoi_ch4,rhoi_c2h6,rhoi_c2h2,                          &
     &      mtetach4,mtetac2h6,mtetac2h2       
   
       common/phys/                                                     &
     &      pi,                                                         &
     &      nav,rgp,kbz,                                                &
     &      rtit,g0,                                                    &
     &      mch4,mc2h6,mc2h2,mar,mn2,mair,                              &
     &      rhol,rhoi_ch4,rhoi_c2h6,rhoi_c2h2,                          &
     &      mtetach4,mtetac2h6,mtetac2h2

