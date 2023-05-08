!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/ocean_slab_mod.F90,v 1.3 2008-02-04 16:24:28 fairhead Exp $
!
MODULE surf_heat_transp_mod

IMPLICIT NONE  

  ! Variables copied over from dyn3d dynamics:
  REAL,SAVE,ALLOCATABLE :: fext(:)
  REAL,SAVE,ALLOCATABLE :: unsairez(:)
  REAL,SAVE,ALLOCATABLE :: unsaire(:)
  REAL,SAVE,ALLOCATABLE :: cu(:)
  REAL,SAVE,ALLOCATABLE :: cv(:)
  REAL,SAVE,ALLOCATABLE :: cuvsurcv(:)
  REAL,SAVE,ALLOCATABLE :: cvusurcu(:)
  REAL,SAVE,ALLOCATABLE :: aire(:)
  REAL,SAVE :: apoln
  REAL,SAVE :: apols
  REAL,SAVE,ALLOCATABLE :: aireu(:)
  REAL,SAVE,ALLOCATABLE :: airev(:) 

  LOGICAL,SAVE :: alpha_var
  LOGICAL,SAVE :: slab_upstream
  REAL,SAVE,ALLOCATABLE :: zmasqu(:)
  REAL,SAVE,ALLOCATABLE :: zmasqv(:)
  REAL,SAVE,ALLOCATABLE :: unsfv(:)
  REAL,SAVE,ALLOCATABLE :: unsev(:)
  REAL,SAVE,ALLOCATABLE :: unsfu(:)
  REAL,SAVE,ALLOCATABLE :: unseu(:)

  ! Routines usable only by routines within this module:
  PRIVATE :: gr_fi_dyn, gr_dyn_fi
  ! Routines from dyn3d, valid on global dynamics grid only:
  PRIVATE :: grad,diverg,gr_v_scal,gr_scal_v,gr_scal_u
  
CONTAINS
  
  SUBROUTINE ini_surf_heat_transp(ip1jm,ip1jmp1,unsairez_,fext_,unsaire_,&
                                  cu_,cuvsurcv_,cv_,cvusurcu_, &
                                  aire_,apoln_,apols_, &
                                  aireu_,airev_)
    USE mod_grid_phy_lmdz, only: nbp_lon, nbp_lat
    IMPLICIT NONE
    ! Transfer some variables from dyn3d dynamics
    INTEGER,INTENT(IN) :: ip1jm
    INTEGER,INTENT(IN) :: ip1jmp1
    REAL,INTENT(IN) :: unsairez_(ip1jm)
    REAL,INTENT(IN) :: fext_(ip1jm)
    REAL,INTENT(IN) :: unsaire_(ip1jmp1)
    REAL,INTENT(IN) :: cu_(ip1jmp1)
    REAL,INTENT(IN) :: cuvsurcv_(ip1jm)
    REAL,INTENT(IN) :: cv_(ip1jm)
    REAL,INTENT(IN) :: cvusurcu_(ip1jmp1)
    REAL,INTENT(IN) :: aire_(ip1jmp1)
    REAL,INTENT(IN) :: apoln_
    REAL,INTENT(IN) :: apols_
    REAL,INTENT(IN) :: aireu_(ip1jmp1)
    REAL,INTENT(IN) :: airev_(ip1jm)
    
    ! Sanity check
    if ((ip1jm.ne.((nbp_lon+1)*(nbp_lat-1))).or. &
        (ip1jmp1.ne.((nbp_lon+1)*nbp_lat))) then
      write(*,*) "ini_surf_heat_transp Error: wrong array sizes"
      stop
    endif
    
    allocate(unsairez(ip1jm))
    unsairez(:)=unsairez_(:)
    allocate(fext(ip1jm))
    fext(:)=fext_(:)
    allocate(unsaire(ip1jmp1))
    unsaire(:)=unsaire_(:)
    allocate(cu(ip1jmp1))
    cu(:)=cu_(:)
    allocate(cuvsurcv(ip1jm))
    cuvsurcv(:)=cuvsurcv_(:)
    allocate(cv(ip1jm))
    cv(:)=cv_(:)
    allocate(cvusurcu(ip1jmp1))
    cvusurcu(:)=cvusurcu_(:)
    allocate(aire(ip1jmp1))
    aire(:)=aire_(:)
    apoln=apoln_
    apols=apols_
    allocate(aireu(ip1jmp1))
    aireu(:)=aireu_(:)
    allocate(airev(ip1jm))
    airev(:)=airev_(:)
    
  END SUBROUTINE ini_surf_heat_transp
  
  SUBROUTINE ini_surf_heat_transp_mod
    USE mod_grid_phy_lmdz, only: nbp_lon, nbp_lat
    IMPLICIT NONE
    INTEGER :: ip1jm, ip1jmp1
    
    ip1jm=(nbp_lon+1)*(nbp_lat-1)
    ip1jmp1=(nbp_lon+1)*nbp_lat
  
    allocate(zmasqu(ip1jmp1))
    allocate(zmasqv(ip1jm))
    allocate(unsfv(ip1jm))
    allocate(unsev(ip1jm))
    allocate(unsfu(ip1jmp1))
    allocate(unseu(ip1jmp1))

  END SUBROUTINE ini_surf_heat_transp_mod

  SUBROUTINE divgrad_phy(ngrid,nlevs,temp,delta)

      USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat
      IMPLICIT NONE

      INTEGER,INTENT(IN) :: ngrid, nlevs
      REAL,INTENT(IN) :: temp(ngrid,nlevs)
      REAL,INTENT(OUT) :: delta(ngrid,nlevs)
      REAL delta_2d((nbp_lon+1)*nbp_lat,nlevs)
      INTEGER :: ll
      REAL ghx((nbp_lon+1)*nbp_lat,nlevs)
      REAL ghy((nbp_lon+1)*(nbp_lat-1),nlevs)
      INTEGER :: iip1,jjp1
      
      iip1=nbp_lon+1
      jjp1=nbp_lat

      CALL gr_fi_dyn(nlevs,ngrid,iip1,jjp1,temp,delta_2d)
      CALL grad(nlevs,delta_2d,ghx,ghy)
      DO ll=1,nlevs
          ghx(:,ll)=ghx(:,ll)*zmasqu
! pas de diffusion zonale  
!          ghx(:,ll)=0.
          ghy(:,ll)=ghy(:,ll)*zmasqv
      END DO
      CALL diverg(nlevs,ghx,ghy,delta_2d)
      CALL gr_dyn_fi(nlevs,iip1,jjp1,ngrid,delta_2d,delta)


  END SUBROUTINE divgrad_phy



  SUBROUTINE init_masquv(ngrid,zmasq)

      USE mod_grid_phy_lmdz, only: nbp_lon, nbp_lat
      IMPLICIT NONE


      INTEGER,INTENT(IN) :: ngrid
      REAL zmasq(ngrid)
      REAL zmasq_2d((nbp_lon+1)*nbp_lat)
      REAL ff((nbp_lon+1)*(nbp_lat-1))
      REAL eps
      INTEGER i
      INTEGER :: iim,iip1,jjp1,ip1jm,ip1jmp1
      
      iim=nbp_lon
      iip1=nbp_lon+1
      jjp1=nbp_lat
      ip1jm=(nbp_lon+1)*(nbp_lat-1)
      ip1jmp1=(nbp_lon+1)*nbp_lat

! Masques u,v
      zmasqu=1.
      zmasqv=1.

      CALL gr_fi_dyn(1,ngrid,iip1,jjp1,zmasq,zmasq_2d)

      DO i=1,ip1jmp1-1
        IF (zmasq_2d(i).GT.1e-5 .OR. zmasq_2d(i+1).GT.1e-5) THEN
                zmasqu(i)=0.
        ENDIF
      END DO
      DO i=iip1,ip1jmp1,iip1
        zmasqu(i)=zmasqu(i-iim)
      END DO
      DO i=1,ip1jm
        IF (zmasq_2d(i).GT.1e-5 .OR. zmasq_2d(i+iip1).GT.1e-5) THEN  
                zmasqv(i)=0.
        END IF
      END DO


! Coriolis (pour Ekman transp.)
      eps=1e-5
!      CALL getin('slab_eps',eps)
!      print *,'epsilon=',eps
      ff=fext*unsairez       
      DO i=1,ip1jm
        unsev(i)=eps/(ff(i)*ff(i)+eps**2)
        unsfv(i)=ff(i)/(ff(i)*ff(i)+eps**2)
      ENDDO
      CALL gr_v_scal(1,unsfv,unsfu)
      CALL gr_v_scal(1,unsev,unseu)
! Alpha variable?
!      alpha_var=.FALSE.
!      CALL getin('slab_alphav',alpha_var)


      
  END SUBROUTINE init_masquv



  SUBROUTINE slab_ekman2(ngrid,tx_phy,ty_phy,ts_phy,dt_phy)
  
      USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat
      USE slab_ice_h, ONLY: noceanmx

      IMPLICIT NONE
      
      INTEGER,INTENT(IN) :: ngrid
      INTEGER ij
      REAL txv((nbp_lon+1)*(nbp_lat-1)),fluxm((nbp_lon+1)*(nbp_lat-1))
      REAL tyv((nbp_lon+1)*(nbp_lat-1))
      REAL fluxtm((nbp_lon+1)*(nbp_lat-1),noceanmx)
      REAL fluxtz((nbp_lon+1)*nbp_lat,noceanmx)
      REAL tyu((nbp_lon+1)*nbp_lat),txu((nbp_lon+1)*nbp_lat)
      REAL fluxz((nbp_lon+1)*nbp_lat),fluxv((nbp_lon+1)*nbp_lat)
      REAL dt((nbp_lon+1)*nbp_lat,noceanmx),ts((nbp_lon+1)*nbp_lat,noceanmx)
      REAL tx_phy(ngrid),ty_phy(ngrid)
      REAL dt_phy(ngrid,noceanmx),ts_phy(ngrid,noceanmx)

      INTEGER iim,iip1,iip2,jjp1,ip1jm,ip1jmi1,ip1jmp1

      iim=nbp_lon
      iip1=nbp_lon+1
      iip2=nbp_lon+2
      jjp1=nbp_lat
      ip1jm=(nbp_lon+1)*(nbp_lat-1) ! = iip1*jjm
      ip1jmi1=(nbp_lon+1)*(nbp_lat-1)-(nbp_lon+1) ! = ip1jm - iip1
      ip1jmp1=(nbp_lon+1)*nbp_lat ! = iip1*jjp1

! Convert taux,y to 2D  scalar grid
      ! north and south poles tx,ty no meaning
      tx_phy(1)=0.
      tx_phy(ngrid)=0.
      ty_phy(1)=0.
      ty_phy(ngrid)=0.
      CALL gr_fi_dyn(1,ngrid,iip1,jjp1,tx_phy,txu)
      CALL gr_fi_dyn(1,ngrid,iip1,jjp1,ty_phy,tyu)
          
! Divide taux,y by f or eps, and convert to 2D u,v grids
! (Arakawa C grid)
      CALL gr_scal_v(1,txu,txv) ! wind stress at v points
      CALL gr_scal_v(1,tyu,tyv)
      fluxm=tyv*unsev-txv*unsfv
!      fluxm=-txv*unsfv
! Zonal flux
      CALL gr_scal_u(1,txu,txu) ! wind stress at u points
      CALL gr_scal_u(1,tyu,tyu)
      fluxz=tyu*unsfu+txu*unseu
!      fluxz=tyu*unsfu
            
! Convert temperature to 2D grid
      CALL gr_fi_dyn(2,ngrid,iip1,jjp1,ts_phy,ts)
       
! Flux de masse
      fluxm=fluxm*cv*cuvsurcv*zmasqv
      fluxz=fluxz*cu*cvusurcu*zmasqu
! Flux de masse vertical
      DO ij=iip2,ip1jm
        fluxv(ij)=fluxz(ij)-fluxz(ij-1)-fluxm(ij)+fluxm(ij-iip1)
      ENDDO
      DO ij=iip1,ip1jmi1,iip1
         fluxv(ij+1)=fluxv(ij+iip1)
      END DO
! Poles
      fluxv(1)=-SUM(fluxm(1:iim))     
      fluxv(ip1jmp1)=SUM(fluxm(ip1jm-iim:ip1jm-1))
      fluxv=fluxv

! Meridional heat fluxes 
      DO ij=1,ip1jm
          ! centered scheme
          fluxtm(ij,1)=fluxm(ij)*(ts(ij+iip1,1)+ts(ij,1))/2.
          fluxtm(ij,2)=-fluxm(ij)*(ts(ij+iip1,2)+ts(ij,2))/2.
      END DO

! Zonal heat fluxes
! Schema upstream      
      fluxtz(1:iip1,:)=0 ! no zonal heat flux at north pole
      DO ij=iip2,ip1jm
      IF (fluxz(ij).GE.0.) THEN
             fluxtz(ij,1)=fluxz(ij)*ts(ij,1)
             fluxtz(ij,2)=-fluxz(ij)*ts(ij+1,2)
      ELSE
             fluxtz(ij,1)=fluxz(ij)*ts(ij+1,1)
             fluxtz(ij,2)=-fluxz(ij)*ts(ij,2)
      ENDIF
      END DO
      DO ij=iip1*2,ip1jmp1,iip1
             fluxtz(ij,:)=fluxtz(ij-iim,:)
      END DO
                   
! Calcul de dT
      DO ij=iip2,ip1jm
         dt(ij,:)=fluxtz(ij-1,:)-fluxtz(ij,:)   &
                  +fluxtm(ij,:)-fluxtm(ij-iip1,:)
         IF (fluxv(ij).GT.0.) THEN
           dt(ij,1)=dt(ij,1)+fluxv(ij)*ts(ij,2)
           dt(ij,2)=dt(ij,2)-fluxv(ij)*ts(ij,2)
         ELSE
           dt(ij,1)=dt(ij,1)+fluxv(ij)*ts(ij,1)
           dt(ij,2)=dt(ij,2)-fluxv(ij)*ts(ij,1)
         ENDIF
         dt(ij,:)=dt(ij,:)/aire(ij)
      END DO
      DO ij=iip1,ip1jmi1,iip1
         dt(ij+1,:)=dt(ij+iip1,:) 
      END DO
! Pôles
      dt(1,:)=SUM(fluxtm(1:iim,:),dim=1)
        IF (fluxv(1).GT.0.) THEN
          dt(1,1)=dt(1,1)+fluxv(1)*ts(1,2)
          dt(1,2)=dt(1,2)-fluxv(1)*ts(1,2)
        ELSE
          dt(1,1)=dt(1,1)+fluxv(1)*ts(1,1)
          dt(1,2)=dt(1,2)-fluxv(1)*ts(1,1)
        ENDIF
      dt(1,:)=dt(1,:)/apoln
      dt(ip1jmp1,:)=-SUM(fluxtm(ip1jm-iim:ip1jm-1,:),dim=1)
       IF (fluxv(ip1jmp1).GT.0.) THEN
         dt(ip1jmp1,1)=dt(ip1jmp1,1)+fluxv(ip1jmp1)*ts(ip1jmp1,2)
         dt(ip1jmp1,2)=dt(ip1jmp1,2)-fluxv(ip1jmp1)*ts(ip1jmp1,2)
       ELSE
         dt(ip1jmp1,1)=dt(ip1jmp1,1)+fluxv(ip1jmp1)*ts(ip1jmp1,1)
         dt(ip1jmp1,2)=dt(ip1jmp1,2)-fluxv(ip1jmp1)*ts(ip1jmp1,1)
       ENDIF
      dt(ip1jmp1,:)=dt(ip1jmp1,:)/apols
      dt(2:iip1,1)=dt(1,1)
      dt(2:iip1,2)=dt(1,2)
      dt(ip1jm+1:ip1jmp1,1)=dt(ip1jmp1,1)
      dt(ip1jm+1:ip1jmp1,2)=dt(ip1jmp1,2)

! Retour grille physique
      CALL gr_dyn_fi(2,iip1,jjp1,ngrid,dt,dt_phy)


  END SUBROUTINE slab_ekman2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE gr_fi_dyn(nfield,ngrid,im,jm,pfi,pdyn)
  ! Transfer a variable on global "physics" grid to global "dynamics" grid
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: im,jm,ngrid,nfield
  REAL,INTENT(IN) :: pfi(ngrid,nfield)
  REAL,INTENT(OUT) :: pdyn(im,jm,nfield)
  
  INTEGER :: i,j,ifield,ig
  
  DO ifield=1,nfield
    ! Handle poles
    DO i=1,im
      pdyn(i,1,ifield)=pfi(1,ifield)
      pdyn(i,jm,ifield)=pfi(ngrid,ifield)
    ENDDO
    ! Other points
    DO j=2,jm-1
      ig=2+(j-2)*(im-1)
      CALL SCOPY(im-1,pfi(ig,ifield),1,pdyn(1,j,ifield),1)
      pdyn(im,j,ifield)=pdyn(1,j,ifield)
    ENDDO
  ENDDO ! of DO ifield=1,nfield
  
  END SUBROUTINE gr_fi_dyn
  


  SUBROUTINE gr_dyn_fi(nfield,im,jm,ngrid,pdyn,pfi)
  ! Transfer a variable on global "dynamics" grid to global "physics" grid
  IMPLICIT NONE
  
  INTEGER,INTENT(IN) :: im,jm,ngrid,nfield
  REAL,INTENT(IN) :: pdyn(im,jm,nfield)
  REAL,INTENT(OUT) :: pfi(ngrid,nfield)
  
  INTEGER j,ifield,ig

  ! Sanity check:
  IF(ngrid.NE.2+(jm-2)*(im-1)) THEN
    WRITE(*,*) "gr_dyn_fi error, wrong sizes"
    STOP
  ENDIF
  
  ! Handle poles
  CALL SCOPY(nfield,pdyn,im*jm,pfi,ngrid)
  CALL SCOPY(nfield,pdyn(1,jm,1),im*jm,pfi(ngrid,1),ngrid)
  ! Other points
  DO ifield=1,nfield
    DO j=2,jm-1
      ig=2+(j-2)*(im-1)
      CALL SCOPY(im-1,pdyn(1,j,ifield),1,pfi(ig,ifield),1)
    ENDDO
  ENDDO

  END SUBROUTINE gr_dyn_fi



  SUBROUTINE  grad(klevel,pg,pgx,pgy)
  ! compute the covariant components x,y of the gradient of pg
  USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat
  IMPLICIT NONE
  
  INTEGER,INTENT(IN) :: klevel
  REAL,INTENT(IN) :: pg((nbp_lon+1)*nbp_lat,klevel)
  REAL,INTENT(OUT) :: pgx((nbp_lon+1)*nbp_lat,klevel)
  REAL,INTENT(OUT) :: pgy((nbp_lon+1)*(nbp_lat-1),klevel)

  INTEGER :: l,ij
  INTEGER :: iim,iip1,ip1jm,ip1jmp1
  
  iim=nbp_lon
  iip1=nbp_lon+1
  ip1jm=(nbp_lon+1)*(nbp_lat-1) ! = iip1*jjm
  ip1jmp1=(nbp_lon+1)*nbp_lat ! = iip1*jjp1
  
  DO l=1,klevel
    DO ij=1,ip1jmp1-1
      pgx(ij,l)=pg(ij+1,l)-pg(ij,l)
    ENDDO
    ! correction for pgx(ip1,j,l) ...
    ! ... pgx(iip1,j,l)=pgx(1,j,l) ...
    DO ij=iip1,ip1jmp1,iip1
      pgx(ij,l)=pgx(ij-iim,l)
    ENDDO
    DO ij=1,ip1jm
      pgy(ij,l)=pg(ij,l)-pg(ij+iip1,l)
    ENDDO
  ENDDO
  
  END SUBROUTINE grad



  SUBROUTINE diverg(klevel,x,y,div)
  ! compute the divergence of a vector field of components
  ! x,y. y and y being covriant components
  USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat
  IMPLICIT NONE
  
  INTEGER,INTENT(IN) :: klevel
  REAL,INTENT(IN) :: x((nbp_lon+1)*nbp_lat,klevel)
  REAL,INTENT(IN) :: y((nbp_lon+1)*(nbp_lat-1),klevel)
  REAL,INTENT(OUT) :: div((nbp_lon+1)*nbp_lat,klevel)
  
  INTEGER :: l,ij
  INTEGER :: iim,iip1,iip2,ip1jm,ip1jmp1,ip1jmi1
  
  REAL :: aiy1(nbp_lon+1),aiy2(nbp_lon+1)
  REAL :: sumypn,sumyps
  REAL,EXTERNAL :: SSUM
  
  iim=nbp_lon
  iip1=nbp_lon+1
  iip2=nbp_lon+2
  ip1jm=(nbp_lon+1)*(nbp_lat-1) ! = iip1*jjm
  ip1jmp1=(nbp_lon+1)*nbp_lat ! = iip1*jjp1
  ip1jmi1=(nbp_lon+1)*(nbp_lat-1)-(nbp_lon+1) ! = ip1jm - iip1
  
  DO l=1,klevel
    DO ij=iip2,ip1jm-1
      div(ij+1,l)= &
        cvusurcu(ij+1)*x(ij+1,l)-cvusurcu(ij)*x(ij,l)+ &
        cuvsurcv(ij-iim)*y(ij-iim,l)-cuvsurcv(ij+1)*y(ij+1,l) 
    ENDDO
    ! correction for div(1,j,l) ...
    ! ... div(1,j,l)= div(iip1,j,l) ...
    DO ij=iip2,ip1jm,iip1
      div(ij,l)=div(ij+iim,l)
    ENDDO
    ! at the poles
    DO ij=1,iim
      aiy1(ij)=cuvsurcv(ij)*y(ij,l)
      aiy2(ij)=cuvsurcv(ij+ip1jmi1)*y(ij+ip1jmi1,l)
    ENDDO
    sumypn=SSUM(iim,aiy1,1)/apoln
    sumyps=SSUM(iim,aiy2,1)/apols
    DO ij=1,iip1
      div(ij,l)=-sumypn
      div(ij+ip1jm,l)=sumyps
    ENDDO
  ENDDO ! of DO l=1,klevel
  
  !!! CALL filtreg( div, jjp1, klevel, 2, 2, .TRUE., 1 )
  
  DO l=1,klevel
    DO ij=iip2,ip1jm
      div(ij,l)=div(ij,l)*unsaire(ij) 
    ENDDO
  ENDDO
  
  END SUBROUTINE diverg



  SUBROUTINE gr_v_scal(nx,x_v,x_scal)
  USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat
  IMPLICIT NONE
  
  INTEGER,INTENT(IN) :: nx
  REAL,INTENT(IN) :: x_v((nbp_lon+1)*(nbp_lat-1),nx)
  REAL,INTENT(OUT) :: x_scal((nbp_lon+1)*nbp_lat,nx)
  
  INTEGER :: l,ij
  INTEGER :: iip1,iip2,ip1jm,ip1jmp1
  
  iip1=nbp_lon+1
  iip2=nbp_lon+2
  ip1jm=(nbp_lon+1)*(nbp_lat-1) ! = iip1*jjm
  ip1jmp1=(nbp_lon+1)*nbp_lat ! = iip1*jjp1
  
  DO l=1,nx
    DO ij=iip2,ip1jm
      x_scal(ij,l)= &
                   (airev(ij-iip1)*x_v(ij-iip1,l)+airev(ij)*x_v(ij,l)) &
                  /(airev(ij-iip1)+airev(ij))
    ENDDO
    DO ij=1,iip1
      x_scal(ij,l)=0.
    ENDDO
    DO ij=ip1jm+1,ip1jmp1
      x_scal(ij,l)=0.
    ENDDO
  ENDDO

  END SUBROUTINE gr_v_scal

  SUBROUTINE gr_scal_v(nx,x_scal,x_v)
  ! convert values from scalar points to v points on C-grid
  ! used to compute wind stress at V points
  USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: nx ! number of levels or fields
  REAL,INTENT(OUT) :: x_v((nbp_lon+1)*(nbp_lat-1),nx)
  REAL,INTENT(IN) :: x_scal((nbp_lon+1)*nbp_lat,nx)

  INTEGER :: l,ij
  INTEGER :: iip1,ip1jm

  iip1=nbp_lon+1
  ip1jm=(nbp_lon+1)*(nbp_lat-1) ! = iip1*jjm

      DO l=1,nx
        DO ij=1,ip1jm
          x_v(ij,l)= &
            (cu(ij)*cvusurcu(ij)*x_scal(ij,l)+ &
            cu(ij+iip1)*cvusurcu(ij+iip1)*x_scal(ij+iip1,l)) &
            /(cu(ij)*cvusurcu(ij)+cu(ij+iip1)*cvusurcu(ij+iip1))
        ENDDO
      ENDDO

  END SUBROUTINE gr_scal_v

  SUBROUTINE gr_scal_u(nx,x_scal,x_u)
  ! convert values from scalar points to U points on C-grid
  ! used to compute wind stress at U points
  USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: nx
  REAL,INTENT(OUT) :: x_u((nbp_lon+1)*nbp_lat,nx)
  REAL,INTENT(IN) :: x_scal((nbp_lon+1)*nbp_lat,nx)

  INTEGER :: l,ij
  INTEGER :: iip1,jjp1,ip1jmp1

  iip1=nbp_lon+1
  jjp1=nbp_lat
  ip1jmp1=(nbp_lon+1)*nbp_lat ! = iip1*jjp1

  DO l=1,nx
     DO ij=1,ip1jmp1-1
        x_u(ij,l)= &
         (aire(ij)*x_scal(ij,l)+aire(ij+1)*x_scal(ij+1,l)) &
         /(aire(ij)+aire(ij+1))
     ENDDO
  ENDDO

  CALL SCOPY(nx*jjp1,x_u(1,1),iip1,x_u(iip1,1),iip1)

  END SUBROUTINE gr_scal_u
  
END MODULE surf_heat_transp_mod



