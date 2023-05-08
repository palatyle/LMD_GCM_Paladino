
!     SUBROUTINE     WSA_ROSA_NEW
!     SUBROUTINE     ITERWV          WSA pour un WV
!     SUBROUTINE     BRACWV          Bracket de ITERWV
!     SUBROUTINE     BRACWSA         Bracket de KEEQ
!     FUNCTION       IRFRMWV         Iterative Root Finder Ridder's Method for WV
!     FUNCTION       IRFRMSA         Iterative Root Finder Ridder's Method for SA
!     FUNCTION       KEEQ            Kelvin Equation EQuality
!     FUNCTION       WVCOND          H2O Condensation with WSA, T, P and H2SO4tot


!----------------------------------------------------------------------------
SUBROUTINE  WSA_ROSA_NEW(TAIR,PAIR,RADIUS,WSAS,MSAD)

  !*     This subroutine calculates the acid mass fraction, density, and
  !*     mass of sulfuric acid in a single aerosol droplet of a specified 
  !*     radius in equilibrium with ambient water vapor partial pressure 
  !*     and temperature.
  !*
  !*     The calculation is performed by iteration of
  !*        ln(PPWV) - [(2Mh2o sigma)/(R T r rho) - ln(ph2osa)] = 0
  !*     using the secant method. Vapor pressures by Gmitro and Vermeulen
  !*     (PWVSAS_GV) are used.  
  !* Zeleznik valid only up to 350 K
  !*
  !*     Input/output variables:
  !*     REAL(KIND=4)  RADIUS,TAIR,PPWV,WSAS,RHOSA,MSA
  !*
  !*     Input:       
  !*         RADIUS:  m         Radius of aerosol droplet
  !*         TAIR:    K         Temperature of ambient air 
  !*         PPWV:    Pa        Partial pressure of ambient water vapor 
  !*
  !*     Output:
  !*         WSAS:              mass fraction of sulfuric acid. [0.1;1]
  !*         RHOSA:  kg/m**3   Density of sulfuric acid solution droplet
  !*         MSAD:     kg        Mass of sulfuric acid in droplet
  !*     modified from
  !*     PROGRAM PSC_MODEL_E
  !*     by A. Määttänen & Slimane Bekki
  !*     subroutine for LMDZ+photochemistry VENUS
  !*     by A. Stolzenbach
  !*
  !*     Modified by S.Guilbon for microphysical module to Venus GCM
  !*

  USE donnees
  USE free_param

  IMPLICIT NONE

  ! Inputs
  REAL, intent(in) :: RADIUS, TAIR, PAIR
  ! Outputs
  REAL, intent(out) :: WSAS, MSAD
  ! Auxilary variables:
  REAL :: mrt_wv, mrt_sa
  REAL :: N_H2SO4, N_H2O
  REAL :: H2SO4_liq, H2O_liq
  REAL :: CONCM
  REAL :: MCONDTOT
  REAL :: RMODE
  REAL :: WSAFLAG
  REAL :: power
  !     Ridder's Method variables:
  REAL :: WVMIN, WVMAX, WVACC

  INTEGER :: NBROOT

  INTEGER :: MAXITE
  PARAMETER(MAXITE=20)

  INTEGER :: NBRAC
  PARAMETER(NBRAC=5)

  INTEGER :: FLAG

  ! External functions needed:
  REAL :: IRFRMWV

  ! Physical constants:
  REAL :: MH2O

  ! External functions needed:
  REAL ::  PWVSAS_GV,STSAS,ROSAS
  ! PWVSAS_GV:   Natural logaritm of water vapor pressure over
  !                      sulfuric acid solution
  ! STSAS:       Surface tension of sulfuric acid solution
  ! ROSAS:       Density of sulfuric acid solution

  !     Auxiliary local variables:
  REAL :: DELW,DELLP,C1,C2,W0,W1,W2,F0,F1,WGUESS,LPPWV,RO
  REAL :: psatwv,watact
  INTEGER :: ITERAT 
!  write(*,*)'WSA ROSA NEW', RADIUS
  MH2O=MWV

  C1=2.0D0*MH2O/RGAS
  C2=4.0D0*PI/3.0D0

  mrt_sa=ppsa/pair
  mrt_wv=ppwv/pair

  ! Initialisation des bornes pour WV
  WVMIN=1.D-35
  WVMAX=mrt_wv

  ! Accuracy de WVeq
  WVACC=WVMAX*1.0D-3

  !  BRACWV borne la fonction f(WV) - WV = 0
  !  de WV=0  WV=WVtot on cherche l'intervalle o f(WV) - WV = 0
  !  avec prcisment f(WVliq de WSA<=WVinput) + WVinput - WVtot = 0 
  !  Elle fait appel  la fct/ssrtine ITERWV()
  CALL BRACWV(TAIR,PAIR,WVMIN,WVMAX,NBRAC,RADIUS,mrt_wv,mrt_sa,FLAG,WSAFLAG,NBROOT)

  SELECT CASE(FLAG)

  CASE(1) 
     ! Cas NROOT=1 ou NROOT>1 mais dans un intervalle restreint WVTOT (cas courant)         
     ! IRFRMWV Ridder's method pour trouver, sur [WVmin,WVmax], WVo tel que f(WVo) - WVo = 0 
     ! Elle fait appel  la fct/ssrtine ITERWV()

     WSAS=IRFRMWV(TAIR,PAIR,WVMIN,WVMAX,WVACC,MAXITE,RADIUS,mrt_wv,mrt_sa,NBROOT)
     RHOSA = ROSAS(TAIR,WSAS)
     MSAD = C2*WSAS*RHOSA*RADIUS**3

  CASE(2)
     ! Cas NROOT=0 mais proche de 0 
     WSAS=WSAFLAG      
     RHOSA=ROSAS(TAIR,WSAS)
     MSAD=C2*WSAS*RHOSA*RADIUS**3
     !     ATTENTION ce IF ne sert a rien en fait, juste a retenir une situation
     !     ubuesque dans mon code ou sans ce IF les valeurs de rho_droplets sont
     !     incohrentes avec TT et WH2SO4 (a priori lorsque NTOT=0) 
     !     Juste le fait de METTRE un IF fait que rho_droplet a la bonne valeur
     !     donne par ROSAS (cf test externe en Python), sinon, la valeur est trop
     !     basse (de l'ordre de 1000 kg/m3) et correspond parfois  la valeur avec
     !     WSA=0.1 (pas totalement sr) 
     !     En tous cas, incohrent avec ce qui est attendue pour le WSA et T donn
     !     La version avec le IF (rho<1100 & WSA>0.1) est CORRECTE, rho_droplet a 
     !     la bonne valeur (tests externes Python confirment)

     IF ((RHOSA.LT.1100.0D0).AND. (WSAS.GT.0.1D0))THEN
        PRINT*,'PROBLEM RHO_DROPLET'
        PRINT*,'rho_droplet',RHOSA
        PRINT*,'T',TAIR,'WSA',WSAS
        PRINT*,'ROSAS',ROSAS(TAIR, WSAS)
        PRINT*,'FLAG',FLAG,'NROOT',NBROOT
        STOP
     ENDIF

  CASE(3)
     write(*,*)'Case 0 NROOT'  
     RHOSA=0.0D+0
     WSAS=0.0D+0
     MSAD=0.0D+0

  END SELECT



  RETURN 

END SUBROUTINE WSA_ROSA_NEW


!*****************************************************************************                              
SUBROUTINE ITERWV(TAIR,PAIR,WV,WVLIQ,WVEQOUT,WVTOT,WSAOUT,SATOT,RADIUS)
  !* Cette routine est la solution par itration afin de trouver WSA pour un WV,
  !* et donc LPPWV, donn. Ce qui nous donne egalement le WV correspondant au 
  !* WSA solution   
  !* For VenusGCM by A. Stolzenbach 07/2014 
  !* OUTPUT: WVEQ et WSAOUT

  USE donnees
  USE free_param
  IMPLICIT NONE

  REAL, INTENT(IN) :: TAIR, PAIR
  REAL, INTENT(IN) :: WV, WVTOT, SATOT, RADIUS 
  REAL, INTENT(OUT) :: WVEQOUT, WSAOUT, WVLIQ

  REAL :: LPPWV

  REAL :: WSAMIN, WSAMAX, WSAACC
  PARAMETER(WSAACC=0.01D0)

  INTEGER :: MAXITSA, NBRACSA, NBROOT
  PARAMETER(MAXITSA=20)
  PARAMETER(NBRACSA=5)

  LOGICAl :: FLAG1,FLAG2

  ! External Function      
  REAL :: IRFRMSA, WVCOND

  IF (RADIUS.LT.1D-30) THEN
     PRINT*,'RMODE == 0 FLAG 3', RADIUS
     STOP
  ENDIF

  ! Initialisation WSA=[0.1,1.0]      
  WSAMIN = 0.1D0
  WSAMAX = 1.0D0
  LPPWV = DLOG(PAIR*WV)

  ! Appel Bracket de KEEQ         
  CALL BRACWSA(WSAMIN,WSAMAX,NBRACSA,RADIUS,TAIR,LPPWV,FLAG1,FLAG2,NBROOT)

  IF ((.NOT.FLAG1).AND.(.NOT.FLAG2).AND.(NBROOT.EQ.1)) THEN          

     WSAOUT=IRFRMSA(TAIR,PAIR,WSAMIN,WSAMAX,WSAACC,MAXITSA,RADIUS,LPPWV,NBROOT)

!!$     ! AM uncommented the two following lines to avoid problems with nucleation
!!$     IF (WSAOUT.GT.1.0) WSAOUT=0.999999
!!$     IF (WSAOUT.LT.0.1) WSAOUT=0.1
!!$     write(*,*) 'in 1 wsaout 2', WSAOUT

     ! Si BRACWSA ne trouve aucun ensemble solution KEEQ=0 on fixe WSA a 0.9999 ou 0.1
  ELSE
     IF (FLAG1.AND.(.NOT.FLAG2)) WSAOUT = 0.999999D0
     IF (FLAG2.AND.(.NOT.FLAG1)) WSAOUT = WSAMIN
     IF (FLAG1.AND.FLAG2) THEN
        PRINT*,'FLAGs BARCWSA tous TRUE'
        STOP
     ENDIF
  ENDIF

  !     WVEQ output correspondant a WVliq lie a WSA calcule
  WVLIQ=WVCOND(WSAOUT,TAIR,PAIR,SATOT)
  WVEQOUT=(WVLIQ+WV)/WVTOT-1.0D0

END SUBROUTINE ITERWV


!*****************************************************************************                             
SUBROUTINE BRACWV(TAIR,PAIR,XA,XB,N,RADIUS,WVTOT,SATOT,FLAGWV,WSAFLAG,NROOT)

  !* Bracket de ITERWV
  !* From Numerical Recipes     
  !* Adapted for VenusGCM A. Stolzenbach 07/2014 
  !* X est WVinput
  !* OUTPUT: XA et XB      

  USE donnees
  USE free_param

  IMPLICIT NONE

  REAL, INTENT(IN) :: WVTOT,SATOT,RADIUS,TAIR, PAIR
  INTEGER, INTENT(IN) :: N

  REAL, INTENT(INOUT) :: XA,XB
  REAL, INTENT(OUT) :: WSAFLAG

  INTEGER :: I,J

  INTEGER, INTENT(OUT) :: NROOT

  REAL :: FP, FC, X, WVEQ, WVLIQ, WSAOUT
  REAL :: XMAX,XMIN,WVEQACC

  INTEGER, INTENT(OUT) :: FLAGWV
!  write(*,*)'BRACWV', RADIUS
  ! WVEQACC est le seuil auquel on accorde un WSA correct meme
  ! si il ne fait pas partie d'une borne. Utile quand le modele
  ! s'approche de 0 mais ne l'atteint pas.
  WVEQACC = 1.0D-3  
  FLAGWV = 1
  NROOT = 0

!     25/11/2016
!     On change ordre on va du max au min
  X = XB
  XMAX = XB
  XMIN = XA

  ! CAS 1 On borne la fonction (WVEQ=0)
  CALL ITERWV(TAIR,PAIR,X,WVLIQ,WVEQ,WVTOT,WSAOUT,SATOT,RADIUS)

  FP=WVEQ

  DO I=N-1,1,-1
     X=(1.-DLOG(DBLE(N-I))/DLOG(DBLE(N)))*XMAX

     CALL ITERWV(TAIR,PAIR,X,WVLIQ,WVEQ,WVTOT,WSAOUT,SATOT,RADIUS)

     FC=WVEQ

     IF ((FP*FC).LT.0.D0) THEN  
        NROOT=NROOT+1
        ! Si NROOT>1 on place la borne sup output  la borne min du calcul en i             
        IF (NROOT.GT.1) THEN
           XB=(1.-DLOG(DBLE(I+1))/DLOG(DBLE(N)))*XMAX
        ENDIF

        IF (I.EQ.1) THEN
           XA=XMIN
        ELSE
           XA=X
        ENDIF
        RETURN
     ENDIF
     FP=FC
  ENDDO

  ! CAS 2 on refait la boucle pour tester si WVEQ est proche de 0 
  ! avec le seuil WVEQACC
  IF (NROOT.EQ.0) THEN
     X=XMAX

     CALL ITERWV(TAIR,PAIR,X,WVLIQ,WVEQ,WVTOT,WSAOUT,SATOT,RADIUS)

     DO J=N-1,1,-1
        X=(1.-DLOG(DBLE(N-J))/DLOG(DBLE(N)))*XMAX
        !             write(*,*) 'BRACWV, bf 4th ITERWV (cas 2) '
        CALL ITERWV(TAIR,PAIR,X,WVLIQ,WVEQ,WVTOT,WSAOUT,SATOT,RADIUS)

        IF (ABS(WVEQ).LE.WVEQACC) THEN
           WSAFLAG=WSAOUT
           FLAGWV=2
           RETURN
        ENDIF
     ENDDO

     !     CAS 3 Pas de borne, WVEQ jamais proche de 0          
     FLAGWV=3
     RETURN
  ENDIF

END SUBROUTINE BRACWV


!*****************************************************************************
SUBROUTINE BRACWSA(XA,XB,N,RADIUS,TAIR,LPPWVINP,FLAGH,FLAGL,NROOT)

  !* Bracket de KEEQ
  !* From Numerical Recipes     
  !* Adapted for Venus GCM A. Stolzenbach 07/2014 

  USE donnees
  USE free_param
  IMPLICIT NONE

  !     External functions needed:
  REAL :: KEEQ

  REAL, INTENT(IN) :: RADIUS,TAIR,LPPWVINP 
  INTEGER, INTENT(IN) :: N

  REAL, INTENT(INOUT) :: XA,XB

  INTEGER, INTENT(OUT) ::  NROOT

  INTEGER :: I, J

  REAL :: DX, FP, FC, X

  LOGICAL, INTENT(OUT) :: FLAGH,FLAGL

  FLAGL=.FALSE.
  FLAGH=.FALSE.    
  NROOT=0
  DX=(XB-XA)/N
  X=XB
  FP=KEEQ(TAIR,RADIUS,X,LPPWVINP)

  DO I=N,1,-1
     X=X-DX

     FC=KEEQ(TAIR,RADIUS,X,LPPWVINP)

     IF ((FP*FC).LE.0.) THEN
        NROOT=NROOT+1
        XA=X
        XB=X+DX
        RETURN
     ENDIF

     FP=FC
  ENDDO

  IF (NROOT.EQ.0) THEN
     ! Test determine la tendance globale KEEQ sur [WSAMIN,WSAMAX]        
     IF ((ABS(KEEQ(TAIR,RADIUS,XA,LPPWVINP))- &
          &    ABS(KEEQ(TAIR,RADIUS,XB,LPPWVINP))).GT.0.0) FLAGH=.TRUE.
     ! On fixe flag low TRUE pour WSA = 0.1
     IF ((ABS(KEEQ(TAIR,RADIUS,XA,LPPWVINP))- &
          &    ABS(KEEQ(TAIR,RADIUS,XB,LPPWVINP))).LT.0.0) FLAGL=.TRUE.
  ENDIF

END SUBROUTINE BRACWSA


!*****************************************************************************
FUNCTION IRFRMWV(TAIR,PAIR,X1,X2,XACC,MAXIT,RADIUS,WVTOT,SATOT,NROOT)

  !* Iterative Root Finder Ridder's Method for Water Vapor calculus
  !* From Numerical Recipes
  !* Adapted for VenusGCM A. Stolzenbach 07/2014

  !* Les iterations sur [X1,X2] sont [WV1,WV2]
  !* la variable X est WV
  !* IRFRMWV sort en OUTPUT : WSALOC pour ITERWV=0 (ou WVEQ=0)

  USE donnees
  USE free_param
  IMPLICIT NONE

  REAL, INTENT(IN) :: TAIR, PAIR
  REAL, INTENT(IN) :: X1, X2
  REAL, INTENT(IN) :: XACC
  REAL :: IRFRMWV 
  INTEGER, INTENT(IN) :: MAXIT,NROOT

  ! LOCAL VARIABLES
  REAL :: XL, XH, XM, XNEW, X
  REAL :: WSALOC, WVEQ, WVLIQ
  REAL :: FL, FH, FM, FNEW
  REAL :: ANS, S, FSIGN
  INTEGER i
  LOGICAL :: FLAGH,FLAGL

  ! External variables needed:
  REAL, INTENT(IN) :: WVTOT,SATOT
  REAL, INTENT(IN) :: RADIUS

  ! Initialisation
  X=X1
  CALL ITERWV(TAIR,PAIR,X,WVLIQ,WVEQ,WVTOT,WSALOC,SATOT,RADIUS) 
  FL=WVEQ
  X=X2
  CALL ITERWV(TAIR,PAIR,X,WVLIQ,WVEQ,WVTOT,WSALOC,SATOT,RADIUS)
  FH=WVEQ

  ! Test Bracketed values 
  IF (((FL.LT.0.).AND.(FH.GT.0.)).OR. &
       &   ((FL.GT.0.).AND.(FH.LT.0.))) &
       &  THEN
     XL=X1
     XH=X2
     ANS=-1.D38

     DO i=1, MAXIT
        XM=0.5D0*(XL+XH)
        CALL ITERWV(TAIR,PAIR,XM,WVLIQ,WVEQ,WVTOT,WSALOC,SATOT,RADIUS)
        FM=WVEQ
        S=SQRT(FM*FM-FL*FH)

        IF (S.EQ.0.0) THEN
           IRFRMWV=WSALOC
           RETURN
        ENDIF

        IF (FL.GT.FH) THEN
           FSIGN=1.0D0
        ELSE
           FSIGN=-1.0D0
        ENDIF

        XNEW=XM+(XM-XL)*(FSIGN*FM/S) 

        IF (ABS(XNEW-ANS).LE.XACC) THEN
           IRFRMWV=WSALOC
           RETURN
        ENDIF

        ANS=XNEW
        CALL ITERWV(TAIR,PAIR,ANS,WVLIQ,WVEQ,WVTOT,WSALOC,SATOT,RADIUS)
        FNEW=WVEQ

        IF (FNEW.EQ.0.D0) THEN
           IRFRMWV=WSALOC
           RETURN
        ENDIF

        IF (SIGN(FM, FNEW).NE.FM) THEN
           XL=XM
           FL=FM
           XH=ANS
           FH=FNEW
        ELSEIF (SIGN(FL, FNEW).NE.FL) THEN
           XH=ANS
           FH=FNEW
        ELSEIF (SIGN(FH, FNEW).NE.FH) THEN
           XL=ANS
           FL=FNEW
        ELSE
           PRINT*,'PROBLEM IRFRMWV dans new_cloud_venus'
           PRINT*,'you shall not PAAAAAASS'
           STOP
        ENDIF
     ENDDO
     PRINT*,'Paaaaas bien MAXIT atteint'
     PRINT*,'PROBLEM IRFRMWV dans new_cloud_venus'
     PRINT*,'you shall not PAAAAAASS'
     XL=X1
     XH=X2
     !         ANS=-9.99e99
     ANS=-1.D38

     DO i=1, MAXIT
        XM=0.5D0*(XL+XH)
        CALL ITERWV(TAIR,PAIR,XM,WVLIQ,WVEQ,WVTOT,WSALOC,SATOT,RADIUS)
        FM=WVEQ
        S=SQRT(FM*FM-FL*FH)
        IF (FL.GT.FH) THEN
           FSIGN=1.0D0
        ELSE
           FSIGN=-1.0D0
        ENDIF

        XNEW=XM+(XM-XL)*(FSIGN*FM/S)     

        ANS=XNEW
        CALL ITERWV(TAIR,PAIR,ANS,WVLIQ,WVEQ,WVTOT,WSALOC,SATOT,RADIUS)
        FNEW=WVEQ
        PRINT*,'WVliq',WVLIQ,'WVtot',WVTOT,'WVeq',WVEQ
        PRINT*,'WSA',WSALOC,'SAtot',SATOT
        PRINT*,'T',TAIR,'P',PAIR

        IF (SIGN(FM, FNEW).NE.FM) THEN
           XL=XM
           FL=FM
           XH=ANS
           FH=FNEW
        ELSEIF (SIGN(FL, FNEW).NE.FL) THEN
           XH=ANS
           FH=FNEW
        ELSEIF (SIGN(FH, FNEW).NE.FH) THEN
           XL=ANS
           FL=FNEW
        ELSE
           PRINT*,'PROBLEM IRFRMWV dans new_cloud_venus'
           PRINT*,'you shall not PAAAAAASS TWIIICE???'
           STOP
        ENDIF
     ENDDO
     STOP
  ELSE
     PRINT*,'IRFRMWV must be bracketed'
     PRINT*,'NROOT de BRACWV', NROOT
     IF (ABS(FL).LT.XACC) THEN
        PRINT*,'IRFRMWV FL == 0',FL 
        PRINT*,'X1',X1,'X2',X2,'FH',FH
        CALL ITERWV(TAIR,PAIR,X1,WVLIQ,WVEQ,WVTOT,WSALOC,SATOT,RADIUS)
        IRFRMWV=WSALOC
        RETURN
     ENDIF
     IF (ABS(FH).LT.XACC) THEN
        PRINT*,'IRFRMWV FH == 0',FH
        PRINT*,'X1',X1,'X2',X2,'FL',FL
        CALL ITERWV(TAIR,PAIR,X2,WVLIQ,WVEQ,WVTOT,WSALOC,SATOT,RADIUS)
        IRFRMWV=WSALOC
        RETURN
     ENDIF
     IF ((ABS(FL).GT.XACC).AND.(ABS(FH).GT.XACC)) THEN 
        PRINT*,'STOP dans IRFRMWV avec rien == 0'
        PRINT*,'X1',X1,'X2',X2
        PRINT*,'Fcalc',FL,FH
        PRINT*,'T',TAIR,'P',PAIR,'R',RADIUS
        STOP   
     ENDIF
     IF ((ABS(FL).LT.XACC).AND.(ABS(FH).LT.XACC)) THEN 
        PRINT*,'STOP dans IRFRMWV Trop de solution < WVACC'
        PRINT*,FL,FH
        STOP   
     ENDIF

  end IF
END FUNCTION IRFRMWV

!*****************************************************************************                             
FUNCTION IRFRMSA(TAIR,PAIR,X1,X2,XACC,MAXIT,RADIUS,LPPWV,NB)

  !* Iterative Root Finder Ridder's Method for Sulfuric Acid calculus
  !* From Numerical Recipes
  !* Adapted for VenusGCM A. Stolzenbach 07/2014
  !*
  !* Les iterations sur [X1,X2] sont [WSA1,WSA2]
  !* la variable X est WSA
  !* IRFRMSA sort en OUTPUT : WSA pour KEEQ=0

  use donnees
  use free_param
  IMPLICIT NONE

  REAL, INTENT(IN) :: TAIR, PAIR
  REAL, INTENT(IN) :: X1, X2
  REAL, INTENT(IN) :: XACC 
  INTEGER, INTENT(IN) :: MAXIT, NB

  !     LOCAL VARIABLES
  REAL :: IRFRMSA
  REAL :: XL, XH, XM, XNEW
  REAL :: Fl, FH, FM, FNEW
  REAL :: ANS, S, FSIGN
  INTEGER i

  !     External variables needed:
  REAL, INTENT(IN) :: LPPWV
  REAL, INTENT(IN) :: RADIUS

  !     External functions needed:
  REAL :: KEEQ

  !     Initialisation
  FL=KEEQ(TAIR,RADIUS,X1,LPPWV)
  FH=KEEQ(TAIR,RADIUS,X2,LPPWV)

  !     Test Bracketed values 
  IF (((FL.LT.0.D0).AND.(FH.GT.0.D0)).OR.((FL.GT.0.D0).AND.(FH.LT.0.D0)))  THEN

     XL=X1
     XH=X2

     ANS=-1.D38

     DO i=1, MAXIT
        XM=0.5D0*(XL+XH)
        FM=KEEQ(TAIR,RADIUS,XM,LPPWV)

        S=SQRT(FM*FM-FL*FH)

        IF (S.EQ.0.D0) THEN
           IRFRMSA=ANS
           RETURN
        ENDIF

        IF (FL.GT.FH) THEN
           FSIGN=1.0D0
        ELSE
           FSIGN=-1.0D0
        ENDIF

        XNEW=XM+(XM-XL)*(FSIGN*FM/S) 

        IF (ABS(XNEW-ANS).LE.XACC) THEN
           IRFRMSA=ANS

           RETURN
        ENDIF

        ANS=XNEW
        FNEW=KEEQ(TAIR,RADIUS,ANS,LPPWV)

        IF (FNEW.EQ.0.D0) THEN
           IRFRMSA=ANS
           RETURN
        ENDIF

        IF (SIGN(FM, FNEW).NE.FM) THEN
           XL=XM
           FL=FM
           XH=ANS
           FH=FNEW
        ELSEIF (SIGN(FL, FNEW).NE.FL) THEN
           XH=ANS
           FH=FNEW
        ELSEIF (SIGN(FH, FNEW).NE.FH) THEN
           XL=ANS
           FL=FNEW
        ELSE
           PRINT*,'PROBLEM IRFRMSA dans new_cloud_venus'
           PRINT*,'you shall not PAAAAAASS'
           STOP
        ENDIF

     ENDDO
     PRINT*,'Paaaaas bien MAXIT atteint'
     PRINT*,'PROBLEM IRFRMSA dans new_cloud_venus'
     PRINT*,'you shall not PAAAAAASS'
     XL=X1
     XH=X2
     PRINT*,'Borne XL',XL,'XH',XH

     ANS=-1.D38

     DO i=1, MAXIT
        XM=0.5D0*(XL+XH)
        FM=KEEQ(TAIR,RADIUS,XM,LPPWV)
        S=SQRT(FM*FM-FL*FH)

        IF (FL.GT.FH) THEN
           FSIGN=1.0D0
        ELSE
           FSIGN=-1.0D0
        ENDIF

        XNEW=XM+(XM-XL)*(FSIGN*FM/S)  
        ANS=XNEW
        FNEW=KEEQ(TAIR,RADIUS,ANS,LPPWV)
        PRINT*,'KEEQ result',FNEW,'T',TAIR,'R',RADIUS
        IF (SIGN(FM, FNEW).NE.FM) THEN
           XL=XM
           FL=FM
           XH=ANS
           FH=FNEW
        ELSEIF (SIGN(FL, FNEW).NE.FL) THEN
           XH=ANS
           FH=FNEW
        ELSEIF (SIGN(FH, FNEW).NE.FH) THEN
           XL=ANS
           FL=FNEW
        ELSE
           PRINT*,'PROBLEM IRFRMSA dans new_cloud_venus'
           PRINT*,'you shall not PAAAAAASS'
           STOP
        ENDIF
     ENDDO
     STOP

  ELSE
     PRINT*,'IRFRMSA must be bracketed'
     IF (FL.EQ.0.D0) THEN
        PRINT*,'IRFRMSA FL == 0',Fl
        IRFRMSA=X1
        RETURN
     ENDIF
     IF (FH.EQ.0.D0) THEN
        PRINT*,'IRFRMSA FH == 0',FH
        IRFRMSA=X2
        RETURN
     ENDIF
     IF ((FL.NE.0.).AND.(FH.NE.0.)) THEN
        PRINT*,'IRFRMSA FH and FL neq 0: ', FL, FH
        PRINT*,'X1',X1,'X2',X2 
        PRINT*,'Kind F', KIND(FL), KIND(FH)
        PRINT*,'Kind X', KIND(X1), KIND(X2)
        PRINT*,'Logical: ',(SIGN(FL,FH).NE.FL)
        PRINT*,'Logical: ',(SIGN(FH,FL).NE.FH)
        PRINT*,'nb root BRACWSA',NB
        STOP
     ENDIF

  ENDIF

END function IRFRMSA

!*****************************************************************************
FUNCTION KEEQ(TAIR,RADIUS,WX,LPPWV)

  !* Kelvin Equation EQuality
  !* ln(PPWV_eq) - (2Mh2o sigma)/(R T r rho) - ln(ph2osa) = 0

  use donnees
  use free_param
  IMPLICIT NONE

  REAL, INTENT(IN) :: RADIUS,WX,LPPWV,TAIR

  ! Physical constants:
  REAL :: KEEQ

  !     External functions needed:
  REAL :: PWVSAS_GV, STSAS, ROSAS
  !     PWVSAS_GV:      Natural logaritm of water vapor pressure over
  !                  sulfuric acid solution
  !     STSAS:       Surface tension of sulfuric acid solution
  !     ROSAS:       Density of sulfuric acid solution
  !
  !     Auxiliary local variables:
  REAL :: C1

  C1=2.0D0*MWV/RGAS

  KEEQ=LPPWV-C1*STSAS(TAIR,WX)/(TAIR*RADIUS*ROSAS(TAIR,WX))- &
       &     PWVSAS_GV(TAIR,WX)

END FUNCTION KEEQ


!*****************************************************************************
FUNCTION WVCOND(WX,T,P,SAt)

  !* Condensation de H2O selon WSA, T et P et H2SO4tot

  !* Adapted for VenusGCM A. Stolzenbach 07/2014
  !     INPUT:
  !     SAt  : VMR of total H2SO4
  !     WX: aerosol H2SO4 weight fraction (fraction)
  !     T: temperature (K)
  !     P: pressure (Pa)
  !     OUTPUT: 
  !	WVCOND : VMR H2O condense

  !      USE chemparam_mod

  use donnees
  use free_param

  IMPLICIT NONE

  REAL, INTENT(IN) :: SAt, WX
  REAL, INTENT(IN) :: T, P

  !     working variables
  REAL :: WVCOND
  REAL :: SA, WV, KBOLTZ
  REAL :: DND2,pstand,lpar,acidps
  REAL :: x1, satpacid
  REAL, DIMENSION(2):: act
  REAL :: CONCM
  REAL :: NH2SO4
  REAL :: H2OCOND, H2SO4COND

  KBOLTZ=KBZ
  CONCM= (P)/(KBOLTZ*T) !air number density, molec/m3? CHECK UNITS!

  NH2SO4=SAt*CONCM
  pstand=1.01325D+5 !Pa  1 atm pressure

  x1=(WX/MSA)/(WX/MSA + ((1.-WX)/MWV))

  CALL zeleznik(x1,T,act)

  !pure acid satur vapor pressure
  lpar= -11.695D0 + DLOG(pstand) ! Zeleznik
  acidps = 1/360.15D0 - 1.0/T + 0.38D0/545.D0*(1.0+DLOG(360.15D0/T)-360.15D0/T)
  acidps = 10156.D0*acidps + lpar
  acidps = DEXP(acidps)    !Pa

  !acid sat.vap.PP over mixture (flat surface):
  satpacid=act(2)*acidps ! Pa 

  !       Conversion from Pa to N.D #/m3
  DND2=satpacid/(KBOLTZ*T)
  !	H2SO4COND N.D #/m3 condensee ssi H2SO4>H2SO4sat
  IF (NH2SO4.GT.DND2) THEN
     H2SO4COND=NH2SO4-DND2
     !	calcul de H2O cond correspondant a H2SO4 cond
     H2OCOND=H2SO4COND*MSA*(1.0D0-WX)/(MWV*WX)
     !	Si on a H2SO4<H2SO4sat on ne condense rien, VMR = 1.0E-30 
  ELSE
     H2OCOND=1.0D-30*CONCM
  END IF

  !*****************************************************
  !     ATTENTION: Ici on ne prends pas en compte 
  !                si H2O en defaut!
  !                On veut la situation thorique
  !                 l'equilibre
  !*****************************************************	     	
  !	Test si H2O en defaut H2Ocond>H2O dispo
  !	IF ((H2OCOND.GT.NH2O).AND.(NH2SO4.GE.DND2)) THEN
  !	On peut alors condenser tout le H2O dispo
  !	H2OCOND=NH2O
  !	On met alors egalement a jour le H2SO4 cond correspondant au H2O cond
  !	H2SO4COND=H2OCOND*18.0153*WSA/(98.078*(1.0-WSA))
  !      END IF

  !     Calcul de H2O condense VMR          
  WVCOND=H2OCOND/CONCM

END FUNCTION WVCOND


!*****************************************************************************
FUNCTION PWVSAS_GV(T,W)

  !*     Natural logaritm of saturated water vapor pressure over plane
  !*     sulfuric acid solution.
  !*
  !*     Source: J.I.Gmitro & T.Vermeulen: A.I.Ch.E.J.  10,740,1964.
  !*             W.F.Giauque et al.: J. Amer. Chem. Soc. 82,62,1960.
  !*
  !*     The formula of Gmitro & Vermeulen for saturation pressure
  !*     is used:
  !*                 ln(p) = A ln(298/T) + B/T + C + DT
  !*     with values of A,B,C and D given by Gmitro & Vermeulen,
  !*     and calculated from partial molal properties given by Giauque et al.
  !* 
  !*     Input:  T: Temperature (K)
  !*             W: Weight fraction of H2SO4  [0;1] 
  !*     Output: Natural logaritm of water vapor pressure 
  !*             over sulfuric acid solution   ( ln(Pa) )
  !*
  !*     External functions needed for calculation of partial molal 
  !*     properties of pure components at 25 C as function of W.

  use donnees
      IMPLICIT NONE

      REAL :: CPH2O,ALH2O,FFH2O,LH2O
!     CPH2O:  Partial molal heat capacity of sulfuric acid solution.
!     ALH2O:  Temparature derivative of CPH2O
!     FFH2O:  Partial molal free energy of sulfuric acid solution.
!     LH2O:   Partial molal enthalpy of sulfuric acid
!
!
!
      REAL, INTENT(IN) :: T,W
      REAL :: PWVSAS_GV
      REAL :: ADOT,BDOT,CDOT,DDOT
      REAL :: RGAScal,MMHGPA
      REAL :: K1,K2
      REAL :: A,B,C,Dd,CP,L,F,ALFA
!     Physical constants given by Gmitro & Vermeulen:
      PARAMETER(                   &
              ADOT=-3.67340,       &
              BDOT=-4143.5,        &
              CDOT=10.24353,       &
              DDOT=0.618943d-3)
      PARAMETER(                   &
!     Gas constant (cal/(deg mole)):
           RGAScal=1.98726,        &
!     Natural logarith of conversion factor between atm. and Pa:     
           MMHGPA=11.52608845,     &
           K1=298.15,              &
           K2=K1*K1/2.0)
!
      INTEGER :: KLO,KHI,K,I,NPOINT
      PARAMETER(NPOINT=110)
      REAL, DIMENSION(NPOINT) :: WTAB(NPOINT)
      DATA (WTAB(I),I=1,NPOINT)/   &
     0.00000,0.08932,0.09819,0.10792,0.11980,0.13461,0.15360,0.16525,   &
     0.17882,0.19482,0.21397,0.23728,0.26629,0.27999,0.29517,0.31209,   &
     0.33107,0.35251,0.36430,0.37691,0.39043,0.40495,0.42059,0.43749,   &
     0.44646,0.45580,0.46555,0.47572,0.48634,0.49745,0.50908,0.52126,   &
     0.53405,0.54747,0.56159,0.57646,0.58263,0.58893,0.59537,0.60195,   &
     0.60868,0.61557,0.62261,0.62981,0.63718,0.64472,0.65245,0.66037,   &
     0.66847,0.67678,0.68530,0.69404,0.70300,0.71220,0.72164,0.73133,   &
     0.73628,0.74129,0.74637,0.75152,0.75675,0.76204,0.76741,0.77286,   &
     0.77839,0.78399,0.78968,0.79545,0.80130,0.80724,0.81327,0.81939,   &
     0.82560,0.83191,0.83832,0.84482,0.85143,0.85814,0.86495,0.87188,   &
     0.87892,0.88607,0.89334,0.90073,0.90824,0.91588,0.92365,0.93156,   &
     0.93959,0.94777,0.95610,0.96457,0.97319,0.98196,0.99090,0.99270,   &
     0.99452,0.99634,0.99725,0.99817,0.99835,0.99853,0.99872,0.99890,   &
     0.99908,0.99927,0.99945,0.99963,0.99982,1.0000/
!
      KLO=1
      KHI=NPOINT
 1    IF(KHI-KLO.GT.1) THEN
          K=(KHI+KLO)/2
          IF(WTAB(K).GT.MAX(WTAB(1),W)) THEN
              KHI=K
          ELSE
              KLO=K
          ENDIF
          GOTO 1
      ENDIF
!
      CP=CPH2O(W,KHI,KLO)
      F=-FFH2O(W,KHI,KLO)
      L=-LH2O(W,KHI,KLO)
      ALFA=ALH2O(W,KHI,KLO)
!
      A=ADOT+(CP-K1*ALFA)/RGAScal
      B=BDOT+(L-K1*CP+K2*ALFA)/RGAScal
      C=CDOT+(CP+(F-L)/K1)/RGAScal
      Dd=DDOT-ALFA/(2.0d0*RGAScal)
!
!     WRITE(*,*) 'TAIR= ',T,'  WSA= ',W
!     WRITE(*,*) 'CPH2O(W)= ',CP
!     WRITE(*,*) 'ALFAH2O(W)= ',ALFA
!     WRITE(*,*) 'FFH2O(W)= ',F
!     WRITE(*,*) 'LH2O(W)= ',L
!
      PWVSAS_GV=A*DLOG(K1/T)+B/T+C+Dd*T+MMHGPA
      
      END FUNCTION PWVSAS_GV


!*****************************************************************************
REAL FUNCTION CPH2O(W,khi_in,klo_in)

!     Relative partial molal heat capacity of water (cal/(deg mole) in 
!     sulfuric acid solution, as a function of H2SO4 weight fraction [0;1],
!     calculated by cubic spline fitting.
!
!     Source: Giauque et al.: J. Amer. Chem. Soc. 82,62,1960.

      IMPLICIT NONE

      INTEGER :: NPOINT,I
      PARAMETER(NPOINT=109)
      REAL, DIMENSION(NPOINT) :: WTAB(NPOINT),CPHTAB(NPOINT),  &
                                 Y2(NPOINT),YWORK(NPOINT)
      REAL, INTENT(IN):: W
      INTEGER, INTENT(IN):: khi_in,klo_in
      INTEGER :: khi,klo
      REAL :: CPH
      LOGICAL :: FIRST
      DATA (WTAB(I),I=1,NPOINT)/   &
     0.00000,0.08932,0.09819,0.10792,0.11980,0.13461,0.15360,0.16525,   &
     0.17882,0.19482,0.21397,0.23728,0.26629,0.27999,0.29517,0.31209,   &
     0.33107,0.35251,0.36430,0.37691,0.39043,0.40495,0.42059,0.43749,   &
     0.44646,0.45580,0.46555,0.47572,0.48634,0.49745,0.50908,0.52126,   &
     0.53405,0.54747,0.56159,0.57646,0.58263,0.58893,0.59537,0.60195,   &
     0.60868,0.61557,0.62261,0.62981,0.63718,0.64472,0.65245,0.66037,   &
     0.66847,0.67678,0.68530,0.69404,0.70300,0.71220,0.72164,0.73133,   &
     0.73628,0.74129,0.74637,0.75152,0.75675,0.76204,0.76741,0.77286,   &
     0.77839,0.78399,0.78968,0.79545,0.80130,0.80724,0.81327,0.81939,   &
     0.82560,0.83191,0.83832,0.84482,0.85143,0.85814,0.86495,0.87188,   &
     0.87892,0.88607,0.89334,0.90073,0.90824,0.91588,0.92365,0.93156,   &
     0.93959,0.94777,0.95610,0.96457,0.97319,0.98196,0.99090,0.99270,   &
     0.99452,0.99634,0.99725,0.99817,0.99835,0.99853,0.99872,0.99890,   &
     0.99908,0.99927,0.99945,0.99963,0.99982/
      DATA (CPHTAB(I),I=1,NPOINT)/   &
      17.996, 17.896, 17.875, 17.858, 17.840, 17.820, 17.800, 17.791,   &
      17.783, 17.777, 17.771, 17.769, 17.806, 17.891, 18.057, 18.248,   &
      18.429, 18.567, 18.613, 18.640, 18.660, 18.660, 18.642, 18.592,   &
      18.544, 18.468, 18.348, 18.187, 17.995, 17.782, 17.562, 17.352,   &
      17.162, 16.993, 16.829, 16.657, 16.581, 16.497, 16.405, 16.302,   &
      16.186, 16.053, 15.901, 15.730, 15.540, 15.329, 15.101, 14.853,   &
      14.586, 14.296, 13.980, 13.638, 13.274, 12.896, 12.507, 12.111,   &
      11.911, 11.711, 11.514, 11.320, 11.130, 10.940, 10.760, 10.570,   &
      10.390, 10.200, 10.000, 9.8400, 9.7600, 9.7900, 9.9500, 10.310,   &
      10.950, 11.960, 13.370, 15.060, 16.860, 18.550, 20.000, 21.170,   &
      22.030, 22.570, 22.800, 22.750, 22.420, 21.850, 21.120, 20.280,   &
      19.360, 18.350, 17.220, 15.940, 14.490, 12.840, 10.800, 9.8000,   &
      7.8000, 3.8000,0.20000,-5.4000,-7.0000,-8.8000,-10.900,-13.500,   &
     -17.000,-22.000,-29.000,-40.000,-59.000/
      DATA FIRST/.TRUE./
      SAVE FIRST,WTAB,CPHTAB,Y2
!
      IF(FIRST) THEN
          FIRST=.FALSE.
          CALL SPLINE(WTAB,CPHTAB,NPOINT,YWORK,Y2)
      ENDIF

      if(khi_in.GT.NPOINT) then
         khi=NPOINT
         klo=NPOINT-1
      else
         khi=khi_in
         klo=klo_in
      endif

      CALL SPLINT(WTAB(khi),WTAB(klo),CPHTAB(khi),CPHTAB(klo),Y2(khi),Y2(klo),W,CPH)
      CPH2O=CPH
      
      END FUNCTION CPH2O


!*******************************************************************************
REAL FUNCTION FFH2O(W,khi,klo)

!     Relative partial molal free energy water (cal/mole) in 
!     sulfuric acid solution, as a function of H2SO4 weight fraction [0;1],
!     calculated by cubic spline fitting.

!     Source: Giauque et al.: J. Amer. Chem. Soc. 82,62,1960.

      IMPLICIT NONE

      INTEGER :: NPOINT,I
      PARAMETER(NPOINT=110)
      REAL, DIMENSION(NPOINT) :: WTAB,FFTAB,Y2,YWORK
      REAL, INTENT(IN) :: W
      INTEGER, INTENT(IN):: khi,klo
      REAL :: FF
      LOGICAL :: FIRST
      DATA (WTAB(I),I=1,NPOINT)/   &
     0.00000,0.08932,0.09819,0.10792,0.11980,0.13461,0.15360,0.16525,   &
     0.17882,0.19482,0.21397,0.23728,0.26629,0.27999,0.29517,0.31209,   &
     0.33107,0.35251,0.36430,0.37691,0.39043,0.40495,0.42059,0.43749,   &
     0.44646,0.45580,0.46555,0.47572,0.48634,0.49745,0.50908,0.52126,   &
     0.53405,0.54747,0.56159,0.57646,0.58263,0.58893,0.59537,0.60195,   &
     0.60868,0.61557,0.62261,0.62981,0.63718,0.64472,0.65245,0.66037,   &
     0.66847,0.67678,0.68530,0.69404,0.70300,0.71220,0.72164,0.73133,   &
     0.73628,0.74129,0.74637,0.75152,0.75675,0.76204,0.76741,0.77286,   &
     0.77839,0.78399,0.78968,0.79545,0.80130,0.80724,0.81327,0.81939,   &
     0.82560,0.83191,0.83832,0.84482,0.85143,0.85814,0.86495,0.87188,   &
     0.87892,0.88607,0.89334,0.90073,0.90824,0.91588,0.92365,0.93156,   &
     0.93959,0.94777,0.95610,0.96457,0.97319,0.98196,0.99090,0.99270,   &
     0.99452,0.99634,0.99725,0.99817,0.99835,0.99853,0.99872,0.99890,   &
     0.99908,0.99927,0.99945,0.99963,0.99982, 1.0000/
      DATA (FFTAB(I),I=1,NPOINT)/   &
     0.00000, 22.840, 25.810, 29.250, 33.790, 39.970, 48.690, 54.560,   &
      61.990, 71.790, 85.040, 103.70, 130.70, 145.20, 163.00, 184.50,   &
      211.50, 245.60, 266.40, 290.10, 317.40, 349.00, 385.60, 428.40,   &
      452.50, 478.80, 507.50, 538.80, 573.30, 611.60, 653.70, 700.50,   &
      752.60, 810.60, 875.60, 948.60, 980.60, 1014.3, 1049.7, 1087.1,   &
      1126.7, 1168.7, 1213.5, 1261.2, 1312.0, 1366.2, 1424.3, 1486.0,   &
      1551.8, 1622.3, 1697.8, 1778.5, 1864.9, 1956.8, 2055.8, 2162.0,   &
      2218.0, 2276.0, 2337.0, 2400.0, 2466.0, 2535.0, 2607.0, 2682.0,   &
      2760.0, 2842.0, 2928.0, 3018.0, 3111.0, 3209.0, 3311.0, 3417.0,   &
      3527.0, 3640.0, 3757.0, 3878.0, 4002.0, 4130.0, 4262.0, 4397.0,   &
      4535.0, 4678.0, 4824.0, 4973.0, 5128.0, 5287.0, 5454.0, 5630.0,   &
      5820.0, 6031.0, 6268.0, 6541.0, 6873.0, 7318.0, 8054.0, 8284.0,   &
      8579.0, 8997.0, 9295.0, 9720.0, 9831.0, 9954.0, 10092., 10248.,   &
      10423., 10618., 10838., 11099., 11460., 12014./
      DATA FIRST/.TRUE./
      SAVE FIRST,WTAB,FFTAB,Y2
!
      IF(FIRST) THEN
          FIRST=.FALSE.
          CALL SPLINE(WTAB,FFTAB,NPOINT,YWORK,Y2)
      ENDIF

      CALL SPLINT(WTAB(khi),WTAB(klo),FFTAB(khi),FFTAB(klo),Y2(khi),Y2(klo),W,FF)
      FFH2O=FF
      
      END FUNCTION FFH2O


!*******************************************************************************
REAL FUNCTION LH2O(W,khi,klo)
  
!     Relative partial molal heat content of water (cal/mole) in 
!     sulfuric acid solution, as a function of H2SO4 weight fraction [0;1],
!     calculated by cubic spline fitting.

!     Source: Giauque et al.: J. Amer. Chem. Soc. 82,62,1960.

      IMPLICIT NONE

      INTEGER :: NPOINT,I
      PARAMETER(NPOINT=110)
      REAL, DIMENSION(NPOINT) ::  WTAB,LTAB,Y2,YWORK
      REAL, INTENT(IN) :: W
      INTEGER, INTENT(IN):: khi,klo
      REAL :: L
      LOGICAL :: FIRST
      DATA (WTAB(I),I=1,NPOINT)/   &
     0.00000,0.08932,0.09819,0.10792,0.11980,0.13461,0.15360,0.16525,   &
     0.17882,0.19482,0.21397,0.23728,0.26629,0.27999,0.29517,0.31209,   &
     0.33107,0.35251,0.36430,0.37691,0.39043,0.40495,0.42059,0.43749,   &
     0.44646,0.45580,0.46555,0.47572,0.48634,0.49745,0.50908,0.52126,   &
     0.53405,0.54747,0.56159,0.57646,0.58263,0.58893,0.59537,0.60195,   &
     0.60868,0.61557,0.62261,0.62981,0.63718,0.64472,0.65245,0.66037,   &
     0.66847,0.67678,0.68530,0.69404,0.70300,0.71220,0.72164,0.73133,   &
     0.73628,0.74129,0.74637,0.75152,0.75675,0.76204,0.76741,0.77286,   &
     0.77839,0.78399,0.78968,0.79545,0.80130,0.80724,0.81327,0.81939,   &
     0.82560,0.83191,0.83832,0.84482,0.85143,0.85814,0.86495,0.87188,   &
     0.87892,0.88607,0.89334,0.90073,0.90824,0.91588,0.92365,0.93156,   &
     0.93959,0.94777,0.95610,0.96457,0.97319,0.98196,0.99090,0.99270,   &
     0.99452,0.99634,0.99725,0.99817,0.99835,0.99853,0.99872,0.99890,   &
     0.99908,0.99927,0.99945,0.99963,0.99982, 1.0000/
      DATA (LTAB(I),I=1,NPOINT)/   &
     0.00000, 5.2900, 6.1000, 7.1800, 8.7800, 11.210, 15.290, 18.680,   &
      23.700, 31.180, 42.500, 59.900, 89.200, 106.70, 128.60, 156.00,   &
      190.40, 233.80, 260.10, 290.00, 324.00, 362.50, 406.50, 456.10,   &
      483.20, 512.40, 543.60, 577.40, 613.80, 653.50, 696.70, 744.50,   &
      797.20, 855.80, 921.70, 995.70, 1028.1, 1062.3, 1098.3, 1136.4,   &
      1176.7, 1219.3, 1264.7, 1313.0, 1364.3, 1418.9, 1477.3, 1539.9,   &
      1607.2, 1679.7, 1757.9, 1842.7, 1934.8, 2035.4, 2145.5, 2267.0,   &
      2332.0, 2401.0, 2473.0, 2550.0, 2631.0, 2716.0, 2807.0, 2904.0,   &
      3007.0, 3118.0, 3238.0, 3367.0, 3507.0, 3657.0, 3821.0, 3997.0,   &
      4186.0, 4387.0, 4599.0, 4819.0, 5039.0, 5258.0, 5476.0, 5694.0,   &
      5906.0, 6103.0, 6275.0, 6434.0, 6592.0, 6743.0, 6880.0, 7008.0,   &
      7133.0, 7255.0, 7376.0, 7497.0, 7618.0, 7739.0, 7855.0, 7876.0,   &
      7905.0, 7985.0, 8110.0, 8415.0, 8515.0, 8655.0, 8835.0, 9125.0,   &
      9575.0, 10325., 11575., 13500., 15200., 16125./
      DATA FIRST/.TRUE./
      SAVE FIRST,WTAB,LTAB,Y2
!
      IF(FIRST) THEN
          FIRST=.FALSE.
          CALL SPLINE(WTAB,LTAB,NPOINT,YWORK,Y2)
      ENDIF

      CALL SPLINT(WTAB(khi),WTAB(klo),LTAB(khi),LTAB(klo),Y2(khi),Y2(klo),W,L)
      LH2O=L
      
      END FUNCTION LH2O


!*******************************************************************************
REAL FUNCTION ALH2O(W,khi_in,klo_in)

!     Relative partial molal temperature derivative of heat capacity (water) 
!     in sulfuric acid solution, (cal/deg**2), calculated by 
!     cubic spline fitting.

!     Source: Giauque et al.: J. Amer. Chem. Soc. 82,62,1960.

      IMPLICIT NONE

      INTEGER :: NPOINT,I
      PARAMETER(NPOINT=96)
      REAL, DIMENSION(NPOINT) :: WTAB,ATAB,Y2,YWORK
      REAL, INTENT(IN) :: W
      INTEGER, INTENT(IN):: khi_in,klo_in
      INTEGER :: khi,klo
      REAL :: A
      LOGICAL :: FIRST
      DATA (WTAB(I),I=1,NPOINT)/   &
     0.29517,0.31209,                                                   &
     0.33107,0.35251,0.36430,0.37691,0.39043,0.40495,0.42059,0.43749,   &
     0.44646,0.45580,0.46555,0.47572,0.48634,0.49745,0.50908,0.52126,   &
     0.53405,0.54747,0.56159,0.57646,0.58263,0.58893,0.59537,0.60195,   &
     0.60868,0.61557,0.62261,0.62981,0.63718,0.64472,0.65245,0.66037,   &
     0.66847,0.67678,0.68530,0.69404,0.70300,0.71220,0.72164,0.73133,   &
     0.73628,0.74129,0.74637,0.75152,0.75675,0.76204,0.76741,0.77286,   &
     0.77839,0.78399,0.78968,0.79545,0.80130,0.80724,0.81327,0.81939,   &
     0.82560,0.83191,0.83832,0.84482,0.85143,0.85814,0.86495,0.87188,   &
     0.87892,0.88607,0.89334,0.90073,0.90824,0.91588,0.92365,0.93156,   &
     0.93959,0.94777,0.95610,0.96457,0.97319,0.98196,0.99090,0.99270,   &
     0.99452,0.99634,0.99725,0.99817,0.99835,0.99853,0.99872,0.99890,   &
     0.99908,0.99927,0.99945,0.99963,0.99982, 1.0000/
      DATA (ATAB(I),I=1,NPOINT)/   &
      0.0190, 0.0182, 0.0180, 0.0177, 0.0174, 0.0169, 0.0167, 0.0164,   &
      0.0172, 0.0212, 0.0239, 0.0264, 0.0276, 0.0273, 0.0259, 0.0238,   &
      0.0213, 0.0190, 0.0170, 0.0155, 0.0143, 0.0133, 0.0129, 0.0124,   &
      0.0120, 0.0114, 0.0106, 0.0097, 0.0084, 0.0067, 0.0047, 0.0024,   &
     -0.0002,-0.0031,-0.0063,-0.0097,-0.0136,-0.0178,-0.0221,-0.0263,   &
     -0.0303,-0.0340,-0.0352,-0.0360,-0.0362,-0.0356,-0.0343,-0.0321,   &
     -0.0290,-0.0251,-0.0201,-0.0137,-0.0058, 0.0033, 0.0136, 0.0254,   &
      0.0388, 0.0550, 0.0738, 0.0962, 0.1198, 0.1300, 0.1208, 0.0790,   &
      0.0348, 0.0058,-0.0102,-0.0211,-0.0292,-0.0350,-0.0390,-0.0418,   &
     -0.0432,-0.0436,-0.0429,-0.0411,-0.0384,-0.0346,-0.0292,-0.0220,   &
     -0.0130,-0.0110,-0.0080,-0.0060,-0.0040,-0.0030,-0.0030,-0.0020,   &
     -0.0020,-0.0020,-0.0020,-0.0010,-0.0010, 0.0000, 0.0000, 0.0000/
      DATA FIRST/.TRUE./
      SAVE FIRST,WTAB,ATAB,Y2
!
      IF(FIRST) THEN
          FIRST=.FALSE.
          CALL SPLINE(WTAB,ATAB,NPOINT,YWORK,Y2)
      ENDIF

      if(klo_in.LT.15) then
         khi=2
         klo=1
      else
         khi=khi_in-14
         klo=klo_in-14
      endif

      CALL SPLINT(WTAB(khi),WTAB(klo),ATAB(khi),ATAB(klo),Y2(khi),Y2(klo),W,A)
      ALH2O=A
      
      END FUNCTION ALH2O

!******************************************************************************
SUBROUTINE SPLINE(X,Y,N,WORK,Y2)

!     Routine to calculate 2.nd derivatives of tabulated function
!     Y(i)=Y(Xi), to be used for cubic spline calculation.

  IMPLICIT NONE

  INTEGER N,I
  REAL,intent(in) ::  X(N),Y(N)
  REAL,intent(out) :: WORK(N),Y2(N)
  REAL :: SIG,P,QN,UN,YP1,YPN

  YP1=(Y(2)-Y(1))/(X(2)-X(1))
  YPN=(Y(N)-Y(N-1))/(X(N)-X(N-1))

  IF(YP1.GT.99.0D+30) THEN
     Y2(1)=0.0
     WORK(1)=0.0
  ELSE
     Y2(1)=-0.5D0
     WORK(1)=(3.0D0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
  ENDIF

  DO I=2,N-1
     SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
     P=SIG*Y2(I-1)+2.0D0
     Y2(I)=(SIG-1.0D0)/P
     WORK(I)=(6.0D0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) &
     &  /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*WORK(I-1))/P
  ENDDO

  IF(YPN.GT.99.0D+30) THEN
     QN=0.0
     UN=0.0
  ELSE
     QN=0.5D0
     UN=(3.0D0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
  ENDIF

  Y2(N)=(UN-QN*WORK(N-1))/(QN*Y2(N-1)+1.0D0)

  DO I=N-1,1,-1
     Y2(I)=Y2(I)*Y2(I+1)+WORK(I)
  ENDDO
  
  RETURN
END SUBROUTINE SPLINE


!******************************************************************************
     SUBROUTINE SPLINT(XAhi,XAlo,YAhi,YAlo,Y2Ahi,Y2Alo,X,Y)

!     Cubic spline calculation

      IMPLICIT NONE

      REAL, INTENT(IN) :: XAhi,XAlo,YAhi,YAlo,Y2Ahi,Y2Alo
      REAL, INTENT(IN) :: X
      REAL, INTENT(OUT) :: Y
      REAL :: H,A,B
!
      H=XAhi-XAlo
      A=(XAhi-X)/H
      B=(X-XAlo)/H
      Y=A*YAlo+B*YAhi+((A**3-A)*Y2Alo+(B**3-B)*Y2Ahi)*(H**2)/6.0d0
!

      END SUBROUTINE SPLINT
