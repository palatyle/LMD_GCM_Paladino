!
! $Id: defrun.F 1403 2010-07-01 09:02:53Z fairhead $
!
c
c
      SUBROUTINE defrun( tapedef, etatinit, clesphy0 )
c
! ========================== ATTENTION =============================
! COMMENTAIRE SL : 
! NE SERT PLUS APPAREMMENT
! DONC PAS MIS A JOUR POUR L'UTILISATION AVEC LES PLANETES
! ==================================================================

      USE control_mod
      USE logic_mod, ONLY: purmats,iflag_phys,fxyhypb,ysinus
      USE serre_mod, ONLY: clon,clat,grossismx,grossismy,dzoomx,dzoomy,
     .		alphax,alphay,taux,tauy
 
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     Auteurs :   L. Fairhead , P. Le Van  .
c
c     Arguments :
c
c     tapedef   :
c     etatinit  :     = TRUE   , on ne  compare pas les valeurs des para- 
c     -metres  du zoom  avec  celles lues sur le fichier start .
c      clesphy0 :  sortie  .
c
       LOGICAL etatinit
       INTEGER tapedef

       INTEGER        longcles
       PARAMETER(     longcles = 20 )
       REAL clesphy0( longcles )
c
c   Declarations :
c   --------------
#include "dimensions.h"
#include "paramet.h"
#include "comdissnew.h"
#include "clesph0.h"
c
c
c   local:
c   ------

      CHARACTER ch1*72,ch2*72,ch3*72,ch4*12
      INTEGER   tapeout
      REAL clonn,clatt,grossismxx,grossismyy
      REAL dzoomxx,dzoomyy,tauxx,tauyy
      LOGICAL  fxyhypbb, ysinuss
      INTEGER i
      
c
c  -------------------------------------------------------------------
c
c       .........     Version  du 29/04/97       ..........
c
c   Nouveaux parametres nitergdiv,nitergrot,niterh,tetagdiv,tetagrot,
c      tetatemp   ajoutes  pour la dissipation   .
c
c   Autre parametre ajoute en fin de liste de tapedef : ** fxyhypb ** 
c
c  Si fxyhypb = .TRUE. , choix de la fonction a derivee tangente hyperb.
c    Sinon , choix de fxynew  , a derivee sinusoidale  ..
c
c   ......  etatinit = . TRUE. si defrun  est appele dans ETAT0_LMD  ou
c         LIMIT_LMD  pour l'initialisation de start.dat (dic) et
c                de limit.dat ( dic)                        ...........
c           Sinon  etatinit = . FALSE .
c
c   Donc etatinit = .F.  si on veut comparer les valeurs de  grossismx ,
c    grossismy,clon,clat, fxyhypb  lues sur  le fichier  start  avec
c   celles passees  par run.def ,  au debut du gcm, apres l'appel a 
c    lectba .  
c   Ces parmetres definissant entre autres la grille et doivent etre
c   pareils et coherents , sinon il y aura  divergence du gcm .
c
c-----------------------------------------------------------------------
c   initialisations:
c   ----------------

      tapeout = 6

c-----------------------------------------------------------------------
c  Parametres de controle du run:
c-----------------------------------------------------------------------

      OPEN( tapedef,file ='gcm.def',status='old',form='formatted')


      READ (tapedef,9000) ch1,ch2,ch3
      WRITE(tapeout,9000) ch1,ch2,ch3

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    dayref
      WRITE(tapeout,9001) ch1,'dayref'
      WRITE(tapeout,*)    dayref

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    anneeref
      WRITE(tapeout,9001) ch1,'anneeref'
      WRITE(tapeout,*)    anneeref

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    nday
      WRITE(tapeout,9001) ch1,'nday'
      WRITE(tapeout,*)    nday

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    day_step
      WRITE(tapeout,9001) ch1,'day_step'
      WRITE(tapeout,*)    day_step

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    iperiod
      WRITE(tapeout,9001) ch1,'iperiod'
      WRITE(tapeout,*)    iperiod

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    iapp_tracvl
      WRITE(tapeout,9001) ch1,'iapp_tracvl'
      WRITE(tapeout,*)    iapp_tracvl

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    iconser
      WRITE(tapeout,9001) ch1,'iconser'
      WRITE(tapeout,*)    iconser

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    iecri
      WRITE(tapeout,9001) ch1,'iecri'
      WRITE(tapeout,*)    iecri

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    periodav
      WRITE(tapeout,9001) ch1,'periodav'
      WRITE(tapeout,*)    periodav

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    dissip_period
      WRITE(tapeout,9001) ch1,'dissip_period'
      WRITE(tapeout,*)    dissip_period

ccc  ....   P. Le Van , modif le 29/04/97 .pour la dissipation  ...
ccc
      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    lstardis
      WRITE(tapeout,9001) ch1,'lstardis'
      WRITE(tapeout,*)    lstardis

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    nitergdiv
      WRITE(tapeout,9001) ch1,'nitergdiv'
      WRITE(tapeout,*)    nitergdiv

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    nitergrot
      WRITE(tapeout,9001) ch1,'nitergrot'
      WRITE(tapeout,*)    nitergrot

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    niterh
      WRITE(tapeout,9001) ch1,'niterh'
      WRITE(tapeout,*)    niterh

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    tetagdiv
      WRITE(tapeout,9001) ch1,'tetagdiv'
      WRITE(tapeout,*)    tetagdiv

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    tetagrot
      WRITE(tapeout,9001) ch1,'tetagrot'
      WRITE(tapeout,*)    tetagrot

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    tetatemp
      WRITE(tapeout,9001) ch1,'tetatemp'
      WRITE(tapeout,*)    tetatemp

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    coefdis
      WRITE(tapeout,9001) ch1,'coefdis'
      WRITE(tapeout,*)    coefdis
c
      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    purmats
      WRITE(tapeout,9001) ch1,'purmats'
      WRITE(tapeout,*)    purmats

c    ...............................................................

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    iflag_phys
      WRITE(tapeout,9001) ch1,'iflag_phys'
      WRITE(tapeout,*)    iflag_phys

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    iphysiq
      WRITE(tapeout,9001) ch1,'iphysiq'
      WRITE(tapeout,*)    iphysiq


      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    cycle_diurne
      WRITE(tapeout,9001) ch1,'cycle_diurne'
      WRITE(tapeout,*)    cycle_diurne

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    soil_model
      WRITE(tapeout,9001) ch1,'soil_model'
      WRITE(tapeout,*)    soil_model

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    new_oliq
      WRITE(tapeout,9001) ch1,'new_oliq'
      WRITE(tapeout,*)    new_oliq

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    ok_orodr
      WRITE(tapeout,9001) ch1,'ok_orodr'
      WRITE(tapeout,*)    ok_orodr

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    ok_orolf
      WRITE(tapeout,9001) ch1,'ok_orolf'
      WRITE(tapeout,*)    ok_orolf

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    ok_limitvrai
      WRITE(tapeout,9001) ch1,'ok_limitvrai'
      WRITE(tapeout,*)    ok_limitvrai

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    nbapp_rad
      WRITE(tapeout,9001) ch1,'nbapp_rad'
      WRITE(tapeout,*)    nbapp_rad

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    iflag_con
      WRITE(tapeout,9001) ch1,'iflag_con'
      WRITE(tapeout,*)    iflag_con

      DO i = 1, longcles
       clesphy0(i) = 0.
      ENDDO
                          clesphy0(1) = REAL( iflag_con )
                          clesphy0(2) = REAL( nbapp_rad )

       IF( cycle_diurne  ) clesphy0(3) =  1.
       IF(   soil_model  ) clesphy0(4) =  1.
       IF(     new_oliq  ) clesphy0(5) =  1.
       IF(     ok_orodr  ) clesphy0(6) =  1.
       IF(     ok_orolf  ) clesphy0(7) =  1.
       IF(  ok_limitvrai ) clesphy0(8) =  1.


ccc  ....   P. Le Van , ajout  le 7/03/95 .pour le zoom ...
c     .........   (  modif  le 17/04/96 )   .........
c
      IF( etatinit ) GO TO 100

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    clonn
      WRITE(tapeout,9001) ch1,'clon'
      WRITE(tapeout,*)    clonn
      IF( ABS(clon - clonn).GE. 0.001 )  THEN
       WRITE(tapeout,*) ' La valeur de clon passee par run.def est diffe
     *rente de  celle lue sur le fichier  start '
        STOP
      ENDIF
c
      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    clatt
      WRITE(tapeout,9001) ch1,'clat'
      WRITE(tapeout,*)    clatt

      IF( ABS(clat - clatt).GE. 0.001 )  THEN
       WRITE(tapeout,*) ' La valeur de clat passee par run.def est diffe
     *rente de  celle lue sur le fichier  start '
        STOP
      ENDIF

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    grossismxx
      WRITE(tapeout,9001) ch1,'grossismx'
      WRITE(tapeout,*)    grossismxx

      IF( ABS(grossismx - grossismxx).GE. 0.001 )  THEN
       WRITE(tapeout,*) ' La valeur de grossismx passee par run.def est
     , differente de celle lue sur le fichier  start '
        STOP
      ENDIF

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    grossismyy
      WRITE(tapeout,9001) ch1,'grossismy'
      WRITE(tapeout,*)    grossismyy

      IF( ABS(grossismy - grossismyy).GE. 0.001 )  THEN
       WRITE(tapeout,*) ' La valeur de grossismy passee par run.def est
     , differente de celle lue sur le fichier  start '
        STOP
      ENDIF
      
      IF( grossismx.LT.1. )  THEN
        WRITE(tapeout,*) ' ***  ATTENTION !! grossismx < 1 .   *** '
         STOP
      ELSE
         alphax = 1. - 1./ grossismx
      ENDIF


      IF( grossismy.LT.1. )  THEN
        WRITE(tapeout,*) ' ***  ATTENTION !! grossismy < 1 .   *** '
         STOP
      ELSE
         alphay = 1. - 1./ grossismy
      ENDIF

c
c    alphax et alphay sont les anciennes formulat. des grossissements
c
c
      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    fxyhypbb
      WRITE(tapeout,9001) ch1,'fxyhypbb'
      WRITE(tapeout,*)    fxyhypbb

      IF( .NOT.fxyhypb )  THEN
           IF( fxyhypbb )     THEN
            WRITE(tapeout,*) ' ********  PBS DANS  DEFRUN  ******** '
            WRITE(tapeout,*)' *** fxyhypb lu sur le fichier start est F'
     *,      '                   alors  qu il est  T  sur  run.def  ***'
              STOP
           ENDIF
      ELSE
           IF( .NOT.fxyhypbb )   THEN
            WRITE(tapeout,*) ' ********  PBS DANS  DEFRUN  ******** '
            WRITE(tapeout,*)' *** fxyhypb lu sur le fichier start est t'
     *,      '                   alors  qu il est  F  sur  run.def  ***'
              STOP
           ENDIF
      ENDIF
c
      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    dzoomxx
      WRITE(tapeout,9001) ch1,'dzoomx'
      WRITE(tapeout,*)    dzoomxx

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    dzoomyy
      WRITE(tapeout,9001) ch1,'dzoomy'
      WRITE(tapeout,*)    dzoomyy

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    tauxx
      WRITE(tapeout,9001) ch1,'taux'
      WRITE(tapeout,*)    tauxx

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    tauyy
      WRITE(tapeout,9001) ch1,'tauy'
      WRITE(tapeout,*)    tauyy

      IF( fxyhypb )  THEN

       IF( ABS(dzoomx - dzoomxx).GE. 0.001 )  THEN
        WRITE(tapeout,*)' La valeur de dzoomx passee par run.def est dif
     *ferente de celle lue sur le fichier  start '
        CALL ABORT_gcm("defrun", "", 1)
       ENDIF

       IF( ABS(dzoomy - dzoomyy).GE. 0.001 )  THEN
        WRITE(tapeout,*)' La valeur de dzoomy passee par run.def est dif
     *ferente de celle lue sur le fichier  start '
        CALL ABORT_gcm("defrun", "", 1)
       ENDIF

       IF( ABS(taux - tauxx).GE. 0.001 )  THEN
        WRITE(6,*)' La valeur de taux passee par run.def est differente
     *  de celle lue sur le fichier  start '
        CALL ABORT_gcm("defrun", "", 1)
       ENDIF

       IF( ABS(tauy - tauyy).GE. 0.001 )  THEN
        WRITE(6,*)' La valeur de tauy passee par run.def est differente
     *  de celle lue sur le fichier  start '
        CALL ABORT_gcm("defrun", "", 1)
       ENDIF

      ENDIF
      
cc
      IF( .NOT.fxyhypb  )  THEN
        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    ysinuss
        WRITE(tapeout,9001) ch1,'ysinus'
        WRITE(tapeout,*)    ysinuss


        IF( .NOT.ysinus )  THEN
           IF( ysinuss )     THEN
              WRITE(6,*) ' ********  PBS DANS  DEFRUN  ******** '
              WRITE(tapeout,*)'** ysinus lu sur le fichier start est F',
     *       ' alors  qu il est  T  sur  run.def  ***'
              STOP
           ENDIF
        ELSE
           IF( .NOT.ysinuss )   THEN
              WRITE(6,*) ' ********  PBS DANS  DEFRUN  ******** '
              WRITE(tapeout,*)'** ysinus lu sur le fichier start est T',
     *       ' alors  qu il est  F  sur  run.def  ***'
              STOP
           ENDIF
        ENDIF
      ENDIF
c
      WRITE(6,*) ' alphax alphay defrun ',alphax,alphay

      CLOSE(tapedef)

      RETURN
c   ...............................................
c
100   CONTINUE
c
      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    clon
      WRITE(tapeout,9001) ch1,'clon'
      WRITE(tapeout,*)    clon
c
      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    clat
      WRITE(tapeout,9001) ch1,'clat'
      WRITE(tapeout,*)    clat

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    grossismx
      WRITE(tapeout,9001) ch1,'grossismx'
      WRITE(tapeout,*)    grossismx

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    grossismy
      WRITE(tapeout,9001) ch1,'grossismy'
      WRITE(tapeout,*)    grossismy

      IF( grossismx.LT.1. )  THEN
        WRITE(tapeout,*) '***  ATTENTION !! grossismx < 1 .   *** '
         STOP
      ELSE
         alphax = 1. - 1./ grossismx
      ENDIF

      IF( grossismy.LT.1. )  THEN
        WRITE(tapeout,*) ' ***  ATTENTION !! grossismy < 1 .   *** '
         STOP
      ELSE
         alphay = 1. - 1./ grossismy
      ENDIF

c
      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    fxyhypb
      WRITE(tapeout,9001) ch1,'fxyhypb'
      WRITE(tapeout,*)    fxyhypb

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    dzoomx
      WRITE(tapeout,9001) ch1,'dzoomx'
      WRITE(tapeout,*)    dzoomx

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    dzoomy
      WRITE(tapeout,9001) ch1,'dzoomy'
      WRITE(tapeout,*)    dzoomy

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    taux
      WRITE(tapeout,9001) ch1,'taux'
      WRITE(tapeout,*)    taux
c
      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    tauy
      WRITE(tapeout,9001) ch1,'tauy'
      WRITE(tapeout,*)    tauy

      READ (tapedef,9001) ch1,ch4
      READ (tapedef,*)    ysinus
      WRITE(tapeout,9001) ch1,'ysinus'
      WRITE(tapeout,*)    ysinus
       
      WRITE(tapeout,*) ' alphax alphay defrun ',alphax,alphay
c
9000  FORMAT(3(/,a72))
9001  FORMAT(/,a72,/,a12)
cc
      CLOSE(tapedef)

      RETURN
      END
