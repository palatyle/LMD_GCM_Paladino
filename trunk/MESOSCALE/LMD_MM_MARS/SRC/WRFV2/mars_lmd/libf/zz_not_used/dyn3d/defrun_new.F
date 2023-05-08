      SUBROUTINE defrun_new( tapedef, etatinit )
c
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     Auteurs :   L. Fairhead , P. Le Van  .
c      Modif C. Hourdin F. Forget VERSION MARTIENNE
c
c
c   Declarations :
c   --------------

#include "dimensions.h"
#include "paramet.h"
#include "control.h"
#include "logic.h"
#include "serre.h"
#include "comdissnew.h"
#include "sponge.h"
c
      LOGICAL  etatinit
c
c   local:
c   ------

      CHARACTER ch1*72,ch2*72,ch3*72,ch4*8
      INTEGER tapeout,tapedef

c     REAL clonn,clatt,alphaxx,alphayy
c     LOGICAL  fxyhypbb
      INTEGER ierr
      REAL clonn,clatt,grossismxx,grossismyy
      REAL dzoomxx,dzoomyy,tauxx,tauyy,temp
      LOGICAL  fxyhypbb, ysinuss


c
c  -------------------------------------------------------------------
c
c                    MODIF JUIN 2000 (zoom)
c       .........     Version  du 29/04/97       ..........
c
c   Nouveaux parametres nitergdiv,nitergrot,niterh,tetagdiv,tetagrot,
c   tetatemp   ajoutes  pour la dissipation   .
c
c   Autre parametre ajoute en fin de liste : ** fxyhypb ** 
c
c   Si fxyhypb = .TRUE. , choix de la fonction a derivee tangente hyperb.
c   Sinon , choix de fxynew  , a derivee sinusoidale  ..
c
c   ......  etatinit = . TRUE. si defrun_new  est appele dans NEWSTART
c   ETAT0_LMD  ou  LIMIT_LMD  pour l'initialisation de start.dat (dic) et
c   de limit.dat (dic)  ...........
c   Sinon  etatinit = . FALSE .
c
c   Donc etatinit = .F.  si on veut comparer les valeurs de  alphax ,
c   alphay,clon,clat, fxyhypb  lues sur  le fichier  start  avec
c   celles passees  par run.def ,  au debut du gcm, apres l'appel a 
c   lectba .  
c   Ces parametres definissant entre autres la grille et doivent etre
c   pareils et coherents , sinon il y aura  divergence du gcm .
c
c-----------------------------------------------------------------------
c   initialisations:
c   ----------------
 
      tapeout=6

c-----------------------------------------------------------------------
c  Parametres de controle du run:
c-----------------------------------------------------------------------


!Initialisation des parametres "terrestres", qui ne concernent pas
!le modele martien et ne sont donc plus lues dans "run.def"

        anneeref=0

      OPEN(tapedef,file='run.def',status='old',form='formatted'
     .     ,iostat=ierr)
      IF(ierr.EQ.0) THEN

        PRINT*
        READ (tapedef,9000) ch1,ch2,ch3
        WRITE(tapeout,9000) ch1,ch2,ch3

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    nday
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9005) 'nday       ',nday

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    day_step
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9005) 'day_step   ',day_step

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    iperiod
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9005) 'iperiod    ',iperiod

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    iconser
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9005) 'iconser    ',iconser

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    temp
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9006) 'iecri not used'

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    periodav
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9003) 'periodav   ',periodav

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    idissip
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9005) 'idissip    ',idissip

ccc  ....   P. Le Van , modif le 29/04/97 .pour la dissipation  ...
ccc
        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    lstardis
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9004) 'lstardis   ',lstardis

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    hybrid
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9004) 'hybrid ',hybrid

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    nitergdiv
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9005) 'nitergdiv  ',nitergdiv

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    nitergrot
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9005) 'nitergrot  ',nitergrot

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    niterh
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9005) 'niterh     ',niterh

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    tetagdiv
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9003) 'tetagdiv   ',tetagdiv

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    tetagrot
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9003) 'tetagrot   ',tetagrot

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    tetatemp
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9003) 'tetatemp   ',tetatemp

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    coefdis
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9003) 'coefdis    ',coefdis
c
c    ...............................................................

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    purmats
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9004) 'purmats    ',purmats

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    physic
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9004) 'physic     ',physic

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    iphysiq
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9005) 'iphysiq    ',iphysiq

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    grireg
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9004) 'grireg     ',grireg

ccc   .... P.Le Van, ajout le 03/01/96 pour l'ecriture phys ...
c
        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    ecritphy
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9003) 'ecritphy   ',ecritphy


ccc  ....   P. Le Van , ajout  le 7/03/95 .pour le zoom ...
c     .........   (  modif  le 17/04/96 )   .........
c
        if (.not.etatinit ) then 

           READ (tapedef,9001) ch1,ch4
           READ (tapedef,*)    clonn
           WRITE(tapeout,9002) ch1
           WRITE(tapeout,9003) 'clon       ',clonn
           IF( ABS(clon - clonn).GE. 0.001 )  THEN
             PRINT *,' La valeur de clon passee par run.def est '
     *       ,'differente de celle lue sur le fichier start '
             STOP
           ENDIF
c
c
           READ (tapedef,9001) ch1,ch4
           READ (tapedef,*)    clatt
           WRITE(tapeout,9002) ch1
           WRITE(tapeout,9003) 'clat       ',clatt
  
           IF( ABS(clat - clatt).GE. 0.001 )  THEN
             PRINT *,' La valeur de clat passee par run.def est '
     *       ,'differente de celle lue sur le fichier start '
             STOP
           ENDIF

           READ (tapedef,9001) ch1,ch4
           READ (tapedef,*)    grossismxx
           WRITE(tapeout,9002) ch1
           WRITE(tapeout,9003) 'grossismx  ',grossismxx

           if(grossismxx.eq.0) then  
             write(*,*)
             write(*,*)'ERREUR : dans run.def, grossismx =0'
             write(*,*)'Attention a ne pas utiliser une version de'
             write(*,*)'run.def avant le nouveau zoom LMDZ2.3 (06/2000)'
             write(*,*)'(Il faut ajouter grossismx,dzoomx,etc... a la'
             write(*,*)'place de alphax, alphay. cf. dyn3d). '
             write(*,*)
             stop
           end if

           IF( ABS(grossismx - grossismxx).GE. 0.001 )  THEN
             PRINT *,' La valeur de grossismx passee par run.def est '
     *       ,'differente de celle lue sur le fichier  start =',
     *        grossismx
             if (grossismx.eq.0) then
                  write(*,*) 'OK,Normal : c est un vieux start'
     *             , 'd avant le nouveau zoom LMDZ2.3 (06/2000)'
                 grossismx=grossismxx
             else
                   STOP
             endif
           ENDIF

           READ (tapedef,9001) ch1,ch4
           READ (tapedef,*)    grossismyy
           WRITE(tapeout,9002) ch1
           WRITE(tapeout,9003) 'grossismy  ',grossismyy

           IF( ABS(grossismy - grossismyy).GE. 0.001 )  THEN
             PRINT *,' La valeur de grossismy passee par run.def est '
     *       ,'differente de celle lue sur le fichier  start =',
     *        grossismy
             if (grossismy.eq.0) then
                  write(*,*) 'OK,Normal : c est un vieux start'
     *             , 'd avant le nouveau zoom LMDZ2.3 (06/2000)'
                 grossismy=grossismyy
             else
                   STOP
             endif
           ENDIF


           IF( grossismx.LT.1. )  THEN
             PRINT *,' ***  ATTENTION !! grossismx < 1 .   *** '
             STOP
           ELSE
             alphax = 1. - 1./ grossismx
           ENDIF

           IF( grossismy.LT.1. )  THEN
             PRINT *,' ***  ATTENTION !! grossismy < 1 .   *** '
             STOP
           ELSE
             alphay = 1. - 1./ grossismy
           ENDIF

           PRINT *,' '
           PRINT *,' --> In defrun: alphax alphay  ',alphax,alphay
           PRINT *,' '
c
           READ (tapedef,9001) ch1,ch4
           READ (tapedef,*)    fxyhypbb
           WRITE(tapeout,9002) ch1
           WRITE(tapeout,9004) 'fxyhypbb   ',fxyhypbb
  
           IF( .NOT.fxyhypb )  THEN
             IF( fxyhypbb )     THEN
                PRINT *,' ********  PBS DANS  DEFRUN  ******** '
                PRINT *,' *** fxyhypb lu sur le fichier start est F ',
     *          'alors  qu il est  T  sur  run.def  ***'
                STOP
             ENDIF
           ELSE
             IF( .NOT.fxyhypbb )   THEN
                PRINT *,' ********  PBS DANS  DEFRUN  ******** '
                PRINT *,' ***  fxyhypb lu sur le fichier start est T ',
     *         'alors  qu il est  F  sur  run.def  ****  '
                STOP
             ENDIF
           ENDIF
           READ (tapedef,9001) ch1,ch4
           READ (tapedef,*)    dzoomxx
           WRITE(tapeout,9002) ch1
           WRITE(tapeout,9003) 'dzoomx  ',dzoomxx

           IF( fxyhypb )  THEN
              IF( ABS(dzoomx - dzoomxx).GE. 0.001 )  THEN
                PRINT *,' La valeur de dzoomx passee par run.def est '
     *          ,'differente de celle lue sur le fichier  start '
                STOP
              ENDIF
           ENDIF

           READ (tapedef,9001) ch1,ch4
           READ (tapedef,*)    dzoomyy
           WRITE(tapeout,9002) ch1
           WRITE(tapeout,9003) 'dzoomy  ',dzoomyy

           IF( fxyhypb )  THEN
              IF( ABS(dzoomy - dzoomyy).GE. 0.001 )  THEN
                PRINT *,' La valeur de dzoomy passee par run.def est '
     *          ,'differente de celle lue sur le fichier  start '
                STOP
              ENDIF
           ENDIF

           READ (tapedef,9001) ch1,ch4
           READ (tapedef,*)    tauxx
           WRITE(tapeout,9001) ch1,'taux'
           WRITE(tapeout,*)    tauxx

           READ (tapedef,9001) ch1,ch4
           READ (tapedef,*)    tauyy
           WRITE(tapeout,9001) ch1,'tauy'
           WRITE(tapeout,*)    tauyy

           IF( fxyhypb )  THEN
              IF( ABS(taux - tauxx).GE. 0.001 )  THEN
                WRITE(6,*)' La valeur de taux passee par run.def est', 
     *             'differente de celle lue sur le fichier  start '
                CALL ABORT
              ENDIF

              IF( ABS(tauy - tauyy).GE. 0.001 )  THEN
                WRITE(6,*)' La valeur de tauy passee par run.def est',
     *          'differente de celle lue sur le fichier  start '
                CALL ABORT
              ENDIF
           ENDIF
  
        ELSE    ! Below, case when etainit=.true.

           READ (tapedef,9001) ch1,ch4
           READ (tapedef,*)    clon
           WRITE(tapeout,9002) ch1
           WRITE(tapeout,9003) 'clon       ',clon
c
           READ (tapedef,9001) ch1,ch4
           READ (tapedef,*)    clat
           WRITE(tapeout,9002) ch1
           WRITE(tapeout,9003) 'clat       ',clat

           READ (tapedef,9001) ch1,ch4
           READ (tapedef,*)    grossismx
           WRITE(tapeout,9002) ch1
           WRITE(tapeout,9003) 'grossismx  ',grossismx

           READ (tapedef,9001) ch1,ch4
           READ (tapedef,*)    grossismy
           WRITE(tapeout,9002) ch1
           WRITE(tapeout,9003) 'grossismy  ',grossismy

           IF( grossismx.LT.1. )  THEN
            PRINT *,' ***  ATTENTION !! grossismx < 1 .   *** '
            STOP
           ELSE
             alphax = 1. - 1./ grossismx
           ENDIF

           IF( grossismy.LT.1. )  THEN
             PRINT *,' ***  ATTENTION !! grossismy < 1 .   *** '
             STOP
           ELSE
             alphay = 1. - 1./ grossismy
           ENDIF

           PRINT *,' Defrun  alphax alphay  ',alphax,alphay
c
           READ (tapedef,9001) ch1,ch4
           READ (tapedef,*)    fxyhypb
           WRITE(tapeout,9002) ch1
           WRITE(tapeout,9004) 'fxyhypb  ',fxyhypb

           READ (tapedef,9001) ch1,ch4
           READ (tapedef,*)    dzoomx
           WRITE(tapeout,9002) ch1
           WRITE(tapeout,9003) 'dzoomx   ',dzoomx

           READ (tapedef,9001) ch1,ch4
           READ (tapedef,*)    dzoomy
           WRITE(tapeout,9002) ch1
           WRITE(tapeout,9003) 'dzoomy   ',dzoomy

           READ (tapedef,9001) ch1,ch4
           READ (tapedef,*)    taux
           WRITE(tapeout,9001) ch1,'taux'
           WRITE(tapeout,*)    taux

           READ (tapedef,9001) ch1,ch4
           READ (tapedef,*)    tauy
           WRITE(tapeout,9001) ch1,'tauy'
           WRITE(tapeout,*)    tauy

        END IF

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    ysinus
        WRITE(tapeout,9001) ch1,'ysinus'
        WRITE(tapeout,*)    ysinus

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    callsponge
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9004) 'callsponge ',callsponge

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    mode_sponge
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9005) 'mode_sponge',mode_sponge

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    hsponge
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9003) 'hsponge    ',hsponge

        READ (tapedef,9001) ch1,ch4
        READ (tapedef,*)    tetasponge
        WRITE(tapeout,9002) ch1
        WRITE(tapeout,9003) 'tetasponge ',tetasponge


      WRITE(tapeout,*) '-----------------------------------------------'
      WRITE(tapeout,*) ' '
      WRITE(tapeout,*) ' '
c

c       Unlike on Earth (cf LMDZ2.2) , always a regular grid on Mars :
        ysinus = .false. !Mars Mettre a jour


      WRITE(tapeout,*) '-----------------------------------------------'
      WRITE(tapeout,*) ' '
      WRITE(tapeout,*) ' '
cc
      ELSE
         write(0,*) ' WHERE IS run.def ? WE NEED IT !!!!!!!!!!!!!!!!!'
         stop
      ENDIF

c     Test sur le zoom

      if((grossismx.eq.1).and.(grossismy.eq.1)) then  
c        Pas de zoom :
         write(*,*) 'No zoom ? -> fxyhypb set to False'
     &   ,'           (It runs better that way)'
         fxyhypb = .false.
      else     
c        Avec Zoom
         if (.not.fxyhypb) stop 'With zoom, fxyhypb should be set to T 
     &in run.def for this version... -> STOP ! '     
      end if

      RETURN
c
9000  FORMAT(3(/,a72))
9001  FORMAT(/,a72,/,a9)
9002  FORMAT(/,a72)
9003  FORMAT(t3,a11,t14,f9.3)
9004  FORMAT(t3,a11,t14,l9)
9005  FORMAT(t3,a11,t14,i9)
9006  FORMAT(t3,a14)
      END
