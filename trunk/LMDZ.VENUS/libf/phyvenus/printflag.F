!
! $Header: /home/cvsroot/LMDZ4/libf/phylmd/printflag.F,v 1.1.1.1 2004/05/19 12:53:09 lmdzadmin Exp $
!
       SUBROUTINE  printflag( ok_journe,ok_instan )
c

c
c      Auteur :  P. Le Van 

       IMPLICIT NONE

       LOGICAL cycle_diurn0,soil_model0,ok_orodr0
       LOGICAL ok_orolf0,ok_gw_nonoro0
       LOGICAL ok_journe,ok_instan
       INTEGER radpas0
c
#include "clesphys.h"
#include "tabcontrol.h"
#include "YOMCST.h"
c
c
       PRINT 100
       PRINT *,' *******************************************************
     ,************'
       PRINT *,' ********   Choix  des principales  cles de la physique 
     ,   *********'
       PRINT *,' *******************************************************
     ,************'
       PRINT 100
       PRINT 10, cycle_diurne,  soil_model  
       PRINT 100

       PRINT 11, ok_orodr, ok_orolf, ok_gw_nonoro  
       PRINT 100

       PRINT 12, nbapp_rad
       PRINT 100

       PRINT 8, radpas
       PRINT 100

       PRINT 4,ok_journe,ok_instan
       PRINT 100
       PRINT 100
c
c
        cycle_diurn0  = .FALSE.
        soil_model0   = .FALSE.
        ok_orodr0     = .FALSE.
        ok_orolf0     = .FALSE.
        ok_gw_nonoro0 = .FALSE.

        IF( tabcntr0( 7 ).EQ. 1. )   cycle_diurn0 = .TRUE.
        IF( tabcntr0( 8 ).EQ. 1. )    soil_model0 = .TRUE.
        IF( tabcntr0(10 ).EQ. 1. )      ok_orodr0 = .TRUE.
        IF( tabcntr0(11 ).EQ. 1. )      ok_orolf0 = .TRUE.
        IF( tabcntr0(12 ).EQ. 1. )  ok_gw_nonoro0 = .TRUE.

        PRINT *,' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     ,$$$$$$$$$$$$$'
        PRINT 100
c

       IF( INT( tabcntr0( 6 ) ) .NE. nbapp_rad  )   THEN
        PRINT 21,  INT(tabcntr0(6)), nbapp_rad
        radpas0  = NINT( RDAY/tabcntr0(1)/INT( tabcntr0(6) ) )
        PRINT 100
        PRINT 22, radpas0, radpas
        PRINT 100
       ENDIF

       IF( cycle_diurn0.AND..NOT.cycle_diurne.OR..NOT.cycle_diurn0.AND.
     ,        cycle_diurne )     THEN
        PRINT 13, cycle_diurn0, cycle_diurne
        PRINT 100
       ENDIF

       IF( soil_model0.AND..NOT.soil_model.OR..NOT.soil_model0.AND.
     ,        soil_model )     THEN
        PRINT 14, soil_model0, soil_model
        PRINT 100
       ENDIF

       IF( ok_orodr0.AND..NOT.ok_orodr.OR..NOT.ok_orodr0.AND.
     ,        ok_orodr )     THEN
        PRINT 15, ok_orodr0, ok_orodr
        PRINT 100
       ENDIF

       IF( ok_orolf0.AND..NOT.ok_orolf.OR..NOT.ok_orolf0.AND.
     ,        ok_orolf )     THEN
        PRINT 17, ok_orolf0, ok_orolf
        PRINT 100
       ENDIF

       IF( ok_gw_nonoro0.AND..NOT.ok_gw_nonoro.OR..NOT.ok_gw_nonoro0.
     ,     AND.ok_gw_nonoro )     THEN
        PRINT 18, ok_gw_nonoro0, ok_gw_nonoro
        PRINT 100
       ENDIF

       PRINT 100
       PRINT *,' *******************************************************
     ,************'
       PRINT 100

 4    FORMAT(2x,5(1H*),'      ok_journe = ',l3,6x,' , ok_instan = ',
     , l3,14x,5(1H*) )


 7     FORMAT(2x,5(1H*),15x,'      ok_gw_nonoro   = ',l3,16x,5(1h*) )

 8     FORMAT(2x,'*****             radpas    =                      ' ,
     , i4,6x,' *****')

 10    FORMAT(2x,5(1H*),'    Cycle_diurne = ',l3,4x,', Soil_model = ',
     , l3,12x,6(1H*) )


 11    FORMAT(2x,5(1H*),', Ok_orodr = ',l3,3x,', Ok_orolf = ',l3,3x,
     ,    ' ok_gw_nonoro = ',l3,3x,5(1H*) )


 12    FORMAT(2x,'*****  Nb d appels /jour des routines de rayonn. = ' ,
     , i5,6x,' *****')

 13    FORMAT(2x,'$$$$$$$$   Attention !!  cycle_diurne  different  sur',
     , /1x,10x,' startphy = ',l3,2x,' et  run.def = ',l3)

 14    FORMAT(2x,'$$$$$$$$   Attention !!    soil_model  different  sur',
     , /1x,10x,' startphy = ',l3,2x,' et  run.def = ',l3)

 15    FORMAT(2x,'$$$$$$$$   Attention !!      ok_orodr  different  sur',
     , /1x,10x,' startphy = ',l3,2x,' et  run.def = ',l3)

 17    FORMAT(2x,'$$$$$$$$   Attention !!      ok_orolf  different  sur',
     , /1x,10x,' startphy = ',l3,2x,' et  run.def = ',l3)

 18    FORMAT(2x,'$$$$$$$$   Attention !!  ok_gw_nonoro  different  sur',
     , /1x,10x,' startphy = ',l3,2x,' et  run.def = ',l3)

 20    FORMAT(/2x,'$$$$$$$$   Attention !!    iflag_con  different  sur',
     , /1x,10x,' startphy = ',i3,2x,' et  run.def = ',i3 )

 21    FORMAT(2x,'$$$$$$$$   Attention !!     nbapp_rad  different  sur',
     , /1x,10x,' startphy = ',i3,2x,' et  run.def = ',i3 )

 22    FORMAT(2x,'$$$$$$$$   Attention !!        radpas  different  sur',
     , /1x,10x,' startphy = ',i3,2x,' et  run.def = ',i3 )

 100   FORMAT(/)

       RETURN
       END
