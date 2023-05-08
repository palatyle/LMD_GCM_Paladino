      SUBROUTINE flumass_p(massebx,masseby, vcont, ucont, pbaru, pbarv)
      USE parallel_lmdz
      IMPLICIT NONE

c=======================================================================
c
c   Auteurs:  P. Le Van, F. Hourdin  .
c   -------
c
c   Objet:
c   ------
c
c *********************************************************************
c     .... calcul du flux de masse  aux niveaux s ......
c *********************************************************************
c   massebx,masseby,vcont et ucont sont des argum. d'entree pour le s-pg .
c       pbaru  et pbarv            sont des argum.de sortie pour le s-pg .
c
c=======================================================================


#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"

      REAL massebx( ip1jmp1,llm ),masseby( ip1jm,llm ) ,
     * vcont( ip1jm,llm ),ucont( ip1jmp1,llm ),pbaru( ip1jmp1,llm ),
     * pbarv( ip1jm,llm )

      REAL apbarun( iip1 ),apbarus( iip1 )

      REAL sairen,saireun,saires,saireus,ctn,cts,ctn0,cts0
      INTEGER  l,ij,i
      INTEGER ijb,ije
      
      EXTERNAL   SSUM
      REAL       SSUM
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)    
      DO  5 l = 1,llm

        ijb=ij_begin
        ije=ij_end+iip1
      
        if (pole_nord) ijb=ij_begin+iip1
        if (pole_sud)  ije=ij_end-iip1
        
        DO  1 ij = ijb,ije
          pbaru( ij,l ) = massebx( ij,l ) * ucont( ij,l )
   1    CONTINUE

        ijb=ij_begin-iip1
        ije=ij_end+iip1
      
        if (pole_nord) ijb=ij_begin
        if (pole_sud)  ije=ij_end-iip1
        
        DO 3 ij = ijb,ije
          pbarv( ij,l ) = masseby( ij,l ) * vcont( ij,l )
   3    CONTINUE

   5  CONTINUE
c$OMP END DO NOWAIT
c    ................................................................
c     calcul de la composante du flux de masse en x aux poles .......
c    ................................................................
c     par la resolution d'1 systeme de 2 equations .

c     la premiere equat.decrivant le calcul de la divergence en 1 point i
c     du pole,ce calcul etant itere de i=1 a i=im .
c                 c.a.d   ,
c     ( ( 0.5*pbaru(i)-0.5*pbaru(i-1) - pbarv(i))/aire(i)   =
c                                           - somme de ( pbarv(n) )/aire pole

c     l'autre equat.specifiant que la moyenne du flux de masse au pole est =0.
c     c.a.d    somme de pbaru(n)*aire locale(n) = 0.

c     on en revient ainsi a determiner la constante additive commune aux pbaru
c     qui representait pbaru(0,j,l) dans l'equat.du calcul de la diverg.au pt
c     i=1 .
c     i variant de 1 a im
c     n variant de 1 a im

      IF (pole_nord) THEN
     
        sairen = SSUM( iim,  aire(   1     ), 1 )
        saireun= SSUM( iim, aireu(   1     ), 1 )

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)    
        DO l = 1,llm
 
          ctn =  SSUM( iim, pbarv(    1     ,l),  1 )/ sairen
      
          pbaru(1,l)=pbarv(1,l) - ctn * aire(1)
        
          DO i = 2,iim
            pbaru(i,l) = pbaru(i- 1,l )    +
     *                   pbarv(i,l) - ctn * aire(i )
          ENDDO
        
          DO i = 1,iim
            apbarun(i) = aireu(    i   ) * pbaru(   i    , l)
          ENDDO
      
          ctn0 = -SSUM( iim,apbarun,1 )/saireun
        
          DO i = 1,iim
            pbaru(   i    , l) = 2. * ( pbaru(   i    , l) + ctn0 )
          ENDDO
       
          pbaru(   iip1 ,l ) = pbaru(    1    ,l )
        
        ENDDO
c$OMP END DO NOWAIT              

      ENDIF

      
      IF (pole_sud) THEN
  
        saires = SSUM( iim,  aire( ip1jm+1 ), 1 )
        saireus= SSUM( iim, aireu( ip1jm+1 ), 1 )

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)    
        DO  l = 1,llm
 
          cts =  SSUM( iim, pbarv(ip1jmi1+ 1,l),  1 )/ saires
          pbaru(ip1jm+1,l)= - pbarv(ip1jmi1+1,l) + cts * aire(ip1jm+1)
   
          DO i = 2,iim
            pbaru(i+ ip1jm,l) = pbaru(i+ip1jm-1,l)    -
     *                          pbarv(i+ip1jmi1,l)+cts*aire(i+ip1jm)
          ENDDO
        
          DO i = 1,iim
            apbarus(i) = aireu(i +ip1jm) * pbaru(i +ip1jm, l)
          ENDDO

          cts0 = -SSUM( iim,apbarus,1 )/saireus

          DO i = 1,iim
            pbaru(i+ ip1jm, l) = 2. * ( pbaru(i +ip1jm, l) + cts0 )
          ENDDO

          pbaru( ip1jmp1,l ) = pbaru( ip1jm +1,l )
       
        ENDDO
c$OMP END DO NOWAIT         
      ENDIF
      
      RETURN
      END
