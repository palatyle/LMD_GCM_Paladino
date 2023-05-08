      SUBROUTINE bernoui_p (ngrid,nlay,pphi,pecin,pbern)
      USE parallel_lmdz
      IMPLICIT NONE

c=======================================================================
c
c   Auteur:   P. Le Van
c   -------
c
c   Objet:
c   ------
c     calcul de la fonction de Bernouilli aux niveaux s  .....
c     phi  et  ecin  sont des arguments d'entree pour le s-pg .......
c          bern       est un  argument de sortie pour le s-pg  ......
c
c    fonction de Bernouilli = bern = filtre de( geopotentiel + 
c                              energ.cinet.)
c
c=======================================================================
c
c-----------------------------------------------------------------------
c   Decalrations:
c   -------------
c
#include "dimensions.h"
#include "paramet.h"
c
c   Arguments:
c   ----------
c
      INTEGER nlay,ngrid
      REAL pphi(ngrid,nlay),pecin(ngrid,nlay),pbern(ngrid,nlay)
c
c   Local:
c   ------
c
      INTEGER   ij,l,ijb,ije,jjb,jje
c
c-----------------------------------------------------------------------
c   calcul de Bernouilli:
c   ---------------------
c
      ijb=ij_begin
      ije=ij_end+iip1
      if (pole_sud) ije=ij_end

      jjb=jj_begin
      jje=jj_end+1
      if (pole_sud) jje=jj_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)                
      DO l=1,llm
    
        DO 4 ij = ijb,ije
          pbern( ij,l ) =  pphi( ij,l ) + pecin( ij,l )
   4    CONTINUE
       
       ENDDO
c$OMP END DO NOWAIT
c
c-----------------------------------------------------------------------
c   filtre:
c   -------
c

        
        CALL filtreg_p( pbern,jjb,jje, jjp1, llm, 2,1, .true., 1 )
c
c-----------------------------------------------------------------------
      
      
      RETURN
      END
