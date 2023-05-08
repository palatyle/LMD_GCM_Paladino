       SUBROUTINE heating(dist,rmu0,fract,falbe,sol_htg,zswup,zswdn,
     .                    icld)


c=======================================================================
c
c   Object:  Computation of the solar heating rate 
c                    SOL_HTG(klon,klev)
c
c   Arguments:
c   ----------
c
c      Input:
c      ------
c
c dist-----input-R- distance astronomique terre-soleil
c rmu0-----input-R- cosinus de l'angle zenithal
c fract----input-R- duree d'ensoleillement normalisee
c falbe----input-R- surface albedo
c icld-----input-I- calcul avec nuages.
c        p(klon,nl)    pressure (level)
c
c      Output:
c      -------
c sol_htg-----output-R- echauffement atmospherique (visible) (K/s)
c zswup-------output-R- flux solaire upward  (+ vers le haut)     (W/m2)
c zswdn-------output-R- flux solaire downward (+ vers le bas)     (W/m2)
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

      use dimphy
      use TGMDAT_MOD, ONLY: UBARI,UBARV,UBAR0
      use TGMDAT_MOD, ONLY: CSUBP,F0PI
      IMPLICIT NONE
      include "dimensions.h"

      INTEGER NLEVEL,NLAYER,NSPECV
      PARAMETER(NLAYER=llm,NLEVEL=NLAYER+1)
      PARAMETER (NSPECV=24)
c
c  ASTUCE POUR EVITER klon... EN ATTENDANT MIEUX
      INTEGER   ngrid
      PARAMETER (ngrid=(jjm-1)*iim+2)  ! = klon
c

c   Arguments:
c   ----------


      real dist, rmu0(klon), fract(klon), falbe(klon)
      integer icld

      real sol_htg(klon,klev)
      real zswup(klon,klev+1)
      real zswdn(klon,klev+1)
      
c   Local:
c   ------

      INTEGER I,J,IG,K,IPRINT,ilat,nq
 
c   COMMONS for interface with local subroutines:
c   ---------------------------------------------

      REAL  CH4(NLEVEL),XN2(NLEVEL),H2(NLEVEL),AR(NLEVEL)
      REAL  XMU(NLEVEL),GAS1(NLAYER),COLDEN(NLAYER)
      REAL FNETV(ngrid,NLEVEL),FUPV(ngrid,NLEVEL,NSPECV)  
      REAL FDV(ngrid,NLEVEL,NSPECV),FMNETV(ngrid,NLEVEL)

      COMMON /GASS/ CH4,XN2
     &              ,H2,AR
     &              ,XMU,GAS1
     &              ,COLDEN

      COMMON /FLUXvV/ FNETV,     
     &               FUPV,
     &               FDV, 
     &               FMNETV


c==================================================================

         fnetv  = 0.0
         sol_htg= 0.0
         zswup  = 0.0
         zswdn  = 0.0
c pour sorties dans gfluxv...
         iprint = 0

         DO ig=1,klon
            IF(fract(ig).LT.1.e-5) THEN
               DO j=1,nlayer
                  sol_htg(ig,j)=0.
               ENDDO
            ELSE
               ubar0=rmu0(ig)

               CALL sfluxv(iprint,ig,dist,falbe,icld)      ! #3

               do K=1,NSPECV
                 zswup(ig,:) = zswup(ig,:)+FUPV(ig,:,K)*fract(ig) ! >0 up
                 zswdn(ig,:) = zswdn(ig,:)+FDV(ig,:,K) *fract(ig) ! >0 down
               enddo
               fnetv(ig,:) = fnetv(ig,:) *fract(ig)   ! >0 up

c conversion en W/m2:
               zswup(ig,:) = 1.e-3*zswup(ig,:)
               zswdn(ig,:) = 1.e-3*zswdn(ig,:)
	       
               DO j=1,nlayer
                  sol_htg(ig,j)=                           ! K/s
     s            (fnetv(ig,j+1)-fnetv(ig,j))
     s                              /(colden(j)*csubp)
               ENDDO
            ENDIF
         ENDDO

      RETURN
 191  FORMAT(F8.2,1P10E10.2)
 192  FORMAT(a8,1P10E10.2)
      END
