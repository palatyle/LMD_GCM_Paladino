      subroutine lwdiff (kdlon,kflev
     .         ,pbsur,pbtop,pdbsl
     .         ,tautotal,omegtotal,gtotal
     .         ,pemis,pfluc)

      IMPLICIT NONE
 
#include "dimensions.h"
#include "dimphys.h"
#include "dimradmars.h"
#include "callkeys.h"
#include "comcstfi.h"

#include "yomaer.h"
#include "yomlw.h"
C-----------------------------------------------------------------------
C
c  ABSORPTION ET DIFUSION HORS DE LA BANDE A 15 MICRONS :
c! 1) Dans la bande a 15 micron (CO2), les poussieres 
c! n'interviennent que comme un milieu gris non diffusif avec 
c!                      Q=Qext*(1-Omega)
c! cette bande est decoupee en deux sous bandes (indices 1 et 2)
c! pour lesquelles les parametres optiques des poussieres sont
c! identiques
c! 2)  le reste est decoupe en "nir-2" bandes : une bande qui recouvre toutes
c! les longueurs d'onde inferieures a celles de la bande a 15 microns
c! (indice 3) et nir-3 bandes pour les grandes longueurs d'onde
c! (indices 4...nir) sur chacune de ces  bandes, les poussieres
c! sont supposees diffusantes grises.
c!
C
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
      integer kdlon,kflev
      REAL PBSUR(NDLO2,nir), PBTOP(NDLO2,nir)
     S  ,  PDBSL(NDLO2,nir,KFLEV*2), PEMIS(NDLO2)

      real PFLUC(NDLO2,2,KFLEV+1)
      real tautotal(ndlon,nflev,nir)
      real omegtotal(ndlon,nflev,nir), gtotal(ndlon,nflev,nir)

C
C
C-------------------------------------------------------------------------
C
C*       0.2   LOCAL ARRAYS
C              ------------
C
C
      integer jl, jk, ndd, indd, iir, j1
      integer  j2, j2dd2, j2dd1,j2bot,j2top, j2dd
      REAL ZADJD(NDLON,NFLEV+1), ZADJU(NDLON,NFLEV+1)
     S  ,  ZDBDT(NDLON,nir,NFLEV)
     S  ,  ZDISD(NDLON,NFLEV+1), ZDISU(NDLON,NFLEV+1)
     S  ,  ZFD(NDLON), ZFDN(NDLON,NFLEV+1), ZFU(NDLON)
     S  ,  ZFUP(NDLON,NFLEV+1),ZGLAYD(NDLON),ZGLAYU(NDLON)
     S  ,  ZOMEGADD(NDLON,NFLEV*2),ZGDD(NDLON,NFLEV*2)
     S  ,  ZTAUDD(NDLON,NFLEV*2)
     S  ,  ZBHDD(NDLON,NFLEV*2+1),ZBSDD(NDLON)
     S  ,  ZZBHDD(NDLON,NFLEV*2+1),ZZBSDD(NDLON)
     S  ,  ZFAHDD(NDLON,NFLEV*2+1),ZFDHDD(NDLON,NFLEV*2+1)
     S  ,  ZZFAHDD(NDLON,NFLEV*2+1),ZZFDHDD(NDLON,NFLEV*2+1)
C
C-----------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
 100  CONTINUE
C
C*         1.1     INITIALIZE LAYER CONTRIBUTIONS
C                  ------------------------------
C
 110  CONTINUE
C

      do jl = 1 , kdlon
        do jk = 1 , nlaylte
          PFLUC(jl,1,jk) = 0.
          PFLUC(jl,2,jk) = 0.
        enddo
      enddo

      DO 112 JK = 1 , nlaylte+1
        DO 111 JL = 1 , KDLON
          ZADJD(JL,JK) = 0.
          ZADJU(JL,JK) = 0.
          ZDISD(JL,JK) = 0.
          ZDISU(JL,JK) = 0.
 111    CONTINUE
 112  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         2.      VERTICAL INTEGRATION
C                  --------------------
C
C     ------------------------------------------------------------------
C
C
C  ==================================================================
C*         2.0     contribution des bandes "hors co2"
C  ==================================================================
C
 200  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.0.1   preparation des Planck a chaque hauteur
C                  ----------------------------------
C
c!
c! le nombre de couche pour la diffusion sera le nombre de layer * 2
c! soit NDD=nlaylte*2, donc la taille du vecteur des Planck sera
c! nlaylte*2 + 1. la taille des vecteurs omega / g / tau sera
c! par contre nlaylte*2 (voir dans FLUSV.F).
c!
      NDD=nlaylte*2
      DO indd=1,ndd+1
                                            do jl=1,kdlon
         ZFAHDD(jl,indd)=0.
         ZFDHDD(jl,indd)=0.
         ZBHDD(jl,indd)=0.
                                            enddo
      ENDDO
                                            do jl=1,kdlon
      ZBSDD(jl)=0.
                                            enddo
c!
c! boucle sur les  bandes hors CO2
c!
      DO 10001 iir=3,nir
c!
                                            do jl=1,kdlon
        ZZBHDD(JL,1)=PBTOP(JL,iir)/pi
                                            enddo
        DO J1=2,NDD+1
                                            do jl=1,kdlon
           ZZBHDD(JL,J1)=
     &     ZZBHDD(JL,J1-1)-PDBSL(JL,iir,NDD-J1+2)/pi
                                            enddo
        ENDDO
                                            do jl=1,kdlon
        ZZBSDD(JL)=PBSUR(JL,iir)/pi
                                            enddo

C
C     ------------------------------------------------------------------
C
C*         2.0.2   preparation des coefficients de diffusion
C                  -----------------------------------------
C
c! les omega, g, tau ... boucle de bas en haut
        DO J2=1,nlaylte-1
          J2DD2=(nlaylte-J2+1)*2
          J2DD1=J2DD2-1
          J2BOT=3*J2-2
          J2TOP=3*J2+1
          do jl=1,kdlon
            ZTAUDD(JL,J2DD1)=tautotal(jl,J2,iir)*0.5
            ZTAUDD(JL,J2DD2)=ZTAUDD(JL,J2DD1)
            ZOMEGADD(JL,J2DD1)=omegtotal(jl,J2,iir)
            ZOMEGADD(JL,J2DD2)=omegtotal(jl,J2,iir)
            ZGDD(JL,J2DD1)=gtotal(jl,J2,iir)
            ZGDD(JL,J2DD2)=gtotal(jl,J2,iir)
          enddo
        ENDDO
        J2=nlaylte
        J2DD2=2
        J2DD1=1
        J2BOT=3*J2-2
                                            do jl=1,kdlon
        ZTAUDD(JL,J2DD1)= tautotal(jl,J2,iir)*0.5
        ZTAUDD(JL,J2DD2)= tautotal(jl,J2,iir)*0.5
        ZOMEGADD(JL,J2DD1)= omegtotal(jl,J2,iir)
        ZOMEGADD(JL,J2DD2)= omegtotal(jl,J2,iir)
        ZGDD(JL,J2DD1)= gtotal(jl,J2,iir)
        ZGDD(JL,J2DD2)= gtotal(jl,J2,iir)
                                            enddo
C
C     ------------------------------------------------------------------
C
C*         2.0.3   calcul de la diffusion
C                  ----------------------
C

c-----------------------------------------------------------------------
        CALL flusv(KDLON,0
     &  ,NDD,ZOMEGADD,ZGDD,ZTAUDD,PEMIS
     &  ,ZZBHDD,ZZBSDD
     &  ,ZZFAHDD,ZZFDHDD)
c!
c!  Cumul des flux sur le spectre hors bande du CO2
c!
        DO indd=1,ndd+1
           do jl=1,kdlon
             ZFAHDD(jl,indd)=ZFAHDD(jl,indd)+ZZFAHDD(jl,indd)
             ZFDHDD(jl,indd)=ZFDHDD(jl,indd)+ZZFDHDD(jl,indd)
           enddo
        ENDDO
10001 CONTINUE

      DO J2=1,nlaylte+1
        J2DD=(nlaylte-J2+1)*2+1
        do jl=1,kdlon
          PFLUC(JL,1,J2)=PFLUC(JL,1,J2)+ZFAHDD(JL,J2DD)
          PFLUC(JL,2,J2)=PFLUC(JL,2,J2)-ZFDHDD(JL,J2DD)
        enddo
      ENDDO


      END
