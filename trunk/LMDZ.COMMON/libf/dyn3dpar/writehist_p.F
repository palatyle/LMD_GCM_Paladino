!
! $Id: writehist_p.F 1279 2009-12-10 09:02:56Z fairhead $
!
      subroutine writehist_p( histid, histvid, time, vcov, 
     ,                          ucov,teta,phi,q,masse,ps,phis)

#ifdef CPP_IOIPSL
! This routine needs IOIPSL
      USE ioipsl
#endif
      USE parallel_lmdz
      USE misc_mod
      USE infotrac
      USE temps_mod, ONLY: itau_dyn
      implicit none

C
C   Ecriture du fichier histoire au format IOIPSL
C
C   Appels succesifs des routines: histwrite
C
C   Entree:
C      histid: ID du fichier histoire
C      histvid:ID du fichier histoire pour les vents V (appele a disparaitre)
C      time: temps de l'ecriture
C      vcov: vents v covariants
C      ucov: vents u covariants
C      teta: temperature potentielle
C      phi : geopotentiel instantane
C      q   : traceurs
C      masse: masse
C      ps   :pression au sol
C      phis : geopotentiel au sol
C      
C
C   Sortie:
C      fileid: ID du fichier netcdf cree
C
C   L. Fairhead, LMD, 03/99
C
C =====================================================================
C
C   Declarations
#include "dimensions.h"
#include "paramet.h"
#include "comgeom.h"
#include "iniprint.h"

C
C   Arguments
C

      INTEGER histid, histvid
      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm) 
      REAL teta(ip1jmp1,llm),phi(ip1jmp1,llm)                   
      REAL ps(ip1jmp1),masse(ip1jmp1,llm)                   
      REAL phis(ip1jmp1)                  
      REAL q(ip1jmp1,llm,nqtot)
      integer time

#ifdef CPP_IOIPSL
! This routine needs IOIPSL
C   Variables locales
C
      integer iq, ii, ll
      integer ndexu(ip1jmp1*llm),ndexv(ip1jm*llm),ndex2d(ip1jmp1)
      logical ok_sync
      integer itau_w
      integer :: ijb,ije,jjn
C
C  Initialisations
C
      if (adjust) return
     
    
      ndexu = 0
      ndexv = 0
      ndex2d = 0
      ok_sync =.TRUE.
      itau_w = itau_dyn + time
C
C  Appels a histwrite pour l'ecriture des variables a sauvegarder
C
C  Vents U
C
      ijb=ij_begin
      ije=ij_end
      jjn=jj_nb
          
      call histwrite(histid, 'ucov', itau_w, ucov(ijb:ije,:), 
     .               iip1*jjn*llm, ndexu)

C
C  Vents V
C
      if (pole_sud) ije=ij_end-iip1
      if (pole_sud) jjn=jj_nb-1
      
      call histwrite(histvid, 'vcov', itau_w, vcov(ijb:ije,:), 
     .               iip1*jjn*llm, ndexv)

C
C  Temperature potentielle
C
      ijb=ij_begin
      ije=ij_end
      jjn=jj_nb

      call histwrite(histid, 'teta', itau_w, teta(ijb:ije,:), 
     .                iip1*jjn*llm, ndexu)
C
C  Geopotentiel
C
      call histwrite(histid, 'phi', itau_w, phi(ijb:ije,:), 
     .                iip1*jjn*llm, ndexu)
C
C  Traceurs
C
        DO iq=1,nqtot
          call histwrite(histid, ttext(iq), itau_w, q(ijb:ije,:,iq), 
     .                   iip1*jjn*llm, ndexu)
        enddo
C
C  Masse
C
      call histwrite(histid, 'masse', itau_w, masse(ijb:ije,1),
     .               iip1*jjn, ndex2d)
C
C  Pression au sol
C
      call histwrite(histid, 'ps', itau_w, ps(ijb:ije),
     .               iip1*jjn, ndex2d)
C
C  Geopotentiel au sol
C
      call histwrite(histid, 'phis', itau_w, phis(ijb:ije),
     .               iip1*jjn, ndex2d)
C
C  Fin
C
      if (ok_sync) then
        call histsync(histid)
        call histsync(histvid)
      endif
#else
      write(lunout,*)'writehist_p: Needs IOIPSL to function'
#endif
! #endif of #ifdef CPP_IOIPSL
      return
      end
