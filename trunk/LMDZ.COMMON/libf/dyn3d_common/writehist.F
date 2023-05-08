!
! $Id: writehist.F 1403 2010-07-01 09:02:53Z fairhead $
!
      subroutine writehist(time,vcov,ucov,teta,phi,q,masse,ps,phis)

#ifdef CPP_IOIPSL
      USE ioipsl
#endif
      USE infotrac, ONLY : nqtot, ttext
      use com_io_dyn_mod, only : histid,histvid,histuid
      USE temps_mod, ONLY: itau_dyn
      implicit none

C
C   Ecriture du fichier histoire au format IOIPSL
C
C   Appels succesifs des routines: histwrite
C
C   Entree:
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

      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm) 
      REAL teta(ip1jmp1,llm),phi(ip1jmp1,llm)                   
      REAL ps(ip1jmp1),masse(ip1jmp1,llm)                   
      REAL phis(ip1jmp1)                  
      REAL q(ip1jmp1,llm,nqtot)
      integer time


#ifdef CPP_IOIPSL
! This routine needs IOIPSL to work
C   Variables locales
C
      integer iq, ii, ll
      integer ndexu(ip1jmp1*llm),ndexv(ip1jm*llm),ndex2d(ip1jmp1)
      logical ok_sync
      integer itau_w
      REAL vnat(ip1jm,llm),unat(ip1jmp1,llm)

C
C  Initialisations
C
      ndexu = 0
      ndexv = 0
      ndex2d = 0
      ok_sync =.TRUE.
      itau_w = itau_dyn + time
!  Passage aux composantes naturelles du vent
      call covnat(llm, ucov, vcov, unat, vnat)
C
C  Appels a histwrite pour l'ecriture des variables a sauvegarder
C
C  Vents U
C
      call histwrite(histuid, 'u', itau_w, unat, 
     .               iip1*jjp1*llm, ndexu)
C
C  Vents V
C
      call histwrite(histvid, 'v', itau_w, vnat, 
     .               iip1*jjm*llm, ndexv)

C
C  Temperature potentielle
C
      call histwrite(histid, 'teta', itau_w, teta, 
     .                iip1*jjp1*llm, ndexu)
C
C  Geopotentiel
C
      call histwrite(histid, 'phi', itau_w, phi, 
     .                iip1*jjp1*llm, ndexu)
C
C  Traceurs
C
!        DO iq=1,nqtot
!          call histwrite(histid, ttext(iq), itau_w, q(:,:,iq), 
!     .                   iip1*jjp1*llm, ndexu)
!        enddo
!C
C  Masse
C
      call histwrite(histid,'masse',itau_w, masse,iip1*jjp1*llm,ndexu)
C
C  Pression au sol
C
      call histwrite(histid, 'ps', itau_w, ps, iip1*jjp1, ndex2d)
C
C  Geopotentiel au sol
C
!      call histwrite(histid, 'phis', itau_w, phis, iip1*jjp1, ndex2d)
C
C  Fin
C
      if (ok_sync) then
        call histsync(histid)
        call histsync(histvid)
        call histsync(histuid)
      endif
#else
! tell the user this routine should be run with ioipsl
      write(lunout,*)"writehist: Warning this routine should not be",
     &               " used without ioipsl"
#endif
! of #ifdef CPP_IOIPSL
      return
      end
