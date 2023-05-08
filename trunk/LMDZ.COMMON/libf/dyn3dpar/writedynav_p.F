!
! $Id: writedynav_p.F 1279 2009-12-10 09:02:56Z fairhead $
!
      subroutine writedynav_p( histid, time, vcov, 
     ,                          ucov,teta,ppk,phi,q,masse,ps,phis)

#ifdef CPP_IOIPSL
! This routine needs IOIPSL
      USE ioipsl
#endif
      USE parallel_lmdz
      USE misc_mod
      USE infotrac
      USE comconst_mod, ONLY: cpp
      USE temps_mod, ONLY: itau_dyn
      implicit none

C
C   Ecriture du fichier histoire au format IOIPSL
C
C   Appels succesifs des routines: histwrite
C
C   Entree:
C      histid: ID du fichier histoire
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

      INTEGER histid
      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm) 
      REAL teta(ip1jmp1,llm),phi(ip1jmp1,llm),ppk(ip1jmp1,llm)                  
      REAL ps(ip1jmp1),masse(ip1jmp1,llm)                   
      REAL phis(ip1jmp1)                  
      REAL q(ip1jmp1,llm,nqtot)
      integer time


#ifdef CPP_IOIPSL
! This routine needs IOIPSL
C   Variables locales
C
      integer ndex2d(iip1*jjp1),ndex3d(iip1*jjp1*llm),iq, ii, ll
      real us(ip1jmp1,llm), vs(ip1jmp1,llm)
      real tm(ip1jmp1,llm)
      REAL vnat(ip1jm,llm),unat(ip1jmp1,llm) 
      logical ok_sync
      integer itau_w
      integer :: ijb,ije,jjn
C
C  Initialisations
C
      if (adjust) return
      
      ndex3d = 0
      ndex2d = 0
      ok_sync = .TRUE.
      us = 999.999
      vs = 999.999
      tm = 999.999
      vnat = 999.999
      unat = 999.999
      itau_w = itau_dyn + time

C Passage aux composantes naturelles du vent
      call covnat_p(llm, ucov, vcov, unat, vnat)

C
C  Appels a histwrite pour l'ecriture des variables a sauvegarder
C
C  Vents U scalaire
C
      call gr_u_scal_p(llm, unat, us)
      
      ijb=ij_begin
      ije=ij_end
      jjn=jj_nb
      
      call histwrite(histid, 'u', itau_w, us(ijb:ije,:), 
     .               iip1*jjn*llm, ndex3d)
C
C  Vents V scalaire
C
      
      call gr_v_scal_p(llm, vnat, vs)
      call histwrite(histid, 'v', itau_w, vs(ijb:ije,:), 
     .               iip1*jjn*llm, ndex3d)
C
C  Temperature potentielle moyennee
C
     
      call histwrite(histid, 'theta', itau_w, teta(ijb:ije,:), 
     .                iip1*jjn*llm, ndex3d)
C
C  Temperature moyennee
C
      do ll=1,llm
        do ii = ijb, ije
          tm(ii,ll) = teta(ii,ll) * ppk(ii,ll)/cpp
        enddo
      enddo
      
      call histwrite(histid, 'temp', itau_w, tm(ijb:ije,:), 
     .                iip1*jjn*llm, ndex3d)
C
C  Geopotentiel
C
      call histwrite(histid, 'phi', itau_w, phi(ijb:ije,:), 
     .                iip1*jjn*llm, ndex3d)
C
C  Traceurs
C
        DO iq=1,nqtot
          call histwrite(histid, ttext(iq), itau_w, q(ijb:ije,:,iq), 
     .                   iip1*jjn*llm, ndex3d)
        enddo
C
C  Masse
C
       call histwrite(histid, 'masse', itau_w, masse(ijb:ije,1),
     .                iip1*jjn, ndex2d)
C
C  Pression au sol
C
       call histwrite(histid, 'ps', itau_w, ps(ijb:ije), 
     .                 iip1*jjn, ndex2d)
C
C  Geopotentiel au sol
C
       call histwrite(histid, 'phis', itau_w, phis(ijb:ije),
     .                 iip1*jjn, ndex2d)
C
C  Fin
C
      if (ok_sync) call histsync(histid)
#else
      write(lunout,*)'writedynav_p: Needs IOIPSL to function'
#endif
! #endif of #ifdef CPP_IOIPSL
      return
      end
