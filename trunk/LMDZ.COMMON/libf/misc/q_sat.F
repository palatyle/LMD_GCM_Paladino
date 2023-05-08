!
! $Header$
!
!
!

      subroutine q_sat(np,temp,pres,qsat)
!
      IMPLICIT none
!======================================================================
! Autheur(s): Z.X. Li (LMD/CNRS)
!  reecriture vectorisee par F. Hourdin.
! Objet: calculer la vapeur d'eau saturante (formule Centre Euro.)
!======================================================================
! Arguments:
! kelvin---input-R: temperature en Kelvin
! millibar--input-R: pression en mb
!
! q_sat----output-R: vapeur d'eau saturante en kg/kg
!======================================================================
!
      integer np
      REAL temp(np),pres(np),qsat(np)
!
      REAL r2es
      PARAMETER (r2es=611.14 *18.0153/28.9644)
!
      REAL r3les, r3ies, r3es
      PARAMETER (R3LES=17.269)
      PARAMETER (R3IES=21.875)
!
      REAL r4les, r4ies, r4es
      PARAMETER (R4LES=35.86)
      PARAMETER (R4IES=7.66)
!
      REAL rtt
      PARAMETER (rtt=273.16)
!
      REAL retv
      PARAMETER (retv=28.9644/18.0153 - 1.0)

      real zqsat
      integer ip
!
!     ------------------------------------------------------------------
!
!

      do ip=1,np

!      write(*,*)'kelvin,millibar=',kelvin,millibar
!       write(*,*)'temp,pres=',temp(ip),pres(ip)
!
         IF (temp(ip) .LE. rtt) THEN
            r3es = r3ies
            r4es = r4ies
         ELSE
            r3es = r3les
            r4es = r4les
         ENDIF
!
         zqsat=r2es/pres(ip)*EXP(r3es*(temp(ip)-rtt)/(temp(ip)-r4es))
         zqsat=MIN(0.5,ZQSAT)
         zqsat=zqsat/(1.-retv *zqsat)
!
         qsat(ip)= zqsat
!      write(*,*)'qsat=',qsat(ip)

      enddo
!
      RETURN
      END

