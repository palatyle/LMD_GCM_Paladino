!
! $Id: calltherm.F90 1428 2010-09-13 08:43:37Z fairhead $
!
      subroutine calltherm(dtime  &
     &      ,pplay,paprs,pphi,weak_inversion  &
     &      ,u_seri,v_seri,t_seri,q_seri,zqsat,debut  &
     &      ,d_u_ajs,d_v_ajs,d_t_ajs,d_q_ajs  &
     &      ,fm_therm,entr_therm,detr_therm,zqasc,clwcon0,lmax,ratqscth,  &
     &       ratqsdiff,zqsatth,Ale_bl,Alp_bl,lalim_conv,wght_th, &
     &       zmax0,f0,zw2,fraca,ztv,zpspsk,ztla,zthl)

      USE dimphy
      implicit none
#include "dimensions.h"
!#include "dimphy.h"
#include "thermcell.h"
#include "iniprint.h"

!IM 140508
      INTEGER itap
      REAL dtime
      LOGICAL debut
      LOGICAL logexpr0, logexpr2(klon,klev), logexpr1(klon)
      REAL fact(klon)
      INTEGER nbptspb

      REAL u_seri(klon,klev),v_seri(klon,klev)
      REAL t_seri(klon,klev),q_seri(klon,klev),qmemoire(klon,klev)
      REAL weak_inversion(klon)
      REAL paprs(klon,klev+1)
      REAL pplay(klon,klev)
      REAL pphi(klon,klev)
      real zlev(klon,klev+1) 
!test: on sort lentr et a* pour alimenter KE
      REAL wght_th(klon,klev)
      INTEGER lalim_conv(klon)
      REAL zw2(klon,klev+1),fraca(klon,klev+1)

!FH Update Thermiques
      REAL d_t_ajs(klon,klev), d_q_ajs(klon,klev)
      REAL d_u_ajs(klon,klev),d_v_ajs(klon,klev)
      real fm_therm(klon,klev+1)
      real entr_therm(klon,klev),detr_therm(klon,klev)

!********************************************************
!     declarations
      LOGICAL flag_bidouille_stratocu
      real fmc_therm(klon,klev+1),zqasc(klon,klev)
      real zqla(klon,klev)
      real zqta(klon,klev)
      real ztv(klon,klev)
      real zpspsk(klon,klev)
      real ztla(klon,klev)
      real zthl(klon,klev)
      real wmax_sec(klon)
      real zmax_sec(klon)
      real f_sec(klon)
      real detrc_therm(klon,klev)
! FH WARNING : il semble que ces save ne servent a rien
!     save fmc_therm, detrc_therm
      real clwcon0(klon,klev)
      real zqsat(klon,klev)
      real zw_sec(klon,klev+1)
      integer lmix_sec(klon)
      integer lmax(klon)
      real ratqscth(klon,klev)
      real ratqsdiff(klon,klev)
      real zqsatth(klon,klev)  
!nouvelles variables pour la convection
      real Ale_bl(klon)
      real Alp_bl(klon)
      real Ale(klon)
      real Alp(klon)
!RC
      !on garde le zmax du pas de temps precedent
      real zmax0(klon), f0(klon)
!********************************************************


! variables locales
      REAL d_t_the(klon,klev), d_q_the(klon,klev)
      REAL d_u_the(klon,klev),d_v_the(klon,klev)
!
      real zfm_therm(klon,klev+1),zdt
      real zentr_therm(klon,klev),zdetr_therm(klon,klev)
! FH A VERIFIER : SAVE INUTILES
!      save zentr_therm,zfm_therm

      character (len=20) :: modname='calltherm'
      character (len=80) :: abort_message

      integer i,k
      logical, save :: first=.true.
!$OMP THREADPRIVATE(first)
!********************************************************
      if (first) then
        itap=0
        first=.false.
      endif

! Incrementer le compteur de la physique
     itap   = itap + 1

!  Modele du thermique
!  ===================
!         print*,'thermiques: WARNING on passe t au lieu de t_seri'


! On prend comme valeur initiale des thermiques la valeur du pas
! de temps precedent
         zfm_therm(:,:)=fm_therm(:,:)
         zdetr_therm(:,:)=detr_therm(:,:)
         zentr_therm(:,:)=entr_therm(:,:)

! On reinitialise les flux de masse a zero pour le cumul en
! cas de splitting
         fm_therm(:,:)=0.
         entr_therm(:,:)=0.
         detr_therm(:,:)=0.

         Ale_bl(:)=0.
         Alp_bl(:)=0.
         if (prt_level.ge.10) then
          print*,'thermV4 nsplit: ',nsplit_thermals,' weak_inversion'
         endif

!   tests sur les valeurs negatives de l'eau
         logexpr0=prt_level.ge.10
         nbptspb=0
         do k=1,klev
            do i=1,klon
! Attention teste abderr 19-03-09
!               logexpr2(i,k)=.not.q_seri(i,k).ge.0.
                logexpr2(i,k)=.not.q_seri(i,k).ge.1.e-15
               if (logexpr2(i,k)) then
                q_seri(i,k)=1.e-15
                nbptspb=nbptspb+1
               endif
!               if (logexpr0) &
!    &             print*,'WARN eau<0 avant therm i=',i,'  k=',k  &
!    &         ,' dq,q',d_q_the(i,k),q_seri(i,k)
            enddo
         enddo
         if(nbptspb.GT.0) print*,'Number of points with q_seri(i,k)<=0 ',nbptspb   

         zdt=dtime/REAL(nsplit_thermals)
         do isplit=1,nsplit_thermals

          if (iflag_thermals.eq.1) then
            CALL thermcell_2002(klon,klev,zdt   &
     &      ,pplay,paprs,pphi  &
     &      ,u_seri,v_seri,t_seri,q_seri  &
     &      ,d_u_the,d_v_the,d_t_the,d_q_the  &
     &      ,zfm_therm,zentr_therm  &
     &      ,r_aspect_thermals,30.,w2di_thermals  &
     &      ,tau_thermals,3)
          else if (iflag_thermals.eq.2) then
            CALL thermcell_sec(klon,klev,zdt  &
     &      ,pplay,paprs,pphi,zlev  &
     &      ,u_seri,v_seri,t_seri,q_seri  &
     &      ,d_u_the,d_v_the,d_t_the,d_q_the  &
     &      ,zfm_therm,zentr_therm  &
     &      ,r_aspect_thermals,30.,w2di_thermals  &
     &      ,tau_thermals,3)
          else if (iflag_thermals.eq.3) then
            CALL thermcell(klon,klev,zdt  &
     &      ,pplay,paprs,pphi  &
     &      ,u_seri,v_seri,t_seri,q_seri  &
     &      ,d_u_the,d_v_the,d_t_the,d_q_the  &
     &      ,zfm_therm,zentr_therm  &
     &      ,r_aspect_thermals,l_mix_thermals,w2di_thermals  &
     &      ,tau_thermals,3)
          else if (iflag_thermals.eq.10) then
            CALL thermcell_eau(klon,klev,zdt  &
     &      ,pplay,paprs,pphi  &
     &      ,u_seri,v_seri,t_seri,q_seri  &
     &      ,d_u_the,d_v_the,d_t_the,d_q_the  &
     &      ,zfm_therm,zentr_therm  &
     &      ,r_aspect_thermals,l_mix_thermals,w2di_thermals  &
     &      ,tau_thermals,3)
          else if (iflag_thermals.eq.11) then
              abort_message = 'cas non prevu dans calltherm'
              CALL abort_gcm (modname,abort_message,1)

!           CALL thermcell_pluie(klon,klev,zdt  &
!   &      ,pplay,paprs,pphi,zlev  &
!    &      ,u_seri,v_seri,t_seri,q_seri  &
!    &      ,d_u_the,d_v_the,d_t_the,d_q_the  &
!    &      ,zfm_therm,zentr_therm,zqla  &
!    &      ,r_aspect_thermals,l_mix_thermals,w2di_thermals  &
!    &      ,tau_thermals,3)
          else if (iflag_thermals.eq.12) then
            CALL calcul_sec(klon,klev,zdt  &
     &      ,pplay,paprs,pphi,zlev  &
     &      ,u_seri,v_seri,t_seri,q_seri  &
     &      ,zmax_sec,wmax_sec,zw_sec,lmix_sec  &
     &      ,r_aspect_thermals,l_mix_thermals,w2di_thermals  &
     &      ,tau_thermals)
          else if (iflag_thermals==13.or.iflag_thermals==14) then
            CALL thermcellV0_main(itap,klon,klev,zdt  &
     &      ,pplay,paprs,pphi,debut  &
     &      ,u_seri,v_seri,t_seri,q_seri  &
     &      ,d_u_the,d_v_the,d_t_the,d_q_the  &
     &      ,zfm_therm,zentr_therm,zdetr_therm,zqasc,zqla,lmax  &
     &      ,ratqscth,ratqsdiff,zqsatth  &
     &      ,r_aspect_thermals,l_mix_thermals  &
     &      ,tau_thermals,Ale,Alp,lalim_conv,wght_th &
     &      ,zmax0,f0,zw2,fraca)
          else if (iflag_thermals==15.or.iflag_thermals==16) then

!            print*,'THERM iflag_thermas_ed=',iflag_thermals_ed
            CALL thermcell_main(itap,klon,klev,zdt  &
     &      ,pplay,paprs,pphi,debut  &
     &      ,u_seri,v_seri,t_seri,q_seri  &
     &      ,d_u_the,d_v_the,d_t_the,d_q_the  &
     &      ,zfm_therm,zentr_therm,zdetr_therm,zqasc,zqla,lmax  &
     &      ,ratqscth,ratqsdiff,zqsatth  &
     &      ,r_aspect_thermals,l_mix_thermals &
     &      ,tau_thermals,iflag_thermals_ed,Ale,Alp,lalim_conv,wght_th &
     &      ,zmax0,f0,zw2,fraca,ztv,zpspsk &
     &      ,ztla,zthl)
           if (prt_level.gt.10) write(lunout,*)'Apres thermcell_main OK'
         else
           abort_message = 'Cas des thermiques non prevu'
           CALL abort_gcm (modname,abort_message,1)
         endif

       flag_bidouille_stratocu=iflag_thermals.eq.14.or.iflag_thermals.eq.16

      fact(:)=0.
      DO i=1,klon
       logexpr1(i)=flag_bidouille_stratocu.or.weak_inversion(i).gt.0.5
       IF(logexpr1(i)) fact(i)=1./REAL(nsplit_thermals)
      ENDDO

     DO k=1,klev
!  transformation de la derivee en tendance
            d_t_the(:,k)=d_t_the(:,k)*dtime*fact(:)
            d_u_the(:,k)=d_u_the(:,k)*dtime*fact(:)
            d_v_the(:,k)=d_v_the(:,k)*dtime*fact(:)
            d_q_the(:,k)=d_q_the(:,k)*dtime*fact(:)
            fm_therm(:,k)=fm_therm(:,k)  &
     &      +zfm_therm(:,k)*fact(:)
            entr_therm(:,k)=entr_therm(:,k)  &
     &       +zentr_therm(:,k)*fact(:)
            detr_therm(:,k)=detr_therm(:,k)  &
     &       +zdetr_therm(:,k)*fact(:)
      ENDDO
       fm_therm(:,klev+1)=0.



!  accumulation de la tendance
            d_t_ajs(:,:)=d_t_ajs(:,:)+d_t_the(:,:)
            d_u_ajs(:,:)=d_u_ajs(:,:)+d_u_the(:,:)
            d_v_ajs(:,:)=d_v_ajs(:,:)+d_v_the(:,:)
            d_q_ajs(:,:)=d_q_ajs(:,:)+d_q_the(:,:)

!  incrementation des variables meteo
            t_seri(:,:) = t_seri(:,:) + d_t_the(:,:)
            u_seri(:,:) = u_seri(:,:) + d_u_the(:,:)
            v_seri(:,:) = v_seri(:,:) + d_v_the(:,:)
            qmemoire(:,:)=q_seri(:,:)
            q_seri(:,:) = q_seri(:,:) + d_q_the(:,:)
           if (prt_level.gt.10) write(lunout,*)'Apres apres thermcell_main OK'

       DO i=1,klon
        if(prt_level.GE.10) print*,'calltherm i Alp_bl Alp Ale_bl Ale',i,Alp_bl(i),Alp(i),Ale_bl(i),Ale(i)
            fm_therm(i,klev+1)=0.
            Ale_bl(i)=Ale_bl(i)+Ale(i)/REAL(nsplit_thermals)
!            write(22,*)'ALE CALLTHERM',Ale_bl(i),Ale(i)
            Alp_bl(i)=Alp_bl(i)+Alp(i)/REAL(nsplit_thermals)
!            write(23,*)'ALP CALLTHERM',Alp_bl(i),Alp(i)
       ENDDO

!IM 060508 marche pas comme cela !!!        enddo ! isplit

!   tests sur les valeurs negatives de l'eau
         nbptspb=0
            DO k = 1, klev
            DO i = 1, klon
               logexpr2(i,k)=.not.q_seri(i,k).ge.0.
               if (logexpr2(i,k)) then
                q_seri(i,k)=1.e-15
                nbptspb=nbptspb+1
!                if (prt_level.ge.10) then
!                  print*,'WARN eau<0 apres therm i=',i,'  k=',k  &
!    &         ,' dq,q',d_q_the(i,k),q_seri(i,k),  &
!    &         'fm=',zfm_therm(i,k),'entr=',entr_therm(i,k)
                 endif
            ENDDO
            ENDDO
        IF(nbptspb.GT.0) print*,'Number of points with q_seri(i,k)<=0 ',nbptspb   
! tests sur les valeurs de la temperature
        nbptspb=0
            DO k = 1, klev
            DO i = 1, klon
               logexpr2(i,k)=t_seri(i,k).lt.50..or.t_seri(i,k).gt.370.
               if (logexpr2(i,k)) nbptspb=nbptspb+1
!              if ((t_seri(i,k).lt.50.) .or.  &
!    &              (t_seri(i,k).gt.370.)) then
!                 print*,'WARN temp apres therm i=',i,'  k=',k  &
!    &         ,' t_seri',t_seri(i,k)
!              CALL abort
!              endif
            ENDDO
            ENDDO
        IF(nbptspb.GT.0) print*,'Number of points with q_seri(i,k)<=0 ',nbptspb
         enddo ! isplit

!
!***************************************************************
!     calcul du flux ascencant conservatif
!            print*,'<<<<calcul flux ascendant conservatif'

      fmc_therm=0.
               do k=1,klev
            do i=1,klon
                  if (entr_therm(i,k).gt.0.) then
                     fmc_therm(i,k+1)=fmc_therm(i,k)+entr_therm(i,k)
                  else
                     fmc_therm(i,k+1)=fmc_therm(i,k)
                  endif
                  detrc_therm(i,k)=(fmc_therm(i,k+1)-fm_therm(i,k+1))  &
     &                 -(fmc_therm(i,k)-fm_therm(i,k))
               enddo
            enddo
      
     
!****************************************************************
!     calcul de l'humidite dans l'ascendance
!      print*,'<<<<calcul de lhumidite dans thermique'
!CR:on ne le calcule que pour le cas sec
      if (iflag_thermals.le.11) then      
      do i=1,klon
         zqasc(i,1)=q_seri(i,1)
         do k=2,klev
            if (fmc_therm(i,k+1).gt.1.e-6) then
               zqasc(i,k)=(fmc_therm(i,k)*zqasc(i,k-1)  &
     &              +entr_therm(i,k)*q_seri(i,k))/fmc_therm(i,k+1)
!CR:test on asseche le thermique
!               zqasc(i,k)=zqasc(i,k)/2.
!            else
!               zqasc(i,k)=q_seri(i,k)
            endif
         enddo
       enddo
      

!     calcul de l'eau condensee dans l'ascendance
!             print*,'<<<<calcul de leau condensee dans thermique'
             do i=1,klon
                do k=1,klev
                   clwcon0(i,k)=zqasc(i,k)-zqsat(i,k)
                   if (clwcon0(i,k).lt.0. .or.   &
     &             (fm_therm(i,k+1)+detrc_therm(i,k)).lt.1.e-6) then
                      clwcon0(i,k)=0.
                   endif
                enddo
             enddo
       else
              do i=1,klon
                do k=1,klev
                   clwcon0(i,k)=zqla(i,k)  
                   if (clwcon0(i,k).lt.0. .or.   &
     &             (fm_therm(i,k+1)+detrc_therm(i,k)).lt.1.e-6) then
                   clwcon0(i,k)=0. 
                   endif
                enddo
             enddo
       endif
!*******************************************************************    


!jyg  Protection contre les temperatures nulles
          do i=1,klon
             do k=1,klev
                if (ztla(i,k) .lt. 1.e-10) fraca(i,k) =0.
             enddo
          enddo


      return

      end
