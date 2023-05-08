!
! $Header$
!
      SUBROUTINE calcratqs ( flag_ratqs,
     I            paprs,pplay,q_seri,d_t_con,d_t_ajs
     O           ,ratqs,zpt_conv)
      USE dimphy
      IMPLICIT none
c======================================================================
c
c Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
c
c Objet: Moniteur general de la physique du modele
cAA      Modifications quant aux traceurs :
cAA                  -  uniformisation des parametrisations ds phytrac
cAA                  -  stockage des moyennes des champs necessaires
cAA                     en mode traceur off-line 
c======================================================================
c    modif   ( P. Le Van ,  12/10/98 )
c
c  Arguments:
c
c paprs---input-R-pression pour chaque inter-couche (en Pa)
c pplay---input-R-pression pour le mileu de chaque couche (en Pa)
c presnivs-input_R_pressions approximat. des milieux couches ( en PA)
cym#include "dimensions.h"
cym#include "dimphy.h"
      REAL paprs(klon,klev+1)
      REAL pplay(klon,klev)
      REAL d_t_con(klon,klev)
      REAL d_t_ajs(klon,klev)
      REAL ratqs(klon,klev)
      LOGICAL pt_conv(klon,klev)
      REAL q_seri(klon,klev)

      logical firstcall
      save firstcall
      data firstcall/.true./
c$OMP THREADPRIVATE(firstcall)

      REAL ratqsmin,ratqsmax,zx,epmax
      REAL ratqs1,ratqs2,ratqs3,ratqs4
      REAL ratqsc1,ratqsc2,ratqsc3,ratqsc4
      INTEGER i,k
      INTEGER flag_ratqs
      save ratqsmin,ratqsmax,epmax
      save ratqs1,ratqs2,ratqs3,ratqs4
      save ratqsc1,ratqsc2,ratqsc3,ratqsc4
c$OMP THREADPRIVATE(ratqsmin,ratqsmax,epmax)
c$OMP THREADPRIVATE(ratqs1,ratqs2,ratqs3,ratqs4)
c$OMP THREADPRIVATE(ratqsc1,ratqsc2,ratqsc3,ratqsc4)
      real zpt_conv(klon,klev)

      REAL zx_min
      PARAMETER (zx_min=1.0)
      REAL zx_max
      PARAMETER (zx_max=0.1)

	zpt_conv=0.
c
c Appeler le processus de condensation a grande echelle
c et le processus de precipitation
c
      if (flag_ratqs.eq.0) then

         ratqsmax=0.01
         ratqsmin=0.3

         if (firstcall) print*,'RATQS ANCIEN '
         do k=1,klev
         do i=1,klon
            zx = pplay(i,k)/paprs(i,1)
            zx = (zx_max-zx)/(zx_max-zx_min)
            zx = MIN(MAX(zx,0.0),1.0)
            zx = zx * zx * zx
            ratqs(i,k)= zx * (ratqsmax-ratqsmin) + ratqsmin
         enddo
         enddo

      else

c  On aplique un ratqs "interactif" a toutes les mailles affectees
c  par la convection ou se trouvant "sous" une maille affectee.
         do i=1,klon
            pt_conv(i,klev)=.false.
         enddo
         do k=klev-1,1,-1
            do i=1,klon
               pt_conv(i,k)=pt_conv(i,k+1).or.
     s               (abs(d_t_con(i,k))+abs(d_t_ajs(i,k))).gt.1.e-8
               if(pt_conv(i,k)) then
                  zpt_conv(i,k)=1.
               else
                  zpt_conv(i,k)=0.
               endif
            enddo
         enddo

         if (flag_ratqs.eq.1) then

            ratqsmin=0.4
            ratqsmax=0.99
            if (firstcall) print*,'RATQS INTERACTIF '
            do k=1,klev
                do i=1,klon
                   if (pt_conv(i,k)) then
                      ratqs(i,k)=0.01
     s                +1.5*0.25*(q_seri(i,1)-q_seri(i,k))/q_seri(i,k)
                      ratqs(i,k)=min(ratqs(i,k),ratqsmax)
                      ratqs(i,k)=max(ratqs(i,k),0.1)
                   else
                      ratqs(i,k)=0.01+(ratqsmin-0.01)*
     s             min((paprs(i,1)-pplay(i,k))/(paprs(i,1)-30000.),1.)
                   endif
                enddo
            enddo
         else if (flag_ratqs.eq.2) then
            do k=1,klev
                do i=1,klon
                   ratqs(i,k)=0.001+
     s             (q_seri(i,1)-q_seri(i,k))/q_seri(i,k)
                   if (pt_conv(i,k)) then
                      ratqs(i,k)=min(ratqs(i,k),ratqsmax)
                   else
                      ratqs(i,k)=min(ratqs(i,k),ratqsmin)
                   endif
                enddo
            enddo
         else
            do k=1,klev
               do i=1,klon
                  if (pplay(i,k).ge.95000.) then
                     if (pt_conv(i,k)) then
                        ratqs(i,k)=ratqsc1
                     else
                        ratqs(i,k)=ratqs1
                     endif
                  else if (pplay(i,k).ge.75000.) then
                     if (pt_conv(i,k)) then
                        ratqs(i,k)=ratqsc2
                     else
                        ratqs(i,k)=ratqs2
                     endif
                  else if (pplay(i,k).ge.50000.) then
                     if (pt_conv(i,k)) then
                        ratqs(i,k)=ratqsc3
                     else
                        ratqs(i,k)=ratqs3
                     endif
                  else
                     if (pt_conv(i,k)) then
                        ratqs(i,k)=ratqsc4
                     else
                        ratqs(i,k)=ratqs4
                     endif
                  endif
               enddo
            enddo
         endif

      endif

      firstcall=.false.

      return
      end
