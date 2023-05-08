      subroutine groupe(pext,pbaru,pbarv,pbarum,pbarvm,wm)
      implicit none

c   SEE COMMENTS IN groupeun.F
c
c   sous-programme servant a fitlrer les champs de flux de masse aux
c   poles en "regroupant" les mailles 2 par 2 puis 4 par 4 etc. au fur
c   et a mesure qu'on se rapproche du pole.
c
c   en entree: pext, pbaru et pbarv
c
c   en sortie:  pbarum,pbarvm et wm.
c
c   remarque, le wm est recalcule a partir des pbaru pbarv et on n'a donc
c   pas besoin de w en entree.

#include "dimensions.h"
#include "paramet.h"
#include "comconst.h"
#include "comgeom2.h"
#include "comvert.h"


      real pbaru(iip1,jjp1,llm),pbarv(iip1,jjm,llm)
      real pext(iip1,jjp1,llm)

      real pbarum(iip1,jjp1,llm),pbarvm(iip1,jjm,llm)
      real wm(iip1,jjp1,llm)

      real zconvm(iip1,jjp1,llm),zconvmm(iip1,jjp1,llm)

      real uu

      integer i,j,l

c   Champs 1D

      call convflu(pbaru,pbarv,llm,zconvm)

c
      call scopy(ijp1llm,zconvm,1,zconvmm,1)
      call scopy(ijmllm,pbarv,1,pbarvm,1)

c
      call groupeun(jjp1,llm,zconvmm)
      call groupeun(jjm,llm,pbarvm)

c   Champs 3D

      do l=1,llm
         do j=2,jjm
            uu=pbaru(iim,j,l)
            do i=1,iim
               uu=uu+pbarvm(i,j,l)-pbarvm(i,j-1,l)-zconvmm(i,j,l)
               pbarum(i,j,l)=uu
c     zconvm(i,j,l ) =  xflu(i-1,j,l)-xflu(i,j,l)+
c    *                      yflu(i,j,l)-yflu(i,j-1,l)
            enddo
            pbarum(iip1,j,l)=pbarum(1,j,l)
         enddo
      enddo

c    integration de la convergence de masse de haut  en bas ......
      do l=1,llm
         do j=1,jjp1
            do i=1,iip1
               zconvmm(i,j,l)=zconvmm(i,j,l)
            enddo
         enddo
      enddo
      do  l = llm-1,1,-1
          do j=1,jjp1
             do i=1,iip1
                zconvmm(i,j,l)=zconvmm(i,j,l)+zconvmm(i,j,l+1)
             enddo
          enddo
      enddo

      CALL vitvert(zconvmm,wm)

      return
      end

