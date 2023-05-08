      subroutine groupe_p(pext,pbaru,pbarv,pbarum,pbarvm,wm)
      USE parallel_lmdz
      USE comconst_mod, ONLY: ngroup
      implicit none

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
#include "comgeom2.h"

!      integer ngroup
!      parameter (ngroup=3)


      real pbaru(iip1,jjp1,llm),pbarv(iip1,jjm,llm)
      real pext(iip1,jjp1,llm)

      real pbarum(iip1,jjp1,llm),pbarvm(iip1,jjm,llm)
      real wm(iip1,jjp1,llm)

      real,save :: zconvm(iip1,jjp1,llm)
      real,save :: zconvmm(iip1,jjp1,llm)

      real uu

      integer i,j,l

      logical firstcall,groupe_ok
      save firstcall,groupe_ok
c$OMP THREADPRIVATE(firstcall,groupe_ok)

      data firstcall/.true./
      data groupe_ok/.true./

      integer ijb,ije,jjb,jje
      
      if (iim==1) then
         groupe_ok=.false.
      endif

      if (firstcall) then
         if (groupe_ok) then
           if(mod(iim,2**ngroup).ne.0) stop'probleme du nombre de point'
         endif
         firstcall=.false.
      endif

c   Champs 1D

      call convflu_p(pbaru,pbarv,llm,zconvm)

c
c      call scopy(ijp1llm,zconvm,1,zconvmm,1)
c      call scopy(ijmllm,pbarv,1,pbarvm,1)
      
      jjb=jj_begin
      jje=jj_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      do l=1,llm
        zconvmm(:,jjb:jje,l)=zconvm(:,jjb:jje,l)
      enddo
c$OMP END DO NOWAIT

      if (groupe_ok) then
         call groupeun_p(jjp1,llm,jjb,jje,zconvmm)
      endif
      
      jjb=jj_begin-1
      jje=jj_end
      if (pole_nord) jjb=jj_begin
      if (pole_sud)  jje=jj_end-1
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      do l=1,llm
        pbarvm(:,jjb:jje,l)=pbarv(:,jjb:jje,l)
      enddo
c$OMP END DO NOWAIT

      if (groupe_ok) then
         call groupeun_p(jjm,llm,jjb,jje,pbarvm)
      endif

c   Champs 3D
   
      jjb=jj_begin
      jje=jj_end
      if (pole_nord) jjb=jj_begin+1
      if (pole_sud)  jje=jj_end-1
      
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      do l=1,llm
         do j=jjb,jje
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
c$OMP END DO NOWAIT

c    integration de la convergence de masse de haut  en bas ......
   
      jjb=jj_begin
      jje=jj_end

c$OMP BARRIER
c$OMP MASTER      
      do  l = llm-1,1,-1
          do j=jjb,jje
             do i=1,iip1
                zconvmm(i,j,l)=zconvmm(i,j,l)+zconvmm(i,j,l+1)
             enddo
          enddo
      enddo

      if (.not. pole_sud) then
        zconvmm(:,jj_end+1,:)=0
cym	wm(:,jj_end+1,:)=0
      endif
      
c$OMP END MASTER
c$OMP BARRIER      

      CALL vitvert_p(zconvmm(1,1,1),wm(1,1,1))

      return
      end

