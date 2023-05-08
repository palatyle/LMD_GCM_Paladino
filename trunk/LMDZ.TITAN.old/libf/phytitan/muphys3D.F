         subroutine muphys(ngrid,
     &   plev,play,zlev,zlay,
     &   tpt,tq,gaz1,gaz2,gaz3,
     &   nmicro,ptimestep,
     &   pmu0,pfract,     
*    Sorties diagnostiques
     &   flxesp_i,
     &   tau_drop,tau_aer,
     &   solesp,prec)

c
c
c  CETTE NOUVELLE ROUTINE DE MICROPHYSIQUE GERE 
c  LA MICROPHYSIQUE DES AEROSOLS ET CELLE DES NUAGES 
c  EN UN SEUL APPEL...
c
c  TOUS LES TRACEURS AEROSOL+NOYAUX+2*GLACES SONT CONTENUS DANS
c  LE TABLEAU TQ ET LEUR TENDANCES DANS TDQ. CES TABLEAUX
c  SONT COMPOSES DE TROIS PARTIES:
c
c       TQ(         1,    nmicro/4      pour les aerosols
c        +   nmicro/4+1,  2*nmicro/4      pour les noyaux  
c        + 2*nmicro/4+1,  3*nmicro/4      pour la glace 1
c        + 3*nmicro/4+1,    nmicro   )    pour la glace 2
c
c  pour les aerosols, les noyaux, la glace 1 et la glace 2. la separation
c  puis la concatenation de fait juste avant l'appel aux routines
c  dans lesquels la separation est necessaire. 
c
c                                            _______
c                                            |     |  ____
c  EX:         -------->    QAER()    ---->  | M   |      \
c            /                               | I   |       \
c           /                                | C   |        \
c      TQ()  ---------->    QGLACE1()  --->  | R   |  ------  TDQ()
c           \                                | O   |        /
c            \                               | P   |       /
c              -------->    QGLACE2() ---->  | H   |  ----/
c              \                             | Y   |     /
c               \                            |     |    /
c                 ----->    QNOYAUX()        |     |  _/
c                                            |     |
c                                            -------
c
c
c
c  DANS LA MICROPHYSIQUE, SE TROUVENT IMBRIQUES LES PROCESSUS SUIVANTS 
c
c
c
c
c
c   NUCLEATION/SEDIMENTATION   ---> n_ethane.F / n_methane.F
c
c   MICROPHYSIQUE AEROSOLS     ---> brume.F 
c
c   SEDIMENTATION DES GOUTTES  ---> snuages.F 
c
c
c
c
c
c
c------------------------------------------------------
         use dimphy
c         use radcommon_h, only : volume,rayon,vrat,drayon,dvolume
         USE geometry_mod,  only: latitude ! in radians

         IMPLICIT NONE
#include "dimensions.h"
#include "microtab.h"
#include "varmuphy.h"
#include "clesphys.h"
#include "itemps.h"

         integer ngrid

         integer iq,nmicro
         real    ptimestep
         real    pdpsrf(ngrid)

c a la place de radcommon_h:
         common/part/vaer,raer,vrat,draer,dvaer
         real   vaer(nrad),raer(nrad),vrat,
     &          draer(nrad),dvaer(nrad)

c*************************************
c declaration des variables internes *
c*************************************
 
c      sources     *
c------------------*

         REAL plev(ngrid,klev+1)
         REAL play(ngrid,klev)
         REAL zlev(ngrid,klev+1)
         REAL zlay(ngrid,klev) 
         REAL pu(ngrid),pv(ngrid)
         REAL pmu0(ngrid),pfract(ngrid)
         REAL tpt(ngrid,klev)
         REAL tq(ngrid,klev,nmicro), 
     &        gaz1(ngrid,klev),
     &        gaz2(ngrid,klev),
     &        gaz3(ngrid,klev)

c       OUTPUT !!!!! c
c  note : gaz1,...,gazN sont aussi des outputs, ils ont modifié tout au long de muphys.
c--------------------c
         REAL pdq(ngrid,klev,nmicro)
         REAL flxesp_i(ngrid,klev,3)    ! flx esp GLACE
         REAL solesp(ngrid,klev,3)      ! tx prod glace (puit/source)
         REAL tau_drop(ngrid,klev)
         REAL tau_aer(ngrid,klev,nrad)
         REAL prec(ngrid,5)

c       LOCAL
c--------------------c
         real  q(ngrid,klev,nmicro)
         REAL taused(klev,nrad)
         integer jsup,jinf,h,jalt,ihor,k,im1

c    microphysique    *
c---------------------*
         real c(klev,nrad),  cni(klev,nrad)
         real c1i(klev,nrad),c2i(klev,nrad),c3i(klev,nrad)
         real gazc1(klev),gazc2(klev),gazc3(klev)
         real ddt
 
         real vcl,nuc,r,xgsn,xmsn
         real zz,effg,xlog,rapport

         integer IPREM,i,n,j,l
         integer ibid
         save IPREM 
         data IPREM/0/


c         real ttq(ngrid,klev,nmicro,2)
c         real tttq(ngrid,klev,nmicro,2)


c      			**************************
c        		INITIALISATION DE TABLEAUX
c    		 	**************************
c                         A NE FAIRE QU'UNE FOIS
c     		        **************************
c
         IF (IPREM.eq.0) THEN
 
             IF (ngrid.ne.klon) THEN
               print*,"aLeRte :"
               print*,"microfi, mais ngrid.ne.klon"
               print*,ngrid,klon
               stop "je m'arrete... (muphys3D)"
             ENDIF

c initialisation des constantes de la microphysique :
c ----------------------------------------------
           call inimphycst()
           

c initialisation des c(z,r), c1i(z,r), c2i(z,r)
c ----------------------------------------------
 
           do i=1,nmicro 
             do n=1,ngrid
               do j=klev,1,-1
                 q(n,klev+1-j,i)=tq(n,j,i)                           ! glaces...
                 if(i.le.2*nrad) q(n,klev+1-j,i)=tq(n,j,i)*tcorrect  ! noyaux+aerosol 
                 pdq(n,j,i)=0.
               enddo 
             enddo 
           enddo
c ici, les tableaux definissant la structure des aerosols sont
c remplis: rf,df(nmicro),r(nmicro),v(nmicro)......
           call rdf()
c   ici on recopie la grille dans un common specifique a la microfi... 
c          v_e    = volume
c          r_e    = rayon
c          vrat_e = vrat
c          dr_e   = drayon
c          dv_e   = dvolume
           v_e    = vaer
           r_e    = raer
           vrat_e = vrat
           dr_e   = draer
           dv_e   = dvaer 
c
         ELSE   
  
c les tq() doivent etre en nombre d'aerosols / cases 

           do j=1,klev                        ! j de 1 a 119 
             do n=1,ngrid
               do  i=1,nmicro
                 q(n,j,i)=tq(n,klev-j+1,i)
                 pdq(n,j,i)=0.
               enddo
             enddo
           enddo

         ENDIF  ! FIN IPREM


c-----------------------------------------------------
c      !! La premiere fois, on ne passe pas par 
c      !! q--->c et par pg3.F
c      !! on passe directement au remplissage c-->q
         IF (IPREM.eq.0) goto 102 
c-----------------------------------------------------

c****************************************
c                                       *
c         ADAPTATION GCM > micro        *
c                                       *
c****************************************


c correpondance des couches / sens GCM > microphysique 
c-----------------------------------------------------

c***************************************************************
         do IHOR=1,NGRID ! GRANDE BOUCLE HORIZONTALE / SEPARATION DES COLONNES

         if (IHOR.eq.1) then 
           im1=1
         else
           im1=IHOR-1
         endif

c***************************************************************
c On refait les calculs si on est au premier point
c         OU            si on change de latitude
c         OU            si on calcule la microfi en 3D
c***************************************************************
        if((IHOR.eq.1)    
     & .or.(latitude(IHOR).ne.latitude(im1))
     & .or.(microfi.eq.2)) then
c***************************************************************
  
c  Ici, on initialise la grille verticale et les 
c  variables communes aux routines de microphysique.
c*******************************************************
           call inimuphy(ihor,plev(ihor,:),play(ihor,:),
     &                   zlev(ihor,:),zlay(ihor,:),
     &                   tpt(ihor,:))

c  Ici, on scinde les tableaux aerosols et glaces 
c*******************************************************

         if (clouds.eq.0) then
           if(nrad .ne. nmicro) then 
             print*,"aLeRte : nrad != nmicro"
             print*,'nmicro= ',nmicro 
             print*,'nrad= ',nrad 
             stop "je m'arrete..." 
           endif
         else
           if(nrad .ne. nmicro/ntype) then 
             print*,"aLeRte : nrad != nmicro/ntype"
             print*,'nmicro= ',nmicro 
             print*,'ntype=',ntype
             print*,'nmicro/ntype= ',nmicro/ntype
             print*,'nrad= ',nrad 
             stop "je m'arrete..." 
           endif
         endif

           do i=1,nrad
             do j=1,klev
               c(j,i)  =q(IHOR,j,i       )/dzb(j) ! concentration aerosols/m^3
               if (clouds.eq.1) then
                 cni(j,i)=q(IHOR,j,i+  nrad)/dzb(j) ! concentration noyaux /m^3
                 c1i(j,i)=q(IHOR,j,i+2*nrad)/dzb(j) ! concentration volume glace/m^3
                 c2i(j,i)=q(IHOR,j,i+3*nrad)/dzb(j) ! concentration volume glace /m^3
                 c3i(j,i)=q(IHOR,j,i+4*nrad)/dzb(j) ! concentration volume glace /m^3
               endif
             enddo
           enddo 
           if (clouds.eq.1) then
             do j=1,klev
               gazc1(j)  =gaz1(IHOR,klev-j+1)   ! fraction molaire CH4
               gazc2(j)  =gaz2(IHOR,klev-j+1)   ! fraction molaire C2H6
               gazc3(j)  =gaz3(IHOR,klev-j+1)   ! fraction molaire C2H2
             enddo
           endif

c****************************************
c 
c  FIN DE LA PREPARATION:
c
c
c  ON  APPELLE LES MODELES MICROPHYSIQUES
c
c         - brume (coagulation + sedimentation)
c         - nuages (nucleation + condensation)
c         - sedimentation nuages
c
c
c****************************************


           do j=1,klev
             solesp(ihor,klev+1-j,:) = 0.
             do i=1,nrad
               tau_aer(ihor,klev+1-j,i)=0.
             enddo
           enddo

           ddt=ptimestep

* concerne les aerosols (tableau c):
c           call begintime(tt0)
           call brume(ngrid,c,ddt,klev,nrad,taused,ihor,
     &                pmu0(ihor),pfract(ihor),
     &                prec)
c           call endtime(tt0,tt1)
c           tthaze=tthaze+tt1

* concerne aerosols +  gouttes (tableaux c,cn,c1,c2):

           do j=1,klev
             do i=1,nrad
               tau_aer(ihor,klev+1-j,i)=tau_aer(ihor,klev+1-j,i)
     &                                 +taused(j,i)
             enddo
           enddo

           IF (clouds.eq.1) THEN
          
c             do j=1,klev
c             do i=1,nrad
c               ttq(ihor,klev-j+1,i,1) = c(j,i)*dzb(j)
c               ttq(ihor,klev-j+1,i+nrad,1) = cni(j,i)*dzb(j)
c               ttq(ihor,klev-j+1,i+2*nrad,1) = c1i(j,i)*dzb(j)
c               ttq(ihor,klev-j+1,i+3*nrad,1) = c2i(j,i)*dzb(j)
c               ttq(ihor,klev-j+1,i+4*nrad,1) = c3i(j,i)*dzb(j)
c             enddo
c             enddo

c             call begintime(tt0)
             call cnuages(c,c1i,c2i,c3i,cni,
     &                      gazc1,gazc2,gazc3,ddt)
c             call endtime(tt0,tt1)
c             ttcclds=ttcclds+tt1

c  verification des valeurs de c,cni,c1i,c2i et c3i.
c  Lorsque l'on vide completement une case, on peut avoir des chiffres negatifs :s
           do j=1,klev
             do i=1,nrad
               c(j,i)= MAX(c(j,i),0.)
               cni(j,i)= MAX(cni(j,i),0.)
               c1i(j,i)= MAX(c1i(j,i),0.)
               c2i(j,i)= MAX(c2i(j,i),0.)
               c3i(j,i)= MAX(c3i(j,i),0.)
c               ttq(ihor,klev-j+1,i,2) = c(j,i)*dzb(j)
c               ttq(ihor,klev-j+1,i+nrad,2) = cni(j,i)*dzb(j)
c               ttq(ihor,klev-j+1,i+2*nrad,2) = c1i(j,i)*dzb(j)
c               ttq(ihor,klev-j+1,i+3*nrad,2) = c2i(j,i)*dzb(j)
c               ttq(ihor,klev-j+1,i+4*nrad,2) = c3i(j,i)*dzb(j)
             enddo
           enddo
     

           do j=1,klev
             do i=1,nrad
!              solesp en m3/m3 pour passer en m3/m2 il faut faire :
!              (c1i(j,i)*dzb(j) -q(IHOR,j,i+2*nrad))
               solesp(ihor,klev+1-j,1)=solesp(ihor,klev+1-j,1) +
     &         (c1i(j,i)-q(IHOR,j,i+2*nrad)/dzb(j))
               solesp(ihor,klev+1-j,2)=solesp(ihor,klev+1-j,2) +
     &         (c2i(j,i)-q(IHOR,j,i+3*nrad)/dzb(j))
               solesp(ihor,klev+1-j,3)=solesp(ihor,klev+1-j,3) +
     &         (c3i(j,i)-q(IHOR,j,i+4*nrad)/dzb(j))
             enddo
           enddo

* concerne les gouttes (tableaux cn,c1,c2):

c             do j=1,klev
c             do i=1,nrad
c               tttq(ihor,klev-j+1,i,1) = c(j,i)*dzb(j)
c               tttq(ihor,klev-j+1,i+nrad,1) = cni(j,i)*dzb(j)
c               tttq(ihor,klev-j+1,i+2*nrad,1) = c1i(j,i)*dzb(j)
c               tttq(ihor,klev-j+1,i+3*nrad,1) = c2i(j,i)*dzb(j)
c               tttq(ihor,klev-j+1,i+4*nrad,1) = c3i(j,i)*dzb(j)
c             enddo
c             enddo

c             call begintime(tt0)
             call snuages(ngrid,cni,c1i,c2i,c3i,c,ddt,
     &                    klev,nrad,ihor,
     &                    flxesp_i,tau_drop,prec)
c             call endtime(tt0,tt1)
c             ttsclds=ttsclds+tt1

c           do j=1,klev
c             do i=1,nrad
c               tttq(ihor,klev-j+1,i,2) = c(j,i)*dzb(j)
c               tttq(ihor,klev-j+1,i+nrad,2) = cni(j,i)*dzb(j)
c               tttq(ihor,klev-j+1,i+2*nrad,2) = c1i(j,i)*dzb(j)
c               tttq(ihor,klev-j+1,i+3*nrad,2) = c2i(j,i)*dzb(j)
c               tttq(ihor,klev-j+1,i+4*nrad,2) = c3i(j,i)*dzb(j)
c             enddo
c           enddo
           ENDIF    ! flag nuages :)


* on recompose le tableau de traceurs ici.
*--------------------------------------------------
 
           do i=1,nrad
             do j=1,klev
               q(IHOR,j,i)        =   c(j,i)*dzb(j) ! nombre  aerosols /m^2
               if (clouds.eq.1) then
                 q(IHOR,j,i+  nrad) = cni(j,i)*dzb(j) ! nombre noyaux /m^2
                 q(IHOR,j,i+2*nrad) = c1i(j,i)*dzb(j) ! concentration volume glace/m^2
                 q(IHOR,j,i+3*nrad) = c2i(j,i)*dzb(j) ! concentration volume glace/m^2
                 q(IHOR,j,i+4*nrad) = c3i(j,i)*dzb(j) ! concentration volume glace/m^2
               endif
             enddo
           enddo 
           if (clouds.eq.1) then
             do j=1,klev
               gaz1(IHOR,klev-j+1) = gazc1(j)   ! fraction molaire 
               gaz2(IHOR,klev-j+1) = gazc2(j)   ! fraction molaire
               gaz3(IHOR,klev-j+1) = gazc3(j)   ! fraction molaire
             enddo
           endif

c***************************************************************
         else      ! same latitude, we don't do calculations again
           q(ihor,:,:)       = q(im1,:,:)
           tau_aer(ihor,:,:) = tau_aer(im1,:,:)
           prec(ihor,:)    = prec(im1,:)
           if (clouds.eq.1) then
             solesp(ihor,:,:)   = solesp(im1,:,:)
             flxesp_i(ihor,:,:) = flxesp_i(im1,:,:)
             tau_drop(ihor,:) = tau_drop(im1,:)
             gaz1(ihor,:) = gaz1(im1,:)
             gaz2(ihor,:) = gaz2(im1,:)
             gaz3(ihor,:) = gaz3(im1,:)
           endif
         endif

         ENDDO             ! Fin de la boucle IHOR
c***************************************************************

102      CONTINUE           ! la premiere fois, c'est une boucle vide!


c***************************************************************
c FIN: on renvoie les nouvelles valeurs q(t+dt)=q(t) + dq(t)
c
c Pour les aerosols, les noyaux, les glaces et les vapeurs modifiees...
c
c***************************************************************

         do n=1,ngrid
           do i=1,nmicro
             do j=1,klev                ! j de 1 a 54 
               tq(n,j,i) = q(n,klev+1-j,i)
             enddo
           enddo
         enddo

 
c      do j=1,15
cc       CH4 -- cnuages
c        write(210,'(I4,3(ES24.16,1X))') j,
c     &  sum(ttq(40,j,2*nrad+1:3*nrad,1)),
c     &  sum(ttq(40,j,2*nrad+1:3*nrad,2)),
c     &  sum(ttq(40,j,2*nrad+1:3*nrad,2))-
c     &  sum(ttq(40,j,2*nrad+1:3*nrad,1))
cc       C2H6 -- cnuages
c        write(211,'(I4,3(ES24.16,1X))') j,
c     &  sum(ttq(40,j,3*nrad+1:4*nrad,1)),
c     &  sum(ttq(40,j,3*nrad+1:4*nrad,2)),
c     &  sum(ttq(40,j,3*nrad+1:4*nrad,2))-
c     &  sum(ttq(40,j,3*nrad+1:4*nrad,1))
cc       C2H2 -- cnuages
c        write(212,'(I4,3(ES24.16,1X))') j,
c     &  sum(ttq(40,j,4*nrad+1:5*nrad,1)),
c     &  sum(ttq(40,j,4*nrad+1:5*nrad,2)),
c     &  sum(ttq(40,j,4*nrad+1:5*nrad,2))-
c     & sum(ttq(40,j,4*nrad+1:5*nrad,1))
c  
cc       CH4 -- snuages
c        write(310,'(I4,3(ES24.16,1X))') j,
c     &  sum(tttq(40,j,2*nrad+1:3*nrad,1)),
c     &  sum(tttq(40,j,2*nrad+1:3*nrad,2)),
c     &  sum(tttq(40,j,2*nrad+1:3*nrad,2))-
c     &  sum(tttq(40,j,2*nrad+1:3*nrad,1))
cc       C2H6 -- snuages
c        write(311,'(I4,3(ES24.16,1X))') j,
c     &  sum(tttq(40,j,3*nrad+1:4*nrad,1)),
c     &  sum(tttq(40,j,3*nrad+1:4*nrad,2)),
c     &  sum(tttq(40,j,3*nrad+1:4*nrad,2))-
c     &  sum(tttq(40,j,3*nrad+1:4*nrad,1))
cc       C2H2 -- snuages
c        write(312,'(I4,3(ES24.16,1X))') j,
c     &  sum(tttq(40,j,4*nrad+1:5*nrad,1)),
c     &  sum(tttq(40,j,4*nrad+1:5*nrad,2)),
c     &  sum(tttq(40,j,4*nrad+1:5*nrad,2))-
c     &  sum(tttq(40,j,4*nrad+1:5*nrad,1))
c      enddo
c      write(210,*) "NEWLINE"
c      write(211,*) "NEWLINE"
c      write(212,*) "NEWLINE"
c      write(310,*) "NEWLINE"
c      write(311,*) "NEWLINE"
c      write(312,*) "NEWLINE"


c       do j=1,20
c         write(410,'(I4,3(ES24.16,1X))') j,
c     &   flxesp_i(40,j,1),flxesp_i(40,j,2),flxesp_i(40,j,3)  
c         write(510,'(I4,3(ES24.16,1X))') j,
c     &   solesp(40,j,1),solesp(40,j,2),solesp(40,j,3)  
c       enddo
c       write(410,*) "NEWLINE"
c       write(510,*) "NEWLINE"

         IPREM=1       ! LA PROCHAINE FOIS NE SERA PLUS LA 1ERE 


 16      return  

         end 


c---------------------------------------------------------------------
