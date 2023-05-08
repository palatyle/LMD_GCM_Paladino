!
! $Id: bilan_dyn_p.F 1907 2013-11-26 13:10:46Z lguez $
!
      SUBROUTINE bilan_dyn_p (ntrac,dt_app,dt_cum,
     s  ps,masse,pk,flux_u,flux_v,teta,phi,ucov,vcov,trac)

c   AFAIRE
c   Prevoir en champ nq+1 le diagnostique de l'energie
c   en faisant Qzon=Cv T + L * ...
c             vQ..A=Cp T + L * ...

#ifdef CPP_IOIPSL
      USE IOIPSL
#endif
      USE parallel_lmdz
      USE mod_hallo
      use misc_mod
      use write_field_p
      USE comvert_mod, ONLY: presnivs
      USE comconst_mod, ONLY: cpp,pi
      USE temps_mod, ONLY: annee_ref,day_ref,itau_dyn
      IMPLICIT NONE

#include "dimensions.h"
#include "paramet.h"
#include "comgeom2.h"
#include "iniprint.h"

c====================================================================
c
c   Sous-programme consacre à des diagnostics dynamiques de base
c
c 
c   De facon generale, les moyennes des scalaires Q sont ponderees par
c   la masse.
c
c   Les flux de masse sont eux simplement moyennes.
c
c====================================================================

c   Arguments :
c   ===========

      integer ntrac
      real dt_app,dt_cum
      real ps(iip1,jjp1)
      real masse(iip1,jjp1,llm),pk(iip1,jjp1,llm)
      real flux_u(iip1,jjp1,llm)
      real flux_v(iip1,jjm,llm)
      real teta(iip1,jjp1,llm)
      real phi(iip1,jjp1,llm)
      real ucov(iip1,jjp1,llm)
      real vcov(iip1,jjm,llm)
      real trac(iip1,jjp1,llm,ntrac)

c   Local :
c   =======

      integer,save :: icum,ncum
!$OMP THREADPRIVATE(icum,ncum)
      logical,SAVE :: first=.true.
!$OMP THREADPRIVATE(first) 

      real zz,zqy
      real,save :: zfactv(jjm,llm)

      integer,parameter :: nQ=7


cym      character*6 nom(nQ)
cym      character*6 unites(nQ)
      character(len=6),save :: nom(nQ)
      character(len=6),save :: unites(nQ)

      character(len=10) file
      integer ifile
      parameter (ifile=4)

      integer,PARAMETER :: itemp=1,igeop=2,iecin=3,iang=4,iu=5
      INTEGER,PARAMETER :: iovap=6,iun=7
      integer,PARAMETER :: i_sortie=1

      real,SAVE :: time=0.
      integer,SAVE :: itau=0.
!$OMP THREADPRIVATE(time,itau)

      real ww

c   variables dynamiques intermédiaires
      REAL,save :: vcont(iip1,jjm,llm),ucont(iip1,jjp1,llm)
      REAL,save :: ang(iip1,jjp1,llm),unat(iip1,jjp1,llm)
      REAL,save :: massebx(iip1,jjp1,llm),masseby(iip1,jjm,llm)
      REAL,save :: vorpot(iip1,jjm,llm)
      REAL,save :: w(iip1,jjp1,llm),ecin(iip1,jjp1,llm)
      REAL,save ::convm(iip1,jjp1,llm)
      REAL,save :: bern(iip1,jjp1,llm)

c   champ contenant les scalaires advectés.
      real,save :: Q(iip1,jjp1,llm,nQ)
    
c   champs cumulés
      real,save :: ps_cum(iip1,jjp1)
      real,save :: masse_cum(iip1,jjp1,llm)
      real,save :: flux_u_cum(iip1,jjp1,llm)
      real,save :: flux_v_cum(iip1,jjm,llm)
      real,save :: Q_cum(iip1,jjp1,llm,nQ)
      real,save :: flux_uQ_cum(iip1,jjp1,llm,nQ)
      real,save :: flux_vQ_cum(iip1,jjm,llm,nQ)
      real,save :: flux_wQ_cum(iip1,jjp1,llm,nQ)
      real,save :: dQ(iip1,jjp1,llm,nQ)


c   champs de tansport en moyenne zonale
      integer ntr,itr
      parameter (ntr=5)

cym      character*10 znom(ntr,nQ)
cym      character*20 znoml(ntr,nQ)
cym      character*10 zunites(ntr,nQ)
      character*10,save :: znom(ntr,nQ)
      character*20,save :: znoml(ntr,nQ)
      character*10,save :: zunites(ntr,nQ)

      INTEGER,PARAMETER :: iave=1,itot=2,immc=3,itrs=4,istn=5

      character*3 ctrs(ntr)
      data ctrs/'  ','TOT','MMC','TRS','STN'/

      real,save :: zvQ(jjm,llm,ntr,nQ),zvQtmp(jjm,llm)
      real,save :: zavQ(jjm,ntr,nQ),psiQ(jjm,llm+1,nQ)
      real,save :: zmasse(jjm,llm),zamasse(jjm)

      real,save :: zv(jjm,llm),psi(jjm,llm+1)

      integer i,j,l,iQ


c   Initialisation du fichier contenant les moyennes zonales.
c   ---------------------------------------------------------

      character*10 infile

      integer fileid
      integer thoriid, zvertiid
      save fileid

      integer,save :: ndex3d(jjm*llm)

C   Variables locales
C
      integer tau0
      real zjulian
      character*3 str
      character*10 ctrac
      integer ii,jj
      integer zan, dayref
C
      real,save :: rlong(jjm),rlatg(jjm)
      integer :: jjb,jje,jjn,ijb,ije
      type(Request),SAVE :: Req
!$OMP THREADPRIVATE(Req)

! definition du domaine d'ecriture pour le rebuild

      INTEGER,DIMENSION(1) :: ddid
      INTEGER,DIMENSION(1) :: dsg
      INTEGER,DIMENSION(1) :: dsl
      INTEGER,DIMENSION(1) :: dpf
      INTEGER,DIMENSION(1) :: dpl
      INTEGER,DIMENSION(1) :: dhs
      INTEGER,DIMENSION(1) :: dhe 
      
      INTEGER :: bilan_dyn_domain_id


c=====================================================================
c   Initialisation
c=====================================================================
      if (adjust) return
      
      time=time+dt_app
      itau=itau+1

      if (first) then

        ndex3d=0

        icum=0
c       initialisation des fichiers
        first=.false.
c   ncum est la frequence de stokage en pas de temps
        ncum=dt_cum/dt_app
        if (abs(ncum*dt_app-dt_cum).gt.1.e-5*dt_app) then
           WRITE(lunout,*)
     .            'Pb : le pas de cumule doit etre multiple du pas'
           WRITE(lunout,*)'dt_app=',dt_app
           WRITE(lunout,*)'dt_cum=',dt_cum
           stop
        else
          write(lunout,*) "bilan_dyn_p: ncum=",ncum
        endif

!        if (i_sortie.eq.1) then
!	 file='dynzon'
!         if (mpi_rank==0) then
!	 call inigrads(ifile,1
!     s  ,0.,180./pi,0.,0.,jjm,rlatv,-90.,90.,180./pi
!     s  ,llm,presnivs,1.
!     s  ,dt_cum,file,'dyn_zon ')
!         endif
!        endif

!$OMP MASTER
        nom(itemp)='T'
        nom(igeop)='gz'
        nom(iecin)='K'
        nom(iang)='ang'
        nom(iu)='u'
        nom(iovap)='ovap'
        nom(iun)='un'

        unites(itemp)='K'
        unites(igeop)='m2/s2'
        unites(iecin)='m2/s2'
        unites(iang)='ang'
        unites(iu)='m/s'
        unites(iovap)='kg/kg'
        unites(iun)='un'


c   Initialisation du fichier contenant les moyennes zonales.
c   ---------------------------------------------------------

      infile='dynzon'

      zan = annee_ref
      dayref = day_ref
      CALL ymds2ju(zan, 1, dayref, 0.0, zjulian)
      tau0 = itau_dyn
      
      rlong=0.
      rlatg=rlatv*180./pi

      jjb=jj_begin
      jje=jj_end
      jjn=jj_nb
      IF (pole_sud) THEN
        jjn=jj_nb-1
        jje=jj_end-1
      ENDIF

      ddid=(/ 2 /)
      dsg=(/ jjm /)
      dsl=(/ jjn /)
      dpf=(/ jjb /)
      dpl=(/ jje /)
      dhs=(/ 0 /)
      dhe=(/ 0 /)

      call flio_dom_set(mpi_size,mpi_rank,ddid,dsg,dsl,dpf,dpl,dhs,dhe, 
     .                 'box',bilan_dyn_domain_id)
       
      call histbeg(trim(infile),
     .             1, rlong(jjb:jje), jjn, rlatg(jjb:jje),
     .             1, 1, 1, jjn,
     .             tau0, zjulian, dt_cum, thoriid, fileid,
     .             bilan_dyn_domain_id)

C
C  Appel a histvert pour la grille verticale
C
      call histvert(fileid, 'presnivs', 'Niveaux sigma','mb',
     .              llm, presnivs, zvertiid)
C
C  Appels a histdef pour la definition des variables a sauvegarder
      do iQ=1,nQ
         do itr=1,ntr
            if(itr.eq.1) then
               znom(itr,iQ)=nom(iQ)
               znoml(itr,iQ)=nom(iQ)
               zunites(itr,iQ)=unites(iQ)
            else
               znom(itr,iQ)=ctrs(itr)//'v'//nom(iQ)
               znoml(itr,iQ)='transport : v * '//nom(iQ)//' '//ctrs(itr)
               zunites(itr,iQ)='m/s * '//unites(iQ)
            endif
         enddo
      enddo

c   Declarations des champs avec dimension verticale
c      print*,'1HISTDEF'
      do iQ=1,nQ
         do itr=1,ntr
      IF (prt_level > 5)
     . WRITE(lunout,*)'var ',itr,iQ
     .      ,znom(itr,iQ),znoml(itr,iQ),zunites(itr,iQ)
            call histdef(fileid,znom(itr,iQ),znoml(itr,iQ),
     .        zunites(itr,iQ),1,jjn,thoriid,llm,1,llm,zvertiid,
     .        32,'ave(X)',dt_cum,dt_cum)
         enddo
c   Declarations pour les fonctions de courant
c      print*,'2HISTDEF'
          call histdef(fileid,'psi'//nom(iQ)
     .      ,'stream fn. '//znoml(itot,iQ),
     .      zunites(itot,iQ),1,jjn,thoriid,llm,1,llm,zvertiid,
     .      32,'ave(X)',dt_cum,dt_cum)
      enddo


c   Declarations pour les champs de transport d'air
c      print*,'3HISTDEF'
      call histdef(fileid, 'masse', 'masse',
     .             'kg', 1, jjn, thoriid, llm, 1, llm, zvertiid,
     .             32, 'ave(X)', dt_cum, dt_cum)
      call histdef(fileid, 'v', 'v',
     .             'm/s', 1, jjn, thoriid, llm, 1, llm, zvertiid,
     .             32, 'ave(X)', dt_cum, dt_cum)
c   Declarations pour les fonctions de courant
c      print*,'4HISTDEF'
          call histdef(fileid,'psi','stream fn. MMC ','mega t/s',
     .      1,jjn,thoriid,llm,1,llm,zvertiid,
     .      32,'ave(X)',dt_cum,dt_cum)


c   Declaration des champs 1D de transport en latitude
c      print*,'5HISTDEF'
      do iQ=1,nQ
         do itr=2,ntr
            call histdef(fileid,'a'//znom(itr,iQ),znoml(itr,iQ),
     .        zunites(itr,iQ),1,jjn,thoriid,1,1,1,-99,
     .        32,'ave(X)',dt_cum,dt_cum)
         enddo
      enddo


c      print*,'8HISTDEF'
               CALL histend(fileid)

!$OMP END MASTER
!$OMP BARRIER
      endif


c=====================================================================
c   Calcul des champs dynamiques
c   ----------------------------

      jjb=jj_begin
      jje=jj_end
    
c   énergie cinétique
!      ucont(:,jjb:jje,:)=0

      call Register_Hallo(ucov,ip1jmp1,llm,1,1,1,1,Req)
      call Register_Hallo(vcov,ip1jm,llm,1,1,1,1,Req)
      call SendRequest(Req)
c$OMP BARRIER
      call WaitRequest(Req)
c$OMP BARRIER

      CALL covcont_p(llm,ucov,vcov,ucont,vcont)
      CALL enercin_p(vcov,ucov,vcont,ucont,ecin)

c   moment cinétique
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      do l=1,llm
         ang(:,jjb:jje,l)=ucov(:,jjb:jje,l)+constang(:,jjb:jje)
         unat(:,jjb:jje,l)=ucont(:,jjb:jje,l)*cu(:,jjb:jje)
      enddo
!$OMP END DO

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
        Q(:,jjb:jje,l,itemp)=teta(:,jjb:jje,l)*pk(:,jjb:jje,l)/cpp
        Q(:,jjb:jje,l,igeop)=phi(:,jjb:jje,l)
        Q(:,jjb:jje,l,iecin)=ecin(:,jjb:jje,l)
        Q(:,jjb:jje,l,iang)=ang(:,jjb:jje,l)
        Q(:,jjb:jje,l,iu)=unat(:,jjb:jje,l)
        Q(:,jjb:jje,l,iovap)=trac(:,jjb:jje,l,1)
        Q(:,jjb:jje,l,iun)=1.
      ENDDO
!$OMP END DO NOWAIT

c=====================================================================
c   Cumul
c=====================================================================
c
      if(icum.EQ.0) then
         jjb=jj_begin
         jje=jj_end

!$OMP MASTER
         ps_cum(:,jjb:jje)=0.
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO l=1,llm
          masse_cum(:,jjb:jje,l)=0.
          flux_u_cum(:,jjb:jje,l)=0.
          Q_cum(:,jjb:jje,l,:)=0.
          flux_uQ_cum(:,jjb:jje,l,:)=0.
          if (pole_sud) jje=jj_end-1
          flux_v_cum(:,jjb:jje,l)=0.
          flux_vQ_cum(:,jjb:jje,l,:)=0.
        ENDDO
!$OMP END DO NOWAIT
      endif

      IF (prt_level > 5)
     . WRITE(lunout,*)'dans bilan_dyn ',icum,'->',icum+1
      icum=icum+1

c   accumulation des flux de masse horizontaux
      jjb=jj_begin
      jje=jj_end

!$OMP MASTER
      ps_cum(:,jjb:jje)=ps_cum(:,jjb:jje)+ps(:,jjb:jje)
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
        masse_cum(:,jjb:jje,l)=masse_cum(:,jjb:jje,l)+masse(:,jjb:jje,l)
        flux_u_cum(:,jjb:jje,l)=flux_u_cum(:,jjb:jje,l)
     .                         +flux_u(:,jjb:jje,l)
      ENDDO
!$OMP END DO NOWAIT
      
      if (pole_sud) jje=jj_end-1

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
       flux_v_cum(:,jjb:jje,l)=flux_v_cum(:,jjb:jje,l)
     .                          +flux_v(:,jjb:jje,l)
      ENDDO
!$OMP END DO NOWAIT
      
      jjb=jj_begin
      jje=jj_end

      do iQ=1,nQ
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO l=1,llm
          Q_cum(:,jjb:jje,l,iQ)=Q_cum(:,jjb:jje,l,iQ)
     .                       +Q(:,jjb:jje,l,iQ)*masse(:,jjb:jje,l)
        ENDDO
!$OMP END DO NOWAIT
      enddo

c=====================================================================
c  FLUX ET TENDANCES
c=====================================================================

c   Flux longitudinal
c   -----------------
      do iQ=1,nQ
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
         do l=1,llm
            do j=jjb,jje
               do i=1,iim
                  flux_uQ_cum(i,j,l,iQ)=flux_uQ_cum(i,j,l,iQ)
     s            +flux_u(i,j,l)*0.5*(Q(i,j,l,iQ)+Q(i+1,j,l,iQ))
               enddo
               flux_uQ_cum(iip1,j,l,iQ)=flux_uQ_cum(1,j,l,iQ)
            enddo
         enddo
!$OMP END DO NOWAIT
      enddo

c    flux méridien
c    -------------
      do iQ=1,nQ
        call Register_Hallo(Q(1,1,1,iQ),ip1jmp1,llm,0,1,1,0,Req)
      enddo
      call SendRequest(Req)
!$OMP BARRIER      
      call WaitRequest(Req)
!$OMP BARRIER

      jjb=jj_begin
      jje=jj_end
      if (pole_sud) jje=jj_end-1
      
      do iQ=1,nQ
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
         do l=1,llm
            do j=jjb,jje
               do i=1,iip1
                  flux_vQ_cum(i,j,l,iQ)=flux_vQ_cum(i,j,l,iQ)
     s            +flux_v(i,j,l)*0.5*(Q(i,j,l,iQ)+Q(i,j+1,l,iQ))
               enddo
            enddo
         enddo
!$OMP END DO NOWAIT
      enddo


c    tendances
c    ---------

c   convergence horizontale
      call Register_Hallo(flux_uQ_cum,ip1jmp1,llm,2,2,2,2,Req)
      call Register_Hallo(flux_vQ_cum,ip1jm,llm,2,2,2,2,Req)
      call SendRequest(Req)
!$OMP BARRIER      
      call WaitRequest(Req)
c$OMP BARRIER

      call  convflu_p(flux_uQ_cum,flux_vQ_cum,llm*nQ,dQ)

c   calcul de la vitesse verticale
      call Register_Hallo(flux_u_cum,ip1jmp1,llm,2,2,2,2,Req)
      call Register_Hallo(flux_v_cum,ip1jm,llm,2,2,2,2,Req)
      call SendRequest(Req)
!$OMP BARRIER      
      call WaitRequest(Req)
c$OMP BARRIER

      call convmas_p(flux_u_cum,flux_v_cum,convm)
      CALL vitvert_p(convm,w)
!$OMP BARRIER      

      jjb=jj_begin
      jje=jj_end

      do iQ=1,nQ
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
         do l=1,llm
            IF (l<llm) THEN
              do j=jjb,jje
                 do i=1,iip1
                    ww=-0.5*w(i,j,l+1)*(Q(i,j,l,iQ)+Q(i,j,l+1,iQ))
                    dQ(i,j,l  ,iQ)=dQ(i,j,l  ,iQ)-ww
                    dQ(i,j,l+1,iQ)=dQ(i,j,l+1,iQ)+ww
                 enddo
              enddo
            ENDIF
            IF (l>2) THEN
              do j=jjb,jje
                do i=1,iip1
                  ww=-0.5*w(i,j,l)*(Q(i,j,l-1,iQ)+Q(i,j,l,iQ))
                  dQ(i,j,l,iQ)=dQ(i,j,l,iQ)+ww
                enddo
              enddo
            ENDIF
         enddo
!$OMP ENDDO NOWAIT 
      enddo
      IF (prt_level > 5)
     . WRITE(lunout,*)'Apres les calculs fait a chaque pas'
c=====================================================================
c   PAS DE TEMPS D'ECRITURE
c=====================================================================
      if (icum.eq.ncum) then
c=====================================================================

      IF (prt_level > 5)
     . WRITE(lunout,*)'Pas d ecriture'

      jjb=jj_begin
      jje=jj_end

c   Normalisation
      do iQ=1,nQ
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        do l=1,llm
          Q_cum(:,jjb:jje,l,iQ)=Q_cum(:,jjb:jje,l,iQ) 
     .	                        /masse_cum(:,jjb:jje,l)
        enddo
!$OMP ENDDO NOWAIT 
      enddo

      zz=1./REAL(ncum)

!$OMP MASTER
      ps_cum(:,jjb:jje)=ps_cum(:,jjb:jje)*zz
!$OMP END MASTER

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
        masse_cum(:,jjb:jje,l)=masse_cum(:,jjb:jje,l)*zz
        flux_u_cum(:,jjb:jje,l)=flux_u_cum(:,jjb:jje,l)*zz
        flux_uQ_cum(:,jjb:jje,l,:)=flux_uQ_cum(:,jjb:jje,l,:)*zz
        dQ(:,jjb:jje,l,:)=dQ(:,jjb:jje,l,:)*zz
      ENDDO
!$OMP ENDDO NOWAIT 
         
      
      IF (pole_sud) jje=jj_end-1
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l=1,llm
        flux_v_cum(:,jjb:jje,l)=flux_v_cum(:,jjb:jje,l)*zz
        flux_vQ_cum(:,jjb:jje,l,:)=flux_vQ_cum(:,jjb:jje,l,:)*zz
      ENDDO
!$OMP ENDDO

      jjb=jj_begin
      jje=jj_end


c   A retravailler eventuellement
c   division de dQ par la masse pour revenir aux bonnes grandeurs
      do iQ=1,nQ
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO l=1,llm
           dQ(:,jjb:jje,l,iQ)=dQ(:,jjb:jje,l,iQ)/masse_cum(:,jjb:jje,l)
        ENDDO
!$OMP ENDDO NOWAIT 
      enddo
 
c=====================================================================
c   Transport méridien
c=====================================================================

c   cumul zonal des masses des mailles
c   ----------------------------------
      jjb=jj_begin
      jje=jj_end
      if (pole_sud) jje=jj_end-1

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
        DO l=1,llm
          zv(jjb:jje,l)=0.
          zmasse(jjb:jje,l)=0.
        ENDDO
!$OMP ENDDO NOWAIT 

      call Register_Hallo(masse_cum,ip1jmp1,llm,1,1,1,1,Req)
      do iQ=1,nQ
        call Register_Hallo(Q_cum(1,1,1,iQ),ip1jmp1,llm,0,1,1,0,Req)
      enddo

      call SendRequest(Req)
!$OMP BARRIER
      call WaitRequest(Req)
c$OMP BARRIER

      call massbar_p(masse_cum,massebx,masseby)
      
      jjb=jj_begin
      jje=jj_end
      if (pole_sud) jje=jj_end-1
      
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      do l=1,llm
         do j=jjb,jje
            do i=1,iim
               zmasse(j,l)=zmasse(j,l)+masseby(i,j,l)
               zv(j,l)=zv(j,l)+flux_v_cum(i,j,l)
            enddo
            zfactv(j,l)=cv(1,j)/zmasse(j,l)
         enddo
      enddo
!$OMP ENDDO 

c     print*,'3OK'
c   --------------------------------------------------------------
c   calcul de la moyenne zonale du transport :
c   ------------------------------------------
c
c                                     --
c TOT : la circulation totale       [ vq ]
c
c                                      -     -
c MMC : mean meridional circulation [ v ] [ q ]
c
c                                     ----      --       - -
c TRS : transitoires                [ v'q'] = [ vq ] - [ v q ]
c
c                                     - * - *       - -       -     -
c STT : stationaires                [ v   q   ] = [ v q ] - [ v ] [ q ]
c
c                                              - -
c    on utilise aussi l'intermediaire TMP :  [ v q ]
c
c    la variable zfactv transforme un transport meridien cumule
c    en kg/s * unte-du-champ-transporte en m/s * unite-du-champ-transporte
c
c   --------------------------------------------------------------


c   ----------------------------------------
c   Transport dans le plan latitude-altitude
c   ----------------------------------------

      jjb=jj_begin
      jje=jj_end
      if (pole_sud) jje=jj_end-1
      
      zvQ=0.
      psiQ=0.
      do iQ=1,nQ
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
         do l=1,llm
            zvQtmp(:,l)=0.
            do j=jjb,jje
c              print*,'j,l,iQ=',j,l,iQ
c   Calcul des moyennes zonales du transort total et de zvQtmp
               do i=1,iim
                  zvQ(j,l,itot,iQ)=zvQ(j,l,itot,iQ)
     s                            +flux_vQ_cum(i,j,l,iQ)
                  zqy=      0.5*(Q_cum(i,j,l,iQ)*masse_cum(i,j,l)+
     s                           Q_cum(i,j+1,l,iQ)*masse_cum(i,j+1,l))
                  zvQtmp(j,l)=zvQtmp(j,l)+flux_v_cum(i,j,l)*zqy
     s             /(0.5*(masse_cum(i,j,l)+masse_cum(i,j+1,l)))
                  zvQ(j,l,iave,iQ)=zvQ(j,l,iave,iQ)+zqy
               enddo
c              print*,'aOK'
c   Decomposition
               zvQ(j,l,iave,iQ)=zvQ(j,l,iave,iQ)/zmasse(j,l)
               zvQ(j,l,itot,iQ)=zvQ(j,l,itot,iQ)*zfactv(j,l)
               zvQtmp(j,l)=zvQtmp(j,l)*zfactv(j,l)
               zvQ(j,l,immc,iQ)=zv(j,l)*zvQ(j,l,iave,iQ)*zfactv(j,l)
               zvQ(j,l,itrs,iQ)=zvQ(j,l,itot,iQ)-zvQtmp(j,l)
               zvQ(j,l,istn,iQ)=zvQtmp(j,l)-zvQ(j,l,immc,iQ)
            enddo
         enddo
!$OMP ENDDO NOWAIT 
c   fonction de courant meridienne pour la quantite Q
!$OMP BARRIER
!$OMP MASTER
         do l=llm,1,-1
            do j=jjb,jje
               psiQ(j,l,iQ)=psiQ(j,l+1,iQ)+zvQ(j,l,itot,iQ)
            enddo
         enddo
!$OMP END MASTER
!$OMP BARRIER
      enddo ! of do iQ=1,nQ

c   fonction de courant pour la circulation meridienne moyenne
!$OMP BARRIER
!$OMP MASTER
      psi(jjb:jje,:)=0.
      do l=llm,1,-1
         do j=jjb,jje
            psi(j,l)=psi(j,l+1)+zv(j,l)
            zv(j,l)=zv(j,l)*zfactv(j,l)
         enddo
      enddo
!$OMP END MASTER
!$OMP BARRIER

c     print*,'4OK'
c   sorties proprement dites
!$OMP MASTER      
      if (i_sortie.eq.1) then
      jjb=jj_begin
      jje=jj_end
      jjn=jj_nb
      if (pole_sud) jje=jj_end-1
      if (pole_sud) jjn=jj_nb-1
      
      do iQ=1,nQ
         do itr=1,ntr
            call histwrite(fileid,znom(itr,iQ),itau,
     s                     zvQ(jjb:jje,:,itr,iQ)
     s                     ,jjn*llm,ndex3d)
         enddo
         call histwrite(fileid,'psi'//nom(iQ),
     s                  itau,psiQ(jjb:jje,1:llm,iQ)
     s                  ,jjn*llm,ndex3d)
      enddo
      call histwrite(fileid,'masse',itau,zmasse(jjb:jje,1:llm)
     s   ,jjn*llm,ndex3d)
      call histwrite(fileid,'v',itau,zv(jjb:jje,1:llm)
     s   ,jjn*llm,ndex3d)
      psi(jjb:jje,:)=psi(jjb:jje,:)*1.e-9
      call histwrite(fileid,'psi',itau,psi(jjb:jje,1:llm),
     s               jjn*llm,ndex3d)

      endif


c   -----------------
c   Moyenne verticale
c   -----------------

      zamasse(jjb:jje)=0.
      do l=1,llm
         zamasse(jjb:jje)=zamasse(jjb:jje)+zmasse(jjb:jje,l)
      enddo
     
      zavQ(jjb:jje,:,:)=0.
      do iQ=1,nQ
         do itr=2,ntr
            do l=1,llm
               zavQ(jjb:jje,itr,iQ)=zavQ(jjb:jje,itr,iQ)
     s                             +zvQ(jjb:jje,l,itr,iQ)
     s                             *zmasse(jjb:jje,l)
            enddo
            zavQ(jjb:jje,itr,iQ)=zavQ(jjb:jje,itr,iQ)/zamasse(jjb:jje)
            call histwrite(fileid,'a'//znom(itr,iQ),itau,
     s                     zavQ(jjb:jje,itr,iQ),jjn*llm,ndex3d)
         enddo
      enddo
!$OMP END MASTER
!$OMP BARRIER
c     on doit pouvoir tracer systematiquement la fonction de courant.

c=====================================================================
c/////////////////////////////////////////////////////////////////////
      icum=0                  !///////////////////////////////////////
      endif ! icum.eq.ncum    !///////////////////////////////////////
c/////////////////////////////////////////////////////////////////////
c=====================================================================
      return
      end

