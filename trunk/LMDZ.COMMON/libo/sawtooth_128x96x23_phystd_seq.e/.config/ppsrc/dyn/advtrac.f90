










! $Id: advtrac.F90 1549 2011-07-05 08:41:12Z lguez $

SUBROUTINE advtrac(pbaru,pbarv , p,  masse,q,iapptrac,teta, flxw, pk)
  !     Auteur :  F. Hourdin
  !
  !     Modif. P. Le Van     (20/12/97)
  !            F. Codron     (10/99)
  !            D. Le Croller (07/2001)
  !            M.A Filiberti (04/2002)
  !
  USE infotrac, ONLY: nqtot, iadv, nqperes, ok_iso_verif
  USE control_mod, ONLY: iapp_tracvl, day_step
  USE comconst_mod, ONLY: dtvr


  IMPLICIT NONE
  !
  include "dimensions.h"
  include "paramet.h"
  include "comdissip.h"
  include "comgeom2.h"
  include "iniprint.h"

  !-------------------------------------------------------------------
  !     Arguments
  !-------------------------------------------------------------------
  INTEGER,INTENT(OUT) :: iapptrac
  REAL,INTENT(IN) :: pbaru(ip1jmp1,llm)
  REAL,INTENT(IN) :: pbarv(ip1jm,llm)
  REAL,INTENT(INOUT) :: q(ip1jmp1,llm,nqtot)
  REAL,INTENT(IN) :: masse(ip1jmp1,llm)
  REAL,INTENT(IN) :: p( ip1jmp1,llmp1 )
  REAL,INTENT(IN) :: teta(ip1jmp1,llm)
  REAL,INTENT(IN) :: pk(ip1jmp1,llm)
  REAL,INTENT(OUT) :: flxw(ip1jmp1,llm)
  !-------------------------------------------------------------------
  !     Ajout PPM
  !--------------------------------------------------------
  REAL massebx(ip1jmp1,llm),masseby(ip1jm,llm)
  !-------------------------------------------------------------
  !     Variables locales
  !-------------------------------------------------------------

  REAL pbaruc(ip1jmp1,llm),pbarvc(ip1jm,llm)
  REAL massem(ip1jmp1,llm),zdp(ip1jmp1)
  REAL pbarug(ip1jmp1,llm),pbarvg(ip1jm,llm),wg(ip1jmp1,llm) 
  REAL (kind=kind(1.d0)) :: t_initial, t_final, tps_cpu
  INTEGER iadvtr
  INTEGER ij,l,iq,iiq
  REAL zdpmin, zdpmax
  EXTERNAL  minmax
  SAVE iadvtr, massem, pbaruc, pbarvc
  DATA iadvtr/0/
  !----------------------------------------------------------
  !     Rajouts pour PPM
  !----------------------------------------------------------
  INTEGER indice,n
  REAL dtbon ! Pas de temps adaptatif pour que CFL<1
  REAL CFLmaxz,aaa,bbb ! CFL maximum
  REAL psppm(iim,jjp1) ! pression  au sol
  REAL unatppm(iim,jjp1,llm),vnatppm(iim,jjp1,llm)
  REAL qppm(iim*jjp1,llm,nqtot)
  REAL fluxwppm(iim,jjp1,llm)
  REAL apppm(llmp1), bpppm(llmp1)
  LOGICAL dum,fill
  DATA fill/.true./
  DATA dum/.true./

  integer,save :: countcfl=0
  real cflx(ip1jmp1,llm)
  real cfly(ip1jm,llm)
  real cflz(ip1jmp1,llm)
  real, save :: cflxmax(llm),cflymax(llm),cflzmax(llm)

  IF(iadvtr.EQ.0) THEN
     pbaruc(:,:)=0
     pbarvc(:,:)=0
  ENDIF

  !   accumulation des flux de masse horizontaux
  DO l=1,llm
     DO ij = 1,ip1jmp1
        pbaruc(ij,l) = pbaruc(ij,l) + pbaru(ij,l)
     ENDDO
     DO ij = 1,ip1jm
        pbarvc(ij,l) = pbarvc(ij,l) + pbarv(ij,l)
     ENDDO
  ENDDO

  !   selection de la masse instantannee des mailles avant le transport.
  IF(iadvtr.EQ.0) THEN

     CALL SCOPY(ip1jmp1*llm,masse,1,massem,1)
     !cc         CALL filtreg ( massem ,jjp1, llm,-2, 2, .TRUE., 1 )
     !
  ENDIF

  iadvtr   = iadvtr+1
  iapptrac = iadvtr


  !   Test pour savoir si on advecte a ce pas de temps
  IF ( iadvtr.EQ.iapp_tracvl ) THEN

     !c   ..  Modif P.Le Van  ( 20/12/97 )  ....
     !c

     !   traitement des flux de masse avant advection.
     !     1. calcul de w
     !     2. groupement des mailles pres du pole.

     CALL groupe( massem, pbaruc,pbarvc, pbarug,pbarvg,wg )

     ! ... Flux de masse diaganostiques traceurs
     flxw = wg / REAL(iapp_tracvl)

     !  test sur l'eventuelle creation de valeurs negatives de la masse
     DO l=1,llm-1
        DO ij = iip2+1,ip1jm
           zdp(ij) =    pbarug(ij-1,l)   - pbarug(ij,l) &
                - pbarvg(ij-iip1,l) + pbarvg(ij,l) &
                +       wg(ij,l+1)  - wg(ij,l)
        ENDDO
        CALL SCOPY( jjm -1 ,zdp(iip1+iip1),iip1,zdp(iip2),iip1 )
        DO ij = iip2,ip1jm
           zdp(ij)= zdp(ij)*dtvr/ massem(ij,l) 
        ENDDO


        CALL minmax ( ip1jm-iip1, zdp(iip2), zdpmin,zdpmax )

        IF(MAX(ABS(zdpmin),ABS(zdpmax)).GT.0.5) THEN
           PRINT*,'WARNING DP/P l=',l,'  MIN:',zdpmin, &
                '   MAX:', zdpmax
        ENDIF

     ENDDO


     !-------------------------------------------------------------------
     ! Calcul des criteres CFL en X, Y et Z
     !-------------------------------------------------------------------

     if (countcfl == 0. ) then
        cflxmax(:)=0.
        cflymax(:)=0.
        cflzmax(:)=0.
     endif

     countcfl=countcfl+iapp_tracvl
     cflx(:,:)=0.
     cfly(:,:)=0.
     cflz(:,:)=0.
     do l=1,llm
        do ij=iip2,ip1jm-1
           if (pbarug(ij,l)>=0.) then
              cflx(ij,l)=pbarug(ij,l)*dtvr/masse(ij,l)
           else
              cflx(ij,l)=-pbarug(ij,l)*dtvr/masse(ij+1,l)
           endif
        enddo
     enddo
     do l=1,llm
        do ij=iip2,ip1jm-1,iip1
           cflx(ij+iip1,l)=cflx(ij,l)
        enddo
     enddo

     do l=1,llm
        do ij=1,ip1jm
           if (pbarvg(ij,l)>=0.) then
              cfly(ij,l)=pbarvg(ij,l)*dtvr/masse(ij,l)
           else
              cfly(ij,l)=-pbarvg(ij,l)*dtvr/masse(ij+iip1,l)
           endif
        enddo
     enddo

     do l=2,llm
        do ij=1,ip1jm
           if (wg(ij,l)>=0.) then
              cflz(ij,l)=wg(ij,l)*dtvr/masse(ij,l)
           else
              cflz(ij,l)=-wg(ij,l)*dtvr/masse(ij,l-1)
           endif
        enddo
     enddo

     do l=1,llm
        cflxmax(l)=max(cflxmax(l),maxval(cflx(:,l)))
        cflymax(l)=max(cflymax(l),maxval(cfly(:,l)))
        cflzmax(l)=max(cflzmax(l),maxval(cflz(:,l)))
     enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Par defaut, on sort le diagnostic des CFL tous les jours.
     ! Si on veut le sortir a chaque pas d'advection en cas de plantage 
     !     if (countcfl==iapp_tracvl) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if (countcfl==day_step) then
        do l=1,llm
           write(lunout,*) 'L, CFL[xyz]max:', l, cflxmax(l), cflymax(l), &
                cflzmax(l)
        enddo
        countcfl=0
     endif

     !-------------------------------------------------------------------
     !   Advection proprement dite (Modification Le Croller (07/2001)
     !-------------------------------------------------------------------

     !----------------------------------------------------
     !        Calcul des moyennes basées sur la masse
     !----------------------------------------------------
     call massbar(massem,massebx,masseby)          

     !-----------------------------------------------------------
     !     Appel des sous programmes d'advection
     !-----------------------------------------------------------
     if (ok_iso_verif) then
           write(*,*) 'advtrac 227'
           call check_isotopes_seq(q,ip1jmp1,'advtrac 162')
     endif !if (ok_iso_verif) then

     do iq=1,nqperes
        !        call clock(t_initial)
        if(iadv(iq) == 0) cycle 
        !   ----------------------------------------------------------------
        !   Schema de Van Leer I MUSCL
        !   ----------------------------------------------------------------
        if(iadv(iq).eq.10) THEN
           ! CRisi: on fait passer tout q pour avoir acces aux fils
           
           !write(*,*) 'advtrac 239: iq,q(1721,19,:)=',iq,q(1721,19,:)     
           call vlsplt(q,2.,massem,wg,pbarug,pbarvg,dtvr,iq)

           !   ----------------------------------------------------------------
           !   Schema "pseudo amont" + test sur humidite specifique
           !    pour la vapeur d'eau. F. Codron
           !   ----------------------------------------------------------------
        else if(iadv(iq).eq.14) then
           !
           !write(*,*) 'advtrac 248: iq,q(1721,19,:)=',iq,q(1721,19,:)
           CALL vlspltqs( q, 2., massem, wg , &
                pbarug,pbarvg,dtvr,p,pk,teta,iq)

           !   ----------------------------------------------------------------
           !   Schema de Frederic Hourdin
           !   ----------------------------------------------------------------
        else if(iadv(iq).eq.12) then
           !            Pas de temps adaptatif
           call adaptdt(iadv(iq),dtbon,n,pbarug,massem)
           if (n.GT.1) then
              write(*,*) 'WARNING horizontal dt=',dtbon,'dtvr=', &
                   dtvr,'n=',n
           endif
           do indice=1,n
              call advn(q(1,1,iq),massem,wg,pbarug,pbarvg,dtbon,1)
           end do
        else if(iadv(iq).eq.13) then
           !            Pas de temps adaptatif
           call adaptdt(iadv(iq),dtbon,n,pbarug,massem)
           if (n.GT.1) then
              write(*,*) 'WARNING horizontal dt=',dtbon,'dtvr=', &
                   dtvr,'n=',n
           endif
           do indice=1,n
              call advn(q(1,1,iq),massem,wg,pbarug,pbarvg,dtbon,2)
           end do
           !   ----------------------------------------------------------------
           !   Schema de pente SLOPES
           !   ----------------------------------------------------------------
        else if (iadv(iq).eq.20) then
           call pentes_ini (q(1,1,iq),wg,massem,pbarug,pbarvg,0)

           !   ----------------------------------------------------------------
           !   Schema de Prather
           !   ----------------------------------------------------------------
        else if (iadv(iq).eq.30) then
           !            Pas de temps adaptatif
           call adaptdt(iadv(iq),dtbon,n,pbarug,massem)
           if (n.GT.1) then
              write(*,*) 'WARNING horizontal dt=',dtbon,'dtvr=', &
                   dtvr,'n=',n
           endif
           call  prather(q(1,1,iq),wg,massem,pbarug,pbarvg, &
                n,dtbon)

           !   ----------------------------------------------------------------
           !   Schemas PPM Lin et Rood
           !   ----------------------------------------------------------------
        else if (iadv(iq).eq.11.OR.(iadv(iq).GE.16.AND. &
             iadv(iq).LE.18)) then

           !        Test sur le flux horizontal
           !        Pas de temps adaptatif
           call adaptdt(iadv(iq),dtbon,n,pbarug,massem)
           if (n.GT.1) then
              write(*,*) 'WARNING horizontal dt=',dtbon,'dtvr=', &
                   dtvr,'n=',n
           endif
           !        Test sur le flux vertical
           CFLmaxz=0.
           do l=2,llm
              do ij=iip2,ip1jm
                 aaa=wg(ij,l)*dtvr/massem(ij,l)
                 CFLmaxz=max(CFLmaxz,aaa)
                 bbb=-wg(ij,l)*dtvr/massem(ij,l-1)
                 CFLmaxz=max(CFLmaxz,bbb)
              enddo
           enddo
           if (CFLmaxz.GE.1) then
              write(*,*) 'WARNING vertical','CFLmaxz=', CFLmaxz
           endif

           !-----------------------------------------------------------
           !        Ss-prg interface LMDZ.4->PPM3d
           !-----------------------------------------------------------

           call interpre(q(1,1,iq),qppm(1,1,iq),wg,fluxwppm,massem, &
                apppm,bpppm,massebx,masseby,pbarug,pbarvg, &
                unatppm,vnatppm,psppm)

           do indice=1,n
              !----------------------------------------------------------------
              !                         VL (version PPM) horiz. et PPM vert.
              !----------------------------------------------------------------
              if (iadv(iq).eq.11) then
                 !                  Ss-prg PPM3d de Lin
                 call ppm3d(1,qppm(1,1,iq), &
                      psppm,psppm, &
                      unatppm,vnatppm,fluxwppm,dtbon,2,2,2,1, &
                      iim,jjp1,2,llm,apppm,bpppm,0.01,6400000, &
                      fill,dum,220.)

                 !-------------------------------------------------------------
                 !                           Monotonic PPM
                 !-------------------------------------------------------------
              else if (iadv(iq).eq.16) then
                 !                  Ss-prg PPM3d de Lin
                 call ppm3d(1,qppm(1,1,iq), &
                      psppm,psppm, &
                      unatppm,vnatppm,fluxwppm,dtbon,3,3,3,1, &
                      iim,jjp1,2,llm,apppm,bpppm,0.01,6400000, &
                      fill,dum,220.)
                 !-------------------------------------------------------------

                 !-------------------------------------------------------------
                 !                           Semi Monotonic PPM
                 !-------------------------------------------------------------
              else if (iadv(iq).eq.17) then
                 !                  Ss-prg PPM3d de Lin
                 call ppm3d(1,qppm(1,1,iq), &
                      psppm,psppm, &
                      unatppm,vnatppm,fluxwppm,dtbon,4,4,4,1, &
                      iim,jjp1,2,llm,apppm,bpppm,0.01,6400000, &
                      fill,dum,220.)
                 !-------------------------------------------------------------

                 !-------------------------------------------------------------
                 !                         Positive Definite PPM
                 !-------------------------------------------------------------
              else if (iadv(iq).eq.18) then
                 !                  Ss-prg PPM3d de Lin
                 call ppm3d(1,qppm(1,1,iq), &
                      psppm,psppm, &
                      unatppm,vnatppm,fluxwppm,dtbon,5,5,5,1, &
                      iim,jjp1,2,llm,apppm,bpppm,0.01,6400000, &
                      fill,dum,220.)
                 !-------------------------------------------------------------
              endif
           enddo
           !-----------------------------------------------------------------
           !               Ss-prg interface PPM3d-LMDZ.4
           !-----------------------------------------------------------------
           call interpost(q(1,1,iq),qppm(1,1,iq))
        endif
        !----------------------------------------------------------------------

        !-----------------------------------------------------------------
        ! On impose une seule valeur du traceur au pôle Sud j=jjm+1=jjp1
        ! et Nord j=1
        !-----------------------------------------------------------------

        !                  call traceurpole(q(1,1,iq),massem)

        ! calcul du temps cpu pour un schema donne

        !                  call clock(t_final)
        !ym                  tps_cpu=t_final-t_initial
        !ym                  cpuadv(iq)=cpuadv(iq)+tps_cpu

     end DO

     if (ok_iso_verif) then
           write(*,*) 'advtrac 402'
           call check_isotopes_seq(q,ip1jmp1,'advtrac 397')
     endif !if (ok_iso_verif) then

     !------------------------------------------------------------------
     !   on reinitialise a zero les flux de masse cumules
     !---------------------------------------------------
     iadvtr=0

  ENDIF ! if iadvtr.EQ.iapp_tracvl

END SUBROUTINE advtrac
