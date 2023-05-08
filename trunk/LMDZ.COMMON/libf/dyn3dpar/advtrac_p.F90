! $Id: advtrac_p.F90 1549 2011-07-05 08:41:12Z lguez $

SUBROUTINE advtrac_p(pbaru,pbarv , p,  masse,q,iapptrac,teta, flxw, pk)

  !     Auteur :  F. Hourdin
  !
  !     Modif. P. Le Van     (20/12/97)
  !            F. Codron     (10/99)
  !            D. Le Croller (07/2001)
  !            M.A Filiberti (04/2002)
  !
  USE parallel_lmdz, ONLY: ij_begin,ij_end,OMP_CHUNK,pole_nord,pole_sud,&
                           setdistrib
  USE Write_Field_p, ONLY: WriteField_p
  USE Bands, ONLY: jj_Nb_Caldyn,jj_Nb_vanleer
  USE mod_hallo
  USE Vampir
  USE times
  USE infotrac, ONLY: nqtot, iadv
  USE control_mod, ONLY: iapp_tracvl, day_step, planet_type, &
                         force_conserv_tracer
  USE comconst_mod, ONLY: dtvr
  USE planetary_operations_p, ONLY: planetary_tracer_amount_from_mass_p
  IMPLICIT NONE
  !
  include "dimensions.h"
  include "paramet.h"
  include "comdissip.h"
  include "comgeom2.h"

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

  REAL,SAVE :: pbaruc(ip1jmp1,llm)=0
  REAL,SAVE :: pbarvc(ip1jm,llm)=0
  REAL,SAVE :: massem(ip1jmp1,llm)
  REAL zdp(ip1jmp1)
  REAL,SAVE::pbarug(ip1jmp1,llm),pbarvg(ip1jm,llm),wg(ip1jmp1,llm) 
  REAL (kind=kind(1.d0)) :: t_initial, t_final, tps_cpu
  INTEGER ij,l,iq,iiq
  REAL zdpmin, zdpmax
  INTEGER,SAVE :: iadvtr=0
  !$OMP THREADPRIVATE(iadvtr)
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
  REAL,SAVE :: finmasse(ip1jmp1,llm)
  integer ijb,ije,ijb_u,ijb_v,ije_u,ije_v,j
  type(Request) :: Request_vanleer
  REAL,SAVE :: p_tmp( ip1jmp1,llmp1 )
  REAL,SAVE :: teta_tmp(ip1jmp1,llm)
  REAL,SAVE :: pk_tmp(ip1jmp1,llm)
  REAL :: totaltracer_old(nqtot),totaltracer_new(nqtot)
  REAL :: ratio

! Ehouarn : try to fix tracer conservation issues:
  if (force_conserv_tracer) then
    do iq=1,nqtot
       call planetary_tracer_amount_from_mass_p(masse,q(:,:,iq), &
                                          totaltracer_old(iq))
    enddo
  endif

  ijb_u=ij_begin
  ije_u=ij_end

  ijb_v=ij_begin-iip1
  ije_v=ij_end
  if (pole_nord) ijb_v=ij_begin
  if (pole_sud)  ije_v=ij_end-iip1

  IF(iadvtr.EQ.0) THEN
     !         CALL initial0(ijp1llm,pbaruc)
     !         CALL initial0(ijmllm,pbarvc)
     !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)        
     DO l=1,llm   
        pbaruc(ijb_u:ije_u,l)=0.
        pbarvc(ijb_v:ije_v,l)=0.
     ENDDO
     !$OMP END DO NOWAIT  
  ENDIF

  !   accumulation des flux de masse horizontaux
  !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
  DO l=1,llm
     DO ij = ijb_u,ije_u
        pbaruc(ij,l) = pbaruc(ij,l) + pbaru(ij,l)
     ENDDO
     DO ij = ijb_v,ije_v
        pbarvc(ij,l) = pbarvc(ij,l) + pbarv(ij,l)
     ENDDO
  ENDDO
  !$OMP END DO NOWAIT

  !   selection de la masse instantanee des mailles avant le transport.
  IF(iadvtr.EQ.0) THEN

     !         CALL SCOPY(ip1jmp1*llm,masse,1,massem,1)
     ijb=ij_begin
     ije=ij_end

     !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
     DO l=1,llm
        massem(ijb:ije,l)=masse(ijb:ije,l)
     ENDDO
     !$OMP END DO NOWAIT

     !cc         CALL filtreg ( massem ,jjp1, llm,-2, 2, .TRUE., 1 )
     !
  ENDIF ! of IF(iadvtr.EQ.0)

  iadvtr   = iadvtr+1

  !$OMP MASTER
  iapptrac = iadvtr
  !$OMP END MASTER

  !   Test pour savoir si on advecte a ce pas de temps

  IF ( iadvtr.EQ.iapp_tracvl ) THEN
     !$OMP MASTER
     call suspend_timer(timer_caldyn)
     !$OMP END MASTER

     ijb=ij_begin
     ije=ij_end


     !c   ..  Modif P.Le Van  ( 20/12/97 )  ....
     !c

     !   traitement des flux de masse avant advection.
     !     1. calcul de w
     !     2. groupement des mailles pres du pole.

     CALL groupe_p( massem, pbaruc,pbarvc, pbarug,pbarvg,wg )

     !$OMP BARRIER

     !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
     DO l=1,llmp1
        p_tmp(ijb:ije,l)=p(ijb:ije,l)
     ENDDO
     !$OMP END DO NOWAIT

     !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
     DO l=1,llm
        pk_tmp(ijb:ije,l)=pk(ijb:ije,l)
        teta_tmp(ijb:ije,l)=teta(ijb:ije,l)
     ENDDO
     !$OMP END DO NOWAIT

     !$OMP MASTER
     call VTb(VTHallo)
     !$OMP END MASTER

     call Register_SwapFieldHallo(pbarug,pbarug,ip1jmp1,llm, &
          jj_Nb_vanleer,0,0,Request_vanleer)
     call Register_SwapFieldHallo(pbarvg,pbarvg,ip1jm,llm, &
          jj_Nb_vanleer,1,0,Request_vanleer)
     call Register_SwapFieldHallo(massem,massem,ip1jmp1,llm, &
          jj_Nb_vanleer,0,0,Request_vanleer)
     call Register_SwapFieldHallo(wg,wg,ip1jmp1,llm, &
          jj_Nb_vanleer,0,0,Request_vanleer)
     call Register_SwapFieldHallo(teta_tmp,teta_tmp,ip1jmp1,llm, &
          jj_Nb_vanleer,1,1,Request_vanleer)
     call Register_SwapFieldHallo(p_tmp,p_tmp,ip1jmp1,llmp1, &
          jj_Nb_vanleer,1,1,Request_vanleer)
     call Register_SwapFieldHallo(pk_tmp,pk_tmp,ip1jmp1,llm, &
          jj_Nb_vanleer,1,1,Request_vanleer)
     do j=1,nqtot
        call Register_SwapFieldHallo(q(1,1,j),q(1,1,j),ip1jmp1,llm, &
             jj_nb_vanleer,0,0,Request_vanleer)
     enddo

     call SendRequest(Request_vanleer)
     !$OMP BARRIER
     call WaitRequest(Request_vanleer)


     !$OMP BARRIER
     !$OMP MASTER      
     call SetDistrib(jj_nb_vanleer)
     call VTe(VTHallo)
     call VTb(VTadvection)
     call start_timer(timer_vanleer)
     !$OMP END MASTER
     !$OMP BARRIER

     ! ... Flux de masse diagnostiques traceurs
     ijb=ij_begin
     ije=ij_end
     flxw(ijb:ije,1:llm)=wg(ijb:ije,1:llm)/REAL(iapp_tracvl)

     !  test sur l'eventuelle creation de valeurs negatives de la masse
     ijb=ij_begin
     ije=ij_end
     if (pole_nord) ijb=ij_begin+iip1
     if (pole_sud) ije=ij_end-iip1

     !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)          
     DO l=1,llm-1
        DO ij = ijb+1,ije
           zdp(ij) =    pbarug(ij-1,l)   - pbarug(ij,l) &
                - pbarvg(ij-iip1,l) + pbarvg(ij,l) &
                +       wg(ij,l+1)  - wg(ij,l)
        ENDDO

        !            CALL SCOPY( jjm -1 ,zdp(iip1+iip1),iip1,zdp(iip2),iip1 )
        ! ym  ---> pourquoi jjm-1 et non jjm ? a cause du pole ?

        do ij=ijb,ije-iip1+1,iip1
           zdp(ij)=zdp(ij+iip1-1)
        enddo

        DO ij = ijb,ije
           zdp(ij)= zdp(ij)*dtvr/ massem(ij,l) 
        ENDDO


        !            CALL minmax ( ip1jm-iip1, zdp(iip2), zdpmin,zdpmax )
        !  ym ---> eventuellement a revoir
        CALL minmax ( ije-ijb+1, zdp(ijb), zdpmin,zdpmax )

        IF(MAX(ABS(zdpmin),ABS(zdpmax)).GT.0.5) THEN
           PRINT*,'WARNING DP/P l=',l,'  MIN:',zdpmin, &
                '   MAX:', zdpmax
        ENDIF

     ENDDO
     !$OMP END DO NOWAIT

     !-------------------------------------------------------------------
     !   Advection proprement dite (Modification Le Croller (07/2001)
     !-------------------------------------------------------------------

     !----------------------------------------------------
     !        Calcul des moyennes basées sur la masse
     !----------------------------------------------------

     !ym      ----> Normalement, inutile pour les schémas classiques
     !ym      ----> Revérifier lors de la parallélisation des autres schemas

     !ym          call massbar_p(massem,massebx,masseby)  

     call vlspltgen_p( q,iadv, 2., massem, wg , &
          pbarug,pbarvg,dtvr,p_tmp,pk_tmp,teta_tmp )


!!! ATTENTION !!!! TOUT CE QUI EST ENTRE ICI ET 1234 EST OBSOLETE !!!!!!!
     GOTO 1234
!!! ATTENTION !!!!

     !-----------------------------------------------------------
     !     Appel des sous programmes d'advection
     !-----------------------------------------------------------
     do iq=1,nqtot
        !        call clock(t_initial)
        if(iadv(iq) == 0) cycle 
        !   ----------------------------------------------------------------
        !   Schema de Van Leer I MUSCL
        !   ----------------------------------------------------------------
        if(iadv(iq).eq.10) THEN

           call vlsplt_p(q(1,1,iq),2.,massem,wg,pbarug,pbarvg,dtvr)

           !   ----------------------------------------------------------------
           !   Schema "pseudo amont" + test sur humidite specifique
           !    pour la vapeur d'eau. F. Codron
           !   ----------------------------------------------------------------
        else if(iadv(iq).eq.14) then
           !
           !ym        stop 'advtrac : appel à vlspltqs :schema non parallelise'
           CALL vlspltqs_p( q(1,1,1), 2., massem, wg , &
                pbarug,pbarvg,dtvr,p_tmp,pk_tmp,teta_tmp )
           !   ----------------------------------------------------------------
           !   Schema de Frederic Hourdin
           !   ----------------------------------------------------------------
        else if(iadv(iq).eq.12) then
           stop 'advtrac : schema non parallelise'
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
           stop 'advtrac : schema non parallelise'
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
           stop 'advtrac : schema non parallelise'

           call pentes_ini (q(1,1,iq),wg,massem,pbarug,pbarvg,0)

           !   ----------------------------------------------------------------
           !   Schema de Prather
           !   ----------------------------------------------------------------
        else if (iadv(iq).eq.30) then
           stop 'advtrac : schema non parallelise'
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

           stop 'advtrac : schema non parallelise'

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

!!! ATTENTION !!!!
1234 CONTINUE
!!! ATTENTION !!!! LE CODE REPREND ICI !!!!!!!!

     !$OMP BARRIER

     if (planet_type=="earth") then

        ijb=ij_begin
        ije=ij_end

        !$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
        DO l = 1, llm
           DO ij = ijb, ije
              finmasse(ij,l) =  p(ij,l) - p(ij,l+1) 
           ENDDO
        ENDDO
        !$OMP END DO

        CALL qminimum_p( q, 2, finmasse )

     endif ! of if (planet_type=="earth")

     ! Ehouarn : try to fix tracer conservation after tracer advection
     if (force_conserv_tracer) then
       ijb=ij_begin
       ije=ij_end
       do iq=1,nqtot
         call planetary_tracer_amount_from_mass_p(masse,q(:,:,iq), &
                                     totaltracer_new(iq))
         ratio=totaltracer_old(iq)/totaltracer_new(iq)
         q(ijb:ije,1:llm,iq)=q(ijb:ije,1:llm,iq)*ratio
       enddo
     endif !of if (force_conserv_tracer)

        !------------------------------------------------------------------
        !   on reinitialise a zero les flux de masse cumules
        !---------------------------------------------------
        !          iadvtr=0

        !$OMP MASTER
	call VTe(VTadvection)
        call stop_timer(timer_vanleer)
        call VTb(VThallo)
        !$OMP END MASTER


	do j=1,nqtot
           call Register_SwapFieldHallo(q(1,1,j),q(1,1,j),ip1jmp1,llm, &
                jj_nb_caldyn,0,0,Request_vanleer)
        enddo

        call Register_SwapFieldHallo(flxw,flxw,ip1jmp1,llm, &
             jj_nb_caldyn,0,0,Request_vanleer)

        call SendRequest(Request_vanleer)
        !$OMP BARRIER
        call WaitRequest(Request_vanleer)      

        !$OMP BARRIER
        !$OMP MASTER
        call SetDistrib(jj_nb_caldyn)
	call VTe(VThallo)
	call resume_timer(timer_caldyn)
 !$OMP END MASTER
 !$OMP BARRIER	
        iadvtr=0
  ENDIF ! if iadvtr.EQ.iapp_tracvl

END SUBROUTINE advtrac_p
