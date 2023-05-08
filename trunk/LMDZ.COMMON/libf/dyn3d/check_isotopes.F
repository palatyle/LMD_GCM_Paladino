        subroutine check_isotopes_seq(q,ip1jmp1,err_msg)
        USE infotrac, ONLY: nqtot,niso,nqo,ntraceurs_zone,ntraciso,
     &                      ok_isotopes,ok_isotrac,iqiso,use_iso,
     &                      indnum_fn_num,tnat,index_trac
        implicit none

#include "dimensions.h"

        ! inputs
        integer ip1jmp1
        real q(ip1jmp1,llm,nqtot)
        character*(*) err_msg ! message d''erreur à afficher

        ! locals
        integer ixt,phase,k,i,iq,iiso,izone,ieau,iqeau
        real xtractot,xiiso
        real borne
        real qmin
        real errmax ! erreur maximale en absolu.
        real errmaxrel ! erreur maximale en relatif autorisée
        real deltaDmax,deltaDmin
        real ridicule
        parameter (borne=1e19)
        parameter (errmax=1e-8)
        parameter (errmaxrel=1e-3)
        parameter (qmin=1e-11)
        parameter (deltaDmax=200.0,deltaDmin=-999.9)
        parameter (ridicule=1e-12)
        real deltaD

        if (ok_isotopes) then

        write(*,*) 'check_isotopes 31: err_msg=',err_msg
        ! verifier que rien n'est NaN
        do ixt=1,ntraciso
          do phase=1,nqo
            iq=iqiso(ixt,phase)
            do k=1,llm
              DO i = 1,ip1jmp1
                if ((q(i,k,iq).gt.-borne).and.
     :            (q(i,k,iq).lt.borne)) then
                else !if ((x(ixt,i,j).gt.-borne).and.
                  write(*,*) 'erreur detectee par iso_verif_noNaN:'
                  write(*,*) err_msg
                  write(*,*) 'q,i,k,iq=',q(i,k,iq),i,k,iq
                  write(*,*) 'borne=',borne
                  stop
                endif  !if ((x(ixt,i,j).gt.-borne).and.
              enddo !DO i = 1,ip1jmp1
            enddo !do k=1,llm
          enddo !do phase=1,nqo
        enddo !do ixt=1,ntraciso

        !write(*,*) 'check_isotopes 52'
        ! verifier que l'eau normale est OK
        if (use_iso(1)) then
          ixt=indnum_fn_num(1)
          do phase=1,nqo
            iq=iqiso(ixt,phase)
            do k=1,llm
            DO i = 1,ip1jmp1  
              if ((abs(q(i,k,phase)-q(i,k,iq)).gt.errmax).and.
     :          (abs((q(i,k,phase)-q(i,k,iq))/
     :           max(max(abs(q(i,k,phase)),abs(q(i,k,iq))),1e-18))
     :           .gt.errmaxrel)) then
                  write(*,*) 'erreur detectee par iso_verif_egalite:'
                  write(*,*) err_msg
                  write(*,*) 'ixt,phase=',ixt,phase
                  write(*,*) 'q,iq,i,k=',q(i,k,iq),iq,i,k
                  write(*,*) 'q(i,k,phase)=',q(i,k,phase)
                  stop
              endif !if ((abs(q(i,k,phase)-q(i,k,iq)).gt.errmax).and.
              ! bidouille pour éviter divergence:
              q(i,k,iq)= q(i,k,phase) 
            enddo ! DO i = 1,ip1jmp1
            enddo !do k=1,llm
          enddo ! do phase=1,nqo 
        endif !if (use_iso(1)) then
        
        !write(*,*) 'check_isotopes 78'
        ! verifier que HDO est raisonable
        if (use_iso(2)) then
          ixt=indnum_fn_num(2)
          do phase=1,nqo
            iq=iqiso(ixt,phase)
            do k=1,llm
            DO i = 1,ip1jmp1
            if (q(i,k,iq).gt.qmin) then
             deltaD=(q(i,k,iq)/q(i,k,phase)/tnat(2)-1)*1000
             if ((deltaD.gt.deltaDmax).or.(deltaD.lt.deltaDmin)) then
                  write(*,*) 'erreur detectee par iso_verif_aberrant:'
                  write(*,*) err_msg
                  write(*,*) 'ixt,phase=',ixt,phase
                  write(*,*) 'q,iq,i,k,=',q(i,k,iq),iq,i,k
                  write(*,*) 'q=',q(i,k,:)
                  write(*,*) 'deltaD=',deltaD
                  stop
             endif !if ((deltaD.gt.deltaDmax).or.(deltaD.lt.deltaDmin)) then
            endif !if (q(i,k,iq).gt.qmin) then
            enddo !DO i = 1,ip1jmp1
            enddo !do k=1,llm
          enddo ! do phase=1,nqo 
        endif !if (use_iso(2)) then

        !write(*,*) 'check_isotopes 103'
        ! verifier que O18 est raisonable
        if (use_iso(3)) then
          ixt=indnum_fn_num(3)
          do phase=1,nqo
            iq=iqiso(ixt,phase)
            do k=1,llm
            DO i = 1,ip1jmp1
            if (q(i,k,iq).gt.qmin) then
             deltaD=(q(i,k,iq)/q(i,k,phase)/tnat(3)-1)*1000
             if ((deltaD.gt.deltaDmax).or.(deltaD.lt.deltaDmin)) then
                  write(*,*) 'erreur detectee iso_verif_aberrant O18:'
                  write(*,*) err_msg
                  write(*,*) 'ixt,phase=',ixt,phase
                  write(*,*) 'q,iq,i,k,=',q(i,k,phase),iq,i,k
                  write(*,*) 'xt=',q(i,k,:)
                  write(*,*) 'deltaO18=',deltaD
                  stop
             endif !if ((deltaD.gt.deltaDmax).or.(deltaD.lt.deltaDmin)) then
            endif !if (q(i,k,iq).gt.qmin) then
            enddo !DO i = 1,ip1jmp1
            enddo !do k=1,llm
          enddo ! do phase=1,nqo 
        endif !if (use_iso(2)) then


        !write(*,*) 'check_isotopes 129'
        if (ok_isotrac) then

          if (use_iso(2).and.use_iso(1)) then
            do izone=1,ntraceurs_zone
             ixt=index_trac(izone,indnum_fn_num(2))
             ieau=index_trac(izone,indnum_fn_num(1))
             do phase=1,nqo
               iq=iqiso(ixt,phase)
               iqeau=iqiso(ieau,phase)
               do k=1,llm
                DO i = 1,ip1jmp1
                if (q(i,k,iq).gt.qmin) then
                 deltaD=(q(i,k,iq)/q(i,k,iqeau)/tnat(2)-1)*1000
                 if ((deltaD.gt.deltaDmax).or.
     &                   (deltaD.lt.deltaDmin)) then
                  write(*,*) 'erreur dans iso_verif_aberrant trac:'
                  write(*,*) err_msg
                  write(*,*) 'izone,phase=',izone,phase
                  write(*,*) 'ixt,ieau=',ixt,ieau
                  write(*,*) 'q,iq,i,k,=',q(i,k,iq),iq,i,k
                  write(*,*) 'deltaD=',deltaD
                  stop
                 endif !if ((deltaD.gt.deltaDmax).or.
                endif !if (q(i,k,iq).gt.qmin) then
                enddo !DO i = 1,ip1jmp1
                enddo  ! do k=1,llm
              enddo ! do phase=1,nqo    
            enddo !do izone=1,ntraceurs_zone
          endif !if (use_iso(2).and.use_iso(1)) then

          do iiso=1,niso
           do phase=1,nqo
              iq=iqiso(iiso,phase)
              do k=1,llm
                DO i = 1,ip1jmp1
                   xtractot=0.0
                   xiiso=q(i,k,iq)
                   do izone=1,ntraceurs_zone
                      iq=iqiso(index_trac(izone,iiso),phase)
                      xtractot=xtractot+ q(i,k,iq)
                   enddo !do izone=1,ntraceurs_zone
                   if ((abs(xtractot-xiiso).gt.errmax).and.
     :                  (abs(xtractot-xiiso)/
     :                  max(max(abs(xtractot),abs(xiiso)),1e-18)
     :                  .gt.errmaxrel)) then
                  write(*,*) 'erreur detectee par iso_verif_traceurs:'
                  write(*,*) err_msg
                  write(*,*) 'iiso,phase=',iiso,phase
                  write(*,*) 'i,k,=',i,k
                  write(*,*) 'q(i,k,:)=',q(i,k,:)
                  stop
                 endif !if ((abs(q(i,k,phase)-q(i,k,iq)).gt.errmax).and.
                  
                 ! bidouille pour éviter divergence:
                 if (abs(xtractot).gt.ridicule) then
                   do izone=1,ntraceurs_zone
                     ixt=index_trac(izone,iiso) 
                     q(i,k,iq)=q(i,k,iq)/xtractot*xiiso
                   enddo !do izone=1,ntraceurs_zone                
                  endif !if ((abs(xtractot).gt.ridicule) then
                enddo !DO i = 1,ip1jmp1
              enddo !do k=1,llm
           enddo !do phase=1,nqo
          enddo !do iiso=1,niso

        endif !if (ok_isotrac) then

        endif ! if (ok_isotopes)
        !write(*,*) 'check_isotopes 198'
        
        end

