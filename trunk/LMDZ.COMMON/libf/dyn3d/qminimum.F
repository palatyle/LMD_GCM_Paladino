!
! $Header$
!
      SUBROUTINE qminimum( q,nqtot,deltap )

      USE infotrac, ONLY: ok_isotopes,ntraciso,iqiso,ok_iso_verif
      IMPLICIT none
c
c  -- Objet : Traiter les valeurs trop petites (meme negatives)
c             pour l'eau vapeur et l'eau liquide
c
#include "dimensions.h"
#include "paramet.h"
c
      INTEGER nqtot
      REAL q(ip1jmp1,llm,nqtot), deltap(ip1jmp1,llm)
c
      INTEGER iq_vap, iq_liq
      PARAMETER ( iq_vap = 1 ) ! indice pour l'eau vapeur
      PARAMETER ( iq_liq = 2 ) ! indice pour l'eau liquide
      REAL seuil_vap, seuil_liq
      PARAMETER ( seuil_vap = 1.0e-10 ) ! seuil pour l'eau vapeur
      PARAMETER ( seuil_liq = 1.0e-11 ) ! seuil pour l'eau liquide
c
c  NB. ....( Il est souhaitable mais non obligatoire que les valeurs des
c            parametres seuil_vap, seuil_liq soient pareilles a celles 
c            qui  sont utilisees dans la routine    ADDFI       )
c     .................................................................
c
      INTEGER i, k, iq
      REAL zx_defau, zx_abc, zx_pump(ip1jmp1), pompe

      real zx_defau_diag(ip1jmp1,llm,2) 
      real q_follow(ip1jmp1,llm,2)
c
      REAL SSUM
c
      INTEGER imprim
      SAVE imprim
      DATA imprim /0/
      !INTEGER ijb,ije
      !INTEGER Index_pump(ij_end-ij_begin+1)
      !INTEGER nb_pump
      INTEGER ixt
c
c Quand l'eau liquide est trop petite (ou negative), on prend
c l'eau vapeur de la meme couche et la convertit en eau liquide
c (sans changer la temperature !)
c

        if (ok_iso_verif) then
           call check_isotopes_seq(q,ip1jmp1,'qminimum 52')   
        endif !if (ok_iso_verif) then     

      zx_defau_diag(:,:,:)=0.0
      q_follow(:,:,1:2)=q(:,:,1:2)  
      DO 1000 k = 1, llm
        DO 1040 i = 1, ip1jmp1
          if (seuil_liq - q(i,k,iq_liq) .gt. 0.d0 ) then

              if (ok_isotopes) then
                 zx_defau_diag(i,k,iq_liq)=AMAX1
     :               ( seuil_liq - q(i,k,iq_liq), 0.0 )
              endif !if (ok_isotopes) then

             q(i,k,iq_vap) = q(i,k,iq_vap) + q(i,k,iq_liq) - seuil_liq
             q(i,k,iq_liq) = seuil_liq
           endif
 1040   CONTINUE
 1000 CONTINUE
c
c Quand l'eau vapeur est trop faible (ou negative), on complete
c le defaut en prennant de l'eau vapeur de la couche au-dessous.
c
      iq = iq_vap
c
      DO k = llm, 2, -1
ccc      zx_abc = dpres(k) / dpres(k-1)
        DO i = 1, ip1jmp1
          if ( seuil_vap - q(i,k,iq) .gt. 0.d0 ) then

            if (ok_isotopes) then
              zx_defau_diag(i,k,iq)=AMAX1( seuil_vap - q(i,k,iq), 0.0 )
            endif !if (ok_isotopes) then

            q(i,k-1,iq) =  q(i,k-1,iq) - ( seuil_vap - q(i,k,iq) ) *
     &                     deltap(i,k) / deltap(i,k-1)
            q(i,k,iq)   =  seuil_vap  
          endif
        ENDDO
      ENDDO
c
c Quand il s'agit de la premiere couche au-dessus du sol, on
c doit imprimer un message d'avertissement (saturation possible).
c
      DO i = 1, ip1jmp1
         zx_pump(i) = AMAX1( 0.0, seuil_vap - q(i,1,iq) )
         q(i,1,iq)  = AMAX1( q(i,1,iq), seuil_vap )
      ENDDO
      pompe = SSUM(ip1jmp1,zx_pump,1)
      IF (imprim.LE.500 .AND. pompe.GT.0.0) THEN
         WRITE(6,'(1x,"ATT!:on pompe de l eau au sol",e15.7)') pompe
         DO i = 1, ip1jmp1
            IF (zx_pump(i).GT.0.0) THEN
               imprim = imprim + 1
               PRINT*,'QMINIMUM:  en ',i,zx_pump(i)
            ENDIF
         ENDDO
      ENDIF

      !write(*,*) 'qminimum 128'
      if (ok_isotopes) then
      ! CRisi: traiter de même les traceurs d'eau
      ! Mais il faut les prendre à l'envers pour essayer de conserver la
      ! masse.
      ! 1) pompage dans le sol  
      ! On suppose que ce pompage se fait sans isotopes -> on ne modifie
      ! rien ici et on croise les doigts pour que ça ne soit pas trop
      ! génant
      DO i = 1,ip1jmp1
        if (zx_pump(i).gt.0.0) then
          q_follow(i,1,iq_vap)=q_follow(i,1,iq_vap)+zx_pump(i)
        endif !if (zx_pump(i).gt.0.0) then
      enddo !DO i = 1,ip1jmp1

      ! 2) transfert de vap vers les couches plus hautes
      !write(*,*) 'qminimum 139'
      do k=2,llm
        DO i = 1,ip1jmp1
          if (zx_defau_diag(i,k,iq_vap).gt.0.0) then             
              ! on ajoute la vapeur en k              
              do ixt=1,ntraciso
               q(i,k,iqiso(ixt,iq_vap))=q(i,k,iqiso(ixt,iq_vap))
     :              +zx_defau_diag(i,k,iq_vap)
     :              *q(i,k-1,iqiso(ixt,iq_vap))/q_follow(i,k-1,iq_vap)
                
              ! et on la retranche en k-1
               q(i,k-1,iqiso(ixt,iq_vap))=q(i,k-1,iqiso(ixt,iq_vap))
     :              -zx_defau_diag(i,k,iq_vap)
     :              *deltap(i,k)/deltap(i,k-1)
     :              *q(i,k-1,iqiso(ixt,iq_vap))/q_follow(i,k-1,iq_vap)

              enddo !do ixt=1,niso
              q_follow(i,k,iq_vap)=   q_follow(i,k,iq_vap)
     :               +zx_defau_diag(i,k,iq_vap)
              q_follow(i,k-1,iq_vap)=   q_follow(i,k-1,iq_vap)
     :               -zx_defau_diag(i,k,iq_vap)
     :              *deltap(i,k)/deltap(i,k-1)
          endif !if (zx_defau_diag(i,k,iq_vap).gt.0.0) then
        enddo !DO i = 1, ip1jmp1        
       enddo !do k=2,llm

        if (ok_iso_verif) then     
           call check_isotopes_seq(q,ip1jmp1,'qminimum 168')
        endif !if (ok_iso_verif) then
        
      
        ! 3) transfert d'eau de la vapeur au liquide
        !write(*,*) 'qminimum 164'
        do k=1,llm
        DO i = 1,ip1jmp1
          if (zx_defau_diag(i,k,iq_liq).gt.0.0) then

              ! on ajoute eau liquide en k en k              
              do ixt=1,ntraciso
               q(i,k,iqiso(ixt,iq_liq))=q(i,k,iqiso(ixt,iq_liq))
     :              +zx_defau_diag(i,k,iq_liq)
     :              *q(i,k,iqiso(ixt,iq_vap))/q_follow(i,k,iq_vap)
              ! et on la retranche à la vapeur en k
               q(i,k,iqiso(ixt,iq_vap))=q(i,k,iqiso(ixt,iq_vap))
     :              -zx_defau_diag(i,k,iq_liq)
     :              *q(i,k,iqiso(ixt,iq_vap))/q_follow(i,k,iq_vap)   
              enddo !do ixt=1,niso
              q_follow(i,k,iq_liq)=   q_follow(i,k,iq_liq)
     :               +zx_defau_diag(i,k,iq_liq)
              q_follow(i,k,iq_vap)=   q_follow(i,k,iq_vap)
     :               -zx_defau_diag(i,k,iq_liq)
          endif !if (zx_defau_diag(i,k,iq_vap).gt.0.0) then
        enddo !DO i = 1, ip1jmp1
       enddo !do k=2,llm  

        if (ok_iso_verif) then
           call check_isotopes_seq(q,ip1jmp1,'qminimum 197')
        endif !if (ok_iso_verif) then

      endif !if (ok_isotopes) then
      !write(*,*) 'qminimum 188'
      
c
      RETURN
      END
