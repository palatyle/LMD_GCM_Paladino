










      SUBROUTINE gradiv_p(klevel, xcov, ycov, ld, gdx_out, gdy_out )
c
c    Auteur :   P. Le Van
c
c   ***************************************************************
c
c                                ld
c       calcul  de  (grad (div) )   du vect. v ....
c
c     xcov et ycov etant les composant.covariantes de v
c   ****************************************************************
c    xcov , ycov et ld  sont des arguments  d'entree pour le s-prog
c     gdx   et  gdy     sont des arguments de sortie pour le s-prog
c
c     
      USE parallel_lmdz
      USE times
      IMPLICIT NONE
c
!-----------------------------------------------------------------------
!   INCLUDE 'dimensions.h'
!
!   dimensions.h contient les dimensions du modele
!   ndm est tel que iim=2**ndm
!-----------------------------------------------------------------------

      INTEGER iim,jjm,llm,ndm

      PARAMETER (iim= 128,jjm=96,llm=23,ndm=1)

!-----------------------------------------------------------------------
!
! $Header$
!
!
!  ATTENTION!!!!: ce fichier include est compatible format fixe/format libre
!                 veillez  n'utiliser que des ! pour les commentaires
!                 et  bien positionner les & des lignes de continuation
!                 (les placer en colonne 6 et en colonne 73)
!
!
!-----------------------------------------------------------------------
!   INCLUDE 'paramet.h'

      INTEGER  iip1,iip2,iip3,jjp1,llmp1,llmp2,llmm1
      INTEGER  kftd,ip1jm,ip1jmp1,ip1jmi1,ijp1llm
      INTEGER  ijmllm,mvar
      INTEGER jcfil,jcfllm

      PARAMETER( iip1= iim+1,iip2=iim+2,iip3=iim+3                       &
     &    ,jjp1=jjm+1-1/jjm)
      PARAMETER( llmp1 = llm+1,  llmp2 = llm+2, llmm1 = llm-1 )
      PARAMETER( kftd  = iim/2 -ndm )
      PARAMETER( ip1jm  = iip1*jjm,  ip1jmp1= iip1*jjp1 )
      PARAMETER( ip1jmi1= ip1jm - iip1 )
      PARAMETER( ijp1llm= ip1jmp1 * llm, ijmllm= ip1jm * llm )
      PARAMETER( mvar= ip1jmp1*( 2*llm+1) + ijmllm )
      PARAMETER( jcfil=jjm/2+5, jcfllm=jcfil*llm )

!-----------------------------------------------------------------------
!
! $Header$
!
!  Attention : ce fichier include est compatible format fixe/format libre
!                 veillez à n'utiliser que des ! pour les commentaires
!                 et à bien positionner les & des lignes de continuation 
!                 (les placer en colonne 6 et en colonne 73)
!-----------------------------------------------------------------------
! INCLUDE comdissipn.h

      REAL  tetaudiv, tetaurot, tetah, cdivu, crot, cdivh
!
      COMMON/comdissipn/ tetaudiv(llm),tetaurot(llm),tetah(llm)   ,     &
     &                        cdivu,      crot,         cdivh

!
!    Les parametres de ce common proviennent des calculs effectues dans 
!             Inidissip  .
!
!-----------------------------------------------------------------------

      INTEGER klevel
c
      REAL xcov( ip1jmp1,klevel ), ycov( ip1jm,klevel )
      REAL,SAVE :: gdx( ip1jmp1,llm ),   gdy( ip1jm,llm )

      REAL gdx_out( ip1jmp1,klevel ),   gdy_out( ip1jm,klevel )

      REAL,SAVE ::  div(ip1jmp1,llm)

      INTEGER l,ij,iter,ld
c
      INTEGER ijb,ije,jjb,jje
c
c
c      CALL SCOPY( ip1jmp1*klevel,xcov,1,gdx,1 )
c      CALL SCOPY( ip1jm*klevel,  ycov,1,gdy,1 )
      
      ijb=ij_begin
      ije=ij_end

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l = 1,klevel
        gdx(ijb:ije,l)=xcov(ijb:ije,l)
      ENDDO
c$OMP END DO NOWAIT
      
      ijb=ij_begin
      ije=ij_end
      if(pole_sud) ije=ij_end-iip1

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l = 1,klevel
        gdy(ijb:ije,l)=ycov(ijb:ije,l)
      ENDDO
c$OMP END DO NOWAIT

c
      DO 10 iter = 1,ld

c$OMP BARRIER
c$OMP MASTER      
      call suspend_timer(timer_dissip)
      call exchange_Hallo(gdy,ip1jm,llm,1,0)
      call resume_timer(timer_dissip)
c$OMP END MASTER      
c$OMP BARRIER

      CALL  diverg_p( klevel,  gdx , gdy, div          )
      
      jjb=jj_begin
      jje=jj_end
      CALL filtreg_p( div,jjb,jje, jjp1, klevel, 2,1, .true.,2 )
      
c      call exchange_Hallo(div,ip1jmp1,llm,0,1)

c$OMP BARRIER
c$OMP MASTER       
      call suspend_timer(timer_dissip)
      call exchange_Hallo(div,ip1jmp1,llm,1,1)
      call resume_timer(timer_dissip)
c$OMP END MASTER
c$OMP BARRIER
      
      CALL    grad_p( klevel,  div, gdx, gdy           )
c

c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO 5  l = 1, klevel
      
      if(pole_sud) ije=ij_end
      DO 3 ij = ijb, ije
        gdx_out( ij,l ) = - gdx( ij,l ) * cdivu
   3  CONTINUE
   
      if(pole_sud) ije=ij_end-iip1
      DO 4 ij = ijb, ije
        gdy_out( ij,l ) = - gdy( ij,l ) * cdivu
   4  CONTINUE

   5  CONTINUE
c$OMP END DO NOWAIT
c
  10  CONTINUE
      RETURN
      END
