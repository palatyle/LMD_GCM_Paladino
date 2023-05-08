










      SUBROUTINE nxgrarot_p (klevel,xcov, ycov, lr, grx_out, gry_out )
c   ***********************************************************
c
c    Auteur :  P.Le Van  
c
c                                 lr
c      calcul de  ( nXgrad (rot) )   du vect. v  ....
c
c       xcov et ycov  etant les compos. covariantes de  v
c   ***********************************************************
c     xcov , ycov et lr  sont des arguments  d'entree pour le s-prog
c      grx   et  gry     sont des arguments de sortie pour le s-prog
c
c
      USE parallel_lmdz
      USE times
      USE write_field_p
      IMPLICIT NONE
c
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
c
      INTEGER klevel
      REAL xcov( ip1jmp1,klevel ), ycov( ip1jm,klevel )
      REAL  grx_out( ip1jmp1,klevel ),  gry_out( ip1jm,klevel )
      REAL,SAVE ::  grx( ip1jmp1,llm ),  gry( ip1jm,llm )

c
      REAL,SAVE :: rot(ip1jm,llm)

      INTEGER l,ij,iter,lr
c
      INTEGER ijb,ije,jjb,jje
c
c
c      CALL SCOPY ( ip1jmp1*klevel, xcov, 1, grx, 1 )
c      CALL SCOPY (  ip1jm*klevel, ycov, 1, gry, 1 )
c
      ijb=ij_begin
      ije=ij_end
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)  
      DO l = 1, klevel
        grx(ijb:ije,l)=xcov(ijb:ije,l)
      ENDDO 
c$OMP END DO NOWAIT      

      if(pole_sud) ije=ij_end-iip1
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)      
      DO l = 1, klevel
        gry(ijb:ije,l)=ycov(ijb:ije,l)
      ENDDO
c$OMP END DO NOWAIT
      
      DO 10 iter = 1,lr
c$OMP BARRIER
c$OMP MASTER
      call suspend_timer(timer_dissip)
      call exchange_Hallo(grx,ip1jmp1,llm,0,1)
      call resume_timer(timer_dissip)
c$OMP END MASTER
c$OMP BARRIER

      CALL  rotat_p (klevel,grx, gry, rot )
c      call write_field3d_p('rot',reshape(rot,(/iip1,jjm,llm/)))
      
      jjb=jj_begin
      jje=jj_end
      if (pole_sud) jje=jj_end-1
      CALL filtreg_p( rot,jjb,jje, jjm, klevel, 2,1, .false.,2)

c$OMP BARRIER
c$OMP MASTER
      call suspend_timer(timer_dissip)
      call exchange_Hallo(rot,ip1jm,llm,1,0)
      call resume_timer(timer_dissip)
c$OMP END MASTER
c$OMP BARRIER
      
      CALL nxgrad_p (klevel,rot, grx, gry )
c
c$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO 5  l = 1, klevel
      if(pole_sud) ije=ij_end-iip1
      DO 2 ij = ijb, ije
      gry_out( ij,l ) = - gry( ij,l ) * crot
   2  CONTINUE
      if(pole_sud) ije=ij_end
      DO 3 ij = ijb, ije
      grx_out( ij,l ) = - grx( ij,l ) * crot
   3  CONTINUE
   5  CONTINUE
c$OMP END DO NOWAIT
c      call write_field3d_p('grx',reshape(grx,(/iip1,jjp1,llm/)))
c      call write_field3d_p('gry',reshape(gry,(/iip1,jjm,llm/)))
c      stop
  10  CONTINUE
      RETURN
      END
