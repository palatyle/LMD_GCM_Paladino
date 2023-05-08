










!
! $Header$
!
      SUBROUTINE advect_p(ucov,vcov,teta,w,massebx,masseby,du,dv,dteta)
      USE parallel_lmdz
      USE write_field_p
      USE comconst_mod, ONLY: daysec
      USE logic_mod, ONLY: conser
      IMPLICIT NONE
c=======================================================================
c
c   Auteurs:  P. Le Van , Fr. Hourdin  .
c   -------
c
c   Objet:
c   ------
c
c   *************************************************************
c   .... calcul des termes d'advection vertic.pour u,v,teta,q ...
c   *************************************************************
c        ces termes sont ajoutes a du,dv,dteta et dq .
c  Modif F.Forget 03/94 : on retire q de advect
c
c=======================================================================
c-----------------------------------------------------------------------
c   Declarations:
c   -------------

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
!CDK comgeom
      COMMON/comgeom/                                                   &
     & cu(ip1jmp1),cv(ip1jm),unscu2(ip1jmp1),unscv2(ip1jm),             &
     & aire(ip1jmp1),airesurg(ip1jmp1),aireu(ip1jmp1),                  &
     & airev(ip1jm),unsaire(ip1jmp1),apoln,apols,                       &
     & unsairez(ip1jm),airuscv2(ip1jm),airvscu2(ip1jm),                 &
     & aireij1(ip1jmp1),aireij2(ip1jmp1),aireij3(ip1jmp1),              &
     & aireij4(ip1jmp1),alpha1(ip1jmp1),alpha2(ip1jmp1),                &
     & alpha3(ip1jmp1),alpha4(ip1jmp1),alpha1p2(ip1jmp1),               &
     & alpha1p4(ip1jmp1),alpha2p3(ip1jmp1),alpha3p4(ip1jmp1),           &
     & fext(ip1jm),constang(ip1jmp1),rlatu(jjp1),rlatv(jjm),            &
     & rlonu(iip1),rlonv(iip1),cuvsurcv(ip1jm),cvsurcuv(ip1jm),         &
     & cvusurcu(ip1jmp1),cusurcvu(ip1jmp1),cuvscvgam1(ip1jm),           &
     & cuvscvgam2(ip1jm),cvuscugam1(ip1jmp1),                           &
     & cvuscugam2(ip1jmp1),cvscuvgam(ip1jm),cuscvugam(ip1jmp1),         &
     & unsapolnga1,unsapolnga2,unsapolsga1,unsapolsga2,                 &
     & unsair_gam1(ip1jmp1),unsair_gam2(ip1jmp1),unsairz_gam(ip1jm),    &
     & aivscu2gam(ip1jm),aiuscv2gam(ip1jm),xprimu(iip1),xprimv(iip1)

!
        REAL                                                            &
     & cu,cv,unscu2,unscv2,aire,airesurg,aireu,airev,unsaire,apoln     ,&
     & apols,unsairez,airuscv2,airvscu2,aireij1,aireij2,aireij3,aireij4,&
     & alpha1,alpha2,alpha3,alpha4,alpha1p2,alpha1p4,alpha2p3,alpha3p4 ,&
     & fext,constang,rlatu,rlatv,rlonu,rlonv,cuvscvgam1,cuvscvgam2     ,&
     & cvuscugam1,cvuscugam2,cvscuvgam,cuscvugam,unsapolnga1,unsapolnga2&
     & ,unsapolsga1,unsapolsga2,unsair_gam1,unsair_gam2,unsairz_gam    ,&
     & aivscu2gam ,aiuscv2gam,cuvsurcv,cvsurcuv,cvusurcu,cusurcvu,xprimu&
     & , xprimv
!

c   Arguments:
c   ----------

      REAL vcov(ip1jm,llm),ucov(ip1jmp1,llm),teta(ip1jmp1,llm)
      REAL massebx(ip1jmp1,llm),masseby(ip1jm,llm),w(ip1jmp1,llm)
      REAL dv(ip1jm,llm),du(ip1jmp1,llm),dteta(ip1jmp1,llm)

c   Local:
c   ------

      REAL uav(ip1jmp1,llm),vav(ip1jm,llm),wsur2(ip1jmp1)
      REAL unsaire2(ip1jmp1), ge(ip1jmp1)
      REAL deuxjour, ww, gt, uu, vv

      INTEGER  ij,l,ijb,ije

      EXTERNAL  SSUM
      REAL      SSUM

c-----------------------------------------------------------------------
c   2. Calculs preliminaires:
c   -------------------------

      IF (conser)  THEN
         deuxjour = 2. * daysec

         DO   1  ij   = 1, ip1jmp1
         unsaire2(ij) = unsaire(ij) * unsaire(ij)
   1     CONTINUE
      END IF


c------------------  -yy ----------------------------------------------
c   .  Calcul de     u

      DO  l=1,llm
         
         ijb=ij_begin
         ije=ij_end
         if (pole_nord) ijb=ijb+iip1
         if (pole_sud)  ije=ije-iip1
         
c         DO    ij     = iip2, ip1jmp1
c            uav(ij,l) = 0.25 * ( ucov(ij,l) + ucov(ij-iip1,l) )
c         ENDDO

c         DO    ij     = iip2, ip1jm
c            uav(ij,l) = uav(ij,l) + uav(ij+iip1,l)
c         ENDDO
         
         DO    ij     = ijb, ije
                  
           uav(ij,l)=0.25*(ucov(ij,l)+ucov(ij-iip1,l))
     .	             +0.25*(ucov(ij+iip1,l)+ucov(ij,l))
         ENDDO
         
         if (pole_nord) then
           DO      ij         = 1, iip1
              uav(ij      ,l) = 0.
           ENDDO
         endif
         
         if (pole_sud) then
           DO      ij         = 1, iip1
              uav(ip1jm+ij,l) = 0.
           ENDDO
         endif
         
      ENDDO
      
c      call write_field3d_p('uav',reshape(uav,(/iip1,jjp1,llm/)))
      
c------------------  -xx ----------------------------------------------
c   .  Calcul de     v
      
      ijb=ij_begin
      ije=ij_end
      if (pole_sud)  ije=ij_end-iip1
      
      DO  l=1,llm
         
         DO    ij   = ijb+1, ije
           vav(ij,l) = 0.25 * ( vcov(ij,l) + vcov(ij-1,l) )
         ENDDO
         
         DO    ij   = ijb,ije,iip1
          vav(ij,l) = vav(ij+iim,l)
         ENDDO
         
         
         DO    ij   = ijb, ije-1
          vav(ij,l) = vav(ij,l) + vav(ij+1,l)
         ENDDO
         
         DO    ij       = ijb, ije, iip1
          vav(ij+iim,l) = vav(ij,l)
         ENDDO
         
      ENDDO
c       call write_field3d_p('vav',reshape(vav,(/iip1,jjm,llm/)))
c-----------------------------------------------------------------------


      
      DO 20 l = 1, llmm1


c       ......   calcul de  - w/2.    au niveau  l+1   .......
      ijb=ij_begin
      ije=ij_end+iip1
      if (pole_sud)  ije=ij_end
      
      DO 5   ij   = ijb, ije
      wsur2( ij ) = - 0.5 * w( ij,l+1 )
   5  CONTINUE


c     .....................     calcul pour  du     ..................
      
      ijb=ij_begin
      ije=ij_end
      if (pole_nord) ijb=ijb+iip1
      if (pole_sud)  ije=ije-iip1
         
      DO 6 ij = ijb ,ije-1
      ww        = wsur2 (  ij  )     + wsur2( ij+1 ) 
      uu        = 0.5 * ( ucov(ij,l) + ucov(ij,l+1) )
      du(ij,l)  = du(ij,l)   - ww * ( uu - uav(ij, l ) )/massebx(ij, l )
      du(ij,l+1)= du(ij,l+1) + ww * ( uu - uav(ij,l+1) )/massebx(ij,l+1)
   6  CONTINUE

c     .....  correction pour  du(iip1,j,l)  ........
c     .....     du(iip1,j,l)= du(1,j,l)   .....

CDIR$ IVDEP
      DO   7  ij   = ijb+iip1-1, ije, iip1
      du( ij, l  ) = du( ij -iim, l  )
      du( ij,l+1 ) = du( ij -iim,l+1 )
   7  CONTINUE

c     .................    calcul pour   dv      .....................
      ijb=ij_begin
      ije=ij_end
      if (pole_sud)  ije=ij_end-iip1
      
      DO 8 ij = ijb, ije
      ww        = wsur2( ij+iip1 )   + wsur2( ij )
      vv        = 0.5 * ( vcov(ij,l) + vcov(ij,l+1) )
      dv(ij,l)  = dv(ij, l ) - ww * (vv - vav(ij, l ) )/masseby(ij, l )
      dv(ij,l+1)= dv(ij,l+1) + ww * (vv - vav(ij,l+1) )/masseby(ij,l+1)
   8  CONTINUE

c

c     ............................................................
c     ...............    calcul pour   dh      ...................
c     ............................................................

c                       ---z
c       calcul de  - d( teta  * w )      qu'on ajoute a   dh
c                   ...............
        ijb=ij_begin
        ije=ij_end
        
        DO 15 ij = ijb, ije
         ww            = wsur2(ij) * (teta(ij,l) + teta(ij,l+1) )
         dteta(ij, l ) = dteta(ij, l )  -  ww
         dteta(ij,l+1) = dteta(ij,l+1)  +  ww
  15    CONTINUE

c ym ---> conser a voir plus tard

c      IF( conser)  THEN
c        
c        DO 17 ij = 1,ip1jmp1
c        ge(ij)   = wsur2(ij) * wsur2(ij) * unsaire2(ij)
c  17    CONTINUE
c        gt       = SSUM( ip1jmp1,ge,1 )
c        gtot(l)  = deuxjour * SQRT( gt/ip1jmp1 )
c      END IF

  20  CONTINUE
 
      RETURN
      END
