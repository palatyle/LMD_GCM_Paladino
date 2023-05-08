










!
! $Id: sortvarc0.F 1403 2010-07-01 09:02:53Z fairhead $
!
      SUBROUTINE sortvarc0
     $(itau,ucov,teta,ps,masse,pk,phis,vorpot,phi,bern,dp,time ,
     $ vcov)

      USE comconst_mod, ONLY: daysec,dtvr,rad,g,omeg
      USE ener_mod, ONLY: etot0,ptot0,ztot0,stot0,ang0,rmsv,rmsdpdt

      IMPLICIT NONE

c=======================================================================
c
c   Auteur:    P. Le Van
c   -------
c
c   Objet:
c   ------
c
c   sortie des variables de controle
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

      INTEGER itau
      REAL ucov(ip1jmp1,llm),teta(ip1jmp1,llm),masse(ip1jmp1,llm)
      REAL vcov(ip1jm,llm)
      REAL ps(ip1jmp1),phis(ip1jmp1)
      REAL vorpot(ip1jm,llm)
      REAL phi(ip1jmp1,llm),bern(ip1jmp1,llm)
      REAL dp(ip1jmp1)
      REAL time
      REAL pk(ip1jmp1,llm)

c   Local:
c   ------

      REAL vor(ip1jm),bernf(ip1jmp1,llm),ztotl(llm)
      REAL etotl(llm),stotl(llm),rmsvl(llm),angl(llm),ge(ip1jmp1)
      REAL cosphi(ip1jm),omegcosp(ip1jm)
      REAL dtvrs1j,rjour,heure,radsg,radomeg
      REAL rday, massebxy(ip1jm,llm)
      INTEGER  l, ij, imjmp1

      REAL       SSUM
      integer  ismin,ismax

c-----------------------------------------------------------------------

       dtvrs1j   = dtvr/daysec
       rjour     = REAL( INT( itau * dtvrs1j ))
       heure     = ( itau*dtvrs1j-rjour ) * 24.
       imjmp1    = iim * jjp1
       IF(ABS(heure - 24.).LE.0.0001 ) heure = 0.
c
       CALL massbarxy ( masse, massebxy )

c   .....  Calcul  de  rmsdpdt  .....

       ge=dp*dp

       rmsdpdt = SSUM(ip1jmp1,ge,1) - SSUM(jjp1,ge,iip1)
c
       rmsdpdt = daysec* 1.e-2 * SQRT(rmsdpdt/imjmp1) 

       CALL SCOPY( ijp1llm,bern,1,bernf,1 )
       CALL filtreg(bernf,jjp1,llm,-2,2,.TRUE.,1)

c   .....  Calcul du moment  angulaire   .....

       radsg    = rad /g
       radomeg  = rad * omeg
c
       DO ij=iip2,ip1jm
          cosphi( ij ) = COS(rlatu((ij-1)/iip1+1))
          omegcosp(ij) = radomeg   * cosphi(ij)
       ENDDO

c  ...  Calcul  de l'energie,de l'enstrophie,de l'entropie et de rmsv  .

       DO l=1,llm
          DO ij = 1,ip1jm
             vor(ij)=vorpot(ij,l)*vorpot(ij,l)*massebxy(ij,l)
          ENDDO
          ztotl(l)=(SSUM(ip1jm,vor,1)-SSUM(jjm,vor,iip1))

          DO ij = 1,ip1jmp1
             ge(ij)= masse(ij,l)*(phis(ij)+teta(ij,l)*pk(ij,l)  +
     s        bernf(ij,l)-phi(ij,l))
          ENDDO
          etotl(l) = SSUM(ip1jmp1,ge,1) - SSUM(jjp1,ge,iip1)

          DO   ij   = 1, ip1jmp1
             ge(ij) = masse(ij,l)*teta(ij,l)
          ENDDO
          stotl(l)= SSUM(ip1jmp1,ge,1) - SSUM(jjp1,ge,iip1)

          DO ij=1,ip1jmp1
             ge(ij)=masse(ij,l)*AMAX1(bernf(ij,l)-phi(ij,l),0.)
          ENDDO
          rmsvl(l)=2.*(SSUM(ip1jmp1,ge,1)-SSUM(jjp1,ge,iip1))

          DO ij =iip2,ip1jm
             ge(ij)=(ucov(ij,l)/cu(ij)+omegcosp(ij))*masse(ij,l) *
     *               cosphi(ij)
          ENDDO
          angl(l) = rad *
     s    (SSUM(ip1jm-iip1,ge(iip2),1)-SSUM(jjm-1,ge(iip2),iip1))
      ENDDO

          DO ij=1,ip1jmp1
            ge(ij)= ps(ij)*aire(ij)
          ENDDO
      ptot0  = SSUM(ip1jmp1,ge,1)-SSUM(jjp1,ge,iip1)
      etot0  = SSUM(     llm, etotl, 1 )
      ztot0  = SSUM(     llm, ztotl, 1 )
      stot0  = SSUM(     llm, stotl, 1 )
      rmsv   = SSUM(     llm, rmsvl, 1 )
      ang0   = SSUM(     llm,  angl, 1 )

      rday = REAL(INT (time ))
c
      PRINT 3500, itau, rday, heure, time
      PRINT *, ptot0,etot0,ztot0,stot0,ang0

3500   FORMAT(10("*"),4x,'pas',i7,5x,'jour',f5.0,'heure',f5.1,4x 
     *   ,'date',f10.5,4x,10("*"))
      RETURN
      END

