










!
! $Id: addfi.F 1446 2010-10-22 09:27:25Z emillour $
!
      SUBROUTINE addfi(pdt, leapf, forward,
     S          pucov, pvcov, pteta, pq   , pps ,
     S          pdufi, pdvfi, pdhfi,pdqfi, pdpfi  )

      USE infotrac, ONLY : nqtot
      USE control_mod, ONLY : planet_type
      USE comconst_mod, ONLY: kappa
      IMPLICIT NONE
c
c=======================================================================
c
c    Addition of the physical tendencies
c
c    Interface :
c    -----------
c
c      Input :
c      -------
c      pdt                    time step of integration
c      leapf                  logical
c      forward                logical
c      pucov(ip1jmp1,llm)     first component of the covariant velocity
c      pvcov(ip1ip1jm,llm)    second component of the covariant velocity
c      pteta(ip1jmp1,llm)     potential temperature
c      pts(ip1jmp1,llm)       surface temperature
c      pdufi(ip1jmp1,llm)     |
c      pdvfi(ip1jm,llm)       |   respective
c      pdhfi(ip1jmp1)         |      tendencies  (unit/s)
c      pdtsfi(ip1jmp1)        |
c
c      Output :
c      --------
c      pucov
c      pvcov
c      ph
c      pts
c
c
c=======================================================================
c
c-----------------------------------------------------------------------
c
c    0.  Declarations :
c    ------------------
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
c
c    Arguments :
c    -----------
c
      REAL,INTENT(IN) :: pdt ! time step for the integration (s)
c
      REAL,INTENT(INOUT) :: pvcov(ip1jm,llm) ! covariant meridional wind
      REAL,INTENT(INOUT) :: pucov(ip1jmp1,llm) ! covariant zonal wind
      REAL,INTENT(INOUT) :: pteta(ip1jmp1,llm) ! potential temperature
      REAL,INTENT(INOUT) :: pq(ip1jmp1,llm,nqtot) ! tracers
      REAL,INTENT(INOUT) :: pps(ip1jmp1) ! surface pressure (Pa)
c respective tendencies (.../s) to add
      REAL,INTENT(IN) :: pdvfi(ip1jm,llm)
      REAL,INTENT(IN) :: pdufi(ip1jmp1,llm)
      REAL,INTENT(IN) :: pdqfi(ip1jmp1,llm,nqtot)
      REAL,INTENT(IN) :: pdhfi(ip1jmp1,llm)
      REAL,INTENT(IN) :: pdpfi(ip1jmp1)
c
      LOGICAL,INTENT(IN) :: leapf,forward ! not used
c
c
c    Local variables :
c    -----------------
c
      REAL xpn(iim),xps(iim),tpn,tps
      INTEGER j,k,iq,ij
      REAL,PARAMETER :: qtestw = 1.0e-15
      REAL,PARAMETER :: qtestt = 1.0e-40

      REAL SSUM
c
c-----------------------------------------------------------------------

      DO k = 1,llm
         DO j = 1,ip1jmp1
            pteta(j,k)= pteta(j,k) + pdhfi(j,k) * pdt
         ENDDO
      ENDDO

      DO  k    = 1, llm
       DO  ij   = 1, iim
         xpn(ij) = aire(   ij   ) * pteta(  ij    ,k)
         xps(ij) = aire(ij+ip1jm) * pteta(ij+ip1jm,k)
       ENDDO
       tpn      = SSUM(iim,xpn,1)/ apoln
       tps      = SSUM(iim,xps,1)/ apols

       DO ij   = 1, iip1
         pteta(   ij   ,k)  = tpn
         pteta(ij+ip1jm,k)  = tps
       ENDDO
      ENDDO
!***********************
! Correction on teta due to surface pressure changes
      DO k = 1,llm
        DO j = 1,ip1jmp1
           pteta(j,k)= pteta(j,k)*(1+pdpfi(j)*pdt/pps(j))**kappa
        ENDDO
      ENDDO
!***********************

      DO k = 1,llm
         DO j = iip2,ip1jm
            pucov(j,k)= pucov(j,k) + pdufi(j,k) * pdt
         ENDDO
      ENDDO

      DO k = 1,llm
         DO j = 1,ip1jm
            pvcov(j,k)= pvcov(j,k) + pdvfi(j,k) * pdt
         ENDDO
      ENDDO

c
      DO j = 1,ip1jmp1
         pps(j) = pps(j) + pdpfi(j) * pdt
      ENDDO
 
      if (planet_type=="earth") then
      ! earth case, special treatment for first 2 tracers (water)
       DO iq = 1, 2
         DO k = 1,llm
            DO j = 1,ip1jmp1
               pq(j,k,iq)= pq(j,k,iq) + pdqfi(j,k,iq) * pdt
               pq(j,k,iq)= AMAX1( pq(j,k,iq), qtestw )
            ENDDO
         ENDDO
       ENDDO

       DO iq = 3, nqtot
         DO k = 1,llm
            DO j = 1,ip1jmp1
               pq(j,k,iq)= pq(j,k,iq) + pdqfi(j,k,iq) * pdt
               pq(j,k,iq)= AMAX1( pq(j,k,iq), qtestt )
            ENDDO
         ENDDO
       ENDDO
      else
      ! general case, treat all tracers equally)
       DO iq = 1, nqtot
         DO k = 1,llm
            DO j = 1,ip1jmp1
               pq(j,k,iq)= pq(j,k,iq) + pdqfi(j,k,iq) * pdt
               pq(j,k,iq)= AMAX1( pq(j,k,iq), qtestt )
            ENDDO
         ENDDO
       ENDDO
      endif ! of if (planet_type=="earth")

      DO  ij   = 1, iim
        xpn(ij) = aire(   ij   ) * pps(  ij     )
        xps(ij) = aire(ij+ip1jm) * pps(ij+ip1jm )
      ENDDO
      tpn      = SSUM(iim,xpn,1)/apoln
      tps      = SSUM(iim,xps,1)/apols

      DO ij   = 1, iip1
        pps (   ij     )  = tpn
        pps ( ij+ip1jm )  = tps
      ENDDO


      DO iq = 1, nqtot
        DO  k    = 1, llm
          DO  ij   = 1, iim
            xpn(ij) = aire(   ij   ) * pq(  ij    ,k,iq)
            xps(ij) = aire(ij+ip1jm) * pq(ij+ip1jm,k,iq)
          ENDDO
          tpn      = SSUM(iim,xpn,1)/apoln
          tps      = SSUM(iim,xps,1)/apols

          DO ij   = 1, iip1
            pq (   ij   ,k,iq)  = tpn
            pq (ij+ip1jm,k,iq)  = tps
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END
