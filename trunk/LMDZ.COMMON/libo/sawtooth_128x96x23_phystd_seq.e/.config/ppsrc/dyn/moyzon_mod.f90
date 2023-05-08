










      MODULE moyzon_mod
! Moyennes zonales pour transmission a la physique 
!======================================================================
! Specifique a Titan
!
!======================================================================
! Declaration des variables

      REAL,ALLOCATABLE,SAVE :: zplevbar(:,:),zplaybar(:,:)
      REAL,ALLOCATABLE,SAVE :: ztfibar(:,:),zqfibar(:,:,:)
      REAL,ALLOCATABLE,SAVE :: zphibar(:,:),zphisbar(:)
      REAL,ALLOCATABLE,SAVE :: zpkbar(:,:),ztetabar(:,:)
      REAL,ALLOCATABLE,SAVE :: zzlevbar(:,:),zzlaybar(:,:)

! pmoy: global averaged pressure...
! tmoy: global averaged temperature...
! put here to be transfered to Titan routines...
! to be changed...
      REAL,ALLOCATABLE,SAVE :: plevmoy(:),playmoy(:)
      REAL,ALLOCATABLE,SAVE :: zlevmoy(:),zlaymoy(:),phimoy(:)
      REAL,ALLOCATABLE,SAVE :: tmoy(:),tetamoy(:),pkmoy(:)
      INTEGER,ALLOCATABLE,SAVE :: klat(:)

CONTAINS
!======================================================================
SUBROUTINE moyzon_init(klon,llm,nqtot)
!! This routine needs physics
IMPLICIT NONE
       INTEGER :: klon,llm,nqtot

      ALLOCATE(zplevbar(klon,llm+1),zplaybar(klon,llm))
      ALLOCATE(zphibar(klon,llm),zphisbar(klon))
      ALLOCATE(ztfibar(klon,llm),zqfibar(klon,llm,nqtot))
      ALLOCATE(zpkbar(klon,llm),ztetabar(klon,llm))
      ALLOCATE(zzlevbar(klon,llm+1),zzlaybar(klon,llm))
END SUBROUTINE moyzon_init

!======================================================================
SUBROUTINE moyzon(nlev,var,varbar)

IMPLICIT NONE
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

      INTEGER :: nlev
      REAL,dimension(iip1,nlev) :: var
      REAL,dimension(nlev)      :: varbar

      INTEGER :: i

      varbar(:) = 0.
      do i=1,iim
        varbar(:)=varbar(:)+var(i,:)/iim
      enddo

      return
END SUBROUTINE moyzon

!======================================================================
      END MODULE moyzon_mod
