!
!$Header$
!
MODULE mod_grid_phy_lmdz
  INTEGER,SAVE :: nbp_lon  ! == iim
  INTEGER,SAVE :: nbp_lat  ! == jjmp1
  INTEGER,SAVE :: nbp_lev  ! == llm
  INTEGER,SAVE :: klon_glo

  INTERFACE grid1dTo2d_glo
    MODULE PROCEDURE grid1dTo2d_glo_i,grid1dTo2d_glo_i1,grid1dTo2d_glo_i2,grid1dTo2d_glo_i3, &
                     grid1dTo2d_glo_r,grid1dTo2d_glo_r1,grid1dTo2d_glo_r2,grid1dTo2d_glo_r3, &
		     grid1dTo2d_glo_l,grid1dTo2d_glo_l1,grid1dTo2d_glo_l2,grid1dTo2d_glo_l3
   END INTERFACE 

   INTERFACE grid2dTo1d_glo
    MODULE PROCEDURE grid2dTo1d_glo_i,grid2dTo1d_glo_i1,grid2dTo1d_glo_i2,grid2dTo1d_glo_i3, &
                     grid2dTo1d_glo_r,grid2dTo1d_glo_r1,grid2dTo1d_glo_r2,grid2dTo1d_glo_r3, &
		     grid2dTo1d_glo_l,grid2dTo1d_glo_l1,grid2dTo1d_glo_l2,grid2dTo1d_glo_l3
   END INTERFACE 
 
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINE grid1dTo2d  !!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE Init_grid_phy_lmdz(iim,jjp1,llm)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: iim
  INTEGER, INTENT(in) :: jjp1
  INTEGER, INTENT(in) :: llm
  
    nbp_lon=iim
    nbp_lat=jjp1
    nbp_lev=llm
    klon_glo=(iim*jjp1)-2*(iim-1)
  
  END SUBROUTINE Init_grid_phy_lmdz
  
  
  SUBROUTINE grid1dTo2d_glo_i(VarIn,VarOut)  
  IMPLICIT NONE  
    INTEGER,INTENT(IN),DIMENSION(:)     :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:)  :: VarOut
    
    CALL grid1dTo2d_glo_igen(VarIn,VarOut,1)
  
  END SUBROUTINE grid1dTo2d_glo_i
  

  SUBROUTINE grid1dTo2d_glo_i1(VarIn,VarOut)  
  IMPLICIT NONE  
    INTEGER,INTENT(IN),DIMENSION(:,:)     :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:)  :: VarOut
    
    CALL grid1dTo2d_glo_igen(VarIn,VarOut,size(VarIn,2))
  
  END SUBROUTINE grid1dTo2d_glo_i1

  SUBROUTINE grid1dTo2d_glo_i2(VarIn,VarOut)  
  IMPLICIT NONE  
    INTEGER,INTENT(IN),DIMENSION(:,:,:)     :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:)  :: VarOut
    
    CALL grid1dTo2d_glo_igen(VarIn,VarOut,size(VarIn,2)*size(VarIn,3))
  
  END SUBROUTINE grid1dTo2d_glo_i2
  
  SUBROUTINE grid1dTo2d_glo_i3(VarIn,VarOut)  
  IMPLICIT NONE  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:)     :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:,:)  :: VarOut
    
    CALL grid1dTo2d_glo_igen(VarIn,VarOut,size(VarIn,2)*size(VarIn,3)*size(VarIn,4))
  
  END SUBROUTINE grid1dTo2d_glo_i3


  SUBROUTINE grid1dTo2d_glo_r(VarIn,VarOut)  
  IMPLICIT NONE  
    REAL,INTENT(IN),DIMENSION(:)     :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:)  :: VarOut
    
    CALL grid1dTo2d_glo_rgen(VarIn,VarOut,1)
  
  END SUBROUTINE grid1dTo2d_glo_r
  

  SUBROUTINE grid1dTo2d_glo_r1(VarIn,VarOut)  
  IMPLICIT NONE  
    REAL,INTENT(IN),DIMENSION(:,:)     :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:)  :: VarOut
    
    CALL grid1dTo2d_glo_rgen(VarIn,VarOut,size(VarIn,2))
  
  END SUBROUTINE grid1dTo2d_glo_r1

  SUBROUTINE grid1dTo2d_glo_r2(VarIn,VarOut)  
  IMPLICIT NONE  
    REAL,INTENT(IN),DIMENSION(:,:,:)     :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:)  :: VarOut
    
    CALL grid1dTo2d_glo_rgen(VarIn,VarOut,size(VarIn,2)*size(VarIn,3))
  
  END SUBROUTINE grid1dTo2d_glo_r2
  
  SUBROUTINE grid1dTo2d_glo_r3(VarIn,VarOut)  
  IMPLICIT NONE  
    REAL,INTENT(IN),DIMENSION(:,:,:,:)     :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:,:)  :: VarOut
    
    CALL grid1dTo2d_glo_rgen(VarIn,VarOut,size(VarIn,2)*size(VarIn,3)*size(VarIn,4))
  
  END SUBROUTINE grid1dTo2d_glo_r3
  
  
  
  SUBROUTINE grid1dTo2d_glo_l(VarIn,VarOut)  
  IMPLICIT NONE  
    LOGICAL,INTENT(IN),DIMENSION(:)     :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:)  :: VarOut
    
    CALL grid1dTo2d_glo_lgen(VarIn,VarOut,1)
  
  END SUBROUTINE grid1dTo2d_glo_l
  

  SUBROUTINE grid1dTo2d_glo_l1(VarIn,VarOut)  
  IMPLICIT NONE  
    LOGICAL,INTENT(IN),DIMENSION(:,:)     :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:)  :: VarOut
    
    CALL grid1dTo2d_glo_lgen(VarIn,VarOut,size(VarIn,2))
  
  END SUBROUTINE grid1dTo2d_glo_l1

  SUBROUTINE grid1dTo2d_glo_l2(VarIn,VarOut)  
  IMPLICIT NONE  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:)     :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:)  :: VarOut
    
    CALL grid1dTo2d_glo_lgen(VarIn,VarOut,size(VarIn,2)*size(VarIn,3))
  
  END SUBROUTINE grid1dTo2d_glo_l2
  
  SUBROUTINE grid1dTo2d_glo_l3(VarIn,VarOut)  
  IMPLICIT NONE  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:)     :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:,:)  :: VarOut
    
    CALL grid1dTo2d_glo_lgen(VarIn,VarOut,size(VarIn,2)*size(VarIn,3)*size(VarIn,4))
  
  END SUBROUTINE grid1dTo2d_glo_l3  
  
    SUBROUTINE grid2dTo1d_glo_i(VarIn,VarOut)  
  IMPLICIT NONE  
    INTEGER,INTENT(IN),DIMENSION(:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:)  :: VarOut
    
    CALL grid2dTo1d_glo_igen(VarIn,VarOut,1)
  
  END SUBROUTINE grid2dTo1d_glo_i
  

  SUBROUTINE grid2dTo1d_glo_i1(VarIn,VarOut)  
  IMPLICIT NONE  
    INTEGER,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:)  :: VarOut
    
    CALL grid2dTo1d_glo_igen(VarIn,VarOut,size(VarIn,3))
  
  END SUBROUTINE grid2dTo1d_glo_i1

  SUBROUTINE grid2dTo1d_glo_i2(VarIn,VarOut)  
  IMPLICIT NONE  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:)  :: VarOut
    
    CALL grid2dTo1d_glo_igen(VarIn,VarOut,size(VarIn,3)*size(VarIn,4))
  
  END SUBROUTINE grid2dTo1d_glo_i2
  
  SUBROUTINE grid2dTo1d_glo_i3(VarIn,VarOut)  
  IMPLICIT NONE  
    INTEGER,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(:,:,:,:)  :: VarOut
    
    CALL grid2dTo1d_glo_igen(VarIn,VarOut,size(VarIn,3)*size(VarIn,4)*size(VarIn,5))
  
  END SUBROUTINE grid2dTo1d_glo_i3
 



  SUBROUTINE grid2dTo1d_glo_r(VarIn,VarOut)  
  IMPLICIT NONE  
    REAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:)  :: VarOut
    
    CALL grid2dTo1d_glo_rgen(VarIn,VarOut,1)
  
  END SUBROUTINE grid2dTo1d_glo_r
  

  SUBROUTINE grid2dTo1d_glo_r1(VarIn,VarOut)  
  IMPLICIT NONE  
    REAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:)  :: VarOut
    
    CALL grid2dTo1d_glo_rgen(VarIn,VarOut,size(VarIn,3))
  
  END SUBROUTINE grid2dTo1d_glo_r1

  SUBROUTINE grid2dTo1d_glo_r2(VarIn,VarOut)  
  IMPLICIT NONE  
    REAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:)  :: VarOut
    
    CALL grid2dTo1d_glo_rgen(VarIn,VarOut,size(VarIn,3)*size(VarIn,4))
  
  END SUBROUTINE grid2dTo1d_glo_r2
  
  SUBROUTINE grid2dTo1d_glo_r3(VarIn,VarOut)  
  IMPLICIT NONE  
    REAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    REAL,INTENT(OUT),DIMENSION(:,:,:,:)  :: VarOut
    
    CALL grid2dTo1d_glo_rgen(VarIn,VarOut,size(VarIn,3)*size(VarIn,4)*size(VarIn,5))
  
  END SUBROUTINE grid2dTo1d_glo_r3



  SUBROUTINE grid2dTo1d_glo_l(VarIn,VarOut)  
  IMPLICIT NONE  
    LOGICAL,INTENT(IN),DIMENSION(:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:)  :: VarOut
    
    CALL grid2dTo1d_glo_lgen(VarIn,VarOut,1)
  
  END SUBROUTINE grid2dTo1d_glo_l
  

  SUBROUTINE grid2dTo1d_glo_l1(VarIn,VarOut)  
  IMPLICIT NONE  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:)  :: VarOut
    
    CALL grid2dTo1d_glo_lgen(VarIn,VarOut,size(VarIn,3))
  
  END SUBROUTINE grid2dTo1d_glo_l1

  SUBROUTINE grid2dTo1d_glo_l2(VarIn,VarOut)  
  IMPLICIT NONE  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:)  :: VarOut
    
    CALL grid2dTo1d_glo_lgen(VarIn,VarOut,size(VarIn,3)*size(VarIn,4))
  
  END SUBROUTINE grid2dTo1d_glo_l2
  
  SUBROUTINE grid2dTo1d_glo_l3(VarIn,VarOut)  
  IMPLICIT NONE  
    LOGICAL,INTENT(IN),DIMENSION(:,:,:,:,:) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(:,:,:,:)  :: VarOut
    
    CALL grid2dTo1d_glo_lgen(VarIn,VarOut,size(VarIn,3)*size(VarIn,4)*size(VarIn,5))
  
  END SUBROUTINE grid2dTo1d_glo_l3

END MODULE mod_grid_phy_lmdz


  
  SUBROUTINE grid1dTo2d_glo_igen(VarIn,VarOut,dimsize)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: dimsize
    INTEGER,INTENT(IN) ,DIMENSION(klon_glo,dimsize)       :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(nbp_lon*nbp_lat,dimsize)  :: VarOut
    INTEGER :: i,ij,Offset

    
    Offset=nbp_lon
        
    DO i=1,dimsize
      DO ij=1,klon_glo
        VarOut(ij+offset-1,i)=VarIn(ij,i)
      ENDDO
    ENDDO
    
    
    DO i=1,dimsize
      DO ij=1,nbp_lon
       VarOut(ij,i)=VarIn(1,i)
      ENDDO
    ENDDO
    
    
    DO i=1,dimsize
      DO ij=nbp_lon*(nbp_lat-1)+1,nbp_lat*nbp_lon
       VarOut(ij,i)=VarIn(klon_glo,i)
      ENDDO
    ENDDO

  END SUBROUTINE grid1dTo2d_glo_igen   


  SUBROUTINE grid1dTo2d_glo_rgen(VarIn,VarOut,dimsize)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: dimsize
    REAL,INTENT(IN) ,DIMENSION(klon_glo,dimsize)       :: VarIn
    REAL,INTENT(OUT),DIMENSION(nbp_lon*nbp_lat,dimsize)  :: VarOut
    INTEGER :: i,ij,Offset

   
    Offset=nbp_lon
        
    DO i=1,dimsize
      DO ij=1,klon_glo
        VarOut(ij+offset-1,i)=VarIn(ij,i)
      ENDDO
    ENDDO
    
    
    DO i=1,dimsize
      DO ij=1,nbp_lon
       VarOut(ij,i)=VarIn(1,i)
      ENDDO
    ENDDO
    
    
    DO i=1,dimsize
      DO ij=nbp_lon*(nbp_lat-1)+1,nbp_lat*nbp_lon
       VarOut(ij,i)=VarIn(klon_glo,i)
      ENDDO
    ENDDO

  END SUBROUTINE grid1dTo2d_glo_rgen   

  SUBROUTINE grid1dTo2d_glo_lgen(VarIn,VarOut,dimsize)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: dimsize
    LOGICAL,INTENT(IN) ,DIMENSION(klon_glo,dimsize)       :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(nbp_lon*nbp_lat,dimsize)  :: VarOut
    INTEGER :: i,ij,Offset

    Offset=nbp_lon
        
    DO i=1,dimsize
      DO ij=1,klon_glo
        VarOut(ij+offset-1,i)=VarIn(ij,i)
      ENDDO
    ENDDO
    
    
    DO i=1,dimsize
      DO ij=1,nbp_lon
       VarOut(ij,i)=VarIn(1,i)
      ENDDO
    ENDDO
    
    
    DO i=1,dimsize
      DO ij=nbp_lon*(nbp_lat-1)+1,nbp_lat*nbp_lon
       VarOut(ij,i)=VarIn(klon_glo,i)
      ENDDO
    ENDDO

  END SUBROUTINE grid1dTo2d_glo_lgen     
  
  
  SUBROUTINE grid2dTo1d_glo_igen(VarIn,VarOut,dimsize)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: dimsize
    INTEGER,INTENT(IN) ,DIMENSION(nbp_lon*nbp_lat,dimsize) :: VarIn
    INTEGER,INTENT(OUT),DIMENSION(klon_glo,dimsize)      :: VarOut
    INTEGER :: i,ij,offset

    offset=nbp_lon

    DO i=1,dimsize
      DO ij=1,klon_glo
        VarOut(ij,i)=VarIn(ij+offset-1,i)
      ENDDO
    ENDDO

    DO i=1,dimsize
      VarOut(1,i)=VarIn(1,i)
    ENDDO
    
  END SUBROUTINE grid2dTo1d_glo_igen   
  
  SUBROUTINE grid2dTo1d_glo_rgen(VarIn,VarOut,dimsize)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: dimsize
    REAL,INTENT(IN) ,DIMENSION(nbp_lon*nbp_lat,dimsize) :: VarIn
    REAL,INTENT(OUT),DIMENSION(klon_glo,dimsize)      :: VarOut
    INTEGER :: i,ij,offset

    offset=nbp_lon

    DO i=1,dimsize
      DO ij=1,klon_glo
        VarOut(ij,i)=VarIn(ij+offset-1,i)
      ENDDO
    ENDDO

    DO i=1,dimsize
      VarOut(1,i)=VarIn(1,i)
    ENDDO
    
  END SUBROUTINE grid2dTo1d_glo_rgen 
    
  SUBROUTINE grid2dTo1d_glo_lgen(VarIn,VarOut,dimsize)
    USE mod_grid_phy_lmdz
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: dimsize
    LOGICAL,INTENT(IN) ,DIMENSION(nbp_lon*nbp_lat,dimsize) :: VarIn
    LOGICAL,INTENT(OUT),DIMENSION(klon_glo,dimsize)      :: VarOut
    INTEGER :: i,ij,offset

    offset=nbp_lon

    DO i=1,dimsize
      DO ij=1,klon_glo
        VarOut(ij,i)=VarIn(ij+offset-1,i)
      ENDDO
    ENDDO

    DO i=1,dimsize
      VarOut(1,i)=VarIn(1,i)
    ENDDO
    
  END SUBROUTINE grid2dTo1d_glo_lgen   
