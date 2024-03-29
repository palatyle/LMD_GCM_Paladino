PROGRAM testflio
!-
!$Id: testflio.f90 887 2010-02-08 09:48:39Z bellier $
!-
! This software is governed by the CeCILL license
! See IOIPSL/IOIPSL_License_CeCILL.txt
!---------------------------------------------------------------------
!- This program is an example of how to use "fliocom".
!---------------------------------------------------------------------
  USE ioipsl
  USE defprec
!-
  IMPLICIT NONE
!-
  INTEGER,PARAMETER :: iimf=3, jjmf=2, llmf=3, itmf=5
!-
  REAL,DIMENSION(iimf,jjmf,llmf,itmf) :: forcing
  REAL,DIMENSION(iimf,jjmf)           :: lonf,latf
  REAL,DIMENSION(llmf)                :: levf
  INTEGER,DIMENSION(itmf)             :: itauf
  REAL :: parf=7.
!-
! REAL,DIMENSION(iimf,jjmf,llmf,itmf) :: forcingw
  REAL,DIMENSION(iimf)                :: lonw
  REAL,DIMENSION(jjmf)                :: latw
  REAL,DIMENSION(llmf)                :: levw
  INTEGER,DIMENSION(itmf)             :: itauw
! REAL :: parw
!-
  CHARACTER(LEN=10),DIMENSION(flio_max_dims) :: f_n_d
  INTEGER,DIMENSION(flio_max_dims) :: f_l_d,f_i_d
  INTEGER :: i_t_w,n_d_w,n_a_w
  INTEGER,DIMENSION(flio_max_var_dims) :: l_d_w,i_d_w
  LOGICAL :: exv,exa
  INTEGER :: iim,jjm,llm,ttm
  INTEGER :: i,j,l,it,fid,id_dom
  CHARACTER(LEN=25) :: f_n,v_n,c_n,c_tmp=' '
  INTEGER :: year,month,day
  REAL :: date,datew,dt,dtw,sec
  REAL,DIMENSION(10) :: r_w
  CHARACTER(LEN=30),DIMENSION(:),ALLOCATABLE :: cn_d,cn_v,cn_a
  CHARACTER(LEN=30) :: cn_u,f_n_c
!---------------------------------------------------------------------
!-
! Create a file for the lateral forcing
!-
  DO i=1,iimf
    DO j=1,jjmf
      lonf(i,j) = (REAL(i-1)/REAL(MAX(iimf-1,1)))*180.
      latf(i,j) = (REAL(j-1)/REAL(MAX(jjmf-1,1)))*90.
    ENDDO
  ENDDO
!-
  levf(1:3) = (/ 0.7, 0.5, 0.3 /)
!-
  DO it = 1,itmf
    itauf(it) = it+10
    DO i=1,iimf
      DO j=1,jjmf
        DO l=1,llmf
          forcing(i,j,l,it) = 10.*(10.*((10.*it)+l)+j)+i
        ENDDO
      ENDDO
    ENDDO
  ENDDO
!-
  WRITE (*,'(" -----------------------------")')
  WRITE (*,'(" ------- Using fliocom -------")')
  WRITE (*,'(" -----------------------------")')
  WRITE(*,'(/," ",A,3(1X,I5))') &
 &  'flio_max_files,flio_max_dims,flio_max_var_dims : ', &
 &   flio_max_files,flio_max_dims,flio_max_var_dims
!-
  year = 1997; month = 5; day = 10; sec = 0.;
  CALL ymds2ju (year,month,day,sec,date)
  dt = 86400./2.
  WRITE (*,'(/," ----------------------------------")')
  WRITE (*,'(/," t_init : ",F11.3)') date
!-
! Create a file and add variables and attributes
!-
  WRITE (*,'(/," ----------------------------------")')
!-
  f_n = 'testflio.nc'
  f_l_d(1:5) = (/  iimf,  jjmf,  llmf,    -1,   10   /)
  f_n_d(1:5) = (/ "lon ","lat ","lev ","time","ud1 " /)
!-
  CALL flio_dom_set &
 & (1,1,(/1,2/),(/iimf,jjmf/),(/iimf,jjmf/), &
 &  (/1,1/),(/iimf,jjmf/),(/0,0/),(/0,0/),"orange",id_dom)
  WRITE (*,'(/," --> fliocrfd (",A,",...)")') TRIM(f_n)
  CALL fliocrfd (TRIM(f_n),f_n_d(1:5),f_l_d(1:5),fid, &
 &               id_dom=id_dom,c_f_n=f_n_c)
  WRITE (*,'(/,"   created file name : ",A)') TRIM(f_n_c)
  CALL flio_dom_unset (id_dom)
!-
  WRITE (*,'(/," --> fliopstc")')
  CALL fliopstc (fid, &
 &  x_axis_2d=lonf(:,:),y_axis_2d=latf(:,:),z_axis=levf(:), &
 &  t_axis=itauf(:),t_init=date,t_step=dt,t_calendar="gregorian")
!-
  CALL fliodefv (fid,'forcing',(/ 1,2,3,4 /), &
 &               v_t=flio_r4,units='1',long_name='Forcing field')
  CALL flioputa (fid,'forcing','vecteur_3',(/3.,5.,7./))
  CALL flioputa (fid,'forcing','text_n',"Text attribute")
  CALL fliodefv (fid,'my_var_1',(/ 5 /), &
 &               v_t=flio_r4,units='1',long_name='my_var_1', &
 &               valid_min=-10.,valid_max=+20.,fillvalue=+50.)
  CALL fliodefv (fid,'Var_vr4', &
 &               v_t=flio_r4,units='1',long_name='Var_vr4')
  CALL flioputa (fid,'Var_vr4','att_1',735)
  CALL flioputa (fid,'?','Param_a4',REAL(parf,KIND=4))
  CALL fliodefv (fid,'Var_vr8', &
 &               v_t=flio_r8,units='1',long_name='Var_vr8')
  CALL fliocpya (fid,'Var_vr4','att_1',fid,'Var_vr8')
  CALL flioputa (fid,'?','Param_a8',REAL(parf,KIND=8))
  CALL fliodefv (fid,'Var_vi2', &
 &               v_t=flio_i2,units='1',long_name='Var_vi2')
  CALL flioputv (fid,'forcing',forcing)
  CALL flioputv (fid,'my_var_1',(/3.,4.,5./),start=(/3/))
  CALL flioputv (fid,'Var_vr4',parf)
  CALL flioputv (fid,'Var_vr8',parf)
  CALL flioputv (fid,'Var_vi2',INT(parf,KIND=i_2))
!-
  WRITE (*,'(/," --> flioclo")')
  CALL flioclo (fid)
!-
! Inspect the file
!-
  WRITE (*,'(/," ----------------------------------")')
!-
  WRITE (*,'(/," --> fliodmpf")')
  CALL fliodmpf (TRIM(f_n_c))
!-
! Open the file and obtain information about the file
!-
  WRITE (*,'(/," ----------------------------------")')
!-
  WRITE (*,'(/," --> flioopfd")')
  i = 0; j = 0; l = 0;
  CALL flioopfd (TRIM(f_n_c),fid,nb_dim=i,nb_var=j,nb_gat=l)
  WRITE (*,'(" Number of dimensions        in the file : ",I2)') i
  WRITE (*,'(" Number of variables         in the file : ",I2)') j
  WRITE (*,'(" Number of global attributes in the file : ",I2)') l
!-
  WRITE (*,'(/," --> flioinqf")')
  CALL flioinqf &
 & (fid,id_dim=f_i_d,ln_dim=f_l_d)
  IF (i > 0) THEN
    WRITE (*,'(" Identifiers of the dimensions :",/,(5(1X,I7),:))') &
 &    f_i_d(1:MIN(i,SIZE(f_i_d)))
    WRITE (*,'(" Dimensions :",/,(5(1X,I7),:))') &
 &    f_l_d(1:MIN(i,SIZE(f_l_d)))
  ENDIF
!-
  WRITE (*,'(/," --> flioinqn")')
  ALLOCATE(cn_d(i),cn_v(j),cn_a(l))
  CALL flioinqn (fid,cn_dim=cn_d,cn_var=cn_v,cn_gat=cn_a,cn_uld=cn_u)
  WRITE (*,'(" Names of the dimensions in the file : ")')
  DO it=1,i
    WRITE (*,'("   """,A,"""")') TRIM(cn_d(it))
  ENDDO
  WRITE (*,'(" Names of the variables in the file : ")')
  DO it=1,j
    WRITE (*,'("   """,A,"""")') TRIM(cn_v(it))
  ENDDO
  WRITE (*,'(" Names of the global attributes in the file : ")')
  DO it=1,l
    WRITE (*,'("   """,A,"""")') TRIM(cn_a(it))
  ENDDO
  WRITE (*,'(" Name of the unlimited dimension : ")')
  WRITE (*,'("   """,A,"""")') TRIM(cn_u)
  DEALLOCATE(cn_d,cn_v,cn_a)
!-
  WRITE (*,'(/," --> flioqstc")')
  c_n = "x"
  CALL flioqstc (fid,TRIM(c_n),exv,v_n)
  IF (exv) THEN
    WRITE (*,'(" Name of the """,A,""" coordinate : ",A)') &
 &   TRIM(c_n),TRIM(v_n)
    CALL flioinqv &
 &   (fid,TRIM(v_n),exv,nb_dims=n_d_w,len_dims=l_d_w,nb_atts=n_a_w) 
    IF (n_d_w > 0) THEN
      WRITE (*,'(" Number of dimensions : ",I2)') n_d_w
      WRITE (*,'(" Dimensions :",/,(5(1X,I7),:))') &
 &      l_d_w(1:MIN(n_d_w,SIZE(l_d_w)))
      iim=l_d_w(1)
    ENDIF
    IF (n_a_w > 0) THEN
      WRITE (*,'(" Number of attributes : ",I2)') n_a_w
    ENDIF
  ELSE
    WRITE (*,'(" Coordinate """,A,""" not found")') TRIM(c_n)
  ENDIF
!-
  WRITE (*,'(/," --> flioqstc")')
  c_n = "y"
  CALL flioqstc (fid,TRIM(c_n),exv,v_n)
  IF (exv) THEN
    WRITE (*,'(" Name of the """,A,""" coordinate : ",A)') &
 &   TRIM(c_n),TRIM(v_n)
    CALL flioinqv &
 &   (fid,TRIM(v_n),exv,nb_dims=n_d_w,len_dims=l_d_w,nb_atts=n_a_w) 
    IF (n_d_w > 0) THEN
      WRITE (*,'(" Number of dimensions : ",I2)') n_d_w
      WRITE (*,'(" Dimensions :",/,(5(1X,I7),:))') &
 &      l_d_w(1:MIN(n_d_w,SIZE(l_d_w)))
      jjm=l_d_w(n_d_w)
    ENDIF
    IF (n_a_w > 0) THEN
      WRITE (*,'(" Number of attributes : ",I2)') n_a_w
    ENDIF
  ELSE
    WRITE (*,'(" Coordinate """,A,""" not found")') TRIM(c_n)
  ENDIF
!-
  WRITE (*,'(/," --> flioqstc")')
  c_n = "z"
  CALL flioqstc (fid,TRIM(c_n),exv,v_n)
  IF (exv) THEN
    WRITE (*,'(" Name of the """,A,""" coordinate : ",A)') &
 &   TRIM(c_n),TRIM(v_n)
    CALL flioinqv &
 &   (fid,TRIM(v_n),exv,nb_dims=n_d_w,len_dims=l_d_w,nb_atts=n_a_w) 
    IF (n_d_w > 0) THEN
      WRITE (*,'(" Number of dimensions : ",I2)') n_d_w
      WRITE (*,'(" Dimensions :",/,(5(1X,I7),:))') &
 &      l_d_w(1:MIN(n_d_w,SIZE(l_d_w)))
      llm=l_d_w(1)
    ENDIF
    IF (n_a_w > 0) THEN
      WRITE (*,'(" Number of attributes : ",I2)') n_a_w
    ENDIF
  ELSE
    WRITE (*,'(" Coordinate """,A,""" not found")') TRIM(c_n)
  ENDIF
!-
  WRITE (*,'(/," --> flioqstc")')
  c_n = "t"
  CALL flioqstc (fid,TRIM(c_n),exv,v_n)
  IF (exv) THEN
    WRITE (*,'(" Name of the """,A,""" coordinate : ",A)') &
 &   TRIM(c_n),TRIM(v_n)
    CALL flioinqv &
 &   (fid,TRIM(v_n),exv,nb_dims=n_d_w,len_dims=l_d_w,nb_atts=n_a_w) 
    IF (n_d_w > 0) THEN
      WRITE (*,'(" Number of dimensions : ",I2)') n_d_w
      WRITE (*,'(" Dimensions :",/,(5(1X,I7),:))') &
 &      l_d_w(1:MIN(n_d_w,SIZE(l_d_w)))
      ttm=l_d_w(1)
    ENDIF
    IF (n_a_w > 0) THEN
      WRITE (*,'(" Number of attributes : ",I2)') n_a_w
    ENDIF
  ELSE
    WRITE (*,'(" Coordinate """,A,""" not found")') TRIM(c_n)
  ENDIF
!-
  WRITE (*,'(/," --> fliogstc")')
  CALL fliogstc (fid, &
 &  x_axis=lonw,y_axis=latw,z_axis=levw, &
 &  t_axis=itauw(:),t_init=datew,t_step=dtw,t_calendar=c_tmp, &
 &  t_count=3)
  WRITE (*,'(" x_axis :",/,(5(1X,1PE11.3),:))') lonw(1:iim)
  WRITE (*,'(" y_axis :",/,(5(1X,1PE11.3),:))') latw(1:jjm)
  WRITE (*,'(" z_axis :",/,(5(1X,1PE11.3),:))') levw(1:llm)
  WRITE (*,'(" t_axis(1:3) :",/,(5(1X,I5),:))') itauw(1:3)
  WRITE (*,'(" t_calendar  : ",A)') '"'//TRIM(c_tmp)//'"'
  WRITE (*,'(" t_init,t_step :",(2(1X,1PE15.7),:))') datew,dtw
!-
  v_n = 'forcing'
  WRITE (*,'(/," Variable : """,A,"""")') TRIM(v_n)
  WRITE (*,'(" --> flioinqv(...,""",A,""",...)")') TRIM(v_n)
  CALL flioinqv &
 & (fid,TRIM(v_n),exv,v_t=i_t_w, &
 &  nb_dims=n_d_w,len_dims=l_d_w,id_dims=i_d_w,nb_atts=n_a_w)
  IF (exv) THEN
    WRITE (*,'(" External type        : ",I2)') i_t_w
    WRITE (*,'(" Number of dimensions : ",I2)') n_d_w
    WRITE (*,'(" Dimensions :",/,5(1X,I7,:))') l_d_w(1:n_d_w)
    WRITE (*,'(" Identifiers :",/,5(1X,I7,:))') i_d_w(1:n_d_w)
    WRITE (*,'(" Number of attributes : ",I2)') n_a_w
!---
    IF (n_a_w > 0) THEN
      ALLOCATE(cn_a(5))
      CALL flioinqv (fid,TRIM(v_n),exv,cn_atts=cn_a)
      WRITE (*,'(" Names of the attributes : ")')
      DO it=1,n_a_w
        WRITE (*,'("   """,A,"""")') TRIM(cn_a(it))
      ENDDO
      DEALLOCATE(cn_a)
    ENDIF
!---
    WRITE (*, &
 &   '(" --> flioinqa(...,""",A,""",""text_n"",...)")') TRIM(v_n)
    CALL flioinqa (fid,TRIM(v_n),"text_n",exa,a_l=l)
    IF (exa) THEN
      WRITE (*,'("  len(text_n) : ",I3)') l
!-----
      WRITE (*,'(" --> fliogeta(...,""",A,""",""text_n"",...)")') &
 &     TRIM(v_n)
      CALL fliogeta (fid,TRIM(v_n),"text_n",c_tmp)
      WRITE (*,'("  ""text_n"" : """,A,"""")') TRIM(c_tmp)
    ELSE
      WRITE (*,'(" Attribute not found")')
    ENDIF
  ELSE
    WRITE (*,'(" Variable not found")')
  ENDIF
!-
  v_n = 'my_var_1'
  WRITE (*,'(/," Variable : """,A,"""")') TRIM(v_n)
  WRITE (*,'(" --> flioinqv(...,""",A,""",...)")') TRIM(v_n)
  CALL flioinqv &
 & (fid,TRIM(v_n),exv, &
 &  nb_dims=n_d_w,len_dims=l_d_w,id_dims=i_d_w,nb_atts=n_a_w)
  IF (exv) THEN
    WRITE (*,'(" Number of dimensions : ",I2)') n_d_w
    WRITE (*,'(" Dimensions :",/,5(1X,I7,:))') l_d_w(1:n_d_w)
    WRITE (*,'(" Identifiers :",/,5(1X,I7,:))') i_d_w(1:n_d_w)
    WRITE (*,'(" Number of attributes : ",I2)') n_a_w
    WRITE (*,'(" --> fliogetv(...,""",A,""",...)")') TRIM(v_n)
    CALL fliogetv (fid,TRIM(v_n),r_w)
    WRITE (*, &
 &   '("  Values :",/,(5(1X,1PE11.3),:))') r_w
    i = 2; j = 5;
    WRITE (*,'(" --> fliogetv(...,",A, &
 &             ",start=(/",I2,"/),count=(/",I2,"/))")') &
 &   '"'//TRIM(v_n)//'"',i,j
    CALL fliogetv (fid,TRIM(v_n),r_w,start=(/i/),count=(/j/))
    WRITE (*, &
 &   '("  Values(",I2,":",I2,") :",/,(5(1X,1PE11.3),:))') &
 &   i,i+j-1,r_w(1:j)
  ENDIF
!-
  WRITE (*,'(/," --> flioclo")')
  CALL flioclo (fid)
!-
  WRITE (*,'(/," ----------------------------------")')
!-------------------
END PROGRAM testflio
