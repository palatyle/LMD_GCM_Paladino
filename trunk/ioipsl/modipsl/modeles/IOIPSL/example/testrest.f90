PROGRAM testrest
!-
!$Id: testrest.f90 386 2008-09-04 08:38:48Z bellier $
!-
! This software is governed by the CeCILL license
! See IOIPSL/IOIPSL_License_CeCILL.txt
!---------------------------------------------------------------------
!- This program provide a an example of the basic usage of REST.
!- Here the test the time sampling and averaging. Thus a long
!- time-series is produced and sampled in different ways.
!---------------------------------------------------------------------
  USE ioipsl
!
  IMPLICIT NONE
!
  INTEGER :: iim,jjm,llm
  PARAMETER (iim=12,jjm=10,llm=2)
!
  REAL :: champ1(iim,jjm,llm),champ2(iim,jjm,llm+1),champ3(iim,jjm,llm)
  REAL :: champ4(iim,jjm)
  REAL :: champ_read(iim,jjm,llm)
  REAL :: lon(iim,jjm),lat(iim,jjm),lev(llm)
!
  INTEGER :: i,j,l,fid
  INTEGER :: day=1,month=1,year=1997,itau=1
!
  REAL :: julday,un_mois,un_an
  REAL :: deltat=86400,dt_wrt,dt_op,dt_wrt2,dt_op2
  CHARACTER*20 :: fnamein,fnameout
!
  REAL :: pi=3.1415
!---------------------------------------------------------------------
!-
! 0.0 Choose a  gregorian calendar
!-
  CALL ioconf_calendar ('gregorian')
!-
! 1.0 Define a few variables we will need.
!     These are the coordinates the file name and the date.
!-
  DO i=1,iim
    DO j=1,jjm
      lon(i,j) = &
 &     ((float(iim/2)+0.5)-float(i))*pi/float(iim/2)*(-1.)*180./pi
      lat(i,j) = &
 &     (180./pi)*ASIN(((float(jjm/2)+0.5)-float(j))/float(jjm/2))
    ENDDO
  ENDDO
!-
  DO l=1,llm
    lev(l) = float(l)/llm
  ENDDO
!-
! 1.1 The chosen date is 15 Feb. 1997 as stated above. It has to be 
!     transformed into julian days using the calendar provided by 
!     IOIPSL.
!-
  CALL ymds2ju (year,month,day,0.,julday)
  CALL ioget_calendar (un_an)
  un_mois = un_an/12.
  dt_wrt = un_mois*deltat
  dt_op = deltat
  dt_wrt2 = -1.
  dt_op2 = deltat
!-
  fnamein = 'NONE'
  fnameout = 'restfile'
!-
! 2.0 Create a restart file from nothing !
!-
  CALL restini (fnamein,iim,jjm,lon,lat,llm,lev,fnameout, &
 &              itau,julday,deltat,fid)
!-
  champ1(:,:,:) = ASIN(1.0)
  champ2(:,:,:) = EXP(ASIN(1.0))
!-
  CALL ioconf_setatt ('units','1')
  CALL ioconf_setatt ('long_name','Tests 1 for a real variable')
  CALL restput (fid,'test1',iim,jjm,llm,itau,champ1)
!-
  CALL ioconf_setatt ('units','1')
  CALL ioconf_setatt ('long_name','Tests 2 for a real variable')
  CALL restput (fid,'test2',iim,jjm,llm+1,itau,champ2)
!-
  CALL restclo ()
!-
  WRITE(*,*) '============== FIRST FILE CLOSED =============='
!-
!  3.0 Reopen the restart file and check that the values read are correct
!-
  fnamein = 'restfile'
  fnameout = 'restfilebis'
!-
  CALL restini (fnamein,iim,jjm,lon,lat,llm,lev,fnameout, &
 &              itau,julday,deltat,fid)
!-
  CALL restget (fid,'test1',iim,jjm,llm,itau,.FALSE.,champ_read)
!-
  itau = itau+10
  CALL restput (fid,'test1',iim,jjm,llm,itau,champ_read)
  CALL restput (fid,'test2',iim,jjm,llm+1,itau,champ2)
!-
  itau = itau+10
  champ3(:,:,:) = champ_read(:,:,:)+champ2(:,:,1:llm) 
  CALL restput (fid,'test1',iim,jjm,llm,itau,champ3)
!-
  CALL restclo ()
!-
  WRITE(*,'(a25,e36.30)') 'The input data : ',champ1(1,1,1)
  WRITE(*,'(a25,e36.30)') 'The restart data : ',champ_read(1,1,1)
!-
!  4.0 Reopen the restart file and add another time step
!-
  fnamein = 'restfilebis'
  fnameout = 'restfilebis'
!-
  CALL restini (fnamein,iim,jjm,lon,lat,llm,lev,fnameout, &
 &              itau,julday,deltat,fid)
!-
  itau = itau+10
  CALL restput (fid,'test1',iim,jjm,llm,itau,champ1)
  CALL ioconf_setatt ('units','1')
  CALL ioconf_setatt ('long_name', &
 &                    'Test a variable with another dimension')
  CALL restput (fid,'test4',iim,jjm,0,itau,champ4)
!-
  CALL restclo ()
!-------------------
END PROGRAM testrest
