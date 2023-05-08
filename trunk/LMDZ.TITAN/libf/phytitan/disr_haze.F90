subroutine disr_haze(dz,press,wno,taeros,ssa,cbar)

use datafile_mod, only: datadir

IMPLICIT NONE

! ==========================================================================
!
!  Purpose :
!    Interpolate values of extinction coefficient, single scattering albedo
!    and asymetry factor from hazetable ( Lavvas et al. 2010, mean profile, no
!    detached layer)
!
!    + JVO 2017 : Vertical extension out of table implemented
!
!  Author :
!     Jan Vatant d'Ollone (2016)
!
! ==========================================================================

real*8,intent(in)   :: dz, press, wno
real*8,intent(inout):: taeros, ssa, cbar

!---------------------------
! NB !!
! taeros is the integrated extinction over the layer
! (extinction * thickness of layer)
!---------------------------

integer             :: i, j, iw, ip, ierr
real*8              :: wln, factw, factp
real*8              :: tmp_p, fact_t
integer,parameter   :: nbwl_PL=328, nblev_PL=162
real*8,save         :: ext_PL(nblev_PL,nbwl_PL), ssa_PL(nblev_PL,nbwl_PL), asf_PL(nblev_PL,nbwl_PL)
real*8,save         :: wl_PL(nbwl_PL), press_PL(nblev_PL)
logical,save        :: firstcall=.true.
character(len=15)   :: dummy

!$OMP THREADPRIVATE(ext_PL,ssa_PL,asf_PL,wl_PL,press_PL,firstcall)

if (firstcall) then
   print*,"We use DISR haze mean profile from P.Lavvas"

! read PL table   
! wl_PL in nm
! press_PL in Pa
   open(11,file=TRIM(datadir)//'/hazetable_PL_original.dat',status="old",iostat=ierr)
   read(11,*) dummy,wl_PL
   do i=1,nblev_PL
     read(11,*) press_PL(i),ext_PL(i,:)  ! in cm-1
   enddo
   do i=1,nblev_PL
     read(11,*) press_PL(i),ssa_PL(i,:)
   enddo
   do i=1,nblev_PL
     read(11,*) press_PL(i),asf_PL(i,:)
   enddo
   close(11)
   ! convert press_PL into millibar for comparison to press in the generic
   press_PL(:)=press_PL(:)*1E-2
   
   firstcall=.false.
endif

! convert wno (in cm-1) into wln (nm)
wln=1E7/wno

! interpolate the needed values from the table

iw=1
do i=2,nbwl_PL
  if(wln.gt.wl_PL(i)) then
    iw=i
  endif
enddo

ip=1
do j=2,nblev_PL
  if(press.lt.press_PL(j)) then
    ip=j
  endif
enddo

!----------------- Interpolate values from the hazetable  --------------------
if (iw.ne. nbwl_PL) then
  factw = (wln-wl_PL(iw)) / (wl_PL(iw+1)-wl_PL(iw))
endif

if (ip .ne. nblev_PL) then
  factp = (press-press_PL(ip)) / (press_PL(ip+1)-press_PL(ip))
endif

! Lin-Log interpolation : linear on wln, logarithmic on press

if((ip.ne.nblev_PL) .and. (iw.ne.nbwl_PL)) then 

taeros = ( ext_PL(ip,iw)*(1.-factw)   + ext_PL(ip,iw+1)  *factw ) ** (1.-factp) &
        *( ext_PL(ip+1,iw)*(1.-factw) + ext_PL(ip+1,iw+1)*factw ) ** factp

ssa    = ( ssa_PL(ip,iw)*(1.-factw)   + ssa_PL(ip,iw+1)  *factw ) ** (1.-factp) &
        *( ssa_PL(ip+1,iw)*(1.-factw) + ssa_PL(ip+1,iw+1)*factw ) ** factp

cbar   = ( asf_PL(ip,iw)*(1.-factw)   + asf_PL(ip,iw+1)  *factw ) ** (1.-factp) &
        *( asf_PL(ip+1,iw)*(1.-factw) + asf_PL(ip+1,iw+1)*factw ) ** factp

else if ((ip.ne.nblev_PL) .and. (iw.eq.nbwl_PL)) then

taeros =  ext_PL(ip,iw)**(1.-factp) * ext_PL(ip+1,iw)**factp
ssa    =  ssa_PL(ip,iw)**(1.-factp) * ssa_PL(ip+1,iw)**factp
cbar   =  asf_PL(ip,iw)**(1.-factp) * asf_PL(ip+1,iw)**factp


! In case of vertical extension over the max of the table
! We take the scale height on the last 5 levels (more it's not quite log)
! Arbitray threshold pressure value, just to deal with the last level press=0
! We do not touch to ssa and cbar and let them at the value of last level
! (extrap would lead to too dark aerosols)

else if(ip.eq.nblev_PL) then 

  tmp_p = press

  if ( tmp_p .lt. 1.E-15 ) then 
    tmp_p = 1.E-15
  endif

  if(iw.ne.nbwl_PL) then  

     fact_t = log10( ( ext_PL(ip,iw)*(1.-factw) + ext_PL(ip,iw+1)  *factw )   & 
           / ( ext_PL(ip-5,iw)*(1.-factw) + ext_PL(ip-5,iw+1)  *factw ) )
           
     taeros =  ext_PL(ip,iw)*(1.-factw) + ext_PL(ip,iw+1)*factw
     ssa    =  ssa_PL(ip,iw)*(1.-factw) + ssa_PL(ip,iw+1)*factw
     cbar   =  asf_PL(ip,iw)*(1.-factw) + asf_PL(ip,iw+1)*factw

  else if (iw.eq.nbwl_PL) then

     fact_t = log10( ext_PL(ip,iw) / ext_PL(ip-5,iw) )
    
     taeros =  ext_PL(ip,iw)
     ssa    =  ssa_PL(ip,iw)
     cbar   =  asf_PL(ip,iw)

  endif

  fact_t = fact_t / log10( press_PL(ip) / press_PL(ip-5) )

  taeros = taeros * ( tmp_p / press_PL(ip) ) ** fact_t

endif

taeros=taeros*dz*1.E2 ! ext in cm-1 * thickness in m * 1E2

end subroutine disr_haze
