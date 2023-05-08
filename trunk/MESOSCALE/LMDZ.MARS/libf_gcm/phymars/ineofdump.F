      subroutine ineofdump

c._.  subroutine ineofdump(long, lati)

      implicit none
c
c     Initialise dumping of profiles for EOF calculations
c
#include "dimensions.h"
#include "dimphys.h"
#include "comvert.h"
#include "comcstfi.h"
#include "comgeomfi.h"
c
c._.#include "param2"
c._.#include "bats"
c._.#include "outcon"
c._.#include "phycmn"
c
#include "eofdump.h"
c
c._.      real*4 long(cols,rows)
c._.      real*4 lati(cols,rows)
      integer ig,i,j,l    
      logical firstcall
      save firstcall

      data firstcall/.true./

      if (firstcall) then
         firstcall=.false.
      endif

c
c     Set frequency for dumping at once per day
c
      ieofs=nint(daysec/dtphys)
      if (abs(float(ieofs)-daysec/dtphys).gt.1.e-8*daysec)
     :   stop'In ineofdump:  1 day .ne. n physics timesteps'
c
c     Header
c
      open(uehead,file='profiles.hdr',form='formatted')
c._.  write(uehead,*) rntape, kstart, ktotal, itspd, kphysic, refday
      write(uehead,*) 0.E+0,0,0,ieofs,1,0
c._.  write(uehead,*) cols, rows, npgrid, nl
      write(uehead,*) iim,npgrid/iim,npgrid,nlayermx

c._.  do i=1,cols,eofskip
c._.    do j=1+eofskip/2,rows,eofskip
c._.      write(uehead,*) long(i,j)*180./pi, lati(i,j)*180./pi
c._.    enddo
c._.  enddo
      do i=1,iim,eofskip
        do j=1+eofskip/2,jjm+1,eofskip    
          ig = 1+ (j-2)*iim +i
          if(j.eq.1) stop 'Problem in ineofdump.F'
          if(j.eq.jjm+1) stop 'Problem in ineofdump.F'
          write(uehead,*) long(ig)*180./pi, lati(ig)*180./pi
c         write(*,*) 'eof grid j=',j,' lat= ', lati(ig)*180./pi
        enddo
      enddo

      write(uehead,*) aps
      write(uehead,*) bps
      close(uehead)
c
c     Main profile file
c
      open(uedata,file='profiles.dat',form='unformatted')
      end
