      subroutine inigrads(if,im
     s  ,x,fx,xmin,xmax,jm,y,ymin,ymax,fy,lm,z,fz
     s  ,dt,file,titlel)


      implicit none

      integer if,im,jm,lm,i,j,l,lnblnk
      real x(im),y(jm),z(lm),fx,fy,fz,dt
      real xmin,xmax,ymin,ymax

      character file*10,titlel*40

#include "gradsdef.h"

      data unit/66,32,34,36,38,40,42,44,46,48/
      data nf/0/

      if (if.le.nf) stop'verifier les appels a inigrads'

      print*,'Entree dans inigrads'

      nf=if
      title(if)=titlel
      ivar(if)=0

      fichier(if)=file(1:lnblnk(file))

      firsttime(if)=.true.
      dtime(if)=dt

      iid(if)=1
      ifd(if)=im
      imd(if)=im
      do i=1,im
         xd(i,if)=x(i)*fx
         if(xd(i,if).lt.xmin) iid(if)=i+1
         if(xd(i,if).le.xmax) ifd(if)=i
      enddo
      print*,'On stoke du point ',iid(if),'  a ',ifd(if),' en x'

      jid(if)=1
      jfd(if)=jm
      jmd(if)=jm
      do j=1,jm
         yd(j,if)=y(j)*fy
         if(yd(j,if).gt.ymax) jid(if)=j+1
         if(yd(j,if).ge.ymin) jfd(if)=j
      enddo
      print*,'On stoke du point ',jid(if),'  a ',jfd(if),' en y'

      print*,'Open de dat'
      print*,'file=',file
      print*,'fichier(if)=',fichier(if)

      print*,4*(ifd(if)-iid(if))*(jfd(if)-jid(if))
      print*,file(1:lnblnk(file))//'.dat'

      OPEN (unit(if)+1,FILE=file(1:lnblnk(file))//'.dat',
     s   FORM='UNFORMATTED',
     s   ACCESS='DIRECT'
     s  ,RECL=4*(ifd(if)-iid(if)+1)*(jfd(if)-jid(if)+1))

      print*,'Open de dat ok'

      lmd(if)=lm
      do l=1,lm
         zd(l,if)=z(l)*fz
      enddo

      irec(if)=0

      print*,if,imd(if),jmd(if),lmd(if)
      print*,'if,imd(if),jmd(if),lmd(if)'

      return
      end
