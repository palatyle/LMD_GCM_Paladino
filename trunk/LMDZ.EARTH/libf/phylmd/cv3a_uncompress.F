      SUBROUTINE cv3a_uncompress(nloc,len,ncum,nd,ntra,idcum
     :         ,iflag,kbas,ktop
     :         ,precip,cbmf,sig,w0,ptop2
     :         ,ft,fq,fu,fv,ftra
     :         ,Ma,mip,Vprecip,upwd,dnwd,dnwd0
     :         ,qcondc,wd,cape,cin
     :         ,tvp
     :         ,ftd,fqd
     :         ,Plim1,Plim2,asupmax,supmax0
     :         ,asupmaxmin
     o         ,iflag1,kbas1,ktop1
     :         ,precip1,cbmf1,sig1,w01,ptop21
     :         ,ft1,fq1,fu1,fv1,ftra1
     :         ,Ma1,mip1,Vprecip1,upwd1,dnwd1,dnwd01
     :         ,qcondc1,wd1,cape1,cin1
     :         ,tvp1
     :         ,ftd1,fqd1
     :         ,Plim11,Plim21,asupmax1,supmax01
     :         ,asupmaxmin1     )
***************************************************************
*                                                             *
* CV3A_UNCOMPRESS                                             *
*                                                             *
*                                                             *
* written by   : Sandrine Bony-Lena , 17/05/2003, 11.22.15    *
* modified by  : Jean-Yves Grandpeix, 23/06/2003, 10.36.17    *
***************************************************************
*
      implicit none

#include "cv3param.h"

c inputs:
      integer nloc, len, ncum, nd, ntra
      integer idcum(nloc)
      integer iflag(nloc),kbas(nloc),ktop(nloc)
      real precip(nloc),cbmf(nloc)
      real sig(nloc,nd), w0(nloc,nd),ptop2(nloc)
      real ft(nloc,nd), fq(nloc,nd), fu(nloc,nd), fv(nloc,nd)
      real ftra(nloc,nd,ntra)
      real Ma(nloc,nd),mip(nloc,nd),Vprecip(nloc,nd+1)
      real upwd(nloc,nd),dnwd(nloc,nd),dnwd0(nloc,nd)
      real qcondc(nloc,nd)
      real wd(nloc),cape(nloc),cin(nloc)
      real tvp(nloc,nd)
      real ftd(nloc,nd), fqd(nloc,nd)
      real Plim1(nloc),Plim2(nloc)
      real asupmax(nloc,nd),supmax0(nloc)
      real asupmaxmin(nloc)

c outputs:
      integer iflag1(len),kbas1(len),ktop1(len)
      real precip1(len),cbmf1(len)
      real sig1(len,nd), w01(len,nd),ptop21(len)
      real ft1(len,nd), fq1(len,nd), fu1(len,nd), fv1(len,nd)
      real ftra1(len,nd,ntra)
      real Ma1(len,nd),mip1(len,nd),Vprecip1(len,nd+1)
      real upwd1(len,nd),dnwd1(len,nd),dnwd01(len,nd)
      real qcondc1(len,nd)
      real wd1(len),cape1(len),cin1(len)
      real tvp1(len,nd)
      real ftd1(len,nd), fqd1(len,nd)
      real Plim11(len),Plim21(len)
      real asupmax1(len,nd),supmax01(len)
      real asupmaxmin1(len)
c
c local variables:
      integer i,k,j,k1,k2

        do 2000 i=1,ncum
         ptop21(idcum(i))=ptop2(i)
         precip1(idcum(i))=precip(i)
         cbmf1(idcum(i))=cbmf(i)
         iflag1(idcum(i))=iflag(i)
         kbas1(idcum(i))=kbas(i)
         ktop1(idcum(i))=ktop(i)
         wd1(idcum(i))=wd(i)
         cape1(idcum(i))=cape(i)
         cin1(idcum(i))=cin(i)
         Plim11(idcum(i))=Plim1(i)
         Plim21(idcum(i))=Plim2(i)
         supmax01(idcum(i))=supmax0(i)
         asupmaxmin1(idcum(i))=asupmaxmin(i)
 2000   continue

        do 2020 k=1,nd
          do 2010 i=1,ncum
            sig1(idcum(i),k)=sig(i,k)
            w01(idcum(i),k)=w0(i,k)
            ft1(idcum(i),k)=ft(i,k)
            fq1(idcum(i),k)=fq(i,k)
            fu1(idcum(i),k)=fu(i,k)
            fv1(idcum(i),k)=fv(i,k)
            Ma1(idcum(i),k)=Ma(i,k)
            mip1(idcum(i),k)=mip(i,k)
            Vprecip1(idcum(i),k)=Vprecip(i,k)
            upwd1(idcum(i),k)=upwd(i,k)
            dnwd1(idcum(i),k)=dnwd(i,k)
            dnwd01(idcum(i),k)=dnwd0(i,k)
            qcondc1(idcum(i),k)=qcondc(i,k)
            tvp1(idcum(i),k)=tvp(i,k)
            ftd1(idcum(i),k)=ftd(i,k)
            fqd1(idcum(i),k)=fqd(i,k)
            asupmax1(idcum(i),k)=asupmax(i,k)
 2010     continue
 2020   continue

        do 2040 i=1,ncum
          sig1(idcum(i),nd)=sig(i,nd)
2040    continue


        do 2100 j=1,ntra
c oct3         do 2110 k=1,nl
         do 2110 k=1,nd ! oct3
          do 2120 i=1,ncum
            ftra1(idcum(i),k,j)=ftra(i,k,j)
 2120     continue
 2110    continue
 2100   continue
c
c        do 2220 k2=1,nd
c         do 2210 k1=1,nd
c          do 2200 i=1,ncum
c            ment1(idcum(i),k1,k2) = ment(i,k1,k2)
c            sij1(idcum(i),k1,k2) = sij(i,k1,k2)
c2200      enddo
c2210     enddo
c2220    enddo

       RETURN
      END

