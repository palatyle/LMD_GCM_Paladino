      SUBROUTINE cv3a_compress( len,nloc,ncum,nd,ntra
     :    ,iflag1,nk1,icb1,icbs1
     :    ,plcl1,tnk1,qnk1,gznk1,hnk1,unk1,vnk1
     :    ,wghti1,pbase1,buoybase1
     :    ,t1,q1,qs1,t1_wake,q1_wake,qs1_wake,s1_wake
     :    ,u1,v1,gz1,th1,th1_wake
     :    ,tra1
     :    ,h1     ,lv1     ,cpn1   ,p1,ph1,tv1    ,tp1,tvp1,clw1
     :    ,h1_wake,lv1_wake,cpn1_wake     ,tv1_wake
     :    ,sig1,w01,ptop21
     :    ,Ale1,Alp1
     o    ,iflag,nk,icb,icbs
     o    ,plcl,tnk,qnk,gznk,hnk,unk,vnk
     o    ,wghti,pbase,buoybase
     o    ,t,q,qs,t_wake,q_wake,qs_wake,s_wake
     o    ,u,v,gz,th,th_wake
     o    ,tra
     o    ,h     ,lv     ,cpn    ,p,ph,tv    ,tp,tvp,clw
     o    ,h_wake,lv_wake,cpn_wake    ,tv_wake
     o    ,sig,w0,ptop2
     o    ,Ale,Alp  )
***************************************************************
*                                                             *
* CV3A_COMPRESS                                               *
*                                                             *
*                                                             *
* written by   : Sandrine Bony-Lena , 17/05/2003, 11.22.15    *
* modified by  : Jean-Yves Grandpeix, 23/06/2003, 10.28.09    *
***************************************************************
*
      implicit none

#include "cv3param.h"

c inputs:
      integer len,nloc,ncum,nd,ntra
      integer iflag1(len),nk1(len),icb1(len),icbs1(len)
      real plcl1(len),tnk1(len),qnk1(len),gznk1(len)
      real hnk1(len),unk1(len),vnk1(len)
      real wghti1(len,nd),pbase1(len),buoybase1(len)
      real t1(len,nd),q1(len,nd),qs1(len,nd)
      real t1_wake(len,nd),q1_wake(len,nd),qs1_wake(len,nd)
      real s1_wake(len)
      real u1(len,nd),v1(len,nd)
      real gz1(len,nd),th1(len,nd),th1_wake(len,nd)
      real tra1(len,nd,ntra)
      real h1(len,nd),lv1(len,nd),cpn1(len,nd)
      real p1(len,nd),ph1(len,nd+1),tv1(len,nd),tp1(len,nd)
      real tvp1(len,nd),clw1(len,nd)
      real h1_wake(len,nd),lv1_wake(len,nd),cpn1_wake(len,nd)
      real tv1_wake(len,nd)
      real sig1(len,nd), w01(len,nd), ptop21(len)
      real Ale1(len),Alp1(len)

c outputs:
c en fait, on a nloc=len pour l'instant (cf cv_driver)
      integer iflag(len),nk(len),icb(len),icbs(len)
      real plcl(len),tnk(len),qnk(len),gznk(len)
      real hnk(len),unk(len),vnk(len)
      real wghti(len,nd),pbase(len),buoybase(len)
      real t(len,nd),q(len,nd),qs(len,nd)
      real t_wake(len,nd),q_wake(len,nd),qs_wake(len,nd)
      real s_wake(len)
      real u(len,nd),v(len,nd)
      real gz(len,nd),th(len,nd),th_wake(len,nd)
      real tra(len,nd,ntra)
      real h(len,nd),lv(len,nd),cpn(len,nd)
      real p(len,nd),ph(len,nd+1),tv(len,nd),tp(len,nd)
      real tvp(len,nd),clw(len,nd)
      real h_wake(len,nd),lv_wake(len,nd),cpn_wake(len,nd)
      real tv_wake(len,nd)
      real sig(len,nd), w0(len,nd), ptop2(len)
      real Ale(len),Alp(len)

c local variables:
      integer i,k,nn,j

      CHARACTER (LEN=20) :: modname='cv3a_compress'
      CHARACTER (LEN=80) :: abort_message


      do 110 k=1,nl+1
       nn=0
      do 100 i=1,len
      if(iflag1(i).eq.0)then
        nn=nn+1
        wghti(nn,k)=wghti1(i,k)
        t(nn,k)=t1(i,k)
        q(nn,k)=q1(i,k)
        qs(nn,k)=qs1(i,k)
        t_wake(nn,k)=t1_wake(i,k)
        q_wake(nn,k)=q1_wake(i,k)
        qs_wake(nn,k)=qs1_wake(i,k)
        u(nn,k)=u1(i,k)
        v(nn,k)=v1(i,k)
        gz(nn,k)=gz1(i,k)
        th(nn,k)=th1(i,k)
        th_wake(nn,k)=th1_wake(i,k)
        h(nn,k)=h1(i,k)
        lv(nn,k)=lv1(i,k)
        cpn(nn,k)=cpn1(i,k)
        p(nn,k)=p1(i,k)
        ph(nn,k)=ph1(i,k)
        tv(nn,k)=tv1(i,k)
        tp(nn,k)=tp1(i,k)
        tvp(nn,k)=tvp1(i,k)
        clw(nn,k)=clw1(i,k)
        h_wake(nn,k)=h1_wake(i,k)
        lv_wake(nn,k)=lv1_wake(i,k)
        cpn_wake(nn,k)=cpn1_wake(i,k)
        tv_wake(nn,k)=tv1_wake(i,k)
        sig(nn,k)=sig1(i,k)
        w0(nn,k)=w01(i,k)
      endif
 100    continue
 110  continue

      do 121 j=1,ntra
ccccc      do 111 k=1,nl+1
      do 111 k=1,nd
       nn=0
      do 101 i=1,len
      if(iflag1(i).eq.0)then
       nn=nn+1
       tra(nn,k,j)=tra1(i,k,j)
      endif
 101  continue
 111  continue
 121  continue

      if (nn.ne.ncum) then
        print*,'WARNING nn not equal to ncum: ',nn,ncum
        abort_message = ''
        CALL abort_gcm (modname,abort_message,1)
      endif

      nn=0
      do 150 i=1,len
      if(iflag1(i).eq.0)then
      nn=nn+1
      s_wake(nn)=s1_wake(i)
      iflag(nn)=iflag1(i)
      nk(nn)=nk1(i)
      icb(nn)=icb1(i)
      icbs(nn)=icbs1(i)
      plcl(nn)=plcl1(i)
      tnk(nn)=tnk1(i)
      qnk(nn)=qnk1(i)
      gznk(nn)=gznk1(i)
      hnk(nn)=hnk1(i)
      unk(nn)=unk1(i)
      vnk(nn)=vnk1(i)
      pbase(nn)=pbase1(i)
      buoybase(nn)=buoybase1(i)
      ptop2(nn)=ptop2(i)
      ale(nn) = ale1(i)
      alp(nn) = alp1(i)
      endif
 150  continue

      if (nn.ne.ncum) then
         print*,'WARNING nn not equal to ncum: ',nn,ncum
         abort_message = ''
         CALL abort_gcm (modname,abort_message,1)
      endif

      RETURN
      END
