c 	************************************************************
c
c		chapman function for grazing incidence
c
c 	************************************************************

	function ch(x,z) 

	integer ierr
	real ch,x

	double precision ch1,sublim,z,error,aerr,rerr
	double precision dcadre,fch
	external fch,dcadre 

	common/integ/ax,az
	double precision ax,az 

c

	ax=x
	az=z

	aerr= 0.0d0
	rerr= 0.0001d0 
	sublim = 1.0d-4
	ch1 = dcadre( fch, sublim, az, aerr, rerr, error, ierr )
	ch = ax * sin(az) * ch1

	return
	end

c 	************************************************************
c
c		integrand for the chapman function 
c
c 	************************************************************

	double precision function fch(gi)

	double precision gi

	common/integ/x,z
	double precision x,z 

!	real x,z,a

	fch = (1./sin(gi))**2. * exp( x *(1.- sin(z)/sin(gi)) )

	return 
	end




******************************************************************************

      double precision function dcadre (f,a,b,aerr,rerr,error,ier)              
c                                  specifications for arguments                 
      integer            ier                                                    
      double precision   f,a,b,aerr,rerr,error                                  
c                                  specifications for local variables           
      integer            ibegs(30),maxts,maxtbl,mxstge,ibeg,ii,nnleft           
      integer            i,n2,iii,istep2,iend,istep,l,lm1,it,istage,n           
      double precision   t(10,10),r(10),ait(10),dif(10),rn(4),ts(2049)          
      double precision   begin(30),finis(30),est(30)                            
      double precision   h2tol,aittol,length,jumptl,zero,p1,half,one            
      double precision   two,four,fourp5,ten,hun,cadre,aitlow                   
      double precision   stepmn,stepnm,stage,curest,fnsize,hrerr                
      double precision   prever,beg,fbeg,edn,fend,step,astep,tabs,hovn          
      double precision   fn,sum,sumabs,vint,tabtlm,ergl,ergoal                  
      double precision   erra,errr,fextrp,errer,diff,sing,fextm1                
      double precision   h2next,singnx,slope,fbeg2,erret,h2tfex,fi              
      logical            h2conv,aitken,right,reglar,reglsv(30)                  
      data               aitlow,h2tol,aittol,jumptl,maxts,maxtbl,mxstge         
     1                   /1.1d0,.15d0,.1d0,.01d0,2049,10,30/                    
      data               rn(1),rn(2),rn(3),rn(4)/                               
     1                   .7142005d0,.3466282d0,.843751d0,.1263305d0/            
      data               zero,p1,half,one,two,four,fourp5,ten,hun               
     1                   /0.0d0,0.1d0,0.5d0,1.0d0,2.0d0,4.0d0,                  
     2                   4.5d0,10.0d0,100.0d0/                                  
c                                  first executable statement                   
      ier = 0                                                                   
      cadre = zero                                                              
      error = zero                                                              
      curest = zero                                                             
      vint = zero                                                               
      length = dabs(b-a)

      if (length .eq. zero) go to 215                                           
      if (rerr .gt. p1 .or. rerr .lt. zero) go to 210                           
      hrerr = rerr+hun                                                          
      if (aerr .eq. zero .and. hrerr .le. hun) go to 210                        
      errr = rerr                                                               
      erra = dabs(aerr)                                                         
      stepmn = length/(two**mxstge)                                             
      stepnm = dmax1(length,dabs(a),dabs(b))*ten                                
      stage = half                                                              
      istage = 1                                                                
      fnsize = zero                                                             
      prever = zero                                                             
      reglar = .false.                                                          
c                                  the given interval of integration            
c                                    is the first interval considered.          
      beg = a                                                                   
      fbeg = f(beg)*half                                                        
      ts(1) = fbeg                                                              
      ibeg = 1                                                                  
      edn = b                                                                   
      fend = f(edn)*half                                                        
      ts(2) = fend                                                              
      iend = 2                                                                  
    5 right = .false.                                                           
c                                  investigation of a particular                
c                                    subinterval begins at this point.          
   10 step = edn - beg                                                          
      astep =  dabs(step)                                                       
      if (astep .lt. stepmn) go to 205                                          
      hrerr = stepnm+astep                                                      
      if (hrerr .eq. stepnm) go to 205                                          
      t(1,1) = fbeg + fend                                                      
      tabs = dabs(fbeg) + dabs(fend)                                            
      l = 1                                                                     
      n = 1                                                                     
      h2conv = .false.                                                          
      aitken = .false.                                                          
   15 lm1 = l                                                                   
      l = l + 1                                                                 
c                                  calculate the next trapezoid sum,            
c                                    t(l,1), which is based on *n2* + 1         
c                                    equispaced points. here,                   
c                                    n2 = n*2 = 2**(l-1).                       
      n2 = n+n                                                                  
      fn = n2                                                                   
      istep = (iend - ibeg)/n                                                   
      if (istep .gt. 1) go to 25                                                
      ii = iend                                                                 
      iend = iend + n                                                           
      if (iend .gt. maxts) go to 200                                            
      hovn = step/fn                                                            
      iii = iend                                                                
      fi = one                                                                  
      do 20 i=1,n2,2                                                            
         ts(iii) = ts(ii)                                                       
         ts(iii-1) = f(edn - fi * hovn)                                         
         fi = fi+two                                                            
         iii = iii-2                                                            
         ii = ii-1                                                              
   20 continue                                                                  
      istep = 2                                                                 
   25 istep2 = ibeg + istep/2                                                   
      sum = zero                                                                
      sumabs = zero                                                             
      do 30 i=istep2,iend,istep                                                 
         sum = sum + ts(i)                                                      
         sumabs = sumabs + dabs(ts(i))                                          
   30 continue                                                                  
      t(l,1) = t(l-1,1)*half+sum/fn                                             
      tabs = tabs*half+sumabs/fn                                                
      n = n2                                                                    
c                                  get preliminary value for *vint*             
c                                    from last trapezoid sum and update         
c                                    the error requirement *ergoal*             
c                                    for this subinterval.                      
      it = 1                                                                    
      vint = step*t(l,1)                                                        
      tabtlm = tabs*ten                                                         
      fnsize = dmax1(fnsize,dabs(t(l,1)))                                       
      ergl = astep*fnsize*ten                                                   
      ergoal = stage*dmax1(erra,errr*dabs(curest+vint))                         
c                                  complete row l and column l of *t*           
c                                    array.                                     
      fextrp = one                                                              
      do 35 i=1,lm1                                                             
         fextrp = fextrp*four                                                   
         t(i,l) = t(l,i) - t(l-1,i)                                             
         t(l,i+1) = t(l,i) + t(i,l)/(fextrp-one)                                
   35 continue                                                                  
      errer = astep*dabs(t(1,l))                                                
c                                  preliminary decision procedure               
c                                    if l = 2 and t(2,1) = t(1,1),              
c                                    go to 135 to follow up the                 
c                                    impression that intergrand is              
c                                    straight line.                             
      if (l .gt. 2) go to 40                                                    
      hrerr = tabs+p1*dabs(t(1,2))                                              
      if (hrerr .eq. tabs) go to 135                                            
      go to 15                                                                  
c                                  caculate next ratios for                     
c                                    columns 1,...,l-2 of t-table               
c                                    ratio is set to zero if difference         
c                                    in last two entries of column is           
c                                    about zero                                 
   40 do 45 i=2,lm1                                                             
         diff = zero                                                            
         hrerr = tabtlm+dabs(t(i-1,l))                                          
         if (hrerr .ne. tabtlm) diff = t(i-1,lm1)/t(i-1,l)                      
         t(i-1,lm1) = diff                                                      
   45 continue                                                                  
      if (dabs(four-t(1,lm1)) .le. h2tol) go to 60                              
      if (t(1,lm1) .eq. zero) go to 55                                          
      if (dabs(two-dabs(t(1,lm1))) .lt. jumptl) go to 130                       
      if (l .eq. 3) go to 15                                                    
      h2conv = .false.                                                          
      if (dabs((t(1,lm1)-t(1,l-2))/t(1,lm1)) .le. aittol) go to 75              
   50 if (reglar) go to 55                                                      
      if (l .eq. 4) go to 15                                                    
      hrerr = ergl+errer                                                        
   55 if (errer .gt. ergoal .and. hrerr .ne. ergl) go to 175                    
      go to 145                                                                 
c                                  cautious romberg extrapolation               
   60 if (h2conv) go to 65                                                      
      aitken = .false.                                                          
      h2conv = .true.                                                           
   65 fextrp = four                                                             
   70 it = it + 1                                                               
      vint = step*t(l,it)                                                       
      errer = dabs(step/(fextrp-one)*t(it-1,l))                                 
      if (errer .le. ergoal) go to 160                                          
      hrerr = ergl+errer                                                        
      if (hrerr .eq. ergl) go to 160                                            
      if (it .eq. lm1) go to 125                                                
      if (t(it,lm1) .eq. zero) go to 70                                         
      if (t(it,lm1) .le. fextrp) go to 125                                      
      if (dabs(t(it,lm1)/four-fextrp)/fextrp .lt. aittol)                       
     1       fextrp = fextrp*four                                               
      go to 70                                                                  
c                                  integrand may have x**alpha type             
c                                    singularity                                
c                                    resulting in a ratio of *sing*  =          
c                                    2**(alpha + 1)                             
   75 if (t(1,lm1) .lt. aitlow) go to 175                                       
      if (aitken) go to 80                                                      
      h2conv = .false.                                                          
      aitken = .true.                                                           
   80 fextrp = t(l-2,lm1)                                                       
      if (fextrp .gt. fourp5) go to 65                                          
      if (fextrp .lt. aitlow) go to 175                                         
      if (dabs(fextrp-t(l-3,lm1))/t(1,lm1) .gt. h2tol) go to 175                
      sing = fextrp                                                             
      fextm1 = one/(fextrp - one)                                               
      ait(1) = zero                                                             
      do 85 i=2,l                                                               
         ait(i) = t(i,1) + (t(i,1)-t(i-1,1))*fextm1                             
         r(i) = t(1,i-1)                                                        
         dif(i) = ait(i) - ait(i-1)                                             
   85 continue                                                                  
      it = 2                                                                    
   90 vint = step*ait(l)                                                        
      errer = errer*fextm1                                                      
      hrerr = ergl+errer                                                        
      if (errer .gt. ergoal .and. hrerr .ne. ergl) go to 95                     
      ier = max0(ier,65)                                                        
      go to 160                                                                 
   95 it = it + 1                                                               
      if (it .eq. lm1) go to 125                                                
      if (it .gt. 3) go to 100                                                  
      h2next = four                                                             
      singnx = sing+sing                                                        
  100 if (h2next .lt. singnx) go to 105                                         
      fextrp = singnx                                                           
      singnx = singnx+singnx                                                    
      go to 110                                                                 
  105 fextrp = h2next                                                           
      h2next = four*h2next                                                      
  110 do 115 i=it,lm1                                                           
         r(i+1) = zero                                                          
         hrerr = tabtlm+dabs(dif(i+1))                                          
         if (hrerr .ne. tabtlm) r(i+1) = dif(i)/dif(i+1)                        
  115 continue                                                                  
      h2tfex = -h2tol*fextrp                                                    
      if (r(l) - fextrp .lt. h2tfex) go to 125                                  
      if (r(l-1)-fextrp .lt. h2tfex) go to 125                                  
      errer = astep*dabs(dif(l))                                                
      fextm1 = one/(fextrp - one)                                               
      do 120 i=it,l                                                             
         ait(i) = ait(i) + dif(i)*fextm1                                        
         dif(i) = ait(i) - ait(i-1)                                             
  120 continue                                                                  
      go to 90                                                                  
c                                  current trapezoid sum and resulting          
c                                    extrapolated values did not give           
c                                    a small enough *errer*.                    
c                                    note -- having prever .lt. errer           
c                                    is an almost certain sign of               
c                                    beginning trouble with in the func-        
c                                    tion values. hence, a watch for,           
c                                    and control of, noise should               
c                                    begin here.                                
  125 fextrp = dmax1(prever/errer,aitlow)                                       
      prever = errer                                                            
      if (l .lt. 5) go to 15                                                    
      if (l-it .gt. 2 .and. istage .lt. mxstge) go to 170                       
      erret = errer/(fextrp**(maxtbl-l))                                        
      hrerr = ergl+erret                                                        
      if (erret .gt. ergoal .and. hrerr .ne. ergl) go to 170                    
      go to 15                                                                  
c                                  integrand has jump (see notes)               
  130 hrerr = ergl+errer                                                        
      if (errer .gt. ergoal .and. hrerr .ne. ergl) go to 170                    
c                                  note that  2*fn = 2**l                       
      diff = dabs(t(1,l))*(fn+fn)                                               
      go to 160                                                                 
c                                  integrand is straight line                   
c                                    test this assumption by comparing          
c                                    the value of the integrand at              
c                                    four *randomly chosen* points with         
c                                    the value of the straight line             
c                                    interpolating the integrand at the         
c                                    two end points of the sub-interval.        
c                                    if test is passed, accept *vint*           
  135 slope = (fend-fbeg)*two                                                   
      fbeg2 = fbeg+fbeg                                                         
      do 140 i=1,4                                                              
         diff = dabs(f(beg+rn(i)*step) - fbeg2-rn(i)*slope)                     
         hrerr = tabtlm+diff                                                    
         if(hrerr .ne. tabtlm) go to 155                                        
  140 continue                                                                  
      go to 160                                                                 
c                                  noise may be dominant feature                
c                                    estimate noise level by comparing          
c                                    the value of the integrand at              
c                                    four *randomly chosen* points with         
c                                    the value of the straight line             
c                                    interpolating the integrand at the         
c                                    two endpoints. if small enough,            
c                                    accept *vint*                              
  145 slope = (fend-fbeg)*two                                                   
      fbeg2 = fbeg+fbeg                                                         
      i = 1                                                                     
  150 diff = dabs(f(beg+rn(i)*step) - fbeg2-rn(i)*slope)                        
  155 errer = dmax1(errer,astep*diff)                                           
      hrerr = ergl+errer                                                        
      if (errer .gt. ergoal .and. hrerr .ne. ergl) go to 175                    
      i = i+1                                                                   
      if (i .le. 4) go to 150                                                   
      ier = 66                                                                  
c                                  intergration over current sub-               
c                                    interval successful                        
c                                    add *vint* to *dcadre* and *errer*         
c                                    to *error*, then set up next sub-          
c                                    interval, if any.                          
  160 cadre = cadre + vint                                                      
      error = error + errer                                                     
      if (right) go to 165                                                      
      istage = istage - 1                                                       
      if (istage .eq. 0) go to 220                                              
      reglar = reglsv(istage)                                                   
      beg = begin(istage)                                                       
      edn = finis(istage)                                                       
      curest = curest - est(istage+1) + vint                                    
      iend = ibeg - 1                                                           
      fend = ts(iend)                                                           
      ibeg = ibegs(istage)                                                      
      go to 180                                                                 
  165 curest = curest + vint                                                    
      stage = stage+stage                                                       
      iend = ibeg                                                               
      ibeg = ibegs(istage)                                                      
      edn = beg                                                                 
      beg = begin(istage)                                                       
      fend = fbeg                                                               
      fbeg = ts(ibeg)                                                           
      go to 5                                                                   
c                                  integration over current subinterval         
c                                    is unsuccessful. mark subinterval          
c                                    for further subdivision. set up            
c                                    next subinterval.                          
  170 reglar = .true.                                                           
  175 if (istage .eq. mxstge) go to 205                                         
      if (right) go to 185                                                      
      reglsv(istage+1) = reglar                                                 
      begin(istage) = beg                                                       
      ibegs(istage) = ibeg                                                      
      stage = stage*half                                                        
  180 right = .true.                                                            
      beg = (beg+edn)*half                                                      
      ibeg = (ibeg+iend)/2                                                      
      ts(ibeg) = ts(ibeg)*half                                                  
      fbeg = ts(ibeg)                                                           
      go to 10                                                                  
  185 nnleft = ibeg - ibegs(istage)                                             
      if (iend+nnleft .ge. maxts) go to 200                                     
      iii = ibegs(istage)                                                       
      ii = iend                                                                 
      do 190 i=iii,ibeg                                                         
         ii = ii + 1                                                            
         ts(ii) = ts(i)                                                         
  190 continue                                                                  
      do 195 i=ibeg,ii                                                          
         ts(iii) = ts(i)                                                        
         iii = iii + 1                                                          
  195 continue                                                                  
      iend = iend + 1                                                           
      ibeg = iend - nnleft                                                      
      fend = fbeg                                                               
      fbeg = ts(ibeg)                                                           
      finis(istage) = edn                                                       
      edn = beg                                                                 
      beg = begin(istage)                                                       
      begin(istage) = edn                                                       
      reglsv(istage) = reglar                                                   
      istage = istage + 1                                                       
      reglar = reglsv(istage)                                                   
      est(istage) = vint                                                        
      curest = curest + est(istage)                                             
      go to 5                                                                   
c                                  failure to handle given integra-             
c                                    tion problem                               
  200 ier = 131                                                                 
      go to 215                                                                 
  205 ier = 132                                                                 
      go to 215                                                                 
  210 ier = 133                                                                 
  215 cadre = curest + vint                                                     
  220 dcadre = cadre                                                            
 9000 continue                                                                  
      !if (ier .ne. 0) call uertst (ier,6hdcadre)                                
 9005 return                                                                    
      end     


!******************************************************************************
!
!	subroutine uertst (ier,name)                                              
!c                                  specifications for arguments                 
!      integer            ier                                                    
!      integer            name(2)                                                
!c                                  specifications for local variables           
!      integer            i,ieq,ieqdf,iounit,level,levold,nameq(6),              
!     *                   namset(6),namupk(6),nin,nmtb                           
!      data               namset/1hu,1he,1hr,1hs,1he,1ht/                        
!      data               nameq/6*1h /                                           
!      data               level/4/,ieqdf/0/,ieq/1h=/                             
!c                                  unpack name into namupk                      
!c                                  first executable statement                   
!      call uspkd (name,6,namupk,nmtb)                                           
!c                                  get output unit number                       
!      call ugetio(1,nin,iounit)                                                 
!c                                  check ier                                    
!      if (ier.gt.999) go to 25                                                  
!      if (ier.lt.-32) go to 55                                                  
!      if (ier.le.128) go to 5                                                   
!      if (level.lt.1) go to 30                                                  
!c                                  print terminal message                       
!      if (ieqdf.eq.1) write(iounit,35) ier,nameq,ieq,namupk                     
!      if (ieqdf.eq.0) write(iounit,35) ier,namupk                               
!      go to 30                                                                  
!    5 if (ier.le.64) go to 10                                                   
!      if (level.lt.2) go to 30                                                  
!c                                  print warning with fix message               
!c      if (ieqdf.eq.1) write(iounit,40) ier,nameq,ieq,namupk                     
!c      if (ieqdf.eq.0) write(iounit,40) ier,namupk                               
!      if (ieqdf.eq.1) continue 
!      if (ieqdf.eq.0) continue 
!      go to 30                                                                  
!   10 if (ier.le.32) go to 15                                                   
!c                                  print warning message                        
!      if (level.lt.3) go to 30                                                  
!      if (ieqdf.eq.1) write(iounit,45) ier,nameq,ieq,namupk                     
!      if (ieqdf.eq.0) write(iounit,45) ier,namupk                               
!      go to 30                                                                  
!   15 continue                                                                  
!c                                  check for uerset call                        
!      do 20 i=1,6                                                               
!         if (namupk(i).ne.namset(i)) go to 25                                   
!   20 continue                                                                  
!      levold = level                                                            
!      level = ier                                                               
!      ier = levold                                                              
!      if (level.lt.0) level = 4                                                 
!      if (level.gt.4) level = 4                                                 
!      go to 30                                                                  
!   25 continue                                                                  
!      if (level.lt.4) go to 30                                                  
!c                                  print non-defined message                    
!      if (ieqdf.eq.1) write(iounit,50) ier,nameq,ieq,namupk                     
!      if (ieqdf.eq.0) write(iounit,50) ier,namupk                               
!   30 ieqdf = 0                                                                 
!      return                                                                    
!   35 format(19h *** terminal error,10x,7h(ier = ,i3,                           
!     1       20h) from imsl routine ,6a1,a1,6a1)                                
!   40 format(27h *** warning with fix error,2x,7h(ier = ,i3,                    
!     1       20h) from imsl routine ,6a1,a1,6a1)                                
!   45 format(18h *** warning error,11x,7h(ier = ,i3,                            
!     1       20h) from imsl routine ,6a1,a1,6a1)                                
!   50 format(20h *** undefined error,9x,7h(ier = ,i5,                           
!     1       20h) from imsl routine ,6a1,a1,6a1)                                
!c                                                                               
!c                                  save p for p = r case                        
!c                                    p is the page namupk                       
!c                                    r is the routine namupk                    
!   55 ieqdf = 1                                                                 
!      do 60 i=1,6                                                               
!   60 nameq(i) = namupk(i)                                                      
!   65 return                                                                    
!      end                          


************************************************************************

      subroutine uspkd  (packed,nchars,unpakd,nchmtb)                           
c                                  specifications for arguments                 
      integer            nc,nchars,nchmtb                                       
c                                                                               
c      integer            unpakd(1),iblank                                       
      integer            unpakd(nchars),iblank                                       
      character          packed(1)           ! Modificado el 29-Ago-2001 porque
!      integer            packed(1)            esto deberia ser CHARACTER, no??
      data               iblank /1h /                                           
c                                  initialize nchmtb                            
      nchmtb = 0                                                                
c                                  return if nchars is le zero                  
      if(nchars.le.0) return                                                    
c                                  set nc=number of chars to be decoded         
      nc = min0 (129,nchars)                                                   
!      decode (nc,150,packed) (unpakd(i),i=1,nc)                                
      read (packed,150) (unpakd(i),i=1,nc)                                
  150 format (129a1)                                                           
c                                  check unpakd array and set nchmtb            
c                                  based on trailing blanks found               
      do 200 n = 1,nc                                                           
         nn = nc - n + 1                                                        
         if(unpakd(nn) .ne. 536870912) go to 210                                
  200 continue                                                                  
  210 nchmtb = nn                                                               
      return                                                                    
      end            


****************************************************************************

      subroutine ugetio(iopt,nin,nout)                                          
c                                  specifications for arguments                 
      integer            iopt,nin,nout                                          
c                                  specifications for local variables           
      integer            nind,noutd                                             
      data               nind/5/,noutd/6/                                       
c                                  first executable statement                   
      if (iopt.eq.3) go to 10                                                   
      if (iopt.eq.2) go to 5                                                    
      if (iopt.ne.1) go to 9005                                                 
      nin = nind                                                                
      nout = noutd                                                              
      go to 9005                                                                
    5 nind = nin                                                                
      go to 9005                                                                
   10 noutd = nout                                                              
 9005 return                                                                    
      end       
