!
! $Id: juldate.F 2197 2015-02-09 07:13:05Z emillour $
!
	subroutine juldate(ian,imoi,ijou,oh,om,os,tjd,tjdsec)
c	Sous-routine de changement de date:
c	gregorien>>>date julienne
c	En entree:an,mois,jour,heure,min.,sec.
c	En sortie:tjd
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: ian,imoi,ijou,oh,om,os
        REAL,INTENT(OUT) :: tjd,tjdsec
        
        REAL frac,year,rmon,cf,a,b
        INTEGER ojou
        
	frac=((os/60.+om)/60.+oh)/24.
	ojou=dble(ijou)+frac
	    year=dble(ian)
	    rmon=dble(imoi)
	if (imoi .le. 2) then
	    year=year-1.
	    rmon=rmon+12.
	endif
	cf=year+(rmon/100.)+(ojou/10000.)
	if (cf .ge. 1582.1015) then
	    a=int(year/100)
	    b=2-a+int(a/4)
	else
	    b=0
	endif
	tjd=int(365.25*year)+int(30.6001*(rmon+1))+int(ojou)
     +   +1720994.5+b
        tjdsec=(ojou-int(ojou))+(tjd-int(tjd))
        tjd=int(tjd)+int(tjdsec)
	tjdsec=tjdsec-int(tjdsec)
	return
	end


