load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "visu-utils.ncl"

begin

;---------------------
;---- load the file
;---------------------

infile=addfile(inFile,"r")

;---------------------
;---- get variable list and dimension names
;---------------------

liste=getfilevarnames(infile)
nvars=dimsizes(liste)

;---- finding the dimensions 
tname="dummy"
pname="dummy"
lname="dummy"
Lname="dummy"
do i=0,nvars-1
  if ((liste(i).eq."time").or.(liste(i).eq."time_counter").or.(liste(i).eq."Time")) then
    tname=liste(i)
  end if
  if ((liste(i).eq."lev").or.(liste(i).eq."presnivs").or.(liste(i).eq."levgrnd")) then
    pname=liste(i)
  end if
  if ((liste(i).eq."lat").or.(liste(i).eq."latitude")) then
    lname=liste(i)
  end if
  if ((liste(i).eq."lon").or.(liste(i).eq."longitude")) then
    Lname=liste(i)
  end if
end do
if ((tname.eq."dummy").or.(pname.eq."dummy").or.(tname.eq."dummy").or.(tname.eq."dummy")) then
  print("One dimension was not found... See -h for help.")
  exit
else
  names=(/tname,pname,lname,Lname/)
end if

;---- temporal axis : for Titan, Ls ; for Venus, Vdays.
Vday=1.0087e7      ; Venus day
Tday=1.37889e6     ; Titan day
axels=False
timeax=infile->$tname$
timeaxis=tofloat(timeax-timeax(0))/Vday
copy_VarMeta(timeax,timeaxis)
timeaxis@units="Local days"
if (isfilevar(infile,"ls")) then
   delete(timeaxis)
   timeaxis=infile->ls
   axels=True
end if
if (isfilevar(infile,"Ls")) then
   delete(timeaxis)
   timeaxis=infile->Ls
   axels=True
end if
ndls=dimsizes(dimsizes(timeaxis))
if (ndls.ne.1) then
   delete(timeax)
   timeax=timeaxis(:,0,0)
   delete(timeaxis)
   timeaxis=timeax
end if

;---------------------
;---- print only list of variables
;---------------------
if (var.eq."liste") then
  print(liste)
  exit
end if

;==============================
; Prepare variable(s) for the plot
;==============================

;---------------------
;---- load variable(s)
;---------------------

  if (var.ne."custom") then
    variable=infile->$var$
    nbdim=dimsizes(filevardimsizes(infile,var))
  else
    variable=customVar(infile,labelvar)
    nbdim=dimsizes(dimsizes(variable))
  end if
  
  overplot=False
  diffdims=False
  if (var2.ne."dummy") then
     variable2=infile->$var2$
     nbdim2=dimsizes(filevardimsizes(infile,var2))
     if (nbdim2.ne.nbdim) then
       diffdims=True
     end if
     overplot=True
  end if
  
;---------------------
;---- trim variable to desired region
;---------------------
     dimname=new(nbdim,string)
     mindimval=new(nbdim,float)
     maxdimval=new(nbdim,float)
     axetim=(/0./)          ; dummy
     dimavg=new(4,logical)
     do i=0,3
       dimavg(i)=False
     end do

; locate dimensions of the variable
     flag=positionDims(variable,nbdim,names)
     if (diffdims) then
       flag2=positionDims(variable2,nbdim2,names)
       dimname2=new(nbdim2,string)
       mindimval2=new(nbdim2,float)
       maxdimval2=new(nbdim2,float)
       dimavg2=new(4,logical)
     end if

; store dim names, associated limits and flags
     index=0
     if (flag(0).ne.10) then
        dimname(index)=tname
        if (ntcutmin.ne.-999) then
          mindimval(index)=tofloat(variable&$tname$(ntcutmin))
          maxdimval(index)=tofloat(variable&$tname$(ntcutmax))
; A REVOIR !!!
     ; !! stupid bug, due to precision when time axis in double in file...
           if ((ntcutmin.eq.0).and.(ntcutmax.eq.0)) then
             mindimval(index)=mindimval(index)*0
             maxdimval(index)=maxdimval(index)*1.1
           end if
     ; !! stupid bug, due to precision when time axis in double in file...
        else
          mindimval(index)=-999
          maxdimval(index)=-999
        end if
        if (tavg.eq.-999) then
          dimavg(index)=True
        else
     ; trim time axis
          delete(axetim)
          axetim=timeaxis({$dimname(0)$|mindimval(0):maxdimval(0)})
        end if
        index=index+1
     end if
     if (flag(1).ne.10) then
        dimname(index)=pname
        mindimval(index)=pcutmin
        maxdimval(index)=pcutmax
        if (pavg.eq.-999) then
          dimavg(index)=True
        end if
        index=index+1
     end if
     if (flag(2).ne.10) then
        dimname(index)=lname
        mindimval(index)=lcutmin
        maxdimval(index)=lcutmax
        if (lavg.eq.-999) then
          dimavg(index)=True
        end if
        index=index+1
     end if
     if (flag(3).ne.10) then
        dimname(index)=Lname
        mindimval(index)=Lcutmin
        maxdimval(index)=Lcutmax
        if (Lavg.eq.-999) then
          dimavg(index)=True
        end if
     end if

; if global average, prepare limits
     do i=0,nbdim-1
        if (mindimval(i).eq.-999) then
          mindimval(i)=tofloat(variable&$dimname(i)$(0))
          maxdimval(i)=tofloat(variable&$dimname(i)$(dimsizes(variable&$dimname(i)$)-1))
        end if 
     end do

; if orthographic projection, default limits for lat and lon
     if (typeplot.eq."ortho") then
      do i=0,nbdim-1
        if ((dimname(i).eq.lname).or.(dimname(i).eq.Lname)) then
          mindimval(i)=tofloat(variable&$dimname(i)$(0))
          maxdimval(i)=tofloat(variable&$dimname(i)$(dimsizes(variable&$dimname(i)$)-1))
        end if 
      end do
      print("Latitude and Longitude are defaulted to full map.")
     end if
     
; if stereographic projection, default limits for lon
     if (typeplot.eq."stereo") then
      do i=0,nbdim-1
        if (dimname(i).eq.Lname) then
          mindimval(i)=tofloat(variable&$dimname(i)$(0))
          maxdimval(i)=tofloat(variable&$dimname(i)$(dimsizes(variable&$dimname(i)$)-1))
        end if 
      end do
      print("Longitude is defaulted to full map.")
     end if
     
; case of second variable
     if (.not.diffdims) then
        dimname2  =dimname
	mindimval2=mindimval
	maxdimval2=maxdimval
	dimavg2   =dimavg
     else
      index=0
      if (flag2(0).ne.10) then
        dimname2(index)=tname
        index=index+1
      end if
      if (flag2(1).ne.10) then
        dimname2(index)=pname
        index=index+1
      end if
      if (flag2(2).ne.10) then
        dimname2(index)=lname
        index=index+1
      end if
      if (flag2(3).ne.10) then
        dimname2(index)=Lname
      end if
      do i=0,nbdim2-1
        do j=0,nbdim-1
	  if (dimname2(i).eq.dimname(j)) then
	    mindimval2(i)=mindimval(j)
	    maxdimval2(i)=maxdimval(j)
	    dimavg2(i)   =dimavg(j)
	  end if
	end do
      end do
     end if
     
; trim variable(s)
     vartrim=trimVar(variable,nbdim,dimname,mindimval,maxdimval)

     if (overplot) then
       vartrim2=trimVar(variable2,nbdim2,dimname2,mindimval2,maxdimval2)
     end if

;---------------------
;---- make averages  
; (if limits are identical, average need to be made to trim the dimension)
;---------------------
     varavg=avgVar(vartrim,nbdim,dimavg)

     if (overplot) then
       varavg2=avgVar(vartrim2,nbdim2,dimavg2)
     end if

; after trim, redo dim sizes
     nbdimavg=dimsizes(dimsizes(varavg))
     if (overplot) then
       nbdimavg2=dimsizes(dimsizes(varavg2))
       if (nbdimavg2.ne.nbdimavg) then
         print("This case is not yet possible: the second variable must be plotted similarly to the first one...")
         overplot=False
       end if
     end if

;==============================
; Single point
;==============================

if ((nbdimavg.eq.1).and.(dimsizes(varavg).eq.1)) then
  print(mindimval)
  print(varavg)
  if (overplot) then
    print(varavg2)
  end if
  exit
end if

;==============================
;       TRACES
;==============================

; after trim, redo flag for dims
     flag=positionDims(varavg,nbdimavg,names)

; 2d maps (e.g. for ortho proj):
; revert longitudes if decreasing
  if ( (nbdimavg.eq.2).and.(varavg!0.eq.lname).and.(varavg!1.eq.Lname).and.(varavg&$Lname$(0).gt.varavg&$Lname$(1)) ) then
     newvaravg=revertLon_2dmap(varavg)
     delete(varavg)
     varavg=newvaravg
     if (overplot) then
        newvaravg=revertLon_2dmap(varavg2)
        delete(varavg2)
        varavg2=newvaravg
     end if
  end if
  
; Pressure axis used ?
  axep=False
  if (flag(1).ne.10) then
    axep=True
  end if

; Time axis used ?
  axet=False
  if (flag(0).ne.10) then
    axet=True
  end if

; limits for contours (2d) or variable axis (1d) ?
  limits=new(3,float)
  limits(0)=valmin
  limits(1)=valmax
  limits(2)=step
  if (overplot) then
    limits2=new(3,float)
    limits2(0)=valmin2
    limits2(1)=valmax2
    limits2(2)=step2
  end if

; output workstation
  wks   = gsn_open_wks(media,outputname)
  resbase          = True 
  resbase@gsnDraw  = False
  resbase@gsnFrame = False

;---------------------
if (typeplot.eq."coupe") then
;---------------------

; Too many dims
  if (nbdimavg.gt.2) then
    print("Too many dimensions for a simple 1D/2D plot...")
    exit
  end if

; 2D MAPS
  if (nbdimavg.eq.2) then
    res2d = resbase 
  ;Time axis used with lat or pressure: it must be on X 
    if ((axet).and.(flag(3).eq.10)) then
      revvar=inversDim_2d(varavg)
      delete(varavg)
      varavg=revvar
    end if
  ;defaults and options for plot
    res2d=deflimits_2d(varavg,res2d)
    res2d=contours_2d(res2d,limits)
  ;pressure axis
    if (axep) then
      res2d=axePress_2d(res2d)
    end if
  ;time axis
    if ((axet).and.(flag(3).eq.10)) then
      if (axels) then
        res2d=axeTimeXls_2d(res2d,axetim)
      else
        res2d=axeTimeX_2d(res2d,axetim)
      end if
    end if
    if ((axet).and.(flag(3).ne.10)) then
      if (axels) then
        res2d=axeTimeYls_2d(res2d,axetim)
      else
        res2d=axeTimeY_2d(res2d,axetim)
      end if
    end if
  ;creating map
    plotvar=Map_2d(varavg,wks,res2d)
  ;overlay
    if (overplot) then
      res2d_op=res2d
      res2d_op=contours_2d(res2d_op,limits2)
      res2d_op=optionsOverPlot_2d(res2d_op)
      plotvar2=OPlot_2d(varavg2,wks,res2d_op)
      overlay(plotvar,plotvar2)
    end if
  end if

; 1D PLOT
  if (nbdimavg.eq.1) then
    res1d = resbase 
    axename=getvardims(varavg)
    axe1d=varavg&$axename(0)$
  ;time axis
    if (axet) then
      delete(axe1d)
      axe1d=axetim
      if (axels) then
        res1d=axeTimels_1d(res1d,axe1d)
      else
        res1d=axeTime_1d(res1d,axe1d)
      end if
    end if
  ;creating plot
    plotvar=Plot_1d(varavg,wks,res1d,axe1d,axep,limits)
  end if

end if

;---------------------
if (typeplot.eq."ortho") then
;---------------------

; Need 2D
  if ((nbdimavg.ne.2).and.((flag(1).eq.10).or.(flag(2).eq.10))) then
    print("You need 2D lat/lon for orthographic projection...")
    exit
  end if

; options
  center  = (/Lcenter,lcenter/)
  grd     = grid
  spacing = (/gLspc,glspc/)
  gridcol = gcol
  
  res2d = resbase 
  res2d = contours_2d(res2d,limits)
  res2d = optionsOrtho(res2d,center,grd,spacing,gridcol)

;creating plot
  plotvar=Ortho_2d(varavg,wks,res2d)

end if

;---------------------
if (typeplot.eq."stereo") then
;---------------------

; Need 2D
  if ((nbdimavg.ne.2).and.((flag(1).eq.10).or.(flag(2).eq.10))) then
    print("You need 2D lat/lon for stereographic projection...")
    exit
  end if

; options
  if (lcenter.eq.-90) then
    hemisph="SH"
  else
    hemisph="NH"
  end if
  grd     = grid
  spacing = (/gLspc,glspc/)
  gridcol = gcol
  
  res2d = resbase 
  res2d = contours_2d(res2d,limits)
  res2d = optionsStereo(res2d,hemisph,grd,spacing,gridcol)
  
;creating plot
  plotvar=Stereo_2d(varavg,wks,res2d)

end if

;---------------------
  draw(plotvar)                   
  frame(wks)

end

