;;;;;;;;;;;; positionDims ;;;;;;;;;;;;;;;;

undef("positionDims")
function positionDims(var,nbdim,names)
;============================

; var   : variable
; nbdim : number of dimensions for var
; names : names of the 4 dimensions: time,levs,lat,lon

begin

  flag=(/10,10,10,10/)
  do i=0,nbdim-1
     if (var!i.eq.names(0)) then
        flag(0)=i
     end if
     if (var!i.eq.names(1)) then
        flag(1)=i
     end if
     if (var!i.eq.names(2)) then
        flag(2)=i
     end if
     if (var!i.eq.names(3)) then
        flag(3)=i
     end if
  end do

  return(flag)

end

;;;;;;;;;;;; trimVar ;;;;;;;;;;;;;;;;

undef("trimVar")
function trimVar(var,nbdim,dimname,mindimval,maxdimval)
;============================

; var      : variable
; nbdim    : number of dimensions for var
; dimname  : names of the dimensions for this var
; mindimval: min values for the dimensions
; maxdimval: max values for the dimensions

begin

  if (nbdim.eq.1) then
     vartrim=var({$dimname(0)$|mindimval(0):maxdimval(0)})
  end if
  if (nbdim.eq.2) then
       vartrim=var({$dimname(0)$|mindimval(0):maxdimval(0)},{$dimname(1)$|mindimval(1):maxdimval(1)})
  end if
  if (nbdim.eq.3) then
       vartrim=var({$dimname(0)$|mindimval(0):maxdimval(0)},{$dimname(1)$|mindimval(1):maxdimval(1)},{$dimname(2)$|mindimval(2):maxdimval(2)})
  end if
  if (nbdim.eq.4) then
       vartrim=var({$dimname(0)$|mindimval(0):maxdimval(0)},{$dimname(1)$|mindimval(1):maxdimval(1)},{$dimname(2)$|mindimval(2):maxdimval(2)},{$dimname(3)$|mindimval(3):maxdimval(3)})
;       vartrim=var($dimname(0)$|0:0,{$dimname(1)$|mindimval(1):maxdimval(1)},{$dimname(2)$|mindimval(2):maxdimval(2)},{$dimname(3)$|mindimval(3):maxdimval(3)})
  end if

  return(vartrim)

end

;;;;;;;;;;;; avgVar ;;;;;;;;;;;;;;;;

undef("avgVar")
function avgVar(vartrim,nbdim,dimavg)
;============================

; vartrim : variable
; nbdim   : number of dimensions for var
; dimavg  : logicals for averaging the dimensions

begin

  if ((nbdim.eq.4).and.(dimavg(3))) then
       vartrim_3=dim_avg_n_Wrap(vartrim,3)
  else
       vartrim_3=vartrim
  end if
  if ((nbdim.ge.3).and.(dimavg(2))) then
       vartrim_2=dim_avg_n_Wrap(vartrim_3,2)
  else
       vartrim_2=vartrim_3
  end if
  if ((nbdim.ge.2).and.(dimavg(1))) then
       vartrim_1=dim_avg_n_Wrap(vartrim_2,1)
  else
       vartrim_1=vartrim_2
  end if
  if ((nbdim.ge.1).and.(dimavg(0))) then
       varavg=dim_avg_n_Wrap(vartrim_1,0)
  else
       varavg=vartrim_1
  end if

  return(varavg)

end

;;;;;;;;;;;; deflimits_2d ;;;;;;;;;;;;;;;;

undef("deflimits_2d")
function deflimits_2d(var,res)
;============================

; var : 2D variable
; res : ressource

begin

  dims=dimsizes(var)
  dimname=getvardims(var)
  res@trXMinF = tofloat(var&$dimname(1)$(0))
  res@trXMaxF = tofloat(var&$dimname(1)$(dims(1)-1))
  res@trYMinF = tofloat(var&$dimname(0)$(0))
  res@trYMaxF = tofloat(var&$dimname(0)$(dims(0)-1))

  return(res)

end

;;;;;;;;;;;; contours_2d ;;;;;;;;;;;;;;;;

undef("contours_2d")
function contours_2d(res,limits)
;============================

; res    : modified ressource for the plot
; limits : min, max and values for contours (-888 => automatic)

begin

  if (limits(0).ne.-888) then
     res@cnLevelSelectionMode = "ManualLevels"
     res@cnMinLevelValF       = limits(0)
     res@cnMaxLevelValF       = limits(1)
     if (limits(2).ne.-888) then
        res@cnLevelSpacingF   = limits(2)
     end if
  end if

  return(res)

end

;;;;;;;;;;;; optionsOrtho ;;;;;;;;;;;;;;;;

undef("optionsOrtho")
function optionsOrtho(res,center,grid,spacing,gridcol)
;============================

; res    : modified ressource for the plot
; center : lon,lat coordinates for center of projection
; grid   : logical for printing grid
; spacing: lon,lat spacing for grid
; gridcol: color for grid

begin

  res@mpCenterLonF      = center(0)  ; choose center lon
  res@mpCenterLatF      = center(1)  ; choose center lat
  res@mpGridAndLimbOn   = grid
  res@mpGridLonSpacingF = spacing(0)
  res@mpGridLatSpacingF = spacing(1)
  res@mpGridLineColor   = gridcol

  return(res)

end

;;;;;;;;;;;; optionsStereo ;;;;;;;;;;;;;;;;

undef("optionsStereo")
function optionsStereo(res,hm,grid,spacing,gridcol)
;============================

; res    : modified ressource for the plot
; hm     : hemisphere to be plotted
; grid   : logical for printing grid
; spacing: lon,lat spacing for grid
; gridcol: color for grid

begin

  res@gsnPolar          = hm

  res@mpGridAndLimbOn   = grid
  res@mpGridLonSpacingF = spacing(0)
  res@mpGridLatSpacingF = spacing(1)
  res@mpGridLineColor   = gridcol

  return(res)

end

;;;;;;;;;;;; optionsOverPlot_2d ;;;;;;;;;;;;;;;;

undef("optionsOverPlot_2d")
function optionsOverPlot_2d(res)
;============================

; res    : modified ressource for the plot

begin

  res@cnLineLabelFormat            = "@5.1f"
  res@cnLineLabelsOn               = True
  res@cnLabelMasking               = True		
  res@cnLineLabelBackgroundColor   = -1
  res@cnInfoLabelOn                = False
  res@gsnContourNegLineDashPattern = 1 

  return(res)

end

;;;;;;;;;;;; axePress_2d ;;;;;;;;;;;;;;;;

undef("axePress_2d")
function axePress_2d(res)
;============================

; res : modified ressource for the plot

begin

  res@tiYAxisString  = "Pressure"
  res@trYReverse     = True
  res@gsnYAxisIrregular2Log = True

  return(res)

end

;;;;;;;;;;;; axePress_1d ;;;;;;;;;;;;;;;;

undef("axePress_1d")
function axePress_1d(res,axe,limits)
;============================

; res    : modified ressource for the plot
; axe    : Y axis for plot (pressure)
; limits : min and max values for variable axis

begin

  res@tiYAxisString  = "Pressure"
  res@trYReverse     = True
  res@xyYStyle       = "Log"
  res@trYMinF         = axe(0)
  res@trYMaxF         = axe(dimsizes(axe)-1)
  if (limits(0).ne.-888) then
    res@trXMinF       = limits(0)
    res@trXMaxF       = limits(1)
  end if

  return(res)

end


undef("axeNoPress_1d")
function axeNoPress_1d(res,axe,limits)
;============================

; res    : modified ressource for the plot
; axe    : X axis for plot (not pressure)
; limits : min and max values for variable axis

begin

  res@trXMinF         = axe(0)
  res@trXMaxF         = axe(dimsizes(axe)-1)
  if (limits(0).ne.-888) then
    res@trYMinF       = limits(0)
    res@trYMaxF       = limits(1)
  end if

  return(res)

end

;;;;;;;;;;;; axeTimeX_2d ;;;;;;;;;;;;;;;;

undef("axeTimeX_2d")
function axeTimeX_2d(res,axe)
;============================

; res : modified ressource for the plot
; axe : X axis is time

begin

  res@tiXAxisString        = "Time (local days)"
  res@sfXArray             = axe
  res@tmXBLabelFontHeightF = 0.015
  res@trXMinF              = axe(0)
  res@trXMaxF              = axe(dimsizes(axe)-1)

  return(res)

end

;;;;;;;;;;;; axeTimeXls_2d ;;;;;;;;;;;;;;;;

undef("axeTimeXls_2d")
function axeTimeXls_2d(res,axe)
;============================

; res : modified ressource for the plot
; axe : X axis is Ls

begin

  res@tiXAxisString            = "Solar longitude"
  res@sfXArray                 = axe
  res@gsnXAxisIrregular2Linear = True
  res@tmXBMode                 = "Manual"	
  res@tmXBTickStartF           = 0
  res@tmXBTickEndF             = 360
  res@tmXBTickSpacingF         = 30
  res@tmXBMinorPerMajor        = 2
  res@tmXBLabelFontHeightF     = 0.015
  res@trXMinF                  = axe(0)
  res@trXMaxF                  = axe(dimsizes(axe)-1)

  return(res)

end

;;;;;;;;;;;; axeTimeY_2d ;;;;;;;;;;;;;;;;

undef("axeTimeY_2d")
function axeTimeY_2d(res,axe)
;============================

; res : modified ressource for the plot
; axe : Y axis is time

begin

  res@tiYAxisString        = "Time (local days)"
  res@sfYArray             = axe
  res@tmYLLabelFontHeightF = 0.015
  res@trYMinF              = axe(0)
  res@trYMaxF              = axe(dimsizes(axe)-1)

  return(res)

end

;;;;;;;;;;;; axeTimeYls_2d ;;;;;;;;;;;;;;;;

undef("axeTimeYls_2d")
function axeTimeYls_2d(res,axe)
;============================

; res : modified ressource for the plot
; axe : Y axis is Ls

begin

  res@tiYAxisString            = "Solar longitude"
  res@sfYArray                 = axe
  res@gsnYAxisIrregular2Linear = True
  res@tmYLMode                 = "Manual"	
  res@tmYLTickStartF           = 0
  res@tmYLTickEndF             = 360
  res@tmYLTickSpacingF         = 30
  res@tmYLMinorPerMajor        = 2
  res@tmYLLabelFontHeightF     = 0.015
  res@trYMinF                  = axe(0)
  res@trYMaxF                  = axe(dimsizes(axe)-1)

  return(res)

end

;;;;;;;;;;;; axeTime_1d ;;;;;;;;;;;;;;;;

undef("axeTime_1d")
function axeTime_1d(res,axe)
;============================

; res : modified ressource for the plot
; axe : X axis is time

begin

  res@tiXAxisString        = "Time (local days)"
  res@tmXBLabelFontHeightF = 0.015

  return(res)

end

;;;;;;;;;;;; axeTimels_1d ;;;;;;;;;;;;;;;;

undef("axeTimels_1d")
function axeTimels_1d(res,axe)
;============================

; res : modified ressource for the plot
; axe : X axis is Ls

begin

  res@tiXAxisString            = "Solar longitude"
  res@tmXBMode                 = "Manual"	
  res@tmXBTickStartF           = 0
  res@tmXBTickEndF             = 360
  res@tmXBTickSpacingF         = 30
  res@tmXBMinorPerMajor        = 2
  res@tmXBLabelFontHeightF     = 0.015

  return(res)

end

;;;;;;;;;;;; inversDim_2d ;;;;;;;;;;;;;;;;

undef("inversDim_2d")
function inversDim_2d(var)
;============================

; var is a 2d variable

begin

  dimname=getvardims(var)
  tmpvar = var($dimname(1)$|:,$dimname(0)$|:)

  return(tmpvar)

end

;;;;;;;;;;;; revertLon_2dmap ;;;;;;;;;;;;;;;;

undef("revertLon_2dmap")
function revertLon_2dmap(var)
;============================

; var is a 2d lat-lon variable

begin

  tmpvar  = var
  dimname = getvardims(var)
  Lname   = dimname(1)
  nbL     = dimsizes(var&$Lname$)
  tmpaxL  = new(nbL,float)

;revert lon axis:
  do i=0,nbL-1
    tmpaxL(i)=var&$Lname$(nbL-1-i)
  end do
  tmpvar&$Lname$ = tmpaxL

;revert values along lon axis:
  do i=0,nbL-1
    tmpvar(:,i) = var(:,nbL-1-i)
  end do

  return(tmpvar)

end

;;;;;;;;;;;; Map_2d ;;;;;;;;;;;;;;;;

undef("Map_2d")
function Map_2d(var,wks,res)
;============================

; var : 2D variable
; wks : workstation
; res : ressource

begin

; main plot characteristics

  res@gsnSpreadColors = True
  res@tmXTOn          = False
  res@cnFillOn        = True
  res@cnLineLabelsOn  = True
  res@cnLabelMasking  = True
  res@cnLineLabelBackgroundColor = "white"
  res@lbOrientation   = "vertical"

; Create plot

  plotvar  = gsn_csm_contour(wks, var, res )

  return(plotvar)

end

;;;;;;;;;;;; OPlot_2d ;;;;;;;;;;;;;;;;

undef("OPlot_2d")
function OPlot_2d(var,wks,res)
;============================

; var : 2D variable
; wks : workstation
; res : ressource

begin

; main plot characteristics

  res@tmXTOn          = False
  res@tmYROn          = False
  res@cnFillOn        = False

  res@gsnLeftString   = ""
  res@gsnRightString  = ""

; Create plot

  plotvar  = gsn_csm_contour(wks, var, res )

  return(plotvar)

end

;;;;;;;;;;;; Plot_1d ;;;;;;;;;;;;;;;;

undef("Plot_1d")
function Plot_1d(var,wks,res,axe,axep,limits)
;============================

; var            : 1D variable
; wks            : workstation
; res            : ressource
; axe            : axis for plot
; axep           : logical indication for pressure axis
; limits         : min and max values for variable axis

begin

; main plot characteristics

; Create plot
; pressure axis ?
  if (axep) then
    res=axePress_1d(res,axe,limits)
    plotvar  = gsn_csm_xy(wks, var, axe, res )
  else
    res1d=axeNoPress_1d(res,axe,limits)
    plotvar  = gsn_csm_xy(wks, axe, var, res )
  end if

  return(plotvar)

end

;;;;;;;;;;;; Ortho_2d ;;;;;;;;;;;;;;;;

undef("Ortho_2d")
function Ortho_2d(var,wks,res)
;============================

; var : 2D variable
; wks : workstation
; res : ressource

begin

; main plot characteristics

  res@gsnSpreadColors = True
  res@tmXTOn          = False
  res@cnFillOn        = True
  res@cnLineLabelsOn  = True
  res@cnLabelMasking  = True
  res@cnLineLabelBackgroundColor = "white"
  res@lbOrientation   = "vertical"

  res@mpProjection      = "Orthographic"       ; choose projection
  res@mpOutlineOn       = False
  res@mpPerimOn         = False             ; turn off box around plot
  res@mpFillOn          = False
  res@mpGridLineDashPattern = 2
  
; Create plot

  plotvar  = gsn_csm_contour_map(wks, var, res )

  return(plotvar)

end

;;;;;;;;;;;; Stereo_2d ;;;;;;;;;;;;;;;;

undef("Stereo_2d")
function Stereo_2d(var,wks,res)
;============================

; var : 2D variable
; wks : workstation
; res : ressource

begin

; main plot characteristics

  res@gsnSpreadColors = True
  res@tmXTOn          = False
  res@cnFillOn        = True
  res@cnLineLabelsOn  = True
  res@cnLabelMasking  = True
  res@cnLineLabelBackgroundColor = "white"
  res@lbOrientation   = "vertical"

  res@mpProjection      = "Stereographic"       ; choose projection
  res@mpOutlineOn       = False
  res@mpPerimOn         = False             ; turn off box around plot
  res@mpFillOn          = False
  
; Create plot

  plotvar  = gsn_csm_contour_map_polar(wks, var, res )

  return(plotvar)

end

;;;;;;;;;;;; customVar ;;;;;;;;;;;;;;;;

undef("customVar")
function customVar(infile,labelvar)
;============================

; infile   : input file
; labelvar : input label for variable

; THIS IS TO BE CUSTOMIZED EACH TIME ACCORDING TO YOUR NEEDS 
; because it is impossible to make it automatic...

begin

prepared = "dummy"    ; DEFAULT

if (labelvar.eq."ustar") then
  if (isfilevar(infile,"vitu")) then
    u = infile->vitu
  else
    u = infile->U
  end if
  myvar = dim_rmvmean_n_Wrap(u,3)
  myvar@long_name="Deviation from zonal-mean zonal wind"
  prepared = labelvar
end if

if (labelvar.eq."uprime") then
  if (isfilevar(infile,"vitu")) then
    u = infile->vitu
  else
    u = infile->U
  end if
  myvar = dim_rmvmean_n_Wrap(u,0)
  myvar@long_name="Deviation from time-mean zonal wind"
  prepared = labelvar
end if

if (labelvar.eq."dTdyn") then
  if (isfilevar(infile,"vitu")) then
    dtdyn = infile->dtdyn
  else
    dtdyn = infile->DTCORE
  end if
  myvar = dtdyn*1e7
  copy_VarMeta(dtdyn,myvar)
  myvar@units    ="K/Vd"
  myvar@long_name="Dynamics heating rate"
  prepared = labelvar
end if

if (labelvar.eq."dTadj") then
  if (isfilevar(infile,"vitu")) then
    dtadj = infile->dtajs
  else
    dtadj = infile->DADJ
  end if
  myvar = dtadj*1e7
  copy_VarMeta(dtadj,myvar)
  myvar@units    ="K/Vd"
  myvar@long_name="Dry convection heating rate"
  prepared = labelvar
end if

if (labelvar.eq."dTvdf") then
  if (isfilevar(infile,"vitu")) then
    dtvdf = infile->dtvdf
  else
    dtvdf = infile->DTV
  end if
  myvar = dtvdf*1e7
  copy_VarMeta(dtvdf,myvar)
  myvar@units    ="K/Vd"
  myvar@long_name="PBL heating rate"
  prepared = labelvar
end if

if (labelvar.eq."dTrad") then
  if (isfilevar(infile,"vitu")) then
    dtlwr = infile->dtlwr
    dtswr = infile->dtswr
  else
    dtlwr = infile->QRL
    dtswr = infile->QRS
  end if
  myvar = (dtlwr+dtswr)*1e7
  copy_VarMeta(dtlwr,myvar)
  myvar@units    ="K/Vd"
  myvar@long_name="Radiative balance"
  prepared = labelvar
end if

if (labelvar.eq."Tbal") then
  if (isfilevar(infile,"vitu")) then
    dtlwr = infile->dtlwr
    dtswr = infile->dtswr
    dtajs = infile->dtajs
    dtdyn = infile->dtdyn
    dtvdf = infile->dtvdf
  else
    dtlwr = infile->QRL
    dtswr = infile->QRS
    dtajs = infile->DADJ
    dtdyn = infile->DTCORE
    dtvdf = infile->DTV
  end if
  myvar = (dtlwr+dtswr+dtajs+dtdyn+dtvdf)*1e7
  copy_VarMeta(dtlwr,myvar)
  myvar@units    ="K/Vd"
  myvar@long_name="Thermal balance"
  prepared = labelvar
end if

if (labelvar.eq."Tforc") then
  if (isfilevar(infile,"vitu")) then
    dtlwr = infile->dtlwr
    dtswr = infile->dtswr
    dtajs = infile->dtajs
    if (isfilevar(infile,"dtvdf")) then
          dtvdf = infile->dtvdf
    else
          dtvdf = dtajs*0.
    end if
  else
    dtlwr = infile->QRL
    dtswr = infile->QRS
    dtajs = infile->DADJ
    if (isfilevar(infile,"DTV")) then
          dtvdf = infile->DTV
    else
          dtvdf = dtajs*0.
    end if
  end if
  myvar = (dtlwr+dtswr+dtajs+dtvdf)*1e7
  copy_VarMeta(dtlwr,myvar)
  myvar@units    ="K/Vd"
  myvar@long_name="Thermal physical forcing"
  prepared = labelvar
end if

if (labelvar.eq."ske") then
  if (isfilevar(infile,"vitu")) then
    u = infile->vitu
    v = infile->vitv
  else
    u = infile->U
    v = infile->V
  end if
  myvar = u*u+v*v
  copy_VarMeta(u,myvar)
  myvar@units    ="m2/s2"
  myvar@long_name="Specific kinetic energy"
  prepared = labelvar
end if

if (prepared.eq."dummy") then   ; DEFAULT: LEAVE AS IS
  print("You chose a customized variable !")
  print("Modify the function customVar in visu-utils.ncl to prepare this variable")
  exit
end if

  return(myvar)

end

