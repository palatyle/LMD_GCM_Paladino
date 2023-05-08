# LMD_GCM_Paladino
Modified version of the Generic LMD GCM (Now known as the Generic PCM) used in Paladino et al. 2023

Based on the main LMD SVN revision: 2289.

## Modified codes
`trunk/LMDZ.COMMON/libf/phystd/volcano.F90` - Fully rewritten to interface with ATHAM

`trunk/LMDZ.COMMON/libf/phystd/writediagfi.F` - Added `callkeys_mod` call

`trunk/LMDZ.COMMON/libf/phystd/tracer_h.F90` - Added volcano tracer variables

`trunk/LMDZ.COMMON/libf/phystd/initracer.F` - Added volcano tracer qualities (density, radius, etc.)

`trunk/LMDZ.COMMON/libf/phystd/callkeys_mod.F90` - Added variables to be read in from `callphys.def` including volcano tracer density (`rho_volc`), volcano name (`volc_name`), keyword (i.e. tharis, no tharsis - for labeling purposes) (`input_key`), atmosphere type (again mainly for labeling) (`atmos_type`), output directory for netcdf file (`output_dir`).

`trunk/LMDZ.COMMON/libf/phystd/inifis_mod.F90` - Read in variables defined in `callkeys_mod.F90`

`trunk/LMDZ.COMMON/libf/phystd/physiq_mod.F90`- Call `volcano.F90`. Read in correct ash profile depending on season. Also use volcano name to find where volcano should be erupting in model as well as which processor is in charge of this zone. 

`trunk/LMDZ.COMMON/libf/phystd/datafile_mod.F90` - Changed datadir variable to reflect where my particular datadir is. 
