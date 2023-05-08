#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Definitions de Macros pour Make
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Repertoires :

GCM     = /home/palatyle/LMD_gen/trunk/LMDZ.COMMON
LIBF    = $(GCM)/libf
LIBO    = $(LIBOGCM)/$(MACHINE)
LOCAL_DIR=
MISC    = $(LIBF)/misc
FILTRE   = filtre
PHYS  = 
LIBPHY = $(LIBO)/libphy$(PHYS).a
LIBDYN_PHY = $(LIBO)/libdynlmdz_phy$(PHYS).a
DIRMAIN=dyn$(DIM)d$(FLAG_PARA)
RM=rm

OPLINK = 

# Les differentes librairies pour l"edition des liens:
L_DYN      = -ldyn$(DIM)d$(FLAG_PARA) -ldyn3d_common
L_FILTRE   = -l$(FILTRE)
L_DYN_PHY = -ldynlmdz_phy$(PHYS) 
L_PHY = -lphy$(PHYS) -lphy_common
L_MISC    = -lmisc
L_ADJNT    =
L_COSP     = -lcosp
L_AERONOMARS = -laeronomars
L_AERONOSTD = -laeronostd

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Option de compilation FORTRAN
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COMPILE = $(F77) $(OPTIM) $(INCLUDE) -c
COMPILE90 = $(F90) $(OPTIM90) $(INCLUDE) -c
COMPTRU90 = $(F90) $(OPTIMTRU90) $(INCLUDE) -c
LINK    = 
AR      = 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Option de compilation C
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COMPILEC = $(CCC) $(OPTIMC) $(INCLUDEC) -c


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Creation des differents executables
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Executables:
# ------------

PROG = code
DYN  = dyn$(DIM)d

main : $(DYN) misc phys phy_common dyn_phy $(OPTION_DEP) 
	cd $(LIBO) ; $(RANLIB) lib*.a ; cd $(GCM) ;\
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/$(DIRMAIN)/$(SOURCE) -o $(PROG).o ; \
	$(LINK) $(PROG).o -L$(LIBO) $(L_DYN) $(L_DYN_PHY) $(L_ADJNT) $(L_COSP) $(L_AERONOMARS) $(L_AERONOSTD) $(L_PHY) $(L_DYN) $(L_DYN3D_COMMON) $(L_MISC) $(L_DYN3D_COMMON) $(L_PHY) $(L_DYN_PHY) $(L_DYN) $(L_FILTRE) $(OPLINK) $(OPTION_LINK) -o $(LOCAL_DIR)/$(PROG).e ; $(RM) $(PROG).o 

phys : $(LIBPHY)
phy_common : $(LIBPHY_COMMON)
dyn3d : $(LIBO)/libdyn$(DIM)d$(FLAG_PARA).a $(LIBO)/lib$(FILTRE).a $(LIBO)/libdyn3d_common.a
dyn2d : $(LIBO)/libdyn2d.a
dyn1d :

dyn_phy: $(LIBDYN_PHYS)


misc : $(LIBO)/libmisc.a

adjnt : $(LIBO)/libadjnt.a

cosp : $(LIBO)/libcosp.a

chimtitan : $(LIBO)/libchimtitan.a

aeronomars : $(LIBO)/libaeronomars.a

aeronostd : $(LIBO)/libaeronostd.a

$(FILTRE)3d : $(LIBO)/lib$(FILTRE).a

$(FILTRE)2d :

$(FILTRE)1d :

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Contenu des differentes bibliotheques
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory misc
#=======================================================================

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(cbrt.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(chfev.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(cray.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(fdump.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(formcoord.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(i1mach.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(ismax.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(ismin.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(j4save.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(juldate.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(minmax2.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(minmax.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(pchdf.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(pchfe.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(pchsp.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(q_sat.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(ran1.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(sort.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(xercnt.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(xerhlt.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(xermsg.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(xerprn.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(xersve.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(xgetua.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(arth_m.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(assert_eq_m.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(assert_m.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(coefpoly_m.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(handle_err_m.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(interpolation.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(ioipsl_errioipsl.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(ioipsl_getincom.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(ioipsl_stringop.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(netcdf95.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(new_unit_m.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(nf95_def_var_m.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(nf95_get_att_m.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(nf95_get_var_m.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(nf95_gw_var_m.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(nf95_put_att_m.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(nf95_put_var_m.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(nrtype.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(pchfe_95_m.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(pchsp_95_m.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(pres2lev_mod.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(regr1_conserv_m.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(regr1_lint_m.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(regr3_lint_m.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(simple.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(slopes_m.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(sort_mod.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(vampir.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(write_field.o)

$(LIBO)/libmisc.a :  $(LIBO)/libmisc.a(wxios.o)

.PRECIOUS	: $(LIBO)/libmisc.a


# Compilation of elements in misc of library libmisc.a

$(LIBO)/libmisc.a(cbrt.o) : $(LIBF)/misc/cbrt.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a cbrt.o ; \
	$(COMPILE) $(LIBF)/misc/cbrt.F ; \
	$(AR) r $(LIBO)/libmisc.a cbrt.o ; $(RM) cbrt.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(chfev.o) : $(LIBF)/misc/chfev.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a chfev.o ; \
	$(COMPILE) $(LIBF)/misc/chfev.F ; \
	$(AR) r $(LIBO)/libmisc.a chfev.o ; $(RM) chfev.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(cray.o) : $(LIBF)/misc/cray.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a cray.o ; \
	$(COMPILE) $(LIBF)/misc/cray.F ; \
	$(AR) r $(LIBO)/libmisc.a cray.o ; $(RM) cray.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(fdump.o) : $(LIBF)/misc/fdump.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a fdump.o ; \
	$(COMPILE) $(LIBF)/misc/fdump.F ; \
	$(AR) r $(LIBO)/libmisc.a fdump.o ; $(RM) fdump.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(formcoord.o) : $(LIBF)/misc/formcoord.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a formcoord.o ; \
	$(COMPILE) $(LIBF)/misc/formcoord.F ; \
	$(AR) r $(LIBO)/libmisc.a formcoord.o ; $(RM) formcoord.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(i1mach.o) : $(LIBF)/misc/i1mach.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a i1mach.o ; \
	$(COMPILE) $(LIBF)/misc/i1mach.F ; \
	$(AR) r $(LIBO)/libmisc.a i1mach.o ; $(RM) i1mach.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(ismax.o) : $(LIBF)/misc/ismax.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a ismax.o ; \
	$(COMPILE) $(LIBF)/misc/ismax.F ; \
	$(AR) r $(LIBO)/libmisc.a ismax.o ; $(RM) ismax.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(ismin.o) : $(LIBF)/misc/ismin.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a ismin.o ; \
	$(COMPILE) $(LIBF)/misc/ismin.F ; \
	$(AR) r $(LIBO)/libmisc.a ismin.o ; $(RM) ismin.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(j4save.o) : $(LIBF)/misc/j4save.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a j4save.o ; \
	$(COMPILE) $(LIBF)/misc/j4save.F ; \
	$(AR) r $(LIBO)/libmisc.a j4save.o ; $(RM) j4save.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(juldate.o) : $(LIBF)/misc/juldate.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a juldate.o ; \
	$(COMPILE) $(LIBF)/misc/juldate.F ; \
	$(AR) r $(LIBO)/libmisc.a juldate.o ; $(RM) juldate.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(minmax2.o) : $(LIBF)/misc/minmax2.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a minmax2.o ; \
	$(COMPILE) $(LIBF)/misc/minmax2.F ; \
	$(AR) r $(LIBO)/libmisc.a minmax2.o ; $(RM) minmax2.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(minmax.o) : $(LIBF)/misc/minmax.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a minmax.o ; \
	$(COMPILE) $(LIBF)/misc/minmax.F ; \
	$(AR) r $(LIBO)/libmisc.a minmax.o ; $(RM) minmax.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(pchdf.o) : $(LIBF)/misc/pchdf.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a pchdf.o ; \
	$(COMPILE) $(LIBF)/misc/pchdf.F ; \
	$(AR) r $(LIBO)/libmisc.a pchdf.o ; $(RM) pchdf.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(pchfe.o) : $(LIBF)/misc/pchfe.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a pchfe.o ; \
	$(COMPILE) $(LIBF)/misc/pchfe.F ; \
	$(AR) r $(LIBO)/libmisc.a pchfe.o ; $(RM) pchfe.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(pchsp.o) : $(LIBF)/misc/pchsp.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a pchsp.o ; \
	$(COMPILE) $(LIBF)/misc/pchsp.F ; \
	$(AR) r $(LIBO)/libmisc.a pchsp.o ; $(RM) pchsp.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(q_sat.o) : $(LIBF)/misc/q_sat.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a q_sat.o ; \
	$(COMPILE) $(LIBF)/misc/q_sat.F ; \
	$(AR) r $(LIBO)/libmisc.a q_sat.o ; $(RM) q_sat.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(ran1.o) : $(LIBF)/misc/ran1.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a ran1.o ; \
	$(COMPILE) $(LIBF)/misc/ran1.F ; \
	$(AR) r $(LIBO)/libmisc.a ran1.o ; $(RM) ran1.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(sort.o) : $(LIBF)/misc/sort.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a sort.o ; \
	$(COMPILE) $(LIBF)/misc/sort.F ; \
	$(AR) r $(LIBO)/libmisc.a sort.o ; $(RM) sort.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(xercnt.o) : $(LIBF)/misc/xercnt.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a xercnt.o ; \
	$(COMPILE) $(LIBF)/misc/xercnt.F ; \
	$(AR) r $(LIBO)/libmisc.a xercnt.o ; $(RM) xercnt.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(xerhlt.o) : $(LIBF)/misc/xerhlt.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a xerhlt.o ; \
	$(COMPILE) $(LIBF)/misc/xerhlt.F ; \
	$(AR) r $(LIBO)/libmisc.a xerhlt.o ; $(RM) xerhlt.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(xermsg.o) : $(LIBF)/misc/xermsg.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a xermsg.o ; \
	$(COMPILE) $(LIBF)/misc/xermsg.F ; \
	$(AR) r $(LIBO)/libmisc.a xermsg.o ; $(RM) xermsg.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(xerprn.o) : $(LIBF)/misc/xerprn.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a xerprn.o ; \
	$(COMPILE) $(LIBF)/misc/xerprn.F ; \
	$(AR) r $(LIBO)/libmisc.a xerprn.o ; $(RM) xerprn.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(xersve.o) : $(LIBF)/misc/xersve.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a xersve.o ; \
	$(COMPILE) $(LIBF)/misc/xersve.F ; \
	$(AR) r $(LIBO)/libmisc.a xersve.o ; $(RM) xersve.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(xgetua.o) : $(LIBF)/misc/xgetua.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a xgetua.o ; \
	$(COMPILE) $(LIBF)/misc/xgetua.F ; \
	$(AR) r $(LIBO)/libmisc.a xgetua.o ; $(RM) xgetua.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(arth_m.o) : $(LIBF)/misc/arth_m.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a arth_m.o ; \
	$(COMPILE90) $(LIBF)/misc/arth_m.F90 ; \
	$(AR) r $(LIBO)/libmisc.a arth_m.o ; $(RM) arth_m.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(assert_eq_m.o) : $(LIBF)/misc/assert_eq_m.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a assert_eq_m.o ; \
	$(COMPILE90) $(LIBF)/misc/assert_eq_m.F90 ; \
	$(AR) r $(LIBO)/libmisc.a assert_eq_m.o ; $(RM) assert_eq_m.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(assert_m.o) : $(LIBF)/misc/assert_m.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a assert_m.o ; \
	$(COMPILE90) $(LIBF)/misc/assert_m.F90 ; \
	$(AR) r $(LIBO)/libmisc.a assert_m.o ; $(RM) assert_m.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(coefpoly_m.o) : $(LIBF)/misc/coefpoly_m.F90 \
$(LIBO)/libmisc.a(nrtype.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a coefpoly_m.o ; \
	$(COMPILE90) $(LIBF)/misc/coefpoly_m.F90 ; \
	$(AR) r $(LIBO)/libmisc.a coefpoly_m.o ; $(RM) coefpoly_m.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(handle_err_m.o) : $(LIBF)/misc/handle_err_m.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a handle_err_m.o ; \
	$(COMPILE90) $(LIBF)/misc/handle_err_m.F90 ; \
	$(AR) r $(LIBO)/libmisc.a handle_err_m.o ; $(RM) handle_err_m.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(interpolation.o) : $(LIBF)/misc/interpolation.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a interpolation.o ; \
	$(COMPILE90) $(LIBF)/misc/interpolation.F90 ; \
	$(AR) r $(LIBO)/libmisc.a interpolation.o ; $(RM) interpolation.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(ioipsl_errioipsl.o) : $(LIBF)/misc/ioipsl_errioipsl.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a ioipsl_errioipsl.o ; \
	$(COMPILE90) $(LIBF)/misc/ioipsl_errioipsl.F90 ; \
	$(AR) r $(LIBO)/libmisc.a ioipsl_errioipsl.o ; $(RM) ioipsl_errioipsl.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(ioipsl_getincom.o) : $(LIBF)/misc/ioipsl_getincom.F90 \
$(LIBO)/libmisc.a(ioipsl_errioipsl.o) \
$(LIBO)/libmisc.a(ioipsl_stringop.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a ioipsl_getincom.o ; \
	$(COMPILE90) $(LIBF)/misc/ioipsl_getincom.F90 ; \
	$(AR) r $(LIBO)/libmisc.a ioipsl_getincom.o ; $(RM) ioipsl_getincom.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(ioipsl_stringop.o) : $(LIBF)/misc/ioipsl_stringop.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a ioipsl_stringop.o ; \
	$(COMPILE90) $(LIBF)/misc/ioipsl_stringop.F90 ; \
	$(AR) r $(LIBO)/libmisc.a ioipsl_stringop.o ; $(RM) ioipsl_stringop.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(netcdf95.o) : $(LIBF)/misc/netcdf95.F90 \
$(LIBO)/libmisc.a(handle_err_m.o) \
$(LIBO)/libmisc.a(nf95_def_var_m.o) \
$(LIBO)/libmisc.a(nf95_get_att_m.o) \
$(LIBO)/libmisc.a(nf95_get_var_m.o) \
$(LIBO)/libmisc.a(nf95_gw_var_m.o) \
$(LIBO)/libmisc.a(nf95_put_att_m.o) \
$(LIBO)/libmisc.a(nf95_put_var_m.o) \
$(LIBO)/libmisc.a(simple.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a netcdf95.o ; \
	$(COMPILE90) $(LIBF)/misc/netcdf95.F90 ; \
	$(AR) r $(LIBO)/libmisc.a netcdf95.o ; $(RM) netcdf95.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(new_unit_m.o) : $(LIBF)/misc/new_unit_m.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a new_unit_m.o ; \
	$(COMPILE90) $(LIBF)/misc/new_unit_m.F90 ; \
	$(AR) r $(LIBO)/libmisc.a new_unit_m.o ; $(RM) new_unit_m.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(nf95_def_var_m.o) : $(LIBF)/misc/nf95_def_var_m.F90 \
$(LIBO)/libmisc.a(handle_err_m.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a nf95_def_var_m.o ; \
	$(COMPILE90) $(LIBF)/misc/nf95_def_var_m.F90 ; \
	$(AR) r $(LIBO)/libmisc.a nf95_def_var_m.o ; $(RM) nf95_def_var_m.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(nf95_get_att_m.o) : $(LIBF)/misc/nf95_get_att_m.F90 \
$(LIBO)/libmisc.a(handle_err_m.o) \
$(LIBO)/libmisc.a(simple.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a nf95_get_att_m.o ; \
	$(COMPILE90) $(LIBF)/misc/nf95_get_att_m.F90 ; \
	$(AR) r $(LIBO)/libmisc.a nf95_get_att_m.o ; $(RM) nf95_get_att_m.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(nf95_get_var_m.o) : $(LIBF)/misc/nf95_get_var_m.F90 \
$(LIBO)/libmisc.a(handle_err_m.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a nf95_get_var_m.o ; \
	$(COMPILE90) $(LIBF)/misc/nf95_get_var_m.F90 ; \
	$(AR) r $(LIBO)/libmisc.a nf95_get_var_m.o ; $(RM) nf95_get_var_m.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(nf95_gw_var_m.o) : $(LIBF)/misc/nf95_gw_var_m.F90 \
$(LIBO)/libmisc.a(nf95_get_var_m.o) \
$(LIBO)/libmisc.a(simple.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a nf95_gw_var_m.o ; \
	$(COMPILE90) $(LIBF)/misc/nf95_gw_var_m.F90 ; \
	$(AR) r $(LIBO)/libmisc.a nf95_gw_var_m.o ; $(RM) nf95_gw_var_m.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(nf95_put_att_m.o) : $(LIBF)/misc/nf95_put_att_m.F90 \
$(LIBO)/libmisc.a(handle_err_m.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a nf95_put_att_m.o ; \
	$(COMPILE90) $(LIBF)/misc/nf95_put_att_m.F90 ; \
	$(AR) r $(LIBO)/libmisc.a nf95_put_att_m.o ; $(RM) nf95_put_att_m.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(nf95_put_var_m.o) : $(LIBF)/misc/nf95_put_var_m.F90 \
$(LIBO)/libmisc.a(handle_err_m.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a nf95_put_var_m.o ; \
	$(COMPILE90) $(LIBF)/misc/nf95_put_var_m.F90 ; \
	$(AR) r $(LIBO)/libmisc.a nf95_put_var_m.o ; $(RM) nf95_put_var_m.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(nrtype.o) : $(LIBF)/misc/nrtype.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a nrtype.o ; \
	$(COMPILE90) $(LIBF)/misc/nrtype.F90 ; \
	$(AR) r $(LIBO)/libmisc.a nrtype.o ; $(RM) nrtype.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(pchfe_95_m.o) : $(LIBF)/misc/pchfe_95_m.F90 \
$(LIBO)/libmisc.a(assert_eq_m.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a pchfe_95_m.o ; \
	$(COMPILE90) $(LIBF)/misc/pchfe_95_m.F90 ; \
	$(AR) r $(LIBO)/libmisc.a pchfe_95_m.o ; $(RM) pchfe_95_m.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(pchsp_95_m.o) : $(LIBF)/misc/pchsp_95_m.F90 \
$(LIBO)/libmisc.a(assert_eq_m.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a pchsp_95_m.o ; \
	$(COMPILE90) $(LIBF)/misc/pchsp_95_m.F90 ; \
	$(AR) r $(LIBO)/libmisc.a pchsp_95_m.o ; $(RM) pchsp_95_m.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(pres2lev_mod.o) : $(LIBF)/misc/pres2lev_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a pres2lev_mod.o ; \
	$(COMPILE90) $(LIBF)/misc/pres2lev_mod.F90 ; \
	$(AR) r $(LIBO)/libmisc.a pres2lev_mod.o ; $(RM) pres2lev_mod.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(regr1_conserv_m.o) : $(LIBF)/misc/regr1_conserv_m.F90 \
$(LIBO)/libmisc.a(assert_eq_m.o) \
$(LIBO)/libmisc.a(assert_m.o) \
$(LIBO)/libmisc.a(interpolation.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a regr1_conserv_m.o ; \
	$(COMPILE90) $(LIBF)/misc/regr1_conserv_m.F90 ; \
	$(AR) r $(LIBO)/libmisc.a regr1_conserv_m.o ; $(RM) regr1_conserv_m.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(regr1_lint_m.o) : $(LIBF)/misc/regr1_lint_m.F90 \
$(LIBO)/libmisc.a(assert_eq_m.o) \
$(LIBO)/libmisc.a(interpolation.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a regr1_lint_m.o ; \
	$(COMPILE90) $(LIBF)/misc/regr1_lint_m.F90 ; \
	$(AR) r $(LIBO)/libmisc.a regr1_lint_m.o ; $(RM) regr1_lint_m.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(regr3_lint_m.o) : $(LIBF)/misc/regr3_lint_m.F90 \
$(LIBO)/libmisc.a(assert_eq_m.o) \
$(LIBO)/libmisc.a(interpolation.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a regr3_lint_m.o ; \
	$(COMPILE90) $(LIBF)/misc/regr3_lint_m.F90 ; \
	$(AR) r $(LIBO)/libmisc.a regr3_lint_m.o ; $(RM) regr3_lint_m.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(simple.o) : $(LIBF)/misc/simple.F90 \
$(LIBO)/libmisc.a(handle_err_m.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a simple.o ; \
	$(COMPILE90) $(LIBF)/misc/simple.F90 ; \
	$(AR) r $(LIBO)/libmisc.a simple.o ; $(RM) simple.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(slopes_m.o) : $(LIBF)/misc/slopes_m.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a slopes_m.o ; \
	$(COMPILE90) $(LIBF)/misc/slopes_m.F90 ; \
	$(AR) r $(LIBO)/libmisc.a slopes_m.o ; $(RM) slopes_m.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(sort_mod.o) : $(LIBF)/misc/sort_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a sort_mod.o ; \
	$(COMPILE90) $(LIBF)/misc/sort_mod.F90 ; \
	$(AR) r $(LIBO)/libmisc.a sort_mod.o ; $(RM) sort_mod.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(vampir.o) : $(LIBF)/misc/vampir.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a vampir.o ; \
	$(COMPILE90) $(LIBF)/misc/vampir.F90 ; \
	$(AR) r $(LIBO)/libmisc.a vampir.o ; $(RM) vampir.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(write_field.o) : $(LIBF)/misc/write_field.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a write_field.o ; \
	$(COMPILE90) $(LIBF)/misc/write_field.F90 ; \
	$(AR) r $(LIBO)/libmisc.a write_field.o ; $(RM) write_field.o ; \
	cd $(GCM)

$(LIBO)/libmisc.a(wxios.o) : $(LIBF)/misc/wxios.F90 \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphy_common.a(geometry_mod.o) \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_transfert_para.o) \
$(LIBO)/libmisc.a(nrtype.o) \
$(LIBO)/libphy_common.a(print_control_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libmisc.a wxios.o ; \
	$(COMPILE90) $(LIBF)/misc/wxios.F90 ; \
	$(AR) r $(LIBO)/libmisc.a wxios.o ; $(RM) wxios.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory phy_common
#=======================================================================

$(LIBO)/libphy_common.a :  $(LIBO)/libphy_common.a(abort_physic.o)

$(LIBO)/libphy_common.a :  $(LIBO)/libphy_common.a(geometry_mod.o)

$(LIBO)/libphy_common.a :  $(LIBO)/libphy_common.a(init_print_control_mod.o)

$(LIBO)/libphy_common.a :  $(LIBO)/libphy_common.a(ioipsl_getin_p_mod.o)

$(LIBO)/libphy_common.a :  $(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o)

$(LIBO)/libphy_common.a :  $(LIBO)/libphy_common.a(mod_phys_lmdz_mpi_data.o)

$(LIBO)/libphy_common.a :  $(LIBO)/libphy_common.a(mod_phys_lmdz_mpi_transfert.o)

$(LIBO)/libphy_common.a :  $(LIBO)/libphy_common.a(mod_phys_lmdz_omp_data.o)

$(LIBO)/libphy_common.a :  $(LIBO)/libphy_common.a(mod_phys_lmdz_omp_transfert.o)

$(LIBO)/libphy_common.a :  $(LIBO)/libphy_common.a(mod_phys_lmdz_para.o)

$(LIBO)/libphy_common.a :  $(LIBO)/libphy_common.a(mod_phys_lmdz_transfert_para.o)

$(LIBO)/libphy_common.a :  $(LIBO)/libphy_common.a(physics_distribution_mod.o)

$(LIBO)/libphy_common.a :  $(LIBO)/libphy_common.a(print_control_mod.o)

$(LIBO)/libphy_common.a :  $(LIBO)/libphy_common.a(regular_lonlat_mod.o)

$(LIBO)/libphy_common.a :  $(LIBO)/libphy_common.a(vertical_layers_mod.o)

$(LIBO)/libphy_common.a :  $(LIBO)/libphy_common.a(write_field_phy.o)

$(LIBO)/libphy_common.a :  $(LIBO)/libphy_common.a(xios_writefield.o)

.PRECIOUS	: $(LIBO)/libphy_common.a


# Compilation of elements in phy_common of library libphy_common.a

$(LIBO)/libphy_common.a(abort_physic.o) : $(LIBF)/phy_common/abort_physic.F90 \
$(LIBO)/libmisc.a(ioipsl_getincom.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphy_common.a(print_control_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphy_common.a abort_physic.o ; \
	$(COMPILE90) $(LIBF)/phy_common/abort_physic.F90 ; \
	$(AR) r $(LIBO)/libphy_common.a abort_physic.o ; $(RM) abort_physic.o ; \
	cd $(GCM)

$(LIBO)/libphy_common.a(geometry_mod.o) : $(LIBF)/phy_common/geometry_mod.F90 \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libmisc.a(nrtype.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphy_common.a geometry_mod.o ; \
	$(COMPILE90) $(LIBF)/phy_common/geometry_mod.F90 ; \
	$(AR) r $(LIBO)/libphy_common.a geometry_mod.o ; $(RM) geometry_mod.o ; \
	cd $(GCM)

$(LIBO)/libphy_common.a(init_print_control_mod.o) : $(LIBF)/phy_common/init_print_control_mod.F90 \
$(LIBO)/libphy_common.a(ioipsl_getin_p_mod.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphy_common.a(print_control_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphy_common.a init_print_control_mod.o ; \
	$(COMPILE90) $(LIBF)/phy_common/init_print_control_mod.F90 ; \
	$(AR) r $(LIBO)/libphy_common.a init_print_control_mod.o ; $(RM) init_print_control_mod.o ; \
	cd $(GCM)

$(LIBO)/libphy_common.a(ioipsl_getin_p_mod.o) : $(LIBF)/phy_common/ioipsl_getin_p_mod.F90 \
$(LIBO)/libmisc.a(ioipsl_getincom.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_transfert_para.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphy_common.a ioipsl_getin_p_mod.o ; \
	$(COMPILE90) $(LIBF)/phy_common/ioipsl_getin_p_mod.F90 ; \
	$(AR) r $(LIBO)/libphy_common.a ioipsl_getin_p_mod.o ; $(RM) ioipsl_getin_p_mod.o ; \
	cd $(GCM)

$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) : $(LIBF)/phy_common/mod_grid_phy_lmdz.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphy_common.a mod_grid_phy_lmdz.o ; \
	$(COMPILE90) $(LIBF)/phy_common/mod_grid_phy_lmdz.F90 ; \
	$(AR) r $(LIBO)/libphy_common.a mod_grid_phy_lmdz.o ; $(RM) mod_grid_phy_lmdz.o ; \
	cd $(GCM)

$(LIBO)/libphy_common.a(mod_phys_lmdz_mpi_data.o) : $(LIBF)/phy_common/mod_phys_lmdz_mpi_data.F90 \
$(LIBO)/libphy_common.a(print_control_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphy_common.a mod_phys_lmdz_mpi_data.o ; \
	$(COMPILE90) $(LIBF)/phy_common/mod_phys_lmdz_mpi_data.F90 ; \
	$(AR) r $(LIBO)/libphy_common.a mod_phys_lmdz_mpi_data.o ; $(RM) mod_phys_lmdz_mpi_data.o ; \
	cd $(GCM)

$(LIBO)/libphy_common.a(mod_phys_lmdz_mpi_transfert.o) : $(LIBF)/phy_common/mod_phys_lmdz_mpi_transfert.F90 \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_mpi_data.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphy_common.a mod_phys_lmdz_mpi_transfert.o ; \
	$(COMPILE90) $(LIBF)/phy_common/mod_phys_lmdz_mpi_transfert.F90 ; \
	$(AR) r $(LIBO)/libphy_common.a mod_phys_lmdz_mpi_transfert.o ; $(RM) mod_phys_lmdz_mpi_transfert.o ; \
	cd $(GCM)

$(LIBO)/libphy_common.a(mod_phys_lmdz_omp_data.o) : $(LIBF)/phy_common/mod_phys_lmdz_omp_data.F90 \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphy_common.a(print_control_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphy_common.a mod_phys_lmdz_omp_data.o ; \
	$(COMPILE90) $(LIBF)/phy_common/mod_phys_lmdz_omp_data.F90 ; \
	$(AR) r $(LIBO)/libphy_common.a mod_phys_lmdz_omp_data.o ; $(RM) mod_phys_lmdz_omp_data.o ; \
	cd $(GCM)

$(LIBO)/libphy_common.a(mod_phys_lmdz_omp_transfert.o) : $(LIBF)/phy_common/mod_phys_lmdz_omp_transfert.F90 \
$(LIBO)/libphy_common.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_omp_data.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphy_common.a mod_phys_lmdz_omp_transfert.o ; \
	$(COMPILE90) $(LIBF)/phy_common/mod_phys_lmdz_omp_transfert.F90 ; \
	$(AR) r $(LIBO)/libphy_common.a mod_phys_lmdz_omp_transfert.o ; $(RM) mod_phys_lmdz_omp_transfert.o ; \
	cd $(GCM)

$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) : $(LIBF)/phy_common/mod_phys_lmdz_para.F90 \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_transfert_para.o) \
$(LIBO)/libphy_common.a(print_control_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphy_common.a mod_phys_lmdz_para.o ; \
	$(COMPILE90) $(LIBF)/phy_common/mod_phys_lmdz_para.F90 ; \
	$(AR) r $(LIBO)/libphy_common.a mod_phys_lmdz_para.o ; $(RM) mod_phys_lmdz_para.o ; \
	cd $(GCM)

$(LIBO)/libphy_common.a(mod_phys_lmdz_transfert_para.o) : $(LIBF)/phy_common/mod_phys_lmdz_transfert_para.F90 \
$(LIBO)/libphy_common.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_mpi_transfert.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_omp_transfert.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphy_common.a mod_phys_lmdz_transfert_para.o ; \
	$(COMPILE90) $(LIBF)/phy_common/mod_phys_lmdz_transfert_para.F90 ; \
	$(AR) r $(LIBO)/libphy_common.a mod_phys_lmdz_transfert_para.o ; $(RM) mod_phys_lmdz_transfert_para.o ; \
	cd $(GCM)

$(LIBO)/libphy_common.a(physics_distribution_mod.o) : $(LIBF)/phy_common/physics_distribution_mod.F90 \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphy_common.a physics_distribution_mod.o ; \
	$(COMPILE90) $(LIBF)/phy_common/physics_distribution_mod.F90 ; \
	$(AR) r $(LIBO)/libphy_common.a physics_distribution_mod.o ; $(RM) physics_distribution_mod.o ; \
	cd $(GCM)

$(LIBO)/libphy_common.a(print_control_mod.o) : $(LIBF)/phy_common/print_control_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphy_common.a print_control_mod.o ; \
	$(COMPILE90) $(LIBF)/phy_common/print_control_mod.F90 ; \
	$(AR) r $(LIBO)/libphy_common.a print_control_mod.o ; $(RM) print_control_mod.o ; \
	cd $(GCM)

$(LIBO)/libphy_common.a(regular_lonlat_mod.o) : $(LIBF)/phy_common/regular_lonlat_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphy_common.a regular_lonlat_mod.o ; \
	$(COMPILE90) $(LIBF)/phy_common/regular_lonlat_mod.F90 ; \
	$(AR) r $(LIBO)/libphy_common.a regular_lonlat_mod.o ; $(RM) regular_lonlat_mod.o ; \
	cd $(GCM)

$(LIBO)/libphy_common.a(vertical_layers_mod.o) : $(LIBF)/phy_common/vertical_layers_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphy_common.a vertical_layers_mod.o ; \
	$(COMPILE90) $(LIBF)/phy_common/vertical_layers_mod.F90 ; \
	$(AR) r $(LIBO)/libphy_common.a vertical_layers_mod.o ; $(RM) vertical_layers_mod.o ; \
	cd $(GCM)

$(LIBO)/libphy_common.a(write_field_phy.o) : $(LIBF)/phy_common/write_field_phy.F90 \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libmisc.a(write_field.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphy_common.a write_field_phy.o ; \
	$(COMPILE90) $(LIBF)/phy_common/write_field_phy.F90 ; \
	$(AR) r $(LIBO)/libphy_common.a write_field_phy.o ; $(RM) write_field_phy.o ; \
	cd $(GCM)

$(LIBO)/libphy_common.a(xios_writefield.o) : $(LIBF)/phy_common/xios_writefield.F90 \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphy_common.a(print_control_mod.o) \
$(LIBO)/libmisc.a(wxios.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphy_common.a xios_writefield.o ; \
	$(COMPILE90) $(LIBF)/phy_common/xios_writefield.F90 ; \
	$(AR) r $(LIBO)/libphy_common.a xios_writefield.o ; $(RM) xios_writefield.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory phystd
#=======================================================================

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(aerave.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(aerave_new.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(blackl.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(callsedim.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(convadj.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(cvmgp.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(cvmgt.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(dsolver.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(dtridgl.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(evap.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(forceWCfn.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(gfluxi.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(gfluxv.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(iniaerosol.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(iniorbit.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(inistats.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(initracer.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(iniwrite.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(iniwrite_specIR.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(iniwrite_specVI.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(interp_line.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(lagrange.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(mucorr.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(mufract.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(multipl.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(newsedim.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(orbite.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(profile.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(scopyi.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(sfluxi.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(sfluxv.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(soil.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(soil_settings.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(stelang.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(stellarlong.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(surface_nature.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(surfini.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(tpindex.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(vdif_cd.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(vdifc.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(vdif_kc.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(vlz_fi.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(writediagfi.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(writediagspecIR.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(writediagspecVI.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(writeg1d.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(aeropacity.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(aeroptproperties.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(aerosol_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(albedo_calcv.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(albedo_snow_calc.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(ave_stelspec.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(bilinearbig.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(bilinear.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(calc_cpp3d.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(calc_cpp_mugaz.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(calc_rayleigh.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(callcorrk.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(callkeys_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(call_rings.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(comcstfi_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(comdiurn_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(comg1d_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(comgeomfi_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(comm_wrf.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(comsaison_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(comsoil_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(conc_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(condense_co2.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(cpdet_phy_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(cp_neutral.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(datafile_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(def_var.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(dimphy.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(gases_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(hydrol.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(inifis_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(iniwritesoil.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(interpolateH2H2.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(interpolateH2He.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(interpolateH2Ocont_CKD.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(interpolateH2Ocont_PPC.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(interpolateN2H2.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(interpolateN2N2.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(iostart.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(largescale.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(mass_redistribution.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(mkstat.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(moistadj.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(newtrelax.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(nonoro_gwd_ran_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(ocean_slab_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(optci.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(optcv.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(params_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(phyetat0_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(phyredem.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(physiq_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(phys_state_var_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(planete_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(planetwide_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(radcommon_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(radii_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(radinc_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(rain.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(rings.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(setspi.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(setspv.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(slab_ice_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(statto_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(stokes.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(suaer_corrk.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(sugas_corrk.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(su_gases.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(surfdat_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(surf_heat_transp_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(tabfi_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(thermcell_closure.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(thermcell_dq.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(thermcell_dv2.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(thermcell_env.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(thermcell_flux.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(thermcell_height.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(thermcell_main.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(thermcell_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(thermcell_plume.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(time_phylmdz_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(totalcloudfrac.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(tracer_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(turbdiff.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(turb_mod.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(watercommon_h.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(writediagsoil.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(wstats.o)

$(LIBO)/libphystd.a :  $(LIBO)/libphystd.a(xios_output_mod.o)

.PRECIOUS	: $(LIBO)/libphystd.a


# Compilation of elements in phystd of library libphystd.a

$(LIBO)/libphystd.a(aerave.o) : $(LIBF)/phystd/aerave.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a aerave.o ; \
	$(COMPILE) $(LIBF)/phystd/aerave.F ; \
	$(AR) r $(LIBO)/libphystd.a aerave.o ; $(RM) aerave.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(aerave_new.o) : $(LIBF)/phystd/aerave_new.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a aerave_new.o ; \
	$(COMPILE) $(LIBF)/phystd/aerave_new.F ; \
	$(AR) r $(LIBO)/libphystd.a aerave_new.o ; $(RM) aerave_new.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(blackl.o) : $(LIBF)/phystd/blackl.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a blackl.o ; \
	$(COMPILE) $(LIBF)/phystd/blackl.F ; \
	$(AR) r $(LIBO)/libphystd.a blackl.o ; $(RM) blackl.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(callsedim.o) : $(LIBF)/phystd/callsedim.F \
$(LIBO)/libphystd.a(aerosol_mod.o) \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(radii_mod.o) \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a callsedim.o ; \
	$(COMPILE90) $(LIBF)/phystd/callsedim.F ; \
	$(AR) r $(LIBO)/libphystd.a callsedim.o ; $(RM) callsedim.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(convadj.o) : $(LIBF)/phystd/convadj.F \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a convadj.o ; \
	$(COMPILE90) $(LIBF)/phystd/convadj.F ; \
	$(AR) r $(LIBO)/libphystd.a convadj.o ; $(RM) convadj.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(cvmgp.o) : $(LIBF)/phystd/cvmgp.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a cvmgp.o ; \
	$(COMPILE) $(LIBF)/phystd/cvmgp.F ; \
	$(AR) r $(LIBO)/libphystd.a cvmgp.o ; $(RM) cvmgp.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(cvmgt.o) : $(LIBF)/phystd/cvmgt.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a cvmgt.o ; \
	$(COMPILE) $(LIBF)/phystd/cvmgt.F ; \
	$(AR) r $(LIBO)/libphystd.a cvmgt.o ; $(RM) cvmgt.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(dsolver.o) : $(LIBF)/phystd/dsolver.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a dsolver.o ; \
	$(COMPILE) $(LIBF)/phystd/dsolver.F ; \
	$(AR) r $(LIBO)/libphystd.a dsolver.o ; $(RM) dsolver.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(dtridgl.o) : $(LIBF)/phystd/dtridgl.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a dtridgl.o ; \
	$(COMPILE) $(LIBF)/phystd/dtridgl.F ; \
	$(AR) r $(LIBO)/libphystd.a dtridgl.o ; $(RM) dtridgl.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(evap.o) : $(LIBF)/phystd/evap.F \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libphystd.a(watercommon_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a evap.o ; \
	$(COMPILE90) $(LIBF)/phystd/evap.F ; \
	$(AR) r $(LIBO)/libphystd.a evap.o ; $(RM) evap.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(forceWCfn.o) : $(LIBF)/phystd/forceWCfn.F \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a forceWCfn.o ; \
	$(COMPILE90) $(LIBF)/phystd/forceWCfn.F ; \
	$(AR) r $(LIBO)/libphystd.a forceWCfn.o ; $(RM) forceWCfn.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(gfluxi.o) : $(LIBF)/phystd/gfluxi.F \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a gfluxi.o ; \
	$(COMPILE90) $(LIBF)/phystd/gfluxi.F ; \
	$(AR) r $(LIBO)/libphystd.a gfluxi.o ; $(RM) gfluxi.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(gfluxv.o) : $(LIBF)/phystd/gfluxv.F \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a gfluxv.o ; \
	$(COMPILE90) $(LIBF)/phystd/gfluxv.F ; \
	$(AR) r $(LIBO)/libphystd.a gfluxv.o ; $(RM) gfluxv.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(iniaerosol.o) : $(LIBF)/phystd/iniaerosol.F \
$(LIBO)/libphystd.a(aerosol_mod.o) \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a iniaerosol.o ; \
	$(COMPILE90) $(LIBF)/phystd/iniaerosol.F ; \
	$(AR) r $(LIBO)/libphystd.a iniaerosol.o ; $(RM) iniaerosol.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(iniorbit.o) : $(LIBF)/phystd/iniorbit.F \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(planete_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a iniorbit.o ; \
	$(COMPILE90) $(LIBF)/phystd/iniorbit.F ; \
	$(AR) r $(LIBO)/libphystd.a iniorbit.o ; $(RM) iniorbit.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(inistats.o) : $(LIBF)/phystd/inistats.F \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libmisc.a(nrtype.o) \
$(LIBO)/libphy_common.a(regular_lonlat_mod.o) \
$(LIBO)/libphystd.a(statto_mod.o) \
$(LIBO)/libphystd.a(time_phylmdz_mod.o) \
$(LIBO)/libphy_common.a(vertical_layers_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a inistats.o ; \
	$(COMPILE90) $(LIBF)/phystd/inistats.F ; \
	$(AR) r $(LIBO)/libphystd.a inistats.o ; $(RM) inistats.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(initracer.o) : $(LIBF)/phystd/initracer.F \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a initracer.o ; \
	$(COMPILE90) $(LIBF)/phystd/initracer.F ; \
	$(AR) r $(LIBO)/libphystd.a initracer.o ; $(RM) initracer.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(iniwrite.o) : $(LIBF)/phystd/iniwrite.F \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphy_common.a(regular_lonlat_mod.o) \
$(LIBO)/libphystd.a(time_phylmdz_mod.o) \
$(LIBO)/libphy_common.a(vertical_layers_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a iniwrite.o ; \
	$(COMPILE90) $(LIBF)/phystd/iniwrite.F ; \
	$(AR) r $(LIBO)/libphystd.a iniwrite.o ; $(RM) iniwrite.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(iniwrite_specIR.o) : $(LIBF)/phystd/iniwrite_specIR.F \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphy_common.a(regular_lonlat_mod.o) \
$(LIBO)/libphystd.a(time_phylmdz_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a iniwrite_specIR.o ; \
	$(COMPILE90) $(LIBF)/phystd/iniwrite_specIR.F ; \
	$(AR) r $(LIBO)/libphystd.a iniwrite_specIR.o ; $(RM) iniwrite_specIR.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(iniwrite_specVI.o) : $(LIBF)/phystd/iniwrite_specVI.F \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphy_common.a(regular_lonlat_mod.o) \
$(LIBO)/libphystd.a(time_phylmdz_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a iniwrite_specVI.o ; \
	$(COMPILE90) $(LIBF)/phystd/iniwrite_specVI.F ; \
	$(AR) r $(LIBO)/libphystd.a iniwrite_specVI.o ; $(RM) iniwrite_specVI.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(interp_line.o) : $(LIBF)/phystd/interp_line.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a interp_line.o ; \
	$(COMPILE) $(LIBF)/phystd/interp_line.F ; \
	$(AR) r $(LIBO)/libphystd.a interp_line.o ; $(RM) interp_line.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(lagrange.o) : $(LIBF)/phystd/lagrange.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a lagrange.o ; \
	$(COMPILE) $(LIBF)/phystd/lagrange.F ; \
	$(AR) r $(LIBO)/libphystd.a lagrange.o ; $(RM) lagrange.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(mucorr.o) : $(LIBF)/phystd/mucorr.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a mucorr.o ; \
	$(COMPILE) $(LIBF)/phystd/mucorr.F ; \
	$(AR) r $(LIBO)/libphystd.a mucorr.o ; $(RM) mucorr.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(mufract.o) : $(LIBF)/phystd/mufract.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a mufract.o ; \
	$(COMPILE) $(LIBF)/phystd/mufract.F ; \
	$(AR) r $(LIBO)/libphystd.a mufract.o ; $(RM) mufract.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(multipl.o) : $(LIBF)/phystd/multipl.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a multipl.o ; \
	$(COMPILE) $(LIBF)/phystd/multipl.F ; \
	$(AR) r $(LIBO)/libphystd.a multipl.o ; $(RM) multipl.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(newsedim.o) : $(LIBF)/phystd/newsedim.F \
$(LIBO)/libphystd.a(comcstfi_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a newsedim.o ; \
	$(COMPILE90) $(LIBF)/phystd/newsedim.F ; \
	$(AR) r $(LIBO)/libphystd.a newsedim.o ; $(RM) newsedim.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(orbite.o) : $(LIBF)/phystd/orbite.F \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(planete_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a orbite.o ; \
	$(COMPILE90) $(LIBF)/phystd/orbite.F ; \
	$(AR) r $(LIBO)/libphystd.a orbite.o ; $(RM) orbite.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(profile.o) : $(LIBF)/phystd/profile.F \
$(LIBO)/libphy_common.a(ioipsl_getin_p_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a profile.o ; \
	$(COMPILE90) $(LIBF)/phystd/profile.F ; \
	$(AR) r $(LIBO)/libphystd.a profile.o ; $(RM) profile.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(scopyi.o) : $(LIBF)/phystd/scopyi.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a scopyi.o ; \
	$(COMPILE) $(LIBF)/phystd/scopyi.F ; \
	$(AR) r $(LIBO)/libphystd.a scopyi.o ; $(RM) scopyi.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(sfluxi.o) : $(LIBF)/phystd/sfluxi.F \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a sfluxi.o ; \
	$(COMPILE90) $(LIBF)/phystd/sfluxi.F ; \
	$(AR) r $(LIBO)/libphystd.a sfluxi.o ; $(RM) sfluxi.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(sfluxv.o) : $(LIBF)/phystd/sfluxv.F \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a sfluxv.o ; \
	$(COMPILE90) $(LIBF)/phystd/sfluxv.F ; \
	$(AR) r $(LIBO)/libphystd.a sfluxv.o ; $(RM) sfluxv.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(soil.o) : $(LIBF)/phystd/soil.F \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphy_common.a(geometry_mod.o) \
$(LIBO)/libphystd.a(planete_mod.o) \
$(LIBO)/libphystd.a(time_phylmdz_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a soil.o ; \
	$(COMPILE90) $(LIBF)/phystd/soil.F ; \
	$(AR) r $(LIBO)/libphystd.a soil.o ; $(RM) soil.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(soil_settings.o) : $(LIBF)/phystd/soil_settings.F \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphystd.a(iostart.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a soil_settings.o ; \
	$(COMPILE90) $(LIBF)/phystd/soil_settings.F ; \
	$(AR) r $(LIBO)/libphystd.a soil_settings.o ; $(RM) soil_settings.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(stelang.o) : $(LIBF)/phystd/stelang.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a stelang.o ; \
	$(COMPILE) $(LIBF)/phystd/stelang.F ; \
	$(AR) r $(LIBO)/libphystd.a stelang.o ; $(RM) stelang.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(stellarlong.o) : $(LIBF)/phystd/stellarlong.F \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(planete_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a stellarlong.o ; \
	$(COMPILE90) $(LIBF)/phystd/stellarlong.F ; \
	$(AR) r $(LIBO)/libphystd.a stellarlong.o ; $(RM) stellarlong.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(surface_nature.o) : $(LIBF)/phystd/surface_nature.F \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphy_common.a(geometry_mod.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a surface_nature.o ; \
	$(COMPILE90) $(LIBF)/phystd/surface_nature.F ; \
	$(AR) r $(LIBO)/libphystd.a surface_nature.o ; $(RM) surface_nature.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(surfini.o) : $(LIBF)/phystd/surfini.F \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(planetwide_mod.o) \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a surfini.o ; \
	$(COMPILE90) $(LIBF)/phystd/surfini.F ; \
	$(AR) r $(LIBO)/libphystd.a surfini.o ; $(RM) surfini.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(tpindex.o) : $(LIBF)/phystd/tpindex.F \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a tpindex.o ; \
	$(COMPILE90) $(LIBF)/phystd/tpindex.F ; \
	$(AR) r $(LIBO)/libphystd.a tpindex.o ; $(RM) tpindex.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(vdif_cd.o) : $(LIBF)/phystd/vdif_cd.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a vdif_cd.o ; \
	$(COMPILE) $(LIBF)/phystd/vdif_cd.F ; \
	$(AR) r $(LIBO)/libphystd.a vdif_cd.o ; $(RM) vdif_cd.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(vdifc.o) : $(LIBF)/phystd/vdifc.F \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libphystd.a(watercommon_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a vdifc.o ; \
	$(COMPILE90) $(LIBF)/phystd/vdifc.F ; \
	$(AR) r $(LIBO)/libphystd.a vdifc.o ; $(RM) vdifc.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(vdif_kc.o) : $(LIBF)/phystd/vdif_kc.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a vdif_kc.o ; \
	$(COMPILE) $(LIBF)/phystd/vdif_kc.F ; \
	$(AR) r $(LIBO)/libphystd.a vdif_kc.o ; $(RM) vdif_kc.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(vlz_fi.o) : $(LIBF)/phystd/vlz_fi.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a vlz_fi.o ; \
	$(COMPILE) $(LIBF)/phystd/vlz_fi.F ; \
	$(AR) r $(LIBO)/libphystd.a vlz_fi.o ; $(RM) vlz_fi.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(writediagfi.o) : $(LIBF)/phystd/writediagfi.F \
$(LIBO)/libphy_common.a(geometry_mod.o) \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(time_phylmdz_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a writediagfi.o ; \
	$(COMPILE90) $(LIBF)/phystd/writediagfi.F ; \
	$(AR) r $(LIBO)/libphystd.a writediagfi.o ; $(RM) writediagfi.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(writediagspecIR.o) : $(LIBF)/phystd/writediagspecIR.F \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphy_common.a(geometry_mod.o) \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(time_phylmdz_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a writediagspecIR.o ; \
	$(COMPILE90) $(LIBF)/phystd/writediagspecIR.F ; \
	$(AR) r $(LIBO)/libphystd.a writediagspecIR.o ; $(RM) writediagspecIR.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(writediagspecVI.o) : $(LIBF)/phystd/writediagspecVI.F \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphy_common.a(geometry_mod.o) \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(time_phylmdz_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a writediagspecVI.o ; \
	$(COMPILE90) $(LIBF)/phystd/writediagspecVI.F ; \
	$(AR) r $(LIBO)/libphystd.a writediagspecVI.o ; $(RM) writediagspecVI.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(writeg1d.o) : $(LIBF)/phystd/writeg1d.F \
$(LIBO)/libphystd.a(comg1d_mod.o) \
$(LIBO)/libphystd.a(time_phylmdz_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a writeg1d.o ; \
	$(COMPILE90) $(LIBF)/phystd/writeg1d.F ; \
	$(AR) r $(LIBO)/libphystd.a writeg1d.o ; $(RM) writeg1d.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(aeropacity.o) : $(LIBF)/phystd/aeropacity.F90 \
$(LIBO)/libphystd.a(aerosol_mod.o) \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphy_common.a(geometry_mod.o) \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a aeropacity.o ; \
	$(COMPILE90) $(LIBF)/phystd/aeropacity.F90 ; \
	$(AR) r $(LIBO)/libphystd.a aeropacity.o ; $(RM) aeropacity.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(aeroptproperties.o) : $(LIBF)/phystd/aeroptproperties.F90 \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a aeroptproperties.o ; \
	$(COMPILE90) $(LIBF)/phystd/aeroptproperties.F90 ; \
	$(AR) r $(LIBO)/libphystd.a aeroptproperties.o ; $(RM) aeroptproperties.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(aerosol_mod.o) : $(LIBF)/phystd/aerosol_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a aerosol_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/aerosol_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a aerosol_mod.o ; $(RM) aerosol_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(albedo_calcv.o) : $(LIBF)/phystd/albedo_calcv.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a albedo_calcv.o ; \
	$(COMPILE90) $(LIBF)/phystd/albedo_calcv.F90 ; \
	$(AR) r $(LIBO)/libphystd.a albedo_calcv.o ; $(RM) albedo_calcv.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(albedo_snow_calc.o) : $(LIBF)/phystd/albedo_snow_calc.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a albedo_snow_calc.o ; \
	$(COMPILE) $(LIBF)/phystd/albedo_snow_calc.F90 ; \
	$(AR) r $(LIBO)/libphystd.a albedo_snow_calc.o ; $(RM) albedo_snow_calc.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(ave_stelspec.o) : $(LIBF)/phystd/ave_stelspec.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(datafile_mod.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a ave_stelspec.o ; \
	$(COMPILE90) $(LIBF)/phystd/ave_stelspec.F90 ; \
	$(AR) r $(LIBO)/libphystd.a ave_stelspec.o ; $(RM) ave_stelspec.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(bilinearbig.o) : $(LIBF)/phystd/bilinearbig.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a bilinearbig.o ; \
	$(COMPILE) $(LIBF)/phystd/bilinearbig.F90 ; \
	$(AR) r $(LIBO)/libphystd.a bilinearbig.o ; $(RM) bilinearbig.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(bilinear.o) : $(LIBF)/phystd/bilinear.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a bilinear.o ; \
	$(COMPILE) $(LIBF)/phystd/bilinear.F90 ; \
	$(AR) r $(LIBO)/libphystd.a bilinear.o ; $(RM) bilinear.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(calc_cpp3d.o) : $(LIBF)/phystd/calc_cpp3d.F90 \
$(LIBO)/libphystd.a(comcstfi_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a calc_cpp3d.o ; \
	$(COMPILE90) $(LIBF)/phystd/calc_cpp3d.F90 ; \
	$(AR) r $(LIBO)/libphystd.a calc_cpp3d.o ; $(RM) calc_cpp3d.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(calc_cpp_mugaz.o) : $(LIBF)/phystd/calc_cpp_mugaz.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(gases_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a calc_cpp_mugaz.o ; \
	$(COMPILE90) $(LIBF)/phystd/calc_cpp_mugaz.F90 ; \
	$(AR) r $(LIBO)/libphystd.a calc_cpp_mugaz.o ; $(RM) calc_cpp_mugaz.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(calc_rayleigh.o) : $(LIBF)/phystd/calc_rayleigh.F90 \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(gases_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a calc_rayleigh.o ; \
	$(COMPILE90) $(LIBF)/phystd/calc_rayleigh.F90 ; \
	$(AR) r $(LIBO)/libphystd.a calc_rayleigh.o ; $(RM) calc_rayleigh.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(callcorrk.o) : $(LIBF)/phystd/callcorrk.F90 \
$(LIBO)/libphystd.a(aerosol_mod.o) \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(datafile_mod.o) \
$(LIBO)/libphystd.a(gases_h.o) \
$(LIBO)/libphy_common.a(ioipsl_getin_p_mod.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radii_mod.o) \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libphystd.a(watercommon_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a callcorrk.o ; \
	$(COMPILE90) $(LIBF)/phystd/callcorrk.F90 ; \
	$(AR) r $(LIBO)/libphystd.a callcorrk.o ; $(RM) callcorrk.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(callkeys_mod.o) : $(LIBF)/phystd/callkeys_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a callkeys_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/callkeys_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a callkeys_mod.o ; $(RM) callkeys_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(call_rings.o) : $(LIBF)/phystd/call_rings.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(comdiurn_h.o) \
$(LIBO)/libphystd.a(comsaison_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a call_rings.o ; \
	$(COMPILE90) $(LIBF)/phystd/call_rings.F90 ; \
	$(AR) r $(LIBO)/libphystd.a call_rings.o ; $(RM) call_rings.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(comcstfi_mod.o) : $(LIBF)/phystd/comcstfi_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a comcstfi_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/comcstfi_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a comcstfi_mod.o ; $(RM) comcstfi_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(comdiurn_h.o) : $(LIBF)/phystd/comdiurn_h.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a comdiurn_h.o ; \
	$(COMPILE90) $(LIBF)/phystd/comdiurn_h.F90 ; \
	$(AR) r $(LIBO)/libphystd.a comdiurn_h.o ; $(RM) comdiurn_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(comg1d_mod.o) : $(LIBF)/phystd/comg1d_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a comg1d_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/comg1d_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a comg1d_mod.o ; $(RM) comg1d_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(comgeomfi_h.o) : $(LIBF)/phystd/comgeomfi_h.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a comgeomfi_h.o ; \
	$(COMPILE90) $(LIBF)/phystd/comgeomfi_h.F90 ; \
	$(AR) r $(LIBO)/libphystd.a comgeomfi_h.o ; $(RM) comgeomfi_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(comm_wrf.o) : $(LIBF)/phystd/comm_wrf.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a comm_wrf.o ; \
	$(COMPILE90) $(LIBF)/phystd/comm_wrf.F90 ; \
	$(AR) r $(LIBO)/libphystd.a comm_wrf.o ; $(RM) comm_wrf.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(comsaison_h.o) : $(LIBF)/phystd/comsaison_h.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a comsaison_h.o ; \
	$(COMPILE90) $(LIBF)/phystd/comsaison_h.F90 ; \
	$(AR) r $(LIBO)/libphystd.a comsaison_h.o ; $(RM) comsaison_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(comsoil_h.o) : $(LIBF)/phystd/comsoil_h.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a comsoil_h.o ; \
	$(COMPILE90) $(LIBF)/phystd/comsoil_h.F90 ; \
	$(AR) r $(LIBO)/libphystd.a comsoil_h.o ; $(RM) comsoil_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(conc_mod.o) : $(LIBF)/phystd/conc_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a conc_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/conc_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a conc_mod.o ; $(RM) conc_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(condense_co2.o) : $(LIBF)/phystd/condense_co2.F90 \
$(LIBO)/libphystd.a(aerosol_mod.o) \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(gases_h.o) \
$(LIBO)/libphy_common.a(geometry_mod.o) \
$(LIBO)/libphystd.a(radii_mod.o) \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a condense_co2.o ; \
	$(COMPILE90) $(LIBF)/phystd/condense_co2.F90 ; \
	$(AR) r $(LIBO)/libphystd.a condense_co2.o ; $(RM) condense_co2.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(cpdet_phy_mod.o) : $(LIBF)/phystd/cpdet_phy_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a cpdet_phy_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/cpdet_phy_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a cpdet_phy_mod.o ; $(RM) cpdet_phy_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(cp_neutral.o) : $(LIBF)/phystd/cp_neutral.F90 \
$(LIBO)/libphystd.a(gases_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a cp_neutral.o ; \
	$(COMPILE90) $(LIBF)/phystd/cp_neutral.F90 ; \
	$(AR) r $(LIBO)/libphystd.a cp_neutral.o ; $(RM) cp_neutral.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(datafile_mod.o) : $(LIBF)/phystd/datafile_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a datafile_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/datafile_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a datafile_mod.o ; $(RM) datafile_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(def_var.o) : $(LIBF)/phystd/def_var.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a def_var.o ; \
	$(COMPILE) $(LIBF)/phystd/def_var.F90 ; \
	$(AR) r $(LIBO)/libphystd.a def_var.o ; $(RM) def_var.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(dimphy.o) : $(LIBF)/phystd/dimphy.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a dimphy.o ; \
	$(COMPILE90) $(LIBF)/phystd/dimphy.F90 ; \
	$(AR) r $(LIBO)/libphystd.a dimphy.o ; $(RM) dimphy.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(gases_h.o) : $(LIBF)/phystd/gases_h.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a gases_h.o ; \
	$(COMPILE90) $(LIBF)/phystd/gases_h.F90 ; \
	$(AR) r $(LIBO)/libphystd.a gases_h.o ; $(RM) gases_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(hydrol.o) : $(LIBF)/phystd/hydrol.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comdiurn_h.o) \
$(LIBO)/libphy_common.a(geometry_mod.o) \
$(LIBO)/libphy_common.a(ioipsl_getin_p_mod.o) \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(slab_ice_h.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libphystd.a(watercommon_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a hydrol.o ; \
	$(COMPILE90) $(LIBF)/phystd/hydrol.F90 ; \
	$(AR) r $(LIBO)/libphystd.a hydrol.o ; $(RM) hydrol.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(inifis_mod.o) : $(LIBF)/phystd/inifis_mod.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(comdiurn_h.o) \
$(LIBO)/libphystd.a(comgeomfi_h.o) \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphystd.a(datafile_mod.o) \
$(LIBO)/libphy_common.a(init_print_control_mod.o) \
$(LIBO)/libphy_common.a(ioipsl_getin_p_mod.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(planete_mod.o) \
$(LIBO)/libphystd.a(planetwide_mod.o) \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(time_phylmdz_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a inifis_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/inifis_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a inifis_mod.o ; $(RM) inifis_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(iniwritesoil.o) : $(LIBF)/phystd/iniwritesoil.F90 \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphy_common.a(regular_lonlat_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a iniwritesoil.o ; \
	$(COMPILE90) $(LIBF)/phystd/iniwritesoil.F90 ; \
	$(AR) r $(LIBO)/libphystd.a iniwritesoil.o ; $(RM) iniwritesoil.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(interpolateH2H2.o) : $(LIBF)/phystd/interpolateH2H2.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(datafile_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a interpolateH2H2.o ; \
	$(COMPILE90) $(LIBF)/phystd/interpolateH2H2.F90 ; \
	$(AR) r $(LIBO)/libphystd.a interpolateH2H2.o ; $(RM) interpolateH2H2.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(interpolateH2He.o) : $(LIBF)/phystd/interpolateH2He.F90 \
$(LIBO)/libphystd.a(datafile_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a interpolateH2He.o ; \
	$(COMPILE90) $(LIBF)/phystd/interpolateH2He.F90 ; \
	$(AR) r $(LIBO)/libphystd.a interpolateH2He.o ; $(RM) interpolateH2He.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(interpolateH2Ocont_CKD.o) : $(LIBF)/phystd/interpolateH2Ocont_CKD.F90 \
$(LIBO)/libphystd.a(datafile_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a interpolateH2Ocont_CKD.o ; \
	$(COMPILE90) $(LIBF)/phystd/interpolateH2Ocont_CKD.F90 ; \
	$(AR) r $(LIBO)/libphystd.a interpolateH2Ocont_CKD.o ; $(RM) interpolateH2Ocont_CKD.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(interpolateH2Ocont_PPC.o) : $(LIBF)/phystd/interpolateH2Ocont_PPC.F90 \
$(LIBO)/libphystd.a(datafile_mod.o) \
$(LIBO)/libphystd.a(watercommon_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a interpolateH2Ocont_PPC.o ; \
	$(COMPILE90) $(LIBF)/phystd/interpolateH2Ocont_PPC.F90 ; \
	$(AR) r $(LIBO)/libphystd.a interpolateH2Ocont_PPC.o ; $(RM) interpolateH2Ocont_PPC.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(interpolateN2H2.o) : $(LIBF)/phystd/interpolateN2H2.F90 \
$(LIBO)/libphystd.a(datafile_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a interpolateN2H2.o ; \
	$(COMPILE90) $(LIBF)/phystd/interpolateN2H2.F90 ; \
	$(AR) r $(LIBO)/libphystd.a interpolateN2H2.o ; $(RM) interpolateN2H2.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(interpolateN2N2.o) : $(LIBF)/phystd/interpolateN2N2.F90 \
$(LIBO)/libphystd.a(datafile_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a interpolateN2N2.o ; \
	$(COMPILE90) $(LIBF)/phystd/interpolateN2N2.F90 ; \
	$(AR) r $(LIBO)/libphystd.a interpolateN2N2.o ; $(RM) interpolateN2N2.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(iostart.o) : $(LIBF)/phystd/iostart.F90 \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(slab_ice_h.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a iostart.o ; \
	$(COMPILE90) $(LIBF)/phystd/iostart.F90 ; \
	$(AR) r $(LIBO)/libphystd.a iostart.o ; $(RM) iostart.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(largescale.o) : $(LIBF)/phystd/largescale.F90 \
$(LIBO)/libphy_common.a(ioipsl_getin_p_mod.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libphystd.a(watercommon_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a largescale.o ; \
	$(COMPILE90) $(LIBF)/phystd/largescale.F90 ; \
	$(AR) r $(LIBO)/libphystd.a largescale.o ; $(RM) largescale.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(mass_redistribution.o) : $(LIBF)/phystd/mass_redistribution.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(planete_mod.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libphystd.a(watercommon_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a mass_redistribution.o ; \
	$(COMPILE90) $(LIBF)/phystd/mass_redistribution.F90 ; \
	$(AR) r $(LIBO)/libphystd.a mass_redistribution.o ; $(RM) mass_redistribution.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(mkstat.o) : $(LIBF)/phystd/mkstat.F90 \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(statto_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a mkstat.o ; \
	$(COMPILE90) $(LIBF)/phystd/mkstat.F90 ; \
	$(AR) r $(LIBO)/libphystd.a mkstat.o ; $(RM) mkstat.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(moistadj.o) : $(LIBF)/phystd/moistadj.F90 \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libphystd.a(watercommon_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a moistadj.o ; \
	$(COMPILE90) $(LIBF)/phystd/moistadj.F90 ; \
	$(AR) r $(LIBO)/libphystd.a moistadj.o ; $(RM) moistadj.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(newtrelax.o) : $(LIBF)/phystd/newtrelax.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a newtrelax.o ; \
	$(COMPILE90) $(LIBF)/phystd/newtrelax.F90 ; \
	$(AR) r $(LIBO)/libphystd.a newtrelax.o ; $(RM) newtrelax.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(nonoro_gwd_ran_mod.o) : $(LIBF)/phystd/nonoro_gwd_ran_mod.F90 \
$(LIBO)/libmisc.a(assert_m.o) \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(inifis_mod.o) \
$(LIBO)/libphy_common.a(ioipsl_getin_p_mod.o) \
$(LIBO)/libphy_common.a(vertical_layers_mod.o) \
$(LIBO)/libphystd.a(xios_output_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a nonoro_gwd_ran_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/nonoro_gwd_ran_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a nonoro_gwd_ran_mod.o ; $(RM) nonoro_gwd_ran_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(ocean_slab_mod.o) : $(LIBF)/phystd/ocean_slab_mod.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(slab_ice_h.o) \
$(LIBO)/libphystd.a(surf_heat_transp_mod.o) \
$(LIBO)/libphystd.a(watercommon_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a ocean_slab_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/ocean_slab_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a ocean_slab_mod.o ; $(RM) ocean_slab_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(optci.o) : $(LIBF)/phystd/optci.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(gases_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a optci.o ; \
	$(COMPILE90) $(LIBF)/phystd/optci.F90 ; \
	$(AR) r $(LIBO)/libphystd.a optci.o ; $(RM) optci.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(optcv.o) : $(LIBF)/phystd/optcv.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(gases_h.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a optcv.o ; \
	$(COMPILE90) $(LIBF)/phystd/optcv.F90 ; \
	$(AR) r $(LIBO)/libphystd.a optcv.o ; $(RM) optcv.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(params_h.o) : $(LIBF)/phystd/params_h.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a params_h.o ; \
	$(COMPILE90) $(LIBF)/phystd/params_h.F90 ; \
	$(AR) r $(LIBO)/libphystd.a params_h.o ; $(RM) params_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(phyetat0_mod.o) : $(LIBF)/phystd/phyetat0_mod.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphy_common.a(ioipsl_getin_p_mod.o) \
$(LIBO)/libphystd.a(iostart.o) \
$(LIBO)/libphystd.a(nonoro_gwd_ran_mod.o) \
$(LIBO)/libphystd.a(slab_ice_h.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(tabfi_mod.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a phyetat0_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/phyetat0_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a phyetat0_mod.o ; $(RM) phyetat0_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(phyredem.o) : $(LIBF)/phystd/phyredem.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphy_common.a(geometry_mod.o) \
$(LIBO)/libphystd.a(iostart.o) \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(nonoro_gwd_ran_mod.o) \
$(LIBO)/libphystd.a(planete_mod.o) \
$(LIBO)/libphystd.a(slab_ice_h.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(time_phylmdz_mod.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a phyredem.o ; \
	$(COMPILE90) $(LIBF)/phystd/phyredem.F90 ; \
	$(AR) r $(LIBO)/libphystd.a phyredem.o ; $(RM) phyredem.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(physiq_mod.o) : $(LIBF)/phystd/physiq_mod.F90 \
$(LIBO)/libphystd.a(aerosol_mod.o) \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(comdiurn_h.o) \
$(LIBO)/libphystd.a(comgeomfi_h.o) \
$(LIBO)/libphystd.a(comm_wrf.o) \
$(LIBO)/libphystd.a(comsaison_h.o) \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphystd.a(conc_mod.o) \
$(LIBO)/libphystd.a(gases_h.o) \
$(LIBO)/libphy_common.a(geometry_mod.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(nonoro_gwd_ran_mod.o) \
$(LIBO)/libphystd.a(ocean_slab_mod.o) \
$(LIBO)/libphystd.a(phyetat0_mod.o) \
$(LIBO)/libphystd.a(phyredem.o) \
$(LIBO)/libphystd.a(phys_state_var_mod.o) \
$(LIBO)/libphystd.a(planete_mod.o) \
$(LIBO)/libphystd.a(planetwide_mod.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radii_mod.o) \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(slab_ice_h.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(surf_heat_transp_mod.o) \
$(LIBO)/libphystd.a(thermcell_mod.o) \
$(LIBO)/libphystd.a(time_phylmdz_mod.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libphystd.a(turb_mod.o) \
$(LIBO)/libphy_common.a(vertical_layers_mod.o) \
$(LIBO)/libphystd.a(watercommon_h.o) \
$(LIBO)/libmisc.a(wxios.o) \
$(LIBO)/libphystd.a(xios_output_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a physiq_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/physiq_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a physiq_mod.o ; $(RM) physiq_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(phys_state_var_mod.o) : $(LIBF)/phystd/phys_state_var_mod.F90 \
$(LIBO)/libphystd.a(comsaison_h.o) \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphystd.a(nonoro_gwd_ran_mod.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(slab_ice_h.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(turb_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a phys_state_var_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/phys_state_var_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a phys_state_var_mod.o ; $(RM) phys_state_var_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(planete_mod.o) : $(LIBF)/phystd/planete_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a planete_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/planete_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a planete_mod.o ; $(RM) planete_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(planetwide_mod.o) : $(LIBF)/phystd/planetwide_mod.F90 \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a planetwide_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/planetwide_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a planetwide_mod.o ; $(RM) planetwide_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(radcommon_h.o) : $(LIBF)/phystd/radcommon_h.F90 \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a radcommon_h.o ; \
	$(COMPILE90) $(LIBF)/phystd/radcommon_h.F90 ; \
	$(AR) r $(LIBO)/libphystd.a radcommon_h.o ; $(RM) radcommon_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(radii_mod.o) : $(LIBF)/phystd/radii_mod.F90 \
$(LIBO)/libphystd.a(aerosol_mod.o) \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphy_common.a(ioipsl_getin_p_mod.o) \
$(LIBO)/libphystd.a(radinc_h.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libphystd.a(watercommon_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a radii_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/radii_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a radii_mod.o ; $(RM) radii_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(radinc_h.o) : $(LIBF)/phystd/radinc_h.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a radinc_h.o ; \
	$(COMPILE90) $(LIBF)/phystd/radinc_h.F90 ; \
	$(AR) r $(LIBO)/libphystd.a radinc_h.o ; $(RM) radinc_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(rain.o) : $(LIBF)/phystd/rain.F90 \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphy_common.a(ioipsl_getin_p_mod.o) \
$(LIBO)/libphystd.a(radii_mod.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libphystd.a(watercommon_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a rain.o ; \
	$(COMPILE90) $(LIBF)/phystd/rain.F90 ; \
	$(AR) r $(LIBO)/libphystd.a rain.o ; $(RM) rain.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(rings.o) : $(LIBF)/phystd/rings.F90 \
$(LIBO)/libphystd.a(comdiurn_h.o) \
$(LIBO)/libphy_common.a(geometry_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a rings.o ; \
	$(COMPILE90) $(LIBF)/phystd/rings.F90 ; \
	$(AR) r $(LIBO)/libphystd.a rings.o ; $(RM) rings.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(setspi.o) : $(LIBF)/phystd/setspi.F90 \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(datafile_mod.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a setspi.o ; \
	$(COMPILE90) $(LIBF)/phystd/setspi.F90 ; \
	$(AR) r $(LIBO)/libphystd.a setspi.o ; $(RM) setspi.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(setspv.o) : $(LIBF)/phystd/setspv.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(datafile_mod.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a setspv.o ; \
	$(COMPILE90) $(LIBF)/phystd/setspv.F90 ; \
	$(AR) r $(LIBO)/libphystd.a setspv.o ; $(RM) setspv.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(slab_ice_h.o) : $(LIBF)/phystd/slab_ice_h.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a slab_ice_h.o ; \
	$(COMPILE90) $(LIBF)/phystd/slab_ice_h.F90 ; \
	$(AR) r $(LIBO)/libphystd.a slab_ice_h.o ; $(RM) slab_ice_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(statto_mod.o) : $(LIBF)/phystd/statto_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a statto_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/statto_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a statto_mod.o ; $(RM) statto_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(stokes.o) : $(LIBF)/phystd/stokes.F90 \
$(LIBO)/libphystd.a(comcstfi_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a stokes.o ; \
	$(COMPILE90) $(LIBF)/phystd/stokes.F90 ; \
	$(AR) r $(LIBO)/libphystd.a stokes.o ; $(RM) stokes.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(suaer_corrk.o) : $(LIBF)/phystd/suaer_corrk.F90 \
$(LIBO)/libphystd.a(aerosol_mod.o) \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(datafile_mod.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a suaer_corrk.o ; \
	$(COMPILE90) $(LIBF)/phystd/suaer_corrk.F90 ; \
	$(AR) r $(LIBO)/libphystd.a suaer_corrk.o ; $(RM) suaer_corrk.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(sugas_corrk.o) : $(LIBF)/phystd/sugas_corrk.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(datafile_mod.o) \
$(LIBO)/libphystd.a(gases_h.o) \
$(LIBO)/libphy_common.a(ioipsl_getin_p_mod.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(radinc_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a sugas_corrk.o ; \
	$(COMPILE90) $(LIBF)/phystd/sugas_corrk.F90 ; \
	$(AR) r $(LIBO)/libphystd.a sugas_corrk.o ; $(RM) sugas_corrk.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(su_gases.o) : $(LIBF)/phystd/su_gases.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a su_gases.o ; \
	$(COMPILE90) $(LIBF)/phystd/su_gases.F90 ; \
	$(AR) r $(LIBO)/libphystd.a su_gases.o ; $(RM) su_gases.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(surfdat_h.o) : $(LIBF)/phystd/surfdat_h.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a surfdat_h.o ; \
	$(COMPILE90) $(LIBF)/phystd/surfdat_h.F90 ; \
	$(AR) r $(LIBO)/libphystd.a surfdat_h.o ; $(RM) surfdat_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(surf_heat_transp_mod.o) : $(LIBF)/phystd/surf_heat_transp_mod.F90 \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(slab_ice_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a surf_heat_transp_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/surf_heat_transp_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a surf_heat_transp_mod.o ; $(RM) surf_heat_transp_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(tabfi_mod.o) : $(LIBF)/phystd/tabfi_mod.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphy_common.a(ioipsl_getin_p_mod.o) \
$(LIBO)/libphystd.a(iostart.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(planete_mod.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(time_phylmdz_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a tabfi_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/tabfi_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a tabfi_mod.o ; $(RM) tabfi_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(thermcell_closure.o) : $(LIBF)/phystd/thermcell_closure.F90 \
$(LIBO)/libphystd.a(thermcell_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a thermcell_closure.o ; \
	$(COMPILE90) $(LIBF)/phystd/thermcell_closure.F90 ; \
	$(AR) r $(LIBO)/libphystd.a thermcell_closure.o ; $(RM) thermcell_closure.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(thermcell_dq.o) : $(LIBF)/phystd/thermcell_dq.F90 \
$(LIBO)/libphy_common.a(print_control_mod.o) \
$(LIBO)/libphystd.a(thermcell_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a thermcell_dq.o ; \
	$(COMPILE90) $(LIBF)/phystd/thermcell_dq.F90 ; \
	$(AR) r $(LIBO)/libphystd.a thermcell_dq.o ; $(RM) thermcell_dq.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(thermcell_dv2.o) : $(LIBF)/phystd/thermcell_dv2.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a thermcell_dv2.o ; \
	$(COMPILE) $(LIBF)/phystd/thermcell_dv2.F90 ; \
	$(AR) r $(LIBO)/libphystd.a thermcell_dv2.o ; $(RM) thermcell_dv2.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(thermcell_env.o) : $(LIBF)/phystd/thermcell_env.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphy_common.a(print_control_mod.o) \
$(LIBO)/libphystd.a(thermcell_mod.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libphystd.a(watercommon_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a thermcell_env.o ; \
	$(COMPILE90) $(LIBF)/phystd/thermcell_env.F90 ; \
	$(AR) r $(LIBO)/libphystd.a thermcell_env.o ; $(RM) thermcell_env.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(thermcell_flux.o) : $(LIBF)/phystd/thermcell_flux.F90 \
$(LIBO)/libphy_common.a(print_control_mod.o) \
$(LIBO)/libphystd.a(thermcell_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a thermcell_flux.o ; \
	$(COMPILE90) $(LIBF)/phystd/thermcell_flux.F90 ; \
	$(AR) r $(LIBO)/libphystd.a thermcell_flux.o ; $(RM) thermcell_flux.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(thermcell_height.o) : $(LIBF)/phystd/thermcell_height.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a thermcell_height.o ; \
	$(COMPILE) $(LIBF)/phystd/thermcell_height.F90 ; \
	$(AR) r $(LIBO)/libphystd.a thermcell_height.o ; $(RM) thermcell_height.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(thermcell_main.o) : $(LIBF)/phystd/thermcell_main.F90 \
$(LIBO)/libphy_common.a(print_control_mod.o) \
$(LIBO)/libphystd.a(thermcell_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a thermcell_main.o ; \
	$(COMPILE90) $(LIBF)/phystd/thermcell_main.F90 ; \
	$(AR) r $(LIBO)/libphystd.a thermcell_main.o ; $(RM) thermcell_main.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(thermcell_mod.o) : $(LIBF)/phystd/thermcell_mod.F90 \
$(LIBO)/libphy_common.a(ioipsl_getin_p_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a thermcell_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/thermcell_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a thermcell_mod.o ; $(RM) thermcell_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(thermcell_plume.o) : $(LIBF)/phystd/thermcell_plume.F90 \
$(LIBO)/libphy_common.a(print_control_mod.o) \
$(LIBO)/libphystd.a(thermcell_mod.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libphystd.a(watercommon_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a thermcell_plume.o ; \
	$(COMPILE90) $(LIBF)/phystd/thermcell_plume.F90 ; \
	$(AR) r $(LIBO)/libphystd.a thermcell_plume.o ; $(RM) thermcell_plume.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(time_phylmdz_mod.o) : $(LIBF)/phystd/time_phylmdz_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a time_phylmdz_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/time_phylmdz_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a time_phylmdz_mod.o ; $(RM) time_phylmdz_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(totalcloudfrac.o) : $(LIBF)/phystd/totalcloudfrac.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comdiurn_h.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libphystd.a(watercommon_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a totalcloudfrac.o ; \
	$(COMPILE90) $(LIBF)/phystd/totalcloudfrac.F90 ; \
	$(AR) r $(LIBO)/libphystd.a totalcloudfrac.o ; $(RM) totalcloudfrac.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(tracer_h.o) : $(LIBF)/phystd/tracer_h.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a tracer_h.o ; \
	$(COMPILE90) $(LIBF)/phystd/tracer_h.F90 ; \
	$(AR) r $(LIBO)/libphystd.a tracer_h.o ; $(RM) tracer_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(turbdiff.o) : $(LIBF)/phystd/turbdiff.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(comm_wrf.o) \
$(LIBO)/libphystd.a(radcommon_h.o) \
$(LIBO)/libphystd.a(surfdat_h.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libphystd.a(turb_mod.o) \
$(LIBO)/libphystd.a(watercommon_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a turbdiff.o ; \
	$(COMPILE90) $(LIBF)/phystd/turbdiff.F90 ; \
	$(AR) r $(LIBO)/libphystd.a turbdiff.o ; $(RM) turbdiff.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(turb_mod.o) : $(LIBF)/phystd/turb_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a turb_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/turb_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a turb_mod.o ; $(RM) turb_mod.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(watercommon_h.o) : $(LIBF)/phystd/watercommon_h.F90 \
$(LIBO)/libphystd.a(comcstfi_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a watercommon_h.o ; \
	$(COMPILE90) $(LIBF)/phystd/watercommon_h.F90 ; \
	$(AR) r $(LIBO)/libphystd.a watercommon_h.o ; $(RM) watercommon_h.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(writediagsoil.o) : $(LIBF)/phystd/writediagsoil.F90 \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphy_common.a(geometry_mod.o) \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(time_phylmdz_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a writediagsoil.o ; \
	$(COMPILE90) $(LIBF)/phystd/writediagsoil.F90 ; \
	$(AR) r $(LIBO)/libphystd.a writediagsoil.o ; $(RM) writediagsoil.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(wstats.o) : $(LIBF)/phystd/wstats.F90 \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(statto_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a wstats.o ; \
	$(COMPILE90) $(LIBF)/phystd/wstats.F90 ; \
	$(AR) r $(LIBO)/libphystd.a wstats.o ; $(RM) wstats.o ; \
	cd $(GCM)

$(LIBO)/libphystd.a(xios_output_mod.o) : $(LIBF)/phystd/xios_output_mod.F90 \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libmisc.a(nrtype.o) \
$(LIBO)/libphy_common.a(print_control_mod.o) \
$(LIBO)/libphy_common.a(regular_lonlat_mod.o) \
$(LIBO)/libmisc.a(wxios.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libphystd.a xios_output_mod.o ; \
	$(COMPILE90) $(LIBF)/phystd/xios_output_mod.F90 ; \
	$(AR) r $(LIBO)/libphystd.a xios_output_mod.o ; $(RM) xios_output_mod.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory aeronostd
#=======================================================================

$(LIBO)/libaeronostd.a :  $(LIBO)/libaeronostd.a(concentrations.o)

$(LIBO)/libaeronostd.a :  $(LIBO)/libaeronostd.a(intrplf.o)

$(LIBO)/libaeronostd.a :  $(LIBO)/libaeronostd.a(inv.o)

$(LIBO)/libaeronostd.a :  $(LIBO)/libaeronostd.a(calchim_asis.o)

$(LIBO)/libaeronostd.a :  $(LIBO)/libaeronostd.a(photochemistry_asis.o)

$(LIBO)/libaeronostd.a :  $(LIBO)/libaeronostd.a(photolysis_asis.o)

$(LIBO)/libaeronostd.a :  $(LIBO)/libaeronostd.a(read_phototable.o)

$(LIBO)/libaeronostd.a :  $(LIBO)/libaeronostd.a(types_asis.o)

.PRECIOUS	: $(LIBO)/libaeronostd.a


# Compilation of elements in aeronostd of library libaeronostd.a

$(LIBO)/libaeronostd.a(concentrations.o) : $(LIBF)/aeronostd/concentrations.F \
$(LIBF)/aeronostd/chimiedata.h \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(conc_mod.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libaeronostd.a concentrations.o ; \
	$(COMPILE90) $(LIBF)/aeronostd/concentrations.F ; \
	$(AR) r $(LIBO)/libaeronostd.a concentrations.o ; $(RM) concentrations.o ; \
	cd $(GCM)

$(LIBO)/libaeronostd.a(intrplf.o) : $(LIBF)/aeronostd/intrplf.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libaeronostd.a intrplf.o ; \
	$(COMPILE) $(LIBF)/aeronostd/intrplf.F ; \
	$(AR) r $(LIBO)/libaeronostd.a intrplf.o ; $(RM) intrplf.o ; \
	cd $(GCM)

$(LIBO)/libaeronostd.a(inv.o) : $(LIBF)/aeronostd/inv.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libaeronostd.a inv.o ; \
	$(COMPILE) $(LIBF)/aeronostd/inv.F ; \
	$(AR) r $(LIBO)/libaeronostd.a inv.o ; $(RM) inv.o ; \
	cd $(GCM)

$(LIBO)/libaeronostd.a(calchim_asis.o) : $(LIBF)/aeronostd/calchim_asis.F90 \
$(LIBF)/aeronostd/chimiedata.h \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(conc_mod.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libaeronostd.a calchim_asis.o ; \
	$(COMPILE90) $(LIBF)/aeronostd/calchim_asis.F90 ; \
	$(AR) r $(LIBO)/libaeronostd.a calchim_asis.o ; $(RM) calchim_asis.o ; \
	cd $(GCM)

$(LIBO)/libaeronostd.a(photochemistry_asis.o) : $(LIBF)/aeronostd/photochemistry_asis.F90 \
$(LIBF)/aeronostd/chimiedata.h \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o) \
$(LIBO)/libphystd.a(tracer_h.o) \
$(LIBO)/libaeronostd.a(types_asis.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libaeronostd.a photochemistry_asis.o ; \
	$(COMPILE90) $(LIBF)/aeronostd/photochemistry_asis.F90 ; \
	$(AR) r $(LIBO)/libaeronostd.a photochemistry_asis.o ; $(RM) photochemistry_asis.o ; \
	cd $(GCM)

$(LIBO)/libaeronostd.a(photolysis_asis.o) : $(LIBF)/aeronostd/photolysis_asis.F90 \
$(LIBF)/aeronostd/chimiedata.h \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libphystd.a(comcstfi_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libaeronostd.a photolysis_asis.o ; \
	$(COMPILE90) $(LIBF)/aeronostd/photolysis_asis.F90 ; \
	$(AR) r $(LIBO)/libaeronostd.a photolysis_asis.o ; $(RM) photolysis_asis.o ; \
	cd $(GCM)

$(LIBO)/libaeronostd.a(read_phototable.o) : $(LIBF)/aeronostd/read_phototable.F90 \
$(LIBF)/aeronostd/chimiedata.h \
$(LIBO)/libphystd.a(datafile_mod.o) \
$(LIBO)/libmisc.a(ioipsl_getincom.o) \
$(LIBO)/libphy_common.a(ioipsl_getin_p_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libaeronostd.a read_phototable.o ; \
	$(COMPILE90) $(LIBF)/aeronostd/read_phototable.F90 ; \
	$(AR) r $(LIBO)/libaeronostd.a read_phototable.o ; $(RM) read_phototable.o ; \
	cd $(GCM)

$(LIBO)/libaeronostd.a(types_asis.o) : $(LIBF)/aeronostd/types_asis.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libaeronostd.a types_asis.o ; \
	$(COMPILE90) $(LIBF)/aeronostd/types_asis.F90 ; \
	$(AR) r $(LIBO)/libaeronostd.a types_asis.o ; $(RM) types_asis.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory filtrez
#=======================================================================

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(acc.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(eigen.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(eigen_sort.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(filtreg.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(inifgn.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(filtreg_mod.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(jacobi.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(mkl_dfti.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(mkl_dft_type.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(mod_fft.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(mod_fft_fftw.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(mod_fft_mathkeisan.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(mod_fft_mkl.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(mod_fft_wrapper.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(mod_filtre_fft.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(mod_filtre_fft_loc.o)

$(LIBO)/libfiltrez.a :  $(LIBO)/libfiltrez.a(timer_filtre.o)

.PRECIOUS	: $(LIBO)/libfiltrez.a


# Compilation of elements in filtrez of library libfiltrez.a

$(LIBO)/libfiltrez.a(acc.o) : $(LIBF)/filtrez/acc.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libfiltrez.a acc.o ; \
	$(COMPILE) $(LIBF)/filtrez/acc.F ; \
	$(AR) r $(LIBO)/libfiltrez.a acc.o ; $(RM) acc.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(eigen.o) : $(LIBF)/filtrez/eigen.F \
$(LIBF)/grid/dimensions.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libfiltrez.a eigen.o ; \
	$(COMPILE) $(LIBF)/filtrez/eigen.F ; \
	$(AR) r $(LIBO)/libfiltrez.a eigen.o ; $(RM) eigen.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(eigen_sort.o) : $(LIBF)/filtrez/eigen_sort.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libfiltrez.a eigen_sort.o ; \
	$(COMPILE) $(LIBF)/filtrez/eigen_sort.F ; \
	$(AR) r $(LIBO)/libfiltrez.a eigen_sort.o ; $(RM) eigen_sort.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(filtreg.o) : $(LIBF)/filtrez/filtreg.F \
$(LIBF)/filtrez/coefils.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libfiltrez.a(filtreg_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libfiltrez.a filtreg.o ; \
	$(COMPILE90) $(LIBF)/filtrez/filtreg.F ; \
	$(AR) r $(LIBO)/libfiltrez.a filtreg.o ; $(RM) filtreg.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(inifgn.o) : $(LIBF)/filtrez/inifgn.F \
$(LIBF)/filtrez/coefils.h \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libfiltrez.a inifgn.o ; \
	$(COMPILE) $(LIBF)/filtrez/inifgn.F ; \
	$(AR) r $(LIBO)/libfiltrez.a inifgn.o ; $(RM) inifgn.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(filtreg_mod.o) : $(LIBF)/filtrez/filtreg_mod.F90 \
$(LIBF)/filtrez/coefils.h \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(logic_mod.o) \
$(LIBO)/libfiltrez.a(mod_filtre_fft.o) \
$(LIBO)/libfiltrez.a(mod_filtre_fft_loc.o) \
$(LIBO)/libdyn3d_common.a(serre_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libfiltrez.a filtreg_mod.o ; \
	$(COMPILE90) $(LIBF)/filtrez/filtreg_mod.F90 ; \
	$(AR) r $(LIBO)/libfiltrez.a filtreg_mod.o ; $(RM) filtreg_mod.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(jacobi.o) : $(LIBF)/filtrez/jacobi.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libfiltrez.a jacobi.o ; \
	$(COMPILE) $(LIBF)/filtrez/jacobi.F90 ; \
	$(AR) r $(LIBO)/libfiltrez.a jacobi.o ; $(RM) jacobi.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(mkl_dfti.o) : $(LIBF)/filtrez/mkl_dfti.f90 \
$(LIBO)/libfiltrez.a(mkl_dft_type.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libfiltrez.a mkl_dfti.o ; \
	$(COMPILE90) $(LIBF)/filtrez/mkl_dfti.f90 ; \
	$(AR) r $(LIBO)/libfiltrez.a mkl_dfti.o ; $(RM) mkl_dfti.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(mkl_dft_type.o) : $(LIBF)/filtrez/mkl_dft_type.f90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libfiltrez.a mkl_dft_type.o ; \
	$(COMPILE90) $(LIBF)/filtrez/mkl_dft_type.f90 ; \
	$(AR) r $(LIBO)/libfiltrez.a mkl_dft_type.o ; $(RM) mkl_dft_type.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(mod_fft.o) : $(LIBF)/filtrez/mod_fft.F90 \
$(LIBO)/libfiltrez.a(mod_fft_fftw.o) \
$(LIBO)/libfiltrez.a(mod_fft_mathkeisan.o) \
$(LIBO)/libfiltrez.a(mod_fft_mkl.o) \
$(LIBO)/libfiltrez.a(mod_fft_wrapper.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libfiltrez.a mod_fft.o ; \
	$(COMPILE90) $(LIBF)/filtrez/mod_fft.F90 ; \
	$(AR) r $(LIBO)/libfiltrez.a mod_fft.o ; $(RM) mod_fft.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(mod_fft_fftw.o) : $(LIBF)/filtrez/mod_fft_fftw.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libfiltrez.a mod_fft_fftw.o ; \
	$(COMPILE90) $(LIBF)/filtrez/mod_fft_fftw.F90 ; \
	$(AR) r $(LIBO)/libfiltrez.a mod_fft_fftw.o ; $(RM) mod_fft_fftw.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(mod_fft_mathkeisan.o) : $(LIBF)/filtrez/mod_fft_mathkeisan.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libfiltrez.a mod_fft_mathkeisan.o ; \
	$(COMPILE90) $(LIBF)/filtrez/mod_fft_mathkeisan.F90 ; \
	$(AR) r $(LIBO)/libfiltrez.a mod_fft_mathkeisan.o ; $(RM) mod_fft_mathkeisan.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(mod_fft_mkl.o) : $(LIBF)/filtrez/mod_fft_mkl.F90 \
$(LIBO)/libfiltrez.a(mkl_dfti.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libfiltrez.a mod_fft_mkl.o ; \
	$(COMPILE90) $(LIBF)/filtrez/mod_fft_mkl.F90 ; \
	$(AR) r $(LIBO)/libfiltrez.a mod_fft_mkl.o ; $(RM) mod_fft_mkl.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(mod_fft_wrapper.o) : $(LIBF)/filtrez/mod_fft_wrapper.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libfiltrez.a mod_fft_wrapper.o ; \
	$(COMPILE90) $(LIBF)/filtrez/mod_fft_wrapper.F90 ; \
	$(AR) r $(LIBO)/libfiltrez.a mod_fft_wrapper.o ; $(RM) mod_fft_wrapper.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(mod_filtre_fft.o) : $(LIBF)/filtrez/mod_filtre_fft.F90 \
$(LIBO)/libfiltrez.a(mod_fft.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libfiltrez.a mod_filtre_fft.o ; \
	$(COMPILE90) $(LIBF)/filtrez/mod_filtre_fft.F90 ; \
	$(AR) r $(LIBO)/libfiltrez.a mod_filtre_fft.o ; $(RM) mod_filtre_fft.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(mod_filtre_fft_loc.o) : $(LIBF)/filtrez/mod_filtre_fft_loc.F90 \
$(LIBO)/libfiltrez.a(mod_fft.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libfiltrez.a mod_filtre_fft_loc.o ; \
	$(COMPILE90) $(LIBF)/filtrez/mod_filtre_fft_loc.F90 ; \
	$(AR) r $(LIBO)/libfiltrez.a mod_filtre_fft_loc.o ; $(RM) mod_filtre_fft_loc.o ; \
	cd $(GCM)

$(LIBO)/libfiltrez.a(timer_filtre.o) : $(LIBF)/filtrez/timer_filtre.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libfiltrez.a timer_filtre.o ; \
	$(COMPILE90) $(LIBF)/filtrez/timer_filtre.F90 ; \
	$(AR) r $(LIBO)/libfiltrez.a timer_filtre.o ; $(RM) timer_filtre.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory dyn3d_common
#=======================================================================

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(adaptdt.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(advn.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(advx.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(advxp.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(advy.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(advyp.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(advz.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(advzp.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(bernoui.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(convflu.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(coordij.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(covnat.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(defrun.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(diagedyn.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(disvert_noterre.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(diverg.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(divergf.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(diverg_gam.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(divergst.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(divgrad2.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(divgrad.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(dump2d.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(extrapol.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(fxy.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(fxyhyper.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(fxysinus.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(geopot.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(grad.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(gradiv2.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(gradiv.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(gr_ecrit_fi.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(gr_int_dyn.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(gr_u_scal.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(gr_v_scal.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(heavyside.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(inigeom.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(inigrads.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(initfluxsto.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(inithist.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(initial0.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(interpost.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(interpre.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(iso_verif_dyn.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(laplacien.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(laplacien_gam.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(laplacien_rot.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(laplacien_rotgam.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(limx.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(limy.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(limz.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(massdair.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(nxgrad.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(nxgrad_gam.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(nxgradst.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(nxgraro2.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(nxgrarot.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(pbar.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(pentes_ini.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(ppm3d.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(prather.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(psextbar.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(relaxetat0.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(rotat.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(rotatf.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(rotat_nfil.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(rotatst.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(sortvarc0.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(sortvarc.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(test_period.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(tidal_forces.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(traceurpole.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(writehist.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(comconst_mod.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(com_io_dyn_mod.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(comuforc_h.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(comvert_mod.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(conf_dat_m.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(conf_planete.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(control_mod.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(convmas.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(covcont.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(cpdet_mod.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(disvert.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(dynetat0.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(dynredem.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(dynredem_mod.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(enercin.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(ener_mod.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(exner_hyb_m.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(exner_milieu_m.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(flumass.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(fxhyp_m.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(fyhyp_m.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(grilles_gcm_netcdf_sub.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(infotrac.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(iniacademic.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(iniconst.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(inidissip.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(initdynav.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(inter_barxy_m.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(invert_lat.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(invert_zoom_x_m.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(massbar.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(massbarxy.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(misc_mod.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(pression.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(principal_cshift_m.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(serre_mod.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(tourpot.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(ugeostr.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(vitvert.o)

$(LIBO)/libdyn3d_common.a :  $(LIBO)/libdyn3d_common.a(writedynav.o)

.PRECIOUS	: $(LIBO)/libdyn3d_common.a


# Compilation of elements in dyn3d_common of library libdyn3d_common.a

$(LIBO)/libdyn3d_common.a(adaptdt.o) : $(LIBF)/dyn3d_common/adaptdt.F \
$(LIBF)/dyn3d_common/comdissip.h \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a adaptdt.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/adaptdt.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a adaptdt.o ; $(RM) adaptdt.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(advn.o) : $(LIBF)/dyn3d_common/advn.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a advn.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/advn.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a advn.o ; $(RM) advn.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(advx.o) : $(LIBF)/dyn3d_common/advx.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a advx.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/advx.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a advx.o ; $(RM) advx.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(advxp.o) : $(LIBF)/dyn3d_common/advxp.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a advxp.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/advxp.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a advxp.o ; $(RM) advxp.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(advy.o) : $(LIBF)/dyn3d_common/advy.F \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a advy.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/advy.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a advy.o ; $(RM) advy.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(advyp.o) : $(LIBF)/dyn3d_common/advyp.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a advyp.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/advyp.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a advyp.o ; $(RM) advyp.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(advz.o) : $(LIBF)/dyn3d_common/advz.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a advz.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/advz.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a advz.o ; $(RM) advz.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(advzp.o) : $(LIBF)/dyn3d_common/advzp.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a advzp.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/advzp.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a advzp.o ; $(RM) advzp.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(bernoui.o) : $(LIBF)/dyn3d_common/bernoui.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a bernoui.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/bernoui.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a bernoui.o ; $(RM) bernoui.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(convflu.o) : $(LIBF)/dyn3d_common/convflu.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a convflu.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/convflu.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a convflu.o ; $(RM) convflu.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(coordij.o) : $(LIBF)/dyn3d_common/coordij.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a coordij.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/coordij.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a coordij.o ; $(RM) coordij.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(covnat.o) : $(LIBF)/dyn3d_common/covnat.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a covnat.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/covnat.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a covnat.o ; $(RM) covnat.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(defrun.o) : $(LIBF)/dyn3d_common/defrun.F \
$(LIBF)/dyn3d_common/clesph0.h \
$(LIBF)/dyn3d_common/comdissnew.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libdyn3dpar.a(logic_mod.o) \
$(LIBO)/libdyn3d_common.a(serre_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a defrun.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/defrun.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a defrun.o ; $(RM) defrun.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(diagedyn.o) : $(LIBF)/dyn3d_common/diagedyn.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libdyn3d_common.a(cpdet_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a diagedyn.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/diagedyn.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a diagedyn.o ; $(RM) diagedyn.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(disvert_noterre.o) : $(LIBF)/dyn3d_common/disvert_noterre.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libmisc.a(ioipsl_getincom.o) \
$(LIBO)/libdyn3dpar.a(logic_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a disvert_noterre.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/disvert_noterre.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a disvert_noterre.o ; $(RM) disvert_noterre.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(diverg.o) : $(LIBF)/dyn3d_common/diverg.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a diverg.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/diverg.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a diverg.o ; $(RM) diverg.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(divergf.o) : $(LIBF)/dyn3d_common/divergf.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a divergf.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/divergf.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a divergf.o ; $(RM) divergf.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(diverg_gam.o) : $(LIBF)/dyn3d_common/diverg_gam.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a diverg_gam.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/diverg_gam.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a diverg_gam.o ; $(RM) diverg_gam.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(divergst.o) : $(LIBF)/dyn3d_common/divergst.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a divergst.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/divergst.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a divergst.o ; $(RM) divergst.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(divgrad2.o) : $(LIBF)/dyn3d_common/divgrad2.F \
$(LIBF)/dyn3d_common/comdissipn.h \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a divgrad2.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/divgrad2.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a divgrad2.o ; $(RM) divgrad2.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(divgrad.o) : $(LIBF)/dyn3d_common/divgrad.F \
$(LIBF)/dyn3d_common/comdissipn.h \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a divgrad.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/divgrad.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a divgrad.o ; $(RM) divgrad.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(dump2d.o) : $(LIBF)/dyn3d_common/dump2d.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a dump2d.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/dump2d.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a dump2d.o ; $(RM) dump2d.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(extrapol.o) : $(LIBF)/dyn3d_common/extrapol.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a extrapol.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/extrapol.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a extrapol.o ; $(RM) extrapol.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(fxy.o) : $(LIBF)/dyn3d_common/fxy.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/grid/fxy_new.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(serre_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a fxy.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/fxy.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a fxy.o ; $(RM) fxy.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(fxyhyper.o) : $(LIBF)/dyn3d_common/fxyhyper.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a fxyhyper.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/fxyhyper.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a fxyhyper.o ; $(RM) fxyhyper.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(fxysinus.o) : $(LIBF)/dyn3d_common/fxysinus.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/grid/fxy_sin.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a fxysinus.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/fxysinus.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a fxysinus.o ; $(RM) fxysinus.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(geopot.o) : $(LIBF)/dyn3d_common/geopot.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a geopot.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/geopot.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a geopot.o ; $(RM) geopot.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(grad.o) : $(LIBF)/dyn3d_common/grad.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a grad.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/grad.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a grad.o ; $(RM) grad.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(gradiv2.o) : $(LIBF)/dyn3d_common/gradiv2.F \
$(LIBF)/dyn3d_common/comdissipn.h \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a gradiv2.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/gradiv2.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a gradiv2.o ; $(RM) gradiv2.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(gradiv.o) : $(LIBF)/dyn3d_common/gradiv.F \
$(LIBF)/dyn3d_common/comdissipn.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a gradiv.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/gradiv.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a gradiv.o ; $(RM) gradiv.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(gr_ecrit_fi.o) : $(LIBF)/dyn3d_common/gr_ecrit_fi.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a gr_ecrit_fi.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/gr_ecrit_fi.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a gr_ecrit_fi.o ; $(RM) gr_ecrit_fi.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(gr_int_dyn.o) : $(LIBF)/dyn3d_common/gr_int_dyn.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a gr_int_dyn.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/gr_int_dyn.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a gr_int_dyn.o ; $(RM) gr_int_dyn.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(gr_u_scal.o) : $(LIBF)/dyn3d_common/gr_u_scal.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a gr_u_scal.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/gr_u_scal.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a gr_u_scal.o ; $(RM) gr_u_scal.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(gr_v_scal.o) : $(LIBF)/dyn3d_common/gr_v_scal.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a gr_v_scal.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/gr_v_scal.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a gr_v_scal.o ; $(RM) gr_v_scal.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(heavyside.o) : $(LIBF)/dyn3d_common/heavyside.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a heavyside.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/heavyside.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a heavyside.o ; $(RM) heavyside.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(inigeom.o) : $(LIBF)/dyn3d_common/inigeom.F \
$(LIBF)/dyn3d_common/comdissnew.h \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(fxhyp_m.o) \
$(LIBO)/libdyn3d_common.a(fyhyp_m.o) \
$(LIBO)/libdyn3dpar.a(logic_mod.o) \
$(LIBO)/libdyn3d_common.a(serre_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a inigeom.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/inigeom.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a inigeom.o ; $(RM) inigeom.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(inigrads.o) : $(LIBF)/dyn3d_common/inigrads.F \
$(LIBF)/dyn3d_common/gradsdef.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a inigrads.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/inigrads.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a inigrads.o ; $(RM) inigrads.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(initfluxsto.o) : $(LIBF)/dyn3d_common/initfluxsto.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3dpar.a(temps_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a initfluxsto.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/initfluxsto.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a initfluxsto.o ; $(RM) initfluxsto.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(inithist.o) : $(LIBF)/dyn3d_common/inithist.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(com_io_dyn_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libdyn3dpar.a(temps_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a inithist.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/inithist.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a inithist.o ; $(RM) inithist.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(initial0.o) : $(LIBF)/dyn3d_common/initial0.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a initial0.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/initial0.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a initial0.o ; $(RM) initial0.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(interpost.o) : $(LIBF)/dyn3d_common/interpost.F \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a interpost.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/interpost.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a interpost.o ; $(RM) interpost.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(interpre.o) : $(LIBF)/dyn3d_common/interpre.F \
$(LIBF)/dyn3d_common/comdissip.h \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a interpre.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/interpre.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a interpre.o ; $(RM) interpre.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(iso_verif_dyn.o) : $(LIBF)/dyn3d_common/iso_verif_dyn.F \
$(LIBO)/libdyn3d_common.a(infotrac.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a iso_verif_dyn.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/iso_verif_dyn.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a iso_verif_dyn.o ; $(RM) iso_verif_dyn.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(laplacien.o) : $(LIBF)/dyn3d_common/laplacien.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a laplacien.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/laplacien.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a laplacien.o ; $(RM) laplacien.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(laplacien_gam.o) : $(LIBF)/dyn3d_common/laplacien_gam.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a laplacien_gam.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/laplacien_gam.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a laplacien_gam.o ; $(RM) laplacien_gam.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(laplacien_rot.o) : $(LIBF)/dyn3d_common/laplacien_rot.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a laplacien_rot.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/laplacien_rot.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a laplacien_rot.o ; $(RM) laplacien_rot.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(laplacien_rotgam.o) : $(LIBF)/dyn3d_common/laplacien_rotgam.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a laplacien_rotgam.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/laplacien_rotgam.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a laplacien_rotgam.o ; $(RM) laplacien_rotgam.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(limx.o) : $(LIBF)/dyn3d_common/limx.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a limx.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/limx.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a limx.o ; $(RM) limx.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(limy.o) : $(LIBF)/dyn3d_common/limy.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a limy.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/limy.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a limy.o ; $(RM) limy.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(limz.o) : $(LIBF)/dyn3d_common/limz.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a limz.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/limz.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a limz.o ; $(RM) limz.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(massdair.o) : $(LIBF)/dyn3d_common/massdair.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a massdair.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/massdair.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a massdair.o ; $(RM) massdair.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(nxgrad.o) : $(LIBF)/dyn3d_common/nxgrad.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a nxgrad.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/nxgrad.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a nxgrad.o ; $(RM) nxgrad.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(nxgrad_gam.o) : $(LIBF)/dyn3d_common/nxgrad_gam.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a nxgrad_gam.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/nxgrad_gam.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a nxgrad_gam.o ; $(RM) nxgrad_gam.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(nxgradst.o) : $(LIBF)/dyn3d_common/nxgradst.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a nxgradst.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/nxgradst.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a nxgradst.o ; $(RM) nxgradst.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(nxgraro2.o) : $(LIBF)/dyn3d_common/nxgraro2.F \
$(LIBF)/dyn3d_common/comdissipn.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a nxgraro2.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/nxgraro2.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a nxgraro2.o ; $(RM) nxgraro2.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(nxgrarot.o) : $(LIBF)/dyn3d_common/nxgrarot.F \
$(LIBF)/dyn3d_common/comdissipn.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a nxgrarot.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/nxgrarot.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a nxgrarot.o ; $(RM) nxgrarot.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(pbar.o) : $(LIBF)/dyn3d_common/pbar.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a pbar.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/pbar.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a pbar.o ; $(RM) pbar.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(pentes_ini.o) : $(LIBF)/dyn3d_common/pentes_ini.F \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a pentes_ini.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/pentes_ini.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a pentes_ini.o ; $(RM) pentes_ini.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(ppm3d.o) : $(LIBF)/dyn3d_common/ppm3d.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a ppm3d.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/ppm3d.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a ppm3d.o ; $(RM) ppm3d.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(prather.o) : $(LIBF)/dyn3d_common/prather.F \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a prather.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/prather.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a prather.o ; $(RM) prather.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(psextbar.o) : $(LIBF)/dyn3d_common/psextbar.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a psextbar.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/psextbar.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a psextbar.o ; $(RM) psextbar.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(relaxetat0.o) : $(LIBF)/dyn3d_common/relaxetat0.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comuforc_h.o) \
$(LIBO)/libmisc.a(ioipsl_getincom.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a relaxetat0.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/relaxetat0.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a relaxetat0.o ; $(RM) relaxetat0.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(rotat.o) : $(LIBF)/dyn3d_common/rotat.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a rotat.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/rotat.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a rotat.o ; $(RM) rotat.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(rotatf.o) : $(LIBF)/dyn3d_common/rotatf.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a rotatf.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/rotatf.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a rotatf.o ; $(RM) rotatf.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(rotat_nfil.o) : $(LIBF)/dyn3d_common/rotat_nfil.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a rotat_nfil.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/rotat_nfil.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a rotat_nfil.o ; $(RM) rotat_nfil.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(rotatst.o) : $(LIBF)/dyn3d_common/rotatst.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a rotatst.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/rotatst.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a rotatst.o ; $(RM) rotatst.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(sortvarc0.o) : $(LIBF)/dyn3d_common/sortvarc0.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(ener_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a sortvarc0.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/sortvarc0.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a sortvarc0.o ; $(RM) sortvarc0.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(sortvarc.o) : $(LIBF)/dyn3d_common/sortvarc.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libdyn3d_common.a(ener_mod.o) \
$(LIBO)/libdyn3dpar.a(logic_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a sortvarc.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/sortvarc.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a sortvarc.o ; $(RM) sortvarc.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(test_period.o) : $(LIBF)/dyn3d_common/test_period.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(infotrac.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a test_period.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/test_period.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a test_period.o ; $(RM) test_period.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(tidal_forces.o) : $(LIBF)/dyn3d_common/tidal_forces.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a tidal_forces.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/tidal_forces.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a tidal_forces.o ; $(RM) tidal_forces.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(traceurpole.o) : $(LIBF)/dyn3d_common/traceurpole.F \
$(LIBF)/dyn3d_common/comdissip.h \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(control_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a traceurpole.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/traceurpole.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a traceurpole.o ; $(RM) traceurpole.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(writehist.o) : $(LIBF)/dyn3d_common/writehist.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(com_io_dyn_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libdyn3dpar.a(temps_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a writehist.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/writehist.F ; \
	$(AR) r $(LIBO)/libdyn3d_common.a writehist.o ; $(RM) writehist.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(comconst_mod.o) : $(LIBF)/dyn3d_common/comconst_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a comconst_mod.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/comconst_mod.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a comconst_mod.o ; $(RM) comconst_mod.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(com_io_dyn_mod.o) : $(LIBF)/dyn3d_common/com_io_dyn_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a com_io_dyn_mod.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/com_io_dyn_mod.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a com_io_dyn_mod.o ; $(RM) com_io_dyn_mod.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(comuforc_h.o) : $(LIBF)/dyn3d_common/comuforc_h.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a comuforc_h.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/comuforc_h.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a comuforc_h.o ; $(RM) comuforc_h.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(comvert_mod.o) : $(LIBF)/dyn3d_common/comvert_mod.F90 \
$(LIBF)/grid/dimensions.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a comvert_mod.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/comvert_mod.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a comvert_mod.o ; $(RM) comvert_mod.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(conf_dat_m.o) : $(LIBF)/dyn3d_common/conf_dat_m.F90 \
$(LIBO)/libmisc.a(assert_eq_m.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a conf_dat_m.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/conf_dat_m.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a conf_dat_m.o ; $(RM) conf_dat_m.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(conf_planete.o) : $(LIBF)/dyn3d_common/conf_planete.F90 \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libmisc.a(ioipsl_getincom.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a conf_planete.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/conf_planete.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a conf_planete.o ; $(RM) conf_planete.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(control_mod.o) : $(LIBF)/dyn3d_common/control_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a control_mod.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/control_mod.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a control_mod.o ; $(RM) control_mod.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(convmas.o) : $(LIBF)/dyn3d_common/convmas.F90 \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a convmas.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/convmas.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a convmas.o ; $(RM) convmas.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(covcont.o) : $(LIBF)/dyn3d_common/covcont.F90 \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a covcont.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/covcont.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a covcont.o ; $(RM) covcont.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(cpdet_mod.o) : $(LIBF)/dyn3d_common/cpdet_mod.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a cpdet_mod.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/cpdet_mod.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a cpdet_mod.o ; $(RM) cpdet_mod.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(disvert.o) : $(LIBF)/dyn3d_common/disvert.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libmisc.a(assert_m.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libmisc.a(ioipsl_getincom.o) \
$(LIBO)/libdyn3dpar.a(logic_mod.o) \
$(LIBO)/libmisc.a(new_unit_m.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a disvert.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/disvert.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a disvert.o ; $(RM) disvert.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(dynetat0.o) : $(LIBF)/dyn3d_common/dynetat0.F90 \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libdyn3d_common.a(ener_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libdyn3dpar.a(logic_mod.o) \
$(LIBO)/libdyn3d_common.a(serre_mod.o) \
$(LIBO)/libdyn3dpar.a(temps_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a dynetat0.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/dynetat0.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a dynetat0.o ; $(RM) dynetat0.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(dynredem.o) : $(LIBF)/dyn3d_common/dynredem.F90 \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libdyn3d_common.a(dynredem_mod.o) \
$(LIBO)/libdyn3d_common.a(ener_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libdyn3dpar.a(logic_mod.o) \
$(LIBO)/libmisc.a(netcdf95.o) \
$(LIBO)/libdyn3d_common.a(serre_mod.o) \
$(LIBO)/libdyn3dpar.a(temps_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a dynredem.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/dynredem.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a dynredem.o ; $(RM) dynredem.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(dynredem_mod.o) : $(LIBF)/dyn3d_common/dynredem_mod.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a dynredem_mod.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/dynredem_mod.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a dynredem_mod.o ; $(RM) dynredem_mod.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(enercin.o) : $(LIBF)/dyn3d_common/enercin.F90 \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a enercin.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/enercin.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a enercin.o ; $(RM) enercin.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(ener_mod.o) : $(LIBF)/dyn3d_common/ener_mod.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a ener_mod.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/ener_mod.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a ener_mod.o ; $(RM) ener_mod.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(exner_hyb_m.o) : $(LIBF)/dyn3d_common/exner_hyb_m.F90 \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a exner_hyb_m.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/exner_hyb_m.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a exner_hyb_m.o ; $(RM) exner_hyb_m.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(exner_milieu_m.o) : $(LIBF)/dyn3d_common/exner_milieu_m.F90 \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a exner_milieu_m.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/exner_milieu_m.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a exner_milieu_m.o ; $(RM) exner_milieu_m.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(flumass.o) : $(LIBF)/dyn3d_common/flumass.F90 \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a flumass.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/flumass.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a flumass.o ; $(RM) flumass.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(fxhyp_m.o) : $(LIBF)/dyn3d_common/fxhyp_m.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBO)/libmisc.a(arth_m.o) \
$(LIBO)/libdyn3d_common.a(invert_zoom_x_m.o) \
$(LIBO)/libmisc.a(nrtype.o) \
$(LIBO)/libdyn3d_common.a(principal_cshift_m.o) \
$(LIBO)/libdyn3d_common.a(serre_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a fxhyp_m.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/fxhyp_m.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a fxhyp_m.o ; $(RM) fxhyp_m.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(fyhyp_m.o) : $(LIBF)/dyn3d_common/fyhyp_m.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBO)/libmisc.a(coefpoly_m.o) \
$(LIBO)/libmisc.a(nrtype.o) \
$(LIBO)/libdyn3d_common.a(serre_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a fyhyp_m.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/fyhyp_m.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a fyhyp_m.o ; $(RM) fyhyp_m.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(grilles_gcm_netcdf_sub.o) : $(LIBF)/dyn3d_common/grilles_gcm_netcdf_sub.F90 \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a grilles_gcm_netcdf_sub.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/grilles_gcm_netcdf_sub.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a grilles_gcm_netcdf_sub.o ; $(RM) grilles_gcm_netcdf_sub.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(infotrac.o) : $(LIBF)/dyn3d_common/infotrac.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libmisc.a(ioipsl_getincom.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a infotrac.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/infotrac.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a infotrac.o ; $(RM) infotrac.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(iniacademic.o) : $(LIBF)/dyn3d_common/iniacademic.F90 \
$(LIBF)/dyn3d_common/academic.h \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libdyn3d_common.a(ener_mod.o) \
$(LIBO)/libdyn3d_common.a(exner_hyb_m.o) \
$(LIBO)/libdyn3d_common.a(exner_milieu_m.o) \
$(LIBO)/libfiltrez.a(filtreg_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libmisc.a(ioipsl_getincom.o) \
$(LIBO)/libdyn3dpar.a(logic_mod.o) \
$(LIBO)/libdyn3dpar.a(temps_mod.o) \
$(LIBO)/libmisc.a(write_field.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a iniacademic.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/iniacademic.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a iniacademic.o ; $(RM) iniacademic.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(iniconst.o) : $(LIBF)/dyn3d_common/iniconst.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libmisc.a(ioipsl_getincom.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a iniconst.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/iniconst.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a iniconst.o ; $(RM) iniconst.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(inidissip.o) : $(LIBF)/dyn3d_common/inidissip.F90 \
$(LIBF)/dyn3d_common/comdissipn.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libdyn3dpar.a(logic_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a inidissip.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/inidissip.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a inidissip.o ; $(RM) inidissip.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(initdynav.o) : $(LIBF)/dyn3d_common/initdynav.F90 \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(com_io_dyn_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libdyn3dpar.a(temps_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a initdynav.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/initdynav.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a initdynav.o ; $(RM) initdynav.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(inter_barxy_m.o) : $(LIBF)/dyn3d_common/inter_barxy_m.F90 \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libmisc.a(assert_eq_m.o) \
$(LIBO)/libmisc.a(assert_m.o) \
$(LIBO)/libdyn3d_common.a(comconst_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a inter_barxy_m.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/inter_barxy_m.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a inter_barxy_m.o ; $(RM) inter_barxy_m.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(invert_lat.o) : $(LIBF)/dyn3d_common/invert_lat.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a invert_lat.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/invert_lat.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a invert_lat.o ; $(RM) invert_lat.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(invert_zoom_x_m.o) : $(LIBF)/dyn3d_common/invert_zoom_x_m.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBO)/libmisc.a(coefpoly_m.o) \
$(LIBO)/libmisc.a(nrtype.o) \
$(LIBO)/libdyn3d_common.a(serre_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a invert_zoom_x_m.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/invert_zoom_x_m.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a invert_zoom_x_m.o ; $(RM) invert_zoom_x_m.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(massbar.o) : $(LIBF)/dyn3d_common/massbar.F90 \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a massbar.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/massbar.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a massbar.o ; $(RM) massbar.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(massbarxy.o) : $(LIBF)/dyn3d_common/massbarxy.F90 \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a massbarxy.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/massbarxy.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a massbarxy.o ; $(RM) massbarxy.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(misc_mod.o) : $(LIBF)/dyn3d_common/misc_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a misc_mod.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/misc_mod.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a misc_mod.o ; $(RM) misc_mod.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(pression.o) : $(LIBF)/dyn3d_common/pression.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a pression.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/pression.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a pression.o ; $(RM) pression.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(principal_cshift_m.o) : $(LIBF)/dyn3d_common/principal_cshift_m.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBO)/libmisc.a(nrtype.o) \
$(LIBO)/libdyn3d_common.a(serre_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a principal_cshift_m.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/principal_cshift_m.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a principal_cshift_m.o ; $(RM) principal_cshift_m.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(serre_mod.o) : $(LIBF)/dyn3d_common/serre_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a serre_mod.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/serre_mod.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a serre_mod.o ; $(RM) serre_mod.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(tourpot.o) : $(LIBF)/dyn3d_common/tourpot.F90 \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a tourpot.o ; \
	$(COMPILE) $(LIBF)/dyn3d_common/tourpot.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a tourpot.o ; $(RM) tourpot.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(ugeostr.o) : $(LIBF)/dyn3d_common/ugeostr.F90 \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a ugeostr.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/ugeostr.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a ugeostr.o ; $(RM) ugeostr.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(vitvert.o) : $(LIBF)/dyn3d_common/vitvert.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comvert_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a vitvert.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/vitvert.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a vitvert.o ; $(RM) vitvert.o ; \
	cd $(GCM)

$(LIBO)/libdyn3d_common.a(writedynav.o) : $(LIBF)/dyn3d_common/writedynav.F90 \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(com_io_dyn_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libdyn3dpar.a(temps_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3d_common.a writedynav.o ; \
	$(COMPILE90) $(LIBF)/dyn3d_common/writedynav.F90 ; \
	$(AR) r $(LIBO)/libdyn3d_common.a writedynav.o ; $(RM) writedynav.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory dyn3dpar
#=======================================================================

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(abort_gcm.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(addfi_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(advect_new_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(advect_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(bernoui_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(bilan_dyn_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(caldyn_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(check_isotopes_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(convflu_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(convmas1_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(convmas2_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(convmas_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(covcont_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(covnat_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(dissip_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(divergf_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(diverg_gam_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(diverg_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(divgrad2_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(divgrad_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(dteta1_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(dudv1_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(dudv2_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(enercin_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(filtreg_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(flumass_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(fluxstokenc_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(friction_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(geopot_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(gradiv2_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(gradiv_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(grad_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(groupe_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(groupeun_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(gr_u_scal_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(gr_v_scal_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(initdynav_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(initfluxsto_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(inithist_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(integrd_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(laplacien_gam_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(laplacien_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(laplacien_rotgam_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(laplacien_rot_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(leapfrog_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(massbar_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(massbarxy_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(massdair_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(nxgrad_gam_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(nxgrad_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(nxgraro2_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(nxgrarot_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(pression_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(qminimum_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(rotatf_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(rotat_nfil_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(rotat_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(sw_case_williamson91_6.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(top_bound_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(tourpot_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(vitvert_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(vlspltgen_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(vlsplt_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(vlspltqs_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(wrgrads.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(writedynav_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(writehist_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(advtrac_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(bands.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(conf_gcm.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(dynredem_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(exner_hyb_p_m.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(exner_milieu_p_m.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(getparam.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(guide_p_mod.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(logic_mod.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(mod_const_mpi.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(mod_hallo.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(moyzon_mod.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(parallel_lmdz.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(planetary_operations_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(sponge_mod_p.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(temps_mod.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(times.o)

$(LIBO)/libdyn3dpar.a :  $(LIBO)/libdyn3dpar.a(write_field_p.o)

.PRECIOUS	: $(LIBO)/libdyn3dpar.a


# Compilation of elements in dyn3dpar of library libdyn3dpar.a

$(LIBO)/libdyn3dpar.a(abort_gcm.o) : $(LIBF)/dyn3dpar/abort_gcm.F \
$(LIBF)/misc/iniprint.h \
$(LIBO)/libmisc.a(ioipsl_getincom.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a abort_gcm.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/abort_gcm.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a abort_gcm.o ; $(RM) abort_gcm.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(addfi_p.o) : $(LIBF)/dyn3dpar/addfi_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a addfi_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/addfi_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a addfi_p.o ; $(RM) addfi_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(advect_new_p.o) : $(LIBF)/dyn3dpar/advect_new_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3dpar.a(logic_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(write_field_p.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a advect_new_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/advect_new_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a advect_new_p.o ; $(RM) advect_new_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(advect_p.o) : $(LIBF)/dyn3dpar/advect_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3dpar.a(logic_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(write_field_p.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a advect_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/advect_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a advect_p.o ; $(RM) advect_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(bernoui_p.o) : $(LIBF)/dyn3dpar/bernoui_p.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a bernoui_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/bernoui_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a bernoui_p.o ; $(RM) bernoui_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(bilan_dyn_p.o) : $(LIBF)/dyn3dpar/bilan_dyn_p.F \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(misc_mod.o) \
$(LIBO)/libdyn3dpar.a(mod_hallo.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(temps_mod.o) \
$(LIBO)/libdyn3dpar.a(write_field_p.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a bilan_dyn_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/bilan_dyn_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a bilan_dyn_p.o ; $(RM) bilan_dyn_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(caldyn_p.o) : $(LIBF)/dyn3dpar/caldyn_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(write_field_p.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a caldyn_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/caldyn_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a caldyn_p.o ; $(RM) caldyn_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(check_isotopes_p.o) : $(LIBF)/dyn3dpar/check_isotopes_p.F \
$(LIBF)/grid/dimensions.h \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a check_isotopes_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/check_isotopes_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a check_isotopes_p.o ; $(RM) check_isotopes_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(convflu_p.o) : $(LIBF)/dyn3dpar/convflu_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a convflu_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/convflu_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a convflu_p.o ; $(RM) convflu_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(convmas1_p.o) : $(LIBF)/dyn3dpar/convmas1_p.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a convmas1_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/convmas1_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a convmas1_p.o ; $(RM) convmas1_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(convmas2_p.o) : $(LIBF)/dyn3dpar/convmas2_p.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a convmas2_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/convmas2_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a convmas2_p.o ; $(RM) convmas2_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(convmas_p.o) : $(LIBF)/dyn3dpar/convmas_p.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a convmas_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/convmas_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a convmas_p.o ; $(RM) convmas_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(covcont_p.o) : $(LIBF)/dyn3dpar/covcont_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a covcont_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/covcont_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a covcont_p.o ; $(RM) covcont_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(covnat_p.o) : $(LIBF)/dyn3dpar/covnat_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a covnat_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/covnat_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a covnat_p.o ; $(RM) covnat_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(dissip_p.o) : $(LIBF)/dyn3dpar/dissip_p.F \
$(LIBF)/dyn3d_common/comdissipn.h \
$(LIBF)/dyn3d_common/comdissnew.h \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(write_field_p.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a dissip_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/dissip_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a dissip_p.o ; $(RM) dissip_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(divergf_p.o) : $(LIBF)/dyn3dpar/divergf_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a divergf_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/divergf_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a divergf_p.o ; $(RM) divergf_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(diverg_gam_p.o) : $(LIBF)/dyn3dpar/diverg_gam_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a diverg_gam_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/diverg_gam_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a diverg_gam_p.o ; $(RM) diverg_gam_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(diverg_p.o) : $(LIBF)/dyn3dpar/diverg_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a diverg_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/diverg_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a diverg_p.o ; $(RM) diverg_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(divgrad2_p.o) : $(LIBF)/dyn3dpar/divgrad2_p.F \
$(LIBF)/dyn3d_common/comdissipn.h \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(mod_hallo.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(times.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a divgrad2_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/divgrad2_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a divgrad2_p.o ; $(RM) divgrad2_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(divgrad_p.o) : $(LIBF)/dyn3dpar/divgrad_p.F \
$(LIBF)/dyn3d_common/comdissipn.h \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(times.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a divgrad_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/divgrad_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a divgrad_p.o ; $(RM) divgrad_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(dteta1_p.o) : $(LIBF)/dyn3dpar/dteta1_p.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(write_field_p.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a dteta1_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/dteta1_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a dteta1_p.o ; $(RM) dteta1_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(dudv1_p.o) : $(LIBF)/dyn3dpar/dudv1_p.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a dudv1_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/dudv1_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a dudv1_p.o ; $(RM) dudv1_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(dudv2_p.o) : $(LIBF)/dyn3dpar/dudv2_p.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a dudv2_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/dudv2_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a dudv2_p.o ; $(RM) dudv2_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(enercin_p.o) : $(LIBF)/dyn3dpar/enercin_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a enercin_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/enercin_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a enercin_p.o ; $(RM) enercin_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(filtreg_p.o) : $(LIBF)/dyn3dpar/filtreg_p.F \
$(LIBF)/filtrez/coefils.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libfiltrez.a(filtreg_mod.o) \
$(LIBO)/libfiltrez.a(mod_filtre_fft.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libfiltrez.a(timer_filtre.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a filtreg_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/filtreg_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a filtreg_p.o ; $(RM) filtreg_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(flumass_p.o) : $(LIBF)/dyn3dpar/flumass_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a flumass_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/flumass_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a flumass_p.o ; $(RM) flumass_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(fluxstokenc_p.o) : $(LIBF)/dyn3dpar/fluxstokenc_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBF)/dyn3d_common/tracstoke.h \
$(LIBO)/libdyn3d_common.a(misc_mod.o) \
$(LIBO)/libdyn3dpar.a(mod_hallo.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a fluxstokenc_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/fluxstokenc_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a fluxstokenc_p.o ; $(RM) fluxstokenc_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(friction_p.o) : $(LIBF)/dyn3dpar/friction_p.F \
$(LIBF)/dyn3d_common/academic.h \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libmisc.a(ioipsl_getincom.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a friction_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/friction_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a friction_p.o ; $(RM) friction_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(geopot_p.o) : $(LIBF)/dyn3dpar/geopot_p.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a geopot_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/geopot_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a geopot_p.o ; $(RM) geopot_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(gradiv2_p.o) : $(LIBF)/dyn3dpar/gradiv2_p.F \
$(LIBF)/dyn3d_common/comdissipn.h \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(mod_hallo.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(times.o) \
$(LIBO)/libdyn3dpar.a(write_field_p.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a gradiv2_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/gradiv2_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a gradiv2_p.o ; $(RM) gradiv2_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(gradiv_p.o) : $(LIBF)/dyn3dpar/gradiv_p.F \
$(LIBF)/dyn3d_common/comdissipn.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(times.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a gradiv_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/gradiv_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a gradiv_p.o ; $(RM) gradiv_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(grad_p.o) : $(LIBF)/dyn3dpar/grad_p.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a grad_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/grad_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a grad_p.o ; $(RM) grad_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(groupe_p.o) : $(LIBF)/dyn3dpar/groupe_p.F \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a groupe_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/groupe_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a groupe_p.o ; $(RM) groupe_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(groupeun_p.o) : $(LIBF)/dyn3dpar/groupeun_p.F \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(write_field_p.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a groupeun_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/groupeun_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a groupeun_p.o ; $(RM) groupeun_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(gr_u_scal_p.o) : $(LIBF)/dyn3dpar/gr_u_scal_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a gr_u_scal_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/gr_u_scal_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a gr_u_scal_p.o ; $(RM) gr_u_scal_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(gr_v_scal_p.o) : $(LIBF)/dyn3dpar/gr_v_scal_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a gr_v_scal_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/gr_v_scal_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a gr_v_scal_p.o ; $(RM) gr_v_scal_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(initdynav_p.o) : $(LIBF)/dyn3dpar/initdynav_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libdyn3d_common.a(misc_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(temps_mod.o) \
$(LIBO)/libmisc.a(write_field.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a initdynav_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/initdynav_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a initdynav_p.o ; $(RM) initdynav_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(initfluxsto_p.o) : $(LIBF)/dyn3dpar/initfluxsto_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(misc_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(temps_mod.o) \
$(LIBO)/libmisc.a(write_field.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a initfluxsto_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/initfluxsto_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a initfluxsto_p.o ; $(RM) initfluxsto_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(inithist_p.o) : $(LIBF)/dyn3dpar/inithist_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libdyn3d_common.a(misc_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(temps_mod.o) \
$(LIBO)/libmisc.a(write_field.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a inithist_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/inithist_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a inithist_p.o ; $(RM) inithist_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(integrd_p.o) : $(LIBF)/dyn3dpar/integrd_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libdyn3dpar.a(logic_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(temps_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a integrd_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/integrd_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a integrd_p.o ; $(RM) integrd_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(laplacien_gam_p.o) : $(LIBF)/dyn3dpar/laplacien_gam_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a laplacien_gam_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/laplacien_gam_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a laplacien_gam_p.o ; $(RM) laplacien_gam_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(laplacien_p.o) : $(LIBF)/dyn3dpar/laplacien_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a laplacien_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/laplacien_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a laplacien_p.o ; $(RM) laplacien_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(laplacien_rotgam_p.o) : $(LIBF)/dyn3dpar/laplacien_rotgam_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a laplacien_rotgam_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/laplacien_rotgam_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a laplacien_rotgam_p.o ; $(RM) laplacien_rotgam_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(laplacien_rot_p.o) : $(LIBF)/dyn3dpar/laplacien_rot_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a laplacien_rot_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/laplacien_rot_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a laplacien_rot_p.o ; $(RM) laplacien_rot_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(leapfrog_p.o) : $(LIBF)/dyn3dpar/leapfrog_p.F \
$(LIBF)/dyn3d_common/academic.h \
$(LIBF)/dyn3d_common/comdissnew.h \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBF)/dyn3d_common/write_grads_dyn.h \
$(LIBO)/libdyn3dpar.a(bands.o) \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comuforc_h.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libdyn3d_common.a(cpdet_mod.o) \
$(LIBO)/libdyn3d_common.a(exner_hyb_m.o) \
$(LIBO)/libdyn3dpar.a(exner_hyb_p_m.o) \
$(LIBO)/libdyn3d_common.a(exner_milieu_m.o) \
$(LIBO)/libdyn3dpar.a(exner_milieu_p_m.o) \
$(LIBO)/libdyn3dpar.a(getparam.o) \
$(LIBO)/libdyn3dpar.a(guide_p_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libdyn3dpar.a(logic_mod.o) \
$(LIBO)/libdyn3d_common.a(misc_mod.o) \
$(LIBO)/libdyn3dpar.a(mod_hallo.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(sponge_mod_p.o) \
$(LIBO)/libdyn3dpar.a(temps_mod.o) \
$(LIBO)/libfiltrez.a(timer_filtre.o) \
$(LIBO)/libdyn3dpar.a(times.o) \
$(LIBO)/libmisc.a(vampir.o) \
$(LIBO)/libmisc.a(write_field.o) \
$(LIBO)/libdyn3dpar.a(write_field_p.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a leapfrog_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/leapfrog_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a leapfrog_p.o ; $(RM) leapfrog_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(massbar_p.o) : $(LIBF)/dyn3dpar/massbar_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a massbar_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/massbar_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a massbar_p.o ; $(RM) massbar_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(massbarxy_p.o) : $(LIBF)/dyn3dpar/massbarxy_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a massbarxy_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/massbarxy_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a massbarxy_p.o ; $(RM) massbarxy_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(massdair_p.o) : $(LIBF)/dyn3dpar/massdair_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a massdair_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/massdair_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a massdair_p.o ; $(RM) massdair_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(nxgrad_gam_p.o) : $(LIBF)/dyn3dpar/nxgrad_gam_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a nxgrad_gam_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/nxgrad_gam_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a nxgrad_gam_p.o ; $(RM) nxgrad_gam_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(nxgrad_p.o) : $(LIBF)/dyn3dpar/nxgrad_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a nxgrad_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/nxgrad_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a nxgrad_p.o ; $(RM) nxgrad_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(nxgraro2_p.o) : $(LIBF)/dyn3dpar/nxgraro2_p.F \
$(LIBF)/dyn3d_common/comdissipn.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(mod_hallo.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(times.o) \
$(LIBO)/libdyn3dpar.a(write_field_p.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a nxgraro2_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/nxgraro2_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a nxgraro2_p.o ; $(RM) nxgraro2_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(nxgrarot_p.o) : $(LIBF)/dyn3dpar/nxgrarot_p.F \
$(LIBF)/dyn3d_common/comdissipn.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(times.o) \
$(LIBO)/libdyn3dpar.a(write_field_p.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a nxgrarot_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/nxgrarot_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a nxgrarot_p.o ; $(RM) nxgrarot_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(pression_p.o) : $(LIBF)/dyn3dpar/pression_p.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a pression_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/pression_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a pression_p.o ; $(RM) pression_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(qminimum_p.o) : $(LIBF)/dyn3dpar/qminimum_p.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a qminimum_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/qminimum_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a qminimum_p.o ; $(RM) qminimum_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(rotatf_p.o) : $(LIBF)/dyn3dpar/rotatf_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a rotatf_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/rotatf_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a rotatf_p.o ; $(RM) rotatf_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(rotat_nfil_p.o) : $(LIBF)/dyn3dpar/rotat_nfil_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a rotat_nfil_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/rotat_nfil_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a rotat_nfil_p.o ; $(RM) rotat_nfil_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(rotat_p.o) : $(LIBF)/dyn3dpar/rotat_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a rotat_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/rotat_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a rotat_p.o ; $(RM) rotat_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(sw_case_williamson91_6.o) : $(LIBF)/dyn3dpar/sw_case_williamson91_6.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a sw_case_williamson91_6.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/sw_case_williamson91_6.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a sw_case_williamson91_6.o ; $(RM) sw_case_williamson91_6.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(top_bound_p.o) : $(LIBF)/dyn3dpar/top_bound_p.F \
$(LIBF)/dyn3d_common/comdissipn.h \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a top_bound_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/top_bound_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a top_bound_p.o ; $(RM) top_bound_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(tourpot_p.o) : $(LIBF)/dyn3dpar/tourpot_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a tourpot_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/tourpot_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a tourpot_p.o ; $(RM) tourpot_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(vitvert_p.o) : $(LIBF)/dyn3dpar/vitvert_p.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a vitvert_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/vitvert_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a vitvert_p.o ; $(RM) vitvert_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(vlspltgen_p.o) : $(LIBF)/dyn3dpar/vlspltgen_p.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libdyn3dpar.a(mod_hallo.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libmisc.a(vampir.o) \
$(LIBO)/libdyn3dpar.a(write_field_p.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a vlspltgen_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/vlspltgen_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a vlspltgen_p.o ; $(RM) vlspltgen_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(vlsplt_p.o) : $(LIBF)/dyn3dpar/vlsplt_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libdyn3dpar.a(mod_hallo.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libmisc.a(vampir.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a vlsplt_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/vlsplt_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a vlsplt_p.o ; $(RM) vlsplt_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(vlspltqs_p.o) : $(LIBF)/dyn3dpar/vlspltqs_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(cpdet_mod.o) \
$(LIBO)/libdyn3dpar.a(mod_hallo.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libmisc.a(vampir.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a vlspltqs_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/vlspltqs_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a vlspltqs_p.o ; $(RM) vlspltqs_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(wrgrads.o) : $(LIBF)/dyn3dpar/wrgrads.F \
$(LIBF)/dyn3d_common/gradsdef.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a wrgrads.o ; \
	$(COMPILE) $(LIBF)/dyn3dpar/wrgrads.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a wrgrads.o ; $(RM) wrgrads.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(writedynav_p.o) : $(LIBF)/dyn3dpar/writedynav_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libdyn3d_common.a(misc_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(temps_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a writedynav_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/writedynav_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a writedynav_p.o ; $(RM) writedynav_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(writehist_p.o) : $(LIBF)/dyn3dpar/writehist_p.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libdyn3d_common.a(misc_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(temps_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a writehist_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/writehist_p.F ; \
	$(AR) r $(LIBO)/libdyn3dpar.a writehist_p.o ; $(RM) writehist_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(advtrac_p.o) : $(LIBF)/dyn3dpar/advtrac_p.F90 \
$(LIBF)/dyn3d_common/comdissip.h \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(bands.o) \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libdyn3dpar.a(mod_hallo.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(planetary_operations_p.o) \
$(LIBO)/libdyn3dpar.a(times.o) \
$(LIBO)/libmisc.a(vampir.o) \
$(LIBO)/libdyn3dpar.a(write_field_p.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a advtrac_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/advtrac_p.F90 ; \
	$(AR) r $(LIBO)/libdyn3dpar.a advtrac_p.o ; $(RM) advtrac_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(bands.o) : $(LIBF)/dyn3dpar/bands.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(times.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a bands.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/bands.F90 ; \
	$(AR) r $(LIBO)/libdyn3dpar.a bands.o ; $(RM) bands.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(conf_gcm.o) : $(LIBF)/dyn3dpar/conf_gcm.F90 \
$(LIBF)/dyn3d_common/comdissnew.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libmisc.a(assert_m.o) \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libmisc.a(ioipsl_getincom.o) \
$(LIBO)/libdyn3dpar.a(logic_mod.o) \
$(LIBO)/libdyn3d_common.a(misc_mod.o) \
$(LIBO)/libfiltrez.a(mod_filtre_fft.o) \
$(LIBO)/libdyn3dpar.a(mod_hallo.o) \
$(LIBO)/libdyn3d_common.a(serre_mod.o) \
$(LIBO)/libdyn3dpar.a(sponge_mod_p.o) \
$(LIBO)/libdyn3dpar.a(temps_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a conf_gcm.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/conf_gcm.F90 ; \
	$(AR) r $(LIBO)/libdyn3dpar.a conf_gcm.o ; $(RM) conf_gcm.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(dynredem_p.o) : $(LIBF)/dyn3dpar/dynredem_p.F90 \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libdyn3d_common.a(dynredem_mod.o) \
$(LIBO)/libdyn3d_common.a(ener_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libdyn3dpar.a(logic_mod.o) \
$(LIBO)/libmisc.a(netcdf95.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3d_common.a(serre_mod.o) \
$(LIBO)/libdyn3dpar.a(temps_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a dynredem_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/dynredem_p.F90 ; \
	$(AR) r $(LIBO)/libdyn3dpar.a dynredem_p.o ; $(RM) dynredem_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(exner_hyb_p_m.o) : $(LIBF)/dyn3dpar/exner_hyb_p_m.F90 \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a exner_hyb_p_m.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/exner_hyb_p_m.F90 ; \
	$(AR) r $(LIBO)/libdyn3dpar.a exner_hyb_p_m.o ; $(RM) exner_hyb_p_m.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(exner_milieu_p_m.o) : $(LIBF)/dyn3dpar/exner_milieu_p_m.F90 \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a exner_milieu_p_m.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/exner_milieu_p_m.F90 ; \
	$(AR) r $(LIBO)/libdyn3dpar.a exner_milieu_p_m.o ; $(RM) exner_milieu_p_m.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(getparam.o) : $(LIBF)/dyn3dpar/getparam.F90 \
$(LIBO)/libmisc.a(ioipsl_getincom.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a getparam.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/getparam.F90 ; \
	$(AR) r $(LIBO)/libdyn3dpar.a getparam.o ; $(RM) getparam.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(guide_p_mod.o) : $(LIBF)/dyn3dpar/guide_p_mod.F90 \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(bands.o) \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libdyn3dpar.a(exner_hyb_p_m.o) \
$(LIBO)/libdyn3dpar.a(exner_milieu_p_m.o) \
$(LIBO)/libdyn3dpar.a(getparam.o) \
$(LIBO)/libdyn3dpar.a(mod_hallo.o) \
$(LIBO)/libmisc.a(netcdf95.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libmisc.a(pres2lev_mod.o) \
$(LIBO)/libdyn3d_common.a(serre_mod.o) \
$(LIBO)/libdyn3dpar.a(write_field_p.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a guide_p_mod.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/guide_p_mod.F90 ; \
	$(AR) r $(LIBO)/libdyn3dpar.a guide_p_mod.o ; $(RM) guide_p_mod.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(logic_mod.o) : $(LIBF)/dyn3dpar/logic_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a logic_mod.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/logic_mod.F90 ; \
	$(AR) r $(LIBO)/libdyn3dpar.a logic_mod.o ; $(RM) logic_mod.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(mod_const_mpi.o) : $(LIBF)/dyn3dpar/mod_const_mpi.F90 \
$(LIBO)/libmisc.a(ioipsl_getincom.o) \
$(LIBO)/libmisc.a(wxios.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a mod_const_mpi.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/mod_const_mpi.F90 ; \
	$(AR) r $(LIBO)/libdyn3dpar.a mod_const_mpi.o ; $(RM) mod_const_mpi.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(mod_hallo.o) : $(LIBF)/dyn3dpar/mod_hallo.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(mod_const_mpi.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a mod_hallo.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/mod_hallo.F90 ; \
	$(AR) r $(LIBO)/libdyn3dpar.a mod_hallo.o ; $(RM) mod_hallo.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(moyzon_mod.o) : $(LIBF)/dyn3dpar/moyzon_mod.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a moyzon_mod.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/moyzon_mod.F90 ; \
	$(AR) r $(LIBO)/libdyn3dpar.a moyzon_mod.o ; $(RM) moyzon_mod.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) : $(LIBF)/dyn3dpar/parallel_lmdz.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libmisc.a(ioipsl_getincom.o) \
$(LIBO)/libdyn3dpar.a(mod_const_mpi.o) \
$(LIBO)/libmisc.a(vampir.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a parallel_lmdz.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/parallel_lmdz.F90 ; \
	$(AR) r $(LIBO)/libdyn3dpar.a parallel_lmdz.o ; $(RM) parallel_lmdz.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(planetary_operations_p.o) : $(LIBF)/dyn3dpar/planetary_operations_p.F90 \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(mod_const_mpi.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a planetary_operations_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/planetary_operations_p.F90 ; \
	$(AR) r $(LIBO)/libdyn3dpar.a planetary_operations_p.o ; $(RM) planetary_operations_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(sponge_mod_p.o) : $(LIBF)/dyn3dpar/sponge_mod_p.F90 \
$(LIBF)/dyn3d_common/comdissip.h \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(write_field_p.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a sponge_mod_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/sponge_mod_p.F90 ; \
	$(AR) r $(LIBO)/libdyn3dpar.a sponge_mod_p.o ; $(RM) sponge_mod_p.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(temps_mod.o) : $(LIBF)/dyn3dpar/temps_mod.F90
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a temps_mod.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/temps_mod.F90 ; \
	$(AR) r $(LIBO)/libdyn3dpar.a temps_mod.o ; $(RM) temps_mod.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(times.o) : $(LIBF)/dyn3dpar/times.F90 \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a times.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/times.F90 ; \
	$(AR) r $(LIBO)/libdyn3dpar.a times.o ; $(RM) times.o ; \
	cd $(GCM)

$(LIBO)/libdyn3dpar.a(write_field_p.o) : $(LIBF)/dyn3dpar/write_field_p.F90 \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libmisc.a(write_field.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdyn3dpar.a write_field_p.o ; \
	$(COMPILE90) $(LIBF)/dyn3dpar/write_field_p.F90 ; \
	$(AR) r $(LIBO)/libdyn3dpar.a write_field_p.o ; $(RM) write_field_p.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory dynphy_lonlat
#=======================================================================

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(calfis.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(calfis_p.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(gr_dyn_fi.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(gr_dyn_fi_p.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(gr_fi_dyn.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(gr_fi_dyn_p.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(grid_atob_m.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(inigeomphy_mod.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(mod_interface_dyn_phys.o)

.PRECIOUS	: $(LIBO)/libdynphy_lonlat.a


# Compilation of elements in dynphy_lonlat of library libdynphy_lonlat.a

$(LIBO)/libdynphy_lonlat.a(calfis.o) : $(LIBF)/dynphy_lonlat/calfis.F \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdynphy_lonlat.a(callphysiq_mod.o) \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libdyn3d_common.a(cpdet_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libdyn3dpar.a(logic_mod.o) \
$(LIBO)/libdyn3dpar.a(moyzon_mod.o) \
$(LIBO)/libmisc.a(write_field.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a calfis.o ; \
	$(COMPILE90) $(LIBF)/dynphy_lonlat/calfis.F ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a calfis.o ; $(RM) calfis.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(calfis_p.o) : $(LIBF)/dynphy_lonlat/calfis_p.F \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdynphy_lonlat.a(callphysiq_mod.o) \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libdyn3d_common.a(cpdet_mod.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libdyn3dpar.a(logic_mod.o) \
$(LIBO)/libdyn3dpar.a(mod_const_mpi.o) \
$(LIBO)/libdynphy_lonlat.a(mod_interface_dyn_phys.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_mpi_data.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_omp_data.o) \
$(LIBO)/libdyn3dpar.a(moyzon_mod.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o) \
$(LIBO)/libdyn3dpar.a(times.o) \
$(LIBO)/libmisc.a(write_field.o) \
$(LIBO)/libdyn3dpar.a(write_field_p.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a calfis_p.o ; \
	$(COMPILE90) $(LIBF)/dynphy_lonlat/calfis_p.F ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a calfis_p.o ; $(RM) calfis_p.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(gr_dyn_fi.o) : $(LIBF)/dynphy_lonlat/gr_dyn_fi.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a gr_dyn_fi.o ; \
	$(COMPILE) $(LIBF)/dynphy_lonlat/gr_dyn_fi.F ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a gr_dyn_fi.o ; $(RM) gr_dyn_fi.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(gr_dyn_fi_p.o) : $(LIBF)/dynphy_lonlat/gr_dyn_fi_p.F \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libdynphy_lonlat.a(mod_interface_dyn_phys.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a gr_dyn_fi_p.o ; \
	$(COMPILE90) $(LIBF)/dynphy_lonlat/gr_dyn_fi_p.F ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a gr_dyn_fi_p.o ; $(RM) gr_dyn_fi_p.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(gr_fi_dyn.o) : $(LIBF)/dynphy_lonlat/gr_fi_dyn.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a gr_fi_dyn.o ; \
	$(COMPILE) $(LIBF)/dynphy_lonlat/gr_fi_dyn.F ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a gr_fi_dyn.o ; $(RM) gr_fi_dyn.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(gr_fi_dyn_p.o) : $(LIBF)/dynphy_lonlat/gr_fi_dyn_p.F \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libdynphy_lonlat.a(mod_interface_dyn_phys.o) \
$(LIBO)/libdyn3dpar.a(parallel_lmdz.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a gr_fi_dyn_p.o ; \
	$(COMPILE90) $(LIBF)/dynphy_lonlat/gr_fi_dyn_p.F ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a gr_fi_dyn_p.o ; $(RM) gr_fi_dyn_p.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(grid_atob_m.o) : $(LIBF)/dynphy_lonlat/grid_atob_m.F90 \
$(LIBO)/libmisc.a(assert_eq_m.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a grid_atob_m.o ; \
	$(COMPILE90) $(LIBF)/dynphy_lonlat/grid_atob_m.F90 ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a grid_atob_m.o ; $(RM) grid_atob_m.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(inigeomphy_mod.o) : $(LIBF)/dynphy_lonlat/inigeomphy_mod.F90 \
$(LIBF)/misc/iniprint.h \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libphy_common.a(geometry_mod.o) \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libdynphy_lonlat.a(mod_interface_dyn_phys.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libmisc.a(nrtype.o) \
$(LIBO)/libphy_common.a(physics_distribution_mod.o) \
$(LIBO)/libphy_common.a(regular_lonlat_mod.o) \
$(LIBO)/libphy_common.a(vertical_layers_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a inigeomphy_mod.o ; \
	$(COMPILE90) $(LIBF)/dynphy_lonlat/inigeomphy_mod.F90 ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a inigeomphy_mod.o ; $(RM) inigeomphy_mod.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(mod_interface_dyn_phys.o) : $(LIBF)/dynphy_lonlat/mod_interface_dyn_phys.F90 \
$(LIBO)/libphy_common.a(mod_phys_lmdz_mpi_data.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a mod_interface_dyn_phys.o ; \
	$(COMPILE90) $(LIBF)/dynphy_lonlat/mod_interface_dyn_phys.F90 ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a mod_interface_dyn_phys.o ; $(RM) mod_interface_dyn_phys.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory dynphy_lonlat/phystd
#=======================================================================

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(caldyn0.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(datareadnc.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(defrun_new.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(exner_hyb.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(grid_noro1.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(ini_archive.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(interp_vert.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(lect_start_archive.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(scal_wind.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(wind_scal.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(write_archive.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(callphysiq_mod.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(inichim_newstart.o)

$(LIBO)/libdynphy_lonlat.a :  $(LIBO)/libdynphy_lonlat.a(iniphysiq_mod.o)

.PRECIOUS	: $(LIBO)/libdynphy_lonlat.a


# Compilation of elements in dynphy_lonlat/phystd of library libdynphy_lonlat.a

$(LIBO)/libdynphy_lonlat.a(caldyn0.o) : $(LIBF)/dynphy_lonlat/phystd/caldyn0.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comvert_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a caldyn0.o ; \
	$(COMPILE90) $(LIBF)/dynphy_lonlat/phystd/caldyn0.F ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a caldyn0.o ; $(RM) caldyn0.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(datareadnc.o) : $(LIBF)/dynphy_lonlat/phystd/datareadnc.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libphystd.a(datafile_mod.o) \
$(LIBO)/libmisc.a(ioipsl_getincom.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a datareadnc.o ; \
	$(COMPILE90) $(LIBF)/dynphy_lonlat/phystd/datareadnc.F ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a datareadnc.o ; $(RM) datareadnc.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(defrun_new.o) : $(LIBF)/dynphy_lonlat/phystd/defrun_new.F \
$(LIBF)/dyn3d_common/comdissnew.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libmisc.a(ioipsl_getincom.o) \
$(LIBO)/libdyn3dpar.a(logic_mod.o) \
$(LIBO)/libdyn3d_common.a(serre_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a defrun_new.o ; \
	$(COMPILE90) $(LIBF)/dynphy_lonlat/phystd/defrun_new.F ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a defrun_new.o ; $(RM) defrun_new.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(exner_hyb.o) : $(LIBF)/dynphy_lonlat/phystd/exner_hyb.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a exner_hyb.o ; \
	$(COMPILE90) $(LIBF)/dynphy_lonlat/phystd/exner_hyb.F ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a exner_hyb.o ; $(RM) exner_hyb.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(grid_noro1.o) : $(LIBF)/dynphy_lonlat/phystd/grid_noro1.F \
$(LIBF)/grid/dimensions.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a grid_noro1.o ; \
	$(COMPILE90) $(LIBF)/dynphy_lonlat/phystd/grid_noro1.F ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a grid_noro1.o ; $(RM) grid_noro1.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(ini_archive.o) : $(LIBF)/dynphy_lonlat/phystd/ini_archive.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libphystd.a(slab_ice_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a ini_archive.o ; \
	$(COMPILE90) $(LIBF)/dynphy_lonlat/phystd/ini_archive.F ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a ini_archive.o ; $(RM) ini_archive.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(interp_vert.o) : $(LIBF)/dynphy_lonlat/phystd/interp_vert.F
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a interp_vert.o ; \
	$(COMPILE) $(LIBF)/dynphy_lonlat/phystd/interp_vert.F ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a interp_vert.o ; $(RM) interp_vert.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(lect_start_archive.o) : $(LIBF)/dynphy_lonlat/phystd/lect_start_archive.F \
$(LIBF)/dyn3d_common/comgeom2.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libdyn3d_common.a(comconst_mod.o) \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libphystd.a(slab_ice_h.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a lect_start_archive.o ; \
	$(COMPILE90) $(LIBF)/dynphy_lonlat/phystd/lect_start_archive.F ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a lect_start_archive.o ; $(RM) lect_start_archive.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(scal_wind.o) : $(LIBF)/dynphy_lonlat/phystd/scal_wind.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a scal_wind.o ; \
	$(COMPILE) $(LIBF)/dynphy_lonlat/phystd/scal_wind.F ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a scal_wind.o ; $(RM) scal_wind.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(wind_scal.o) : $(LIBF)/dynphy_lonlat/phystd/wind_scal.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a wind_scal.o ; \
	$(COMPILE) $(LIBF)/dynphy_lonlat/phystd/wind_scal.F ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a wind_scal.o ; $(RM) wind_scal.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(write_archive.o) : $(LIBF)/dynphy_lonlat/phystd/write_archive.F \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libphystd.a(comsoil_h.o) \
$(LIBO)/libphystd.a(slab_ice_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a write_archive.o ; \
	$(COMPILE90) $(LIBF)/dynphy_lonlat/phystd/write_archive.F ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a write_archive.o ; $(RM) write_archive.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(callphysiq_mod.o) : $(LIBF)/dynphy_lonlat/phystd/callphysiq_mod.F90 \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(physiq_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a callphysiq_mod.o ; \
	$(COMPILE90) $(LIBF)/dynphy_lonlat/phystd/callphysiq_mod.F90 ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a callphysiq_mod.o ; $(RM) callphysiq_mod.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(inichim_newstart.o) : $(LIBF)/dynphy_lonlat/phystd/inichim_newstart.F90 \
$(LIBO)/libphystd.a(callkeys_mod.o) \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libphystd.a(datafile_mod.o) \
$(LIBO)/libphy_common.a(mod_grid_phy_lmdz.o) \
$(LIBO)/libphystd.a(tracer_h.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a inichim_newstart.o ; \
	$(COMPILE90) $(LIBF)/dynphy_lonlat/phystd/inichim_newstart.F90 ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a inichim_newstart.o ; $(RM) inichim_newstart.o ; \
	cd $(GCM)

$(LIBO)/libdynphy_lonlat.a(iniphysiq_mod.o) : $(LIBF)/dynphy_lonlat/phystd/iniphysiq_mod.F90 \
$(LIBF)/dyn3d_common/comgeom.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/misc/iniprint.h \
$(LIBF)/dyn3d_common/paramet.h \
$(LIBO)/libdyn3d_common.a(comvert_mod.o) \
$(LIBO)/libdyn3d_common.a(control_mod.o) \
$(LIBO)/libphystd.a(dimphy.o) \
$(LIBO)/libphy_common.a(geometry_mod.o) \
$(LIBO)/libdyn3d_common.a(infotrac.o) \
$(LIBO)/libphystd.a(inifis_mod.o) \
$(LIBO)/libdynphy_lonlat.a(inigeomphy_mod.o) \
$(LIBO)/libphy_common.a(ioipsl_getin_p_mod.o) \
$(LIBO)/libphy_common.a(mod_phys_lmdz_para.o) \
$(LIBO)/libphystd.a(planete_mod.o) \
$(LIBO)/libphystd.a(surf_heat_transp_mod.o)
	cd $(LIBO); \
	$(AR) d $(LIBO)/libdynphy_lonlat.a iniphysiq_mod.o ; \
	$(COMPILE90) $(LIBF)/dynphy_lonlat/phystd/iniphysiq_mod.F90 ; \
	$(AR) r $(LIBO)/libdynphy_lonlat.a iniphysiq_mod.o ; $(RM) iniphysiq_mod.o ; \
	cd $(GCM)


