#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Definitions de Macros pour Make
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Repertoires :

GCM     = /donnees/aslmd/WRF/MARS_WRF/64bits_pgf90/WRFV2/mars_lmd
MACHINE = sun4
LIBF    = $(GCM)/libf
LIBO    = $(LIBOGCM)/$(MACHINE)
LOCAL_DIR=/donnees/aslmd/WRF/MARS_WRF/64bits_pgf90/WRFV2/mars_lmd
BIBIO    = $(LIBF)/bibio
AERONOMARS    = $(LIBF)/aeronomars
FILTRE   = filtre
PHYS  = 
DYN  = dyn 
LIBPHY = $(LIBO)/libphy$(PHYS).a
DIRMAIN=dyn$(DIM)d
RM=rm

OPLINK = 

# Les differentes librairies pour l"edition des liens:

dyn3d      = $(LIBO)/libdyn3d.a $(LIBO)/lib$(FILTRE).a
dyn2d      = $(LIBO)/libdyn2d.a
dyn1d      = $(LIBO)/libdyn1d.a
L_DYN      = -ldyn$(DIM)d
L_FILTRE   = -l$(FILTRE)
L_PHY = -lphy$(PHYS) 
L_BIBIO    = -lbibio
L_AERONOMARS    = -laeronomars
L_ADJNT    =

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Option de compilation FORTRAN
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COMPILE90 = $(F90) $(OPTIM90) $(INCLUDE) -c
LINK    = 
AR      = 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Creation des differents executables
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Executables:
# ------------

PROG = code

main : $(DYN) bibio phys aeronomars $(OPTION_DEP) 
	cd $(LIBO) ; $(RANLIB) lib*.a ; cd $(GCM) ;\
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/$(DIRMAIN)/$(PROG).F -o $(PROG).o ; \
	$(LINK) $(PROG).o -L$(LIBO) $(L_DYN) $(L_ADJNT) $(L_FILTRE) $(L_AERONOMARS) $(L_PHY) $(L_BIBIO) $(L_DYN) $(L_AERONOMARS) $(OPLINK) $(OPTION_LINK) -o $(LOCAL_DIR)/$(PROG).e ; $(RM) $(PROG).o 

dyn : $(LIBO)/libdyn$(DIM)d.a $(FILTRE)$(DIM)d

phys : $(LIBPHY)

bibio : $(LIBO)/libbibio.a

aeronomars : $(LIBO)/libaeronomars.a

adjnt : $(LIBO)/libadjnt.a

$(FILTRE)3d : $(LIBO)/lib$(FILTRE).a

$(FILTRE)2d :

$(FILTRE)1d :

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Contenu des differentes bibliotheques
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory aeronomars
#=======================================================================

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(calchim.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(chemthermos.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(ch.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(column.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(concentrations.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(conduction.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(dtridgl.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(euvheat.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(flujo.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(getch.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(hrtherm.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(inichim.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(inichim_newstart.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(inichim_readcallphys.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(init_chimie_B.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(interpfast.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(intrplf.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(inv.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(jthermcalc.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(moldiffcoeff.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(moldiff.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(molvis.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(paramfoto.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(param_read.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(perosat.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(phdisrate.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(photochemist_B.o)

$(LIBO)/libaeronomars.a :  $(LIBO)/libaeronomars.a(thermosphere.o)

.PRECIOUS	: $(LIBO)/libaeronomars.a


# Compilation des membres de la bibliotheque libaeronomars.a

$(LIBO)/libaeronomars.a(calchim.o) : $(LIBF)/aeronomars/calchim.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/fisice.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/calchim.F ; \
	$(AR) r $(LIBO)/libaeronomars.a calchim.o ; $(RM) calchim.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(chemthermos.o) : $(LIBF)/aeronomars/chemthermos.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/aeronomars/param.h \
$(LIBF)/aeronomars/param_v3.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/chemthermos.F ; \
	$(AR) r $(LIBO)/libaeronomars.a chemthermos.o ; $(RM) chemthermos.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(ch.o) : $(LIBF)/aeronomars/ch.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/ch.F ; \
	$(AR) r $(LIBO)/libaeronomars.a ch.o ; $(RM) ch.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(column.o) : $(LIBF)/aeronomars/column.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/column.F ; \
	$(AR) r $(LIBO)/libaeronomars.a column.o ; $(RM) column.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(concentrations.o) : $(LIBF)/aeronomars/concentrations.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/concentrations.F ; \
	$(AR) r $(LIBO)/libaeronomars.a concentrations.o ; $(RM) concentrations.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(conduction.o) : $(LIBF)/aeronomars/conduction.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/conduction.F ; \
	$(AR) r $(LIBO)/libaeronomars.a conduction.o ; $(RM) conduction.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(dtridgl.o) : $(LIBF)/aeronomars/dtridgl.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/dtridgl.F ; \
	$(AR) r $(LIBO)/libaeronomars.a dtridgl.o ; $(RM) dtridgl.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(euvheat.o) : $(LIBF)/aeronomars/euvheat.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/aeronomars/param.h \
$(LIBF)/aeronomars/param_v3.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/euvheat.F ; \
	$(AR) r $(LIBO)/libaeronomars.a euvheat.o ; $(RM) euvheat.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(flujo.o) : $(LIBF)/aeronomars/flujo.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/flujo.F ; \
	$(AR) r $(LIBO)/libaeronomars.a flujo.o ; $(RM) flujo.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(getch.o) : $(LIBF)/aeronomars/getch.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/getch.F ; \
	$(AR) r $(LIBO)/libaeronomars.a getch.o ; $(RM) getch.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(hrtherm.o) : $(LIBF)/aeronomars/hrtherm.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/hrtherm.F ; \
	$(AR) r $(LIBO)/libaeronomars.a hrtherm.o ; $(RM) hrtherm.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(inichim.o) : $(LIBF)/aeronomars/inichim.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/datafile.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/inichim.F ; \
	$(AR) r $(LIBO)/libaeronomars.a inichim.o ; $(RM) inichim.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(inichim_newstart.o) : $(LIBF)/aeronomars/inichim_newstart.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/phymars/datafile.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/inichim_newstart.F ; \
	$(AR) r $(LIBO)/libaeronomars.a inichim_newstart.o ; $(RM) inichim_newstart.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(inichim_readcallphys.o) : $(LIBF)/aeronomars/inichim_readcallphys.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/planete.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/comsaison.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/surfdat.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/inichim_readcallphys.F ; \
	$(AR) r $(LIBO)/libaeronomars.a inichim_readcallphys.o ; $(RM) inichim_readcallphys.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(init_chimie_B.o) : $(LIBF)/aeronomars/init_chimie_B.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/datafile.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/init_chimie_B.F ; \
	$(AR) r $(LIBO)/libaeronomars.a init_chimie_B.o ; $(RM) init_chimie_B.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(interpfast.o) : $(LIBF)/aeronomars/interpfast.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/interpfast.F ; \
	$(AR) r $(LIBO)/libaeronomars.a interpfast.o ; $(RM) interpfast.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(intrplf.o) : $(LIBF)/aeronomars/intrplf.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/intrplf.F ; \
	$(AR) r $(LIBO)/libaeronomars.a intrplf.o ; $(RM) intrplf.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(inv.o) : $(LIBF)/aeronomars/inv.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/inv.F ; \
	$(AR) r $(LIBO)/libaeronomars.a inv.o ; $(RM) inv.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(jthermcalc.o) : $(LIBF)/aeronomars/jthermcalc.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/jthermcalc.F ; \
	$(AR) r $(LIBO)/libaeronomars.a jthermcalc.o ; $(RM) jthermcalc.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(moldiffcoeff.o) : $(LIBF)/aeronomars/moldiffcoeff.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/moldiffcoeff.F ; \
	$(AR) r $(LIBO)/libaeronomars.a moldiffcoeff.o ; $(RM) moldiffcoeff.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(moldiff.o) : $(LIBF)/aeronomars/moldiff.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/moldiff.F ; \
	$(AR) r $(LIBO)/libaeronomars.a moldiff.o ; $(RM) moldiff.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(molvis.o) : $(LIBF)/aeronomars/molvis.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/molvis.F ; \
	$(AR) r $(LIBO)/libaeronomars.a molvis.o ; $(RM) molvis.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(paramfoto.o) : $(LIBF)/aeronomars/paramfoto.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/paramfoto.F ; \
	$(AR) r $(LIBO)/libaeronomars.a paramfoto.o ; $(RM) paramfoto.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(param_read.o) : $(LIBF)/aeronomars/param_read.F \
$(LIBF)/aeronomars/param.h \
$(LIBF)/aeronomars/param_v3.h \
$(LIBF)/phymars/datafile.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/param_read.F ; \
	$(AR) r $(LIBO)/libaeronomars.a param_read.o ; $(RM) param_read.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(perosat.o) : $(LIBF)/aeronomars/perosat.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/perosat.F ; \
	$(AR) r $(LIBO)/libaeronomars.a perosat.o ; $(RM) perosat.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(phdisrate.o) : $(LIBF)/aeronomars/phdisrate.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/phdisrate.F ; \
	$(AR) r $(LIBO)/libaeronomars.a phdisrate.o ; $(RM) phdisrate.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(photochemist_B.o) : $(LIBF)/aeronomars/photochemist_B.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/photochemist_B.F ; \
	$(AR) r $(LIBO)/libaeronomars.a photochemist_B.o ; $(RM) photochemist_B.o ; \
	cd $(GCM)

$(LIBO)/libaeronomars.a(thermosphere.o) : $(LIBF)/aeronomars/thermosphere.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/aeronomars/param.h \
$(LIBF)/aeronomars/param_v3.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/aeronomars/thermosphere.F ; \
	$(AR) r $(LIBO)/libaeronomars.a thermosphere.o ; $(RM) thermosphere.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory bibio
#=======================================================================

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(cbrt.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(cvmgp.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(cvmgt.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(from36.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(gather.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ilsum.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ini36.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(isamax.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(isitor.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ismax.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ismin.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(lnblnk.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(mxm.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(mxva.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(mxv.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(mywrite.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(readcoord.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(readstd.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(scatter.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(scopy.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(scopyi.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(sdot.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(ssum.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(to36.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(uniplot.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(wheneq.o)

$(LIBO)/libbibio.a :  $(LIBO)/libbibio.a(writestd.o)

.PRECIOUS	: $(LIBO)/libbibio.a


# Compilation des membres de la bibliotheque libbibio.a

$(LIBO)/libbibio.a(cbrt.o) : $(LIBF)/bibio/cbrt.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/cbrt.F ; \
	$(AR) r $(LIBO)/libbibio.a cbrt.o ; $(RM) cbrt.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(cvmgp.o) : $(LIBF)/bibio/cvmgp.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/cvmgp.F ; \
	$(AR) r $(LIBO)/libbibio.a cvmgp.o ; $(RM) cvmgp.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(cvmgt.o) : $(LIBF)/bibio/cvmgt.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/cvmgt.F ; \
	$(AR) r $(LIBO)/libbibio.a cvmgt.o ; $(RM) cvmgt.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(from36.o) : $(LIBF)/bibio/from36.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/from36.F ; \
	$(AR) r $(LIBO)/libbibio.a from36.o ; $(RM) from36.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(gather.o) : $(LIBF)/bibio/gather.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/gather.F ; \
	$(AR) r $(LIBO)/libbibio.a gather.o ; $(RM) gather.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ilsum.o) : $(LIBF)/bibio/ilsum.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/ilsum.F ; \
	$(AR) r $(LIBO)/libbibio.a ilsum.o ; $(RM) ilsum.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ini36.o) : $(LIBF)/bibio/ini36.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/ini36.F ; \
	$(AR) r $(LIBO)/libbibio.a ini36.o ; $(RM) ini36.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(isamax.o) : $(LIBF)/bibio/isamax.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/isamax.F ; \
	$(AR) r $(LIBO)/libbibio.a isamax.o ; $(RM) isamax.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(isitor.o) : $(LIBF)/bibio/isitor.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/isitor.F ; \
	$(AR) r $(LIBO)/libbibio.a isitor.o ; $(RM) isitor.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ismax.o) : $(LIBF)/bibio/ismax.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/ismax.F ; \
	$(AR) r $(LIBO)/libbibio.a ismax.o ; $(RM) ismax.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ismin.o) : $(LIBF)/bibio/ismin.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/ismin.F ; \
	$(AR) r $(LIBO)/libbibio.a ismin.o ; $(RM) ismin.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(lnblnk.o) : $(LIBF)/bibio/lnblnk.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/lnblnk.F ; \
	$(AR) r $(LIBO)/libbibio.a lnblnk.o ; $(RM) lnblnk.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(mxm.o) : $(LIBF)/bibio/mxm.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/mxm.F ; \
	$(AR) r $(LIBO)/libbibio.a mxm.o ; $(RM) mxm.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(mxva.o) : $(LIBF)/bibio/mxva.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/mxva.F ; \
	$(AR) r $(LIBO)/libbibio.a mxva.o ; $(RM) mxva.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(mxv.o) : $(LIBF)/bibio/mxv.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/mxv.F ; \
	$(AR) r $(LIBO)/libbibio.a mxv.o ; $(RM) mxv.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(mywrite.o) : $(LIBF)/bibio/mywrite.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/mywrite.F ; \
	$(AR) r $(LIBO)/libbibio.a mywrite.o ; $(RM) mywrite.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(readcoord.o) : $(LIBF)/bibio/readcoord.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/readcoord.F ; \
	$(AR) r $(LIBO)/libbibio.a readcoord.o ; $(RM) readcoord.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(readstd.o) : $(LIBF)/bibio/readstd.F \
$(LIBF)/bibio/lmdstd.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/readstd.F ; \
	$(AR) r $(LIBO)/libbibio.a readstd.o ; $(RM) readstd.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(scatter.o) : $(LIBF)/bibio/scatter.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/scatter.F ; \
	$(AR) r $(LIBO)/libbibio.a scatter.o ; $(RM) scatter.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(scopy.o) : $(LIBF)/bibio/scopy.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/scopy.F ; \
	$(AR) r $(LIBO)/libbibio.a scopy.o ; $(RM) scopy.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(scopyi.o) : $(LIBF)/bibio/scopyi.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/scopyi.F ; \
	$(AR) r $(LIBO)/libbibio.a scopyi.o ; $(RM) scopyi.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(sdot.o) : $(LIBF)/bibio/sdot.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/sdot.F ; \
	$(AR) r $(LIBO)/libbibio.a sdot.o ; $(RM) sdot.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(ssum.o) : $(LIBF)/bibio/ssum.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/ssum.F ; \
	$(AR) r $(LIBO)/libbibio.a ssum.o ; $(RM) ssum.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(to36.o) : $(LIBF)/bibio/to36.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/to36.F ; \
	$(AR) r $(LIBO)/libbibio.a to36.o ; $(RM) to36.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(uniplot.o) : $(LIBF)/bibio/uniplot.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/uniplot.F ; \
	$(AR) r $(LIBO)/libbibio.a uniplot.o ; $(RM) uniplot.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(wheneq.o) : $(LIBF)/bibio/wheneq.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/wheneq.F ; \
	$(AR) r $(LIBO)/libbibio.a wheneq.o ; $(RM) wheneq.o ; \
	cd $(GCM)

$(LIBO)/libbibio.a(writestd.o) : $(LIBF)/bibio/writestd.F \
$(LIBF)/bibio/lmdstd.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/bibio/writestd.F ; \
	$(AR) r $(LIBO)/libbibio.a writestd.o ; $(RM) writestd.o ; \
	cd $(GCM)





#=======================================================================
# Contenu de la bibliotheque correspondant au Directory phymars
#=======================================================================

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(aerave.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(blackl.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(blendrad.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(calldrag_noro.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(callradite.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(callsedim2q.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(callsedim.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(co2snow.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(convadj.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(cvmgp.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(cvmgt.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(datareadnc.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(dedd.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(drag_noro.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(dustdevil.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(dustlift.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(eofdump.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(flusv.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(growthrate.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(gwprofil.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(gwstress.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(ineofdump.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(iniorbit.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(inistats.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(initracer.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(ismax.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(ismin.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwb.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwdiff.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwflux.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwi.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwmain.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwtt.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwu.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwxb.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwxd.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(lwxn.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(meso_dustopacity.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(meso_inifis.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(meso_physiq.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(meso_writediagfi.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(mucorr.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(mufract.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(multipl.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(newcondens.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(newsedim.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(nirco2abs.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(nltecool.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(nlthermeq.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(orbite.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(orodrag.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(orosetup.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(phyetat0.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(physdem1.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(profile.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(scatter.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(scopyi.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(soil.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(solang.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(solarlong.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(suaer.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(sugwd.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(sulw.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(surfini.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(swmain.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(swrayleigh.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(swr.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(tabfi.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(uniplot.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(vdif_cd.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(vdifc.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(vdif_kc.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(vlz_fi.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(watercloud.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(watersat.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(writeg1d.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(zerophys.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(def_var.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(meso_readtesassim.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(mkstat.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(readtesassim.o)

$(LIBO)/libphymars.a :  $(LIBO)/libphymars.a(wstats.o)

.PRECIOUS	: $(LIBO)/libphymars.a


# Compilation des membres de la bibliotheque libphymars.a

$(LIBO)/libphymars.a(aerave.o) : $(LIBF)/phymars/aerave.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/aerave.F ; \
	$(AR) r $(LIBO)/libphymars.a aerave.o ; $(RM) aerave.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(blackl.o) : $(LIBF)/phymars/blackl.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/blackl.F ; \
	$(AR) r $(LIBO)/libphymars.a blackl.o ; $(RM) blackl.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(blendrad.o) : $(LIBF)/phymars/blendrad.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/nlteparams.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/blendrad.F ; \
	$(AR) r $(LIBO)/libphymars.a blendrad.o ; $(RM) blendrad.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(calldrag_noro.o) : $(LIBF)/phymars/calldrag_noro.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/surfdat.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/calldrag_noro.F ; \
	$(AR) r $(LIBO)/libphymars.a calldrag_noro.o ; $(RM) calldrag_noro.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(callradite.o) : $(LIBF)/phymars/callradite.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/callradite.F ; \
	$(AR) r $(LIBO)/libphymars.a callradite.o ; $(RM) callradite.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(callsedim2q.o) : $(LIBF)/phymars/callsedim2q.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/tracer.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/callsedim2q.F ; \
	$(AR) r $(LIBO)/libphymars.a callsedim2q.o ; $(RM) callsedim2q.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(callsedim.o) : $(LIBF)/phymars/callsedim.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/fisice.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/callsedim.F ; \
	$(AR) r $(LIBO)/libphymars.a callsedim.o ; $(RM) callsedim.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(co2snow.o) : $(LIBF)/phymars/co2snow.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/co2snow.F ; \
	$(AR) r $(LIBO)/libphymars.a co2snow.o ; $(RM) co2snow.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(convadj.o) : $(LIBF)/phymars/convadj.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/tracer.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/convadj.F ; \
	$(AR) r $(LIBO)/libphymars.a convadj.o ; $(RM) convadj.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(cvmgp.o) : $(LIBF)/phymars/cvmgp.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/cvmgp.F ; \
	$(AR) r $(LIBO)/libphymars.a cvmgp.o ; $(RM) cvmgp.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(cvmgt.o) : $(LIBF)/phymars/cvmgt.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/cvmgt.F ; \
	$(AR) r $(LIBO)/libphymars.a cvmgt.o ; $(RM) cvmgt.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(datareadnc.o) : $(LIBF)/phymars/datareadnc.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/phymars/datafile.h \
$(LIBF)/bibio/lmdstd.h \
$(LIBF)/grid/fxyprim.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/datareadnc.F ; \
	$(AR) r $(LIBO)/libphymars.a datareadnc.o ; $(RM) datareadnc.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(dedd.o) : $(LIBF)/phymars/dedd.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/dedd.F ; \
	$(AR) r $(LIBO)/libphymars.a dedd.o ; $(RM) dedd.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(drag_noro.o) : $(LIBF)/phymars/drag_noro.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/comcstfi.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/drag_noro.F ; \
	$(AR) r $(LIBO)/libphymars.a drag_noro.o ; $(RM) drag_noro.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(dustdevil.o) : $(LIBF)/phymars/dustdevil.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/phymars/planete.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/tracer.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/dustdevil.F ; \
	$(AR) r $(LIBO)/libphymars.a dustdevil.o ; $(RM) dustdevil.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(dustlift.o) : $(LIBF)/phymars/dustlift.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/tracer.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/dustlift.F ; \
	$(AR) r $(LIBO)/libphymars.a dustlift.o ; $(RM) dustlift.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(eofdump.o) : $(LIBF)/phymars/eofdump.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/eofdump.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/eofdump.F ; \
	$(AR) r $(LIBO)/libphymars.a eofdump.o ; $(RM) eofdump.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(flusv.o) : $(LIBF)/phymars/flusv.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/flusv.F ; \
	$(AR) r $(LIBO)/libphymars.a flusv.o ; $(RM) flusv.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(growthrate.o) : $(LIBF)/phymars/growthrate.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/growthrate.F ; \
	$(AR) r $(LIBO)/libphymars.a growthrate.o ; $(RM) growthrate.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(gwprofil.o) : $(LIBF)/phymars/gwprofil.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yoegwd.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/gwprofil.F ; \
	$(AR) r $(LIBO)/libphymars.a gwprofil.o ; $(RM) gwprofil.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(gwstress.o) : $(LIBF)/phymars/gwstress.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yoegwd.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/gwstress.F ; \
	$(AR) r $(LIBO)/libphymars.a gwstress.o ; $(RM) gwstress.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(ineofdump.o) : $(LIBF)/phymars/ineofdump.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/eofdump.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/ineofdump.F ; \
	$(AR) r $(LIBO)/libphymars.a ineofdump.o ; $(RM) ineofdump.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(iniorbit.o) : $(LIBF)/phymars/iniorbit.F \
$(LIBF)/phymars/planete.h \
$(LIBF)/phymars/comcstfi.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/iniorbit.F ; \
	$(AR) r $(LIBO)/libphymars.a iniorbit.o ; $(RM) iniorbit.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(inistats.o) : $(LIBF)/phymars/inistats.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comconst.h \
$(LIBF)/phymars/statto.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/inistats.F ; \
	$(AR) r $(LIBO)/libphymars.a inistats.o ; $(RM) inistats.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(initracer.o) : $(LIBF)/phymars/initracer.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/watercap.h \
$(LIBF)/phymars/aerice.h \
$(LIBF)/phymars/fisice.h \
$(LIBF)/aeronomars/chimiedata.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/initracer.F ; \
	$(AR) r $(LIBO)/libphymars.a initracer.o ; $(RM) initracer.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(ismax.o) : $(LIBF)/phymars/ismax.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/ismax.F ; \
	$(AR) r $(LIBO)/libphymars.a ismax.o ; $(RM) ismax.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(ismin.o) : $(LIBF)/phymars/ismin.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/ismin.F ; \
	$(AR) r $(LIBO)/libphymars.a ismin.o ; $(RM) ismin.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwb.o) : $(LIBF)/phymars/lwb.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/lwb.F ; \
	$(AR) r $(LIBO)/libphymars.a lwb.o ; $(RM) lwb.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwdiff.o) : $(LIBF)/phymars/lwdiff.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/yomaer.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/lwdiff.F ; \
	$(AR) r $(LIBO)/libphymars.a lwdiff.o ; $(RM) lwdiff.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwflux.o) : $(LIBF)/phymars/lwflux.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comg1d.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/lwflux.F ; \
	$(AR) r $(LIBO)/libphymars.a lwflux.o ; $(RM) lwflux.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwi.o) : $(LIBF)/phymars/lwi.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/comg1d.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/lwi.F ; \
	$(AR) r $(LIBO)/libphymars.a lwi.o ; $(RM) lwi.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwmain.o) : $(LIBF)/phymars/lwmain.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comg1d.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/lwmain.F ; \
	$(AR) r $(LIBO)/libphymars.a lwmain.o ; $(RM) lwmain.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwtt.o) : $(LIBF)/phymars/lwtt.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/lwtt.F ; \
	$(AR) r $(LIBO)/libphymars.a lwtt.o ; $(RM) lwtt.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwu.o) : $(LIBF)/phymars/lwu.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/yomaer.h \
$(LIBF)/phymars/yomlw.h \
$(LIBF)/phymars/fisice.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/aerice.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/lwu.F ; \
	$(AR) r $(LIBO)/libphymars.a lwu.o ; $(RM) lwu.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwxb.o) : $(LIBF)/phymars/lwxb.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/lwxb.F ; \
	$(AR) r $(LIBO)/libphymars.a lwxb.o ; $(RM) lwxb.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwxd.o) : $(LIBF)/phymars/lwxd.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomlw.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/lwxd.F ; \
	$(AR) r $(LIBO)/libphymars.a lwxd.o ; $(RM) lwxd.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(lwxn.o) : $(LIBF)/phymars/lwxn.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomlw.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/lwxn.F ; \
	$(AR) r $(LIBO)/libphymars.a lwxn.o ; $(RM) lwxn.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(meso_dustopacity.o) : $(LIBF)/phymars/meso_dustopacity.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomaer.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/planete.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/meso_dustopacity.F ; \
	$(AR) r $(LIBO)/libphymars.a meso_dustopacity.o ; $(RM) meso_dustopacity.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(meso_inifis.o) : $(LIBF)/phymars/meso_inifis.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/planete.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/comsaison.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/surfdat.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/meso_inifis.F ; \
	$(AR) r $(LIBO)/libphymars.a meso_inifis.o ; $(RM) meso_inifis.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(meso_physiq.o) : $(LIBF)/phymars/meso_physiq.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/comdiurn.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/planete.h \
$(LIBF)/phymars/comsaison.h \
$(LIBF)/dyn3d/control.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/comg1d.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/nlteparams.h \
$(LIBF)/aeronomars/chimiedata.h \
$(LIBF)/phymars/watercap.h \
$(LIBF)/phymars/fisice.h \
$(LIBF)/aeronomars/param.h \
$(LIBF)/aeronomars/param_v3.h \
$(LIBF)/aeronomars/conc.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/meso_physiq.F ; \
	$(AR) r $(LIBO)/libphymars.a meso_physiq.o ; $(RM) meso_physiq.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(meso_writediagfi.o) : $(LIBF)/phymars/meso_writediagfi.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/control.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom.h \
$(LIBF)/dyn3d/description.h \
$(LIBF)/dyn3d/temps.h \
$(LIBF)/phymars/surfdat.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/meso_writediagfi.F ; \
	$(AR) r $(LIBO)/libphymars.a meso_writediagfi.o ; $(RM) meso_writediagfi.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(mucorr.o) : $(LIBF)/phymars/mucorr.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/mucorr.F ; \
	$(AR) r $(LIBO)/libphymars.a mucorr.o ; $(RM) mucorr.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(mufract.o) : $(LIBF)/phymars/mufract.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/mufract.F ; \
	$(AR) r $(LIBO)/libphymars.a mufract.o ; $(RM) mufract.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(multipl.o) : $(LIBF)/phymars/multipl.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/multipl.F ; \
	$(AR) r $(LIBO)/libphymars.a multipl.o ; $(RM) multipl.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(newcondens.o) : $(LIBF)/phymars/newcondens.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/fisice.h \
$(LIBF)/grid/dimensions.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/newcondens.F ; \
	$(AR) r $(LIBO)/libphymars.a newcondens.o ; $(RM) newcondens.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(newsedim.o) : $(LIBF)/phymars/newsedim.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/newsedim.F ; \
	$(AR) r $(LIBO)/libphymars.a newsedim.o ; $(RM) newsedim.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(nirco2abs.o) : $(LIBF)/phymars/nirco2abs.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/comdiurn.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/nirco2abs.F ; \
	$(AR) r $(LIBO)/libphymars.a nirco2abs.o ; $(RM) nirco2abs.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(nltecool.o) : $(LIBF)/phymars/nltecool.F \
$(LIBF)/phymars/nltedata.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/nltecool.F ; \
	$(AR) r $(LIBO)/libphymars.a nltecool.o ; $(RM) nltecool.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(nlthermeq.o) : $(LIBF)/phymars/nlthermeq.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/nlteparams.h \
$(LIBF)/phymars/yomlw.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/nlthermeq.F ; \
	$(AR) r $(LIBO)/libphymars.a nlthermeq.o ; $(RM) nlthermeq.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(orbite.o) : $(LIBF)/phymars/orbite.F \
$(LIBF)/phymars/planete.h \
$(LIBF)/phymars/comcstfi.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/orbite.F ; \
	$(AR) r $(LIBO)/libphymars.a orbite.o ; $(RM) orbite.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(orodrag.o) : $(LIBF)/phymars/orodrag.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/yoegwd.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/orodrag.F ; \
	$(AR) r $(LIBO)/libphymars.a orodrag.o ; $(RM) orodrag.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(orosetup.o) : $(LIBF)/phymars/orosetup.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/yoegwd.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/orosetup.F ; \
	$(AR) r $(LIBO)/libphymars.a orosetup.o ; $(RM) orosetup.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(phyetat0.o) : $(LIBF)/phymars/phyetat0.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/planete.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomaer.h \
$(LIBF)/phymars/comcstfi.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/phyetat0.F ; \
	$(AR) r $(LIBO)/libphymars.a phyetat0.o ; $(RM) phyetat0.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(physdem1.o) : $(LIBF)/phymars/physdem1.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/dyn3d/paramet.h \
$(LIBF)/dyn3d/comvert.h \
$(LIBF)/dyn3d/comgeom2.h \
$(LIBF)/dyn3d/control.h \
$(LIBF)/dyn3d/comdissnew.h \
$(LIBF)/dyn3d/logic.h \
$(LIBF)/dyn3d/ener.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/dyn3d/serre.h \
$(LIBF)/dyn3d/clesph0.h \
$(LIBF)/grid/fxyprim.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/planete.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomaer.h \
$(LIBF)/phymars/comcstfi.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/physdem1.F ; \
	$(AR) r $(LIBO)/libphymars.a physdem1.o ; $(RM) physdem1.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(profile.o) : $(LIBF)/phymars/profile.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/profile.F ; \
	$(AR) r $(LIBO)/libphymars.a profile.o ; $(RM) profile.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(scatter.o) : $(LIBF)/phymars/scatter.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/scatter.F ; \
	$(AR) r $(LIBO)/libphymars.a scatter.o ; $(RM) scatter.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(scopyi.o) : $(LIBF)/phymars/scopyi.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/scopyi.F ; \
	$(AR) r $(LIBO)/libphymars.a scopyi.o ; $(RM) scopyi.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(soil.o) : $(LIBF)/phymars/soil.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/soil.F ; \
	$(AR) r $(LIBO)/libphymars.a soil.o ; $(RM) soil.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(solang.o) : $(LIBF)/phymars/solang.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/solang.F ; \
	$(AR) r $(LIBO)/libphymars.a solang.o ; $(RM) solang.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(solarlong.o) : $(LIBF)/phymars/solarlong.F \
$(LIBF)/phymars/planete.h \
$(LIBF)/phymars/comcstfi.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/solarlong.F ; \
	$(AR) r $(LIBO)/libphymars.a solarlong.o ; $(RM) solarlong.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(suaer.o) : $(LIBF)/phymars/suaer.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomaer.h \
$(LIBF)/phymars/aerdust.h \
$(LIBF)/phymars/aerice.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/suaer.F ; \
	$(AR) r $(LIBO)/libphymars.a suaer.o ; $(RM) suaer.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(sugwd.o) : $(LIBF)/phymars/sugwd.F \
$(LIBF)/phymars/yoegwd.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/sugwd.F ; \
	$(AR) r $(LIBO)/libphymars.a sugwd.o ; $(RM) sugwd.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(sulw.o) : $(LIBF)/phymars/sulw.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/yomlw.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/sulw.F ; \
	$(AR) r $(LIBO)/libphymars.a sulw.o ; $(RM) sulw.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(surfini.o) : $(LIBF)/phymars/surfini.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/tracer.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/surfini.F ; \
	$(AR) r $(LIBO)/libphymars.a surfini.o ; $(RM) surfini.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(swmain.o) : $(LIBF)/phymars/swmain.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomaer.h \
$(LIBF)/phymars/yomlw.h \
$(LIBF)/phymars/callkeys.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/swmain.F ; \
	$(AR) r $(LIBO)/libphymars.a swmain.o ; $(RM) swmain.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(swrayleigh.o) : $(LIBF)/phymars/swrayleigh.F \
$(LIBF)/phymars/comcstfi.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/swrayleigh.F ; \
	$(AR) r $(LIBO)/libphymars.a swrayleigh.o ; $(RM) swrayleigh.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(swr.o) : $(LIBF)/phymars/swr.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/yomaer.h \
$(LIBF)/phymars/yomlw.h \
$(LIBF)/phymars/fisice.h \
$(LIBF)/phymars/aerice.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/swr.F ; \
	$(AR) r $(LIBO)/libphymars.a swr.o ; $(RM) swr.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(tabfi.o) : $(LIBF)/phymars/tabfi.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/planete.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/dimradmars.h \
$(LIBF)/phymars/yomaer.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/tabfi.F ; \
	$(AR) r $(LIBO)/libphymars.a tabfi.o ; $(RM) tabfi.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(uniplot.o) : $(LIBF)/phymars/uniplot.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/uniplot.F ; \
	$(AR) r $(LIBO)/libphymars.a uniplot.o ; $(RM) uniplot.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(vdif_cd.o) : $(LIBF)/phymars/vdif_cd.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/vdif_cd.F ; \
	$(AR) r $(LIBO)/libphymars.a vdif_cd.o ; $(RM) vdif_cd.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(vdifc.o) : $(LIBF)/phymars/vdifc.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/surfdat.h \
$(LIBF)/phymars/comgeomfi.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/watercap.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/vdifc.F ; \
	$(AR) r $(LIBO)/libphymars.a vdifc.o ; $(RM) vdifc.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(vdif_kc.o) : $(LIBF)/phymars/vdif_kc.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/vdif_kc.F ; \
	$(AR) r $(LIBO)/libphymars.a vdif_kc.o ; $(RM) vdif_kc.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(vlz_fi.o) : $(LIBF)/phymars/vlz_fi.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/vlz_fi.F ; \
	$(AR) r $(LIBO)/libphymars.a vlz_fi.o ; $(RM) vlz_fi.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(watercloud.o) : $(LIBF)/phymars/watercloud.F \
$(LIBF)/grid/dimensions.h \
$(LIBF)/phymars/dimphys.h \
$(LIBF)/phymars/comcstfi.h \
$(LIBF)/phymars/callkeys.h \
$(LIBF)/phymars/tracer.h \
$(LIBF)/phymars/fisice.h \
$(LIBF)/phymars/comgeomfi.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/watercloud.F ; \
	$(AR) r $(LIBO)/libphymars.a watercloud.o ; $(RM) watercloud.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(watersat.o) : $(LIBF)/phymars/watersat.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/watersat.F ; \
	$(AR) r $(LIBO)/libphymars.a watersat.o ; $(RM) watersat.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(writeg1d.o) : $(LIBF)/phymars/writeg1d.F \
$(LIBF)/phymars/comg1d.h \
$(LIBF)/phymars/comg1d.h \
$(LIBF)/phymars/comcstfi.h
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/writeg1d.F ; \
	$(AR) r $(LIBO)/libphymars.a writeg1d.o ; $(RM) writeg1d.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(zerophys.o) : $(LIBF)/phymars/zerophys.F
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/zerophys.F ; \
	$(AR) r $(LIBO)/libphymars.a zerophys.o ; $(RM) zerophys.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(def_var.o) : $(LIBF)/phymars/def_var.F90
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/def_var.F90 ; \
	$(AR) r $(LIBO)/libphymars.a def_var.o ; $(RM) def_var.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(meso_readtesassim.o) : $(LIBF)/phymars/meso_readtesassim.F90
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/meso_readtesassim.F90 ; \
	$(AR) r $(LIBO)/libphymars.a meso_readtesassim.o ; $(RM) meso_readtesassim.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(mkstat.o) : $(LIBF)/phymars/mkstat.F90
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/mkstat.F90 ; \
	$(AR) r $(LIBO)/libphymars.a mkstat.o ; $(RM) mkstat.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(readtesassim.o) : $(LIBF)/phymars/readtesassim.F90
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/readtesassim.F90 ; \
	$(AR) r $(LIBO)/libphymars.a readtesassim.o ; $(RM) readtesassim.o ; \
	cd $(GCM)

$(LIBO)/libphymars.a(wstats.o) : $(LIBF)/phymars/wstats.F90
	cd $(LOCAL_DIR); \
	$(COMPILE90) $(LIBF)/phymars/wstats.F90 ; \
	$(AR) r $(LIBO)/libphymars.a wstats.o ; $(RM) wstats.o ; \
	cd $(GCM)


