#! /bin/bash

wheregen=$MOD/LMDZ.GENERIC/

mkdir phygeneric
cd phygeneric

mkdir bands
cd bands
ln -sf $wheregen/libf/grid/dimension/makbands .
cd ..


ln -sf $wheregen/libf/phystd/* .
\rm for_lmdz5
\rm ismin*
\rm ismax*
\rm rcm1d*
\rm kcm1d*

ln -sf $wheregen/libf/bibio/cbrt.F .
ln -sf $wheregen/libf/bibio/lmdstd.h .

ln -sf $wheregen/libf/dyn3d/control.h .

ln -sf $wheregen/libf/phystd/for_lmdz5/generatedoth.sh .
ln -sf $wheregen/libf/phystd/for_lmdz5/iniphysiq.F . 
ln -sf $wheregen/libf/phystd/for_lmdz5/dimphys.h .

ln -sf ../phylmd/comgeomphy.F90 .
ln -sf ../phylmd/dimphy.F90 .
ln -sf ../phylmd/init_phys_lmdz.F90 .
ln -sf ../phylmd/iophy.F90 .
ln -sf ../phylmd/iostart.F90 .
ln -sf ../phylmd/mod_grid_phy_lmdz.F90 .
ln -sf ../phylmd/mod_phys_lmdz_mpi_data.F90 .
ln -sf ../phylmd/mod_phys_lmdz_mpi_transfert.F90 .
ln -sf ../phylmd/mod_phys_lmdz_omp_data.F90 .
ln -sf ../phylmd/mod_phys_lmdz_omp_transfert.F90 .
ln -sf ../phylmd/mod_phys_lmdz_para.F90 .
ln -sf ../phylmd/mod_phys_lmdz_transfert_para.F90 .
ln -sf ../phylmd/tetalevel.F .
ln -sf ../phylmd/write_field_phy.F90 .
