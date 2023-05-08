! Copyright 2013-2015,2017 Université de Reims Champagne-Ardenne 
! Contributor: J. Burgalat (GSMA, URCA)
! email of the author : jeremie.burgalat@univ-reims.fr
! 
! This software is a computer program whose purpose is to compute
! microphysics processes using a two-moments scheme.
! 
! This library is governed by the CeCILL-B license under French law and
! abiding by the rules of distribution of free software.  You can  use, 
! modify and/ or redistribute the software under the terms of the CeCILL-B
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info". 
! 
! As a counterpart to the access to the source code and  rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty  and the software's author,  the holder of the
! economic rights,  and the successive licensors  have only  limited
! liability. 
! 
! In this respect, the user's attention is drawn to the risks associated
! with loading,  using,  modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean  that it is complicated to manipulate,  and  that  also
! therefore means  that it is reserved for developers  and  experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or 
! data to be ensured and,  more generally, to use and operate it in the 
! same conditions as regards security. 
! 
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL-B license and that you accept its terms.

!! file: mm_mp2m.f90
!! summary: MP2M library interface module.
!! author: J. Burgalat
!! date: 2013-2015,2017

MODULE MM_LIB
  !! MP2M library interface module.
  !!
  !! This module is only intended to get a overall acces to all library module. It contains no 
  !! definitions and just __uses__ all others modules of the library.
  USE MM_MPREC
  USE MM_GLOBALS
  USE MM_INTERFACES
  USE MM_METHODS
  USE MM_HAZE
  USE MM_CLOUDS
  USE MM_MICROPHYSIC
END MODULE MM_LIB

