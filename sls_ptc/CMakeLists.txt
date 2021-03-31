cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(INC_DIRS 
)

set (SRC_DIRS
  code
  modules
)

set (EXE_SPECS
  #cmake.test_lsq_h
  #cmake.ptc_et
  #cmake.ptc_et_dis1
  cmake.ptc_mod_driver
  #cmake.ptc_co
  #cmake.ptc_emit
  cmake.ptc_map
  cmake.ptc_c_map
  #cmake.ptc_shifts
  cmake.ptc_make_flat
  #cmake.ptc_jac
  #cmake.ptc_h
  #cmake.ptc_ndts
  #cmake.ptc_response_sext
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
