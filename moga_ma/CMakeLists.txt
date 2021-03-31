set (LIBNAME moga_ma_lib)
cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(moga_ma)

set(INC_DIRS 
)

set (SRC_DIRS
  modules
)

set (EXE_SPECS
  cmake.moga_ma
  cmake.error_tool
  cmake.error_tool_special
  cmake.aperture_and_lifetime
  cmake.breeder
)

#set (FFLAGS "-pg")  # For profiling
#set (CFLAGS "-pg")

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
