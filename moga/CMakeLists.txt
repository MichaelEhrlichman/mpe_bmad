set (LIBNAME moga_lib)
cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(INC_DIRS 
)

set (SRC_DIRS
  modules
)

set (EXE_SPECS
#  cmake.moga
  cmake.ansga2
  cmake.aspea2
)

#set (FFLAGS "-pg")  # For profiling
#set (CFLAGS "-pg")

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
