cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(sls_lib)

set(LIBNAME sls_lib)

set(CREATE_SHARED false)

set(SRC_DIRS
  code
  modules
  interfaces
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)





