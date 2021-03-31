cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

project(custom_tao)

set (EXE_SPECS
  cmake.custom_tao
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
