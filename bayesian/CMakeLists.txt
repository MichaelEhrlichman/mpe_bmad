cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})
project(bayes)

SET (EXE_SPECS
  cmake.bayes
  )

SET (SHARED_DEPS
  bmad
  sim_utils
  ${PLOT_LINK_LIBS}
  $ENV{ACC_BMAD_LINK_LIBS}
  gomp
  X11
  readline
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
