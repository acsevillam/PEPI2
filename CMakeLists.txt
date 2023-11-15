#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(PEPI2)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(PEPI2 PEPI2.cc ${sources} ${headers})
target_link_libraries(PEPI2 ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build PEPI2. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(PEPI2_SCRIPTS
  mac/init_vis.mac
  mac/run1.mac
  mac/run2.mac
  mac/vis.mac
  mac/tsg_offscreen.mac
  detmask_config.in
  gps.in
  geometry_config.in
  data/Energy.txt
  data/Graphite_delta.txt
  data/Nylon_delta.txt
  data/PMMA_delta.txt
  data/Polycharbonate_delta.txt
  data/Silicon_delta.txt
  data/Water_delta.txt
  spectra/40kV42deg1Air0.5Be0.5Al.in
  spectra/50kV42deg1Air0.5Be0.5Al.in
  spectra/60kV42deg1Air0.5Be0.5Al.in
  spectra/70kV42deg1Air0.5Be0.5Al.in
  spectra/80kV42deg1Air0.5Be0.5Al.in
  spectra/90kV42deg1Air0.5Be0.5Al.in
  spectra/100kV42deg1Air0.5Be0.5Al.in
  )

foreach(_script ${PEPI2_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
#add_custom_target(PEPI2 DEPENDS PEPI2)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS PEPI2 DESTINATION bin)
