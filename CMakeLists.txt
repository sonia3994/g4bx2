# $Id: CMakeLists.txt,v 1.3 2015/04/15 16:21:04 marcocci Exp $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(g4bx2)

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
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(G4BX2_SCRIPTS
init_vis.mac
vis.mac
gun.mac
AmBeSource.mac
antineutrino.mac
rdm.mac
rdmchain.mac
Li9He8.mac
please_source_me
neutron.mac
neutrino_C13.mac
ops.mac
snu.mac
multi.mac
scs.mac
PostgreSQLConfig.cmake
.rootrc
)

foreach(_script ${G4BX2_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT REQUIRED)
if(NOT ROOT_FOUND)
  message(STATUS "G4 Examples: ROOT package not found. --> ")
  return()
endif()

#------------------------------------------------------------------------------
# PostgreSQL - required package to connect to bxdb's bx_runvalidation database
#

find_package(PostgreSQL REQUIRED PATHS ${PROJECT_SOURCE_DIR})

if( PostgreSQL_FOUND )
include_directories(${PostgreSQL_INCLUDE_DIRS})
link_directories(${PostgreSQL_LIBRARY_DIRS})
endif( PostgreSQL_FOUND )

 if(NOT PostgreSQL_FOUND)
	message(STATUS "PostgreSQL not found. --> ")
	return()
 endif()
#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
include_directories(${ROOT_INCLUDE_DIR})
#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${PROJECT_SOURCE_DIR}/../shared/include
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)




#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(g4bx2 g4bx2.cc ${sources} ${headers})
target_link_libraries(g4bx2 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES}
	${POSTGRESQL_LIBRARIES})

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(BX2 DEPENDS g4bx2)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS g4bx2 DESTINATION bin)


