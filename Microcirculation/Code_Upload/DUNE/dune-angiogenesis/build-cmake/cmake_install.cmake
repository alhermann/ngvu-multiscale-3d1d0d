# Install script for directory: /home/hermann/Schreibtisch/DVC_Paper/Code/Microcirculation/Code_Upload/DUNE/dune-angiogenesis

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/dunecontrol/dune-angiogenesis" TYPE FILE FILES "/home/hermann/Schreibtisch/DVC_Paper/Code/Microcirculation/Code_Upload/DUNE/dune-angiogenesis/dune.module")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/dune-angiogenesis" TYPE FILE FILES
    "/home/hermann/Schreibtisch/DVC_Paper/Code/Microcirculation/Code_Upload/DUNE/dune-angiogenesis/build-cmake/cmake/pkg/dune-angiogenesis-config.cmake"
    "/home/hermann/Schreibtisch/DVC_Paper/Code/Microcirculation/Code_Upload/DUNE/dune-angiogenesis/build-cmake/dune-angiogenesis-config-version.cmake"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/dune-angiogenesis" TYPE FILE FILES "/home/hermann/Schreibtisch/DVC_Paper/Code/Microcirculation/Code_Upload/DUNE/dune-angiogenesis/config.h.cmake")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/home/hermann/Schreibtisch/DVC_Paper/Code/Microcirculation/Code_Upload/DUNE/dune-angiogenesis/build-cmake/dune-angiogenesis.pc")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/hermann/Schreibtisch/DVC_Paper/Code/Microcirculation/Code_Upload/DUNE/dune-angiogenesis/build-cmake/src/cmake_install.cmake")
  include("/home/hermann/Schreibtisch/DVC_Paper/Code/Microcirculation/Code_Upload/DUNE/dune-angiogenesis/build-cmake/dune/cmake_install.cmake")
  include("/home/hermann/Schreibtisch/DVC_Paper/Code/Microcirculation/Code_Upload/DUNE/dune-angiogenesis/build-cmake/doc/cmake_install.cmake")
  include("/home/hermann/Schreibtisch/DVC_Paper/Code/Microcirculation/Code_Upload/DUNE/dune-angiogenesis/build-cmake/cmake/modules/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/hermann/Schreibtisch/DVC_Paper/Code/Microcirculation/Code_Upload/DUNE/dune-angiogenesis/build-cmake/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
