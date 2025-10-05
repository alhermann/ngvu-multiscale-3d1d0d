if(NOT dune-angiogenesis_FOUND)
# Whether this module is installed or not
set(dune-angiogenesis_INSTALLED OFF)

# Settings specific to the module

# Package initialization
# Set prefix to source dir
set(PACKAGE_PREFIX_DIR /home/hermann/Schreibtisch/DVC_Paper/Code/Microcirculation/Code_Upload/DUNE/dune-angiogenesis)
macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

#report other information
set_and_check(dune-angiogenesis_PREFIX "${PACKAGE_PREFIX_DIR}")
set_and_check(dune-angiogenesis_INCLUDE_DIRS "/home/hermann/Schreibtisch/DVC_Paper/Code/Microcirculation/Code_Upload/DUNE/dune-angiogenesis")
set(dune-angiogenesis_CXX_FLAGS "-std=c++17 -Wno-deprecated -Wno-inline -Wno-gnu-zero-variadic-macro-arguments -Wno-deprecated-declarations -Wno-c++17-extensions -std=c++14 ")
set(dune-angiogenesis_CXX_FLAGS_DEBUG "-g")
set(dune-angiogenesis_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")
set(dune-angiogenesis_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(dune-angiogenesis_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")
set(dune-angiogenesis_DEPENDS "dune-common;dune-istl;dune-geometry;dune-grid;dune-foamgrid;dune-localfunctions;dune-typetree;dune-functions")
set(dune-angiogenesis_SUGGESTS "")
set(dune-angiogenesis_MODULE_PATH "/home/hermann/Schreibtisch/DVC_Paper/Code/Microcirculation/Code_Upload/DUNE/dune-angiogenesis/cmake/modules")
set(dune-angiogenesis_LIBRARIES "")

# Lines that are set by the CMake build system via the variable DUNE_CUSTOM_PKG_CONFIG_SECTION


#import the target
if(dune-angiogenesis_LIBRARIES)
  get_filename_component(_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)
  include("${_dir}/dune-angiogenesis-targets.cmake")
endif()
endif()
