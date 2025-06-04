# Preset that turns on all existing packages. Using the combination
# of this preset followed by the nolib.cmake preset should configure
# a LAMMPS binary, with as many packages included, that can be compiled
# with just a working C++ compiler and an MPI library.

set(ALL_PACKAGES
  BODY
  EXTRA-FIX
  EXTRA-COMPUTE
  MANYBODY
  MOLECULE
  VORONOI)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()
