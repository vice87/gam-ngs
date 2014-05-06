# ==========================
# BamTools CMakeLists.txt
# (c) 2010 Derek Barnett
#
# top-level
# ==========================

project( BamTools )

# BamTools version information
set (BamTools_VERSION_MAJOR 2)
set (BamTools_VERSION_MINOR 3)
set (BamTools_VERSION_BUILD 0)

# add our includes root path
include_directories( bamtools-2.3.0/src )

# list subdirectories to build in
add_subdirectory( bamtools-2.3.0/src/api )
