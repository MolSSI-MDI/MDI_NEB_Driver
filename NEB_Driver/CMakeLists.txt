# Locate MPI

find_package(MPI)
if(MPI_FOUND)
   include_directories(${MPI_INCLUDE_PATH})
else()
   configure_file(${CMAKE_CURRENT_SOURCE_DIR}/STUBS_MPI/mpi.h ${CMAKE_CURRENT_BINARY_DIR}/STUBS_MPI/mpi.h COPYONLY)
endif()



# Link to MDI

set( MDI_LOCATION ${CMAKE_BINARY_DIR}/lib/mdi/MDI_Library/ )
link_directories( ${MDI_LOCATION} )
include_directories(${MDI_LOCATION})



# Compile the driver

add_executable(NEB_Driver
               NEB_Driver.cpp
	       NEB_utils.cpp
	       NEB_utils.h)
target_link_libraries(NEB_Driver mdi
                      ${MPI_LIBRARIES})



# Ensure that MPI is properly linked

if(NOT MPI_FOUND)
   target_include_directories(NEB_Driver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/STUBS_MPI/)
endif()
if(MPI_COMPILE_FLAGS)
   set_target_properties(NEB_Driver PROPERTIES
      COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()
if(MPI_LINK_FLAGS)
   set_target_properties(NEB_Driver PROPERTIES
      LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
