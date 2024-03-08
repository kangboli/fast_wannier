cmake_minimum_required(VERSION 3.14)
project(wannieroracles LANGUAGES Fortran)

find_package(MKL)

set(SRCS src/param_parser.f90 src/oracles.f90 src/optimizer.f90 src/zgpadm.f90)
add_library(wannieroracles SHARED ${SRCS})
add_executable(wannier.x src/main.f90 ${SRCS})

if(${MKL_FOUND})
    
  #   add_custom_command(TARGET wannieroracles POST_BUILD
  # COMMAND ${CMAKE_COMMAND} -E echo 
  # "MKL_COMPILE = $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>")
    target_link_libraries(wannier.x PRIVATE $<LINK_ONLY:MKL::MKL>)
    target_compile_options(wannier.x PRIVATE $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    target_link_libraries(wannieroracles PRIVATE $<LINK_ONLY:MKL::MKL>)
else()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    target_link_libraries(wannier.x PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    target_compile_options(wannier.x PRIVATE -Wall -pedantic -march=native -fdefault-integer-8)
    target_link_libraries(wannieroracles PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()


install(TARGETS wannieroracles
	EXPORT wannieroraclesTargets
	LIBRARY DESTINATION lib
	INCLUDES DESTINATION include
	)
