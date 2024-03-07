cmake_minimum_required(VERSION 3.14)
project(wannier LANGUAGES Fortran)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
set(SRCS src/param_parser.f90 src/oracles.f90 src/optimizer.f90 src/zgpadm.f90)
add_library(libwannier.a STATIC ${SRCS})
add_executable(wannier.x src/main.f90 ${SRCS})
target_link_libraries(wannier.x PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_compile_options(wannier.x PRIVATE -Wall -pedantic -march=native)
target_link_libraries(libwannier.a PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

