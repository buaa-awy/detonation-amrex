cmake_minimum_required(VERSION 3.24)
project(detonation2d VERSION 0.1.4 LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)

set(
  DETONATION_CUDA_ARCHITECTURES "89"
  CACHE STRING "CUDA compute capability without the decimal point"
)
set(CMAKE_CUDA_ARCHITECTURES "${DETONATION_CUDA_ARCHITECTURES}")

find_package(AMReX CONFIG REQUIRED)

add_executable(detonation2d
  src/main.cpp
  src/Detonation.cpp
  src/FieldExport.cpp
  src/Detonation.H
  src/FieldExport.H
  src/Physics.H
  src/MechanismData.H
)

set_source_files_properties(
  src/main.cpp
  src/Detonation.cpp
  src/FieldExport.cpp
  PROPERTIES LANGUAGE CUDA
)
target_link_libraries(detonation2d PRIVATE AMReX::amrex)

add_executable(chemtest
  src/ChemTest.cpp
  src/Detonation.H
  src/Physics.H
  src/MechanismData.H
)
set_source_files_properties(
  src/ChemTest.cpp
  PROPERTIES LANGUAGE CUDA
)
target_link_libraries(chemtest PRIVATE AMReX::amrex)
