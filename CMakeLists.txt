cmake_minimum_required(VERSION 3.16)

project(project_1_fys4411)

# Set CXX standard to C++14
set(CMAKE_CXX_STANDARD 14)

# Check if the user has specified a build type via the command line, e.g., (in the
# build-directory)
#
# $cmake .. -DCMAKE_BUILD_TYPE=Release
#
# if not, set Debug as the default build type.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

# Add all compiler warnings. Remember, a warning in C++ can very often lead to a bug
# later on so _fix your warnings!_
set(CMAKE_CXX_FLAGS "-Wall -Wextra")

# Add debugging flag for tools such as gdb and valgrind.
set(CMAKE_CXX_FLAGS_DEBUG "-g")

# Add full compiler optimization when in the release build type. Other options include
# -O0, -O1, -O2.
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Add directories containing header-files
include_directories(
  .
  Hamiltonians
  InitialStates
  Math
  WaveFunctions
  Solvers
)

# Add sources from all directories
file(
  GLOB LIB_SOURCES "particle.cpp"
  GLOB LIB_SOURCES "sampler.cpp"
  GLOB LIB_SOURCES "system.cpp"
  GLOB LIB_SOURCES "Hamiltonians/*.cpp"
  GLOB LIB_SOURCES "InitialStates/*.cpp"
  GLOB LIB_SOURCES "Math/*.cpp"
  GLOB LIB_SOURCES "Solvers/*.cpp"
  GLOB LIB_SOURCES "WaveFunctions/*.cpp"
)

# Add libary to compile path
add_library(Vmclib ${LIB_SOURCES})

# Create different executables
add_executable(vmc main.cpp)
target_link_libraries(vmc PRIVATE Vmclib)

add_executable(interact maininteracting.cpp)
target_link_libraries(interact PRIVATE Vmclib)

add_executable(timing timing.cpp)
target_link_libraries(timing PRIVATE Vmclib)

add_executable(gradient gradientdescent.cpp)
target_link_libraries(gradient PRIVATE Vmclib)