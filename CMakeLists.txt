cmake_minimum_required(VERSION 3.16)
project(ccc2026 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Binaries land in build/bin. solve6b.py looks for them there (and in the
# repo root for backwards compatibility).
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Search helpers driven by solve6b.py / run.sh.
foreach(tool compress cycle2 cyclek pathk)
  add_executable(${tool} ${tool}.cpp)
  target_compile_options(${tool} PRIVATE -Wall -Wextra)
endforeach()
