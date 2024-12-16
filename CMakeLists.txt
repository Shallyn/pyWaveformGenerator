# CMakeLists.txt

cmake_minimum_required(VERSION 3.10)
project(pySEOBNREPHM VERSION 1.0 LANGUAGES C)

# Set C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

# Specify the directories containing header files
include_directories(include)

# Collect all C source files
file(GLOB_RECURSE SOURCES "src/*.c")

# Define the shared library name as EOBLib
add_library(EOBLib SHARED ${SOURCES})

# Remove the 'lib' prefix to name the library EOBLib.so instead of libEOBLib.so
set_target_properties(EOBLib PROPERTIES PREFIX "")
set_target_properties(EOBLib PROPERTIES OUTPUT_NAME "EOBLib")

# Specify the version of the library
set_target_properties(EOBLib PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION 1)

# Find the GSL library
find_package(GSL REQUIRED)
if (GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
    target_link_libraries(EOBLib PRIVATE ${GSL_LIBRARIES})
else()
    message(FATAL_ERROR "GSL library not found. Please install GSL.")
endif()

# **Set CMAKE_INSTALL_PREFIX to the project root directory**
# This ensures that installation paths are relative to the project directory
set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}" CACHE PATH "Install prefix" FORCE)

# Specify where to install the library relative to CMAKE_INSTALL_PREFIX
install(TARGETS EOBLib
        LIBRARY DESTINATION pySEOBNREPHM/lib)

# Optional: Add an executable (e.g., main)
# Uncomment and modify if you have a main executable
# add_executable(pySEOBNREPHM_exec src/main.c)
# target_link_libraries(pySEOBNREPHM_exec PRIVATE EOBLib)
# install(TARGETS pySEOBNREPHM_exec
#         RUNTIME DESTINATION pySEOBNREPHM/bin)
