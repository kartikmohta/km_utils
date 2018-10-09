cmake_minimum_required(VERSION 3.8)

project(km_utils VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_TESTS "Build the tests." OFF)

find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} SYSTEM INTERFACE
  ${EIGEN3_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
#target_link_libraries(${PROJECT_NAME} INTERFACE Boost::boost)

install(TARGETS ${PROJECT_NAME} EXPORT "${PROJECT_NAME}Targets"
  LIBRARY DESTINATION "lib"
  ARCHIVE DESTINATION "lib"
  RUNTIME DESTINATION "bin"
  INCLUDES DESTINATION "include"
  )

install(DIRECTORY "include/${PROJECT_NAME}"
  DESTINATION "include"
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE
  )

include(cmake/GenerateAndInstallConfig.cmake)
generate_and_install_config_file(
  INCLUDE_DIRS "include"
  LIBRARIES
  DEPS Boost Eigen3
  DEPS_INCLUDE_DIRS Boost EIGEN3
  DEPS_LIBRARIES Boost
  )

#############
## Testing ##
#############

if(BUILD_TESTS)
  add_executable(test_math tests/test_math.cpp)
  target_include_directories(test_math PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
  target_link_libraries(test_math ${PROJECT_NAME})
  enable_testing()
  add_test(Math test_math)
endif()
