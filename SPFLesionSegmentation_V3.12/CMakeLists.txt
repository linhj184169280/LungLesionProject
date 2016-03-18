cmake_minimum_required(VERSION 2.8.7 FATAL_ERROR)
project(LesionSegmentation)

find_package(VTK REQUIRED)

find_package(ITK REQUIRED)
include(${VTK_USE_FILE})
include(${ITK_USE_FILE})
add_executable(LesionSegmentation LesionSegmentation.cxx LesionSegmentationCLI.h)
target_link_libraries(LesionSegmentation ${VTK_LIBRARIES} ${ITK_LIBRARIES})