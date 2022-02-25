



########## BACKUP
cmake_minimum_required(VERSION 3.10.2)

project(skel_itk)

## Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(skel_itk skel_itk.cxx )

target_link_libraries(skel_itk ${ITK_LIBRARIES})
