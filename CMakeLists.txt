# cmake solution referenced from
# stackoverflow.com/questions/17511496/how-to-create-a-shared-library-with-cmake
cmake_minimum_required(VERSION 3.13)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(TPP_VMath VERSION 0.1)

add_library(TPP_VMath STATIC
    src/.cpp/Matrix.cpp
    src/.cpp/Matrix3.cpp
    src/.cpp/VWrap.cpp
    src/.cpp/Vect.cpp
    src/.cpp/Vect1D.cpp
    src/.cpp/Vect2D.cpp
    src/.cpp/Vect3D.cpp
    src/.cpp/Vect4D.cpp
    src/.cpp/VectND.cpp
)

set_target_properties(TPP_VMath PROPERTIES VERSION ${PROJECT_VERSION})

set_target_properties(TPP_VMath PROPERTIES PUBLIC_HEADER src/.h/TPP_VMath.h)

target_include_directories(TPP_VMath PRIVATE .h)
target_include_directories(TPP_VMath PRIVATE .cpp)

include(GNUInstallDirs)

install(TARGETS TPP_VMath
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

configure_file(TPP_VMath.pc.in TPP_VMath.pc @ONLY)

install(FILES ${CMAKE_BINARY_DIR}/TPP_VMath.pc 
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)