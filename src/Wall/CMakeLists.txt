project(MapWall)

# Add some libraries
add_library(CureMapWall src/MapWall.cc  src/PosedWall.cc src/WallHelper.cc)


#add_executable(demo src/demo.cc)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(CureMapWall PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    PRIVATE )

# Link the libraries

target_link_libraries(CureMapWall PUBLIC CureBase CureFilters CureAddressBank CureSensorData)
# CureTransformation CureMath CureSensorData)

set(CURE_MSPACE_MAPWALL_LIBS CureMapWall)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/../../../install/include)

# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS ${CURE_MSPACE_MAPWALL_LIBS} EXPORT CureMapWallConfig
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into 'cmake'.

install(EXPORT CureMapWallConfig DESTINATION share/include/cmake)


# This makes the project importable from the build directory

export(TARGETS ${CURE_MSPACE_MAPWALL_LIBS} CureBase CureFilters CureAddressBank CureSensorData FILE CureMapWallConfig.cmake )


#target_link_libraries(demo CureFilters CureAddressBank CureSensory CureSensorData CureMath CureTransformation CureGeometry CureUtils -lpthread)

