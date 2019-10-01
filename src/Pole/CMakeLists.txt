project(MapPole)

# Add some libraries
add_library(CureMapPole src/MapPole.cc src/PosedPole.cc src/PoleHelper.cc)


#add_executable(demo src/demo.cc)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(CureMapPole PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    PRIVATE )

# Link the libraries

target_link_libraries(CureMapPole PUBLIC CureBase CureFilters CureAddressBank CureSensorData)
# CureTransformation CureMath CureSensorData)

set(CURE_MSPACE_MAPPOLE_LIBS CureMapPole)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/../../../install/include)

# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS ${CURE_MSPACE_MAPPOLE_LIBS} EXPORT CureMapPoleConfig
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into 'cmake'.

install(EXPORT CureMapPoleConfig DESTINATION share/include/cmake)


# This makes the project importable from the build directory

export(TARGETS ${CURE_MSPACE_MAPPOLE_LIBS} CureBase CureFilters CureAddressBank CureSensorData FILE CureMapPoleConfig.cmake )


#target_link_libraries(demo CureFilters CureAddressBank CureSensory CureSensorData CureMath CureTransformation CureGeometry CureUtils -lpthread)

