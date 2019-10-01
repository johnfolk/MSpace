project(MapHline)

# Add some libraries
add_library(CureMapPoint src/MapHLine.cc  src/PosedHLine.cc  src/HLineHelper.cc)


#add_executable(demo src/demo.cc)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(CureMapHLine PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    PRIVATE )

# Link the libraries

target_link_libraries(CureMapHLine PUBLIC CureBase CureFilters CureAddressBank CureSensorData)
# CureTransformation CureMath CureSensorData)

set(CURE_MSPACE_MAPHLINE_LIBS CureMapHLine)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/../../../install/include)

# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS ${CURE_MSPACE_MAPHLINE_LIBS} EXPORT CureMapHLineConfig
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into 'cmake'.

install(EXPORT CureMapHLineConfig DESTINATION share/include/cmake)


# This makes the project importable from the build directory

export(TARGETS ${CURE_MSPACE_MAPHLINE_LIBS} CureBase CureFilters CureAddressBank CureSensorData FILE CureMapHLineConfig.cmake )


#target_link_libraries(demo CureFilters CureAddressBank CureSensory CureSensorData CureMath CureTransformation CureGeometry CureUtils -lpthread)

