project(MSpace)

# Add some libraries
add_library(CureBase src/MapObject.cc  src/MapObjectList.cc src/MapBank.cc src/MapFeatureList.cc  src/MapPoint.cc src/MapPointList.cc  src/MapFeature.cc  src/FeatureFun.cc src/PosedFeature.cc src/PosedFeatureList.cc src/MapHelper.cc src/EKF.cc src/EKLocate.cc src/Match.cc )


#add_executable(demo src/demo.cc)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(CureBase PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    PRIVATE )



# Link the libraries

target_link_libraries(CureBase PUBLIC CureFilters CureAddressBank CureSensorData)
# CureTransformation CureMath CureSensorData)

set(CURE_MSPACE_BASE_LIBS CureBase)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/../../../install/include)

# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS ${CURE_MSPACE_BASE_LIBS} EXPORT CureBaseConfig
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into 'cmake'.

install(EXPORT CureBaseConfig DESTINATION share/include/cmake)


# This makes the project importable from the build directory

export(TARGETS ${CURE_MSPACE_BASE_LIBS} CureFilters CureAddressBank CureSensorData FILE CureBaseConfig.cmake )


#target_link_libraries(demo CureFilters CureAddressBank CureSensory CureSensorData CureMath CureTransformation CureGeometry CureUtils -lpthread)

