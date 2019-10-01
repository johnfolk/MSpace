project(RelativePointFeature)

# Add some libraries
add_library(CureRelativePointFeature src/RelativePointFeature.cc src/FeatureDescriptors.cc src/PosedRelativePointFeature.cc src/RelativePointFeatureHelper.cc src/MapPoseTree.cc)


#add_executable(demo src/demo.cc)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(CureRelativePointFeature PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    PRIVATE )

# Link the libraries

target_link_libraries(CureRelativePointFeature PUBLIC CureBase CureFilters CureAddressBank CureSensorData)
# CureTransformation CureMath CureSensorData)

set(CURE_MSPACE_RELATIVEPOINTFEATURE_LIBS CureRelativePointFeature)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/../../../install/include)

# 'make install' to the correct locations (provided by GNUInstallDirs).
install(TARGETS ${CURE_MSPACE_RELATIVEPOINTFEATURE_LIBS} EXPORT CureRelativePointFeatureConfig
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})  # This is for Windows
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# This makes the project importable from the install directory
# Put config file in per-project dir (name MUST match), can also
# just go into 'cmake'.

install(EXPORT CureRelativePointFeatureConfig DESTINATION share/include/cmake)


# This makes the project importable from the build directory

export(TARGETS ${CURE_MSPACE_RELATIVEPOINTFEATURE_LIBS} CureBase CureFilters CureAddressBank CureSensorData FILE CureRelativePointFeatureConfig.cmake )


#target_link_libraries(demo CureFilters CureAddressBank CureSensory CureSensorData CureMath CureTransformation CureGeometry CureUtils -lpthread)

