# MSpace
A module for Cure that used to be part of the toolbox.  This has basic mapping classes.

You should clone this in the Cure directory then add the line to 

Cure/CMakeLists.txt:

add_subdirectory(MSpace)

Then 

cd build
cmake -DCMAKE_INSTALL_PREFIX=../install ..
make -j4
make install

