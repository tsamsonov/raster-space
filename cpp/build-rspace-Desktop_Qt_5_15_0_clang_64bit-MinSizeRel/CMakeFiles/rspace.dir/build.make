# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/Qt/Tools/CMake/CMake.app/Contents/bin/cmake

# The command to remove a file.
RM = /Applications/Qt/Tools/CMake/CMake.app/Contents/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/tsamsonov/GitHub/raster-space/cpp/rspace

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/tsamsonov/GitHub/raster-space/cpp/build-rspace-Desktop_Qt_5_15_0_clang_64bit-MinSizeRel

# Include any dependencies generated for this target.
include CMakeFiles/rspace.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/rspace.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/rspace.dir/flags.make

CMakeFiles/rspace.dir/rspace.cpp.o: CMakeFiles/rspace.dir/flags.make
CMakeFiles/rspace.dir/rspace.cpp.o: /Users/tsamsonov/GitHub/raster-space/cpp/rspace/rspace.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/tsamsonov/GitHub/raster-space/cpp/build-rspace-Desktop_Qt_5_15_0_clang_64bit-MinSizeRel/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/rspace.dir/rspace.cpp.o"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/rspace.dir/rspace.cpp.o -c /Users/tsamsonov/GitHub/raster-space/cpp/rspace/rspace.cpp

CMakeFiles/rspace.dir/rspace.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rspace.dir/rspace.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/tsamsonov/GitHub/raster-space/cpp/rspace/rspace.cpp > CMakeFiles/rspace.dir/rspace.cpp.i

CMakeFiles/rspace.dir/rspace.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rspace.dir/rspace.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/tsamsonov/GitHub/raster-space/cpp/rspace/rspace.cpp -o CMakeFiles/rspace.dir/rspace.cpp.s

# Object files for target rspace
rspace_OBJECTS = \
"CMakeFiles/rspace.dir/rspace.cpp.o"

# External object files for target rspace
rspace_EXTERNAL_OBJECTS =

rspace.cpython-37m-darwin.so: CMakeFiles/rspace.dir/rspace.cpp.o
rspace.cpython-37m-darwin.so: CMakeFiles/rspace.dir/build.make
rspace.cpython-37m-darwin.so: CMakeFiles/rspace.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/tsamsonov/GitHub/raster-space/cpp/build-rspace-Desktop_Qt_5_15_0_clang_64bit-MinSizeRel/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared module rspace.cpython-37m-darwin.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/rspace.dir/link.txt --verbose=$(VERBOSE)
	/usr/bin/strip -x /Users/tsamsonov/GitHub/raster-space/cpp/build-rspace-Desktop_Qt_5_15_0_clang_64bit-MinSizeRel/rspace.cpython-37m-darwin.so

# Rule to build all files generated by this target.
CMakeFiles/rspace.dir/build: rspace.cpython-37m-darwin.so

.PHONY : CMakeFiles/rspace.dir/build

CMakeFiles/rspace.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/rspace.dir/cmake_clean.cmake
.PHONY : CMakeFiles/rspace.dir/clean

CMakeFiles/rspace.dir/depend:
	cd /Users/tsamsonov/GitHub/raster-space/cpp/build-rspace-Desktop_Qt_5_15_0_clang_64bit-MinSizeRel && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/tsamsonov/GitHub/raster-space/cpp/rspace /Users/tsamsonov/GitHub/raster-space/cpp/rspace /Users/tsamsonov/GitHub/raster-space/cpp/build-rspace-Desktop_Qt_5_15_0_clang_64bit-MinSizeRel /Users/tsamsonov/GitHub/raster-space/cpp/build-rspace-Desktop_Qt_5_15_0_clang_64bit-MinSizeRel /Users/tsamsonov/GitHub/raster-space/cpp/build-rspace-Desktop_Qt_5_15_0_clang_64bit-MinSizeRel/CMakeFiles/rspace.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/rspace.dir/depend

