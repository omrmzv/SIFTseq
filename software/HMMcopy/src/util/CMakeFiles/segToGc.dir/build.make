# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /workdir/apc88/KTxBS/Bin/HMMcopy

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /workdir/apc88/KTxBS/Bin/HMMcopy

# Include any dependencies generated for this target.
include src/util/CMakeFiles/segToGc.dir/depend.make

# Include the progress variables for this target.
include src/util/CMakeFiles/segToGc.dir/progress.make

# Include the compile flags for this target's objects.
include src/util/CMakeFiles/segToGc.dir/flags.make

src/util/CMakeFiles/segToGc.dir/seg/segToGc.cpp.o: src/util/CMakeFiles/segToGc.dir/flags.make
src/util/CMakeFiles/segToGc.dir/seg/segToGc.cpp.o: src/util/seg/segToGc.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /workdir/apc88/KTxBS/Bin/HMMcopy/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/util/CMakeFiles/segToGc.dir/seg/segToGc.cpp.o"
	cd /workdir/apc88/KTxBS/Bin/HMMcopy/src/util && /bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/segToGc.dir/seg/segToGc.cpp.o -c /workdir/apc88/KTxBS/Bin/HMMcopy/src/util/seg/segToGc.cpp

src/util/CMakeFiles/segToGc.dir/seg/segToGc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/segToGc.dir/seg/segToGc.cpp.i"
	cd /workdir/apc88/KTxBS/Bin/HMMcopy/src/util && /bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /workdir/apc88/KTxBS/Bin/HMMcopy/src/util/seg/segToGc.cpp > CMakeFiles/segToGc.dir/seg/segToGc.cpp.i

src/util/CMakeFiles/segToGc.dir/seg/segToGc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/segToGc.dir/seg/segToGc.cpp.s"
	cd /workdir/apc88/KTxBS/Bin/HMMcopy/src/util && /bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /workdir/apc88/KTxBS/Bin/HMMcopy/src/util/seg/segToGc.cpp -o CMakeFiles/segToGc.dir/seg/segToGc.cpp.s

src/util/CMakeFiles/segToGc.dir/seg/segToGc.cpp.o.requires:
.PHONY : src/util/CMakeFiles/segToGc.dir/seg/segToGc.cpp.o.requires

src/util/CMakeFiles/segToGc.dir/seg/segToGc.cpp.o.provides: src/util/CMakeFiles/segToGc.dir/seg/segToGc.cpp.o.requires
	$(MAKE) -f src/util/CMakeFiles/segToGc.dir/build.make src/util/CMakeFiles/segToGc.dir/seg/segToGc.cpp.o.provides.build
.PHONY : src/util/CMakeFiles/segToGc.dir/seg/segToGc.cpp.o.provides

src/util/CMakeFiles/segToGc.dir/seg/segToGc.cpp.o.provides.build: src/util/CMakeFiles/segToGc.dir/seg/segToGc.cpp.o

# Object files for target segToGc
segToGc_OBJECTS = \
"CMakeFiles/segToGc.dir/seg/segToGc.cpp.o"

# External object files for target segToGc
segToGc_EXTERNAL_OBJECTS =

util/seg/segToGc: src/util/CMakeFiles/segToGc.dir/seg/segToGc.cpp.o
util/seg/segToGc: src/util/CMakeFiles/segToGc.dir/build.make
util/seg/segToGc: lib/libfastahack.a
util/seg/segToGc: lib/libsplit.a
util/seg/segToGc: src/util/CMakeFiles/segToGc.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../util/seg/segToGc"
	cd /workdir/apc88/KTxBS/Bin/HMMcopy/src/util && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/segToGc.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/util/CMakeFiles/segToGc.dir/build: util/seg/segToGc
.PHONY : src/util/CMakeFiles/segToGc.dir/build

src/util/CMakeFiles/segToGc.dir/requires: src/util/CMakeFiles/segToGc.dir/seg/segToGc.cpp.o.requires
.PHONY : src/util/CMakeFiles/segToGc.dir/requires

src/util/CMakeFiles/segToGc.dir/clean:
	cd /workdir/apc88/KTxBS/Bin/HMMcopy/src/util && $(CMAKE_COMMAND) -P CMakeFiles/segToGc.dir/cmake_clean.cmake
.PHONY : src/util/CMakeFiles/segToGc.dir/clean

src/util/CMakeFiles/segToGc.dir/depend:
	cd /workdir/apc88/KTxBS/Bin/HMMcopy && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /workdir/apc88/KTxBS/Bin/HMMcopy /workdir/apc88/KTxBS/Bin/HMMcopy/src/util /workdir/apc88/KTxBS/Bin/HMMcopy /workdir/apc88/KTxBS/Bin/HMMcopy/src/util /workdir/apc88/KTxBS/Bin/HMMcopy/src/util/CMakeFiles/segToGc.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/util/CMakeFiles/segToGc.dir/depend

