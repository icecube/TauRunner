# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/austin/PROPOSALsrc

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/austin/PROPOSALsrc/build

# Include any dependencies generated for this target.
include CMakeFiles/PROPOSALtest.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/PROPOSALtest.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/PROPOSALtest.dir/flags.make

CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.o: CMakeFiles/PROPOSALtest.dir/flags.make
CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.o: ../private/test/PROPOSAL.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/austin/PROPOSALsrc/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.o -c /home/austin/PROPOSALsrc/private/test/PROPOSAL.cxx

CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/austin/PROPOSALsrc/private/test/PROPOSAL.cxx > CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.i

CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/austin/PROPOSALsrc/private/test/PROPOSAL.cxx -o CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.s

CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.o.requires:

.PHONY : CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.o.requires

CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.o.provides: CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.o.requires
	$(MAKE) -f CMakeFiles/PROPOSALtest.dir/build.make CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.o.provides.build
.PHONY : CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.o.provides

CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.o.provides.build: CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.o


# Object files for target PROPOSALtest
PROPOSALtest_OBJECTS = \
"CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.o"

# External object files for target PROPOSALtest
PROPOSALtest_EXTERNAL_OBJECTS =

bin/PROPOSALtest: CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.o
bin/PROPOSALtest: CMakeFiles/PROPOSALtest.dir/build.make
bin/PROPOSALtest: lib/libPROPOSAL.a
bin/PROPOSALtest: /usr/lib/x86_64-linux-gnu/libboost_program_options.so
bin/PROPOSALtest: /usr/local/lib/liblog4cplus.so
bin/PROPOSALtest: CMakeFiles/PROPOSALtest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/austin/PROPOSALsrc/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable bin/PROPOSALtest"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/PROPOSALtest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/PROPOSALtest.dir/build: bin/PROPOSALtest

.PHONY : CMakeFiles/PROPOSALtest.dir/build

CMakeFiles/PROPOSALtest.dir/requires: CMakeFiles/PROPOSALtest.dir/private/test/PROPOSAL.cxx.o.requires

.PHONY : CMakeFiles/PROPOSALtest.dir/requires

CMakeFiles/PROPOSALtest.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/PROPOSALtest.dir/cmake_clean.cmake
.PHONY : CMakeFiles/PROPOSALtest.dir/clean

CMakeFiles/PROPOSALtest.dir/depend:
	cd /home/austin/PROPOSALsrc/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/austin/PROPOSALsrc /home/austin/PROPOSALsrc /home/austin/PROPOSALsrc/build /home/austin/PROPOSALsrc/build /home/austin/PROPOSALsrc/build/CMakeFiles/PROPOSALtest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/PROPOSALtest.dir/depend
