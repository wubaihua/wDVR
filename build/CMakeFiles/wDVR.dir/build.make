# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/public/wbh/github/wDVR

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/public/wbh/github/wDVR/build

# Include any dependencies generated for this target.
include CMakeFiles/wDVR.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/wDVR.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/wDVR.dir/flags.make

CMakeFiles/wDVR.dir/src/constant.f90.o: CMakeFiles/wDVR.dir/flags.make
CMakeFiles/wDVR.dir/src/constant.f90.o: ../src/constant.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/public/wbh/github/wDVR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/wDVR.dir/src/constant.f90.o"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/public/wbh/github/wDVR/src/constant.f90 -o CMakeFiles/wDVR.dir/src/constant.f90.o

CMakeFiles/wDVR.dir/src/constant.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/wDVR.dir/src/constant.f90.i"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/public/wbh/github/wDVR/src/constant.f90 > CMakeFiles/wDVR.dir/src/constant.f90.i

CMakeFiles/wDVR.dir/src/constant.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/wDVR.dir/src/constant.f90.s"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/public/wbh/github/wDVR/src/constant.f90 -o CMakeFiles/wDVR.dir/src/constant.f90.s

CMakeFiles/wDVR.dir/src/def.f90.o: CMakeFiles/wDVR.dir/flags.make
CMakeFiles/wDVR.dir/src/def.f90.o: ../src/def.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/public/wbh/github/wDVR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object CMakeFiles/wDVR.dir/src/def.f90.o"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/public/wbh/github/wDVR/src/def.f90 -o CMakeFiles/wDVR.dir/src/def.f90.o

CMakeFiles/wDVR.dir/src/def.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/wDVR.dir/src/def.f90.i"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/public/wbh/github/wDVR/src/def.f90 > CMakeFiles/wDVR.dir/src/def.f90.i

CMakeFiles/wDVR.dir/src/def.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/wDVR.dir/src/def.f90.s"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/public/wbh/github/wDVR/src/def.f90 -o CMakeFiles/wDVR.dir/src/def.f90.s

CMakeFiles/wDVR.dir/src/fileio.f90.o: CMakeFiles/wDVR.dir/flags.make
CMakeFiles/wDVR.dir/src/fileio.f90.o: ../src/fileio.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/public/wbh/github/wDVR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building Fortran object CMakeFiles/wDVR.dir/src/fileio.f90.o"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/public/wbh/github/wDVR/src/fileio.f90 -o CMakeFiles/wDVR.dir/src/fileio.f90.o

CMakeFiles/wDVR.dir/src/fileio.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/wDVR.dir/src/fileio.f90.i"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/public/wbh/github/wDVR/src/fileio.f90 > CMakeFiles/wDVR.dir/src/fileio.f90.i

CMakeFiles/wDVR.dir/src/fileio.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/wDVR.dir/src/fileio.f90.s"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/public/wbh/github/wDVR/src/fileio.f90 -o CMakeFiles/wDVR.dir/src/fileio.f90.s

CMakeFiles/wDVR.dir/src/finiteDVR.f90.o: CMakeFiles/wDVR.dir/flags.make
CMakeFiles/wDVR.dir/src/finiteDVR.f90.o: ../src/finiteDVR.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/public/wbh/github/wDVR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building Fortran object CMakeFiles/wDVR.dir/src/finiteDVR.f90.o"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/public/wbh/github/wDVR/src/finiteDVR.f90 -o CMakeFiles/wDVR.dir/src/finiteDVR.f90.o

CMakeFiles/wDVR.dir/src/finiteDVR.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/wDVR.dir/src/finiteDVR.f90.i"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/public/wbh/github/wDVR/src/finiteDVR.f90 > CMakeFiles/wDVR.dir/src/finiteDVR.f90.i

CMakeFiles/wDVR.dir/src/finiteDVR.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/wDVR.dir/src/finiteDVR.f90.s"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/public/wbh/github/wDVR/src/finiteDVR.f90 -o CMakeFiles/wDVR.dir/src/finiteDVR.f90.s

CMakeFiles/wDVR.dir/src/math.f90.o: CMakeFiles/wDVR.dir/flags.make
CMakeFiles/wDVR.dir/src/math.f90.o: ../src/math.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/public/wbh/github/wDVR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building Fortran object CMakeFiles/wDVR.dir/src/math.f90.o"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/public/wbh/github/wDVR/src/math.f90 -o CMakeFiles/wDVR.dir/src/math.f90.o

CMakeFiles/wDVR.dir/src/math.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/wDVR.dir/src/math.f90.i"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/public/wbh/github/wDVR/src/math.f90 > CMakeFiles/wDVR.dir/src/math.f90.i

CMakeFiles/wDVR.dir/src/math.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/wDVR.dir/src/math.f90.s"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/public/wbh/github/wDVR/src/math.f90 -o CMakeFiles/wDVR.dir/src/math.f90.s

CMakeFiles/wDVR.dir/src/potential.f90.o: CMakeFiles/wDVR.dir/flags.make
CMakeFiles/wDVR.dir/src/potential.f90.o: ../src/potential.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/public/wbh/github/wDVR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building Fortran object CMakeFiles/wDVR.dir/src/potential.f90.o"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/public/wbh/github/wDVR/src/potential.f90 -o CMakeFiles/wDVR.dir/src/potential.f90.o

CMakeFiles/wDVR.dir/src/potential.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/wDVR.dir/src/potential.f90.i"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/public/wbh/github/wDVR/src/potential.f90 > CMakeFiles/wDVR.dir/src/potential.f90.i

CMakeFiles/wDVR.dir/src/potential.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/wDVR.dir/src/potential.f90.s"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/public/wbh/github/wDVR/src/potential.f90 -o CMakeFiles/wDVR.dir/src/potential.f90.s

CMakeFiles/wDVR.dir/src/wDVR.f90.o: CMakeFiles/wDVR.dir/flags.make
CMakeFiles/wDVR.dir/src/wDVR.f90.o: ../src/wDVR.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/public/wbh/github/wDVR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building Fortran object CMakeFiles/wDVR.dir/src/wDVR.f90.o"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/public/wbh/github/wDVR/src/wDVR.f90 -o CMakeFiles/wDVR.dir/src/wDVR.f90.o

CMakeFiles/wDVR.dir/src/wDVR.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/wDVR.dir/src/wDVR.f90.i"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/public/wbh/github/wDVR/src/wDVR.f90 > CMakeFiles/wDVR.dir/src/wDVR.f90.i

CMakeFiles/wDVR.dir/src/wDVR.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/wDVR.dir/src/wDVR.f90.s"
	ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/public/wbh/github/wDVR/src/wDVR.f90 -o CMakeFiles/wDVR.dir/src/wDVR.f90.s

# Object files for target wDVR
wDVR_OBJECTS = \
"CMakeFiles/wDVR.dir/src/constant.f90.o" \
"CMakeFiles/wDVR.dir/src/def.f90.o" \
"CMakeFiles/wDVR.dir/src/fileio.f90.o" \
"CMakeFiles/wDVR.dir/src/finiteDVR.f90.o" \
"CMakeFiles/wDVR.dir/src/math.f90.o" \
"CMakeFiles/wDVR.dir/src/potential.f90.o" \
"CMakeFiles/wDVR.dir/src/wDVR.f90.o"

# External object files for target wDVR
wDVR_EXTERNAL_OBJECTS =

wDVR: CMakeFiles/wDVR.dir/src/constant.f90.o
wDVR: CMakeFiles/wDVR.dir/src/def.f90.o
wDVR: CMakeFiles/wDVR.dir/src/fileio.f90.o
wDVR: CMakeFiles/wDVR.dir/src/finiteDVR.f90.o
wDVR: CMakeFiles/wDVR.dir/src/math.f90.o
wDVR: CMakeFiles/wDVR.dir/src/potential.f90.o
wDVR: CMakeFiles/wDVR.dir/src/wDVR.f90.o
wDVR: CMakeFiles/wDVR.dir/build.make
wDVR: CMakeFiles/wDVR.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/public/wbh/github/wDVR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking Fortran executable wDVR"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/wDVR.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/wDVR.dir/build: wDVR

.PHONY : CMakeFiles/wDVR.dir/build

CMakeFiles/wDVR.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/wDVR.dir/cmake_clean.cmake
.PHONY : CMakeFiles/wDVR.dir/clean

CMakeFiles/wDVR.dir/depend:
	cd /home/public/wbh/github/wDVR/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/public/wbh/github/wDVR /home/public/wbh/github/wDVR /home/public/wbh/github/wDVR/build /home/public/wbh/github/wDVR/build /home/public/wbh/github/wDVR/build/CMakeFiles/wDVR.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/wDVR.dir/depend

