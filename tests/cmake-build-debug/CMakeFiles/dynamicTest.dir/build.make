# CMAKE generated file: DO NOT EDIT!
# Generated by "NMake Makefiles" Generator, CMake Version 3.14

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE
NULL=nul
!ENDIF
SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = D:\MainDrive\University\Computing\leapfrog_BarnesHut\tests

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = D:\MainDrive\University\Computing\leapfrog_BarnesHut\tests\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles\dynamicTest.dir\depend.make

# Include the progress variables for this target.
include CMakeFiles\dynamicTest.dir\progress.make

# Include the compile flags for this target's objects.
include CMakeFiles\dynamicTest.dir\flags.make

CMakeFiles\dynamicTest.dir\dynamicTest.cpp.obj: CMakeFiles\dynamicTest.dir\flags.make
CMakeFiles\dynamicTest.dir\dynamicTest.cpp.obj: ..\dynamicTest.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\MainDrive\University\Computing\leapfrog_BarnesHut\tests\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/dynamicTest.dir/dynamicTest.cpp.obj"
	C:\PROGRA~2\MIB055~1\2019\COMMUN~1\VC\Tools\MSVC\1422~1.279\bin\Hostx64\x64\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoCMakeFiles\dynamicTest.dir\dynamicTest.cpp.obj /FdCMakeFiles\dynamicTest.dir\ /FS -c D:\MainDrive\University\Computing\leapfrog_BarnesHut\tests\dynamicTest.cpp
<<

CMakeFiles\dynamicTest.dir\dynamicTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dynamicTest.dir/dynamicTest.cpp.i"
	C:\PROGRA~2\MIB055~1\2019\COMMUN~1\VC\Tools\MSVC\1422~1.279\bin\Hostx64\x64\cl.exe > CMakeFiles\dynamicTest.dir\dynamicTest.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\MainDrive\University\Computing\leapfrog_BarnesHut\tests\dynamicTest.cpp
<<

CMakeFiles\dynamicTest.dir\dynamicTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dynamicTest.dir/dynamicTest.cpp.s"
	C:\PROGRA~2\MIB055~1\2019\COMMUN~1\VC\Tools\MSVC\1422~1.279\bin\Hostx64\x64\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\dynamicTest.dir\dynamicTest.cpp.s /c D:\MainDrive\University\Computing\leapfrog_BarnesHut\tests\dynamicTest.cpp
<<

# Object files for target dynamicTest
dynamicTest_OBJECTS = \
"CMakeFiles\dynamicTest.dir\dynamicTest.cpp.obj"

# External object files for target dynamicTest
dynamicTest_EXTERNAL_OBJECTS =

dynamicTest.exe: CMakeFiles\dynamicTest.dir\dynamicTest.cpp.obj
dynamicTest.exe: CMakeFiles\dynamicTest.dir\build.make
dynamicTest.exe: CMakeFiles\dynamicTest.dir\objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=D:\MainDrive\University\Computing\leapfrog_BarnesHut\tests\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable dynamicTest.exe"
	"C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" -E vs_link_exe --intdir=CMakeFiles\dynamicTest.dir --rc=C:\PROGRA~2\WI3CF2~1\10\bin\100183~1.0\x64\rc.exe --mt=C:\PROGRA~2\WI3CF2~1\10\bin\100183~1.0\x64\mt.exe --manifests  -- C:\PROGRA~2\MIB055~1\2019\COMMUN~1\VC\Tools\MSVC\1422~1.279\bin\Hostx64\x64\link.exe /nologo @CMakeFiles\dynamicTest.dir\objects1.rsp @<<
 /out:dynamicTest.exe /implib:dynamicTest.lib /pdb:D:\MainDrive\University\Computing\leapfrog_BarnesHut\tests\cmake-build-debug\dynamicTest.pdb /version:0.0  /machine:x64 /debug /INCREMENTAL /subsystem:console kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib 
<<

# Rule to build all files generated by this target.
CMakeFiles\dynamicTest.dir\build: dynamicTest.exe

.PHONY : CMakeFiles\dynamicTest.dir\build

CMakeFiles\dynamicTest.dir\clean:
	$(CMAKE_COMMAND) -P CMakeFiles\dynamicTest.dir\cmake_clean.cmake
.PHONY : CMakeFiles\dynamicTest.dir\clean

CMakeFiles\dynamicTest.dir\depend:
	$(CMAKE_COMMAND) -E cmake_depends "NMake Makefiles" D:\MainDrive\University\Computing\leapfrog_BarnesHut\tests D:\MainDrive\University\Computing\leapfrog_BarnesHut\tests D:\MainDrive\University\Computing\leapfrog_BarnesHut\tests\cmake-build-debug D:\MainDrive\University\Computing\leapfrog_BarnesHut\tests\cmake-build-debug D:\MainDrive\University\Computing\leapfrog_BarnesHut\tests\cmake-build-debug\CMakeFiles\dynamicTest.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles\dynamicTest.dir\depend

