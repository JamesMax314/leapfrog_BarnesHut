# CMAKE generated file: DO NOT EDIT!
# Generated by "NMake Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2019.2.5\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2019.2.5\bin\cmake\win\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles\barnesHut.dir\depend.make

# Include the progress variables for this target.
include CMakeFiles\barnesHut.dir\progress.make

# Include the compile flags for this target's objects.
include CMakeFiles\barnesHut.dir\flags.make

CMakeFiles\barnesHut.dir\main.cpp.obj: CMakeFiles\barnesHut.dir\flags.make
CMakeFiles\barnesHut.dir\main.cpp.obj: ..\main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/barnesHut.dir/main.cpp.obj"
	C:\PROGRA~2\MICROS~1\2019\COMMUN~1\VC\Tools\MSVC\1423~1.281\bin\Hostx86\x86\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoCMakeFiles\barnesHut.dir\main.cpp.obj /FdCMakeFiles\barnesHut.dir\ /FS -c C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\main.cpp
<<

CMakeFiles\barnesHut.dir\main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/barnesHut.dir/main.cpp.i"
	C:\PROGRA~2\MICROS~1\2019\COMMUN~1\VC\Tools\MSVC\1423~1.281\bin\Hostx86\x86\cl.exe > CMakeFiles\barnesHut.dir\main.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\main.cpp
<<

CMakeFiles\barnesHut.dir\main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/barnesHut.dir/main.cpp.s"
	C:\PROGRA~2\MICROS~1\2019\COMMUN~1\VC\Tools\MSVC\1423~1.281\bin\Hostx86\x86\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\barnesHut.dir\main.cpp.s /c C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\main.cpp
<<

# Object files for target barnesHut
barnesHut_OBJECTS = \
"CMakeFiles\barnesHut.dir\main.cpp.obj"

# External object files for target barnesHut
barnesHut_EXTERNAL_OBJECTS =

barnesHut.exe: CMakeFiles\barnesHut.dir\main.cpp.obj
barnesHut.exe: CMakeFiles\barnesHut.dir\build.make
barnesHut.exe: trees.lib
barnesHut.exe: vecMaths.lib
barnesHut.exe: leapfrog.lib
barnesHut.exe: bodies.lib
barnesHut.exe: CMakeFiles\barnesHut.dir\objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable barnesHut.exe"
	"C:\Program Files\JetBrains\CLion 2019.2.5\bin\cmake\win\bin\cmake.exe" -E vs_link_exe --intdir=CMakeFiles\barnesHut.dir --rc=C:\PROGRA~2\WI3CF2~1\10\bin\100183~1.0\x86\rc.exe --mt=C:\PROGRA~2\WI3CF2~1\10\bin\100183~1.0\x86\mt.exe --manifests  -- C:\PROGRA~2\MICROS~1\2019\COMMUN~1\VC\Tools\MSVC\1423~1.281\bin\Hostx86\x86\link.exe /nologo @CMakeFiles\barnesHut.dir\objects1.rsp @<<
 /out:barnesHut.exe /implib:barnesHut.lib /pdb:C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\cmake-build-debug\barnesHut.pdb /version:0.0  /machine:X86 /debug /INCREMENTAL /subsystem:console trees.lib vecMaths.lib leapfrog.lib bodies.lib kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib 
<<

# Rule to build all files generated by this target.
CMakeFiles\barnesHut.dir\build: barnesHut.exe

.PHONY : CMakeFiles\barnesHut.dir\build

CMakeFiles\barnesHut.dir\clean:
	$(CMAKE_COMMAND) -P CMakeFiles\barnesHut.dir\cmake_clean.cmake
.PHONY : CMakeFiles\barnesHut.dir\clean

CMakeFiles\barnesHut.dir\depend:
	$(CMAKE_COMMAND) -E cmake_depends "NMake Makefiles" C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\cmake-build-debug C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\cmake-build-debug C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\cmake-build-debug\CMakeFiles\barnesHut.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles\barnesHut.dir\depend
