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
CMAKE_BINARY_DIR = C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\cmake-build-debug-x64

# Include any dependencies generated for this target.
include CMakeFiles\treecode.dir\depend.make

# Include the progress variables for this target.
include CMakeFiles\treecode.dir\progress.make

# Include the compile flags for this target's objects.
include CMakeFiles\treecode.dir\flags.make

CMakeFiles\treecode.dir\pyInterface.cpp.obj: CMakeFiles\treecode.dir\flags.make
CMakeFiles\treecode.dir\pyInterface.cpp.obj: ..\pyInterface.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\cmake-build-debug-x64\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/treecode.dir/pyInterface.cpp.obj"
	C:\PROGRA~2\MICROS~1\2019\COMMUN~1\VC\Tools\MSVC\1423~1.281\bin\Hostx64\x64\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoCMakeFiles\treecode.dir\pyInterface.cpp.obj /FdCMakeFiles\treecode.dir\ /FS -c C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\pyInterface.cpp
<<

CMakeFiles\treecode.dir\pyInterface.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/treecode.dir/pyInterface.cpp.i"
	C:\PROGRA~2\MICROS~1\2019\COMMUN~1\VC\Tools\MSVC\1423~1.281\bin\Hostx64\x64\cl.exe > CMakeFiles\treecode.dir\pyInterface.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\pyInterface.cpp
<<

CMakeFiles\treecode.dir\pyInterface.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/treecode.dir/pyInterface.cpp.s"
	C:\PROGRA~2\MICROS~1\2019\COMMUN~1\VC\Tools\MSVC\1423~1.281\bin\Hostx64\x64\cl.exe @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\treecode.dir\pyInterface.cpp.s /c C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\pyInterface.cpp
<<

# Object files for target treecode
treecode_OBJECTS = \
"CMakeFiles\treecode.dir\pyInterface.cpp.obj"

# External object files for target treecode
treecode_EXTERNAL_OBJECTS =

treecode.cp37-win_amd64.pyd: CMakeFiles\treecode.dir\pyInterface.cpp.obj
treecode.cp37-win_amd64.pyd: CMakeFiles\treecode.dir\build.make
treecode.cp37-win_amd64.pyd: C:\ProgramData\Anaconda3\libs\Python37.lib
treecode.cp37-win_amd64.pyd: trees.lib
treecode.cp37-win_amd64.pyd: vecMaths.lib
treecode.cp37-win_amd64.pyd: leapfrog.lib
treecode.cp37-win_amd64.pyd: bodies.lib
treecode.cp37-win_amd64.pyd: treeShow.lib
treecode.cp37-win_amd64.pyd: CMakeFiles\treecode.dir\objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\cmake-build-debug-x64\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared module treecode.cp37-win_amd64.pyd"
	"C:\Program Files\JetBrains\CLion 2019.2.5\bin\cmake\win\bin\cmake.exe" -E vs_link_dll --intdir=CMakeFiles\treecode.dir --rc=C:\PROGRA~2\WI3CF2~1\10\bin\100183~1.0\x64\rc.exe --mt=C:\PROGRA~2\WI3CF2~1\10\bin\100183~1.0\x64\mt.exe --manifests  -- C:\PROGRA~2\MICROS~1\2019\COMMUN~1\VC\Tools\MSVC\1423~1.281\bin\Hostx64\x64\link.exe /nologo @CMakeFiles\treecode.dir\objects1.rsp @<<
 /out:treecode.cp37-win_amd64.pyd /implib:treecode.lib /pdb:C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\cmake-build-debug-x64\treecode.pdb /dll /version:0.0 /machine:x64 /debug /INCREMENTAL C:\ProgramData\Anaconda3\libs\Python37.lib trees.lib vecMaths.lib leapfrog.lib bodies.lib treeShow.lib kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib  
<<

# Rule to build all files generated by this target.
CMakeFiles\treecode.dir\build: treecode.cp37-win_amd64.pyd

.PHONY : CMakeFiles\treecode.dir\build

CMakeFiles\treecode.dir\clean:
	$(CMAKE_COMMAND) -P CMakeFiles\treecode.dir\cmake_clean.cmake
.PHONY : CMakeFiles\treecode.dir\clean

CMakeFiles\treecode.dir\depend:
	$(CMAKE_COMMAND) -E cmake_depends "NMake Makefiles" C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\cmake-build-debug-x64 C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\cmake-build-debug-x64 C:\Users\10jam\Documents\ComputingL3\Barnes-Hut\leapfrog_BarnesHut\src\cmake-build-debug-x64\CMakeFiles\treecode.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles\treecode.dir\depend
