cmake_minimum_required(VERSION 3.1.0)
project(geoEM CXX C)

###############################################################################
# End-user options
###############################################################################
option(CREATE_GUI "Create the GUI programs and options" ON)
option(KERNEL_INEXACT "Use the inexact geometric kernel" ON)
option(USE_PRECOMPILED_HEADERS "Use precompiled headers" ON)

if (${CREATE_GUI})
	add_definitions(-DCREATE_GUI)
endif (${CREATE_GUI})
if (${KERNEL_INEXACT})
	add_definitions(-DKERNEL_INEXACT)
endif (${KERNEL_INEXACT})


###############################################################################
# Load external packages that we need
###############################################################################

# Load Boost
find_package(Boost COMPONENTS timer chrono system thread)
set(LIBRARIES ${Boost_LIBRARIES})
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

if (CREATE_GUI)
	# Load CGAL (with Qt)
	find_package(CGAL REQUIRED COMPONENTS Qt5 Eigen3 MPFR GMP)
	include(${CGAL_USE_FILE})
	set(LIBRARIES ${LIBRARIES} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

	# Load OpenGL
	find_package(OpenGL REQUIRED)
	set(LIBRARIES ${LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})

	# Load Qt package including Core, Gui, Xml, and OpenGL.
	find_package(Qt5 REQUIRED COMPONENTS Core Gui OpenGL Xml)
	include_directories(${QtCore_INCLUDE_DIRS})
	include_directories(${QtGui_INCLUDE_DIRS})
	include_directories(${Qt5OpenGL_INCLUDE_DIRS})
	include_directories(${Qt5Xml_INCLUDE_DIRS})
	add_definitions(${Qt5Core_DEFINITIONS})
	add_definitions(${Qt5Gui_DEFINITIONS})
	add_definitions(${Qt5OpenGL_DEFINITIONS})
	add_definitions(${Qt5Xml_DEFINITIONS})
	set(LIBRARIES ${LIBRARIES} ${Qt5Core_LIBRARIES} ${Qt5Gui_LIBRARIES} ${Qt5OpenGL_LIBRARIES} ${Qt5Xml_LIBRARIES})

	# Load QGlViewer package
	find_package(QGLViewer REQUIRED)
	set(LIBRARIES ${LIBRARIES} ${QGLVIEWER_LIBRARIES})
	include_directories(${QGLVIEWER_INCLUDE_DIR})
else (CREATE_GUI)
	# Load CGAL (without Qt)
	find_package(CGAL REQUIRED COMPONENTS Eigen3 MPFR GMP)
	include(${CGAL_USE_FILE})
	set(LIBRARIES ${LIBRARIES} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
endif (CREATE_GUI)
#separate_arguments(LIBRARIES)
#list(REMOVE_DUPLICATES LIBRARIES)


###############################################################################
# Adjust the compiler/linker flags
###############################################################################

add_definitions(-DCGAL_EIGEN3_ENABLED)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_DEBUG")
set(STATIC_LIBRARY_FLAGS_DEBUG "${STATIC_LIBRARY_FLAGS_DEBUG} -D_DEBUG")
if (CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS}      -std=c++11")
	set(STATIC_LIBRARY_FLAGS           "${STATIC_LIBRARY_FLAGS} -std=c++11")
	set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE}      -O3") #-O -flto -fno-fat-lto-objects
	set(STATIC_LIBRARY_FLAGS_RELEASE   "${STATIC_LIBRARY_FLAGS_RELEASE} -O3") #-O -flto -fno-fat-lto-objects
	# TODO: link-time optimization (need stuff above as well)
	# See: http://hubicka.blogspot.com/2014/04/linktime-optimization-in-gcc-2-firefox.html
	#set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -flto -O3")
	#set(_CMAKE_TOOLCHAIN_PREFIX gcc-) # we want to use gcc-ar/gcc-nm/gcc-ld which allows us to use static libraries and flto
	#SET (CMAKE_AR      "/usr/bin/gcc-ar")
	#SET (CMAKE_LINKER  "/usr/bin/gcc-ld")
	#SET (CMAKE_NM      "/usr/bin/gcc-nm")
	#SET (CMAKE_OBJDUMP "/usr/bin/gcc-objdump")
	#SET (CMAKE_RANLIB  "/usr/bin/gcc-ranlib")
elseif (MSVC)
	string(REPLACE "/GR" "" CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS}")
	string(REPLACE "/GR" "" STATIC_LIBRARY_FLAGS "${STATIC_LIBRARY_FLAGS}")
	#set(CMAKE_CXX_FLAGS_DEBUG              "${CMAKE_CXX_FLAGS_DEBUG}      ")
	#set(STATIC_LIBRARY_FLAGS_DEBUG         "${STATIC_LIBRARY_FLAGS_DEBUG} ")
	set(CMAKE_CXX_FLAGS_RELWITHDEBINFO      "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}      /MP /Oi")
	set(STATIC_LIBRARY_FLAGS_RELWITHDEBINFO "${STATIC_LIBRARY_FLAGS_RELWITHDEBINFO} /MP /Oi")
	set(CMAKE_CXX_FLAGS_MINSIZEREL          "${CMAKE_CXX_FLAGS_MINSIZEREL}          /MP")
	set(STATIC_LIBRARY_FLAGS_MINSIZEREL     "${STATIC_LIBRARY_FLAGS_MINSIZEREL}     /MP")
	set(CMAKE_CXX_FLAGS_RELEASE             "${CMAKE_CXX_FLAGS_RELEASE}             /MP /Oi /GL")
	set(STATIC_LIBRARY_FLAGS_RELEASE        "${STATIC_LIBRARY_FLAGS_RELEASE}        /MP /Oi /GL")
	set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO} /OPT:REF /OPT:ICF")
	set(CMAKE_EXE_LINKER_FLAGS_MINSIZEREL     "${CMAKE_EXE_LINKER_FLAGS_MINSIZEREL}     /OPT:REF /OPT:ICF")
	set(CMAKE_EXE_LINKER_FLAGS_RELEASE        "${CMAKE_EXE_LINKER_FLAGS_RELEASE}        /OPT:REF /OPT:ICF /LTCG")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS}      -std=c++11")
	set(STATIC_LIBRARY_FLAGS           "${STATIC_LIBRARY_FLAGS} -std=c++11")
	# TODO: any other options?
endif ()


###############################################################################
# Setup precompiled header
###############################################################################
if (USE_PRECOMPILED_HEADERS)
	set(PCH PrecompiledHeader)
	set(PCH_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/${PCH}.hpp")

	if (MSVC)
		set(PCH_OUTDIR "${PROJECT_BINARY_DIR}\\${PCH}.dir\\$(Configuration)")
		string(REPLACE "/" "\\" PCH_OUTDIR "${PCH_OUTDIR}")
		set(PCH_OUTPUT "${PCH_OUTDIR}\\${PCH}.hpp.pch")
		set(PCH_OUTOBJ "${PCH_OUTDIR}\\${PCH}.obj")
		set(PCH_OUTPDB "${PCH_OUTDIR}\\vc$(PlatformToolsetVersion).pdb")

		add_library(${PCH} OBJECT ${PCH}.hpp ${PCH}.cpp) # TODO: OBJECT here should NOT produce a .lib but it does
		set_source_files_properties(${PCH}.cpp PROPERTIES COMPILE_FLAGS "/Yc\"${PCH_SOURCE}\" /Fd\"${PCH_OUTPDB}\" /Fo\"${PCH_OUTOBJ}\" /Fp\"${PCH_OUTPUT}\" /FI\"${PCH_SOURCE}\" /Zm140")
		set_target_properties(${PCH} PROPERTIES STATIC_LIBRARY_FLAGS_RELEASE "/LTCG")

		function(add_pch) # TARGETS...
			foreach(TARGET ${ARGN})
				set_target_properties(${TARGET} PROPERTIES COMPILE_FLAGS "/Yu\"${PCH_SOURCE}\" /Fp\"${PCH_OUTPUT}\" /FI\"${PCH_SOURCE}\" /Zm140")
				set_target_properties(${TARGET} PROPERTIES LINK_FLAGS "/ignore:4042")
				add_dependencies(${TARGET} ${PCH})
				target_link_libraries(${TARGET} "${PCH_OUTOBJ}")
				add_custom_command(TARGET ${TARGET} PRE_BUILD
					COMMAND
						IF EXIST \"${PCH_OUTPDB}\" ( \n
							COPY \"${PCH_OUTPDB}\" \"$(IntDir)\" /Y \n
							IF NOT EXIST \"$(IntDir)vc$(PlatformToolsetVersion).pdb\" ( type NUL > \"$(IntDir)vc$(PlatformToolsetVersion).pdb\" ) \n
						)
				)
			endforeach()
		endfunction()

	elseif (CMAKE_COMPILER_IS_GNUCXX)
		set(PCH_OUTDIR "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PCH}.dir")
		make_directory(${PCH_OUTDIR})
		set(PCH_OUT_H  "${PCH_OUTDIR}/${PCH}.hpp")
		set(PCH_OUTPUT "${PCH_OUTDIR}/${PCH}.hpp.gch")

		string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" FLAGS_VAR_NAME)
		set(PCH_FLAGS ${CMAKE_CXX_FLAGS} ${${FLAGS_VAR_NAME}})
		get_directory_property(DIR_FLAGS INCLUDE_DIRECTORIES)
		foreach (ITEM ${DIR_FLAGS})
			list(APPEND PCH_FLAGS "-I${ITEM}")
		endforeach (ITEM)
		get_directory_property(DIR_FLAGS COMPILE_DEFINITIONS)
		foreach (ITEM ${DIR_FLAGS})
			list(APPEND PCH_FLAGS "-D${ITEM}")
		endforeach (ITEM)
		separate_arguments(PCH_FLAGS)

		add_custom_command(
			OUTPUT ${PCH_OUTPUT}
			COMMAND ${CMAKE_CXX_COMPILER} ${PCH_FLAGS} -x c++-header -o ${PCH_OUTPUT} ${PCH_SOURCE}
			DEPENDS ${PCH_SOURCE}
		)
		add_custom_target(${PCH} DEPENDS ${PCH_OUTPUT})
		#target_link_libraries(${PCH} ${LIBRARIES})

		function(add_pch) # TARGETS...
			foreach(TARGET ${ARGN})
				set_target_properties(${TARGET} PROPERTIES COMPILE_FLAGS "-include \"${PCH_OUT_H}\" -Winvalid-pch")
				add_dependencies(${TARGET} ${PCH})
			endforeach()
		endfunction()

	# elseif Xcode?
	else (MSVC)
		message("Precompiled headers not available for your generator.")
		set(USE_PRECOMPILED_HEADERS FALSE)
	endif (MSVC)
endif (USE_PRECOMPILED_HEADERS)


###############################################################################
# Add static libraries to be compiled
###############################################################################

add_library(GeometryTypes STATIC GeometryTypes.hpp Polygons_3.hpp SkeletonGraph_3.hpp GeometryTypes.cpp)
add_library(Utilities STATIC Flags.hpp Heap.hpp Iterators.hpp GeometryUtils.hpp GeometryUtils.cpp Segments2Cylinders.hpp Segments2Cylinders.cpp Points2Spheres.hpp Points2Spheres.cpp Progress.hpp Progress.cpp Strings.hpp)
add_library(IO STATIC IO.hpp IO.cpp IO_core.hpp IO_core.cpp IO_OBJ.hpp IO_OBJ.cpp TriangulatePolyhedron.hpp TriangulatePolyhedron.cpp)
add_library(Intersection STATIC Intersection.hpp Intersection.cpp PolyBuilder.hpp)
add_library(Skeleton STATIC Skeleton.hpp Skeleton.cpp)
set(MY_LIBRARIES IO Intersection Skeleton GeometryTypes Utilities)
if (CREATE_GUI)
	add_library(Viewer STATIC Viewer.hpp Viewer.cpp GLRender.hpp GLRender.cpp glext.h PolygonSorter.hpp PolygonSorter.cpp)
	list(APPEND MY_LIBRARIES Viewer)
	if (CMAKE_COMPILER_IS_GNUCXX)
		set_target_properties(Viewer PROPERTIES COMPILE_FLAGS "-fPIC")
	endif (CMAKE_COMPILER_IS_GNUCXX)
endif (CREATE_GUI)

if (USE_PRECOMPILED_HEADERS)
	add_pch(IO ${MY_LIBRARIES})
endif (USE_PRECOMPILED_HEADERS)

if (MSVC)
	foreach(TARGET ${MY_LIBRARIES})
		set_target_properties(${TARGET} PROPERTIES STATIC_LIBRARY_FLAGS_RELEASE "/LTCG")
	endforeach()
endif (MSVC)


###############################################################################
# Add programs to be compiled
###############################################################################

add_executable(testing-app main.cpp)
target_link_libraries(testing-app ${MY_LIBRARIES} ${LIBRARIES})
#message("X: ${MY_LIBRARIES} ${LIBRARIES}")
if (USE_PRECOMPILED_HEADERS)
	add_pch(testing-app)
endif (USE_PRECOMPILED_HEADERS)


add_executable(mesh-simplifier mesh-simplifier.cpp)
target_link_libraries(mesh-simplifier ${MY_LIBRARIES} ${LIBRARIES})
if (USE_PRECOMPILED_HEADERS)
	add_pch(mesh-simplifier)
endif (USE_PRECOMPILED_HEADERS)
