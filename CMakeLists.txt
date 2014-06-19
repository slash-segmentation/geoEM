cmake_minimum_required(VERSION 2.8.6)
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
find_package(Boost COMPONENTS chrono system thread)
set(LIBRARIES ${Boost_LIBRARIES})

if (CREATE_GUI)
	# Load CGAL (with Qt4)
	find_package(CGAL REQUIRED COMPONENTS Qt4 MPFR GMP)
	include(${CGAL_USE_FILE})
	set(LIBRARIES ${LIBRARIES} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

	# Load OpenGL
	find_package(OpenGL REQUIRED)
	set(LIBRARIES ${LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})

	# Load Qt4 package including QtCore (default), QtGui (default), QtScript, QtXml, and QtOpenGL.
	set(QT_USE_QTSCRIPT TRUE)
	set(QT_USE_QTXML    TRUE)
	set(QT_USE_QTOPENGL TRUE)
	find_package(Qt4 REQUIRED)
	include(${QT_USE_FILE})
	set(LIBRARIES ${LIBRARIES} ${QT_LIBRARIES})

	# Load QGlViewer package
	find_package(QGLViewer REQUIRED)
	set(LIBRARIES ${LIBRARIES} ${QGLVIEWER_LIBRARIES})
	include_directories(${QGLVIEWER_INCLUDE_DIR})
else (CREATE_GUI)
	# Load CGAL (without Qt4)
	find_package(CGAL REQUIRED MPFR GMP)
	include(${CGAL_USE_FILE})
	set(LIBRARIES ${LIBRARIES} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
endif (CREATE_GUI)
#separate_arguments(LIBRARIES)
#list(REMOVE_DUPLICATES LIBRARIES)


###############################################################################
# Adjust the compiler/linker flags
###############################################################################

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_DEBUG")
set(STATIC_LIBRARY_FLAGS_DEBUG "${STATIC_LIBRARY_FLAGS_DEBUG} -D_DEBUG")
if (CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS}      -std=c++11")
	set(STATIC_LIBRARY_FLAGS           "${STATIC_LIBRARY_FLAGS} -std=c++11")
	set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE}      -O3") #-O -flto -fno-fat-lto-objects
	set(STATIC_LIBRARY_FLAGS_RELEASE   "${STATIC_LIBRARY_FLAGS_RELEASE} -O3") #-O -flto -fno-fat-lto-objects
	# TODO: link-time optimization (need stuff above as well)	
	#set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -flto -O3")
	#set(_CMAKE_TOOLCHAIN_PREFIX gcc-) # we want to use gcc-ar/gcc-nm/gcc-ld which allows us to use static libraries and flto
	#SET (CMAKE_AR      "/usr/bin/gcc-ar")
	#SET (CMAKE_LINKER  "/usr/bin/gcc-ld")
	#SET (CMAKE_NM      "/usr/bin/gcc-nm")
	#SET (CMAKE_OBJDUMP "/usr/bin/gcc-objdump")
	#SET (CMAKE_RANLIB  "/usr/bin/gcc-ranlib")
elseif (MSVC)
	set(CMAKE_CXX_FLAGS                     "${CMAKE_CXX_FLAGS}      /GR-")
	set(STATIC_LIBRARY_FLAGS                "${STATIC_LIBRARY_FLAGS} /GR-")
	#set(CMAKE_CXX_FLAGS_DEBUG              "${CMAKE_CXX_FLAGS_DEBUG}      ")
	#set(STATIC_LIBRARY_FLAGS_DEBUG         "${STATIC_LIBRARY_FLAGS_DEBUG} ")
	set(CMAKE_CXX_FLAGS_RELWITHDEBINFO      "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}      /MP /OPT:REF /OPT:ICF /Oi")
	set(STATIC_LIBRARY_FLAGS_RELWITHDEBINFO "${STATIC_LIBRARY_FLAGS_RELWITHDEBINFO} /MP /OPT:REF /OPT:ICF /Oi")
	set(CMAKE_CXX_FLAGS_MINSIZEREL          "${CMAKE_CXX_FLAGS_MINSIZEREL}          /MP /OPT:REF /OPT:ICF")
	set(STATIC_LIBRARY_FLAGS_MINSIZEREL     "${STATIC_LIBRARY_FLAGS_MINSIZEREL}     /MP /OPT:REF /OPT:ICF")
	set(CMAKE_CXX_FLAGS_RELEASE             "${CMAKE_CXX_FLAGS_RELEASE}             /MP /OPT:REF /OPT:ICF /Oi /GL")
	set(STATIC_LIBRARY_FLAGS_RELEASE        "${STATIC_LIBRARY_FLAGS_RELEASE}        /MP /OPT:REF /OPT:ICF /Oi /GL")
	set(CMAKE_EXE_LINKER_FLAGS_RELEASE      "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")
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
	set(PCH_OUTDIR "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PCH}.dir")
	make_directory(${PCH_OUTDIR})

	string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" FLAGS_VAR_NAME)
	set(PCH_FLAGS ${CMAKE_CXX_FLAGS} ${${FLAGS_VAR_NAME}})
	get_directory_property(DIR_FLAGS INCLUDE_DIRECTORIES)
	foreach (ITEM ${DIR_FLAGS})
		list(APPEND PCH_FLAGS "-I${ITEM}")
	endforeach (ITEM)
	get_directory_property(DIR_FLAGS DEFINITIONS)
	list(APPEND PCH_FLAGS ${DIR_FLAGS})
	separate_arguments(PCH_FLAGS)

	if (MSVC)
		# TODO: debug vs release? (above with ${PROJECT_BINARY_DIR} and ${CMAKE_BUILD_TYPE})
		# TODO: need /Zm136 ?
		set(PCH_OUTPUT "${PCH_OUTDIR}/${PCH}.hpp.pch")

		add_custom_command(
			OUTPUT ${PCH_OUTPUT}
			COMMAND ${CMAKE_CXX_COMPILER} ${PCH_FLAGS} /Yc\"${PCH_SOURCE}\" /Fp\"${PCH_OUTPUT}\" /FI\"${PCH_SOURCE}\" \"${CMAKE_CURRENT_SOURCE_DIR}/${PCH}.cpp\"
			DEPENDS ${PCH_SOURCE}
		)
		add_custom_target(${PCH} DEPENDS ${PCH_OUTPUT})

		set(CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS}      /Yu\"${PCH_SOURCE}\" /Fp\"${PCH_OUTPUT}\" /FI\"${PCH_SOURCE}\"")
		set(STATIC_LIBRARY_FLAGS "${STATIC_LIBRARY_FLAGS} /Yu\"${PCH_SOURCE}\" /Fp\"${PCH_OUTPUT}\" /FI\"${PCH_SOURCE}\"")

	elseif (CMAKE_COMPILER_IS_GNUCXX)
		set(PCH_OUT_H  "${PCH_OUTDIR}/${PCH}.hpp")
		set(PCH_OUTPUT "${PCH_OUTDIR}/${PCH}.hpp.gch")

		add_custom_command(
			OUTPUT ${PCH_OUTPUT}
			COMMAND ${CMAKE_CXX_COMPILER} ${PCH_FLAGS} -x c++-header -o ${PCH_OUTPUT} ${PCH_SOURCE}
			DEPENDS ${PCH_SOURCE}
		)
		add_custom_target(${PCH} DEPENDS ${PCH_OUTPUT})

		set(CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS}      -include \"${PCH_OUT_H}\" -Winvalid-pch")
		set(STATIC_LIBRARY_FLAGS "${STATIC_LIBRARY_FLAGS} -include \"${PCH_OUT_H}\" -Winvalid-pch")
	# elseif Xcode?
	endif (MSVC)
endif (USE_PRECOMPILED_HEADERS)


###############################################################################
# Add static libraries to be compiled
###############################################################################

add_library(GeometryTypes STATIC GeometryTypes.hpp Polyhedron_3.hpp Polygons_3.hpp Skeleton_3.hpp SkeletonGraph_3.hpp GeometryTypes.cpp)
add_library(Utilities STATIC Flags.hpp Flags-MSVC.hpp Heap.hpp Iterators.hpp GeometryUtils.hpp GeometryUtils.cpp Progress.hpp Progress.cpp)
add_library(IO STATIC IO.hpp IO.cpp TriangulatePolyhedron.hpp TriangulatePolyhedron.cpp)
add_library(Intersection STATIC Intersection.hpp Intersection.cpp)
add_library(MedialAxisTransform STATIC MedialAxisTransform.hpp MedialAxisTransform.cpp MedialAxisTransform_Types_MAT.hpp MedialAxisTransform_Types_Triangulation.hpp MedialAxisTransform_Init.cpp MedialAxisTransform_Types_Geodesic.hpp MedialAxisTransform_Geodesic.cpp MedialAxisTransform_IO.hpp MedialAxisTransform_IO.cpp)
add_library(Skeleton STATIC Skeleton.hpp Skeleton.cpp Skeleton_Compute.cpp)
set(MY_LIBRARIES IO Intersection MedialAxisTransform Skeleton GeometryTypes Utilities)
if (CREATE_GUI)
	add_library(Viewer STATIC Viewer.hpp Viewer.cpp GLRender.hpp GLRender.cpp glext.h PolygonSorter.hpp PolygonSorter.cpp)
	list(APPEND MY_LIBRARIES Viewer)
endif (CREATE_GUI)

if (USE_PRECOMPILED_HEADERS)
	add_dependencies(GeometryTypes       ${PCH})
	add_dependencies(Utilities           ${PCH})
	add_dependencies(IO                  ${PCH})
	add_dependencies(Intersection        ${PCH})
	add_dependencies(MedialAxisTransform ${PCH})
	add_dependencies(Skeleton            ${PCH})
	if (CREATE_GUI)
		add_dependencies(Viewer          ${PCH})
	endif (CREATE_GUI)
endif (USE_PRECOMPILED_HEADERS)


###############################################################################
# Add programs to be compiled
###############################################################################

add_executable(testing-app main.cpp)
target_link_libraries(testing-app ${MY_LIBRARIES} ${LIBRARIES})
if (USE_PRECOMPILED_HEADERS)
	add_dependencies(testing-app ${PCH})
endif (USE_PRECOMPILED_HEADERS)

