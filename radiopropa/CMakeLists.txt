if(APPLE)
	# rpath specific patches
	cmake_minimum_required(VERSION 2.8.12)
else(APPLE)
	# require > 2.8.8 for FILE DOWNLOAD fixes
	# allow < 2.8.12 for debian backports
	cmake_minimum_required(VERSION 2.8.11)
endif(APPLE)

cmake_policy(SET CMP0048 NEW)

project(RadioPropa VERSION 1.0.0 LANGUAGES C CXX)
message(STATUS "RadioPropa version: ${PROJECT_VERSION}")

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(RADIOPROPA_EXTRA_SOURCES)
set(RADIOPROPA_EXTRA_INCLUDES)
set(RADIOPROPA_EXTRA_LIBRARIES)
set(RADIOPROPA_SWIG_DEFINES)

macro(USE_CXX11)
  if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
    endif ()
  else ()
    set (CMAKE_CXX_STANDARD 11)
  endif ()
endmacro(USE_CXX11)
USE_CXX11()

if(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)
	set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--as-needed")
	set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Wl,--as-needed")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--as-needed")
	message(STATUS "Use --as-needed linker flags!")
endif(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)

# Add build type for profiling
SET(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS} -ggdb -fno-omit-frame-pointer")

# Set default build-type to release to enable performance improvements
if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")


# ----------------------------------------------------------------------------
# Version info from Git
# ----------------------------------------------------------------------------
# Embed the information about RadioPropa version into the binaries for easier
# debugging with a reproducibility in mind
option(ENABLE_GIT "Embedding information about RadioPropa version from git" ON)

if(ENABLE_GIT)
	include(GetGitRevisionDescription)
	get_git_head_revision(GIT_REFSPEC GIT_SHA1)
	git_describe(GIT_DESC --tags)
else(ENABLE_GIT)
	set(GIT_REFSPEC "")
	set(GIT_SHA1 "")
	set(GIT_DESC "no-git")
endif(ENABLE_GIT)

message(STATUS "RadioPropa git version: ${GIT_DESC} ${GIT_SHA1} ${GIT_REFSPEC}")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/Version.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/src/Version.cpp" @ONLY)
list(APPEND RADIOPROPA_EXTRA_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/src/Version.cpp")




# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------
# googletest (provided, see https://code.google.com/p/googletest/wiki/FAQ
#     Why is it not recommended use a pre-compiled copy of Google Test?)
option(ENABLE_TESTING "Build tests and enable test target" ON)
if(ENABLE_TESTING)
	include_directories(libs/gtest/include)
	add_subdirectory(libs/gtest)
	if(APPLE)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGTEST_USE_OWN_TR1_TUPLE=1")
	endif(APPLE)
endif(ENABLE_TESTING)

# kiss (provided)
add_subdirectory(libs/kiss)
list(APPEND RADIOPROPA_EXTRA_LIBRARIES kiss)
list(APPEND RADIOPROPA_EXTRA_INCLUDES libs/kiss/include)

# healpix redux (provided)
add_subdirectory(libs/healpix_base)
list(APPEND RADIOPROPA_EXTRA_LIBRARIES healpix_base)
list(APPEND RADIOPROPA_EXTRA_INCLUDES libs/healpix_base/include)
install(DIRECTORY libs/healpix_base/include/ DESTINATION include FILES_MATCHING PATTERN "*.h")

# OpenMP (optional for shared memory multiprocessing)
option(ENABLE_OPENMP "OpenMP for multithreading" ON)
if(ENABLE_OPENMP)
include(FindOpenMP)
	if(OPENMP_FOUND)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	endif(OPENMP_FOUND)
endif(ENABLE_OPENMP)

# Google Performance Tools (optional as possible performance tweak for OpenMP)
find_package(GooglePerfTools)
set(TCMALLOC)
if(GOOGLE_PERFTOOLS_FOUND AND NOT NO_TCMALLOC)
	set(TCMALLOC ${TCMALLOC_LIBRARY})
	list(APPEND RADIOPROPA_EXTRA_LIBRARIES profiler)
endif(GOOGLE_PERFTOOLS_FOUND AND NOT NO_TCMALLOC)

# zlib / zstream-cpp (optional for compressed output files)
find_package(ZLIB)
if(ZLIB_FOUND)
	list(APPEND RADIOPROPA_EXTRA_INCLUDES ${ZLIB_INCLUDE_DIRS})
	list(APPEND RADIOPROPA_EXTRA_INCLUDES "libs/zstream-cpp")
	list(APPEND RADIOPROPA_EXTRA_LIBRARIES ${ZLIB_LIBRARIES})
	add_definitions (-DRADIOPROPA_HAVE_ZLIB)
	list(APPEND RADIOPROPA_SWIG_DEFINES -DRADIOPROPA_HAVE_ZLIB)
	list(APPEND SWIG_INCLUDE_DIRECTORIES ${ZLIB_INCLUDE_DIRS})
endif(ZLIB_FOUND)

# HDF5 (optional for HDF5 output files)
option(ENABLE_HDF5 "HDF5 Support" ON)
if(ENABLE_HDF5)
	find_package( HDF5 COMPONENTS C )
	if(HDF5_FOUND)
		# parallel version requires MPI which is not supported in RadioPropa
		if(NOT HDF5_IS_PARALLEL)
			list(APPEND RADIOPROPA_EXTRA_INCLUDES ${HDF5_INCLUDE_DIRS})
			list(APPEND RADIOPROPA_EXTRA_LIBRARIES ${HDF5_LIBRARIES})
			add_definitions (-DRADIOPROPA_HAVE_HDF5)
			list(APPEND RADIOPROPA_SWIG_DEFINES -DRADIOPROPA_HAVE_HDF5)
			list(APPEND SWIG_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIRS})
			#string(REPLACE " " " -I" HDF5_INCLUDE_DIRS_SWIG ${HDF5_INCLUDE_DIRS})
			#list(APPEND RADIOPROPA_SWIG_DEFINES -I${HDF5_INCLUDE_DIRS_SWIG})
		endif(NOT HDF5_IS_PARALLEL)
	endif(HDF5_FOUND)
endif(ENABLE_HDF5)


# ----------------------------------------------------------------------------
# Fix Apple RPATH
# ----------------------------------------------------------------------------
set(CMAKE_MACOSX_RPATH 1)
option(USE_ABSOLUTE_RPATH "Add absolute rpath to all libraries and executables" ON)
if(APPLE OR USE_ABSOLUTE_RPATH)
	set(CMAKE_SKIP_BUILD_RPATH  FALSE)
	set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
	set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
	set(ABSOLUTE_RPATH "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")
	if(NOT IS_ABSOLUTE ${ABSOLUTE_RPATH})
		set(ABSOLUTE_RPATH ${CMAKE_BINARY_DIR}/${ABSOLUTE_RPATH})
	endif()

	list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${ABSOLUTE_RPATH}" isSystemDir)
	if("${isSystemDir}" STREQUAL "-1")
		message(STATUS "Use absolute RPATH ${ABSOLUTE_RPATH}")
		set(CMAKE_INSTALL_RPATH "${ABSOLUTE_RPATH}")
	endif("${isSystemDir}" STREQUAL "-1")
endif()

# Fix for linker error on mac, Issue: #147
if(APPLE)
	add_definitions(-arch x86_64)
endif(APPLE)

# ----------------------------------------------------------------------------
# Library and Binary
# ----------------------------------------------------------------------------
file(GLOB_RECURSE RADIOPROPA_INCLUDES RELATIVE ${CMAKE_SOURCE_DIR} include/*.h)
include_directories(include ${RADIOPROPA_EXTRA_INCLUDES})

add_library(radiopropa SHARED
	src/Candidate.cpp
	src/Clock.cpp
	src/Common.cpp
	src/Cosmology.cpp
	src/EmissionMap.cpp
	src/GridTools.cpp
	src/Geometry.cpp
	src/Module.cpp
	src/ModuleList.cpp
	src/ParticleState.cpp
	src/ProgressBar.cpp
	src/Random.cpp
	src/Source.cpp
  	src/ScalarField.cpp
  	src/IceModel.cpp
	src/Variant.cpp
	src/module/Boundary.cpp
	src/module/BreakCondition.cpp
	src/module/Discontinuity.cpp
	src/module/RandomN.cpp	
    src/module/RandomScattering.cpp
	# src/module/RandomPointScattering.cpp
	src/module/HorizontalSurface.cpp
	src/module/DensityPerturbation.cpp
	src/module/HDF5Output.cpp
	src/module/Observer.cpp
	src/module/Output.cpp
	src/module/OutputShell.cpp
	src/module/ParticleCollector.cpp
	src/module/PropagationCK.cpp
	src/module/RestrictToRegion.cpp
	src/module/SimplePropagation.cpp
	src/module/TextOutput.cpp
	src/module/Tools.cpp


	${RADIOPROPA_EXTRA_SOURCES}
)
target_link_libraries(radiopropa ${RADIOPROPA_EXTRA_LIBRARIES})

# ----------------------------------------------------------------------------
# Python
# ----------------------------------------------------------------------------
option(ENABLE_PYTHON "Create python library via SWIG" ON)
find_package(PythonInterp)
find_package(PythonLibs)

if(ENABLE_PYTHON AND PYTHONLIBS_FOUND)

	find_package(SWIG 3.0 REQUIRED)
	include(python/Python.cmake)
	include_directories(${PYTHON_INCLUDE_PATH})

	# use SWIG's builtin (maybe should be abandoned if everything will be tested and working with -builtin)
	option(ENABLE_SWIG_BUILTIN "Use SWIG's -builtin option" OFF)

	if(ENABLE_SWIG_BUILTIN)
		set(BUILTIN "-builtin")
		set(PY3 "-py3")
	else(ENABLE_SWIG_BUILTIN)
		set(BUILTIN "")
		set(PY3 "")
	endif(ENABLE_SWIG_BUILTIN)

  if(PYTHON_VERSION_STRING VERSION_GREATER 3.0)
  	list(APPEND RADIOPROPA_SWIG_DEFINES -DSWIG_PYTHON3)
  endif(PYTHON_VERSION_STRING VERSION_GREATER 3.0)

	# tries to import numpy
	execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/checkNumpy.py" OUTPUT_VARIABLE numpyIncludePath)
	if(numpyIncludePath)
		MESSAGE(STATUS "Found numpy headers in " ${numpyIncludePath})
		SET(CMAKE_SWIG_FLAGS -DWITHNUMPY ${CRP})
		list(APPEND RADIOPROPA_SWIG_DEFINES -DWITHNUMPY)
	include_directories(${numpyIncludePath})
	else(numpyIncludePath)
		MESSAGE(STATUS "Numpy not found.")
	endif(numpyIncludePath)

	file(GLOB_RECURSE RADIOPROPA_SWIG_INPUTS python/*.i)
	set_source_files_properties( ${CMAKE_CURRENT_BINARY_DIR}/radiopropa_wrap.cxx PROPERTIES GENERATED true )
	add_custom_command( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/radiopropa_wrap.cxx
		COMMAND swig ${BUILTIN} -c++ -python ${PY3} -I${CMAKE_SOURCE_DIR}/include ${SWIG_INCLUDES} ${RADIOPROPA_SWIG_DEFINES} -dirprot -o ${CMAKE_CURRENT_BINARY_DIR}/radiopropa_wrap.cxx -outdir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR}/python/radiopropa${BUILTIN}.i
		DEPENDS ${RADIOPROPA_SWIG_INPUTS} ${RADIOPROPA_INCLUDES} )

	add_library(radiopropa-swig MODULE ${CMAKE_CURRENT_BINARY_DIR}/radiopropa_wrap.cxx)
	set_target_properties(radiopropa-swig PROPERTIES PREFIX "")
	set_target_properties(radiopropa-swig PROPERTIES COMPILE_FLAGS "-w")
	set_target_properties(radiopropa-swig PROPERTIES OUTPUT_NAME "_radiopropa")
	target_link_libraries(radiopropa-swig radiopropa ${PYTHON_LIBRARIES})

	install(DIRECTORY "${CMAKE_SOURCE_DIR}/python/radiopropa" DESTINATION "${PYTHON_SITE_PACKAGES}")
	install(FILES "${CMAKE_CURRENT_BINARY_DIR}/radiopropa.py" DESTINATION "${PYTHON_SITE_PACKAGES}/radiopropa")
	install(TARGETS radiopropa-swig LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/radiopropa")
  install(FILES ${RADIOPROPA_SWIG_INPUTS} DESTINATION share/radiopropa/swig_interface)
endif(ENABLE_PYTHON AND PYTHONLIBS_FOUND)


# ----------------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------------
add_definitions(-DRADIOPROPA_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
install(TARGETS radiopropa DESTINATION lib)
install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_BINARY_DIR}/include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
#install(DIRECTORY ${CMAKE_BINARY_DIR}/data/ DESTINATION share/radiopropa/ PATTERN ".git" EXCLUDE)


# ----------------------------------------------------------------------------
# Testing
# ----------------------------------------------------------------------------

if(ENABLE_TESTING)
	enable_testing()
	add_executable(testCore test/testCore.cpp)
	target_link_libraries(testCore  radiopropa gtest gtest_main pthread ${COVERAGE_LIBS})
	add_test(testCore testCore)

	add_executable(testVector3 test/testVector3.cpp)
	target_link_libraries(testVector3 radiopropa gtest gtest_main pthread ${COVERAGE_LIBS})
	add_test(testVector3 testVector3)

	add_executable(testModuleList test/testModuleList.cpp)
	target_link_libraries(testModuleList radiopropa gtest gtest_main pthread ${COVERAGE_LIBS})
	add_test(testModuleList testModuleList)

	add_executable(testPropagation test/testPropagation.cpp)
	target_link_libraries(testPropagation radiopropa gtest gtest_main pthread ${COVERAGE_LIBS})
	add_test(testPropagation testPropagation)

	add_executable(testBreakCondition test/testBreakCondition.cpp)
	target_link_libraries(testBreakCondition radiopropa gtest gtest_main pthread ${COVERAGE_LIBS})
	add_test(testBreakCondition testBreakCondition)

	add_executable(testDiscontinuity test/testDiscontinuity.cpp)
	target_link_libraries(testDiscontinuity  radiopropa gtest gtest_main pthread ${COVERAGE_LIBS})
	add_test(testDiscontinuity testDiscontinuity)

	add_executable(testSource test/testSource.cpp)
	target_link_libraries(testSource radiopropa gtest gtest_main pthread ${COVERAGE_LIBS})
	add_test(testSource testSource)

	add_executable(testOutput test/testOutput.cpp)
	target_link_libraries(testOutput radiopropa gtest gtest_main pthread ${COVERAGE_LIBS})
	add_test(testOutput testOutput)

	if(ENABLE_PYTHON AND PYTHONLIBS_FOUND)
		CONFIGURE_FILE(test/testPythonExtension.py.in testPythonExtension.py)
		add_test(testPythonExtension ${PYTHON_EXECUTABLE} testPythonExtension.py)
	endif(ENABLE_PYTHON AND PYTHONLIBS_FOUND)

endif(ENABLE_TESTING)
