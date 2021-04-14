cmake_minimum_required(VERSION 2.8.12)

# Change name of the library to something meaningful!
set(PLUGINNAME "myPlugin")
project(${PLUGINNAME} C CXX)

find_package(PythonInterp)
find_package(PythonLibs)
find_package(SWIG REQUIRED)

# get information on radiopropa
execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/findRadioPropa.py" swig_interface OUTPUT_VARIABLE radiopropaSwigPath)
execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/findRadioPropa.py" install_prefix OUTPUT_VARIABLE radiopropaInstallPrefix)
message(STATUS "RadioPropa includes : ${radiopropaInstallPrefix}")

find_path(RADIOPROPA_SWIG_INTERFACE_FILES
  radiopropa.i
  HINTS radiopropa/ share/radiopropa python/ ${radiopropaSwigPath} 
  )

find_path(RADIOPROPA_INCLUDES
  RadioPropa.h
  HINTS radiopropa/ include/radiopropa ${radiopropaInstallPrefix}/include/
  )

message(STATUS "RadioPropa SWIG interface files: ${RADIOPROPA_SWIG_INTERFACE_FILES}")
message(STATUS "RadioPropa includes : ${RADIOPROPA_INCLUDES}")

# swig plugin
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx PROPERTIES GENERATED true )
add_custom_command( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx
  COMMAND swig -builtin -c++ -python ${PY3} -I${RADIOPROPA_SWIG_INTERFACE_FILES} -I${RADIOPROPA_INCLUDES} -I${CMAKE_SOURCE_DIR}/include -I${PYTHON_INCLUDE_DIRS} ${RADIOPROPA_SWIG_DEFINES} -dirprot -o ${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx -outdir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${PLUGINNAME}.i)

include_directories(include ${CMAKE_CURRENT_BINARY_DIR}/include ${RADIOPROPA_INCLUDES} ${PYTHON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})

# build interface and source in one library
add_library(${PLUGINNAME} MODULE 
  ${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}_wrap.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/${PLUGINNAME}.cc)

set_target_properties(${PLUGINNAME} PROPERTIES PREFIX "")
set_target_properties(${PLUGINNAME} PROPERTIES OUTPUT_NAME "_${PLUGINNAME}")

execute_process(
        COMMAND ${PYTHON_EXECUTABLE} -c "import sys; from distutils import sysconfig; sys.stdout.write(
sysconfig.get_python_lib(1,0,prefix='${CMAKE_INSTALL_PREFIX}'))"
        OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )
# Uncomment these lines for installation of the plugin

# install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PLUGINNAME}.py" DESTINATION "${PYTHON_SITE_PACKAGES}/${PLUGINNAME}")
# install(TARGETS ${PLUGINNAME} LIBRARY DESTINATION "${PYTHON_SITE_PACKAGES}/${PLUGINNAME}")


