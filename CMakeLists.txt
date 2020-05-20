# ZDEngine/Core/CorePointsRegistration/CMakeLists.txt

project(CorePointsRegistration)

#-----------------------------------------------------------------------------
# Add cmake utils and documentation utils
include( ${CMAKE_CURRENT_LIST_DIR}/config/cmake/ZD_utils.cmake ) 
include(${CMAKE_CURRENT_LIST_DIR}/config/cmake/doxygen/doxygen.cmake) 


doxygen_add_library(Core/CorePointsRegistration
  PACKAGE Core
  DESCRIPTION "ZDEngine's Core Points Registration Library"
  )

  include_directories( ${EXTERNAL_LIBS_INCLUDE_DIR} )
  file( GLOB_RECURSE header_list_CorePointsRegistration 	"${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/*.inl" )
  file( GLOB_RECURSE source_list_CorePointsRegistration   "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx" "${CMAKE_CURRENT_SOURCE_DIR}/*.c" )
		
  set( CorePointsRegistration_sources
		${header_list_CorePointsRegistration}
		${source_list_CorePointsRegistration}
	 )
		
source_group_custom( ${CMAKE_CURRENT_SOURCE_DIR} "Header Files" ${header_list_CorePointsRegistration} )
source_group_custom( ${CMAKE_CURRENT_SOURCE_DIR} "Source Files" ${source_list_CorePointsRegistration} )

SOURCE_GROUP_BY_FOLDER (CorePointsRegistration)
	 
zd_add_library(LIBRARY_NAME ${ZD_LIB_PREFIX}CorePointsRegistration
  LIBRARY_SOURCES ${CorePointsRegistration_sources}
  HEADER_INSTALL_DIR CorePointsRegistration)
  
target_compile_options(${ZD_LIB_PREFIX}CorePointsRegistration PUBLIC ${active_compile_opts})
target_link_libraries( ${ZD_LIB_PREFIX}CorePointsRegistration ${ZD_LIB_PREFIX}CoreMaths ${ZD_LIB_PREFIX}CoreCommon  ${CorePointsRegistration_needed_libraries})
set_target_properties (CorePointsRegistration PROPERTIES FOLDER ZDEngine)
	