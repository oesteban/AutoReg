cmake_minimum_required(VERSION 2.8)

PROJECT(AutoRegistration)
OPTION( BUILD_SHARED_LIBS OFF )


#INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)

#The following lines are required to use Dart
#ENABLE_TESTING()
#INCLUDE(Dart)

IF (WIN32)
	#SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "C:/win32libs/boost")
	#SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "C:/win32libs/boost/lib")
	SET( Boost_NO_SYSTEM_PATHS ON )
	SET ( BOOST_ROOT "C:/Program Files (x86)/boost/boost_1_47" )
	SET ( BOOST_INCLUDEDIR "C:/Program Files (x86)/boost/boost_1_47" )
	#SET ( Boost_ADDITIONAL_VERSIONS "1.46" "1.46.1")
	#MESSAGE( "BOOST_ROOT=${BOOST_ROOT}")
ENDIF(WIN32)


SET( Boost_USE_STATIC_LIBS ON )
FIND_PACKAGE ( Boost COMPONENTS system filesystem program_options REQUIRED)

# Find ITK.
SET( ITK_USE_STATIC_LIBS ON )
FIND_PACKAGE(ITK REQUIRED)

INCLUDE(${ITK_USE_FILE})

SET(REVIEW_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../Review" "${CMAKE_CURRENT_SOURCE_DIR}/../include" "${CMAKE_CURRENT_SOURCE_DIR}/../txx" "${CMAKE_CURRENT_SOURCE_DIR}/../../GIBUBTransformTest/include" )

INCLUDE_DIRECTORIES( ${REVIEW_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} )
LINK_DIRECTORIES   ( ${Boost_LIBRARY_DIRS} )

FILE(GLOB REVIEW_FILES_H "${REVIEW_INCLUDE_DIR}/*.h")
FILE(GLOB REVIEW_FILES_TXX "${REVIEW_INCLUDE_DIR}/*.txx")

SOURCE_GROUP("Review Files" FILES ${REVIEW_FILES_H} ${REVIEW_FILES_TXX})

IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)



#Set any libraries that your project depends on.
#examples: ITKCommon, VTKRendering, etc
SET(Libraries
ITKCommon
ITKIO
ITKStatistics
ITKNumerics
)
SET(Boost_LIBS
boost_program_options-vc90-mt-1_47
#boost_system-vc90-mt-1_47
#boost_filesystem-vc90-mt-1_47
)

SET(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" )

SET(autoreg_files "${SOURCE_DIR}/autoreg.h" "${SOURCE_DIR}/itkRobbinsMonroOptimizer.cxx" "${SOURCE_DIR}/OptimizerHelper.cxx" ${REVIEW_FILES_H} ${REVIEW_FILES_TXX})
SET(autoreg "autoreg")
ADD_EXECUTABLE(${autoreg} "${SOURCE_DIR}/autoreg.cxx" ${autoreg_files} ) 
TARGET_LINK_LIBRARIES(${autoreg} ${Libraries} ${Boost_PROGRAM_OPTIONS_LIBRARY})


#SET(vph_reg "vph_reg")
#ADD_EXECUTABLE(${vph_reg} "${SOURCE_DIR}/vph_reg.cxx" ${autoreg_files} ) 
#TARGET_LINK_LIBRARIES(${vph_reg} ${Libraries} ${Boost_LIBS} )

#SET(vph_spect_mask "vph_spect_mask")
#ADD_EXECUTABLE(${vph_spect_mask} "${SOURCE_DIR}/maskSpect.cxx" ${autoreg_files} ) 
#TARGET_LINK_LIBRARIES(${vph_spect_mask} ${Libraries} ${Boost_LIBS} )

#ADD_EXECUTABLE(vph_siscom_test "${SOURCE_DIR}/vphSISCOMTest.cxx" )
#TARGET_LINK_LIBRARIES(vph_siscom_test ITKCommon ITKIO ${Boost_LIBS})

#ADD_EXECUTABLE(vph_transform_images "${SOURCE_DIR}/vphTransformImages.cxx"  "fcmax2.cxx")
#TARGET_LINK_LIBRARIES(vph_transform_images ITKCommon ITKIO ${Boost_LIBS})

#ADD_EXECUTABLE(batch_transform "${SOURCE_DIR}/cumulativeTransform.cxx" )
#TARGET_LINK_LIBRARIES(batch_transform ITKCommon ITKIO ${Boost_SYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})

#ADD_EXECUTABLE(transform "${SOURCE_DIR}/simpleTransformImage.cxx" )
#TARGET_LINK_LIBRARIES(transform ITKCommon ITKIO ${Boost_SYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})

#SET(rire_reg "rire_reg")
#ADD_EXECUTABLE(${rire_reg} "${SOURCE_DIR}/rire_reg.cxx" ${autoreg_files} ) 
#TARGET_LINK_LIBRARIES(${rire_reg} ${Libraries} )

#ADD_EXECUTABLE(TransformTest "${SOURCE_DIR}/transformTests.cxx" )
#TARGET_LINK_LIBRARIES(TransformTest ITKCommon ITKIO )

#ADD_EXECUTABLE(transform-points "${SOURCE_DIR}/transform-points.cxx" )
#TARGET_LINK_LIBRARIES(transform-points ITKCommon ITKIO )

#SET(metricMap "metric-map")
#ADD_EXECUTABLE(${metricMap} "${SOURCE_DIR}/metric-map.cxx" ${autoreg_files}) 
#TARGET_LINK_LIBRARIES(${metricMap} ${Libraries} )

#SET(simpleMetricMap "simple-metric-map")
#ADD_EXECUTABLE(${simpleMetricMap} "${SOURCE_DIR}/simple-metric-map.cxx" ${autoreg_files}) 
#TARGET_LINK_LIBRARIES(${simpleMetricMap} ${Libraries} )

# ADD_EXECUTABLE(TransformTest "${SOURCE_DIR}/transformTests.cxx" )
# TARGET_LINK_LIBRARIES(TransformTest ${Libraries} )


SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_CRT_SECURE_NO_DEPRECATE -D_SCL_SECURE_NO_WARNINGS")
SET( CPACK_DEBIAN_PACKAGE_DEPENDS "libc6 (>= 2.13-ubuntu13), libgcc1 (>= 1:4.5.2-8ubuntu4), libboost-system1.42.0 (>= 1.42.0-4ubuntu2), libboost-filesystem1.42.0 (>= 1.42.0-4ubuntu2), libboost-program-options1.42.0 (>= 1.42.0-4ubuntu2)")

#IF(CMAKE_COMPILER_IS_GNUCC)
#  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
#ENDIF(CMAKE_COMPILER_IS_GNUCC)
#IF(CMAKE_COMPILER_IS_GNUCXX)
#  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
#ENDIF(CMAKE_COMPILER_IS_GNUCXX)
#SET(CMAKE_VERBOSE_MAKEFILE ON)
