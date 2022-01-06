cmake_minimum_required ( VERSION 3.22 FATAL_ERROR )

set ( PROJECT_NAME "TEFEN" )
set ( PROJECT_VERSION 0.0.1 )
set ( PROJECT_DESCRIPTION "Ultimate teftel engine!")
project ( ${PROJECT_NAME} VERSION ${PROJECT_VERSION} DESCRIPTION ${PROJECT_DESCRIPTION} LANGUAGES CXX )

set ( CMAKE_CXX_STANDARD 20 )
set ( CMAKE_POSITION_INDEPENDENT_CODE ON )
set ( BUILD_TYPE_DIR "release" )



if ( WIN32 ) 
    set ( PLATFORM_DIR "win32" )
    add_definitions( -DWINDOWS=1 )
elseif ( WIN64 )
    set ( PLATFORM_DIR "win64" )
    add_definitions( -DWINDOWS=1 )
endif ( WIN32 )

set ( CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON )


set ( BIN_DIR "${CMAKE_SOURCE_DIR}/bin/${PLATFORM_DIR}/${BUILD_TYPE_DIR}" )

set ( TEFEN_FOLDER "src/lib" )

file(GLOB_RECURSE TEFEN_PUBLIC_HEADERS
    "${TEFEN_FOLDER}/public/*.hxx"
    "${TEFEN_FOLDER}/public/*.h"
    "${TEFEN_FOLDER}/public/*.hpp")

file(GLOB_RECURSE TEFEN_HEADERS
    "${TEFEN_FOLDER}/*.hxx"
    "${TEFEN_FOLDER}/*.h"
    "${TEFEN_FOLDER}/*.hpp")

file(GLOB_RECURSE TEFEN_SOURCES
    "${TEFEN_FOLDER}/*.cxx")

set ( TEST_TEFEN_FOLDER "src/exe" )

file(GLOB_RECURSE TEST_TEFEN_HEADERS
    "${TEST_TEFEN_FOLDER}/*.hxx" )

file(GLOB_RECURSE TEST_TEFEN_SOURCES
    "${TEST_TEFEN_FOLDER}/*.cxx" )



add_library ( GLM INTERFACE)
target_include_directories(GLM INTERFACE "lib/glm/include")

add_library ( MINGW_THREADING INTERFACE)
target_include_directories(MINGW_THREADING INTERFACE "lib/mingw_threading/include")

add_library ( TEFEN SHARED ${TEFEN_SOURCES} )
add_executable ( TEST_TEFEN ${TEST_TEFEN_HEADERS} ${TEST_TEFEN_SOURCES} )

set_target_properties ( TEFEN PROPERTIES
    PUBLIC_HEADER TEFEN_PUBLIC_HEADERS
    RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR}
    VERSION ${PROJECT_VERSION} )

set_target_properties ( TEST_TEFEN PROPERTIES
	OUTPUT_NAME "Test TEFEN"
	RUNTIME_OUTPUT_DIRECTORY "${BIN_DIR}" )

target_include_directories( TEST_TEFEN PRIVATE "${TEFEN_FOLDER}/public" )

set(SDL2_DIR "${CMAKE_SOURCE_DIR}/lib/SDL2/x86_64-w64-mingw32/")
set(SDL2_LIB_DIR ${SDL2_DIR}/lib)

include_directories(${SDL2_DIR}/include)

find_package(OpenCL REQUIRED)

target_link_libraries ( TEST_TEFEN TEFEN OpenCL::OpenCL GLM MINGW_THREADING ${SDL2_LIB_DIR}/libSDL2.dll.a ${SDL2_LIB_DIR}/libSDL2main.a)
target_link_libraries ( TEFEN OpenCL::OpenCL GLM MINGW_THREADING ${SDL2_LIB_DIR}/libSDL2.dll.a ${SDL2_LIB_DIR}/libSDL2main.a)

add_definitions(-DSDL_MAIN_HANDLED)