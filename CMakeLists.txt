cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

message("Checking files for CPLEX.")
project(cplex CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-fdiagnostics-color=always -W -Wshadow -Wall -Wno-unused-result")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lemon")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_DEBUG")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GENERATE")
#set(CMAKE_CXX_CLANG_TIDY clang-tidy -checks=-*,readability-*)
#set(DCMAKE_EXPORT_COMPILE_COMMANDS 1)
#set(CMAKE_CXX_CLANG_TIDY clang-tidy -checks=cppcoreguidelines-*,readability-*, -header-filter=.*)
#set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  MESSAGE("Clang")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  MESSAGE("GNU")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  MESSAGE("Intel")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  MESSAGE("MSVC")
endif()


if(MSVC)
	#----------- Since we are using CPLEX >=12.8, only 64 bit is allowed -------------------
	if(NOT "${CMAKE_SIZEOF_VOID_P}" STREQUAL "8" )
		message(FATAL_ERROR "CPLEX >=12.8 only supports 64 bit. Please choose the 'x64' version of the generator (or, in newer CMake versions, 'x64' as optional platform for generator).")
	endif()
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}")

find_package(CPLEX)


if(${CPLEX_FOUND})

	# On Windows, CPLEX comes as dlls - on other systems as static libraries
	if(WIN32)
		add_library(cplex-lib SHARED IMPORTED GLOBAL)
		set_target_properties(cplex-lib  PROPERTIES  IMPORTED_IMPLIB_RELEASE ${CPLEX_LIBRARY})
		set_target_properties(cplex-lib  PROPERTIES  IMPORTED_IMPLIB_GENERAL ${CPLEX_LIBRARY})
		set_target_properties(cplex-lib  PROPERTIES  IMPORTED_IMPLIB_DEBUG ${CPLEX_LIBRARY_D})

		add_library(cplex-concert SHARED IMPORTED GLOBAL)
		set_target_properties(cplex-concert  PROPERTIES  INTERFACE_INCLUDE_DIRECTORIES ${CPLEX_CONCERT_INCLUDE_DIR})
		set_target_properties(cplex-concert  PROPERTIES  IMPORTED_IMPLIB_RELEASE   ${CPLEX_CONCERT_LIBRARY})
		set_target_properties(cplex-concert  PROPERTIES  IMPORTED_IMPLIB_GENERAL   ${CPLEX_CONCERT_LIBRARY})
		set_target_properties(cplex-concert  PROPERTIES  IMPORTED_IMPLIB_DEBUG   ${CPLEX_CONCERT_LIBRARY_D})

		add_library(cplex-ilo SHARED IMPORTED GLOBAL)
		set_target_properties(cplex-ilo  PROPERTIES  IMPORTED_IMPLIB_RELEASE   ${CPLEX_ILOCPLEX_LIBRARY})
		set_target_properties(cplex-ilo  PROPERTIES  IMPORTED_IMPLIB_GENERAL   ${CPLEX_ILOCPLEX_LIBRARY})
		set_target_properties(cplex-ilo  PROPERTIES  IMPORTED_IMPLIB_DEBUG   ${CPLEX_ILOCPLEX_LIBRARY_D})
	else()
		add_library(cplex-lib STATIC IMPORTED GLOBAL)
		set_target_properties(cplex-lib  PROPERTIES  IMPORTED_LOCATION ${CPLEX_LIBRARY})
		target_link_libraries(cplex-lib INTERFACE m pthread dl)

		add_library(cplex-concert STATIC IMPORTED GLOBAL)
		set_target_properties(cplex-concert  PROPERTIES  IMPORTED_LOCATION   ${CPLEX_CONCERT_LIBRARY})

		add_library(cplex-ilo STATIC IMPORTED GLOBAL)
		set_target_properties(cplex-ilo  PROPERTIES  IMPORTED_LOCATION   ${CPLEX_ILOCPLEX_LIBRARY})
	endif()

	set_target_properties(cplex-lib  PROPERTIES  INTERFACE_INCLUDE_DIRECTORIES ${CPLEX_INCLUDE_DIR})
	set_target_properties(cplex-lib  PROPERTIES  INTERFACE_COMPILE_DEFINITIONS IL_STD)
	set_target_properties(cplex-concert  PROPERTIES  INTERFACE_INCLUDE_DIRECTORIES ${CPLEX_CONCERT_INCLUDE_DIR})
	set_target_properties(cplex-concert  PROPERTIES  INTERFACE_COMPILE_DEFINITIONS IL_STD)

	add_library(cplex-available INTERFACE)
	set_target_properties(cplex-available  PROPERTIES  INTERFACE_COMPILE_DEFINITIONS HAVE_CPLEX)

	add_library(cplex INTERFACE)
	target_link_libraries(cplex INTERFACE cplex-available cplex-ilo cplex-concert cplex-lib)
	target_compile_options(cplex
        INTERFACE
            $<$<CXX_COMPILER_ID:AppleClang>: -Wno-deprecated-declarations>
            $<$<CXX_COMPILER_ID:GNU>:-Wno-deprecated-declarations>
        )

else()

	add_library(cplex-lib INTERFACE)
	add_library(cplex-concert INTERFACE)
	add_library(cplex-ilo INTERFACE)
	add_library(cplex-available INTERFACE)
	add_library(cplex INTERFACE)
	message("CPLEX could not be found. A dummy target will be used instead and the flag HAVE_CPLEX will not be defined.")
	message("This may be OK depending on the application (e.g., in MAiNGO), but you will not be able to actually use CPLEX.")

endif()


if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
	if(${CPLEX_FOUND})
		include_directories(includes)
		file(GLOB SOURCES src/*.cpp)
		add_executable(robust_path_tariff ${SOURCES})
	target_link_libraries(robust_path_tariff cplex)
	else()
		message(FATAL_ERROR "CPLEX could not be found. Cannot build CPLEX test.")
	endif()
	if (MSVC)
		set_target_properties(cplex-test PROPERTIES LINK_FLAGS /ignore:4099)
	endif()
endif()
