# This module finds cplex.
#
# The variables CPLEX_ROOT_DIR can serve as a hint stored in the cmake cache.


#--------------------------------------------------
# Get a few paths and infos that help find the right version
if(WIN32)

    # Attempt to find paths using the environment variable CPLEX_STUDIO_DIR<VERSION>
    execute_process(COMMAND cmd /C set CPLEX_STUDIO_DIR OUTPUT_VARIABLE CPLEX_STUDIO_DIR_VAR ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(CPLEX_STUDIO_DIR_VAR)

        # We want to use versions 12.8, 12.9, 12.10, 20.1, 22.1 only for now, where the latter ones are preferred
        string(REGEX MATCH "CPLEX_STUDIO_DIR221[^\n]*" CPLEX_STUDIO_DIR_VAR_CORRECT_VERSION ${CPLEX_STUDIO_DIR_VAR})
        if (NOT(CPLEX_STUDIO_DIR_VAR_CORRECT_VERSION))
            string(REGEX MATCH "CPLEX_STUDIO_DIR201[^\n]*" CPLEX_STUDIO_DIR_VAR_CORRECT_VERSION ${CPLEX_STUDIO_DIR_VAR})
        endif()
        if (NOT(CPLEX_STUDIO_DIR_VAR_CORRECT_VERSION))
            string(REGEX MATCH "CPLEX_STUDIO_DIR1210[^\n]*" CPLEX_STUDIO_DIR_VAR_CORRECT_VERSION ${CPLEX_STUDIO_DIR_VAR})
        endif()
        if (NOT(CPLEX_STUDIO_DIR_VAR_CORRECT_VERSION))
            string(REGEX MATCH "CPLEX_STUDIO_DIR129[^\n]*" CPLEX_STUDIO_DIR_VAR_CORRECT_VERSION ${CPLEX_STUDIO_DIR_VAR})
        endif()
        if (NOT(CPLEX_STUDIO_DIR_VAR_CORRECT_VERSION))
            string(REGEX MATCH "CPLEX_STUDIO_DIR128[^\n]*" CPLEX_STUDIO_DIR_VAR_CORRECT_VERSION ${CPLEX_STUDIO_DIR_VAR})
        endif()
        if (NOT(CPLEX_STUDIO_DIR_VAR_CORRECT_VERSION))
            message("Unable to find CPLEX versions 12.8, 12.9, 12.10, 20.1, or 22.1 based on the following environment variable, ignoring it:\n${CPLEX_STUDIO_DIR_VAR}")
            message("You can manually set the CPLEX_ROOT_DIR variable if you know where your CPLEX version is installed.")
        else()
            string(REGEX REPLACE "CPLEX_STUDIO_DIR" "" CPLEX_STUDIO_DIR_VAR_CORRECT_VERSION ${CPLEX_STUDIO_DIR_VAR_CORRECT_VERSION})
            string(REGEX REPLACE "\"" "" CPLEX_STUDIO_DIR_VAR_CORRECT_VERSION ${CPLEX_STUDIO_DIR_VAR_CORRECT_VERSION})
            string(REGEX MATCH "^[0-9]+" CPLEX_WIN_VERSION ${CPLEX_STUDIO_DIR_VAR_CORRECT_VERSION})
            string(REGEX REPLACE "[0-9]+=" "" CPLEX_STUDIO_DIR_VAR_CORRECT_VERSION ${CPLEX_STUDIO_DIR_VAR_CORRECT_VERSION})
            file(TO_CMAKE_PATH "${CPLEX_STUDIO_DIR_VAR_CORRECT_VERSION}" CPLEX_ROOT_DIR_GUESS)
            message(STATUS "Found CPLEX version ${CPLEX_WIN_VERSION} at '${CPLEX_ROOT_DIR_GUESS}'")
            set(CPLEX_ROOT_DIR "${CPLEX_ROOT_DIR_GUESS}" CACHE PATH "CPLEX root directory.")
            set(CPLEX_VERSION "${CPLEX_WIN_VERSION}" CACHE STRING "CPLEX version (e.g., 128, 129, 1210, 201, or 221)")
        endif()

    endif()

    string(REGEX REPLACE "/VC/bin/.*" "" VISUAL_STUDIO_PATH ${CMAKE_CXX_COMPILER})
    string(REGEX MATCH "Studio/[0-9]+/" CPLEX_WIN_VS_VERSION ${VISUAL_STUDIO_PATH})
    string(REGEX REPLACE "Studio/" "" CPLEX_WIN_VS_VERSION ${CPLEX_WIN_VS_VERSION})
    string(REGEX REPLACE "/" "" CPLEX_WIN_VS_VERSION ${CPLEX_WIN_VS_VERSION})
    set(CPLEX_WIN_VS_VERSION ${CPLEX_WIN_VS_VERSION} CACHE STRING "Visual Studio Version")
    set(CPLEX_WIN_BITNESS "x64" CACHE STRING "x86 or x64 (32bit or 64bit)")
    message(STATUS "CPLEX: searching version for Visual Studio ${CPLEX_WIN_VS_VERSION} ${CPLEX_WIN_BITNESS}")

    set(CPLEX_WIN_LINKAGE mda CACHE STRING "CPLEX linkage variant. One of these: mda (dll, release), mdd (dll, debug), mta (static, release), mtd (static, debug)")
    set(CPLEX_WIN_LINKAGE_D mdd CACHE STRING "CPLEX linkage variant. One of these: mda (dll, release), mdd (dll, debug), mta (static, release), mtd (static, debug)")

    # now, generate platform strings
    set(CPLEX_WIN_PLATFORM "${CPLEX_WIN_BITNESS}_windows_vs${CPLEX_WIN_VS_VERSION}/stat_${CPLEX_WIN_LINKAGE}")
    set(CPLEX_WIN_PLATFORM_D "${CPLEX_WIN_BITNESS}_windows_vs${CPLEX_WIN_VS_VERSION}/stat_${CPLEX_WIN_LINKAGE_D}")
    set(CPLEX_WIN_PLATFORM2 "${CPLEX_WIN_BITNESS}_windows_msvc14/stat_${CPLEX_WIN_LINKAGE}")
    set(CPLEX_WIN_PLATFORM2_D "${CPLEX_WIN_BITNESS}_windows_msvc14/stat_${CPLEX_WIN_LINKAGE_D}")
    if(CPLEX_WIN_BITNESS STREQUAL "x64")
        set(CPLEX_WIN_PLATFORM3 "x64_win64")
    elseif(CPLEX_WIN_BITNESS STREQUAL "x86")
        set(CPLEX_WIN_PLATFORM3 "x86_win32")
    endif()

else()

  set(CPLEX_ROOT_DIR1 "/Applications/CPLEX_Studio_Community221")
  set(CPLEX_ROOT_DIR2 "/usr/local/ILOG/CPLEX_Studio221")
  set(CPLEX_ROOT_DIR3 "~/IBM/ILOG/CPLEX_Studio221")
  set(CPLEX_ROOT_DIR4 "/opt/ibm/ILOG/CPLEX_Studio221")
  set(CPLEX_ROOT_DIR5 "/APPLICATIONS/CPLEX_Studio221")
  set(CPLEX_ROOT_DIR6 "/usr/local/ILOG/CPLEX_Studio201")
  set(CPLEX_ROOT_DIR7 "~/IBM/ILOG/CPLEX_Studio201")
  set(CPLEX_ROOT_DIR8 "/opt/ibm/ILOG/CPLEX_Studio201")
  set(CPLEX_ROOT_DIR9 "/APPLICATIONS/CPLEX_Studio201")
  set(CPLEX_ROOT_DIR10 "/usr/local/ILOG/CPLEX_Studio1210")
  set(CPLEX_ROOT_DIR11 "~/IBM/ILOG/CPLEX_Studio1210")
  set(CPLEX_ROOT_DIR12 "/opt/ibm/ILOG/CPLEX_Studio1210")
  set(CPLEX_ROOT_DIR13 "/APPLICATIONS/CPLEX_Studio1210")
  set(CPLEX_ROOT_DIR14 "/usr/local/ILOG/CPLEX_Studio129")
  set(CPLEX_ROOT_DIR15 "~/IBM/ILOG/CPLEX_Studio129")
  set(CPLEX_ROOT_DIR16 "/opt/ibm/ILOG/CPLEX_Studio129")
  set(CPLEX_ROOT_DIR17 "/APPLICATIONS/CPLEX_Studio129")
  set(CPLEX_ROOT_DIR18 "/usr/local/ILOG/CPLEX_Studio128")
  set(CPLEX_ROOT_DIR19 "~/IBM/ILOG/CPLEX_Studio128")
  set(CPLEX_ROOT_DIR20 "/opt/ibm/ILOG/CPLEX_Studio128")
  set(CPLEX_ROOT_DIR21 "/APPLICATIONS/CPLEX_Studio128")

endif()


#--------------------------------------------------
# Introduce settings in case we haven't found a default yet
if (NOT(CPLEX_ROOT_DIR))
    set(CPLEX_ROOT_DIR "" CACHE PATH "CPLEX root directory.")
endif()

if (NOT(CPLEX_VERSION))
    set(CPLEX_VERSION "" CACHE STRING "CPLEX version (e.g., 128, 129, 1210, 201, or 221)")
endif()


#--------------------------------------------------
# Find include directories & libraries
find_path(CPLEX_INCLUDE_DIR
  ilcplex/cplex.h
  HINTS ${CPLEX_ROOT_DIR}
        ${CPLEX_ROOT_DIR1}
        ${CPLEX_ROOT_DIR2}
        ${CPLEX_ROOT_DIR3}
        ${CPLEX_ROOT_DIR4}
        ${CPLEX_ROOT_DIR5}
        ${CPLEX_ROOT_DIR6}
        ${CPLEX_ROOT_DIR7}
        ${CPLEX_ROOT_DIR8}
        ${CPLEX_ROOT_DIR9}
        ${CPLEX_ROOT_DIR10}
        ${CPLEX_ROOT_DIR11}
        ${CPLEX_ROOT_DIR12}
        ${CPLEX_ROOT_DIR13}
        ${CPLEX_ROOT_DIR14}
        ${CPLEX_ROOT_DIR15}
        ${CPLEX_ROOT_DIR16}
        ${CPLEX_ROOT_DIR17}
        ${CPLEX_ROOT_DIR18}
        ${CPLEX_ROOT_DIR19}
        ${CPLEX_ROOT_DIR20}
        ${CPLEX_ROOT_DIR21}
  NO_DEFAULT_PATH
  PATH_SUFFIXES
    CPLEX_Studio${CPLEX_VERSION}/include
    CPLEX_Studio${CPLEX_VERSION}/cplex/include
    include
	  cplex/include
  )

find_path(CPLEX_CONCERT_INCLUDE_DIR
  ilconcert/iloenv.h
  HINTS ${CPLEX_ROOT_DIR}
        ${CPLEX_ROOT_DIR1}
        ${CPLEX_ROOT_DIR2}
        ${CPLEX_ROOT_DIR3}
        ${CPLEX_ROOT_DIR4}
        ${CPLEX_ROOT_DIR5}
        ${CPLEX_ROOT_DIR6}
        ${CPLEX_ROOT_DIR7}
        ${CPLEX_ROOT_DIR8}
        ${CPLEX_ROOT_DIR9}
        ${CPLEX_ROOT_DIR10}
        ${CPLEX_ROOT_DIR11}
        ${CPLEX_ROOT_DIR12}
        ${CPLEX_ROOT_DIR13}
        ${CPLEX_ROOT_DIR14}
        ${CPLEX_ROOT_DIR15}
        ${CPLEX_ROOT_DIR16}
        ${CPLEX_ROOT_DIR17}
        ${CPLEX_ROOT_DIR18}
        ${CPLEX_ROOT_DIR19}
        ${CPLEX_ROOT_DIR20}
        ${CPLEX_ROOT_DIR21}
  PATH_SUFFIXES
    CPLEX_Studio${CPLEX_VERSION}/include
    CPLEX_Studio${CPLEX_VERSION}/concert/include
    include
	  concert/include
  )

find_library(CPLEX_LIBRARY
  NAMES cplex${CPLEX_VERSION} cplex${CPLEX_VERSION}0 cplex${CPLEX_WIN_VERSION} cplex${CPLEX_WIN_VERSION}0 cplex
  HINTS ${CPLEX_ROOT_DIR}
        ${CPLEX_ROOT_DIR1}
        ${CPLEX_ROOT_DIR2}
        ${CPLEX_ROOT_DIR3}
        ${CPLEX_ROOT_DIR4}
        ${CPLEX_ROOT_DIR5}
        ${CPLEX_ROOT_DIR6}
        ${CPLEX_ROOT_DIR7}
        ${CPLEX_ROOT_DIR8}
        ${CPLEX_ROOT_DIR9}
        ${CPLEX_ROOT_DIR10}
        ${CPLEX_ROOT_DIR11}
        ${CPLEX_ROOT_DIR12}
        ${CPLEX_ROOT_DIR13}
        ${CPLEX_ROOT_DIR14}
        ${CPLEX_ROOT_DIR15}
        ${CPLEX_ROOT_DIR16}
        ${CPLEX_ROOT_DIR17}
        ${CPLEX_ROOT_DIR18}
        ${CPLEX_ROOT_DIR19}
        ${CPLEX_ROOT_DIR20}
        ${CPLEX_ROOT_DIR21}
  PATH_SUFFIXES
     CPLEX_Studio${CPLEX_VERSION}/cplex/lib/x86-64_debian4.0_4.1/static_pic
	   CPLEX_Studio${CPLEX_VERSION}/cplex/lib/x86-64_sles10_4.1/static_pic
	   CPLEX_Studio${CPLEX_VERSION}/cplex/lib/x86-64_linux/static_pic
	   CPLEX_Studio${CPLEX_VERSION}/cplex/lib/x86-64_osx/static_pic
	   CPLEX_Studio${CPLEX_VERSION}/cplex/lib/x86-64_darwin/static_pic
	   CPLEX_Studio${CPLEX_VERSION}/cplex/lib/${CPLEX_WIN_PLATFORM}
	   CPLEX_Studio${CPLEX_VERSION}/cplex/lib/${CPLEX_WIN_PLATFORM2}
     cplex/lib/x86-64_debian4.0_4.1/static_pic
	   cplex/lib/x86-64_sles10_4.1/static_pic
	   cplex/lib/x86-64_linux/static_pic
	   cplex/lib/x86-64_osx/static_pic
	   cplex/lib/x86-64_darwin/static_pic
	   cplex/lib/${CPLEX_WIN_PLATFORM}
	   cplex/lib/${CPLEX_WIN_PLATFORM2}
)


find_library(CPLEX_ILOCPLEX_LIBRARY
  ilocplex
  HINTS ${CPLEX_ROOT_DIR}
        ${CPLEX_ROOT_DIR1}
        ${CPLEX_ROOT_DIR2}
        ${CPLEX_ROOT_DIR3}
        ${CPLEX_ROOT_DIR4}
        ${CPLEX_ROOT_DIR5}
        ${CPLEX_ROOT_DIR6}
        ${CPLEX_ROOT_DIR7}
        ${CPLEX_ROOT_DIR8}
        ${CPLEX_ROOT_DIR9}
        ${CPLEX_ROOT_DIR10}
        ${CPLEX_ROOT_DIR11}
        ${CPLEX_ROOT_DIR12}
        ${CPLEX_ROOT_DIR13}
        ${CPLEX_ROOT_DIR14}
        ${CPLEX_ROOT_DIR15}
        ${CPLEX_ROOT_DIR16}
        ${CPLEX_ROOT_DIR17}
        ${CPLEX_ROOT_DIR18}
        ${CPLEX_ROOT_DIR19}
        ${CPLEX_ROOT_DIR20}
        ${CPLEX_ROOT_DIR21}
  PATH_SUFFIXES
      CPLEX_Studio${CPLEX_VERSION}/cplex/lib/x86-64_debian4.0_4.1/static_pic
	    CPLEX_Studio${CPLEX_VERSION}/cplex/lib/x86-64_sles10_4.1/static_pic
	    CPLEX_Studio${CPLEX_VERSION}/cplex/lib/x86-64_linux/static_pic
	    CPLEX_Studio${CPLEX_VERSION}/cplex/lib/x86-64_osx/static_pic
	    CPLEX_Studio${CPLEX_VERSION}/cplex/lib/x86-64_darwin/static_pic
	    CPLEX_Studio${CPLEX_VERSION}/cplex/lib/${CPLEX_WIN_PLATFORM}
	    CPLEX_Studio${CPLEX_VERSION}/cplex/lib/${CPLEX_WIN_PLATFORM2}
      cplex/lib/x86-64_debian4.0_4.1/static_pic
	    cplex/lib/x86-64_sles10_4.1/static_pic
	    cplex/lib/x86-64_linux/static_pic
	    cplex/lib/x86-64_osx/static_pic
	    cplex/lib/x86-64_darwin/static_pic
	    cplex/lib/${CPLEX_WIN_PLATFORM}
	    cplex/lib/${CPLEX_WIN_PLATFORM2}
  )
message(STATUS "ILOCPLEX Library: ${CPLEX_ILOCPLEX_LIBRARY}")


find_library(CPLEX_CONCERT_LIBRARY
  concert
  HINTS ${CPLEX_ROOT_DIR}
        ${CPLEX_ROOT_DIR1}
        ${CPLEX_ROOT_DIR2}
        ${CPLEX_ROOT_DIR3}
        ${CPLEX_ROOT_DIR4}
        ${CPLEX_ROOT_DIR5}
        ${CPLEX_ROOT_DIR6}
        ${CPLEX_ROOT_DIR7}
        ${CPLEX_ROOT_DIR8}
        ${CPLEX_ROOT_DIR9}
        ${CPLEX_ROOT_DIR10}
        ${CPLEX_ROOT_DIR11}
        ${CPLEX_ROOT_DIR12}
        ${CPLEX_ROOT_DIR13}
        ${CPLEX_ROOT_DIR14}
        ${CPLEX_ROOT_DIR15}
        ${CPLEX_ROOT_DIR16}
        ${CPLEX_ROOT_DIR17}
        ${CPLEX_ROOT_DIR18}
        ${CPLEX_ROOT_DIR19}
        ${CPLEX_ROOT_DIR20}
        ${CPLEX_ROOT_DIR21}
  PATH_SUFFIXES
      CPLEX_Studio${CPLEX_VERSION}/concert/lib/x86-64_debian4.0_4.1/static_pic
  	  CPLEX_Studio${CPLEX_VERSION}/concert/lib/x86-64_sles10_4.1/static_pic
  	  CPLEX_Studio${CPLEX_VERSION}/concert/lib/x86-64_linux/static_pic
  	  CPLEX_Studio${CPLEX_VERSION}/concert/lib/x86-64_osx/static_pic
	    CPLEX_Studio${CPLEX_VERSION}/concert/lib/x86-64_darwin/static_pic
	    CPLEX_Studio${CPLEX_VERSION}/concert/lib/${CPLEX_WIN_PLATFORM}
	    CPLEX_Studio${CPLEX_VERSION}/concert/lib/${CPLEX_WIN_PLATFORM2}
      concert/lib/x86-64_debian4.0_4.1/static_pic
	    concert/lib/x86-64_sles10_4.1/static_pic
	    concert/lib/x86-64_linux/static_pic
	    concert/lib/x86-64_osx/static_pic
	    concert/lib/x86-64_darwin/static_pic
	    concert/lib/${CPLEX_WIN_PLATFORM}
	    concert/lib/${CPLEX_WIN_PLATFORM2}
  )
message(STATUS "CONCERT Library: ${CPLEX_CONCERT_LIBRARY}")

#Debug libraries (available for Windows only)
if(WIN32)
    find_library(CPLEX_LIBRARY_D
      NAMES cplex${CPLEX_VERSION} cplex${CPLEX_VERSION}0 cplex${CPLEX_WIN_VERSION} cplex${CPLEX_WIN_VERSION}0 cplex
      HINTS ${CPLEX_ROOT_DIR}/cplex/lib/${CPLEX_WIN_PLATFORM_D}
            ${CPLEX_ROOT_DIR}/cplex/lib/${CPLEX_WIN_PLATFORM2_D}
            ${CPLEX_ROOT_DIR}/CPLEX_Studio${CPLEX_VERSION}/cplex/lib/${CPLEX_WIN_PLATFORM_D}
            ${CPLEX_ROOT_DIR}/CPLEX_Studio${CPLEX_VERSION}/cplex/lib/${CPLEX_WIN_PLATFORM2_D}
      NO_DEFAULT_PATH
      )
    message(STATUS "CPLEX Library (Debug): ${CPLEX_LIBRARY_D}")

    find_library(CPLEX_ILOCPLEX_LIBRARY_D
      ilocplex
      HINTS ${CPLEX_ROOT_DIR}/cplex/lib/${CPLEX_WIN_PLATFORM_D}
            ${CPLEX_ROOT_DIR}/cplex/lib/${CPLEX_WIN_PLATFORM2_D}
            ${CPLEX_ROOT_DIR}/CPLEX_Studio${CPLEX_VERSION}/cplex/lib/${CPLEX_WIN_PLATFORM_D}
            ${CPLEX_ROOT_DIR}/CPLEX_Studio${CPLEX_VERSION}/cplex/lib/${CPLEX_WIN_PLATFORM2_D}
      NO_DEFAULT_PATH
      )
    message(STATUS "ILOCPLEX Library (Debug): ${CPLEX_ILOCPLEX_LIBRARY_D}")

    find_library(CPLEX_CONCERT_LIBRARY_D
      concert
      HINTS ${CPLEX_ROOT_DIR}/concert/lib/${CPLEX_WIN_PLATFORM_D}
            ${CPLEX_ROOT_DIR}/concert/lib/${CPLEX_WIN_PLATFORM2_D}
            ${CPLEX_ROOT_DIR}/CPLEX_Studio${CPLEX_VERSION}/cplex/lib/${CPLEX_WIN_PLATFORM_D}
            ${CPLEX_ROOT_DIR}/CPLEX_Studio${CPLEX_VERSION}/cplex/lib/${CPLEX_WIN_PLATFORM2_D}
      NO_DEFAULT_PATH
      )
    message(STATUS "CONCERT Library (Debug): ${CPLEX_CONCERT_LIBRARY_D}")
endif()

# Binaries
if(WIN32)
    find_path(CPLEX_BIN_DIR
        cplex${CPLEX_VERSION}.dll cplex${CPLEX_VERSION}0.dll cplex${CPLEX_WIN_VERSION}.dll cplex${CPLEX_WIN_VERSION}0.dll
        HINTS ${CPLEX_ROOT_DIR}/cplex/bin/${CPLEX_WIN_PLATFORM}
              ${CPLEX_ROOT_DIR}/cplex/bin/${CPLEX_WIN_PLATFORM2}
              ${CPLEX_ROOT_DIR}/cplex/bin/${CPLEX_WIN_PLATFORM3}
              ${CPLEX_ROOT_DIR}/CPLEX_Studio${CPLEX_VERSION}/cplex/bin/${CPLEX_WIN_PLATFORM}
              ${CPLEX_ROOT_DIR}/CPLEX_Studio${CPLEX_VERSION}/cplex/bin/${CPLEX_WIN_PLATFORM2}
              ${CPLEX_ROOT_DIR}/CPLEX_Studio${CPLEX_VERSION}/cplex/bin/${CPLEX_WIN_PLATFORM3}
        NO_DEFAULT_PATH
    )
else()
    find_path(CPLEX_BIN_DIR
      cplex
	  HINTS 
      ${CPLEX_ROOT_DIR}
      ${CPLEX_ROOT_DIR1}
			${CPLEX_ROOT_DIR2}
			${CPLEX_ROOT_DIR3}
			${CPLEX_ROOT_DIR4}
			${CPLEX_ROOT_DIR5}
			${CPLEX_ROOT_DIR6}
			${CPLEX_ROOT_DIR7}
			${CPLEX_ROOT_DIR8}
      ${CPLEX_ROOT_DIR9}
      ${CPLEX_ROOT_DIR10}
      ${CPLEX_ROOT_DIR11}
      ${CPLEX_ROOT_DIR12}
      ${CPLEX_ROOT_DIR13}
      ${CPLEX_ROOT_DIR14}
      ${CPLEX_ROOT_DIR15}
      ${CPLEX_ROOT_DIR16}
      ${CPLEX_ROOT_DIR17}
      ${CPLEX_ROOT_DIR18}
      ${CPLEX_ROOT_DIR19}
      ${CPLEX_ROOT_DIR20}
      ${CPLEX_ROOT_DIR21}
	  PATH_SUFFIXES
		  CPLEX_Studio${CPLEX_VERSION}/cplex/bin/x86-64_debian4.0_4.1/static_pic
		  CPLEX_Studio${CPLEX_VERSION}/cplex/bin/x86-64_sles10_4.1/static_pic
		  CPLEX_Studio${CPLEX_VERSION}/cplex/bin/x86-64_linux/static_pic
		  CPLEX_Studio${CPLEX_VERSION}/cplex/bin/x86-64_osx/static_pic
		  CPLEX_Studio${CPLEX_VERSION}/cplex/bin/x86-64_darwin/static_pic
		  cplex/bin/x86-64_debian4.0_4.1/static_pic
		  cplex/bin/x86-64_sles10_4.1/static_pic
		  cplex/bin/x86-64_linux/static_pic
		  cplex/bin/x86-64_osx/static_pic
		  cplex/bin/x86-64_darwin/static_pic
    )
endif()
message(STATUS "CPLEX Bin Dir: ${CPLEX_BIN_DIR}")



#--------------------------------------------------
# Cleanup

include(FindPackageHandleStandardArgs)
if(WIN32)
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(CPLEX DEFAULT_MSG
        CPLEX_LIBRARY CPLEX_LIBRARY_D CPLEX_INCLUDE_DIR CPLEX_ILOCPLEX_LIBRARY CPLEX_ILOCPLEX_LIBRARY_D CPLEX_CONCERT_LIBRARY CPLEX_CONCERT_LIBRARY_D CPLEX_CONCERT_INCLUDE_DIR)
else()
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CPLEX DEFAULT_MSG
        CPLEX_LIBRARY CPLEX_INCLUDE_DIR CPLEX_ILOCPLEX_LIBRARY CPLEX_CONCERT_LIBRARY CPLEX_CONCERT_INCLUDE_DIR)
endif()

if(CPLEX_FOUND)
  set(CPLEX_INCLUDE_DIRS ${CPLEX_INCLUDE_DIR} ${CPLEX_CONCERT_INCLUDE_DIR})
  set(CPLEX_LIBRARIES ${CPLEX_CONCERT_LIBRARY} ${CPLEX_CONCERT_LIBRARY_D} ${CPLEX_ILOCPLEX_LIBRARY} ${CPLEX_ILOCPLEX_LIBRARY_D} ${CPLEX_LIBRARY} ${CPLEX_LIBRARY_D} )
endif(CPLEX_FOUND)

mark_as_advanced(CPLEX_LIBRARY CPLEX_LIBRARY_D CPLEX_INCLUDE_DIR CPLEX_ILOCPLEX_LIBRARY CPLEX_ILOCPLEX_LIBRARY_D CPLEX_CONCERT_INCLUDE_DIR CPLEX_CONCERT_LIBRARY CPLEX_CONCERT_LIBRARY_D)