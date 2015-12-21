cmake_minimum_required(VERSION 2.8.7)

project(libwebp)

# Options for coder / decoder executables.
option(WEBP_BUILD_CWEBP "Build the cwebp command line tool." OFF)
option(WEBP_BUILD_DWEBP "Build the dwebp command line tool." OFF)
option(WEBP_EXPERIMENTAL_FEATURES "Build with experimental features." OFF)
option(WEBP_FORCE_ALIGNED "Force aligned memory operations." OFF)

set(WEBP_DEP_LIBRARIES)
set(WEBP_DEP_INCLUDE_DIRS)

if(WEBP_BUILD_CWEBP OR WEBP_BUILD_DWEBP)
endif()

################################################################################
# Generate the config.h to compile with secific intrinsics / libs.

## Check for compiler options.
include(CheckCSourceCompiles)
check_c_source_compiles("
    in main(void) {
      __builtin_bswap16(0);
      return 0;
    }
  "
  HAVE_BUILTIN_BSWAP16
)
check_c_source_compiles("
    in main(void) {
      __builtin_bswap32(0);
      return 0;
    }
  "
  HAVE_BUILTIN_BSWAP32
)
check_c_source_compiles("
    in main(void) {
      __builtin_bswap64(0);
      return 0;
    }
  "
  HAVE_BUILTIN_BSWAP64
)

## Check for libraries.
find_package(Threads)
if(Threads_FOUND)
  if(CMAKE_USE_PTHREADS_INIT)
    set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} "-pthread")
  endif()
  foreach(PTHREAD_TEST HAVE_PTHREAD_PRIO_INHERIT PTHREAD_CREATE_UNDETACHED)
    check_c_source_compiles("
        #include <pthread.h>
        int main ()
        {
          int attr = ${PTHREAD_TEST}; return attr;
          return 0;
        }
      " ${PTHREAD_TEST}
    )
  endforeach()
  set(WEBP_DEP_LIBRARIES ${WEBP_DEP_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endif()
set(WEBP_USE_THREAD ${Threads_FOUND})

# TODO: this seems unused, check with autotools.
set(LT_OBJDIR ".libs/")

# Only useful for vwebp, so useless for now.
# find_package(OpenGL)
# set(WEBP_HAVE_GL ${OPENGL_FOUND})
# set(WEBP_DEP_INCLUDE_DIRS ${WEBP_DEP_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIRS})
# set(WEBP_DEP_LIBRARIES ${WEBP_DEP_LIBRARIES} ${OPENGL_LIBRARIES})

# Find the standard C math library.
find_library(MATH_LIBRARY NAMES m)
if(MATH_LIBRARY)
  set(WEBP_DEP_LIBRARIES ${WEBP_DEP_LIBRARIES} ${MATH_LIBRARY})
endif()

# Find the standard image libraries.
foreach(I_LIB PNG JPEG TIFF GIF)
  find_package(${I_LIB})
  set(WEBP_HAVE_${I_LIB} ${${I_LIB}_FOUND})
  set(WEBP_DEP_LIBRARIES
    ${WEBP_DEP_LIBRARIES}
    ${${I_LIB}_LIBRARIES}
  )
  set(WEBP_DEP_INCLUDE_DIRS
    ${WEBP_DEP_INCLUDE_DIRS}
    ${${I_LIB}_INCLUDE_DIRS}
  )
endforeach()

## Check for specific headers.
include(CheckIncludeFiles)
check_include_files("stdlib.h;stdarg.h;string.h;float.h" STDC_HEADERS)
check_include_files(dlfcn.h HAVE_DLFCN_H)
check_include_files(GLUT/glut.h HAVE_GLUT_GLUT_H)
check_include_files(GL/glut.h HAVE_GL_GLUT_H)
check_include_files(inttypes.h HAVE_INTTYPES_H)
check_include_files(memory.h HAVE_MEMORY_H)
check_include_files(OpenGL/glut.h HAVE_OPENGL_GLUT_H)
check_include_files(shlwapi.h HAVE_SHLWAPI_H)
check_include_files(stdint.h HAVE_STDINT_H)
check_include_files(stdlib.h HAVE_STDLIB_H)
check_include_files(strings.h HAVE_STRINGS_H)
check_include_files(string.h HAVE_STRING_H)
check_include_files(sys/stat.h HAVE_SYS_STAT_H)
check_include_files(sys/types.h HAVE_SYS_TYPES_H)
check_include_files(unistd.h HAVE_UNISTD_H)
check_include_files(wincodec.h HAVE_WINCODEC_H)
check_include_files(windows.h HAVE_WINDOWS_H)

## Check for CPU architectures.
foreach(I_TEST RANGE 2)
  set(WEBP_FLAGS "SSE2;SSE41;AVX2")
  list(GET WEBP_FLAGS ${I_TEST} WEBP_FLAG)
  if(MSVC)
    set(SIMD_COMPILE_FLAGS "/arch:SSE2;/arch:SSE41;/arch:AVX2")
    list(GET SIMD_COMPILE_FLAGS ${I_TEST} SIMD_COMPILE_FLAG)
  else()
    set(SIMD_COMPILE_FLAGS "-msse2;-msse4.1;-mavx2")
    list(GET SIMD_COMPILE_FLAGS ${I_TEST} SIMD_COMPILE_FLAG)
  endif()
  set(CMAKE_REQUIRED_FLAGS ${SIMD_COMPILE_FLAG})
  check_c_source_compiles("
      #include \"${CMAKE_CURRENT_LIST_DIR}/src/dsp/dsp.h\"
      in main(void) {
        #if !defined(WEBP_USE_${WEBP_FLAG})
        this is not valid code
        #endif
        return 0;
      }
    "
    WEBP_HAVE_${WEBP_FLAG}
  )
  if(WEBP_HAVE_${WEBP_FLAG})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SIMD_COMPILE_FLAG}")
  endif()
endforeach()

## Define extra info.
set(PACKAGE ${PROJECT_NAME})
set(PACKAGE_NAME ${PROJECT_NAME})

# Read from configure.ac.
file(READ ./configure.ac CONFIGURE_AC)
string(REGEX MATCHALL "\\[([0-9a-z\\.:/]*)\\]"
  CONFIGURE_AC_PACKAGE_INFO ${CONFIGURE_AC}
)
function(strip_bracket VAR)
  string(LENGTH ${${VAR}} TMP_LEN)
  math(EXPR TMP_LEN ${TMP_LEN}-2)
  string(SUBSTRING ${${VAR}} 1 ${TMP_LEN} TMP_SUB)
  set(${VAR} ${TMP_SUB} PARENT_SCOPE)
endfunction()

list(GET CONFIGURE_AC_PACKAGE_INFO 1 PACKAGE_VERSION)
strip_bracket(PACKAGE_VERSION)
list(GET CONFIGURE_AC_PACKAGE_INFO 2 PACKAGE_BUGREPORT)
strip_bracket(PACKAGE_BUGREPORT)
list(GET CONFIGURE_AC_PACKAGE_INFO 3 PACKAGE_URL)
strip_bracket(PACKAGE_URL)

# Build more info.
set(PACKAGE_STRING "${PACKAGE_NAME} ${PACKAGE_VERSION}")
set(PACKAGE_TARNAME ${PACKAGE_NAME})
set(VERSION ${PACKAGE_VERSION})

## Generate the config.h header in sources (ideally, we would copy it to the
# CMAKE_CURRENT_BINARY_DIR build dir, but we reference that file relatively
# with ../webp/config.h).
configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/config.h.in
  ${CMAKE_CURRENT_SOURCE_DIR}/src/webp/config.h)
add_definitions(-DHAVE_CONFIG_H)


################################################################################
# WebP source files.
# Read the Makefile.am to get the source files.
set(WEBP_SRCS "")

function(parse_Makefile_am FOLDER WEBP_SRCS)
  file(READ ${FOLDER}/Makefile.am MAKEFILE_AM)
  string(REGEX MATCHALL "_SOURCES \\+= [^\n]*"
    FILES_PER_LINE ${MAKEFILE_AM}
  )
  set(SRCS ${WEBP_SRCS})
  foreach(FILES ${FILES_PER_LINE})
    string(SUBSTRING ${FILES} 12 -1 FILES)
    string(REGEX MATCHALL "[0-9a-z\\._]+"
      FILES ${FILES}
    )
    foreach(FILE ${FILES})
      list(APPEND SRCS ${FOLDER}/${FILE})
    endforeach()
  endforeach()
  set(WEBP_SRCS ${SRCS} PARENT_SCOPE)
endfunction()

parse_Makefile_am(./src/dec "${WEBP_SRCS}")
parse_Makefile_am(./src/demux "${WEBP_SRCS}")
parse_Makefile_am(./src/dsp "${WEBP_SRCS}")
parse_Makefile_am(./src/enc "${WEBP_SRCS}")
parse_Makefile_am(./src/utils "${WEBP_SRCS}")

# Build the library.
add_definitions(-Wall)
include_directories(./src/ ${WEBP_DEP_INCLUDE_DIRS})
add_library(webp ${WEBP_SRCS})
target_link_libraries(webp ${WEBP_DEP_LIBRARIES})

# Build the executables if asked for.
if(WEBP_BUILD_CWEBP OR WEBP_BUILD_DWEBP)
  # Example utility library.
  set(exampleutil_SRCS
    ./examples/example_util.c
    ./examples/example_util.h
    ./examples/stopwatch.h)
  add_library(exampleutil ${exampleutil_SRCS})
  target_link_libraries(exampleutil webp ${WEBP_DEP_LIBRARIES})
endif()

if(WEBP_BUILD_DWEBP)
  # dwebp
  add_executable(dwebp
    ./examples/dwebp.c
    ./examples/stopwatch.h
  )
  target_link_libraries(dwebp webp exampleutil ${WEBP_DEP_LIBRARIES})
endif()

if(WEBP_BUILD_CWEBP)
  # cwebp
  add_executable(cwebp
    ./examples/cwebp.c
    ./examples/metadata.c
    ./examples/metadata.h
    ./examples/stopwatch.h
    ./examples/jpegdec.c
    ./examples/jpegdec.h
    ./examples/pngdec.c
    ./examples/pngdec.h
    ./examples/tiffdec.c
    ./examples/tiffdec.h
    ./examples/webpdec.c
    ./examples/webpdec.h
    ./examples/wicdec.c
    ./examples/wicdec.h)
  target_link_libraries(cwebp webp exampleutil ${WEBP_DEP_LIBRARIES})
endif()
