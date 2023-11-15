cmake_minimum_required (VERSION 2.6)

### Project Configuration
project(sep_topo)

include(CheckCCompilerFlag)
check_c_compiler_flag(-std=gnu99 HAS_STD_GNU99)

if(HAS_STD_GNU99)
  add_definitions(-std=gnu99)
endif()
file(MAKE_DIRECTORY ./bin)
set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

include_directories(
  ./inc
)

set(sources
    src/sep_header.c
    src/convert_topo.c
    src/parse_sep.c)

add_executable(
  create_ascii_bathy_from_sep_vs
  src/create_ascii_bathy_from_sep_vs.c
  ${sources}
)
