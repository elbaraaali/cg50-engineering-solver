cmake_minimum_required(VERSION 3.18)
project(CG50EngineeringSolver LANGUAGES C ASM)

include(GenerateG3A)
include(Fxconv)
find_package(Gint 2.9 REQUIRED)

set(SOURCES
  src/main.c
  src/engine.c
)

add_executable(engsolve ${SOURCES})
target_compile_options(engsolve PRIVATE -Wall -Wextra -Werror -Os)
target_link_libraries(engsolve Gint::Gint)

if("${FXSDK_PLATFORM_LONG}" STREQUAL fxCG50)
  # Empty internal/menu overlay name is intentional for CG50 launch compatibility.
  generate_g3a(TARGET engsolve OUTPUT "ENGSOLV.g3a"
    NAME "" ICONS assets-cg/icon-uns.png assets-cg/icon-sel.png)
endif()
