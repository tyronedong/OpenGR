project (OpenGR-utilsLib)

set(utils_ROOT "${SRC_DIR}/super4pcs/utils")

#############################################

set(utils_INCLUDE
    ${utils_ROOT}/timer.h
    ${utils_ROOT}/geometry.h
    ${utils_ROOT}/disablewarnings.h
    ${utils_ROOT}/logger.h)


if(CMAKE_VERSION VERSION_GREATER 3.1)
    # use library target, to make accessible to test routines
    add_library(opengr_utils INTERFACE)
    target_sources(opengr_utils INTERFACE ${utils_INCLUDE})
    #
    # HACK: have the files showing in the IDE, under the name 'opengr_utils_IDE'
    # Source: https://stackoverflow.com/questions/27039019/how-to-have-cmake-show-headers-that-are-not-part-of-any-binary-target-in-the-ide/29214327#29214327
    #
    add_custom_target(opengr_utils_IDE SOURCES ${utils_INCLUDE})
else ()
    add_library(opengr_utils STATIC ${utils_INCLUDE})
    set_target_properties(opengr_utils PROPERTIES LINKER_LANGUAGE CXX)
endif(CMAKE_VERSION VERSION_GREATER 3.1)

