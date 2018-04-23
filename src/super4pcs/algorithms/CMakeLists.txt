project (OpenGR-AlgoLib)

set(algo_ROOT "${SRC_DIR}/super4pcs/algorithms")

file(GLOB_RECURSE algo_sources ${algo_ROOT}/*.cc)
file(GLOB_RECURSE algo_headers ${algo_ROOT}/*.h ${algo_ROOT}/*.hpp)

add_library(opengr_algo ${algo_headers} ${algo_sources})
# we need that since opengr_accel might be header only library
if(OpenGR_USE_CHEALPIX)
    target_link_libraries(opengr_algo ${Chealpix_LIBS} )
endif(OpenGR_USE_CHEALPIX)
target_link_libraries(opengr_algo opengr_io opengr_accel)

set(targets_export_name "${PROJECT_NAME}Targets")
install( TARGETS opengr_algo
    EXPORT "${targets_export_name}"
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    INCLUDES DESTINATION include  )
