if(NOT TRIBOL_USE_MPI)
    message(FATAL_ERROR "TRIBOL_ENABLE_FUTURE requires MPI")
endif()

if(NOT BUILD_REDECOMP)
    message(FATAL_ERROR "TRIBOL_ENABLE_FUTURE requires BUILD_REDECOMP=ON")
endif()

set(tribol_future_headers
    README.md
    Config.hpp
    Concepts.hpp
    ContactDomain.hpp
    ContactRestrictionOperator.hpp
    ElementBlockAssembler.hpp
    SurfaceRestrictionOperator.hpp
    Execution.hpp
    InteractionBatches.hpp
    MfemBridge.hpp
    MfemContactOperator.hpp
    Policies.hpp
    Results.hpp
    RuntimeFactory.hpp
    Search.hpp
    SurfaceDiscretizations.hpp
    Types.hpp
    algorithms/CommonPlane.hpp
    algorithms/EnergyMortar.hpp
    algorithms/SingleMortar.hpp
    detail/AutomaticDifferentiation.hpp
    detail/ElementGeometry.hpp
    detail/Enzyme.hpp
    TribolFuture.hpp
    )

set(tribol_future_sources
    ContactDomain.cpp
    ContactRestrictionOperator.cpp
    ElementBlockAssembler.cpp
    SurfaceRestrictionOperator.cpp
    MfemBridge.cpp
    )

set(tribol_future_depends
    tribol_shared
    redecomp
    mfem
    axom::core
    axom::primal
    axom::spin
    blt::mpi
    )
blt_list_append(TO tribol_future_depends ELEMENTS umpire IF TRIBOL_USE_UMPIRE)
blt_list_append(TO tribol_future_depends ELEMENTS RAJA IF TRIBOL_USE_RAJA)
blt_list_append(TO tribol_future_depends ELEMENTS ClangEnzymeFlags IF TRIBOL_USE_ENZYME)

blt_add_library(
    NAME tribol_future
    SOURCES ${tribol_future_sources}
    HEADERS ${tribol_future_headers}
    DEPENDS_ON ${tribol_future_depends} ${tribol_device_depends}
    FOLDER tribol/future
    )

add_library(tribol::future ALIAS tribol_future)
set_target_properties(tribol_future PROPERTIES EXPORT_NAME tribol::future)
target_compile_features(tribol_future PUBLIC cxx_std_20)
target_compile_definitions(tribol_future PUBLIC
    TRIBOL_FUTURE_NATIVE_HO_MAX_ORDER=${TRIBOL_FUTURE_NATIVE_HO_MAX_ORDER})
target_include_directories(tribol_future PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>
    )

if(ENABLE_HIP)
    target_compile_options(tribol_future PRIVATE $<$<CONFIG:Debug>:-ggdb>)
endif()

install(TARGETS tribol_future
        EXPORT tribol-targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include)

foreach(tribol_future_header ${tribol_future_headers})
    get_filename_component(tribol_future_base_dir ${tribol_future_header} DIRECTORY)
    install(FILES ${tribol_future_header}
            DESTINATION "include/tribol/future/${tribol_future_base_dir}")
endforeach()
