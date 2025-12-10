# FIND AND LINK EXTERNAL PACKAGE DEPENDENCIES

# ------------------------------------------------------------------------------
# MPI (always C)
# ------------------------------------------------------------------------------
function(find_and_link_mpi_for_target tgt SCOPE)
    find_package(MPI REQUIRED)

    get_target_property(sources ${tgt} SOURCES)
    set(has_cpp FALSE)
    set(has_c FALSE)
    foreach(src IN LISTS sources)
        if(src MATCHES "\\.cpp$|\\.cxx$|\\.cc$")
            set(HAS_CXX TRUE)
        endif()
        if(src MATCHES "\\.c$")
            set(HAS_C TRUE)
        endif ()
    endforeach()

    if(HAS_CXX)
        colored_message(${Blue}  "Link MPI::MPI_CXX to target ${tgt}")
        target_link_libraries(${tgt} ${SCOPE} MPI::MPI_CXX)
    endif()
    if(HAS_C)
        colored_message(${Blue}  "Link MPI::MPI_C to target ${tgt}")
        target_link_libraries(${tgt} ${SCOPE} MPI::MPI_C)
    endif()
endfunction()


# ------------------------------------------------------------------------------
# OpenMP (C/C++ detection)
# ------------------------------------------------------------------------------
# Provide helper function to link OpenMP per-target based on language (C/C++)
function(find_and_link_openmp_for_target tgt SCOPE)
    # Try normal find_package first
    find_package(OpenMP)

    # Create interface targets for Apple fallback if needed
    if(APPLE AND (NOT OpenMP_C_FOUND OR NOT OpenMP_CXX_FOUND))
        message(STATUS "OpenMP not found via find_package, trying Homebrew libomp")

        execute_process(
                COMMAND brew --prefix libomp
                OUTPUT_VARIABLE HOMEBREW_LIBOMP_PREFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
        )

        if(HOMEBREW_LIBOMP_PREFIX)
            message(STATUS "Using Homebrew libomp at ${HOMEBREW_LIBOMP_PREFIX}")

            add_library(OpenMP::OpenMP_C INTERFACE)
            add_library(OpenMP::OpenMP_CXX INTERFACE)

            foreach(target IN ITEMS OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
                target_include_directories(${target} INTERFACE ${HOMEBREW_LIBOMP_PREFIX}/include)
                target_compile_options(${target} INTERFACE -Xpreprocessor -fopenmp -Wall -Wextra -Wpedantic)
                target_link_options(${target} INTERFACE -L${HOMEBREW_LIBOMP_PREFIX}/lib -lomp)
            endforeach()
        else()
            message(FATAL_ERROR "OpenMP not found and Homebrew libomp not detected")
        endif()

    elseif(!APPLE)
        message(FATAL_ERROR "OpenMP not found on this system")
    endif()

    get_target_property(sources ${tgt} SOURCES)
    set(has_cpp FALSE)
    set(has_c FALSE)
    foreach(src IN LISTS sources)
        if(src MATCHES "\\.cpp$|\\.cxx$|\\.cc$")
            set(HAS_CXX TRUE)
        endif()
        if(src MATCHES "\\.c$")
            set(HAS_C TRUE)
        endif ()
    endforeach()

    if(HAS_CXX)
        colored_message(${Blue} "Link OpenMP::OpenMP_CXX to target ${tgt}")
        target_link_libraries(${tgt} ${SCOPE} OpenMP::OpenMP_CXX)
    endif()
    if(HAS_C)
        colored_message(${Blue}  "Link OpenMP::OpenMP_C to target ${tgt}")
        target_link_libraries(${tgt} ${SCOPE} OpenMP::OpenMP_C)
    endif()

endfunction()
