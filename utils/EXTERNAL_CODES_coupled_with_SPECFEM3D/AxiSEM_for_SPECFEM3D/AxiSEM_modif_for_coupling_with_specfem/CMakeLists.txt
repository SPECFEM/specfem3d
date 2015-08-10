CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (AxiSEM)
ENABLE_LANGUAGE (Fortran)

# Add a module search path for custom modules
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)

# Enable testing
INCLUDE (CTest)

SET(AXISEM_SOURCE_DIR ${CMAKE_SOURCE_DIR})

# Find MPI installation
FIND_PACKAGE(MPI)
IF(DEFINED MPI_Fortran_FOUND)
    ADD_DEFINITIONS(${MPI_Fortran_COMPILE_FLAGS})
    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
ENDIF(DEFINED MPI_Fortran_FOUND)

# Determine if compiler supports OpenMP
FIND_PACKAGE(OpenMP_Fortran)

# Check if NetCDF for Fortran90 is installed
SET (NETCDF_F90 "YES")
FIND_PACKAGE(NetCDF)
IF (NETCDF_FOUND)
    INCLUDE_DIRECTORIES(${NETCDF_INCLUDES})
ENDIF (NETCDF_FOUND)

ADD_SUBDIRECTORY("MESHER")
ADD_SUBDIRECTORY("SOLVER")

SET(TEST_NAMES test_01_explosion test_02_quadpole)
IF (NETCDF_FOUND)
    SET(TEST_NAMES test_01_explosion test_02_quadpole test_03_netcdf)
ENDIF (NETCDF_FOUND)

FOREACH(TEST_NAME ${TEST_NAMES})
    SET(TEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/tests/${TEST_NAME}/)
    SET(TEST_SCRIPT_DIR ${CMAKE_SOURCE_DIR}/cmake/)
    SET(TEST_DATA_DIR ${CMAKE_SOURCE_DIR}/TESTING/automated/${TEST_NAME}/)
    FILE(MAKE_DIRECTORY ${TEST_DIR})
    FILE(MAKE_DIRECTORY ${TEST_DIR}/Diags/)
    FILE(COPY ${CMAKE_SOURCE_DIR}/SOLVER/CMTSOLUTION.TEMPLATE DESTINATION ${TEST_DIR})

    # Run the test input file setup
    ADD_TEST(NAME ${TEST_NAME}_setup
        WORKING_DIRECTORY ${TEST_DIR}
        COMMAND ${TEST_SCRIPT_DIR}/ctest_parse_params.py ${TEST_DATA_DIR}/inpython.cfg ${TEST_DATA_DIR}/STATIONS)

    # Run the mesher
    ADD_TEST(NAME ${TEST_NAME}_mesh
        WORKING_DIRECTORY ${TEST_DIR}
        COMMAND ${CMAKE_CURRENT_BINARY_DIR}/MESHER/xmesh)
    SET_TESTS_PROPERTIES(${TEST_NAME}_mesh PROPERTIES DEPENDS ${TEST_NAME}_setup)

    # Move the meshes to the appropriate directory
    ADD_TEST(NAME ${TEST_NAME}_move_mesh
        WORKING_DIRECTORY ${TEST_DIR}
        COMMAND ${TEST_SCRIPT_DIR}/ctest_move_mesh.csh ${TEST_DIR})
    SET_TESTS_PROPERTIES(${TEST_NAME}_move_mesh PROPERTIES DEPENDS ${TEST_NAME}_mesh)

    # Run the axisem program
    # Determine the number of processors to use
    EXECUTE_PROCESS(COMMAND ${TEST_SCRIPT_DIR}/ctest_calc_nprocs.py ${TEST_DATA_DIR} OUTPUT_VARIABLE NPROCS)
    ADD_TEST(NAME ${TEST_NAME}_solver
        WORKING_DIRECTORY ${TEST_DIR}
        COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NPROCS} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/SOLVER/axisem ${MPIEXEC_POSTFLAGS})
    SET_TESTS_PROPERTIES(${TEST_NAME}_solver PROPERTIES DEPENDS ${TEST_NAME}_move_mesh)

    # Run the post processor param file editor
    ADD_TEST(NAME ${TEST_NAME}_post_proc_setup
        WORKING_DIRECTORY ${TEST_DIR}
        COMMAND ${TEST_SCRIPT_DIR}/ctest_prepare_ppp.py ${TEST_DATA_DIR}/inpython.cfg)
    SET_TESTS_PROPERTIES(${TEST_NAME}_post_proc_setup PROPERTIES DEPENDS ${TEST_NAME}_solver)

    # Run the xpost processor
    ADD_TEST(NAME ${TEST_NAME}_post_proc
        WORKING_DIRECTORY ${TEST_DIR}
        COMMAND ${TEST_SCRIPT_DIR}/ctest_post_proc.csh)
    SET_TESTS_PROPERTIES(${TEST_NAME}_post_proc PROPERTIES DEPENDS ${TEST_NAME}_post_proc_setup)

    # Check the results with the baseline values
    ADD_TEST(NAME ${TEST_NAME}_check_results
        WORKING_DIRECTORY ${TEST_DIR}
        COMMAND ${TEST_SCRIPT_DIR}/ctest_check_results.py ${TEST_DATA_DIR}/inpython.cfg ${TEST_DATA_DIR}/ref_data/)
    SET_TESTS_PROPERTIES(${TEST_NAME}_check_results PROPERTIES DEPENDS ${TEST_NAME}_post_proc)
ENDFOREACH()

