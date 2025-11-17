#!/bin/zsh

cd "${LOCAL_WORKING_DIRECTORY}" || { echo "[local_run.sh] Failed to change directory to ${LOCAL_WORKING_DIRECTORY}" >> ~/local_run_error.log; exit 1; }
source "${LOCAL_CONFIG_FILE}" || { echo "[local_run.sh] Failed to source environment file ${LOCAL_ENV_FILE}" >> ~/local_run_error.log; exit 1; }

mpiexec -n "${LOCAL_MPI_PROCESSES}" "${CMAKE_BUILD_LOCAL_DEBUG}/${PROJECT_BINARY}"
