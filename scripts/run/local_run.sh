#!/bin/zsh

cd "${LOCAL_WORKING_DIRECTORY}" || { echo "[local_run.sh] Failed to change directory to ${LOCAL_WORKING_DIRECTORY}" >> ~/local_run_error.log; exit 1; }
source "${LOCAL_CONFIG_FILE}" || { echo "[local_run.sh] Failed to source environment file ${LOCAL_ENV_FILE}" >> ~/local_run_error.log; exit 1; }

export JOB_ID="$(($(date +%s%N)/1000000))"

if [[ "${SANDBOX:-0}" == 1 ]]; then
    BINARY=$SANDBOX_BINARY
elif [[ "${BENCHMARK:-0}" == 1 ]]; then
    BINARY=$BENCHMARK_BINARY
else
    BINARY=$PROJECT_BINARY
fi

if [[ "${MPI_RUN:-0}" == 1 ]]; then
    mpiexec -n "${LOCAL_MPI_PROCESSES}" "${CMAKE_BUILD_LOCAL_DEBUG}/${BINARY}"
else
    "./${CMAKE_BUILD_LOCAL_DEBUG}/${BINARY}"
fi