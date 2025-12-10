#!/bin/bash

cd "${REMOTE_WORKING_DIRECTORY}" || { echo "[hpc_run.sh] Failed to change directory to ${REMOTE_WORKING_DIRECTORY}" >> ~/remote_run_error.log; exit 1; }


source "${HPC_CONFIG_FILE}" || { echo "[hpc_run.sh] Failed to source HPC config file ${HPC_CONFIG_FILE}" >> ~/remote_run_error.log; exit 1; }


LOGS_DIR="${HPC_JOB_LOGS_DIR}"
if [[ "${SANDBOX:-0}" == 1 ]]; then
    BINARY=$SANDBOX_BINARY
elif [[ "${BENCHMARK:-0}" == 1 ]]; then
    BINARY=$BENCHMARK_BINARY
    LOGS_DIR="${HPC_BENCHMARK_DIR}"
else
    BINARY=$PROJECT_BINARY
fi

echo "BINARY: ${BINARY}"

# Copy the application binary is a temporary location to avoid issues with concurrent runs modifying the binary while another job is running (or still in queue) f
TMP_DIR="${HPC_TMP_BINARIES_DIR}/$(uuidgen)"
mkdir -p "${TMP_DIR}" || { echo "[hpc_run.sh] Failed to create temporary directory ${TMP_DIR}" >> ~/remote_run_error.log; exit 1; }
cp "${CMAKE_BUILD_HPC_DEBUG}/${BINARY}" "${TMP_DIR}/" || { echo "[hpc_run.sh] Failed to copy binary to ${TMP_DIR}" >> ~/remote_run_error.log; exit 1; }
TMP_BINARY_PATH="${TMP_DIR}/$(basename -- "${BINARY}")"
export TMP_BINARY_PATH
chmod +x "${TMP_BINARY_PATH}"
echo "TMP_BINARY: ${TMP_BINARY_PATH}"

JOB_ARRAY_COMMAND=""
# Prepare job array command if needed
if [[ "${BENCHMARK:-0}" == 1 && "${JOB_ARRAY_SPEC:-}" != "" ]]; then
     JOB_ARRAY_COMMAND="-J ${JOB_ARRAY_SPEC}"
    echo "JOB ARRAY: '${JOB_ARRAY_SPEC}'"
fi

# Submit the job in held state (-h) to allow for post-processing of the job submission (e.g., logging the job ID)
# Preserving environment variables with -V (for example, REMOTE_WORKING_DIRECTORY)
job_id=$(qsub $JOB_ARRAY_COMMAND -h -V -o /dev/null -e /dev/null -l select="${CHUNKS}":ncpus="${N_CPUS}":mem="${MEMORY}":mpiprocs="${MPI_PROCESSES_PER_CHUNK}" -l walltime="${WALLTIME}" -l place="${PLACING}" -q "${QUEUE_TYPE}" "${REMOTE_WORKING_DIRECTORY}/scripts/run/run.pbs") || { echo "[hpc_run.sh] Failed to submit job" >> ~/remote_run_error.log; exit 1; }
echo "JOB_ID: ${job_id}"

numeric_job_id="${job_id%%.*}"
numeric_job_id="${numeric_job_id%%\[*}"

# Log the job configuration for future reference
mkdir -p "${LOGS_DIR}/${numeric_job_id}" || { echo "[hpc_run.sh] Failed to create job log dir ${LOGS_DIR}/${numeric_job_id}" >> ~/remote_run_error.log; exit 1; }
cp "${HPC_CONFIG_FILE}" "${LOGS_DIR}/${numeric_job_id}/job_config.log" || { echo "[hpc_run.sh] Failed to copy config to ${LOGS_DIR}/${numeric_job_id}" >> ~/remote_run_error.log; exit 1; }
chmod 444 "${LOGS_DIR}/${numeric_job_id}/job_config.log"

#FIXME: if benchmark execution, HPC helpers script don't work
# TODO: store dir of last job for easy access
mkdir -p "${HPC_JOB_LOGS_DIR}/last_job"
echo "${numeric_job_id}" > "${HPC_JOB_LOGS_DIR}/last_job/job_id"


# Release the queued job now that this submission script has completed
qrls "${job_id}" || { echo "[hpc_run.sh] Failed to release job ${job_id}" >> ~/remote_run_error.log; exit 1; }
echo "JOB_RELEASED: ${job_id}"

#TODO cleanup job after the previous job is done
## Cleanup temporary binary
#if [[ -n "${TMP_BINARY_PATH}" ]]; then
#    parent_dir="$(dirname "${TMP_BINARY_PATH}")"
#    rm -f "${TMP_BINARY_PATH}"
#    if [[ -n "$parent_dir" && "$parent_dir" != "/" && "$parent_dir" != "." ]]; then
#        rm -rf "$parent_dir"
#    fi
#fi

# TODO: if benchmark, launch data aggregation script after job completion (using qsub -W depend=afterok:<jobid>)