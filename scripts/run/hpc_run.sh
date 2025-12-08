#!/bin/bash

cd "${REMOTE_WORKING_DIRECTORY}" || { echo "[hpc_run.sh] Failed to change directory to ${REMOTE_WORKING_DIRECTORY}" >> ~/remote_run_error.log; exit 1; }


source "${HPC_CONFIG_FILE}" || { echo "[hpc_run.sh] Failed to source HPC config file ${HPC_CONFIG_FILE}" >> ~/remote_run_error.log; exit 1; }


if [[ "${SANDBOX:-0}" == 1 ]]; then
    BINARY=$SANDBOX_BINARY
elif [[ "${BENCHMARK:-0}" == 1 ]]; then
    BINARY=$BENCHMARK_BINARY
else
    BINARY=$PROJECT_BINARY
fi

# Copy the application binary is a temporary location to avoid issues with concurrent runs modifying the binary while another job is running (or still in queue)

TMP_DIR="${HPC_TMP_BINARIES_DIR}/$(uuidgen)"
mkdir -p "${TMP_DIR}" || { echo "[hpc_run.sh] Failed to create temporary directory ${TMP_DIR}" >> ~/remote_run_error.log; exit 1; }
cp "${CMAKE_BUILD_HPC_DEBUG}/${BINARY}" "${TMP_DIR}/" || { echo "[hpc_run.sh] Failed to copy binary to ${TMP_DIR}" >> ~/remote_run_error.log; exit 1; }
TMP_BINARY_PATH="${TMP_DIR}/$(basename -- "${BINARY}")"
export TMP_BINARY_PATH
chmod +x "${TMP_BINARY_PATH}"
echo "TMP_BINARY: ${TMP_BINARY_PATH}"

# Preserving environment variables with -V (for example, REMOTE_WORKING_DIRECTORY)
job_id=$(qsub -V  -o /dev/null -e /dev/null -l select="${CHUNKS}":ncpus="${N_CPUS}":mem="${MEMORY}":mpiprocs="${MPI_PROCESSES_PER_CHUNK}" -l walltime="${WALLTIME}" -l place="${PLACING}" -q "${QUEUE_TYPE}" "${REMOTE_WORKING_DIRECTORY}/scripts/run/run.pbs")  || { echo "[hpc_run.sh] Failed to submit job" >> ~/remote_run_error.log; exit 1; }
echo "JOB_ID: ${job_id}"

numeric_job_id="${job_id%%.*}"

mkdir -p "${HPC_JOB_LOGS_DIR}/${numeric_job_id}"
cp "${HPC_CONFIG_FILE}" "${HPC_JOB_LOGS_DIR}/${numeric_job_id}/job_config.log"
chmod 444 "${HPC_JOB_LOGS_DIR}/${numeric_job_id}/job_config.log"

mkdir -p "${HPC_JOB_LOGS_DIR}/last_job"
echo "${numeric_job_id}" > "${HPC_JOB_LOGS_DIR}/last_job/job_id"
