#!/bin/bash

if [[ -n "${STARTUP_HPC_HELPER_SCRIPTS:-}" ]]; then
    echo "You have to source the startup.sh script before running this script."
    exit 1
fi

cd "${REMOTE_WORKING_DIRECTORY}" || { echo "[$(basename "$0")] Failed to change directory to ${REMOTE_WORKING_DIRECTORY}"; exit 1; }

LAST_JOB_ID="$(cat "${HPC_JOB_LOGS_DIR}/last_job/job_id")"
LAST_BENCHMARK_FILE="${HPC_JOB_LOGS_DIR}/${LAST_JOB_ID}/benchmark.log"

printf "\n==========================================================\n"
printf "BENCHMARK OF JOB: %s" "${LAST_JOB_ID}"
printf "\n==========================================================\n\n"

cat "${LAST_BENCHMARK_FILE}"
