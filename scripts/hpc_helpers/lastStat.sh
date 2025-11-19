#!/bin/bash

if [[ -n "${STARTUP_HPC_HELPER_SCRIPTS:-}" ]]; then
    echo "You have to source the startup.sh script before running this script."
    exit 1
fi

cd "${REMOTE_WORKING_DIRECTORY}" || { echo "[lastStat.sh] Failed to change directory to ${REMOTE_WORKING_DIRECTORY}"; exit 1; }

LAST_JOB_ID="$(cat "${HPC_JOB_LOGS_DIR}/last_job/job_id")"

qstat -H "${LAST_JOB_ID}"