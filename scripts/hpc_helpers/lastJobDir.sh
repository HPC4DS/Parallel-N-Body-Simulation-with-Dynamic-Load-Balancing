#!/bin/bash

echo "Changing directory to the last job's log directory."
if [[ -n "${STARTUP_HPC_HELPER_SCRIPTS:-}" ]]; then
    echo "You have to source the startup.sh script before running this script."
    exit 1
fi

cd "${REMOTE_WORKING_DIRECTORY}" || { echo "[lastLogDir.sh] Failed to change directory to ${REMOTE_WORKING_DIRECTORY}"; exit 1; }


LAST_JOB_ID="$(cat ./HPC/logs/last_job/job_id)"
cd "./HPC/logs/${LAST_JOB_ID}" || { echo "[lastLogDir.sh] Failed to change directory to ./HPC/logs/${LAST_JOB_ID}"; exit 1; }

echo "Changed directory to: $(pwd)"