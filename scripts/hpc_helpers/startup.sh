#!/bin/bash

# Get the directory of the current script
#if bash use ${BASH_SOURCE[0]} else use ${0} (APPLE: zsh)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-${0}}")" && pwd)"
source "${SCRIPT_DIR}/../../configs/project.env"

cd "${REMOTE_WORKING_DIRECTORY}" || exit

alias lastout="'${REMOTE_WORKING_DIRECTORY}/scripts/hpc_helpers/lastOutput.sh'"
alias lasterr="'${REMOTE_WORKING_DIRECTORY}/scripts/hpc_helpers/lastError.sh'"
alias lastbench="'${REMOTE_WORKING_DIRECTORY}/scripts/hpc_helpers/lastBenchmark.sh'"
alias laststat="'${REMOTE_WORKING_DIRECTORY}/scripts/hpc_helpers/lastStat.sh'"
alias lastjobdir="source '${REMOTE_WORKING_DIRECTORY}/scripts/hpc_helpers/lastJobDir.sh'"
