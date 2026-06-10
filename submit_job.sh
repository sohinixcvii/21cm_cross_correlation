#!/bin/bash
# =============================================================================
# submit_job.sh — Python run script with sendmail email notification
# =============================================================================
# Usage:
#   bash submit_job.sh
# =============================================================================

# ── User settings ─────────────────────────────────────────────────────────────
EMAIL_TO="sohinidutta97@gmail.com"
JOB_NAME="21cm_lightcone"
PYTHON_SCRIPT="run_simulation.py"
CONDA_ENV="21cmfast"

# ── Create outputs directory ──────────────────────────────────────────────────
mkdir -p outputs

LOG_FILE="outputs/${JOB_NAME}_$(date +%Y%m%d_%H%M%S).log"

# ── Start timing ──────────────────────────────────────────────────────────────
START_EPOCH=$(date +%s)
START_TIME=$(date)

echo "Starting ${JOB_NAME} at ${START_TIME}" | tee -a "$LOG_FILE"
echo "Running on host: $(hostname)" | tee -a "$LOG_FILE"
echo "Python script: ${PYTHON_SCRIPT}" | tee -a "$LOG_FILE"
echo "Log file: ${LOG_FILE}" | tee -a "$LOG_FILE"
echo "----------------------------------------" | tee -a "$LOG_FILE"

# ── Activate conda environment ────────────────────────────────────────────────
# Uncomment this if conda activate fails:
# source "$(conda info --base)/etc/profile.d/conda.sh"

conda activate "$CONDA_ENV"

# ── Run simulation and collect CPU usage ──────────────────────────────────────
/usr/bin/time -p python "$PYTHON_SCRIPT" >> "$LOG_FILE" 2>&1

EXIT_CODE=$?

# ── Finish timing ─────────────────────────────────────────────────────────────
END_EPOCH=$(date +%s)
END_TIME=$(date)

RUNTIME_SECONDS=$((END_EPOCH - START_EPOCH))
RUNTIME_HMS=$(printf '%02d:%02d:%02d' \
    $((RUNTIME_SECONDS / 3600)) \
    $(((RUNTIME_SECONDS % 3600) / 60)) \
    $((RUNTIME_SECONDS % 60)))

USER_CPU=$(grep "^user " "$LOG_FILE" | tail -1 | awk '{print $2}')
SYS_CPU=$(grep "^sys " "$LOG_FILE" | tail -1 | awk '{print $2}')

if [ -z "$USER_CPU" ]; then USER_CPU=0; fi
if [ -z "$SYS_CPU" ]; then SYS_CPU=0; fi

CPU_SECONDS=$(awk -v u="$USER_CPU" -v s="$SYS_CPU" 'BEGIN {print u+s}')
CPU_HOURS=$(awk -v c="$CPU_SECONDS" 'BEGIN {printf "%.4f", c/3600}')

if [ $EXIT_CODE -eq 0 ]; then
    STATUS="COMPLETE"
    SUBJECT="[${JOB_NAME}] COMPLETE on $(hostname)"
else
    STATUS="FAILED"
    SUBJECT="[${JOB_NAME}] FAILED on $(hostname)"
fi

# ── Email body ────────────────────────────────────────────────────────────────
EMAIL_BODY=$(cat <<EOF
Job name:      ${JOB_NAME}
Status:        ${STATUS}
Host:          $(hostname)

Start time:    ${START_TIME}
Finish time:   ${END_TIME}
Total runtime: ${RUNTIME_HMS}
CPU hours:     ${CPU_HOURS}

Python script: ${PYTHON_SCRIPT}
Log file:      ${LOG_FILE}
Exit code:     ${EXIT_CODE}

Last 30 lines of log:
------------------------------------------------------------
$(tail -30 "$LOG_FILE")
EOF
)

# ── Send email notification using sendmail ────────────────────────────────────
{
cat <<EOF
To: ${EMAIL_TO}
Subject: ${SUBJECT}
Content-Type: text/plain; charset=UTF-8

${EMAIL_BODY}
EOF
} | sendmail -t

# ── Final console/log summary ─────────────────────────────────────────────────
echo "----------------------------------------" | tee -a "$LOG_FILE"
echo "Finished with status: ${STATUS}" | tee -a "$LOG_FILE"
echo "Finished at ${END_TIME}" | tee -a "$LOG_FILE"
echo "Total runtime: ${RUNTIME_HMS}" | tee -a "$LOG_FILE"
echo "CPU hours: ${CPU_HOURS}" | tee -a "$LOG_FILE"

exit $EXIT_CODE