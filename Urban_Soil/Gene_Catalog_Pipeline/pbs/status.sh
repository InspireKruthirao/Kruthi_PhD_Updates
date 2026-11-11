#!/bin/bash
# Tell Snakemake whether a PBS job is running, finished OK, or failed.
info=$(qstat -xf "$1" 2>/dev/null) || { echo running; exit 0; }
state=$(awk '/job_state =/ {print $3}' <<< "$info")
code=$(awk '/Exit_status =/ {print $3}' <<< "$info")
if [[ $state != F ]]; then echo running
elif [[ $code == 0 ]]; then echo success
else echo failed; fi
