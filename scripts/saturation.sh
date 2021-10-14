#!/bin/bash

# cf ../data/list_ms

# Root dir
ROOT_DIR=$(dirname $(dirname $(readlink -f $0)));
# Subdirs: Exe/Data/logs
EXE_DIR="${ROOT_DIR}/scripts";
DATA_DIR="${ROOT_DIR}/data"; 
LOGS_DIR="${ROOT_DIR}/logs";


# Just check that parent folders are indeed where they should be
[[ ! -d ${DATA_DIR} ]] && {
    echo -e "\x1b[31m[Err]\x1b[0m Data directory ${DATA_DIR} does not exist.";
    exit 1;
};

[[ ! -d ${LOGS_DIR} ]] && {
    echo -e "\x1b[31m[Err]\x1b[0m Logs directory ${LOGS_DIR} does not exist.";
    exit 1;
};


# Saturate from URS sets
for m in "$@"; do
    nf="z$m";
    d=0; until ! [[ -f "${DATA_DIR}/${nf}/${nf}_d$(($d+1)).urs" ]]; do d=$(($d+1)); done;
    echo "2-Saturating real/Stickelberger generators for Q(z$m) [#orbs=$d]";
    ${EXE_DIR}/saturation.sage ${DATA_DIR} $nf $d 1>${LOGS_DIR}/${nf}.satlog 2>&1 &
done;

exit 0;
