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


for m in "$@"; do
    nf="z$m";
    # Determine the max number of orbits.
    # NB: for now we assume that if one type is present (urs/sat/su), it is present for all orbits.
    d=0; until ! [[ -f "${DATA_DIR}/${nf}/${nf}_d$(($d+1)).urs" ]]; do d=$(($d+1)); done;
    # Determine whether we have precomputed S-units
    if [[ -f "${DATA_DIR}/${nf}/${nf}_d$d.su" ]];  then su="su=true";   else su="su=false"; fi;
    # Determine whether we have saturated elements
    if [[ -f "${DATA_DIR}/${nf}/${nf}_d$d.sat" ]]; then sat="sat=true"; else sat="sat=false"; fi;

    for di in `seq 1 1 $d`; do
        echo "Computing GSO's of Tw-PHS lattices for Q(z$m) [orbs=#$di,$sat,$su] (raw/lll/bkz-40)";
        ${EXE_DIR}/gso.sage ${DATA_DIR} $nf $di $sat $su 1>${LOGS_DIR}/${nf}_${di}.gsolog 2>&1 &
    done;
done;

exit 0;
