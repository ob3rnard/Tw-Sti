#!/bin/bash

# cf ../data/list_ms_pcmp

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


# Multiprocess param
NCORES=4; # This run will be launched on NCORES (for all targets, all fields successively)

# Compute approached log-S-Unit lattices
for m in "$@"; do
    nf="z$m";
    # Determine the max number of orbits.
    durs=0; until ! [[ -f "${DATA_DIR}/${nf}/${nf}_d$(($durs+1)).urs" ]]; do durs=$(($durs+1)); done;
    dsat=0; until ! [[ -f "${DATA_DIR}/${nf}/${nf}_d$(($dsat+1)).sat" ]]; do dsat=$(($dsat+1)); done;
    dsu=0;  until ! [[ -f "${DATA_DIR}/${nf}/${nf}_d$(($dsu+1)).su"   ]]; do dsu=$(($dsu+1));   done;

    # Experience gives best results for durs=dsat=1, dsu=dmax
    durs=1; dsat=1;
    
    echo "Simulate IdSVP Solve using Tw-PHS for Q(z$m) [durs=$durs,dsat=$dsat,dsu=$dsu]";
    for suset in "urs" "sat" "su"; do
        dname="d${suset}"; dm=${!dname}; # Now, $dm=$durs or $dsat or $dsu according to $suset
        for di in `seq 1 1 $dm`; do
            sage ${EXE_DIR}/approx_factor.sage ${DATA_DIR} ${NCORES} ${nf} ${di} set=${suset} 1>${LOGS_DIR}/${nf}_d${di}.aflog_${suset} 2>&1 ;
        done;
    done &
done ;

exit 0;
