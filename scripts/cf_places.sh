#!/bin/bash

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


# Places
for m in "$@"; do
    nf="z$m";
    folder="${DATA_DIR}/$nf";
    echo "Field 'Q(z$m)':";
    echo -e "\tCreating $folder"; mkdir -p "$folder";
    echo -e "\tGet places";
    ${EXE_DIR}/cf_places.sage ${DATA_DIR} $nf 1>${LOGS_DIR}/$nf.fblog 2>&1 &
done;

exit 0;
