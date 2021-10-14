#!/bin/bash

# Sorted by degree to simplify workflow
# 23  35  39  45  52  56  72  84  29  31  51  64  68  80  96  120 # n=32
# 37  57  63  76  108 41  55  75  88  100 132 43  49  69  92  47  # n=46
# 65  104 105 112 140 144 156 168 180 53  81  87  116 59  61  77  # n=60
# 93  99  124 85  128 136 160 192 204 240 67  71                  # n=70
# 73  91  95  111 117 135 148 152 216 228 252                     # n=72 (all)
# 79                                                              # n=78  <-----  SUnits OK (.5 day)

# A priori, the following conductors fail:
# 123 164 165 176 200 220 264 300                                 # n=80 (all)
# etc


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


# S-Units (Raw)
for m in "$@"; do
    nf="z$m";
    d=0; until ! [[ -f "${DATA_DIR}/${nf}/${nf}_d$(($d+1)).fb" ]]; do d=$(($d+1)); done;
    echo "Launching S-Units for Q(z$m) [#orbs=$d]";
    magma -b data:=${DATA_DIR} nf:=$nf d:=$d sunits.magma 1>${LOGS_DIR}/${nf}.sulog 2>&1 & # Don't use ${EXE_DIR}/sunits.magma as this script doesn't work if not launched from scripts/
done;

exit 0;
