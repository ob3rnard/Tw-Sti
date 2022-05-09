#!/usr/bin/env sage

# Code in support of ePrint:2021/1384
# Copyright 2021, Olivier Bernard
# GPL-3.0-only (see LICENSE file)

# Small convenient hack to simulate as if we were always executing from '<trunk>/'
import os
import sys
os.chdir(os.path.dirname(__file__) + "/../");
sys.path.append("./src/");
# --------------------------------------------------------------------------------------

from sage.all    import *
from pcmp_io     import *


W_PREC = 500; # Must be sufficient to handle discriminant


# --------------------------------------------------------------------------------------
# Data folder
data_dir = os.getcwd() + "/data";
conductors_file = data_dir + "/list_ms_pcmp";
su_sets  = ["cdw", "urs", "sat", "su"];
measures = ["afsup", "gf", "afinf", "hf"];


# --------------------------------------------------------------------------------------
# List of conductors
if (len(sys.argv) == 1):
    list_m = [];
    for line in open(conductors_file, 'r'):
        for m in line.split():
            list_m.append(int(m));
else:
    list_m = list(map(int, sys.argv[1:]));


# --------------------------------------------------------------------------------------
# Helping I/O functions
# Compute the mean value of AFs in the file for m, d orbs, set "urs"/... and measure "gf/..."
def file_get_mean_min(m, d, su_set, measure):
    af_file_name = data_dir + "/z{}/z{}_d{}_{}.{}".format(m, m, d, su_set, measure)
    if (not os.path.exists(af_file_name)):
        return RDF(+infinity);

    af_file  = open(af_file_name, "r");
    min_vals = [ min( map(float, _line.rstrip().split()[:4]) ) for _line in af_file
                 if not _line.startswith("#") ];

    return mean(min_vals);


# Find the max nb of orbits for a given S-units family
def find_max_nb_orbits(m, su_set):
    assert (su_set in su_sets);
    suffix = "urs" if su_set == "cdw" else su_set;
    
    d = 1;
    while True:
        suset_file = data_dir + "/z{}/z{}_d{}.{}".format(m, m, d, suffix);
        if os.path.exists(suset_file):
            d += 1;
        else:
            break;
    dmax = (d-1) if d > 1 else None;
    
    return dmax;


# Find, for a given S-units family, the number of orbits that minimizes the AF
def find_optimal_nb_orbits(m, su_set, dmax):
    assert (su_set in su_sets);
    if (dmax == None):
        return None;
    
    mean_afs = [ file_get_mean_min(m, _d, su_set, "gf") for _d in range(1, dmax+1) ];
    best_af, d_opt = min( (_af, _idx+1) for (_idx, _af) in enumerate(mean_afs) );

    return d_opt;


def af_write(m, measure, d_opt, dpw, afs):
    out_file_name = data_dir + "/z{}/z{}.{}".format(m, m, measure);
    with open(out_file_name, "w") as out_file:
        # Header
        head  = "# Measure:'{}' nf:'z{}'\n".format(measure, m);
        head += ("# Orbits:"+" d{}={}"*len(su_sets)+"\n").format(*sum(([_su, d_opt[_su]] for _su in su_sets), []));
        head += ("# n\tlow_dpw"+"\t\t{}"*len(su_sets)+"\n").format(*su_sets);
        out_file.write(head);

        # Data
        line = list(map(lambda x: str(x).replace("+infinity","NaN"),
                        [euler_phi(m), dpw[measure]] + afs));
        out_file.write(("{}"+"\t{:<13.13s}"*(len(line)-1)+"\n").format(*line));
 
    return;


# -------------------------------------------------------------------------------------------------
# Hermite Factor Lower bound from [DPW19]
def get_dpw_lower_bound(m):
    Re = RealField(W_PREC);
    K  = CyclotomicField(m);
    n  = euler_phi(m);      
    abs_disc = K.discriminant().abs();
    p  = K.next_split_prime();

    hf = Re(abs_disc)^(-1/Re(2*n)) * exp( ln(Re(p))/Re(n) * Re(0.02927)*Re(n)^(Re(3/2))
                                          + Re(0.1839)*sqrt(Re(n)) );

    dpw          = {};
    dpw["gf"]    = float(hf * sqrt(Re(2)*Re(pi)*Re(e)/Re(n)));
    dpw["afinf"] = float(hf / sqrt(Re(n)));
    dpw["afsup"] = float(hf / sqrt(Re(n)) * Re(abs_disc)^(1/Re(2*n)));
    dpw["hf"]    = float(hf^(1/Re(n)));
    assert(dpw.keys() == set(measures));

    return dpw;


# --------------------------------------------------------------------------------------------------
# To prepare the drawing data, we collect the minimal number of orbits giving the better average GF
# for each experiment, and compute the mean over all experiences (for a given conductor).
# [NB] In the paper's XP, this systematically translates to using d=1 for urs/sat and d=dmax for su.

for m in list_m:
    print("Preprocessing data for m={}".format(m));

    # Find dmax
    d_max = {};
    for _su in su_sets:
        d_max[_su] = find_max_nb_orbits(m, _su);
    print("\tMax nb orbits:\t{}".format(d_max));
        
    # Find optimal number of orbits for sets cdw/urs/sat/su
    d_opt = {};
    for _su in su_sets:
        d_opt[_su] = find_optimal_nb_orbits(m, _su, d_max[_su]);
    print("\tOptimal d:\t{}".format(d_opt));
    if (d_opt != { "cdw":1, "urs":1, "sat":1, "su":d_max["su"] }):
        print ("\t\x1b[31m[Warn] Unexpected optimal d_opt\x1b[0m");

    # DPW Lower bound (Hermite Factor)
    dpw = get_dpw_lower_bound(m);
    
    # Get Approximation factors, write resulting data
    for _mes in measures:
        afs     = [ file_get_mean_min(m, d_opt[_su], _su, _mes) for _su in su_sets ];
        af_write(m, _mes, d_opt, dpw, afs);


sys.exit(0);
