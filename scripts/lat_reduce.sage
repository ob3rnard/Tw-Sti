#!/usr/bin/env sage

# Code in support of ePrint:2021/xxxx
# Copyright 2021, Olivier Bernard, Andrea Lesavourey
# GPL-3.0-only (see LICENSE file)

# Small convenient hack to simulate as if we were always executing from '<trunk>/'
import os
import sys
os.chdir(os.path.dirname(__file__) + "/../");
sys.path.append("./src/");
# --------------------------------------------------------------------------------------

import itertools

from sage.all    import *
from pcmp_io     import * 
from lattice	 import *

if len(sys.argv) != 6:
    print(("Usage: {:s} <data_root> z<m> <d> sat=[true/false] su=[true/false]\n"+
           "\tfor Cyclotomic Field of conductor <m>\n"+
           "\tfor Factor Bases consisting of a maximum of <d> (split prime) orbits\n"+
	   "\tsat/su: precise whether we have saturated elements and/or full S-units <b>")
          .format(sys.argv[0]));
    sys.exit(2);

data_root = sys.argv[1];
tag   = sys.argv[2];
dmax  = ZZ(sys.argv[3]);
b_sat = True if sys.argv[4] == "sat=true" else False; 
b_su  = True if sys.argv[5] == "su=true"  else False; 


# ----------------------------------------------------------------------------------
# Reduction parameters
BLOCK_SZ  = 40;
MAX_LOOPS = 300;
W_PREC    = 500; # fplll performance drops drastically beyond 500

opt_red  = { "lll": lll, "bkz-{}".format(BLOCK_SZ): bkz };


# ----------------------------------------------------------------------------------
# For each precomputed sets, generate the corresponding lattice (with iso/noiso:exp/tw)
opt_iso = { "iso": True,       "noiso": False };
opt_inf = { "exp": "EXPANDED", "tw"   : "TWISTED" };


# ----------------------------------------------------------------------------------
# List of lattices
data_dir  = data_root + "/{}/".format(tag);
sets      = ["urs"] + (["sat"] if b_sat == True else []) + (["su"] if b_su == True else []);
#list_lat  = [ data_dir + "{}_d{}_{}_{}_{}.lat".format(tag,_d+1,_s,_iso,_inf) for
#              _d,_s,_iso,_inf in itertools.product(range(dmax), sets, opt_iso.keys(), opt_inf.keys()) ];
list_lat  = [ data_dir + "{}_d{}_{}_{}_{}.lat".format(tag,_d+1,_s,_iso,_inf) for
              _d,_s,_iso,_inf in itertools.product([dmax-1], sets, opt_iso.keys(), opt_inf.keys()) ];
out_lll   = [ _f.replace(".lat",".lll") for _f in list_lat ];
out_lll   = [ (_f, _f + "_U")           for _f in out_lll ];
out_bkz   = [ _f.replace(".lat",".bkz-{}".format(BLOCK_SZ)) for _f in list_lat ];
out_bkz   = [ (_f, _f + "_U")                               for _f in out_bkz];
assert (all(os.path.exists(_f) for _f in list_lat));


# ----------------------------------------------------------------------------------
# For each lattice, compute LLL / BKZ red and output
for _f_lat, _f_lll, _f_bkz in zip(list_lat, out_lll, out_bkz):
    print("Lat '{}'".format(_f_lat), flush=True);

    # Input lattice
    t = cputime(); B = lattice_read_data(_f_lat); t = cputime(t);
    print("\tReading t={:.2f}".format(t));
    
    # Launching lll
    print("\tlll prec={} ...\t\t".format(W_PREC), end='', flush=True);
    t = walltime();
    B_lll, U_lll = lll(B, work_prec=W_PREC); # B_lll=U_lll*B
    t = walltime(t);
    print("\t[done] t={:.2f}".format(t), flush=True);
    
    # Launching bkz-k
    print("\tbkz prec={} bk={} nloops={} ...".format(W_PREC, BLOCK_SZ, MAX_LOOPS), end='', flush=True);
    t = walltime();
    B_bkz, U_bkz = bkz(B, work_prec=W_PREC, block_size=BLOCK_SZ, bkzmaxloops=MAX_LOOPS); # B_bkz=U_bkz*B
    t = walltime(t);
    print("\t[done] t={:.2f}".format(t), flush=True);

    # Output
    # Output also the U_lll, U_bkz, this will help for CVP in the form y.B_lll = (y*U_lll).B
    lattice_write_data(_f_lll[0], B_lll);
    lattice_write_data(_f_lll[1], U_lll);    
    print("\t--> output in '{}', '{}'".format(*_f_lll));
    lattice_write_data(_f_bkz[0], B_bkz);
    lattice_write_data(_f_bkz[1], U_bkz);    
    print("\t--> output in '{}', '{}'".format(*_f_bkz));
    

exit;
