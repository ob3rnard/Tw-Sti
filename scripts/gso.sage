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

import itertools

from sage.all import *
import fp
from pcmp_io  import *
from ZR_mat   import *


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
W_PREC    = 500;


# ----------------------------------------------------------------------------------
# For each precomputed sets, generate the corresponding lattice (with iso/noiso:exp/tw)
opt_sets = ["urs"] + (["sat"] if b_sat == True else []) + (["su"] if b_su == True else []);
opt_iso  = { "iso": True,       "noiso": False };
opt_inf  = { "exp": "EXPANDED", "tw"   : "TWISTED" };


# ----------------------------------------------------------------------------------
# List of lattices
data_dir  = data_root + "/{}/".format(tag);


# --------------------------------------------------------------------------------------
# Going through all types (#orbs / sets / iso-exp )
#for _d,_s,_iso,_inf in itertools.product(range(dmax),opt_sets,opt_iso.keys(),opt_inf.keys()):
for _d,_s,_iso,_inf in itertools.product([dmax-1],opt_sets,opt_iso.keys(),opt_inf.keys()):
    print ("\n"+"-"*80+"\n#orbs={} typ='{}' [{}/{}]".format(_d+1,_s,_iso,_inf));

    # Lattices
    label = "d{}_{}_{}_{}".format(_d+1,_s,_iso,_inf);
    f_lat = data_dir + "{}_{}.lat".format(tag,label);
    f_lll = f_lat.replace(".lat",".lll");
    f_bkz = f_lat.replace(".lat",".bkz-{}".format(BLOCK_SZ));
    assert(all(os.path.exists(_f) for _f in [f_lat,f_lll,f_bkz]));
    
    # Read data
    print("Reading data to b_prec={}".format(W_PREC), flush=True);
    t = cputime(); B_raw = lattice_read_data(f_lat, to_b_prec=W_PREC); t = cputime(t);
    print("--> '{}' \t[done] t={:.2f}".format(f_lat, t), flush=True);
    t = cputime(); B_lll = lattice_read_data(f_lll, to_b_prec=W_PREC); t = cputime(t);
    print("--> '{}' \t[done] t={:.2f}".format(f_lll, t), flush=True);
    t = cputime(); B_bkz = lattice_read_data(f_bkz, to_b_prec=W_PREC); t = cputime(t);
    print("--> '{}' \t[done] t={:.2f}".format(f_bkz, t), flush=True);

    # Precompute GSO (This is really the hardest part with Sage)
    print("Compute raw/lll/bkz GSOs", end='', flush=True);
    t = cputime(); G_raw, _ = gram_schmidt_ortho(B_raw, normalize=False); t=cputime(t);
    print("\traw:t={:.2f}".format(t), end='', flush=True);
    t = cputime(); G_lll, _ = gram_schmidt_ortho(B_lll, normalize=False); t=cputime(t);
    print("\tlll:t={:.2f}".format(t), end='', flush=True);
    t = cputime(); G_bkz, _ = gram_schmidt_ortho(B_bkz, normalize=False); t=cputime(t);
    print("\tbkz:t={:.2f}".format(t), end='', flush=True);
    print("");
    
    # Output GSO's
    gso_lat = f_lat + "_gso";
    gso_lll = f_lll + "_gso";
    gso_bkz = f_bkz + "_gso";
    t = cputime(); lattice_write_data(gso_lat, G_raw); t = cputime(t);
    print("--> Output GSO raw lattice in '{}' t={:.2f}".format(gso_lat,t));
    t = cputime(); lattice_write_data(gso_lll, G_lll); t = cputime(t);
    print("--> Output GSO lll lattice in '{}' t={:.2f}".format(gso_lll,t));
    t = cputime(); lattice_write_data(gso_bkz, G_bkz); t = cputime(t);
    print("--> Output GSO bkz lattice in '{}' t={:.2f}".format(gso_bkz,t));

    
exit;

