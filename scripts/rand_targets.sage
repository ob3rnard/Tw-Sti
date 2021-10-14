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
from idsvp	 import *
import fp
from random_log_elt import *


if len(sys.argv) != 4:
    print(("Usage: {:s} <data_root> z<m> <d> sat=[true/false] su=[true/false]\n"+
           "\tfor Cyclotomic Field of conductor <m>\n"+
           "\tfor Factor Bases consisting of a maximum of <d> (split prime) orbits\n"+
	   "\tsat/su: precise whether we have saturated elements and/or full S-units <b>")
          .format(sys.argv[0]));
    sys.exit(2);

data_root = sys.argv[1];
tag   = sys.argv[2];
dmax  = ZZ(sys.argv[3]);



# ----------------------------------------------------------------------------------
# Parameters

# Precision of output (Reduced lattices have this prec)
W_PREC     = 1000;
# Number of targets needed
NB_TARGETS = 100;
# Bit size of alg. norm of virtual challenge b
CHAL_BSZ   = 100;


# -------------------------------------------------------------------------------------
# Obtain number field
K = nf_set_tag(tag);
m = cf_get_conductor(K);
n = K.degree();
assert (n == euler_phi(ZZ(tag[len("z"):])));


# ------------------------------------------------------------------------------------
# Input data files
data_dir        = data_root + "/{}/".format(tag);
inf_places_file = data_dir  + "{}.inf".format(tag);
fb_files        = [ data_dir + "{}_d{}.fb".format(tag, _d+1)  for _d in range(dmax) ];
target_files    = [ data_dir + "{}_d{}.{}".format(tag, _d+1, 'targets') for _d
                    in range(dmax) ];
assert(all(os.path.exists(_pcmp_file) for _pcmp_file in [inf_places_file] + fb_files));


# ----------------------------------------------------------------------------------
# Read Infinite places
phi = inf_places_read_data(inf_places_file, K);
phi = adapt_inf_places(K, phi, to_prec=W_PREC//2);


# ----------------------------------------------------------------------------------
# Loop on #orbits (fb)
for d in range(dmax):
#for d in [dmax-1]:
    print ("\n"+"-"*80 +"\nTargets for orbit=#{} \n".format(d+1)+"-"*80,
           flush=True);

    # Read Factor Base
    fb  = fb_read_data(fb_files[d], K);

    # Generate NB_TARGETS targets (logarg)
    t = cputime();
    targets = random_targets(K, phi, fb,
                             chal_bsz=CHAL_BSZ, b_prec=W_PREC, nb_targets=NB_TARGETS);
    t = cputime(t);
    print ("\t[done] t={:.2f}".format(t), flush=True);
    
    # Output data file
    logarg_file = target_files[d];
    print ("--> Output to '{}'".format(logarg_file), end='', flush=True);
    t = cputime(); logarg_write_data(logarg_file, K, targets); t = cputime(t);
    print ("\t[done] t={:.2f}\n".format(t), flush=True);
    
exit;
