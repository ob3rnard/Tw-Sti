#!/usr/bin/env sage

# Code in support of ePrint:2021/1384
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

if len(sys.argv) != 6:
    print(("Usage: {:s} <data_root> z<m> <d> sat=[true/false] su=[true/false]\n"+
           "\tfor Cyclotomic Field of conductor <m>\n"+
           "\tfor Factor Bases consisting of a maximum of <d> (split prime) orbits\n"+
	   "\tsat/su: precise whether we have saturated elements and/or full S-units<b>")
          .format(sys.argv[0]));
    sys.exit(2);

data_root = sys.argv[1];
tag   = sys.argv[2];
dmax  = ZZ(sys.argv[3]);
b_sat = True if sys.argv[4] == "sat=true" else False; 
b_su  = True if sys.argv[5] == "su=true"  else False; 


# ----------------------------------------------------------------------------------
# For each precomputed sets, generate the corresponding lattice (with iso/noiso:exp/tw)
opt_iso = { "iso": True,       "noiso": False };
opt_inf = { "exp": "EXPANDED", "tw"   : "TWISTED" };

# opt_iso = {"noiso": False};
# opt_inf = {"exp": "EXPANDED"};

# Precision of output.
HPREC=500;
 

# -------------------------------------------------------------------------------------
# Obtain number field
K = nf_set_tag(tag);
n = K.degree();
assert (n == euler_phi(ZZ(tag[len("z"):])))


# ------------------------------------------------------------------------------------
# Input data files
data_dir        = data_root + "/{}/".format(tag);
inf_places_file = data_dir  + "{}.inf".format(tag);
fb_files        = [data_dir + "{}_d{}.fb".format(tag, _d+1)  for _d in range(dmax) ];
sets            = ["urs"] + (["sat"] if b_sat == True else []) + (["su"] if b_su == True else []);
su_files        = [ [data_dir + "{}_d{}.{}".format(tag, _d+1, _set) for _set in sets] for _d in range(dmax) ];
assert(all(os.path.exists(_pcmp_file) for _pcmp_file in [inf_places_file] + fb_files + sum(su_files, [])));


# ----------------------------------------------------------------------------------
# Read Infinite places
phi = inf_places_read_data(inf_places_file, K);


# ----------------------------------------------------------------------------------
# Loop on #orbits (fb), then on urs/sat/su then on iso/inf parameters
for d in [dmax-1]:
    # ---------------------------------------
    # Read Factor Base
    fb  = fb_read_data(fb_files[d], K);

    # Looping on S-units sets for this orbit
    for su_pcmp in su_files[d]:
        print ("\n"+"-"*80 +"\nTreating orbit=#{} set='{}'\n".format(d+1, su_pcmp)+"-"*80, flush=True);
        print ("Import raw S-units material", end='', flush=True);
        t = cputime(); (yu, ysu), Bsu, Bvp = sunits_raw_read_data(su_pcmp, K); t = cputime(t);
        print ("\t[done] t={:.2f}".format(t), flush=True);
    
        # -----------------------------------------------------------------------------
        # Obtain logarg representation of S-units.
        print ("Logarg(raw)\t\t", end='', flush=True);
        _t = cputime(); la_Bsu = [ logarg_set(_g, phi, fb=fb, vp=_val_g) for _g, _val_g in zip(Bsu,Bvp) ]; _t = cputime(_t);
        print ("\t[done] t={:.2f}".format(_t), flush=True);
        print ("su=mpow(logarg)\t\t", end='', flush=True);
        _t = cputime(); la_su  = [ logarg_mpow(la_Bsu, _y) for _y in yu + ysu ]; _t = cputime(_t);
        print ("\t[done] t={:.2f}\n".format(_t), flush=True);
        # -----------------------------------------------------------------------------

        
        for iso, inf in itertools.product(opt_iso.keys(), opt_inf.keys()):
            print ("[{}/{}]\n".format(iso,inf) + "-"*60, flush=True);
            t = cputime();
            L = twphs_get_matrix_logarg(la_su, phi, fb,
                                        isometry=opt_iso.get(iso), inf_type=opt_inf.get(inf),
                                        b_prec=HPREC);
            t = cputime(t);
            print ("\tt={:.2f}".format(t), flush=True);
            
            # Output data file
            sulat_file = su_pcmp[::-1].replace(".","_",1)[::-1] + "_{}_{}.lat".format(iso, inf);
            print ("--> Output to '{}'".format(sulat_file), end='', flush=True);
            t = cputime(); lattice_write_data(sulat_file, L); t = cputime(t);
            print ("\t[done] t={:.2f}\n".format(t), flush=True);

exit;
