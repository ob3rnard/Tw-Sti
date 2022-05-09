#!/usr/bin/env sage

# Code in support of ePrint:2021/1384
# Copyright 2021, Olivier Bernard, Andrea Lesavourey, Tuong-Huy Nguyen
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
from lattice  import *
from idsvp    import *
from nf       import *
from time import time


if len(sys.argv) != 6:
    print(("Usage: {:s} <data_root> <Ncpus> z<m> <d> set=[urs/sat/su]\n"+
           "\tfor Cyclotomic Field of conductor <m>\n"+
           "\tfor Factor Bases consisting of a maximum of <d> (split prime) orbits\n"+
	   "\turs/sat/su: precise which family of independent S-units to use")
          .format(sys.argv[0]));
    sys.exit(2);

data_root = sys.argv[1];
ncores    = ZZ(sys.argv[2]);
tag       = sys.argv[3];
orb       = ZZ(sys.argv[4]);
assert (sys.argv[5].startswith("set=")), "\x1b[31m[Err]\x1b[0m Wrong arguments.";
su_set    = sys.argv[5][len("set="):];
assert (su_set in ["urs", "sat", "su"]), "\x1b[31m[Err]\x1b[0m Wrong S-unit set.";


# --------------------------------------------------------------------------------------
# Obtain number field
K = nf_set_tag(tag);
r1, r2 = K.signature();
n = K.degree();
abs_disc = K.discriminant().abs();
print ("{}: eval approx factor of log S-unit lattices".format(tag), flush=True);


# ----------------------------------------------------------------------------------
# Reduction parameters
BLOCK_SZ  = 40;
W_PREC    = 500;
NB_ITER   = 100;


# ----------------------------------------------------------------------------------
# For each precomputed sets, generate the corresponding lattice (with iso/noiso:exp/tw)
opt_iso  = { "iso": True,       "noiso": False };
opt_inf  = { "exp": "EXPANDED", "tw"   : "TWISTED" };

measures_set = ["afsup", "gf", "afinf", "hf"];
l_names      = [ "{}/{}/{}".format(_s,_iso,_inf) for _s,_iso,_inf in itertools.product(opt_sets, opt_iso.keys(), opt_inf.keys()) ];


# ----------------------------------------------------------------------------------
# List of lattices
data_dir  = data_root + "/{}/".format(tag);


# Approx factor
# ------------------------------
# Each line is: approx factor by simu. target t wrt different options for lattices
def print_headers(streams, d):
    assert(len(streams) == len(measures_set));
    for _i in range(len(measures_set)):
        _f_out = streams[_i];
        _f_out.write("# Measure:'{}' nf:'{}' orb={} prec={} bkz_sz={}\n# ".format(measures_set[_i], tag, d, W_PREC, BLOCK_SZ));
        _f_out.write(("{}\t"*(len(l_names)-1)+"{}\n").format(*l_names));
        _f_out.flush();
    return;


def approx_factor_sup(ls, Nb):
    Re     = RealField(W_PREC);
    b_inf  = sqrt(Re(n)) * Re(Nb)^(1/Re(n));
    t2_s   = logarg_t2_norm_cf(ls);
    return Re(t2_s/b_inf);


def approx_factor_inf(ls, Nb):
    Re     = RealField(W_PREC);
    b_sup  = sqrt(Re(n)) * Re(Nb)^(1/Re(n)) * abs_disc^(1/Re(2*n));
    t2_s   = logarg_t2_norm_cf(ls);
    return Re(t2_s/b_sup);


# compute root hermite factor
def hermite_factor(ls, Nb):
    Re     = RealField(W_PREC);
    t2_s   = logarg_t2_norm_cf(ls);
    hf     = ( t2_s / Re(Nb)^(1/Re(n)) / abs_disc^(1/Re(2*n)) )^(1/Re(n));
    return hf;


# compute norm of sol. / gaussian heuristic
def gaussian_factor(ls, Nb):
    Re     = RealField(W_PREC);
    gauss  = sqrt(Re(n)/Re(2)/Re(pi)/Re(e)) * Re(Nb)^(1/Re(n)) * abs_disc^(1/Re(2*n));
    t2_s   = logarg_t2_norm_cf(ls);
    return t2_s / gauss ;


def compute_afs(sols, Nb):
    # measures_set = ["afsup", "gf", "afinf", "hf"];    
    af_sup = [ approx_factor_sup(_sol, Nb) for _sol in sols ];
    gf     = [ gaussian_factor  (_sol, Nb) for _sol in sols ];
    af_inf = [ approx_factor_inf(_sol, Nb) for _sol in sols ];
    hf     = [ hermite_factor(_sol, Nb)    for _sol in sols ];
    return [af_sup, gf, af_inf, hf];

def af_write_data(streams, data):
    assert(len(data) == len(measures_set) and all(len(_af) == len(l_names) for _af in data));

    for _f, _data in zip(streams, data):
        _f.write(("{:<13.4f}\t"*(len(l_names)-1)+"{:<13.4f}\n").format(*(map(float, _data))));
        _f.flush();
    return;


# -------------------------------------------------------------------------------------
# Read infinte places
inf_places_file = data_dir + "{}.inf".format(tag);
p_inf = inf_places_read_data(inf_places_file, K);
p_inf = adapt_inf_places(K, p_inf, to_prec=W_PREC);


# --------------------------------------------------------------------------------
# Read precomputations
# Factor base
f_fb  = data_dir + "{}_d{}.fb".format(tag, orb);
fb    = fb_read_data(f_fb, K);
# S-units for each set
f_su = data_dir + "{}_d{}.{}".format(tag, orb, su_set);
print ("Import raw S-units material for '{}' from '{}'".format(su_set,f_su), end='', flush=True);
t = time(); (yu, ysu), Bsu, Bvp = sunits_raw_read_data(f_su, K); 
print ("\t[done] t={:.2f}".format(time()-t), flush=True);

# -----------------------------------------------------------------------------
# Obtain logarg representation of S-units.
print ("Logarg(raw)\t\t", end='', flush=True);
t = time(); la_Bsu = [ logarg_set(_g, p_inf, fb=fb, vp=_val_g) for _g, _val_g in zip(Bsu,Bvp) ]; 
print ("\t[done] t={:.2f}".format(time()-t), flush=True);
print ("su=mpow(logarg)\t\t", end='', flush=True);
t = time(); la_Su  = [ logarg_mpow(la_Bsu, _y) for _y in yu + ysu ]; 
print ("\t[done] t={:.2f}\n".format(time()-t), flush=True);
# // Su Done ------------------------------------------------------------------

# fHcE for all iso/noiso-exp/tw options
print ("Compute fHcE matrices", end='', flush=True);
t = time();
fHcE = { "{}/{}".format(_iso,_inf): get_twfHcE_matrix(r1, r2, len(fb), inf_type=opt_inf.get(_inf), isometry=opt_iso.get(_iso), b_prec=W_PREC)
         for _iso, _inf in itertools.product(opt_iso.keys(), opt_inf.keys()) };
t = time()-t
print ("\t[done] t={:.2f}\n".format(t), flush=True);

# Read lattices: BKZ+GSO+lu
B_bkz   = [];
U_bkz   = [];
G_bkz   = [];
print ("Reading bases BKZ/GSO/LU", flush=True);
i = 0;
for _iso,_inf in itertools.product(opt_iso.keys(),opt_inf.keys()):
    print("\t{}:".format(l_names[i]), end='', flush=True);
    assert(l_names[i] == "{}/{}/{}".format(su_set, _iso, _inf));
    
    # Lattices
    label     = "d{}_{}_{}_{}".format(orb, su_set, _iso, _inf);
    f_bkz     = data_dir + "{}_{}.bkz-{}".format(tag, label, BLOCK_SZ);
    f_bkz_gso = data_dir + "{}_{}.bkz-{}_gso".format(tag, label, BLOCK_SZ);
    f_bkz_u   = data_dir + "{}_{}.bkz-{}_U".format(tag, label, BLOCK_SZ);
    assert(all(os.path.exists(_file) for _file in [f_bkz,f_bkz_u,f_bkz_gso]));
        
    # Read lattice and transformation matrix
    t = time(); B_bkz += [lattice_read_data(f_bkz,     to_b_prec=W_PREC)]; t = time()-t;
    print("\tBKZ:'{}'\tt={:.2f}".format(f_bkz, t),     end='', flush=True);
    t = time(); G_bkz += [lattice_read_data(f_bkz_gso, to_b_prec=W_PREC)]; t = time()-t;
    print("\tGSO:'{}'\tt={:.2f}".format(f_bkz_gso, t), end='', flush=True);
    t = time(); U_bkz += [lattice_ZZ_read_data(f_bkz_u)]; t = time()-t;
    print("\tLU:'{}'\tt={:.2f}".format(f_bkz_u, t), flush=True);
    i += 1;
    
    
# --------------------------------------------------------------------------------
# Results files
out_names = [ data_dir + "{}_d{}_{}.{}".format(tag, orb, su_set, _measure)
              for _measure in measures_set ]; 
out_files = [ open(_file, "w") for _file in out_names ];
# Prepare Headers
print_headers(out_files, orb-1);

    
# --------------------------------------------------------------------------------
# Read targets in logarg rep
f_targets = data_dir + "{}_d{}.targets".format(tag, orb);    
print("Read targets in '{}'".format(f_targets));
targets   = logarg_read_data(f_targets, K);

     
# --------------------------------------------------------------------------------
# Now, the real deal
t_all = time()
for _k in range(len(targets)):
    print ("\n"+"-"*60+"\nChallenge #{}:".format(_k));
    _target = targets[_k];
    
    # norm of b corresponding to target
    b_ln = logarg_lnSnorm_cf(_target, fb);
    Nb   = round(exp(b_ln));
    print("\tN(b):\t\t{}".format(Nb));
    
    count_ops = 0;
    sols  = []; # indexed by measure set
    t_tar = time()
    for _iso,_inf in itertools.product(opt_iso.keys(),opt_inf.keys()):
        # Indices
        assert(l_names[count_ops] == "{}/{}/{}".format(su_set, _iso, _inf));
        
        # Apply Tw-PHS for one target
        print("\tmethod:'{}':".format(l_names[count_ops]), flush=True);
        t = time();
        ls, ns = twphs_random(_target, p_inf, fb, fHcE.get("{}/{}".format(_iso,_inf)),
                              B_bkz[count_ops], U_bkz[count_ops], la_Su,
                              inf_type=opt_inf.get(_inf), set_tag=su_set, G=G_bkz[count_ops],
                              nb_cores=ncores, b_prec=W_PREC);
        t = time()-t;
        print("\t[done] t2_norm={:7.3f} t={:.2f}".format(float(ns),t), flush=True);
        
        # Compute all ratios from sol
        sols      += [ls];
        count_ops += 1;
        
    af_rat = compute_afs(sols, Nb);
    assert(len(af_rat) == len(measures_set));
    assert(all(len(_af) == len(l_names) for _af in af_rat));
    t_tar  = time() - t_tar;
    af_write_data(out_files, af_rat);
    
    print("%% Challenge #{} [done] in t={:.2f} %%".format(_k, t_tar), flush=True);
    
t_all = time() - t_all;
print("\n[ALL DONE] in t={:.2f}".format(t_all));
        
exit;
