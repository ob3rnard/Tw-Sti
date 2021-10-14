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

from sage.all import *
import fp
from pcmp_io  import *
from ZR_mat   import *
from lattice  import *
from idsvp import get_proj_HE, get_twfHcE_matrix

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


# --------------------------------------------------------------------------------------
# Obtain number field
K = nf_set_tag(tag);
r1, r2 = K.signature();
n = K.degree();
print ("{}: eval geometry of log S-unit lattices".format(tag), flush=True);


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


# ----------------------------------------------------------------------------------
# Output file
out_geo = data_dir + "{}.geo".format(tag);


# --------------------------------------------------------------------------------------
# Output stuff

# GSN File: (one per label)
# ------------------------------
# Each file contains 4 columns: i ln(ai*) ln(bi*) ln(ci*)
# ai* = GSO(raw), bi*=GSO(lll) ci*=GSO(bkz)
def __print_gsn_head(f_out, label):
    f_out.write("# nf:'{}' label:'{}'\n".format(tag, label));
    f_out.write("# i\tln(raw*)\tln(lll*)\tln(bkz*)\n");
    f_out.flush();

def __print_gsn_row(f_out, k, data):    
    f_out.write(("{:d}"+"\t{:14.12f}"*3+"\n").format(k, *[float(_d) for _d in data]));
    f_out.flush();

def gsnorms_write_data(filename, label, G_raw, G_lll, G_bkz):
    assert (G_raw.nrows() == G_lll.nrows() == G_bkz.nrows());
    _f_out = open(filename, "w");

    __print_gsn_head(_f_out, label);
    _dim   = G_raw.nrows();
    for _k in range(_dim):
        _gsn_k = [ ln(G_raw[_k].norm()), ln(G_lll[_k].norm()), ln(G_bkz[_k].norm()) ];
        __print_gsn_row(_f_out, _k, _gsn_k);    
    
    _f_out.close();
    return;


# GEO File:
# --------------------------------
# Each line is: label dim rvol h_raw h_lll h_bkz o_raw o_lll o_bkz
def print_geo_head(f_out):
    f_out.write("# nf:'{}'\n".format(tag));
    f_out.write("# label\t\t\tdim\tr_vol\th_raw\th_lll\th_bkz\to_raw\to_lll\to_bkz\n");
    f_out.flush();
    return;
    
def print_geo_line(f_out, label, data):
    f_out.write("{:16s}".format(label) + "\t{:<4d}".format(data[0])
                + ("\t{:<7.3f}"*7).format(*[float(_d) for _d in data[1:]])
                + "\n");
    f_out.flush();
    return;


# --------------------------------------------------------------------------------------
# Going through all types (#orbs / sets / iso-exp )
f_geo = open(out_geo, "w");
print_geo_head(f_geo);
for _d,_s,_iso,_inf in itertools.product(range(dmax),opt_sets,opt_iso.keys(),opt_inf.keys()):
    print ("\n"+"-"*80+"\n#orbs={} typ='{}' [{}/{}]".format(_d+1,_s,_iso,_inf));

    # Lattices
    label = "d{}_{}_{}_{}".format(_d+1,_s,_iso,_inf);
    f_lat = data_dir + "{}_{}.lat".format(tag,label);
    f_lll = f_lat.replace(".lat",".lll");
    f_bkz = f_lat.replace(".lat",".bkz-{}".format(BLOCK_SZ));
    gso_lat = f_lat + "_gso";
    gso_lll = f_lll + "_gso";
    gso_bkz = f_bkz + "_gso";
    assert(all(os.path.exists(_f) for _f in [f_lat,f_lll,f_bkz,gso_lat,gso_lll,gso_bkz]));
    # GSN output
    out_gsn = data_dir + "{}_{}.gsn".format(tag,label);
    
    # Read data
    print("Reading data to b_prec={}".format(W_PREC), flush=True);
    t = cputime(); B_raw = lattice_read_data(f_lat, to_b_prec=W_PREC); t = cputime(t);
    print("--> '{}' \t[done] t={:.2f}".format(f_lat, t), flush=True);
    t = cputime(); B_lll = lattice_read_data(f_lll, to_b_prec=W_PREC); t = cputime(t);
    print("--> '{}' \t[done] t={:.2f}".format(f_lll, t), flush=True);
    t = cputime(); B_bkz = lattice_read_data(f_bkz, to_b_prec=W_PREC); t = cputime(t);
    print("--> '{}' \t[done] t={:.2f}".format(f_bkz, t), flush=True);

    # Computing GSO is really the hardest part with Sage --> Pcmp
    print("Read precomputed raw/lll/bkz GSOs", end='', flush=True);
    t = cputime(); G_raw = lattice_read_data(gso_lat, to_b_prec=W_PREC); t=cputime(t);
    print("\traw:t={:.2f}".format(t), end='', flush=True);
    t = cputime(); G_lll = lattice_read_data(gso_lll, to_b_prec=W_PREC); t=cputime(t);
    print("\tlll:t={:.2f}".format(t), end='', flush=True);
    t = cputime(); G_bkz = lattice_read_data(gso_bkz, to_b_prec=W_PREC); t=cputime(t);
    print("\tbkz:t={:.2f}".format(t), end='', flush=True);
    print("");
    
    # Output GS norms
    print("--> compute / output gs-norms to '{}'".format(out_gsn), end='', flush=True);
    t = cputime(); gsnorms_write_data(out_gsn, label, G_raw, G_lll, G_bkz); t = cputime(t);
    print("\t[done] t={:.2}".format(t), flush=True);
    
    # Dimension
    dim = G_raw.nrows();
    
    # Reduced volume
    print("Reduced volume", end='', flush=True);
    t = cputime(); r_vol = vol_reduced(B_bkz, gso=G_bkz); t = cputime(t);
    print("\t[done] t={:.2f}".format(t));
    assert(fp.fp_check_zero("vol(raw)=vol(bkz)", [r_vol - vol_reduced(B_raw,gso=G_raw)], target=W_PREC, sloppy=True));
    assert(fp.fp_check_zero("vol(raw)=vol(lll)", [r_vol - vol_reduced(B_lll,gso=G_lll)], target=W_PREC, sloppy=True));
    ln_vol = B_bkz.base_ring()(dim)*ln(r_vol);
    
    # Hermite before/after BKZ
    print("Root Hermite factors raw/lll/bkz", end='', flush=True);
    t = cputime(); d0_raw = rankin_reduced(B_raw, 1, ln_V=ln_vol); t = cputime(t);
    print("\traw:t={:.2f}".format(t), end='', flush=True);
    t = cputime(); d0_lll = rankin_reduced(B_lll, 1, ln_V=ln_vol); t=cputime(t);
    print("\tlll:t={:.2f}".format(t), end='', flush=True);
    t = cputime(); d0_bkz = rankin_reduced(B_bkz, 1, ln_V=ln_vol); t=cputime(t);
    print("\tbkz:t={:.2f}".format(t), end='', flush=True);
    print("");
    
    # Ortho defect before/after BKZ
    print("Orthogonality defect raw/lll/bkz", end='', flush=True);
    t = cputime(); dn_raw = rankin_reduced(B_raw, B_raw.nrows(), ln_V=ln_vol); t = cputime(t);
    print("\traw:t={:.2f}".format(t), end='', flush=True);
    t = cputime(); dn_lll = rankin_reduced(B_lll, B_lll.nrows(), ln_V=ln_vol); t = cputime(t);
    print("\tlll:t={:.2f}".format(t), end='', flush=True);
    t = cputime(); dn_bkz = rankin_reduced(B_bkz, B_bkz.nrows(), ln_V=ln_vol); t = cputime(t);
    print("\tbkz:t={:.2f}".format(t), end='', flush=True);
    print("");
    
    # Output
    geo_line = [dim, r_vol,
                d0_raw,  d0_lll, d0_bkz, dn_raw, dn_lll, dn_bkz];
    print_geo_line(f_geo, label, geo_line);


f_geo.close();
exit;

