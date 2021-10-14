#!/usr/bin/env sage

# Code in support of ePrint:2021/xxxx
# Copyright 2021, Olivier Bernard
# GPL-3.0-only (see LICENSE file)

# Small convenient hack to simulate as if we were always executing from '<trunk>/'
import os
import sys
os.chdir(os.path.dirname(__file__) + "/../");
sys.path.append("./src/");
# --------------------------------------------------------------------------------------

from sage.all    import *
from nf          import *
from cyclotomics import *
from idsvp       import *
from pcmp_io     import *


if len(sys.argv) != 3:
    print("Usage: {:s} <data_root> z<m> for Cyclotomic Field of conductor <m>\n".format(sys.argv[0]));
    sys.exit(2);

data_root = sys.argv[1];
tag = sys.argv[2];


# --------------------------------------------------------------------------------------
# Obtain number field
K = nf_set_tag(tag);
print ("{}: get places".format(tag), flush=True);


# --------------------------------------------------------------------------------------
# Output data folder
data_dir = data_root + "/{}/".format(tag);


# --------------------------------------------------------------------------------------
# Compute inf places at high precision
# Only the order really matters, this is very fast.
HPREC = 5000;
p_inf_name  = data_dir + "{}.inf".format(tag);

_t = cputime(); p_inf = get_inf_places(K, b_prec=HPREC); _t = cputime(_t);
print ("Inf_places -> {} prec={} t={:.2f}".format(p_inf_name, HPREC, _t), end='', flush=True);
inf_places_write_data(p_inf_name, K, p_inf);
print ("\t[done]", flush=True);


# --------------------------------------------------------------------------------------
# Compute factor bases for tw and all d-orbits from d=1 to d(tw)

# Evaluate hR
_t = cputime(); hR = analytic_hR(K); _t = cputime(_t);
print ("Analytic hR: hR=exp({:.2f}) t={:.2f}".format(float(ln(hR)), _t), flush=True);

# Find nb orbits for tw-PHS
_t = cputime(); d_tw = tw_get_dimfb_norb(K, hR=hR); _t = cputime(_t);
print ("tw-PHS orbits: d={} t={:.2f}".format(d_tw, _t), flush=True);

# Generate fb files for all orbits
for _d in range(1,d_tw+1):
    fb_d_orb_name = data_dir + "{}_d{}.fb".format(tag,_d);
    _t = cputime(); fb_d_orb = cf_d_orbits(K, d=_d); _t = cputime(_t);
    print ("Orbs #{}->{} k={} t={:.2f}".format(_d, fb_d_orb_name, len(fb_d_orb), _t),
           end='', flush=True);
    fb_write_data(fb_d_orb_name, K, fb_d_orb);
    print ("\t[done]", flush=True);    
    
exit;
