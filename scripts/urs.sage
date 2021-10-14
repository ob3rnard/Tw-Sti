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
from cyclotomics import * # Real Generators / Subfield
from circular    import * # Cyclotomic units
from stim_q      import * # Stickelberger generators
from pcmp_io     import * 


if len(sys.argv) != 4:
    print(("Usage: {:s} <data_root> z<m> <d>\n"+
           "\tfor Cyclotomic Field of conductor <m>\n"+
           "\tfor Factor Bases consisting of a maximum of <d> (split prime) orbits")
          .format(sys.argv[0]));
    sys.exit(2);

data_root = sys.argv[1];
tag  = sys.argv[2];
dmax = ZZ(sys.argv[3]);


# --------------------------------------------------------------------------------------
# Pari stack size for computing Cl_K+
pari.allocatemem(10^10); # Avoid "out of memory" when using Pari on the laptop


# --------------------------------------------------------------------------------------
# Obtain number field
K = nf_set_tag(tag);
m = cf_get_conductor(K);
n = K.degree();
assert(n == euler_phi(m));
print ("-"*100);
print ("{}: Get Stickelberger and Real generators".format(tag), flush=True);
# Class number
_t = cputime(); hK_pm = cf_hK(m); _t = cputime(_t);
assert(hK_pm[1] != None), "Unknown hK+"; 
assert(hK_pm[1] == 1),    "O_K+ is not principal"; # At this point, we don't allow non principal real subfields
hK    = hK_pm[0]*hK_pm[1];
print ("ClK:\t\t\thK={}\tt={:.2f}".format(hK, _t));


# --------------------------------------------------------------------------------------
# Output data folder
data_dir    = data_root + "/{}/".format(tag);
# Filenames
fb_files  = [ data_dir + "{}_d{}.fb".format(tag, _d+1)  for _d in range(dmax) ];
urs_files = [ data_dir + "{}_d{}.urs".format(tag, _d+1) for _d in range(dmax) ];


# --------------------------------------------------------------------------------------
# Factor base
# Obtain list of split primes
all_ps   = cf_d_first_split_primes(K, d=dmax);
all_orbs = cf_d_orbits(K, d=dmax);
# Read files and check consistency (orbits organization / integrity)
for _d in range(dmax):
    _fbi_d = fb_read_data(fb_files[_d], K);
    assert (all(pid_fast_gens_two(_pid_fb) == pid_fast_gens_two(_pid_all)
                for _pid_fb, _pid_all in zip(_fbi_d, all_orbs[:n*(_d+1)]))), "Bad factor bases";
    # Note that testing directly using Sage built-in '_fbi_d == all_orbs[:n*(_d+1)]' is intractable


# --------------------------------------------------------------------------------------
# Maximal Totally Real Subfield
_t = cputime(); Kp = real_maximal_subfield(K); _t = cputime(_t);
print ("Real Maximal Subfield:\tn+={}\tt={:.2f} K+={}".format(Kp.degree(), _t, str(list(Kp.defining_polynomial())).replace(' ','')), flush=True);


# Dealing with too small orbit of z492.
if m == 492:
    print ("/!\ Fine tuning of Pari");
    Kp_pari = pari(Kp);
    c1 = 0.062; # This is c >= max(norm(fb))/ln^2 (disc Kp)
    _t = cputime(); ClKp = Kp_pari.bnfinit(tech=[c1,c1,4]); hKp = 1; _t = cputime(_t);    
else:
    _t = cputime(); ClKp = Kp.class_group(); hKp = Kp.class_number(); _t = cputime(_t);
print ("Cl_K+:\t\t\thK+={}\tt={:.2f}".format(hKp, _t), flush=True);
assert (hKp == 1), "O_K+ is not principal"; # At this point, we don't allow non principal real subfields


# --------------------------------------------------------------------------------------
# Circular units
_t = cputime(); uc = cf_cyclotomic_units(K); _t = cputime(_t);
print ("Circular units:\t\tvu={}\tt={:.2f}".format(len(uc), _t), flush=True);


# --------------------------------------------------------------------------------------
# Real generators
r_vals = real_gens_valp_orbit(m); # Same for each orbit

if m == 492:
    r_gens = [];
    for _d in range(dmax):
        _p          = all_ps[_d];
        _real_orb_p = real_orbit_p(K, _p, M=Kp);
        _t  = cputime(); _fb_pari = [ pari(_pid) for _pid in _real_orb_p ]; _t = cputime(_t);
        print("Convert orb#{} to Pari:\tp={}\tt={:.2f} [{:.2f}/pid]".format(_d+1, _p, _t, float(_t/Kp.degree())), flush=True);
        _t  = cputime();
        _gens_p = [ Kp(ClKp.bnfisprincipal(_pid)[1]) for _pid in _fb_pari ];
        assert(all(_rid == Kp.ideal(_g0) for _rid,_g0 in zip(_real_orb_p,_gens_p)));
        r_gens   += [ [real_lift_elt(K, _g0) for _g0 in _gens_p] ];
        _t = cputime(_t);
        print ("Real Gens orb#{}:\tp={}\tt={:.2f} [{:.2f}/g]".format(_d+1, _p, _t, float(_t/Kp.degree())), flush=True);
else:
    r_gens = [];
    for _d in range(dmax):
        _p = all_ps[_d];
        _t = cputime(); r_gens += [real_gens(K, _p, M=Kp)]; _t = cputime(_t);
        print ("Real Gens orb#{}:\tp={}\tt={:.2f} [{:.2f}/g]".format(_d+1, _p, _t, float(_t/Kp.degree())), flush=True);


# --------------------------------------------------------------------------------------
# Stickelberger valuations / generators
# 1. Obtain generating set
# 2. Matrice vals
# 3. Compute gens for each orbit
print("Stickelberger valuations", flush=True);
t = cputime(); sw_vals = sti_vals_alpha(m, verbose=True).submatrix(1,0); t = cputime(t);
print("\tShort Kucera vals:\tt={:.2f}".format(t), flush=True);


# --------------------------------------------------------------------------------------
# Valuation for each orbit: check volume
rsw_orb_val = block_matrix([[r_vals],[sw_vals]]);
hK_sti  = rsw_orb_val.determinant().abs(); assert(hK.divides(hK_sti));
sat_idx = ZZ(hK_sti / hK);
t = len(m.factor());
b = (sat_idx == 2**(n//2-1 + (0 if t==1 else 2**(t-2)-1)));
print("det(r+sti) = {} = hK * \x1b[{}m{}\x1b[0m".format(hK_sti, "32" if b else "31", factor(sat_idx)), flush=True);


# --------------------------------------------------------------------------------------
# Stickelberger generators: omega(a)-omega(a+1) for a in S
sw_gens = [];
for _d in range(dmax):
    _p  = all_ps[_d];
    print("Sti gens for p={}".format(_p), flush=True);
    _p0 = all_orbs[_d*n];
    # _t = cputime(); sw_gens += [sti_gens_kucera_short_idx(K, _p0, S)]; _t = cputime(_t);
    _t = cputime(); sw_gens += [sti_gens_alpha(K, _p0)[1:]]; _t = cputime(_t); # Remove s[0] = p
    print("Sti Gens orb#{}:\tp={}\tt={:.2f} [{:.2f}/g]".format(_d+1, _p, _t, float(_t/(m//2))), flush=True);


# --------------------------------------------------------------------------------------
# Output for each orbit
assert (len(uc) == n//2-1);
assert (all( len(_r_orb)  == n//2 for _r_orb  in r_gens  ));
assert (all( len(_sw_orb) == n//2 for _sw_orb in sw_gens ));

for _d in range(dmax):
    B_su_d = uc + sum((r_gens[_i]+sw_gens[_i] for _i in range(_d+1)), []);
    B_vp_d = list(block_matrix(ZZ, [[zero_matrix(len(uc), (_d+1)*n)],
                                    [block_diagonal_matrix([rsw_orb_val]*(_d+1))]],
                               subdivide=False));
    y_all  = list(identity_matrix(ZZ, len(B_su_d)));
    y_u    = y_all[:len(uc)];
    y_su   = y_all[len(uc):];

    print("Orb#{} --> {}".format(_d+1, urs_files[_d]), flush=True);
    sunits_raw_write_data(urs_files[_d], K, y_u, y_su, B_su_d, B_vp_d);


print ("\x1b[32mDone\x1b[0m\n");
exit;

