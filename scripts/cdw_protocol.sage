#!/usr/bin/env sage

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

if len(sys.argv) != 3:
    print(("Usage: {:s} <data_root> z<m>\n"+
           "\tfor Cyclotomic Field of conductor <m>\n"+
           "\tfor Factor Bases consisting of a maximum of <d> (split prime) orbits\n")
          .format(sys.argv[0]));
    sys.exit(2);
    
data_root = sys.argv[1];
tag   = sys.argv[2];
# This code a priori won't work as is for greater dmax
dmax  = 1; # ZZ(sys.argv[3]);


# --------------------------------------------------------------------------------------
# Obtain number field
K = nf_set_tag(tag);
r1, r2 = K.signature();
n = K.degree();
abs_disc = K.discriminant().abs();
print ("{}: eval approx factor of log S-unit lattices".format(tag), flush=True);


# ----------------------------------------------------------------------------------
# Reduction parameters
W_PREC    = 500;
NB_ITER   = n;
BLOCK_SZ  = 40;
MAX_LOOPS = 300;

# ----------------------------------------------------------------------------------
# For each precomputed sets, generate the corresponding lattice (with iso/noiso:exp/tw
opt_sets = [ "cdw" ];
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
        for _s, _iso, _inf in itertools.product(opt_sets, opt_iso.keys(), opt_inf.keys()):
            _f_out.write("{}/{}/{}\t".format(_s,_iso,_inf));

        _f_out.write("\n");
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



# for _d in range(dmax):
for _d in [dmax-1]:
    assert(_d == 0); 
    # --------------------------------------------------------------------------------
    # Read precomputations
    # Factor base
    f_fb  = data_dir + "{}_d{}.fb".format(tag, _d+1);
    fb    = fb_read_data(f_fb, K);
    # S-units for each set
    la_SU_all = [];
    _s   = "urs"; # CDW algorithm uses 'urs' family (when h+ = 1).
    f_su = data_dir + "{}_d{}.{}".format(tag, _d+1, _s);
    print ("Import raw S-units material for '{}' from '{}'".format(_s,f_su), end='', flush=True);
    t = cputime(); (yu, ysu), Bsu, Bvp = sunits_raw_read_data(f_su, K); t = cputime(t);
    print ("\t[done] t={:.2f}".format(t), flush=True);
        
    # -----------------------------------------------------------------------------
    # Obtain logarg representation of S-units.
    print ("Logarg(raw)\t\t", end='', flush=True);
    t = cputime(); la_Bsu = [ logarg_set(_g, p_inf, fb=fb, vp=_val_g) for _g, _val_g in zip(Bsu,Bvp) ]; t = cputime(t);
    print ("\t[done] t={:.2f}".format(t), flush=True);
    print ("su=mpow(logarg)\t\t", end='', flush=True);
    t = cputime(); la_su  = [ logarg_mpow(la_Bsu, _y) for _y in yu + ysu ]; t = cputime(t);
    print ("\t[done] t={:.2f}\n".format(t), flush=True);
    # -----------------------------------------------------------------------------
    la_SU_all.append(la_su);
    
    # ----------------------------------------------------------------------------------
    # Separate log args of UC / Real Gens / Sti Gens using output conventions of urs
    Uc_logarg    = la_SU_all[0][:r1+r2-1];             # /!\ WARN: Valid for URS only
    Rgens_logarg = la_SU_all[0][r1+r2-1:r1+r2-1+n//2]; # /!\ WARN: Valid for URS only
    Sgens_logarg = la_SU_all[0][r1+r2-1+n//2:];        # /!\ WARN: Valid for URS only

    # ----------------------------------------------------------------------------------
    # Lattices
    # Stickelberger lattice in S^-
    B_sti = Matrix([proj_Sminus_cdw(l.vp) for l in Sgens_logarg]);
    assert(B_sti.ncols() == B_sti.nrows() == n//2);
    # Sti_lat_proj, U_lll = bkz(Lsti, work_prec=W_PREC, block_size=BLOCK_SZ, bkzmaxloops=MAX_LOOPS);
    # Sgens_logarg_lll = [logarg_mpow(Sgens_logarg, _e) for _e in U_lll];
    u_sti = identity_matrix(ZZ, B_sti.nrows());
    B_sti, u_sti = bkz(B_sti, work_prec=W_PREC, block_size=BLOCK_SZ, bkzmaxloops=MAX_LOOPS);
    G_sti, _ = gram_schmidt_ortho(B_sti, normalize=False, b_prec=W_PREC);

    
    # Log unit lattice: in CDW, this is the only one impacted by the fHcE choices.
    # fHcE for all iso/noiso-exp/tw options
    print ("Compute fHcE matrices for Log U", end='', flush=True);
    t = cputime();
    fHcE_lu = { "{}/{}".format(_iso,_inf): get_twfHcE_matrix(r1, r2, 0, inf_type=opt_inf.get(_inf), isometry=opt_iso.get(_iso), b_prec=W_PREC)
                for _iso, _inf in itertools.product(opt_iso.keys(), opt_inf.keys()) };
    t = cputime(t);
    print ("\t[done] t={:.2f}\n".format(t), flush=True);
    
    # Log-embeddings for units
    B_logU = [];
    u_logU = [];
    G_logU = [];
    print ("Building log-circular unit lattice", flush=True);
    i = 0;
    for _s,_iso,_inf in itertools.product(opt_sets,opt_iso.keys(),opt_inf.keys()):
        print("\t{}:".format(l_names[i]), end='', flush=True);
        assert(l_names[i] == "{}/{}/{}".format(_s,_iso,_inf));

        # Log-unit lattice with correct embedding
        _logU = [ ];
        for _lu in Uc_logarg:
            _vu = logarg_log_embedding(_lu, p_inf, fb=[], inf_type=opt_inf.get(_inf));
            _vu = _vu*fHcE_lu.get("{}/{}".format(_iso,_inf));
            _logU.append(_vu);
        _logU = matrix(_logU);
        _u_logU = identity_matrix(ZZ, r1+r2-1);
        _logU, _u_logU = bkz(_logU, work_prec=W_PREC, block_size=BLOCK_SZ, bkzmaxloops=MAX_LOOPS);
        
        B_logU += [_logU];
        u_logU += [_u_logU];
        G_logU += [gram_schmidt_ortho(_logU, normalize=False, b_prec=W_PREC)[0]];
        i += 1;
    

    # --------------------------------------------------------------------------------
    # Results files
    out_names = [ data_dir + "{}_d{}_{}.{}".format(tag, _d+1, opt_sets[0], _measure) for _measure in measures_set ];
    out_files = [ open(_file, "w") for _file in out_names ];
    # Prepare Headers
    print_headers(out_files, _d);
  
    
    # --------------------------------------------------------------------------------
    # Read targets in logarg rep
    f_targets = data_dir + "{}_d{}.targets".format(tag, _d+1);    
    print("Read targets in '{}'".format(f_targets));
    targs     = logarg_read_data(f_targets, K);

     
    # --------------------------------------------------------------------------------
    # Now, the real deal
    for _k in range(len(targs)):
        print ("Challenge #{}:".format(_k));
        _t = targs[_k];
               
        # norm of b corresponding to target
        b_ln = logarg_lnSnorm_cf(_t, fb);
        Nb   = round(exp(b_ln));
        print("\tN(b):\t\t{}".format(Nb));
        
        count_ops = 0;
        sols = []; # indexed by measure set
        for _s,_iso,_inf in itertools.product(opt_sets, opt_iso.keys(),opt_inf.keys()):
            # Indices
            ind = opt_sets.index(_s);
            assert(l_names[count_ops] == "{}/{}/{}".format(_s, _iso, _inf) and opt_sets[ind] == _s);
            
            # Apply Tw-PHS for one target
            t = cputime();
            print("\tmethod:'{}':".format(l_names[count_ops]), flush=True);
            ls, ns = cdw_protocol(_t, p_inf, fb, fHcE_lu.get("{}/{}".format(_iso,_inf)),
                                  B_logU[count_ops], u_logU[count_ops], Uc_logarg,
                                  B_sti, u_sti, Rgens_logarg, Sgens_logarg,
                                  G_logU = G_logU[count_ops], G_sti = G_sti,
                                  inf_type=opt_inf.get(_inf), b_prec=W_PREC);
            t = cputime(t);
            print("\t[done] t2_norm={:7.3f} t={:.2f}".format(float(ns),t), flush=True);
            
            # Compute all ratios from sol
            sols      += [ls];
            count_ops += 1;
        
        af_rat = compute_afs(sols, Nb);
        assert(len(af_rat) == len(measures_set));
        assert(all(len(_af) == len(l_names) for _af in af_rat));
        af_write_data(out_files, af_rat);
        
exit;
