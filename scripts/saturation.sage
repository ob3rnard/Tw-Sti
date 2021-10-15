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

from sage.all import *
from nf          import *
from pideal      import *
from cyclotomics import *
from circular    import *
from lattice     import * # CVP Babai NP
from idsvp       import * # For reduction in log_U
from pcmp_io     import *

from nfroot      import *


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
# Quadratic character mod pid
# Returns the next prime q > p such that q = 1 mod div_pm1
def __next_suitable_prime(p, div_pm1):
    # Start from _p0 > p st. _p0 == 1 [md]
    p0 = p + div_pm1 - p.mod(div_pm1) + 1;
    while not is_prime(p0):
        p0 += div_pm1;
    return p0;

def log_chip_quadratic(pid, a):
    # O_K / pid O_K = F_q
    (p, g) = pid_fast_gens_two(pid); f = g.degree();
    q      = ZZ(p**f); Fq = GF(q);
    # a mod pid: using 'pid.small_residue(a)' is several times slower
    Fqy = PolynomialRing(Fq, name='y'); # y = Fqy.gen();
    om  = Fq( Fqy(a.polynomial()).mod(Fqy(g)) );
    # If om \in pid, return +Infinity (character must be discarded)
    if (om == 0):
       print ("/!\\ Warning: {} = 0 mod {}, must be discarded".format(a, pid_fast_gens_two()));
       return +Infinity;

    # Quadratic character mod pid (log_1)
    log_chip = 0 if om.is_square() else 1;

    return log_chip;

# --------------------------------------------------------------------------------------
def reduce_mod_logU(G, log_G, U, log_U, log_U_gso):
    Y_cvp = [ cvp_babai_NP(log_U, _log_g, G=log_U_gso)[1] for _log_g in log_G ];
    S     = [ prod(map(pow, U, _y)) for _y in Y_cvp ];
    G_red = [ _g/_s for _g, _s in zip(G, S) ];
    return G_red;



# --------------------------------------------------------------------------------------
# Working precision
# NB: the lattice reduction will use prec = 500, so...
W_PREC       = 500;
# Supplementary characters
CHI_OVERHEAD = 20;
# Transformation matrix
BKZ_BK       = 40;


# --------------------------------------------------------------------------------------
# Working dir
data_dir  = data_root + "/{}/".format(tag);
inf_file  = data_dir + "{}.inf".format(tag);
fb_files  = [ data_dir + "{}_d{}.fb".format(tag, _d+1)  for _d in range(dmax) ];
urs_files = [ data_dir + "{}_d{}.urs".format(tag, _d+1) for _d in range(dmax) ]; # WARN: only works on .urs files at the moment.
sat_files = [ data_dir + "{}_d{}.sat".format(tag, _d+1) for _d in range(dmax) ];
assert(all(os.path.exists(_f) for _f in [inf_file] + fb_files + urs_files));


# --------------------------------------------------------------------------------------
# Obtain number field
K = nf_set_tag(tag);
m = cf_get_conductor(K);
n = K.degree();
assert(n == euler_phi(m));
r1,r2 = K.signature();
assert((r1 == 0) and (r2 == n//2));


# --------------------------------------------------------------------------------------
# Read places
p_inf = inf_places_read_data(inf_file, K);
p_inf = adapt_inf_places(K, p_inf, to_prec=W_PREC);


# Consistency of factor base (sorted as fb[i] = sigma_i(fb[0])), sigma_i : z -> z^i
all_ps   = cf_d_first_split_primes(K, d=dmax);
all_orbs = cf_d_orbits(K, d=dmax);
# Read files and check consistency (orbits organization / integrity)
fb       = [ all_orbs[n*_d:n*(_d+1)] for _d in range(dmax) ];
for _d in range(dmax):
    _fbi_d = fb_read_data(fb_files[_d], K);
    assert (all(pid_fast_gens_two(_pid_fb) == pid_fast_gens_two(_pid_all)
                for _pid_fb, _pid_all in zip(_fbi_d, all_orbs[:n*(_d+1)]))), "Bad factor bases";
    # Note that testing directly using Sage built-in '_fbi_d == all_orbs[:n*(_d+1)]' is intractable


# --------------------------------------------------------------------------------------
# Units
mu = cf_get_torsion_units(K); # [-z]
uc = cf_cyclotomic_units(K);


# --------------------------------------------------------------------------------------
# Read S-Units
print ("Import all raw S-units material\t\t\t", end='', flush=True);
t = cputime();
(yu_all, ysu_all), B_all, Vp_all = sunits_raw_read_data(urs_files[-1], K);
t = cputime(t);
print ("\t[done] t={:.2f}".format(t), flush=True);

print ("Check consistency\t\t\t\t\t", end='', flush=True);
Su = [ B_all[len(uc)+_d*n:len(uc)+(_d+1)*n] for _d in range(dmax) ];
Vp = [ matrix(Vp_all)[len(uc)+_d*n:len(uc)+(_d+1)*n,
                        _d*n:(_d+1)*n] for _d in range(dmax) ]; # Well they should be all equal...
for _d in range(dmax):
    t = cputime();
    (_yu_d, _ysu_d), _B_d, _Vp_d = sunits_raw_read_data(urs_files[_d], K);
    assert(matrix(_yu_d+_ysu_d) == identity_matrix(len(uc)+(_d+1)*n)); # urs files
    assert(_B_d  == uc + sum(Su[:_d+1], []));
    assert(matrix(_Vp_d[len(uc):]) == block_diagonal_matrix(Vp[:_d+1]));
    t = cputime(t);
    print("\x1b[32m[OK]\x1b[0m orb:#{}/t={:.2f} ".format(_d+1,t), end='', flush=True);
print("");


# --------------------------------------------------------------------------------------
# Obtain Log-unit^2 lattice (for reduction)
u2 = [ _u^2 for _u in uc ];
la_u2  = [ logarg_set(_u2, p_inf) for _u2 in u2 ];
log_U2 = matrix([ logarg_log_embedding(_la_u2, p_inf=p_inf, inf_type='TWISTED')
                  for _la_u2 in la_u2 ]);
fH     = get_twfHcE_matrix(r1, r2, 0, inf_type='TWISTED', b_prec=W_PREC);
log_U2 = log_U2*fH;
log_U2_gso, _ = gram_schmidt_ortho(log_U2, normalize=False);


# --------------------------------------------------------------------------------------
# HERE IS LOOP FOR _D in range(dmax)
new_Su  = [0]*dmax;
new_Vp  = [0]*dmax;
new_Ysu = [0]*dmax;
# h_K . 2^{\varphi(m)/2 - 1} . 2^b, with b = 0 if t = 1, 2^{t-2}-1 if t > 1.
exp_rank = (n//2 - 1 + (0 if is_prime_power(m) else (2^(len(factor(m))-2)-1) ) );
for _d in range(dmax):
    print ("\n"+"-"*80 +"\nTreating orbit=#{} set='{}'\n".format(_d+1, urs_files[_d])+"-"*80, flush=True);

    # ----------------------------------------------------
    # Read Factor Base (actually we really only need only norm(fb[-1]))
    fb_d = fb[_d];

    # ----------------------------------------------------
    # Extended (ie. with torsion) set of S-units / valuations
    ext_B = mu + uc + Su[_d];
    Vsu   = Vp[_d];

    # ----------------------------------------------------
    # Set of characters
    # Choose sufficiently many quadratic characters
    div_pm1  = lcm(m, 2);
    l_pZ     = [__next_suitable_prime(pid_fast_smallest_integer(fb_d[-1]), div_pm1)];
    nb_chi_p = len(ext_B) - len(fb_d) + CHI_OVERHEAD; # Overkill but why not
    print ("Choosing N={} chi_p from p >= {}\t\t".format(nb_chi_p, l_pZ[0]), end='', flush=True);
    t = cputime();
    while (len(l_pZ) < nb_chi_p):
        l_pZ.append(__next_suitable_prime(l_pZ[-1], div_pm1));
    l_chi_p  = sum((cf_orbit_p(K, _p)[:1] for _p in l_pZ), []); # Choose one arbitrary ideal over each p
    t = cputime(t); print ("\t[done] t={:.2f}".format(t), flush=True);
        
    # ----------------------------------------------------
    # Compute characters values
    print ("Compute character values in GF(2) on B_ext (#elts={})".format(len(ext_B)), end='', flush=True);
    t = cputime();
    chi_B    = matrix(GF(2), [[log_chip_quadratic(_pid, _a) for _pid in l_chi_p] for _a in ext_B ]);
    t = cputime(t); print ("\t[done] t={:.2f}".format(t), flush=True);
    print ("Extract values on S-units + tors\t\t", end='', flush=True);
    t = cputime();
    val_extB = block_matrix(GF(2), [[ zero_matrix(1 + len(uc), len(fb_d))],
                                    [ Vsu ]], subdivide=False);
    chi_all  = block_matrix(GF(2), [[ val_extB, chi_B ]]);
    t = cputime(t); print ("\t[done] t={:.2f}".format(t), flush=True);

    # ----------------------------------------------------
    # Squares: Kernel
    t = cputime(); H = chi_all.left_kernel().basis_matrix(); t = cputime(t);
    rk = H.rank();
    print("Kernel rk={} t={:.2f}".format(rk, t), flush=True);
    print("Expected rank {} \x1b[{}\x1b[0m".format(exp_rank, "32m[OK]" if (exp_rank == rk) else "31m[NOK]"), flush=True);
    
    # Squares: compute products
    t = cputime();
    H2    = matrix(ZZ, H);
    sq_B = [ prod(ext_B[j]^H2[i][j] for j in range(len(ext_B))) for i in range(H2.nrows()) ];
    t = cputime(t);
    print("Squares t={:.2f}".format(t), flush=True);
    assert(len(sq_B) == rk);
    
    # Squares: reduce mod log U^2
    t = cputime();
    la_sq  = [ logarg_set(_g, p_inf) for _g in sq_B];
    log_sq = [ logarg_log_embedding(_la_sq, p_inf=p_inf, inf_type='TWISTED')*fH for _la_sq in la_sq ];
    sq_red = reduce_mod_logU(sq_B, log_sq, u2, log_U2, log_U2_gso);
    t = cputime(t);
    print("Reduce mod log_U2 t={:.2f}".format(t), flush=True);
    print(("log_2(t2-norms) sq/sq_red: " +"{:.2f}/{:.2f} "*len(sq_B)).format(*(sum(([float(t2_norm(_r).log(2)), float(t2_norm(_r0).log(2))] for _r, _r0 in zip(sq_B, sq_red)), []))), flush=True);

    # -----------------------------------------------------
    # Extracting the square roots
    T = walltime();
    r_dis = cf_sqrtn_fam(sq_red, 2, K, p_inf);
    T = walltime(T);
    print("Square Roots total: t={:.2f}".format(T), flush=True);
    # Someone decided to mix things up, so reorganize this
    r = r_dis;
    t = cputime();
    s_dis = [ _r^2 for _r in r_dis ];
    r     = [ r_dis[s_dis.index(_sq)] for _sq in sq_red ];
    assert(all(_r^2 == _s for _r, _s in zip(r, sq_red)));
    t = cputime(t);
    print("Reorder square roots: t={:.2f}".format(t), flush=True);
    print(("log_2(t2-norms) sq_red/sqrt: " +"{:.2f}/{:.2f} "*len(sq_red)).format(*(sum(([float(t2_norm(_r).log(2)), float(t2_norm(_r0).log(2))] for _r, _r0 in zip(sq_red, r)), []))));

    # -----------------------------------------------------
    # p-adic valuations of the new set
    r_Vp    = matrix(ZZ, (H2[:,1+len(uc):]*Vsu)/2);
    # Finding the new raw representation
    new_Su[_d] = r + Su[_d];
    new_Vp[_d] = block_matrix([[r_Vp],[Vsu]], subdivide=False);
    # Extracting a set of generators
    t     = walltime();
    L, _  = bkz_ZZ(new_Vp[_d],block_size=BKZ_BK);
    # Remove useless rows
    assert(L[:-len(fb_d)] == 0);
    L     = L[-len(fb_d):];
    t     = walltime(t);
    print("BKZ-{}: t={:.2f}".format(BKZ_BK, t), flush=True);

    # Check volume of extracted S-units valuations
    D     = L.determinant().abs();
    ratio = Vsu.determinant().abs()/D;
    print("Vol of image: {}, prev/new: \x1b[{}m{}\x1b[0m".format(D, "32" if (ratio == 2^rk) else "31", ratio.factor()));

    # Bug in FpLLL, the returned transformation matrix is not always correct when coeffs > 2^60
    t     = cputime();
    lu    = new_Vp[_d].solve_left(L, check=False);
    t     = cputime(t);
    assert(L == lu*new_Vp[_d]), "Wrong unimodular matrix"; # This assertion is very important to keep
    print("Find transfo: {:.2f}".format(t), flush=True);
    print("Height(transfo): 2^{:.2f}".format(float(lu.height().log(2))), flush=True);
    new_Ysu[_d] = lu;  

    # -----------------------------------------------------
    # Concatenate things in the correct order for agglomerated output
    B_orb_d   = uc + sum(new_Su[:_d+1], []);
    Vp_orb_d  = block_matrix([[zero_matrix(len(uc),n*(_d+1))],
                              [block_diagonal_matrix(new_Vp[:_d+1])]],
                             subdivide=False).rows();
    Yu_orb_d  = block_matrix([[identity_matrix(len(uc)), zero_matrix(len(uc), len(B_orb_d)-len(uc))]],
                             subdivide=False).rows();
    Ysu_orb_d = block_matrix([[zero_matrix(n*(_d+1),len(uc)),
                               block_diagonal_matrix(new_Ysu[:_d+1])]],
                             subdivide=False).rows();
    
    # -----------------------------------------------------
    # Output !!
    out_sat = data_dir + "z{}_d{}.sat".format(m, _d+1);
    t = cputime(); sunits_raw_write_data(out_sat, K, Yu_orb_d, Ysu_orb_d, B_orb_d, Vp_orb_d); t = cputime(t);
    print("--> Output saturated orb#{} in '{}' t={:.2f}".format(_d+1, out_sat, t));

exit();

