# Code in support of ePrint:2021/xxxx
# Copyright 2021, Olivier Bernard, Andrea Lesavourey
# GPL-3.0-only (see LICENSE file)

import fp
from lattice import *      # svp / bkz
from nf import *           # NF_INFP_BPREC_SCALE
from cyclotomics import *  # cf_orbit_p()


# --------------------------------------------------------------------------------------------------
# This only keeps Tw-PHS strategy, with options 


# --------------------------------------------------------------------------------------------------
# 1. FACTOR BASE CHOICE
#    Dimensioning depends on:
#          ClK (FB must generate),
#          Disc K (phs),
#          Scaling 'c' (PHS/OPT),
#          Targeted root volume (OPT): there is an *unused* __LNROOTVOL_TARGET set to 0.3.
#          RAND (PHS)
#          Class_number hK and regulator RK (TW). Note: for cyclotomics, HR can be found using zeta_K.
#              (for NTRU Prime, it seems more complicated, which is weird)

#  PHS: k = floor( ln |dK| ) - vu )
#       c = n^{3/2} / k
#       random k prime ideals amongst the set of the 3k of smallest norm.

#  TW:  Adds ideals until V_tw^{1 / vu+k} = [ hR.2^{-r2/2}.sqrt(n+k).prod ln Np ]^{1 / vu+k} increases.

# In any case, the following code is designed for *any* cyclotomic field, even of large degree.
# A consequence of this is that we don't guarantee anymore that ideals in fb generate Cl_K.
# Heuristically, this should be ok whp, since:
# - most class groups are cyclic or almost cyclic,
# - it is likely that the first orbit of non principal ideals generate ClK for a substantial fraction
#   of conductors (cf [DPW19, Hyp.7]).


# PHS strategy [PHS19]
# --------------------------------------------------
#     k = floor(ln |dK|) - vu (as in the code of [PHS19])
#     c = n^{3/2} / k         (paper)
# Rk: According to the paper, k should be sthg like ln(D)+n.lnln(D)-rkU (almost twice as large)
def phs_get_dimfb_scale(K, hR='unused', b_prec=fp.BIT_PREC_DEFAULT):
    Re = RealField(b_prec);
    # k = floor(ln |dK|) - vu
    _rk_u = get_rank_units(K);
    _dK   = K.discriminant().abs();
    _k    = ZZ(floor(ln(_dK))) - _rk_u; # Don't do floor( RDF(ln(_dK)) ) --> +Infinity most of the time
    # c = n^{3/2}/k
    _n = K.degree();
    _c = Re(_n)**Re(3/2) / Re(_k);
    return _k, _c;


# Tw-PHS strategy [BR20]
# --------------------------------------------------
#   Adds ideals until the following reduced volume increases:
#           V_tw^{1 / vu+k} = [ hR.2^{-r2/2}.sqrt(n+k).prod ln Np ]^{1 / vu+k} increases.
#   Rk: we don't cut in the middle of iso-norm ideals.
#   We also don't guarantee that fb contains all generators of ClK.
def tw_get_dimfb_norb(K, hR=0, b_prec=fp.BIT_PREC_DEFAULT):
    Re     = RealField(b_prec);
    n      = K.degree();
    r1, r2 = K.signature();
    un_rk  = get_rank_units(K);
    lnhR   = Re(ln(hR)) if hR != 0 else Re(ln(analytic_hR(K)));
    kphs, _= phs_get_dimfb_scale(K);
    # Sufficiently many orbits so that there is kphs ideals (ie. ceil(kphs / K.degree()) orbits)

    # It appears the max density is reached *beyond* kphs for high degrees.
    fb_all = reduce(lambda l,p:(l[0]+[l[1]]*n, K.next_split_prime(l[1])),
                    range((3*kphs-1)//n+1), ([], K.next_split_prime(2)))[0]; 
    S      = Re(log(log(fb_all[0])));
    k      = 1;
    while (Re(fb_all[k].log().log()) < 1/(un_rk+k)*(ln(sqrt(n+k))-ln(2)*r2/2.+lnhR+S) - ln(sqrt(1+1/(n+k)))):
        S = S + Re(log(log(fb_all[k])));
        k = k+1;
        
    # Include all isonorm ideals of the last one.
    pmax = fb_all[k-1];
    while (fb_all[k] == pmax):
        k=k+1;
    
    return ZZ(k / n);


# -------------------------------------------------------------------------------------------------
# 2. GET LOG S-UNIT LATTICE BASIS
#    Depends on:
#        - Chosen logarithmic embeddings (see nf.py for option details)
#        - Isometry/Projection choices: get_twfHcE_matrix()

# Projection / full-rank isometry wizardry 
# ---------------------------------------------

# The method naively applied from [BR20,§4 or §3] is very costly in Sage.
# Instead, this compute the matrix directly from the proof of [BR20,Pr.3.4] (an extended version)
# In the TWISTED case, this is the normalized GSO of (-1 1... ...) with r1+r2+k-1 rows, r1+r2+k columns,
#     for i = 1 .. r1+r2+k-1: bi = (-1/sqrt(i . i+1) {i times}, sqrt(i/(i+1)), 0 ...)
# In the EXPANDED case, this is closer to Pr.3.4:
#     for j = 1 .. r1-1: # occurs only if r1 > 1
#              i=j,        bj   = (-sqrt(1/i/(i+1)) i times,   sqrt(i/(i+1)), 0 ... )     [as above]
#     if r1 = 1:           b_r1 = (-sqrt(2/3), 1/sqrt(6), 1/sqrt(6), ... )
#     if r1 > 1:           b_r1 = (-sqrt(2/r1/(r1+2)) r1 times,sqrt(r1/2/(r1+2)) 2 times, 0 ... )
#                                                                                         [as below]
#     for j = r1+1 to r1+r2-1:
#              i=2j-r1     bj   = (-sqrt(2/i/(i+2)) i times,sqrt(i/2/(i+2)) 2 times, 0 ... )
#     for j = r1+r2 to r1+r2 + k-1,
#              i=j+r2      bj   = as first case
def get_minkH(r1, r2, k, inf_type='EXPANDED', b_prec=fp.BIT_PREC_DEFAULT):
    assert(inf_type in NF_LOG_TYPE and inf_type != 'NONE');
    assert(r1+2*r2 > 1);
    _Re    = RealField(b_prec);

    if (inf_type == 'TWISTED') or (inf_type == 'FLAT'):
        _dim = r1+r2+k;
        _OmH = matrix(_Re, [ vector([-sqrt(1/_Re(_j)/_Re(_j+1))]*_j
                                    + [sqrt(_Re(_j)/_Re(_j+1))]
                                    + [_Re(0)]*(_dim-_j-1))
                             for _j in range(1, r1+r2+k) ]);
    elif (inf_type == 'EXPANDED'):
        _dim = r1+2*r2+k;
        _Br1 = [ vector(_Re, [-sqrt(1/_Re(_i)/_Re(_i+1))]*_i
                        + [sqrt(_Re(_i)/_Re(_i+1))] + [_Re(0)]*(_dim-_i-1))
                 for _i in range(1, r1) ];
        if r1 == 1: # Then r2 > 0
            _Br1 += [ vector([-sqrt(_Re(2)/_Re(3)), _Re(1)/sqrt(_Re(6)), _Re(1)/sqrt(_Re(6))] + [0]*(_dim-3)) ];
        if r1 > 1 and r2 > 0:
            _Br1 += [ vector([-sqrt(_Re(2)/_Re(r1)/_Re(r1+2))]*r1
                             + [sqrt(_Re(r1)/_Re(2)/_Re(r1+2))]*2 + [0]*(_dim-r1-2)) ];
        _Br2 = [ vector(_Re, [-sqrt(_Re(2)/_Re(_i)/_Re(_i+2))]*_i
                        + [sqrt(_Re(_i)/_Re(2)/_Re(_i+2))]*2 + [0]*(_dim-_i-2))
                 for _i in range(r1+2, r1+2*r2, 2) ];
        _Brk = [ vector(_Re, [-sqrt(1/_Re(_i)/_Re(_i+1))]*_i
                        + [sqrt(_Re(_i)/_Re(_i+1))] + [_Re(0)]*(_dim-_i-1))
                 for _i in range(r1+2*r2, r1+2*r2+k) ];
        _OmH = matrix(_Br1+_Br2+_Brk);
    
    # Verify prec
    _OmH   = _OmH.transpose();
    assert (_OmH.base_ring().precision() >= b_prec);
    return _OmH;


# Projection on H = <1,...,1> (k times) inter E (conjugate coordinates are equal)
# ----------------------------------------------
# If the input has coordinates already paired (so, in E), proj_H is sufficient
# Otherwise, we need to project on E also to pair coordinates.
def get_proj_HE(r1, r2, k, inf_type='EXPANDED', b_prec=fp.BIT_PREC_DEFAULT):
    _d     = r1+2*r2+k if (inf_type == 'EXPANDED') else r1+r2+k;
    _pH_E  = (-1/RealField(b_prec)(_d))*ones_matrix(_d);
    _bk_r2 = matrix(RealField(b_prec), [[1/2,1/2],[1/2,1/2]]);
    _pH_E += (block_diagonal_matrix(identity_matrix(r1), *([_bk_r2]*r2), identity_matrix(k), subdivide=False) if (inf_type == 'EXPANDED') else identity_matrix(_d));
    
    assert (_pH_E.base_ring().precision() >= b_prec);
    return _pH_E;
    

# Projection/isometry thing.
# --------------------------------------------------
# a. Full-rank strategy (isometry):
#        True:  Isometry on the whole space, sending 1,...,1 and duplicates to 0 (TW)
#        False: Do nothing.
# b. Log embedding type: 'TWISTED' or 'EXPANDED' or 'FLAT' (only 'EXPANDED' has a different behaviour)
#
# Default values correspond to [BR20]
def get_twfHcE_matrix(r1, r2, k, inf_type='EXPANDED', isometry=True, b_prec=fp.BIT_PREC_DEFAULT):
    assert(inf_type in NF_LOG_TYPE and inf_type != 'NONE');
    n  = r1+2*r2;
    nu = r1+r2-1;
    
    # Domain dimension
    if (inf_type == 'EXPANDED'):
        dim_KR = n+k;
    elif (inf_type == 'TWISTED'):
        dim_KR = (r1+r2+k);

    # Full-rank strategy
    # ----------------------------------------
    # Actually, in Tw-PHS, we have fH(pH(v)) = fH(v) for any vector, so pH*fH = fH
    # So, if we have an isometry: fHcE = fH (projection is implied)
    #     otherwise:              fHcE = pH (no isometry but we still have to project)
    if (isometry == False): # If no isometry: just apply projection
        fHcE = get_proj_HE(r1, r2, k, inf_type=inf_type, b_prec=b_prec);
        assert (fp.fp_check_zero("PrH(1)", vector([1]*dim_KR)*fHcE, target=b_prec));  
    elif (isometry == True): # MK-isometry on n+k (whole space)
        _t = cputime(); fHcE = get_minkH(r1, r2, k, inf_type=inf_type, b_prec=b_prec); _t = cputime(_t);
        # This assertion takes way too long in high dimensions.
        # assert (fp.fp_check_zero("Vol(fH)==1", [vol(fH.transpose())-RealField(b_prec)(1)],
        #                         target=b_prec, sloppy=True)); # /!\ Sloppy
    
    assert (fHcE.base_ring().precision() >= b_prec);
    return fHcE;



# Get Log S-unit lattice basis (compact representation)
# --------------------------------------------------
# la_su:  Log-arg representations of S-Units (a maximal set of independent...)
# p_inf:  Set of infinite places.
#             The aim is to have a consistent order between several calls or Sage sessions
#             [It seems consistent always, but there is no guarantee in the documentation]
#             Note that p_inf[].codomain().precision() should be sufficient to handle sun coefficients.
#             [This is handled internally, but it saves some time.]
# fb:     List of finite places (factor base of prime ideals)
# b_prec: Output bit precision (of course, it should not exceed the precision of p_inf).
def twphs_get_matrix_logarg(la_su, p_inf, fb, inf_type='EXPANDED', fb_type='TWISTED',
                            isometry=True, b_prec=fp.BIT_PREC_DEFAULT):
    assert (all(len(_la_su.vp) == len(fb) for _la_su in la_su)
            and all(len(_la_su.inf) == len(p_inf) for _la_su in la_su)
            and (p_inf[0].codomain().precision() >= b_prec)
            and inf_type in NF_LOG_TYPE and fb_type in NF_LOG_TYPE);
    K      = p_inf[0].domain();
    r1, r2 = K.signature();
    n_inf  = get_nb_inf_places(K);
    k      = len(fb);
    
    # Transformation to full rank
    print ("Get fHcE [inf:{}/iso:{}]".format(inf_type,isometry), flush=True);
    _t = cputime();
    fHcE  = get_twfHcE_matrix(r1, r2, k, inf_type=inf_type, isometry=isometry, b_prec=b_prec);
    _t = cputime(_t);
    print ("[done] t={:.2f}".format(_t), flush=True);
    
    # Log_sun = log_embeddings(la_su)*fHcE;
    Log_sun  = [ ];
    print ("Compute phi(lsu)*fHcE", flush=True); _i = 0;
    for _la_su in la_su:
        print ("\tL#{}".format(_i), end='', flush=True);
        _t = cputime();
        _lsu = logarg_log_embedding(_la_su, p_inf, fb=fb, inf_type=inf_type, fb_type=fb_type);
        _lsu = _lsu*fHcE;
        _t = cputime(_t); print("\t[done] log t={:.2f}, prec={}".format(_t, _lsu.base_ring().precision()), flush=True);
        Log_sun.append(_lsu);
        _i = _i+1;
        
    twphs = matrix(Log_sun);
    assert (twphs.base_ring().precision() >= b_prec);
    return twphs;



# -------------------------------------------------------------------------------------------------
# 3. GET TARGET AFTER CLDL
from random_log_elt import* # for "norm_from_logarg"

# Taken from [BR20] code and adapted to logargs
def twphs_guess_beta(la, Vred, fb, inf='EXPANDED'):
    assert(inf in ['EXPANDED','TWISTED']);
    Re = la.inf[0].base_ring();
    _n = 2*len(la.inf) if inf == 'EXPANDED' else len(la.inf);
    _k = len(fb);
    _sum = sum(ln(Re(pid_fast_norm(_pid))) for _pid in fb);
    
    beta= Re( (_k+_n)/_n*Vred + norm_from_logarg(la, fb)/_n - _sum/_n - ln(pid_fast_norm(fb[0])));
    return beta;


# computes drifted target in log_embedding from logarg la and drift beta
def twphs_target_drift(la, beta, p_inf, fb, fHcE, inf_type='EXPANDED'):
    Re = RealField(p_inf[0].codomain().precision());
        
    log_t   = logarg_log_embedding(la, p_inf, fb, inf_type=inf_type);
    l_drift = vector(Re, [0]*(len(log_t)-len(fb)) + [beta]*len(fb));
    t = (log_t + l_drift)*fHcE;
    return t;


# build logarg of solution from result of CVP in log-embedding
def twphs_build_solution(la, y, U, SU_logarg):
    v = - y*U;
    lcvp = logarg_mpow(SU_logarg, v);
    ls   = logarg_mul2(la, lcvp);
    return ls;


# Tw-query
EXPANSION=1.15;
def twphs_protocol2(la, p_inf, fb, fHcE, lsubkz, u_bkz, SU_logarg, nb_iter, inf_type='EXPANDED', G=0, b_prec=fp.BIT_PREC_DEFAULT):
    Re = RealField(b_prec);
    assert(G!=0);
    
    # Norm challenge
    Na   = round(exp(norm_from_logarg(la, fb)));
    # Reduced volume
    Vred = vol_reduced(lsubkz, gso=G);

    # Set "BEST" at the worst possible element
    best_norm = Na*sqrt(Re(2*len(la.inf))); #(logarg_t2_norm_cf(_ls));
    best_s    = Na;
    s_in_a = False;
    beta_cur = Re(twphs_guess_beta(la, Vred, fb, inf=inf_type))/Re(EXPANSION);  assert(beta_cur > 0);
    while (s_in_a == False):
        beta_cur = beta_cur*Re(EXPANSION);
        print ("\t\tbeta={:8.6f}".format(float(beta_cur)), flush=True, end='');
        _t   = cputime();
        t    = twphs_target_drift(la, beta_cur, p_inf, fb, fHcE, inf_type=inf_type);
        v, y = cvp_babai_NP(lsubkz, t, G=G);
        ls   = twphs_build_solution(la, y, u_bkz, SU_logarg);
        ns   = (logarg_t2_norm_cf(ls));
        _t   = cputime(_t);
        s_in_a  = all(_ls_vp >= 0 for _ls_vp in ls.vp);
        print ("\t[{}]\tl2={:7.3f} t={:.2f}".format(
            s_in_a, float(ns), _t), flush=True);
    
    # Fine grained sampling
    beta_sup = beta_cur;
    beta_inf = beta_cur;
    k=0;
    while (s_in_a == True) and (k<10): # sometimes beta=0 is ok, so avoid infinite loop
        beta_inf = beta_inf/Re(1.4);
        print ("\t\tbeta={:8.6f}".format(float(beta_inf)), flush=True, end='');
        _t   = cputime();
        t    = twphs_target_drift(ls, beta_inf, p_inf, fb, fHcE, inf_type=inf_type);
        v, y = cvp_babai_NP(lsubkz, t, G=G);
        tst  = twphs_build_solution(ls, y, u_bkz, SU_logarg);
        nt   = (logarg_t2_norm_cf(tst));
        _t   = cputime(_t);
        s_in_a  = all(_ls_vp >= 0 for _ls_vp in tst.vp);
        print ("\t[{}]\tl2={:7.3f} t={:.2f}".format(
            s_in_a, float(ns), _t), flush=True);
        if (s_in_a) and (nt < ns): ls = tst; ns = nt;
        k=k+1;
    
    nb_intervals = nb_iter;
    cur_t     = ls;
    best_norm = logarg_t2_norm_cf(ls);
    best_s    = ls;
    for _i in range(1,nb_intervals):
        beta_cur = beta_inf + Re(_i)*(beta_sup-beta_inf)/Re(nb_intervals);
        print ("\t\tbeta={:8.6f}".format(float(beta_cur)), flush=True, end='');
        _t   = cputime();
        t    = twphs_target_drift(cur_t, beta_cur, p_inf, fb, fHcE, inf_type=inf_type);
        v, y = cvp_babai_NP(lsubkz, t, G=G);
        ls   = twphs_build_solution(cur_t, y, u_bkz, SU_logarg);
        ns   = (logarg_t2_norm_cf(ls));
        _t   = cputime(_t);
        s_in_a  = all(_ls_vp >= 0 for _ls_vp in ls.vp);
        print ("\t[{}]\tl2={:7.3f} t={:.2f}".format(
            s_in_a, float(ns), _t), flush=True);
        
        if (s_in_a == True):
            best_s, best_norm = (ls, ns) if ns < best_norm else (best_s, best_norm);
    cur_t = best_s;
        
    return best_s, best_norm;


# //-- END OF FILE
