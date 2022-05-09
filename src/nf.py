# Code in support of ePrint:2021/1384
# Copyright 2021, Olivier Bernard
# GPL-3.0-only (see LICENSE file)

from sage.all import *
import fp
from pideal import *


# Hopefully, always use the (e)GRH hypothesis and useful unproven conjectures for number fields
# We will still try to specify it each time it is needed, for debug and checks
proof.number_field(False);


# -------------------------------------------------------------------------------------------
# Infinite places

# There seemed to exist a precision drift somewhere (?).
# At the time, doubling precision on infinite places corrected the phenomenon.
NF_INFP_BPREC_SCALE = 2.0;


# Returns the number of infinite places of K (r1+r2)
def get_nb_inf_places(K):
    return sum(K.signature());


# Returns infinite places of suitable precision to reach b_prec. [In practice, 2*b_prec]
# Rk: At some point, there was some precision issues with complex embeddings (?).
def get_inf_places(K, b_prec=fp.BIT_PREC_DEFAULT):
    # NB: K.complex_embeddings() is much *MUCH* faster; the ordering is as follows: (real to verify)
    #         complex embeddings: _[0] == _[-1].conjugate(), etc
    #                             decreasing real part, positive imaginary part from _[0] to _[n/2-1]
    #         real embeddings     increasing real part
    # For K.places(), increasing real part, positive imaginary
    # --> all (log( abs(places[i](K.gen())-cembed[n/2-i-1](K.gen())),2) < -prec for i in range(n/2))

    # Attempt...
    # Problem is, complex_embeddings() are not ordered the same way for cyclotomic and generic fields.
    _r1, _r2 = K.signature(); 
    _p_inf   = K.complex_embeddings(prec=NF_INFP_BPREC_SCALE*b_prec);
    assert(len(_p_inf) == K.degree());
    _p_inf   = _p_inf[2*_r2:] + _p_inf[_r2-1::-1];
    assert(len(_p_inf) == _r1+_r2);

    return _p_inf;


# Reduces the precision of infinite places
# This is done by recomputing new places, and verifying the order matches.
def adapt_inf_places(K, p_inf, to_prec=fp.BIT_PREC_DEFAULT):
    assert(len(p_inf) == get_nb_inf_places(K));
    from_prec = p_inf[0].codomain().precision();
    _new_p_inf = get_inf_places(K, b_prec=to_prec); # K.places(prec=to_prec*NF_INFP_BPREC_SCALE);
    assert (fp.fp_check_zero("phi-- - phi", [ _new_p_inf[_k](K.gen())-p_inf[_k](K.gen()) for _k in range(len(p_inf)) ], target=min(from_prec, to_prec)));
    return _new_p_inf;


# T2-norm
# This is the l2-norm of some complex embedding of alpha.
# Unlike Magma, we don't return the square of T2-norm, though it would be an integer (or rational)
def t2_norm(alpha, b_prec=fp.BIT_PREC_DEFAULT):
    # Beware of precision issues with complex_embeddings. Twice should compensate.
    return vector(ComplexField(b_prec),
                  alpha.complex_embeddings(prec=NF_INFP_BPREC_SCALE*b_prec)).norm();



# ------------------------------------------------------------------------------------------
# Units and S-units (generic case)

# Obtain the rank of the (fundamental) unit group
def get_rank_units(K):
    return sum(K.signature())-1;


# -------------------------------------------------------------------------------------------
# Class Number Formula
# /!\ This code apparently only works with cyclotomic fields.
#     Sthg wrong is happening for other fields.
# /!\ Precision issues: this works as is, but experience driven...
#     Nevertheless, it seems sensible that zeta_K should be computed at precision following ln |DK|
#     and that the computation (s-1)*zeta(s) is close enough to the pole for s=1+1/sqrt(|DK|)
def analytic_hR(K):
    dK     = K.discriminant().abs();
    r1, r2 = K.signature();
    # Take precision = ln(dK), distance to 1 eps = 1/sqrt(dK)
    b_prec = ceil(RealField()(log(dK, 2)));# + RealField()(log(exp(150), 2)));
    Re     = RealField(b_prec); # print("zeta:{} {} {}".format(b_prec, b_prec/2, b_prec*0.5));
    zs     = K.zeta_function(prec=b_prec);
    s      = Re(1+2**Re(-0.5*b_prec));# Re(1+exp(-b_prec*0.35));
    res_s  = Re(zs(s)); 
    # Apply Class Number Formula
    hr     = (s-1)*res_s*K.number_of_roots_of_unity()*sqrt(Re(dK))/2**Re(r1)/(2*Re(pi))**Re(r2);

    return hr;



# -------------------------------------------------------------------------------------------
# Number field elements representation.
from collections import namedtuple

# Log-Arg representation
# --------------------------------------------------
# Log-Arg representation: keep ln |s_i(.)| for each infinite place s_i.
#                              arg s_i(.)  same, as elements in [-1/2,1/2] (implicitly mult. by 2*Pi).
#                         (?)  v_p (.)     for p in some FB.

# Log-Arg functions take, whenever needed, the list of infinite/finite (FB) places to consider.
# Going to LogArg is easy, pulling it back is done via polynomial interpolation (to be used wisely).

logarg = namedtuple('logarg',
                    ['inf',  # ln | s_i() |
                     'args', # arg[ s_i() ] / Pi \in [-1,1]
                     'vp',
                    ]);

# K <----> LogArg
# ---------------
def logarg_set(elt, p_inf, fb=[], vp=[]):
    K = elt.parent();
    assert(len(p_inf) == get_nb_inf_places(K));
    #assert(p_inf[0].codomain().precision() >= b_prec*NF_INFP_BPREC_SCALE);
    b_prec = p_inf[0].codomain().precision();
    Re     = RealField(b_prec);
    
    inf_embeddings = [ _phi(elt) for _phi in p_inf ];
    la_inf  = vector(Re, [ _si_elt.abs().log()      for _si_elt in inf_embeddings]);
    la_args = vector(Re, [ _si_elt.arg() / Re(2*pi) for _si_elt in inf_embeddings]);
    # If vp's are given, use them.
    # If factor base is empty as well as vp's, vp field stays empty.
    la_vp   = vector(ZZ, [ elt.valuation(_pid) for _pid in fb ]) if len(vp) == 0 else vector(ZZ, vp);
    la_elt  = logarg(inf  = la_inf, args = la_args, vp = la_vp);
    
    return la_elt;


# compute t2 norm of alpha defined by logarg la
# WARN: Works **only** with cyclotomic fields at this point
def logarg_t2_norm_cf(la):
    Re     = la.inf[0].parent();
    t2_sq  = Re(2) * sum( exp(Re(2)*_la_inf) for _la_inf in la.inf );
    t2_n   = sqrt(t2_sq);
    return t2_n;


# return ln N(b) from logarg of <a> = b . prod_{p in FB} p^vp
# WARN: Works **only** with cyclotomic fields at this point
#       Suppose b is coprime with FB.
def logarg_lnSnorm_cf(la, fb):
    assert (len(fb) == len(la.vp));
    Re     = la.inf[0].parent();
    b_prec = Re.precision();
    ln_Nfb = [Re(pid_fast_norm(_pid).log(prec=b_prec)) for _pid in fb];
    ln_Na  = Re(2)*sum(_la for _la in la.inf) - sum(_vp*_ln_Nfb for _vp, _ln_Nfb in zip(la.vp,ln_Nfb));
    Na     = exp(ln_Na);
    assert(fp.fp_check_zero("exp ln Nb in ZZ", [Na-round(Na)], target=b_prec, sloppy=True));
    return ln_Na;


# /!\ WARN: lifts only in the equation order ZZ[a] (K = Q(a))
def logarg_lift_eqn_order(la, p_inf):
    # Number Field
    K      = p_inf[0].domain(); assert(len(p_inf) == get_nb_inf_places(K));
    r1, r2 = K.signature();
    z      = K.gen();
    # Complex domain: asserting needed precision
    w_prec = p_inf[0].codomain().precision();
    Ce     = ComplexField(w_prec);
    Ce_x   = PolynomialRing(Ce, name='C');
    # 1.443 ~ 1/ln(2), so bit prec is height of la*1.443 + size of the discriminant
    # This should be very safe.
    assert(ceil(1.443*max(map(abs,la.inf)) + RR(K.discriminant().abs().log(2)))  < w_prec), "Precision of infinite places ({}) not sufficient (expected:{})".format(w_prec, ceil(1.443*max(map(abs,la.inf)) + RR(K.discriminant()).abs().log(2)));

    # Evaluation points: phis(z)
    lpts  = [ _phi(z) for _phi in p_inf ] + [ _phi(z).conjugate() for _phi in p_inf[r1:]];
    # Values of the polynomial, computed from the logarg representation
    lvals = (  [ exp(_a+Ce(2*pi*I)*_th) for _a, _th in zip(la.inf[:r1], la.args[:r1]) ]   # Real places
             + [ exp(_a+Ce(2*pi*I)*_th) for _a, _th in zip(la.inf[r1:], la.args[r1:]) ]   # Im   places
             + [ exp(_a-Ce(2*pi*I)*_th) for _a, _th in zip(la.inf[r1:], la.args[r1:]) ]); # Conjugates
    # Interpolation 
    gx   = Ce_x.lagrange_polynomial([(_pt,_v) for _pt, _v in zip(lpts,lvals)]);
    assert(fp.fp_check_zero("gx(pts) == vals", vector(map(gx, lpts))-vector(lvals),
                            target=w_prec, sloppy=True)); # Interpolation is quite stable, but still
    gx   = list(gx);
    
    # Sanity checks and map to ZZ[x]--> O_K
    assert(fp.fp_check_zero("gx in RR[x]", [ _gx_C.imag_part() for _gx_C in gx],
                            target=w_prec, sloppy=True)); # Idem
    gx   = [ _gx_C.real_part() for _gx_C in gx ];
    gx_Z = [ _gx_R.round()     for _gx_R in gx ];
    assert(fp.fp_check_zero("gx in ZZ[x]", [_gx_R - _gx_Z for _gx_R, _gx_Z in zip(gx,gx_Z)],
                            target=w_prec, sloppy=True));
    # Little 'Sagerie': omitting padding only works for cyclotomic fields 
    g    = K(gx_Z + [0]*(K.degree()-len(gx_Z)));
    
    return g;



# Arithmetic in Log-Arg representation
# ------------------------------------
def logarg_reduce_args(c_args):
    # Apply mod [-Pi,Pi] / 2Pi : [n,n+1] --> round [-1/2,1/2]
    return vector(map(lambda _th: _th-_th.round(), c_args));


def logarg_is_equal(la1, la2, sloppy=False):
    w_prec = la1.inf[0].parent().precision();
    inf_eq = fp.fp_check_zero("ln|s1|-ln|s2|", la1.inf - la2.inf, target=w_prec, sloppy=sloppy);
    arg_eq = fp.fp_check_zero("theta_1-theta_2 mod 1", [min(_th.abs(), 1-_th.abs()) for _th in nfc_reduce_args(la1.args) - nfc_reduce_args(la2.args)], target=w_prec, sloppy=sloppy);
    vp_eq  = (la1.vp == la2.vp);
    
    return (inf_eq and arg_eq and vp_eq);


def logarg_mul2(la1, la2):
    mul_inf  = la1.inf  + la2.inf;
    mul_args = logarg_reduce_args(la1.args + la2.args);
    mul_vp   = la1.vp   + la2.vp;
    
    la_mul   = logarg(inf=mul_inf, args=mul_args, vp=mul_vp);
    return la_mul;


# Naive accumulation.
# Use a "product(sum)" tree if there is too much precision drift, but tests seem very fine as is.
def logarg_mul(l_las):
    assert(len(l_las) > 0); # Lazy enough not to define "1".
    
    la_mul = l_las[0];
    for _la in l_las[1:]:
        la_mul = logarg_mul2(la_mul, _la);

    return la_mul;


def logarg_pow(la, k):
    pow_inf  = k*la.inf;
    pow_args = logarg_reduce_args(k*la.args);
    pow_vp   = k*la.vp;
    la_pow   = logarg(inf=pow_inf, args=pow_args, vp=pow_vp);

    return la_pow;


def logarg_mpow(l_las, l_exps):
    assert(len(l_las) == len(l_exps));
    
    l_mpow  = [ logarg_pow(_la, _exp) for _la, _exp in zip(l_las, l_exps) ];
    la_mpow = logarg_mul(l_mpow);

    return la_mpow;


# -------------------------------------------------------------------------------------------
# Logarithmic embedding.
# When conjugates are present:
#     the order of infinite embeddings is as follows : 1...r_1,  ..., j, \conj(j), ...
#     the order of finite embeddings follows the same principle: -vp ln p, ..(\times f=[OK/P:Z/p]), -vq etc
# For finite/infinite parts, there are 3 options: TWISTED, EXPANDED and FLAT: (dimensions given if applied on both parts)
#     TWISTED:     [Ks:Qs] ln |a|_s                                      Dim: r_1+r_2  / k
#     FLAT:        ln |a|_s on infinite places, -vp(a) on finite places  Dim: r_1+r_2  / k
#     EXPANDED:    ln |a|_s, [Ks:Qs] times                               Dim: r_1+2r_2 / sum([OK/p:Z/p])
# For ideals, we add the following on infinite places:
#     NONE:        Put 0 everywhere.
# Note the "FLAT" version is designed to fit more or less PHS's choice on finite valuations;
# we would naturally expect instead ln |a|_s on all places (= vp ln(p) for finite places).
# /!\ p_inf[].codomain().precision() should be large enough to handle coefficients of elt.
NF_LOG_TYPE = ['TWISTED', 'FLAT', 'EXPANDED', 'NONE'];


# This returns the chosen log embedding from a logarg representation 
# We need the factor base norms for fb_type='TWISTED' or 'EXPANDED'. The simpler is probably to input the factor base directly.
# /!\ We NEED pid ideals <p,g(x)> in FB.
# 
# NB: Compute directly at the precision of p_inf.
def logarg_log_embedding(la, p_inf, fb=[], inf_type='TWISTED', fb_type='TWISTED'):
    K      = p_inf[0].domain();
    w_prec = p_inf[0].codomain().precision();
    Re     = RealField(w_prec);
    assert(len(la.inf) == len(p_inf) and len(p_inf) == get_nb_inf_places(K));
    assert(len(la.vp)  >= len(fb));
    assert((inf_type in NF_LOG_TYPE) and (fb_type in NF_LOG_TYPE));
    
    r1, r2 = K.signature();
    # Several remarks:
    # 1. abs_val return |s(eta)|^[Ks:R] on infinite places, N(p)^(-v_p(eta)) on finite places
    # 2. the log function has some weird precision issues (?)
    # 3. going through finite places with abs_val() takes too much time in Sage;
    #    computing it from scratch (valuation, norm) can be up to 3 times faster
    ln_inf_vals = list(la.inf);
    if (inf_type == 'FLAT'):
        ln_inf = ln_inf_vals;
    elif (inf_type == 'TWISTED'):
        ln_inf = ln_inf_vals[:r1] + [ 2*_ln_s for _ln_s in ln_inf_vals[r1:] ];
    elif (inf_type == 'EXPANDED'):
        ln_inf = ln_inf_vals[:r1] + sum((2*[_ln_s] for _ln_s in ln_inf_vals[r1:]), []);

    padic_vals = list(la.vp);
    if (fb_type == 'FLAT'): # PHS's version, but logically this would be -vp(eta) ln(p) (vs. -vp(eta) ln(norm(pid)) for TWISTED)
        ln_fb  = [ - _vp for _vp in padic_vals ];
    elif (fb_type == 'TWISTED'):
        ln_fb  = [ - _vp * pid_fast_norm(_pid).log(prec=w_prec) for _vp, _pid in zip(padic_vals, fb)];
    elif (fb_type == 'EXPANDED'):
        ln_fb  = sum((pid_fast_residue_class_degree(_pid) * [- _vp*pid_fast_smallest_integer(_pid).log(prec=w_prec)] for _vp, _pid in zip(padic_vals,fb)), []);
        
    ln_embedding = vector(ln_inf + ln_fb);
    assert (ln_embedding.base_ring().precision() >= w_prec);
    return ln_embedding;


# //-- END OF FILE
