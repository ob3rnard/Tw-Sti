# Code in support of ePrint:2021/1384
# Copyright 2021, Olivier Bernard, Andrea Lesavourey
# GPL-3.0-only (see LICENSE file)

from sage.all import *

import fp
from nf import *
from pideal import *

from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler


#####################################################################
#                        FONCTIONS                                  #
#####################################################################

# idea for random ln |si| s.t. Sum_i ln |si| = ln |N(a)| is the following:
# compute random element in H = ortho(RR*(1,...,1)) then add (ln |N(a)|)/n*(1,...1)
def proj_log(v):
    Re = v.base_ring();
    s  = sum(v);
    ph = vector([_v - s/Re(len(v)) for _v in v]);
    assert(fp.fp_check_zero("sum(pH u) = 0", [sum(_p for _p in ph)], target=Re.precision()));
    return ph;


# return ln |b| from logarg of target
def norm_from_logarg(la, fb):
    Re     = la.inf[0].parent();
    b_prec = Re.precision();
    ln_Nfb = [Re(pid_fast_norm(_pid).log(prec=b_prec)) for _pid in fb];
    ln_Na  = Re(2)*sum(_la for _la in la.inf) - sum(_vp*_ln_Nfb for _vp, _ln_Nfb in zip(la.vp,ln_Nfb));
    Na     = exp(ln_Na);
    assert(fp.fp_check_zero("exp ln Nb in ZZ", [Na-round(Na)], target=b_prec, sloppy=True)); # Sloppy option barely needed
    return ln_Na;


# return random ClDL (simulated) solution
def random_targets(K, p_inf, fb,
                   chal_bsz = 100, nb_targets = 100, b_prec = fp.BIT_PREC_DEFAULT):

    assert(p_inf[0].codomain().precision() >= b_prec);
    assert(len(p_inf) == get_nb_inf_places(K));

    Re = RealField(b_prec);
    n = K.degree();
    ln_Nfb = vector(Re, [pid_fast_norm(fb[i]).log(prec=b_prec) for i in range(len(fb))]);

    # Sigma = 100*dim (*2^prec is relevant)
    DU   = DiscreteGaussianDistributionIntegerSampler(sigma=100*(n//2-1)*2**b_prec);
    DSti = DiscreteGaussianDistributionIntegerSampler(sigma=100*len(fb));
    list_targets = [];
    for k in range(nb_targets):
        p = random_prime(2**(chal_bsz+3), lbound=2**(chal_bsz-3));
        p = K.next_split_prime(p);
        print ("Chal #{}: p={}".format(k, p), flush=True);
        
        # draw random valuation for fb
        vp = vector(ZZ, [DSti() for i in range(len(fb))]);
        l_inf = vector(Re, [DU()/2**b_prec for i in range(n//2)]);
        l_inf = proj_log(l_inf);

        # Fix infinite valuations to obtain a <t> = p . prod(Ni^vi)
        ln_Nt = Re(p.log(prec=b_prec)) + sum(_vi*_ln_Ni for _vi, _ln_Ni in zip(vp, ln_Nfb));
        l_inf = l_inf + vector([ln_Nt / Re(n)]*(n//2));

        # Set simulated logarg representation
        la = logarg(inf = l_inf, args = vector(Re, [0]*len(l_inf)), vp = vp);
        assert(fp.fp_check_zero("sum a_s - sum a_p = ln Nb",
                                [norm_from_logarg(la, fb) - Re(ln(p))], target=b_prec-50));
        list_targets.append(la);
    
    return(list_targets);


# //-- END OF FILE
