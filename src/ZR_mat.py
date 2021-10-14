# Code in support of ePrint:2021/xxxx
# Copyright 2021, Olivier Bernard
# GPL-3.0-only (see LICENSE file)

from sage.all import *

import fp


# ----------------------------------------------------------------------------------------------
# General "matrix on real" caveats handling
def rc_mat_inverse(M):
    _R      = M.base_ring(); # Real or Complex
    _b_prec = _R.precision();
    _RR     = _R.to_prec(4*_b_prec); # Should be enough, 1 is not sufficient, 2 or 3 might work

    _iM     = ((M.change_ring(_RR)).inverse()).change_ring(_R);

    # Check there is no precision issue (NB: compute the product in high precision)
    _chk_M  = (M.change_ring(_RR)*_iM.change_ring(_RR)).change_ring(_R) - identity_matrix(M.nrows());
    assert (fp.fp_check_zero("M*iM-In", _chk_M.coefficients(), target=_b_prec, sloppy=True)); # /!\ Sloppy

    return _iM;


# Conversions RR <--> ZZ
# Scaling factors are given in log_2 to match "precision" of RR

# Get integer equivalent of (M) [ returns floor(2^prec*M) on R ]
# Returns also the log_2 of the applied scaling factor.
def matvec_real_to_ZZ(Mv, work_prec=0):
    _R         = Mv.base_ring().fraction_field();
    _log_scale = _R.precision() if (work_prec == 0) else work_prec;
    _scale     = Integer(2)**(_log_scale);
    _MZ        = Mv.apply_map(lambda _mij:Integer(floor(_scale*_mij)));

    return _MZ, _log_scale;


# Returns to real numbers
def matvec_ZZ_to_real(Mv, log_scale, to_prec=0):
    assert ((to_prec == 0) or log_scale <= to_prec);

    _prec = log_scale if (to_prec == 0) else to_prec;
    _R    = RealField(_prec);
    _MR   = Mv.apply_map(lambda _mij:_R(_mij)/_R(2)**(log_scale));
    
    assert (_MR.base_ring().precision() >= _prec);
    return _MR;



# ----------------------------------------------------------------------------------------------
# Gram-Schmidt Orthogonalization
# Input: M is m*n matrix on any Ring
# Return G,P: G is GSO of M (normalized if normalize=True), and permutation matrix (of the rows)
#             P is the transition matrix such that M = P*G
# NB: annoyingly, sage does not implement this outside of ZZ/QQ, RDF and CDF.
# NB: b_prec is effective iff M has exact ring AND exact=False
def gram_schmidt_ortho(M, normalize=False, exact=False, b_prec=fp.BIT_PREC_DEFAULT):
    _R     = M.base_ring();
    if (_R.is_exact()):
        _R = _R.fraction_field() if (exact==True) else RealField(b_prec);
        print ("[GSO] **Warn** Coeffs of M in exact ring, moved to '{}'.".format(_R));

    _n = M.nrows();
    _G = M.change_ring(_R);
    _P = identity_matrix(_R, _n);
    
    # Main loop
    # NB: Exchanging _i and _j loop would lead to somewhat "order-independent" algorithm,
    #     allowing to choose the smallest length projection for the next step
    #    for _i in range(1,_n):
    #        _G[_i] = _G[_i] - sum( (_G[_i]*_G[_j])/_G[_j].norm()**2 * _G[_j] for _j in range(_i));
    _G_norm = [_G[0].norm()]*_n;
    for _i in range(1,_n):
        t= cputime();
        # _mu_i       = [(_G[_i]*_G[_j])/_G[_j].norm()**2 for _j in range(_i)] + [1] + [0]*(_n-_i-1);
        _mu_i       = [(_G[_i]*_G[_j])/_G_norm[_j]**2 for _j in range(_i)] + [1] + [0]*(_n-_i-1);
        t1 = cputime(t);
        _G[_i]      = _G[_i] - sum(_mu_i[_j] * _G[_j] for _j in range(_i));
        _G_norm[_i] = _G[_i].norm();
        t2 = cputime(t);
        _P[_i] = vector(_R, _mu_i);#[_mu_i[_j] for _j in range(_i)] + [_R(1)] + [_R(0)]*(_n-_i-1));
        t3 = cputime(t);
        # print ("{}: {:.2f} {:.2f} {:.2f}".format(_i, t1, t2-t1, t3-t2));
        
    # Orthonormalization (not by default)
    if (normalize == True):
        for _i in range(_n):
            _norm_i   = _G_norm[_i];
            _P[:,_i] *= _norm_i;
            _G[_i]   /= _norm_i;
            
    # **Warn** These assertions are not costless
    if not (_R.is_exact()) and (_n < 100):
        assert (fp.fp_check_zero("M-PG", list(M-_P*_G), target=_R.precision()));
    elif (_n < 100):
        assert (M-_P*_G == 0);
    return _G, _P;



# ----------------------------------------------------------------------------------------------
# Volumes
# Return Vol M:
# - if M is on an exact ring (ZZ, Fp, etc), volume is best computed as sqrt |det M*Mt|
# - if M is on a field "with precision" (RR, CC), the above is not numerically stable, so it is
#   better advised to compute it as exp( sum ln ||gi*|| for gi* in M*:=gram_schmidt_ortho(M))
def lnvol(M, gso=0):
    assert (M.nrows() <= M.ncols());
    _gM = gram_schmidt_ortho(M)[0] if (gso == 0) else gso; assert(_gM.nrows() == M.nrows());
    return sum([ ln(_gM[_k].norm()) for _k in range(M.nrows())]);

def vol_exact(M):
    assert(M.base_ring().is_exact());
    if M.is_square():
        _vol = M.determinant().abs();
    else:
        assert (M.nrows() <= M.ncols());        
        _gram_M = M*M.transpose();
        _vol    = _gram_M.determinant().abs().sqrt();
    return _vol;

def vol(M, gso=0):
    if (M.base_ring().is_exact()):
        return vol_exact(M);
    else:
        return exp(lnvol(M, gso=gso));



# Return |Vol M|^(1/n), where n is the matrix number of rows
def lnvol_reduced(M, gso=0):
    return (1/M.base_ring()(M.nrows()))*lnvol(M, gso=gso);

def vol_reduced(M, gso=0):
    if (M.base_ring().is_exact()):
        return vol_exact(M)**(1/M.nrows());
    else:
        return exp(lnvol_reduced(M, gso=gso));

# //-- END OF FILE
