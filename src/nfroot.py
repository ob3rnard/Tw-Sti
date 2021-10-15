# Code in support of ePrint:2021/1384
# Copyright 2021, Andrea Lesavourey
# GPL-3.0-only (see LICENSE file)

from sage.all import *
import fp
from pideal import *
from lattice import *
from nf import *
proof.number_field(False);

# sort family of nf elements according to coefficient norm
# Hv = elements in coeff representation
def nf_coeff_sort(Hv, K):
    Hv.sort(key=norm);
    return [K(list(h)) for h in Hv];
    
# evaluation of precision needed
def prec_eval_sqrtn(g, n):
    return len(g)*ceil(log(g.norm(), 2))//n;

# evaluation of precision needed
def prec_eval_sqrtn_fam(G, n):
    return max([len(g)*ceil(log(g.norm(), 2))//n for g in G]);

# creation of reduced basis for given precision
def reduced_basis_latt(K, Binf, bprec):
    dim = K.degree(); 
    B = [-round( (2**bprec*Binf[i]).real_part()) for i in range(dim//2)];
    coeff = dim // 2;
    M = coeff*identity_matrix(dim//2, ZZ);
    M = M.insert_row(0, B)
    M = M.transpose()
    return lll_ZZ(M)[0];

def init_basis_latt(K):
    dim = K.degree(); 
    B = vector(ZZ, [-1 for i in range(dim//2)]);
    coeff = dim // 2;
    
    M = coeff*identity_matrix(dim//2, ZZ);
    M = M.insert_row(0, B)
    M = M.transpose()
    return M;


# use previous unimodular matrix U to start the reduction
def update_basis_latt(K, U, Binf, bprec):
    dim = K.degree(); 
    # print "done here"
    B = [-round(((2**bprec)*(Binf[i].real_part()))) for i in range(dim//2)];
    coeff = dim // 2;
    V = U/coeff;
    M = coeff*identity_matrix(dim//2, ZZ);
    M = M.insert_row(0, B)
    M = M.transpose()
    return lll_ZZ(V*M)[0];


# decoding by nearest plane and kannan's embedding technique
# for better performances and accuracy replace by proper NP
def test_decode(BL, x_emb, bprec):
    dim = BL.nrows()
    coeff = dim
    r = BL.ncols()-BL.nrows()
    Y = [round((2**(bprec))*x_emb[i]) for i in range(r)];
    m = max(Y)
    Z = list(zero_vector(ZZ, dim))
    Y += Z
    
    M = (BL.transpose()).augment(vector(ZZ, Y)).transpose()

    Z += [2**(bprec//2)];

    M = M.augment(vector(ZZ,Z));

    M = lll_ZZ(M)[0];
    return M[dim][r:M.ncols()-1]/coeff


# d-th root computation of h in K
# assume L=basis lattice for decoding  is precomputed
# assume M=inv relative mink precomputed as well (2x2 matrix)
def cf_sqrtn(h, d, Minf, K, L, pinf, Binf, bprec):
    zeta = K.gen();
    eta = zeta + 1/zeta;
    hinf = pinf[0](h); # put in CC
    ginf = hinf**(1/d);
    gconj = ginf.conjugate();

    # going to rep of g as g = g0 + g1 zeta_m
    ginf  = vector([ginf, gconj]);
    ginf = Minf*(ginf.column()); # 

    bprec_max = pinf[0].codomain().precision()
    
    bprec_new = bprec;
    b = False;
    while (not b):
        
        if bprec_new > bprec_max:
            # print("updating pinf and all subsequent objects");
            bprec_max = 2*bprec_new//3;
            pinf = get_inf_places(K, bprec_max);
            Binf = [pinf[0](eta**k) for k in range(K.degree()//2)];
            zinf = pinf[0](zeta);
            Minf = Matrix([[1, zinf], [1, zinf**(-1)]])**(-1);
            hinf = pinf[0](h); # put in CC
            ginf = hinf**(1/d);
            gconj = ginf.conjugate();
            
            # going to rep of g as g = g0 + g1 zeta_m
            ginf  = vector([ginf, gconj]);
            ginf = Minf*(ginf.column()); #
                
        if bprec_new != bprec:
            U = L[:,1:];
            s = cputime()
            L = update_basis_latt(K, U, Binf, bprec_new).change_ring(ZZ);
            bprec = bprec_new;
            
        g0 = test_decode(L, [ginf[0,0].real_part()], bprec);
        g0 = sum([g0[k]*eta**k for k in range(len(g0))]);
        
        g1 = test_decode(L, [ginf[1,0].real_part()], bprec);
        g1 = sum([g1[k]*eta**k for k in range(len(g1))]);
        
        b = ((g0+g1*zeta)**d==h);
        bprec_new = ceil(sqrt(2)*bprec);
        
    return g0+g1*zeta, L, bprec, pinf, Binf, Minf;


# d-th root computation of family in K
# Warn! the order of output is not the same as input...
def cf_sqrtn_fam(H, d, K, pinf):
    dim = K.degree();
    zeta = K.gen();
    eta = zeta + 1/zeta;
    zinf = pinf[0](zeta);
    Minf = Matrix([[1, zinf], [1, zinf**(-1)]])**(-1);
    N = [vector(eta**k) for k in range(dim//2)] + [vector(zeta*eta**k) for k in range(dim//2)];
    N = Matrix(N);
    
    Hv = [vector(H[i]) for i in range(len(H))];
    Hn = nf_coeff_sort(Hv, K);
    bprec_max = prec_eval_sqrtn_fam(Hv, d) + ceil(log(Minf.norm(), 2))*dim;
    pinf = get_inf_places(K, bprec_max//3)
    
    Binf = [pinf[0](eta**k) for k in range(K.degree()//2)];
    L = init_basis_latt(K);
    R = [];

    bprec = 0;
    
    
    for k in range(len(Hn)):
       
        bprec_new = prec_eval_sqrtn(vector(Hn[k]), d)//5 + ceil(log(N.norm(), 2))*K.degree()//2;
        
        if bprec_new > bprec:
            print("Update prec to {}".format(bprec_new), end='', flush=True);
            bprec = bprec_new;
            t = walltime();
            L = update_basis_latt(K, L[:,1:], Binf, bprec);
            t = walltime(t);
            print("\t[done] t={:.2f}".format(t), flush=True);
            
        print("Computing sqrt(#{})\tlog_2 red norm {:.2f}\t".format(k, float(log(t2_norm(Hn[k]),2))), end='', flush=True);
        t = walltime()
        r, L, bprec, pinf, Binf, Minf = cf_sqrtn(Hn[k], d, Minf, K, L, pinf, Binf, bprec);
        t = walltime(t);
        print("\x1b[32m[OK]\x1b[0m {:.2f} t={:.2f}\n".format(float(log(t2_norm(r),2)),t), end='', flush=True);
        R += [r]
        
    return(R);


# //-- END OF FILE
