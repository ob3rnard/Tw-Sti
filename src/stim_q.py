# Code in support of ePrint:2021/1384
# Copyright 2021, Olivier Bernard
# GPL-3.0-only (see LICENSE file)

from sage.all    import *
from characters  import *
from cyclotomics import * # cf_get_conductor()



# -------------------------------------------------------------------------------------------------
# KUCERA
# Adapting [Kuc92] "On bases of the Stickelberger ideal and of the group of circular units of a cyclotomic field"
# The "alpha" function corresponds to a new paper:
#          [BK21] "A short basis of the Stickelberger ideal of a cyclotomic field"


# Special sets of indices mod m
# ------------------------------------------------------------
# M- from [Kuc92, p.293]. We decompose it as:
# X (no 0) \ res_q + N- + { B_k for k in 1..g, g=#factor(m) }
def kucera_get_Mminus(m):
    # m = prod q_i, where q_i = p_i^e_i
    m_p = [ _f[0]        for _f in factor(m) ]; # Cheap anyway.
    m_q = [ _f[0]**_f[1] for _f in factor(m) ]; # Cheap anyway.
    g   = ZZ(len(m_q));
    # Start from X = all [1, m-1] prime to m, plus multiples of p_i^e_i (not of p_i, ... only)
    Mminus = [ ZZ(_a) for _a in range(1,m) if
               all( not m_p[_i].divides(_a) or m_q[_i].divides(_a) for _i in range(g)) ];
    # Remove N+
    Nminus = { _a for _a in Mminus
               if _a.divides(m) and is_even(len([_q for _q in m_q if not _q.divides(_a)])) };
    Mminus = [ _a for _a in Mminus if _a not in Nminus ];
    # Remove all a/(m,a)=-1 [q_i] for some q_i
    res_q = { _a for _a in Mminus
            if any([ _q.divides(_a+gcd(m,_a)) for _q in m_q if not _q.divides(_a)]) };
    Mminus = [ _a for _a in Mminus if _a not in res_q ];
    # Remove B_k = { q_k \not| a and {a/q_k/(a,m)} > 1/2 and for i > k, a = (m,a) mod [q_i] }
    for i in range(g):
        q_i = m_q[i];
        B_i = { _a for _a in Mminus if
                not q_i.divides(_a)
                and 2 * Integers()(mod(_a/gcd(_a,m), q_i)) > q_i
                and all( _q.divides(_a-gcd(m,_a)) for _q in m_q[i+1:] ) };
        Mminus = [ _a for _a in Mminus if _a not in B_i ];

    assert (len(Mminus) == euler_phi(m)/2);
    return Mminus;


# Implementing the 2021 version:
# M'_m = M_m \ {a st. m/q_i | a} u [[1,\phi(qi)/2]]
def kucera_get_Mp(m):
    m_p = [ _f[0]        for _f in factor(m) ]; # Cheap anyway.
    m_q = [ _f[0]**_f[1] for _f in factor(m) ]; # Cheap anyway.
    g   = ZZ(len(m_q));

    Mm  = kucera_get_Mminus(m);
    Mp  = [ a for a in Mm if all( not ZZ(m/qi).divides(a) for qi in m_q ) ];
    ext = sum( ([ZZ(m*b/qi) for b in range(1, euler_phi(qi)/2+1)] for qi in m_q), []);
    Mp  = sorted(list(set( Mp + ext )));
    
    assert (len(Mp) == euler_phi(m)/2);
    return Mp;


# The Alpha function
# ----------------------------------------------------------
# Computes \alpha_m(b)
# Returns the ZZ triple [a,b,c] in [0,m[ st. m | a+b+c and :
#         \alpha_m(b) = om(a)+om(b)+om(c)+1/2.Nm = th(a)+th(b)-th(-c)  is *short* (all coeffs in {0,1})
# The relation is then (s_a + s_b + s_c) th(1) / p
def kucera_alpha(m, b, verbose=False):
    m_p  = [ _f[0]        for _f in factor(m) ]; # Cheap anyway.
    m_q  = [ _f[0]**_f[1] for _f in factor(m) ]; # Cheap anyway.
    g    = len(m_q);
    Jb_c = [ i for i in range(g) if not m_q[i].divides(b) ]; assert(len(Jb_c) > 0);
    
    if len(Jb_c) == 1:
        qj = m_q[Jb_c[0]]; pj = m_p[Jb_c[0]];
        c  = ZZ(b*qj/m); assert( 0 < c and c <= euler_phi(qj)/2 );
        if c > 1:
            a_idx = [ ZZ(_i) for _i in [-b, b-m/qj, m/qj]]; # - {b} + {b-m/qj} + {m/qj}
        else: # c==1
            a_idx = [ ZZ(_i) for _i in [ m*euler_phi(qj)/2/qj ]*2 + [m/pj]]; # 2*{ m.\phi(qj)/(2.qj) } + {m/pj}
    else:
        # Solving u_b x + v_b y = -b/r_b
        ub      = m_q[min(Jb_c)];
        rb      = prod(m_q[i] for i in range(g) if i not in Jb_c);
        vb      = ZZ(m/rb/ub);
        (d,s,t) = ub.xgcd(vb); assert(d==1); assert(d==s*ub+t*vb);
        # s = s*(-b/rb); t = t*(-b/rb); 
        A     = ZZ(mod(-vb*t*b,m)); B = ZZ(mod(-ub*s*b,m));
        a_idx = [b, A, B];

    # Mapping everything in [0,m[
    a_idx = [ ZZ(mod(_i,m)) for _i in a_idx ];
    assert(0 == mod(sum(a_idx),m));
    if verbose: print("{} +{} +{}".format(*a_idx));

    return a_idx;


# All \alpha_m(b) for b in S.
# Computing \alpha is cheap so there is no trick here
def kucera_alpha_idx(m, S):
    al = [ kucera_alpha(m, _b) for _b in S ];
    return al;

# All \alpha_m(b) for b in M'm
def kucera_alpha_all(m):
    S  = kucera_get_Mp(m);
    al = kucera_alpha_idx(m, S);
    return al;



# --------------------------------------------------------------------------------------------------
# Special Stickelberger elements in QQ[G]

# Theta computations
# ------------------------------------------
# \theta(a) = \sum_{(r,m)=1, 0<r<m} {- ar/m} \s_r^{-1}
def theta_kucera_a(m,a):
    n =  euler_phi(m);
    G         = [ ZZ(i) for i in range(m) if gcd(i,m)==1 ]; # s_a : a st (a,m)=1
    Ginv      = [ b.inverse_mod(m) for b in G ]; # Ginv[i] : s_{1/G[i] mod m}
    Gminv_idx = [ G.index(b) for b in Ginv];

    v = [0]*n;
    for i in range(len(G)):
        r      = G[i];
        idx    = Gminv_idx[i];
        v[idx] = ZZ(mod(-a*r,m))/ m;

    return vector(v);


# Computation of all thetas is quite costly (because of G/Ginv and QQ computations)
# So it is interesting to compute all of them in batch.
# - The output is st. ta[a] = \theta(a), and notably, ta[0] = 0.
# - The last element is \theta(m//2).
def theta_kucera_all(m):
    n         =  euler_phi(m);
    G         = [ ZZ(i) for i in range(m) if gcd(i,m)==1 ]; # s_a : a st (a,m)=1
    Ginv      = [ b.inverse_mod(m) for b in G ]; # Ginv[i] : s_{1/G[i] mod m}
    Gminv_idx = [ G.index(b) for b in Ginv];

    ta = [zero_vector(n)];
    for a in range(1, m//2+1):
        v = [0]*n;
        for i in range(len(G)):
            r      = G[i];
            idx    = Gminv_idx[i];
            v[idx] = ZZ(mod(-a*r,m))/ m;

        ta.append(vector(v));
        
    return ta;


# Same for a subset S, computes \theta_m(S)
def theta_kucera_idx(m, S):
    n         =  euler_phi(m);
    G         = [ ZZ(i) for i in range(m) if gcd(i,m)==1 ]; # s_a : a st (a,m)=1
    Ginv      = [ b.inverse_mod(m) for b in G ]; # Ginv[i] : s_{1/G[i] mod m}
    Gminv_idx = [ G.index(b) for b in Ginv];

    ta = [];
    for a in S:
        v = [0]*n;
        for i in range(len(G)):
            r      = G[i];
            idx    = Gminv_idx[i];
            v[idx] = ZZ(mod(-a*r,m))/ m;

        ta.append(vector(v));

    return ta;


def theta_washington(m): # This is \theta(1) with Kucera notation
    v = theta_kucera_a(m,-1);
    return v;


# Omega computations
# ------------------------------------------
# omega* = w/2m (\theta(-1) + (m+1)\theta(1)), w=2m if m is odd, m if m is even
#        = w/2m (s(G) + m\theta(1))
def omega_kucera_star(m):
    w = 2*m if is_odd(m) else m;
    o_s = w/2/m * (theta_kucera_a(m,-1) + (m+1)*theta_kucera_a(m,1));
    return vector(ZZ,o_s);


# omega(a) = \theta(a) - a*\theta(1) + omega* - s(G)
# The plain version applies the formula above regardless of a
def omega_kucera_a_plain(m,a):
    os = omega_kucera_star(m);
    sG = ones_matrix(1, len(os))[0];
    ta = theta_kucera_a(m,a);
    t1 = theta_kucera_a(m,1);
    return vector(ZZ, ( ta -a*t1 +os -sG));


# The "normal" version implements [Kuc92]:
# - omega(0) = 0
# - omega(a) = -omega(m-a) if m/2<a<m
# - omega(a) = omega(b) whenever a=b [m]
def omega_kucera_a(m,a):
    _b = ZZ(mod(a,m));
    if (_b == 0):
        return zero_vector(ZZ, euler_phi(m));
    elif (_b <= m/2):
        return omega_kucera_a_plain(m,_b);
    else:
        return -omega_kucera_a_plain(m,m-_b);



# --------------------------------------------------------------------------------------------------
# Stickelberger valuations

# [CDW17] Computations of the v_i and w_i
# ------------------------------------------------------------
# Based on [DPW19], with correction -b --> 1/b for result idx
# (The i-th relation is for s_{1/a}, not s_{-a})
# Output: array V, st.:
#     V[0]  is the first relation (2-s_2),
#     V[i]  is relation a - s_a, a=G[i+1]  (convention G[0]=1)
#     V[-1] is relation "1 - s_1", ie. (m+1)-s_{m+1}
def cdw17_va(m,a):
    n         = euler_phi(m);
    G         = [ ZZ(i) for i in range(m) if gcd(i,m)==1 ]; # s_a : a st (a,m)=1
    Ginv      = [ b.inverse_mod(m) for b in G ]; # Ginv[i] : s_{1/G[i] mod m}
    Gminv_idx = [ G.index(b) for b in Ginv];

    v = [0]*n;
    for i in range(len(G)):
        r      = G[i];
        idx    = Gminv_idx[i];
        v[idx] = floor(a * (r / m));

    return vector(v);

def sti_vals_cdw17_v(m):
    n = euler_phi(m);
    G         = [ ZZ(i) for i in range(m) if gcd(i,m)==1 ]; # s_a : a st (a,m)=1
    Ginv      = [ b.inverse_mod(m) for b in G ]; # Ginv[i] : s_{1/G[i] mod m}
    Gminv_idx = [ G.index(b) for b in Ginv];

    vs = [];
    for a in [_a for _a in range(2,m+2) if gcd(_a,m) == 1]:
        v = [0]*n;
        for i in range(len(G)):
            b      = G[i];
            idx    = Gminv_idx[i];
            v[idx] = floor(a * (b / m));
            
        assert(len(v) == n);
        vs.append(vector(v));

    assert(len(vs) == n); # Rank is n/2 + 1, but we have phi(m) relations
    vs = Matrix(ZZ, vs);
    return vs;

# For now, compute the vs and substract v_i+1 - v_i
def cdw17_wa(m,a):
    return cdw17_va(m,a) - cdw17_va(m,a-1);

def sti_vals_cdw17_w(m):
    vs = sti_vals_cdw17_v(m);
    n  = euler_phi(m);
    ws = ([ vs[-1]-vs[-2] ] # sG
          + [ vs[0] ] + [ vs[_i+1] - vs[_i] for _i in range(n//2-1) ]);

    assert(len(ws) == n//2+1);
    ws = Matrix(ZZ, ws);
    return ws;


# Vals of [CDW21] "Mildly short vectors etc"
# ----------------------------------------------------------
def cdw21_va(m,a):
    return a*theta_kucera_a(m,1)-theta_kucera_a(m,a);

# v_a = a*theta(1) - theta(a)
# NB: v_1 = 0, v_m = m*theta(1)
def sti_vals_cdw21_v(m):
    ta = theta_kucera_all(m); assert (len(ta) == m//2+1);
    sG = ones_matrix(1,euler_phi(m))[0];
    vs = [ _a*ta[1] - (ta[_a] if _a <= m//2 else (sG-ta[m-_a]) if _a != m else 0)
           for _a in range(2,m+1) ];
    vs = Matrix(ZZ, vs);
    return vs;

# w_a = v_a - v_{a-1}, 2 <= a <= m
# w_m = v_m - v_{m-1} = s(G)
# w_a = w_{m-a+1} so we just keep s(G) + w_a for 2 <= a <= (m-1)//2 + 1 = ceil(m/2)
#     --> v[-1] - v[-2] = sG
#     --> w_2 = v_2 - v_1 = v[0]
#     --> w_{(m-1)//2 + 1} = v_{(m-1)//2+1}-v_{(m-1)//2} = v[(m-1)//2-1] - v[(m-1)//2-2]
def sti_vals_cdw21_w(m):
    vs = sti_vals_cdw21_v(m); assert (vs.nrows() == m-1);
    ws = ([ vs[-1]-vs[-2] ] # sG
          + [ vs[0] ] + [vs[_i+1] - vs[_i] for _i in range((m-1)//2-1)]);
    # ws = [ vs[0] ] + [vs[i+1] - vs[i] for i in range(m-2)];
    ws = Matrix(ZZ, ws);
    assert(ws.nrows() == ceil(m/2));
    return ws;


# Kucera's (large) Stickelberger basis
# -------------------------------------------------------
# Following Th.4.2: {omega*} + {omega(a) for a in M-}
# We don't use the "omega_a" function to save significant computational effort
def sti_vals_kucera(m, thetas=[]):
    Mminus = kucera_get_Mminus(m);

    # Special vectors (computed once)
    os = omega_kucera_star(m);
    sG = ones_matrix(1, len(os))[0];
    C  = os - sG;
    # Theta(a) is ta[a]
    ta = theta_kucera_all(m) if len(thetas) == 0 else thetas; assert(len(ta) == m//2+1);
    MK = matrix(ZZ, [os] + [ ta[_a] - _a*ta[1] + C if _a <= m/2 else - ta[m-_a] + (m-_a)*ta[1] - C
                             for _a in Mminus ]);
    return MK;


# Short Stickelberger basis (using [BK21] alpha_m function) 
# --------------------------------------------------------
def sti_vals_alpha(m, verbose=False):
    # Get all alpha triples
    t = cputime(); alphas = kucera_alpha_all(m); t = cputime(t);
    if verbose: print("\talpha(b):\tt={:.2f} [{:.2f}/b]".format(t, t/len(alphas))) 

    # Get all thetas (restrict only to necessary ones)
    alphas_idx = list({*sum(alphas,[])}); # list of unique indices in alphas (lot smaller than m)
    t = cputime(); ta = theta_kucera_idx(m, alphas_idx); t = cputime(t);
    if verbose: print("\tth(a):\tt={:.2f}".format(t)); assert (len(ta) == len(alphas_idx));
    
    # s(G): will be split as "real gens"
    sG  = ones_matrix(1,euler_phi(m))[0];
    Ma  = matrix(ZZ, [sG] + [ -sG + ta[alphas_idx.index(_a1)] + ta[alphas_idx.index(_a2)]
                              + ta[alphas_idx.index(_a3)] for _a1,_a2,_a3 in alphas ]);
    
    assert (Ma.height() == 1), "||M||>1"; # Nb: This does not verify that all M_ij > 0 (ie. M_ij != -1)
    return Ma;


# --------------------------------------------------------------------------------------------------
# Test tools for verifying the Stickelberger generators
# Computes the Stickelberger product idp^rel
def sti_id_product(K, p, rel, __idx_start=0):
    orb_p  = cf_orbit_p(K, p, __idx_start=__idx_start);
    id_sti = prod(map(pow, orb_p, rel));
    return id_sti;


# --------------------------------------------------------------------------------------------------
# Stickelberger's generators

# Following [Was97, §6] "i-s_i" relations
# ---------------------------------------------------------
# For relations "(i-s_i)\theta", following [Was97,§6], and especially proof of Th.6.10
# These correspond to cdw17_v

# Returns h = x^{mp}-1, phi_{p}(x^m), phi_{m}(x^p)
def __get_sti_moduli_x(m,p):
    assert (is_prime(p));
    
    Zy = PolynomialRing(ZZ, name='y', sparse=True);  y = Zy.gen();
    Zx = PolynomialRing(ZZ, name='x', sparse=False); x = Zx.gen();
    
    h      = x**(m*p)-1;
    # phip_m = Zx( cyclotomic_polynomial(p)(y^m) ); # Forever, but phi_p = sum x^i, i<p
    phip_m = Zx( { _k*m: 1 for _k in range(p) } );
    phim_p = Zx( cyclotomic_polynomial(m)(y**p) );

    return (h, phip_m, phim_p);


# Stickelberger elt is gamma^_a / \sigma_a(gamma), galois[_a]: zm -> zm^_a
# Use [Was97, pf.Lem.6.4] to obtain \sigma_a(g(chi)) = g( chi^a )
# Use [Was97, Lem.6.1]    to obtain 1/g(chi)  = bar(g(chi))/p
#                         to obtain g(chibar) = chi(-1) bar(g(chi)) # Useless
# ---> Compute for each a: g(chi)^{a-\sigma_a} = g(chi)^a * g(chibar^a) / q
# ---> chi is chi_p0 such that chi(a) = a mod p0
def sti_gens_cdw17_v(K, p0):
    m   = cf_get_conductor(K);
    chi = get_chi(K, m, p0); chi = chi_conj(chi); # [Was97] is the conjugate of the residue character 
    Zx  = PolynomialRing(ZZ, name='x');
    Fq  = chi[0].parent(); q = Fq.cardinality(); p = Fq.characteristic();
    # Compute the polynomial moduli used in the computation
    t = cputime(); (h, phip_m, phim_p) = __get_sti_moduli_x(m,p); t = cputime(t);
    print("h/phip_m/phim_p:\t{:.4f}s".format(t), flush=True);

    # g(chi) is re-used each time
    t = cputime(); g_x = gauss_sum_x(chi); t = cputime(t);
    print("g(chi):\t\t\t{:.4f}s".format(t), flush=True);
    
    print("Sti gens v_i = (i-s_i).theta(-1):", flush=True);
    stirel_idx = range(2, m//2+1); # [_a for _a in range(2,m+2) if gcd(_a,m) == 1];
    sti_gens   = [];
    __prev = 1; g_x_pow_a  = g_x;
    for _a in stirel_idx:
        t = cputime();
        # NB: power_mod(g_x, _a, h) is also a viable option, but less efficient when we want all powers
        # g_x_pow_a = power_mod(g_x, _a, h);
        g_x_pow_a = (g_x_pow_a * power_mod(g_x, _a-__prev, h)).mod(h);
        # chi^-a is conj(chi)^a
        # ga_x is hence q / g(chi)^{s_a} (. chi(-1))
        ga_x      = gauss_sum_chipow_x(chi, -_a);
        # This is g(chi)^{a-s_a}
        # NB: factor q (from ga_x) is removable only after all moduli have been applied
        sti_x     = (g_x_pow_a * ga_x).mod(h).mod(phip_m).mod(phim_p) / q;
        # All non zero monomials in sti_x are of degree p*... [ie, sti_x is in Q(\zeta_m) = Q(x^p)]
        sti_gen   = K(Zx( {_key/p: _val for (_key,_val) in sti_x.dict().items()} ));
        t = cputime(t);
        print("\t{}-s_{}\tt={:.4f}s".format(_a,_a,t), flush=True);
        sti_gens.append(sti_gen);
        __prev = _a;

    return sti_gens;


# Computing directly w_i = v_i / v_{i-1} to avoid integers overflows
# [This way, g(w_i)*bar(g(w_i)) = q , vs "q^i" for v_i]
# Using g(chi)^{d-c} / s_d * s_c = g(chi^c) * g(chi^-d) / q
def sti_gens_cdw17_w(K, p0):
    m   = cf_get_conductor(K);
    n   = euler_phi(m);
    chi = get_chi(K, m, p0); chi = chi_conj(chi); # [Was97] is the conjugate of the residue character 
    Zx  = PolynomialRing(Integers(), name='x');
    Fq  = chi[0].parent(); q = Fq.cardinality(); p = Fq.characteristic();
    # Compute the polynomial moduli used in the computation
    t = cputime(); (h, phip_m, phim_p) = __get_sti_moduli_x(m,p); t = cputime(t);
    print("\th/phip_m/phim_p:\tt={:.4f}s".format(t), flush=True);

    # g(chi) is re-used each time
    t = cputime(); g_x = gauss_sum_x(chi); t = cputime(t);
    print("\tg(chi):\t\t\tt={:.4f}s".format(t), flush=True);

    print("\tSti gens w_i = v_i - v_{i-1}:", flush=True);
    # stirel_idx = range(2, m//2+2); # Gives "all" generators, but not the CDW17 one (which can "jump")
    stirel_idx = [_a for _a in range(2,m+2) if gcd(_a,m) == 1][:n//2];
    sti_gens_w = [K(q)]; # sti_gens_w[i] will be conjugate of gen of omega(i)-omega(i+1), [0] is "m+1"
    __prev = 0; ga_prev_x = Zx(1);
    for _a in stirel_idx:
        # See comments of get_sti_gens_vx()
        t = cputime();
        g_x_pow_diff = power_mod(g_x, _a-__prev, h);
        # Compute g(chi^a)
        chi_pow_a = chi_pow(chi, _a);
        ga_x      = gauss_sum_x(chi_pow_a);
         # g(chi^(-a)) = q / s_a(g(chi))
        gca_x     = gauss_sum_chiconj_x(chi_pow_a, g_chi=ga_x);
        # This is g(chi)^{s_prev-s_a} * q
        g_sig_x   = (ga_prev_x*gca_x).mod(h);
        # g(chi)^{a-pr} / s_a(g(chi)) * s_pr(g(chi)) = g(chi)^{a-pr} * g(chi^pr) * g(chi^(-a)) / q
        # NB: factor q (from gca_x) is removable only after all moduli have been applied
        sti_x     = (g_x_pow_diff * g_sig_x).mod(h).mod(phip_m).mod(phim_p) / q;
        # All non zero monomials in sti_x are of degree p*... [ie, sti_x is in Q(\zeta_m) = Q(x^p)]
        sti_gen   = K(Zx( {_key/p: _val for (_key,_val) in sti_x.dict().items()} ));
        t = cputime(t);
        print("\t1-s_{}+s_{}\tt={:.4f}s".format(_a,_a-1,t), flush=True);
        sti_gens_w.append(sti_gen);
        __prev = _a; ga_prev_x = ga_x;

    return sti_gens_w;



# Using [BK21], generators for short relations write nicely as Jacobi sums
# ---------------------------------------------------------

# Returns generators for set of { [a,b,c] } indexes that correspond to the Stickelberger elt:
#         \th(a) + \th(b) + \th(c) - sG , where m | a+b+c
# Using Jacobi sums
# S is a list of [a_i,b_i,c_i] that verify (a_i+b_i+c_i)=0[m]
def sti_gens_alpha_idx(K, p0, S):
    m = cf_get_conductor(K);

    # Residue norm character
    t = cputime(); res_chi = get_chi(K, m, p0); t = cputime(t);
    Fq_a = res_chi[0]; q = Fq_a.parent().cardinality();
    print("\tchi(p0)\tt={:.4f}s".format(t));
    # Computing DLog tables for Jacobi sums
    t = cputime(); dl = jacobi_sum_dlog_table(Fq_a); t = cputime(t);
    print("\tdlog [{}]\tt={:.4f}s".format(q, t));
    # Generators as Jacobi sums J(i[0],i[1])
    sti_gens_al = [K(q)];
    for _a,_b,_c in S:
        assert(mod(_a+_b+_c,m) == 0);
        _chi_a = chi_pow(res_chi, _a);
        _chi_b = chi_pow(res_chi, _b);

        t = cputime(); _j_ab  = jacobi_sum_x(_chi_a, _chi_b, dlog=dl); t = cputime(t);
        _j_ab  = K(_j_ab);
        print("\t[{}+{}+{}]\tt={:.4f}s".format(_a,_b,_c,t), flush=True);
        sti_gens_al.append(_j_ab);
        
    return sti_gens_al;


# For the whole basis
def sti_gens_alpha(K, p0):
    m  = cf_get_conductor(K);
    S  = kucera_alpha_all(m);
    sk = sti_gens_alpha_idx(K, p0, S);
    return sk;
    

# This induces a simplified computation for the cdw17/cdw21_w elements
# Compute the [CDW21] generators (actually the conjugate version of [CDW17]), but using Jacobi sums
def sti_gens_cdw21_w(K, p0):
    m       = cf_get_conductor(K);
    jac_idx = [ [1,_a-1,m-_a] for _a in range(2,(m-1)//2+2) ];
    sw      = sti_gens_alpha_idx(K, p0, jac_idx);
    return sw;


# //-- END OF FILE
