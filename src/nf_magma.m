/* 
 * Code in support of ePrint:2021/xxxx
 * Copyright 2021, Olivier Bernard
 * GPL-3.0-only (see LICENSE file)
 */
SetClassGroupBounds("GRH");

/* -------------------------------------------------------------------------------------------------- */
/* ----- S-UNITS COMPUTATIONS ----- */

/* ----- S-Units wrpt. some FB ----- */
function get_raw_S_units(K, FB)
    _r1, _r2          := Signature(K);
    _G_SU, _mU, _B_SU := SUnitGroup(FB : Raw:=true);
    assert(NumberOfGenerators(_G_SU) eq _r1 + _r2 + #FB);
    
    // Remove torsion (check it is indeed _G_SU.1)
    _tors_K := K! PowerProduct(_B_SU, _mU(_G_SU.1));
    assert (_tors_K eq _B_SU[1]); // Note that _B_SU[1] *is* _tors_K
    assert (AbsoluteValue(Degree(K)-Length(_tors_K)) lt 2^(-100));
    // First "generators" are integers that are never used anywhere.
    _null := 0; // _null is the number of useless first coordinates (all 0)
    repeat
        _null := _null+1;
    until &or[ _mU(_G_SU.i)[_null] ne 0: i in [2..NumberOfGenerators(_G_SU)]];
    _null := _null-1; // _null should be at least 1 (torsion)
    
    // Adapt the basis accordingly
    _B_SU := [K!_b_su : _b_su in Eltseq(_B_SU)[[_null+1..Degree(_B_SU)]]];
    _y_U  := [ Eltseq(_mU(_G_SU.i))[[_null+1..#_B_SU+_null]]: i in [2.._r1+_r2] ]; // Units are 2..r1+r2
    assert (&and[ AbsoluteValue(Norm(_u)) eq 1: _u in [K| PowerProduct(_B_SU, _u): _u in _y_U] ]);
    _y_SU := [ Eltseq(_mU(_G_SU.i))[_null+1..#_B_SU+_null]: i in [_r1+_r2+1..NumberOfGenerators(_G_SU)] ];
    // Valuations at FB
    t := Cputime(); _B_vp := [ [ Valuation(_su, _pid): _pid in FB ]: _su in _B_SU ]; t := Cputime(t);
    if t gt #_B_SU then
        printf "/!\ Valuations of S-units over FB t=%o [%o/su]\n", t, t/#_B_SU;
    end if;

    return _y_U, _y_SU, _B_SU, _B_vp;
end function;


/* -------------------------------------------------------------------------------------------------- */
/* INPUT/OUTPUT 
   Format SU, FB, NF :--> see Sage "pcmp_io.py" */

/* ----------------------------------------- */
/* Number Field */
/* in: 'z23' for Cyclo(23) */
function nf_set_tag(tag)
    Zx<x>:=PolynomialRing(Integers());
    typ  := tag[1]; assert (typ eq "z");
    m    := StringToInteger(tag[2..#tag]);
    K    := CyclotomicField(m);
    return K;
end function;

function __nf_out_str(K, tag)
    nf_str := Sprintf("# number_field: tag=%o eq=%o", tag,
                      &cat[ss: ss in Split(Sprintf("%o",Eltseq(DefiningPolynomial(K))), " ")]);
    return nf_str;
end function;


/* ----------------------------------------- */
/* Prime ideals: for FB and Challenges */
function __fb_pid_in_str(line, K)
    g2 := Split(line, " \n"); assert(#g2 eq 2);
    p  := StringToInteger(g2[1]);
    a  := Evaluate(Polynomial([StringToInteger(_a): _a in Split(g2[2], "[,]")]), K.1);
    pid := ideal< MaximalOrder(K)| p, a>;
    return pid;
end function;


function fb_in_stream(filename, K)
    F := Open(filename, "r");
    // ignore first line
    _ := Gets(F);
    // read second line to find k
    fb_head := Gets(F); assert(not IsEof(fb_head));
    exp     := "# fb_places: k=";
    assert(fb_head[1..#exp] eq exp);
    k       := StringToInteger(fb_head[#exp+1..#fb_head]);
    
    // read each line and convert
    FB := [];
    for _i in [1..k] do
        _line:= Gets(F); assert(not IsEof(_line));
        _pid := __fb_pid_in_str(_line, K);
        Append(~FB, _pid);
    end for;
    
    delete F;
    return FB;
end function;


/* ----------------------------------------- */
/* Sunits */
__END_TAG := "# --- END ---";

// NB: returns a string
function __su_out_str(K, su)
    _su_str := &cat[_su_c: _su_c in Split(Sprintf("%o",Eltseq(Polynomial(Eltseq(su))))," ")];
    return _su_str;
end function;

function __raw_power_out_str(y)
    _y_str  := &cat[_y_c: _y_c in Split(Sprintf("%o",Eltseq(y))," ")];
    return _y_str;
end function;

function __valp_out_str(vp)
    _vp_str := &cat[_vp_c: _vp_c in Split(Sprintf("%o",Eltseq(vp))," ")];
    return _vp_str;
end function;


procedure sunits_raw_write_data(filename, K, tag, y_u, y_su, B_su, B_vp)
    r1, r2 := Signature(K);
    assert (#y_u  eq r1+r2-1);
    assert (#B_su eq #B_vp);
    assert (&and[ #_yi eq #B_su: _yi in y_u cat y_su ]);
    assert (&and[ #_vp eq #y_su: _vp in B_vp ]);
    
    F := Open(filename, "w");
    nf_str := __nf_out_str(K, tag);

    fprintf F, nf_str cat "\n";
    fprintf F, "# raw sunits(fb): nu=%o k=%o N=%o\n", #y_u, #y_su, #B_su;
    for _yi in y_u cat y_su do
        y_str := __raw_power_out_str(_yi);
        fprintf F, y_str cat "\n";
    end for;
    for _su in B_su do
        su_str := __su_out_str(K, _su);
        fprintf F, su_str cat "\n";
    end for;
    for _vp in B_vp do
        vp_str := __valp_out_str(_vp);
        fprintf F, vp_str cat "\n";
    end for;
    fprintf F, __END_TAG cat "\n";
    
    delete F;
end procedure;


//-- END OF FILE
