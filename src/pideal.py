# Code in support of ePrint:2021/xxxx
# Copyright 2021, Olivier Bernard
# GPL-3.0-only (see LICENSE file)

from sage.all import *

import fp
from ZR_mat import *        # Row-style lower triangular HNF



# --------------------------------------------------------------------------------------------------
# Prime ideals

# In Sagemaths, several ideal functions rely on their HNF on some basis order (as they should)
# It seems that Sage always try to compute the maximal order (which it shouldn't).
# This is the case for (e.g.):
#         .norm(), .basis(), .gens_two(), ...
# and this is quite annoying.
#
# I didn't find a way to tell Sage what is the maximal order, or to obtain a less poor behaviour
# for number field ideals, so these are fast (but maybe generically incorrect) methods that should
# work *in our setting*.

# Should be correct in most common situations *FOR PRIME IDEALS*
def pid_fast_gens_two(pid):
    (p, g) = pid.gens();# [:2]; # Removing [:2] should guarantee the 2-elt representation is correct  
    p = ZZ(p); g = g.polynomial();
    assert(is_prime(p)); # Keep it for now but unnecessary if the [:2] trick works 
    return (p, g);


def pid_fast_norm(pid):
    (p, g) = pid_fast_gens_two(pid);
    f      = g.degree();
    return ZZ(p**f);


def pid_fast_smallest_integer(pid):
    (p, g) = pid_fast_gens_two(pid);
    return p;


def pid_fast_residue_class_degree(pid):
    (p, g) = pid_fast_gens_two(pid);
    f      = g.degree();
    return ZZ(f);


# //-- END OF FILE
