from functools import partial, reduce
from math import inf
from Tools.IntTools import pext
from bitarrays.bitset import BitSet


def is_classic(e):
    # evaluates whether this element could found by parity sets, 3SDP etc
    # This is basically, a single parity set to a single monomial
    return len(e) == 1


def reduce_element(e):
    # we iteratively apply the following reduction rule:
    # 1. all elements (x, y) with the same x become (x, \Delta y)
    # 2. all elements (x, y) with the same y become (\delta x, y)
    l = inf
    s = {(frozenset({x[0], }), frozenset({x[1], })) for x in e}
    while len(s) < l:
        l = len(s)
        snew = set()
        while len(s):
            x, y = s.pop()
            for x2, y2 in list(s):
                if x == x2:
                    y = y.symmetric_difference(y2)
                    s.discard((x2, y2))
            snew.add((x, y))
        s = snew
        snew = set()
        while len(s):
            x, y = s.pop()
            for x2, y2 in list(s):
                if y == y2:
                    x = x.symmetric_difference(x2)
                    s.discard((x2, y2))
            snew.add((x, y))
        s = snew
    return s


def is_input_independent(e):
    # check whether the element is of the form r(F(X))
    return len(reduce_element(e)) == 1


def is_linearly_equivalent(e, invmask):
    # evaluate whether this element is the result of
    # a single parity set to a single monomial, but of a linearly equivalent cipher
    def evaluate(vs):
        # linearity is evaluated by counting the number of times function evaluates to one
        # a product of d linear functions in n variables should have 2^{n-d} ones
        d = max(x.bit_count() for x in vs)
        m = reduce(int.__or__, vs)
        n = m.bit_count()
        s = BitSet(2**n)
        for x in vs:
            s.set(pext(x, m))
        s.XORup()
        return s.count() == 2**(n-d)

    rede = reduce_element(e)
    if len(rede) > 1:
        return False
    rede = rede.pop()

    if not all((x ^ invmask).bit_count() == 1 for x in rede[0]):
        # get all involved bits
        if not evaluate([x ^ invmask for x in rede[0]]):
            return False
    if not all(x.bit_count() == 1 for x in rede[1]):
        if not evaluate(rede[1]):
            return False
    return True

def sort_basis(B, invmask):
    # sort into 2 buckets, input dependent or not
    # then sort the second one into two more buckets, linearly equivalent or not
    # sort the last one into two bucket again, classis or linearly equivalen
    input_independent = set(filter(is_input_independent, B))
    input_dependent = B - input_independent
    linear_equiv = set(filter(partial(is_linearly_equivalent, invmask = invmask), input_independent))
    nonlinear_equiv = input_independent - linear_equiv
    classic = set(filter(is_classic, linear_equiv))
    return classic, linear_equiv - classic, nonlinear_equiv, input_dependent

def nice_print(B, ninputs, noutputs):
    classic, lin, nonlin, inputdep = sort_basis(B, 2**ninputs-1)
    print("classic:", len(classic))
    classic = sorted(list(x)[0] for x in classic)
    for e in classic:
        print(bin(e[0])[2:].zfill(ninputs), "->", bin(e[1])[2:].zfill(noutputs))

    print()
    print("linearly equivalent:", len(lin))
    lin = sorted(sorted(x) for x in lin)
    for c in lin:
        for e in c:
            print(bin(e[0])[2:].zfill(ninputs), "->", bin(e[1])[2:].zfill(noutputs))
        print()

    print()
    print("nonlinearly equivalent:", len(nonlin))
    nonlin = sorted(sorted(x) for x in nonlin)
    for c in nonlin:
        for e in c:
            print(bin(e[0])[2:].zfill(ninputs), "->", bin(e[1])[2:].zfill(noutputs))
        print()

    print()
    print("input dependent:", len(inputdep))
    inputdep = sorted(sorted(x) for x in inputdep)
    for c in inputdep:
        for e in c:
            print(bin(e[0])[2:].zfill(ninputs), "->", bin(e[1])[2:].zfill(noutputs))
        print()
    
def stats(B, ninputs, noutputs):
    classic, lin, nonlin, inputdep = sort_basis(B, 2**ninputs-1)
    print("Total:", len(B))
    print("classic:", len(classic))
    print("linearly equivalent:", len(lin))
    print("nonlinearly equivalent:", len(nonlin))
    print("input dependent:", len(inputdep))