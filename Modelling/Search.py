from functools import reduce
from galois import GF2
from multiprocessing import Pool
from itertools import product

def search_integral_properties(Avec, input_size, output_size, is_permutation, num_threads=1):
    # Avec gets an io pair (u, v) (ints) and returns a boolean and a set of pairs (a, b).
    # The boolean indicates keydepence in the product $(A^{F_1}\otimes{A^{F_3}}^\intercal)_{(a, b), (u, v)}\mathcal{A}(F_2, a, b)$
    # if it is key independent, the set contains the support of $a, b \mapsto (A^{F_1}\otimes{A^{F_3}}^\intercal)_{(a, b), (u, v)}\mathcal{A}(F_2, a, b)$
    invmask = 2**input_size - 1
    input_modifiers = tuple((1<<x) ^ invmask for x in range(input_size))
    output_modifiers = tuple(1<<x for x in range(output_size))
    basis_weight_1 = set()
    # sort candidates into B (already in nullspace), UV(to be used for the nullspace computation) or discard
    Bi = [set()]
    WUV = {}
    with Pool(num_threads) as pool:
        while len(Bi) == 1 or len(Bi[-1]) > 0:
            if len(Bi) == 1 and not is_permutation:
                to_test = [(invmask, x) for x in output_modifiers]
            elif len(Bi) == 1 and is_permutation:
                Bi.append(set())
                to_test = list(product(input_modifiers, output_modifiers))
            else:
                to_test = set()
                for u1, v1 in Bi[-1]:
                    for u2 in input_modifiers:
                        u = u1 & u2
                        if (u ^ invmask).bit_count() + v1.bit_count() == len(Bi):
                            flag = True
                            for u3 in input_modifiers:
                                if ((u | u3) ^ invmask) > 0 and not (u^u3^invmask, v1) in Bi[-1]:
                                    flag = False
                                    break
                            if flag:
                                to_test.add((u, v1))
                    for v2 in output_modifiers:
                        v = v1 | v2
                        if (u1 ^ invmask).bit_count() + v.bit_count() == len(Bi):
                            flag = True
                            for v3 in output_modifiers:
                                if (v & v3) > 0 and not (u1, v^v3) in Bi[-1]:
                                    flag = False
                                    break
                            if flag:
                                to_test.add((u1, v))
                to_test = list(to_test)
            Bi.append(set())
            for uv, (kd, W) in zip(to_test, pool.map(Avec, to_test, 1)):
                if not kd:
                    if len(W) == 0:
                        Bi[-1].add(uv)
                    else:
                        WUV[uv] = W

    for uv in reduce(set.union, Bi, set()):
        basis_weight_1.add(frozenset({uv, }))

    if len(WUV) == 0:
        return basis_weight_1
    
    # find properties
    UV = list(reduce(set.union, WUV.values()))
    WUV = list(WUV.items())
    M = GF2.Zeros((len(UV), len(WUV)))
    for i in range(len(UV)):
        for j in range(len(WUV)):
            if UV[i] in WUV[j][1]:
                M[i][j] = 1
    K = M.null_space().T 

    rest_basis = set()
    for j in range(K.shape[1]):
        basis_element = set() 
        for i in range(len(WUV)):
            if K[i][j]==1:
                basis_element.add(WUV[i][0])
        rest_basis.add(frozenset(basis_element))


    return basis_weight_1 | rest_basis