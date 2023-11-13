from galois import GF2
from functools import reduce
from numpy.linalg import matrix_rank


def in_basis(B, e):
    V = list(B)
    U = set(reduce(frozenset.union, V))
    if len(e - U) > 0:
        return False
    U = list(U)
    M = GF2.Zeros((len(U), len(V)+1))
    for i in range(len(U)):
        if U[i] in e:
            M[i][len(V)] = 1
        for j in range(len(V)):
            if U[i] in V[j]:
                M[i][j] = 1
    if matrix_rank(M) > len(V):
        return False
    return True

def find_complementary_space(O, B):
    # assume U contains B
    V = list(B) + list(O-B)
    U = set(reduce(frozenset.union, B))
    U = list(U) + list(set(reduce(frozenset.union, O)) - U)
    M = GF2.Zeros((len(V), len(U)))
    for i in range(len(V)):
        for j in range(len(U)):
            if U[j] in V[i]:
                M[i][j] = 1
    M.row_reduce()
    M = M[len(B):, :]
    M.row_reduce()
    Bn = set()
    for i in range(len(V)-len(B)):
        basis_element = set()
        for j in range(len(U)):
            if M[i][j] == 1:
                basis_element.add(U[j])
        if len(basis_element) == 0:
            break
        Bn.add(frozenset(basis_element))
    return Bn


def merge_bases(Bs):
    V = list(reduce(set.union, Bs))
    U = list(reduce(frozenset.union, V))
    M = GF2.Zeros((len(U), len(V)))
    for i in range(len(U)):
        for j in range(len(V)):
            if U[i] in V[j]:
                M[i][j] = 1
    M = GF2.Zeros((len(V), len(U)))
    for i in range(len(V)):
        for j in range(len(U)):
            if U[j] in V[i]:
                M[i][j] = 1
    M.row_reduce()
    
    Bn = set()
    for i in range(len(V)):
        basis_element = set()
        for j in range(len(U)):
            if M[i][j] == 1:
                basis_element.add(U[j])
        if len(basis_element) == 0:
            break
        Bn.add(frozenset(basis_element))
    return Bn


def intersect_bases(B1, B2):
    V1 = list(B1)
    V2 = list(B2)
    U = list(reduce(frozenset.union, B1 | B2))
    M = GF2.Zeros((len(U), len(V1) + len(V2)))
    for i in range(len(U)):
        for j in range(len(V1)):
            if U[i] in V1[j]:
                M[i][j] = 1
        for j in range(len(V2)):
            if U[i] in V2[j]:
                M[i][j+len(V1)] = 1
    N = M.null_space().T
    I = (M[:, :len(V1)] @ N[:len(V1), :]).column_space().T
    Bn = set()
    for j in range(I.shape[1]):
        basis_element = set()
        for i in range(len(U)):
            if I[i][j] == 1:
                basis_element.add(U[i])
        Bn.add(frozenset(basis_element))
    return Bn
