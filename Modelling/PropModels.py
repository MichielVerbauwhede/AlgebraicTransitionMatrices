from LogicOptimisation.QMC import QMC_optimise_CNF
from bitarrays.bitarray import BitArray
from bitarrays.bitset import BitSet

def transition_matrix(f, input_size, output_size):
    res = BitArray((2**output_size, 2**input_size))
    for i in range(2**input_size):
        res.set((f(i), i))
    return res

def ANF_matrix(f, input_size, output_size):
    res = transition_matrix(f, input_size, output_size)
    res.XORdown(0, -1)
    res.XORup(1, -1)
    return res

def ANF_prop_table(F, input_mask, output_mask):
    res = ANF_matrix(lambda x: F(x), F.input_size, F.output_size)
    res.ORup(0, output_mask)
    res.ORdown(1, input_mask)
    return res

def compute_parity_propagation_model(F, input_mask, output_mask):
    m = ANF_prop_table(F, input_mask, output_mask)
    dont_care = BitSet(2**(F.input_size + F.output_size))
    return QMC_optimise_CNF(m, dont_care)