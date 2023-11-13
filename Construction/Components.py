from .Function import Function 
from Modelling.PropModels import compute_parity_propagation_model
from abc import abstractmethod
from pysat.card import CardEnc, EncType

class Component(Function):
    def __init__(self, input_size, output_size):
        super().__init__(input_size, output_size)
        self.parity_prop_models = {}

    def compute_parity_propagation_model(self, input_key_mask, output_key_mask):
        return compute_parity_propagation_model(self, input_key_mask, output_key_mask)

    def get_parity_propagation_model(self, input_key_mask, output_key_mask):
        if (input_key_mask, output_key_mask) not in self.parity_prop_models:
            self.parity_prop_models[(input_key_mask, output_key_mask)] = self.compute_parity_propagation_model(input_key_mask, output_key_mask)
        return self.parity_prop_models[(input_key_mask, output_key_mask)]

    @abstractmethod
    def __call__(self, v):
        return 0

class SBox(Component):

    def __init__(self, input_size, output_size, lookup_table: list):
        self.lookup_table = lookup_table
        super().__init__(input_size, output_size)

    def __call__(self, v):
        return self.lookup_table[v]

class __COPYn(Component):
    def __init__(self, n):
        self.n = n
        super().__init__(1, n)

    def __call__(self, v):
        r = 0
        for _ in range(self.n):
            r <<= 1
            r += v & 1
        return r

_COPYs = {}
def get_COPYn(i):
    global _COPYs
    if i not in _COPYs:
        _COPYs[i] = __COPYn(i)
    return _COPYs[i]

class __XORn(Component):

    def __init__(self, n):
        super().__init__(n, 1)
        self.n = n

    def __call__(self, v):
        r = 0
        for i in range(self.n):
            r ^= v & 1
            v >>= 1
        return r
    
    def compute_parity_propagation_model(self, input_key_mask, output_key_mask):
        atmost1 = tuple(tuple(x) for x in CardEnc.atmost(lits = list(range(1, self.n+1)), encoding = EncType.pairwise).clauses)
        if input_key_mask.bit_count() == output_key_mask.bit_count() == 0:
            return atmost1 +  tuple((-v, self.n+1) for v in range(1, self.n+1)) + (tuple(range(1, self.n+1)) + (-self.n+1,),)
        else:
            return atmost1 + tuple((-v, self.n+1) for v in range(1, self.n+1))


_XORs = {}
def get_XORn(i):
    global _XORs
    if i not in _XORs.keys():
        _XORs[i] = __XORn(i)
    return _XORs[i]

XOR = get_XORn(2)

class __Sink(Component):

    def __init__(self):
        super().__init__(1, 0)

    def __call__(self, v):
        return 0

    def compute_parity_propagation_model(self, input_key_mask, output_key_mask):
        return ((-1,),)

Sink = __Sink()

class __IDn(Component):

    def __init__(self, n):
        super().__init__(n, n)
        self.n = n

    def __call__(self, v):
        return v
    
    def compute_parity_propagation_model(self, input_key_mask, output_key_mask):
        clauses = []
        km = input_key_mask | output_key_mask
        for i in range(self.n):
            clauses.append((-i-1, i+1+self.n))
            if km & (1<<i) == 0:
                clauses.append((i+1, -i-1-self.n))
        return tuple(clauses)
    
_IDs = {}
def get_IDn(i):
    global _IDs
    if i not in _IDs.keys():
        _IDs[i] = __IDn(i)
    return _IDs[i]

ID = get_XORn(1)