from .Function import Function
from .Components import Component, XOR
from pysat.formula import IDPool

INPUT_ID = 0
OUTPUT_ID = 2**64
KEY_ID = -1
ERROR_ID = -2


class Null(Function):
    def __init__(self, input_size, output_size):
        super().__init__(input_size, output_size)
        self.n_cnf_vars = 0

    def __call__(self, v):
        return 0

class component_record(object):
    def __init__(self, f):
        self.f = f
        self.input_connections = [(ERROR_ID, ERROR_ID)]*f.input_size
        self.output_connections = [(ERROR_ID, ERROR_ID)]*f.output_size
        self.input_key_mask = 0
        self.output_key_mask = 0
        return


class CompoundFunction(Function):
    def __init__(self, input_size, output_size):
        super().__init__(input_size, output_size)
        # component list contains component, inputs, outputs, input_vars, output_vars
        self.components =\
            {INPUT_ID: component_record(Null(0, self.input_size)) ,
             OUTPUT_ID: component_record(Null(self.output_size, 0))}
        self.n_components = 0
        self.n_key_bits = 0
        self.key = 0
        self.local_model = tuple()
        self.local_input_vars = tuple()
        self.local_output_vars = tuple()
        self.local_key_vars = tuple()
        self.n_vars = 0
        return

    def add_component(self, component):
        """
        Add component to Compound Function.
        Components should be added in order of execution.
        The id of the component is returned
        """
        assert(not isinstance(component, CompoundFunction) or component.n_key_bits == 0)
        self.n_components += 1
        self.components[self.n_components] = component_record(component)
        return self.n_components

    def connect_components(self, from_component, from_wire, to_component, to_wire):
        """
        connect one component to another
        """
        # check existence of components
        assert(-1 < from_component <= self.n_components)
        assert(-1 < to_component <= self.n_components or to_component == OUTPUT_ID)
        # check existence of wires
        assert(from_wire < self.components[from_component].f.output_size)
        assert(to_wire < self.components[to_component].f.input_size)
        # enforce order of computation
        assert(from_component < to_component)
        # update information
        assert(self.components[from_component].output_connections[from_wire] == (ERROR_ID, ERROR_ID))
        self.components[from_component].output_connections[from_wire] = (to_component, to_wire)
        assert(self.components[to_component].input_connections[to_wire] == (ERROR_ID, ERROR_ID))
        self.components[to_component].input_connections[to_wire] = (from_component, from_wire)
        return

    def connect_to_key(self, to_component, to_wire):
        assert(-1 < to_component <= self.n_components)
        assert(to_wire < self.components[to_component].f.input_size)
        assert(self.components[to_component].input_connections[to_wire] == (ERROR_ID, ERROR_ID))
        self.components[to_component].input_connections[to_wire] = (KEY_ID, self.n_key_bits)
        self.n_key_bits += 1
        return

    def __call__(self, v):
        # evaluates with random key based on the given key value
        outputs = [v]
        for i in range(1, self.n_components+1):
            # construct input
            x = 0
            for wire in self.components[i].input_connections[::-1]:
                x <<= 1
                if wire[0] == KEY_ID:
                    x |= ((key >> wire[1]) & 1)
                else:
                    x |= ((outputs[wire[0]] >> wire[1]) & 1)
            # get output
            outputs.append(self.components.get(i).f(x))
        # construct output
        x = 0
        for wire in self.components[OUTPUT_ID].input_connections[::-1]:
            x <<= 1
            if wire[0] == KEY_ID:
                x |= ((key >> wire[1]) & 1)
            else:
                x |= ((outputs[wire[0]] >> wire[1]) & 1)
        return x

    def __build_local_model(self):
        for record in self.components.values():
            assert((ERROR_ID, ERROR_ID) not in record.input_connections)
            assert((ERROR_ID, ERROR_ID) not in record.output_connections)

        pool = IDPool()
        self.local_model = []
        # build model
        all_component_output_vars = {}
        self.local_key_vars = tuple(pool.id() for _ in range(self.n_key_bits))
        for i in range(self.n_key_bits):
            all_component_output_vars[(KEY_ID, i)] = self.local_key_vars[i]
        for i in range(self.n_components+1):
            component = self.components[i]
            f = component.f
            # generate output vars
            output_vars = tuple(pool.id() for _ in range(f.output_size))
            for j in range(f.output_size):
                all_component_output_vars[(i, j)] = output_vars[j]
            # generate input vars
            input_vars = tuple(all_component_output_vars[wire] for wire in component.input_connections)
            # generate clauses
            if isinstance(f, Component):
                vs = [0]*(f.input_size + f.output_size+1)
                for j in range(f.input_size):
                    vs[j+1] = input_vars[j]
                for j in range(f.output_size):
                    vs[f.input_size + j + 1] = output_vars[j]
                vs += [-x for x in vs[1:][::-1]]
                clauses = f.get_parity_propagation_model(component.input_key_mask, component.output_key_mask)
                for clause in f.get_parity_propagation_model(component.input_key_mask, component.output_key_mask):
                    self.local_model.append(tuple(vs[x] for x in clause))
            elif isinstance(f, CompoundFunction):
                self.local_model += f.to_model(pool, input_vars, output_vars)[0]
            else:
                self.local_input_vars = output_vars
            # recover local output vars
        self.local_output_vars = tuple(all_component_output_vars[wire] for wire in self.components[OUTPUT_ID].input_connections)
        self.n_vars = max(map(max, self.local_model))
        return

    def to_model(self, pool=None, input_vars = tuple(), output_vars = tuple(), key_vars = tuple()):
        if len(self.local_model) == 0:
            self.__build_local_model()

        if pool is None:
            return tuple(self.local_model), self.local_input_vars, self.local_output_vars, self.local_key_vars
        
        if len(input_vars) == 0:
            input_vars = tuple(pool.id() for _ in range(self.input_size))
        else:
            input_vars = tuple(input_vars)
        if len(output_vars) == 0:
            output_vars = tuple(pool.id() for _ in range(self.output_size))
        else:
            output_vars = tuple(output_vars)

        vs = [0]*(self.n_vars+1)
        for i in range(self.input_size):
            vs[self.local_input_vars[i]] = input_vars[i]
        for i in range(self.output_size):
            vs[self.local_output_vars[i]] = output_vars[i]
        for i in range(1, len(vs)):
            if vs[i] == 0:
                vs[i] = pool.id()
        vs += [-x for x in vs[1:][::-1]]
        model = []
        for clause in self.local_model:
            model.append(tuple(vs[x] for x in clause))
        
        return model, input_vars, output_vars, key_vars

    def __component_is_key_addition(self, c):
        return (c.f is XOR) and (c.input_connections[0][0] == KEY_ID or c.input_connections[1][0] == KEY_ID)

    def optimized_for_nonzero_trail_detection(self, global_input_key_mask=0, global_output_key_mask=0):
        f = CompoundFunction(self.input_size, self.output_size)
        o2n = {INPUT_ID:INPUT_ID}
        for i in [x for x in range(1, self.n_components+1)]:
            component = self.components[i]
            if not self.__component_is_key_addition(component):
                o2n[i] = f.add_component(component.f)
                for j in range(component.f.input_size):
                    fid, fw = component.input_connections[j]
                    if fid >= 0:
                        fc = self.components[fid]
                        while self.__component_is_key_addition(fc):
                            if fc.input_connections[0][0] == KEY_ID:
                                fid, fw = fc.input_connections[1]
                            else:
                                fid, fw = fc.input_connections[0]
                            fc = self.components[fid]
                        f.connect_components(o2n[fid], fw, o2n[i], j)
                        if fid != component.input_connections[j][0]:
                            f.components[o2n[i]].input_key_mask |= 1<<j
                    else:
                        f.connect_to_key(o2n[i], j)
            else:
                if component.input_connections[0][0] == KEY_ID:
                    f.components[o2n[component.input_connections[1][0]]].output_key_mask |= 1 << component.input_connections[1][1]
                else:
                    fc = self.components[i]
                    while self.__component_is_key_addition(fc):
                        if fc.input_connections[0][0] == KEY_ID:
                            fid, fw = fc.input_connections[1]
                        else:
                            fid, fw = fc.input_connections[0]
                        fc = self.components[fid]
                    f.components[o2n[fid]].output_key_mask |= 1<<fw
        for i in range(self.output_size):
            fid, fw = self.components[OUTPUT_ID].input_connections[i]
            fc = self.components[fid]
            while self.__component_is_key_addition(fc):
                if fc.input_connections[0][0] == KEY_ID:
                    fid, fw = fc.input_connections[1]
                else:
                    fid, fw = fc.input_connections[0]
                fc = self.components[fid]
            f.connect_components(o2n[fid], fw, OUTPUT_ID, i)
        # update key masks based on global masks:
        for i in range(f.input_size):
            fid, fw = f.components[INPUT_ID].output_connections[i]
            f.components[fid].input_key_mask |= ((global_input_key_mask >> i) & 1) << fw
        for i in range(f.output_size):
            fid, fw = f.components[OUTPUT_ID].input_connections[i]
            f.components[fid].output_key_mask |= ((global_output_key_mask >> i) & 1) << fw
        # run through components of f and optimize the compound functions
        for c in f.components.values():
            if isinstance(c.f, CompoundFunction):
                c.f = c.f.optimized_for_nonzero_trail_detection(c.input_key_mask, c.output_key_mask)
        return f