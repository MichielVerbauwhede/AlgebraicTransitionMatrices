from Modelling.Trails import enum_projected_models, is_key_dependent_limited, get_key_independent_sum, is_key_dependent
from math import inf

def Avec_unified_model_with_partial_trail_counting_constant(cachek1, cachek2, cachek3, cachec1, cachec2, cachec3, oracle_call_counts, model, input_vars, intermediate_vars1, intermediate_vars2, output_vars, f1_model, f1_input_vars, f1_output_vars, f1_key_vars, f2_model, f2_input_vars, f2_output_vars, f2_key_vars, f3_model, f3_input_vars, f3_output_vars, f3_key_vars, uv, limit=inf):
    # This implementation of Avec uses a complete model of F to decide which pairs (a->b) over F_2 should be partially trail counted.
    # This is the most efficient method available in this test suit, given that the complete model of F is somewhat easy to evaluate.
    # If this is not the case, use Avec_disjunct_model_with_partial_trail_counting
    local_cachek1 = {}
    local_cachek2 = {}
    local_cachek3 = {}
    local_cachec1 = {}
    local_cachec2 = {}
    local_cachec3 = {}
    local_count = 0

    invmask = 2**len(input_vars) - 1
    W = set()
    assumptions = []
    u, v = uv
    for i in range(len(input_vars)):
        assumptions.append((-1+2*((u >> i) & 1))*input_vars[i])
    for i in range(len(output_vars)):
        assumptions.append((-1+2*((v >> i) & 1))*output_vars[i])
    for wm in enum_projected_models(model, intermediate_vars1+intermediate_vars2, assumptions):
        a, b = 0, 0
        for x in wm[:len(intermediate_vars1)][::-1]:
            a <<= 1
            if x > 0:
                a += 1
        for x in wm[len(intermediate_vars1):][::-1]:
            b <<= 1
            if x > 0:
                b += 1
        if (u, a) in local_cachek1:
            k1 = local_cachek1[(u, a)]
        elif (u, a) in cachek1:
            k1 = cachek1[(u, a)]
            local_cachek1[(u, a)] = k1
        else:
            k1 = is_key_dependent(f1_model, f1_input_vars, f1_output_vars, f1_key_vars, u, a)
            cachek1[(u, a)] = k1
            local_cachek1[(u, a)] = k1
        if not k1:
            if (u, a) in local_cachec1:
                c1 = local_cachec1[(u, a)]
            elif (u, a) in cachec1:
                c1 = cachec1[(u, a)]
                local_cachec1[(u, a)] = c1
            else:
                c1 = get_key_independent_sum(f1_model, f1_input_vars, f1_output_vars, f1_key_vars, u, a)
                cachec1[(u, a)] = c1
                local_cachec1[(u, a)] = c1
        else:
            c1 = 1
        
        if (b, v) in local_cachek3:
            k3 = local_cachek3[(b, v)]
        elif (b, v) in cachek3:
            k3 = cachek3[(b, v)]
            local_cachek3[(b, v)] = k3
        else:
            k3 = is_key_dependent(f3_model, f3_input_vars, f3_output_vars, f3_key_vars, b, v)
            cachek3[(b, v)] = k3
            local_cachek3[(b, v)] = k3
        if not k3:
            if (b, v) in local_cachec3:
                c3 = local_cachec3[(b, v)]
            elif (b, v) in cachec3:
                c3 = cachec3[(b, v)]
                local_cachec3[(b, v)] = c3
            else:
                c3 = get_key_independent_sum(f3_model, f3_input_vars, f3_output_vars, f3_key_vars, b, v)
                cachec3[(b, v)] = c3
                local_cachec3[(b, v)] = c3
        else:
            c3 = 1

        if c1 and c3:
            local_count += 1
            if (a, b) in local_cachek2:
                k2 = local_cachek2[(a, b)]
            elif (a, b) in cachek2:
                k2 = cachek2[(a, b)]
                local_cachek2[(a, b)] = k2
            else:
                k2 = is_key_dependent_limited(f2_model, f2_input_vars, f2_output_vars, f2_key_vars, a, b, limit)
                cachek2[(a, b)] = k2
                local_cachek2[(a, b)] = k2
            if not k2:
                if (a, b) in local_cachec2:
                    c2 = local_cachec2[(a, b)]
                elif (a, b) in cachec2:
                    c2 = cachec2[(a, b)]
                    local_cachec2[(a, b)] = c2
                else:
                    c2 = get_key_independent_sum(f2_model, f2_input_vars, f2_output_vars, f2_key_vars, a, b)
                    cachec2[(a, b)] = c2
                    local_cachec2[(a, b)] = c2
            else:
                c2 = 1

            if c2:
                if k1 or k3:
                    return True, None
                elif k2:
                    W.add((a, b))
        oracle_call_counts.append(local_count)
    return False, W

def Avec_no_trail_counting(cachek1, cachek3, cachec1, cachec3, oracle_call_counts, model, input_vars, intermediate_vars1, intermediate_vars2, output_vars, f1_model, f1_input_vars, f1_output_vars, f1_key_vars, f3_model, f3_input_vars, f3_output_vars, f3_key_vars, uv):
    # This is a simple implementation of Avec without trail enumeration in F2
    local_cachek1 = {}
    local_cachek3 = {}
    local_cachec1 = {}
    local_cachec3 = {}

    invmask = 2**len(input_vars) - 1
    W = set()
    assumptions = []
    u, v = uv
    for i in range(len(input_vars)):
        assumptions.append((-1+2*((u >> i) & 1))*input_vars[i])
    for i in range(len(output_vars)):
        assumptions.append((-1+2*((v >> i) & 1))*output_vars[i])
    for wm in enum_projected_models(model, intermediate_vars1+intermediate_vars2, assumptions):
        a, b = 0, 0
        for x in wm[:len(intermediate_vars1)][::-1]:
            a <<= 1
            if x > 0:
                a += 1
        for x in wm[len(intermediate_vars1):][::-1]:
            b <<= 1
            if x > 0:
                b += 1
        if (u, a) in local_cachek1:
            k1 = local_cachek1[(u, a)]
        elif (u, a) in cachek1:
            k1 = cachek1[(u, a)]
            local_cachek1[(u, a)] = k1
        else:
            k1 = is_key_dependent(f1_model, f1_input_vars, f1_output_vars, f1_key_vars, u, a)
            cachek1[(u, a)] = k1
            local_cachek1[(u, a)] = k1
        if not k1:
            if (u, a) in local_cachec1:
                c1 = local_cachec1[(u, a)]
            elif (u, a) in cachec1:
                c1 = cachec1[(u, a)]
                local_cachec1[(u, a)] = c1
            else:
                c1 = get_key_independent_sum(f1_model, f1_input_vars, f1_output_vars, f1_key_vars, u, a)
                cachec1[(u, a)] = c1
                local_cachec1[(u, a)] = c1
        else:
            c1 = 1
        
        if (b, v) in local_cachek3:
            k3 = local_cachek3[(b, v)]
        elif (b, v) in cachek3:
            k3 = cachek3[(b, v)]
            local_cachek3[(b, v)] = k3
        else:
            k3 = is_key_dependent(f3_model, f3_input_vars, f3_output_vars, f3_key_vars, b, v)
            cachek3[(b, v)] = k3
            local_cachek3[(b, v)] = k3
        if not k3:
            if (b, v) in local_cachec3:
                c3 = local_cachec3[(b, v)]
            elif (b, v) in cachec3:
                c3 = cachec3[(b, v)]
                local_cachec3[(b, v)] = c3
            else:
                c3 = get_key_independent_sum(f3_model, f3_input_vars, f3_output_vars, f3_key_vars, b, v)
                cachec3[(b, v)] = c3
                local_cachec3[(b, v)] = c3
        else:
            c3 = 1

        if c1 and c3:
            W.add((a, b))
        oracle_call_counts.append(len(W))

    return False, W