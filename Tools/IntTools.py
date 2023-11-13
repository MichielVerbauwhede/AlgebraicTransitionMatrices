def pext(x, m):
    j = 0
    r = 0
    for i in range(m.bit_length()):
        if m & (1<<i):
            r |= ((x>>i) & 1)<<j
            j += 1
    return r