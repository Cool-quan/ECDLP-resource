# output format: Toffoli depth, Toffoli count, CNOT depth, CNOT count

from math import ceil, floor
from numpy import log2 as log
import numpy as np
from resource_metric import DELTA, MISC, TABLE, delta, metric, zero

def ghz(k):
    """Integer estimate for the paper's 3k/2 GHZ preparation cost."""
    return k + k//2


def ghz_metric(k, axis=3):
    """
    A GHZ contribution placed on one metric axis.
    """
    values = np.array([0, 0, 0, 0], dtype=object)
    values[axis] = ghz(k)
    return metric(values, MISC)

def swap(n):
    """
    for swap. n: qubits number.
    """
    return zero()

def ctrl_swap(n):
    """
    for ctrl_swap. n: qubits number.
    """
    return metric([1, n, 5, 2*n + ghz(n)], MISC)

def shift(n, k=1):
    """
    for k-bit shift. n: qubits number.
    """
    return zero()

def ctrl_shift(n, k=1):
    """
    for k-bit ctrl_shift. n: qubits number.
    """
    # Only the single-bit shift is specified in the paper.
    toffoli_count = 2*n if k == 1 else n*k + n//2*2
    direct_cnot_count = 2*n if k == 1 else toffoli_count
    return metric([k+2, toffoli_count, k+5, direct_cnot_count + ghz(n+k)], MISC)

def and_gate(n):
    """
    for n-bit and_gate.
    """
    return delta([ceil(log(n)), n, 0, 0])

def table_unit(k, m=1):
    """
    for quantum table unit. k: input qubits number. m: output qubits number.
    """
    if k == 0:
        return zero()
    exp = 2 ** k
    # The paper's average includes data CNOTs, GHZ preparation and traversal.
    cnot_count = exp*(2*m + 1)
    return metric([exp-2, exp-2, exp*5, cnot_count], TABLE)

def unary(n):
    """
    for unary encoding and its inverse. n: qubits number.
    """
    if n == 0:
        return zero()
    exp = 2 ** n
    cnot_count = 5*(exp//2) - 1
    return metric([n, exp-1, n+6, cnot_count], TABLE)


def _untable_split(k, h):
    encoding = unary(h)
    resource = encoding + table_unit(k-h, 2**h)
    # The inverse one-hot conversion still uses direct CNOTs.
    return resource.with_cnot(resource + encoding)

def untable1(k):
    """
    for first method to uncompute a table. k: input qubits number. m: output qubits number.
    """
    if k == 0:
        return zero()
    return _untable_split(k, k-1)

def untable2(k):
    """
    for second method to uncompute a table. k: input qubits number. m: output qubits number.
    """
    # Both balanced splits tie on Toffoli count. Round the unary side up
    # to preserve the existing split and its lower lookup depth.
    return _untable_split(k, k - k//2)

def table(k, m=1, anc_c=0, anc_d=0, mode='count'):
    """
    for improved quantum table. k: input qubits number. m: output qubits number. anc_c: number of clean ancillas. anc_d: number of dirty ancillas. mode: 'depth' or 'count'
    """
    if k == 0:
        return zero()
    # the index of the metric axis to optimize
    if mode == 'depth':
        com_idx = 0
    elif mode == 'count':
        com_idx = 1
    else:
        raise ValueError(f"Invalid mode: {mode}")

    # dirty ancilla variant
    N = 2 ** k
    anc = anc_c + anc_d
    dirty_res = table_unit(k,m)
    for log_lamb in range(1, k):
        lamb = 2 ** log_lamb
        mu = anc // (lamb - 1)
        if mu == 0:
            break
        temp_res = (ceil(m/mu)+1) * (table_unit(k-log_lamb, lamb*mu) + unary(log_lamb))\
            + (lamb-1) * metric([1, m+mu, 0, 0], TABLE)
        if temp_res.only_value(com_idx) < dirty_res.only_value(com_idx):
            dirty_res = temp_res
    # clean ancilla variant
    lamb_temp = anc_c // m + 1
    log_lamb = floor(log(lamb_temp))
    lamb = 2 ** log_lamb
    clean_res = table_unit(k-log_lamb, lamb * m) + untable2(k-log_lamb) + metric([1, lamb * m, 0, 0], TABLE)
    res = dirty_res if dirty_res.only_value(com_idx) < clean_res.only_value(com_idx) else clean_res
    return res

if __name__ == "__main__":
    # print("ghz(3):", ghz(3))
    # print("table(3, 2):", table(3, 2))
    # print("swap(3):", swap(3))
    # print("ctrl_swap(3):", ctrl_swap(3))
    # print("shift(3, 2):", shift(3, 2))
    # print("ctrl_shift(3, 2):", ctrl_shift(3, 2))
    # print("and_gate(3):", and_gate(3))
    # print(swap(3) + ctrl_swap(3))
    # funcs = [
    #     ghz, table, untable1, untable2,
    #     swap, ctrl_swap, shift, ctrl_shift,
    #     and_gate
    # ]
    # n = 4
    # for func in funcs:
    #     resource = func(n)
    #     if hasattr(resource, "format_summary"):
    #         print(f"{func.__name__}({n}):\n  {resource.format_summary()}")
    #     else:
    #         print(f"{func.__name__}({n}):", resource)
    n = 256
    print(table_unit(16, n))
    print(table(16, n, anc_c=0, anc_d=3*n))
