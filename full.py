# output format: Toffoli depth, Toffoli count, CNOT depth, CNOT count

from utils import *
from add import AddClass
from mul import *
from div import *
from qubits import full_qubit_stages
from resource_metric import ResourceMetric, delta, zero
from math import ceil
from numpy import log2 as log
import numpy as np

def _get_add_class(add_class):
    return add_class if add_class is not None else AddClass()


def _point_products(n, add_class):
    # Two raw squarings add 2n CNOTs each and at most 2n CNOT layers each.
    # The full multiplication also makes one bitwise copy.
    squaring_cnot = delta([0, 0, 4*n, 4*n])
    copy_cnot = delta([0, 0, 1, n])
    return 2*div(n, add_class) + 6*mul(n, add_class) + squaring_cnot + copy_cnot


def point_add(n, add_class=None):
    """
    for n-bit point addition. n: qubits number.
    """
    add_class = _get_add_class(add_class)
    addition = add_class.mod_add(n)
    raw = _point_products(n, add_class) + 3*add_class.ctrl_mod_add(n) + 5*addition + add_class.ctrl_neg(n)
    depth = raw - 2*addition
    return ResourceMetric.from_axes([depth, raw, depth, raw])

def full_no_win(n, add_class=None):
    """
    for n-bit full addition without window. n: qubits number.
    """
    add_class = _get_add_class(add_class)
    return 2*n * point_add(n, add_class)

def full_win(n, l, add_class=None):
    """
    for n-bit full addition with window. n: qubits number. l: window size.
    """
    add_class = _get_add_class(add_class)
    return _full_win(n, l, add_class, _window_arithmetic(n, add_class))


def _window_arithmetic(n, add_class):
    addition = add_class.mod_add(n)
    ctrl_negation = add_class.ctrl_neg(n)
    # Fig. 10 contains four controlled negations and one ordinary negation.
    return (
        _point_products(n, add_class) + 7*addition
        + 4*ctrl_negation + add_class.neg(n)
    )


def _full_win(n, l, add_class, arithmetic):
    num_full, num_remain = divmod(2*n, l)
    addition = add_class.mod_add(n)

    def window(width):
        address_bits = width - 1
        table_metric = (
            2*table_unit(address_bits, 2*n) + table_unit(address_bits, n)
            + untable2(address_bits)
        )
        # Complement and restore the address, with serial CNOTs on the sign bit.
        address_cnot = delta([0, 0, 2*address_bits, 2*address_bits])
        raw = arithmetic + table_metric + address_cnot
        depth = raw - 2*addition
        return ResourceMetric.from_axes([depth, raw, depth, raw])

    if num_remain > 1:
        remain = window(num_remain)
    elif num_remain == 1:
        remain = point_add(n, add_class)
    else:
        remain = zero()

    return num_full * window(l) + remain

def full_full(n, l, add_class=None):
    '''
    for n-bit full addition. n: qubits number. l: window size.
    '''
    add_class = _get_add_class(add_class)
    if l == 1: 
        return full_no_win(n, add_class)
    else:
        return full_win(n, l, add_class)
    
def opt_full_win_size(n, add_class=None):
    """Choose the outer window by Toffoli depth or gate count, according to mode."""
    add_class = _get_add_class(add_class)
    axis = {'depth': 0, 'count': 1}.get(add_class.mode)
    if axis is None:
        raise ValueError("Invalid mode. Mode should be 'depth' or 'count'.")
    if n == 1:
        return 1

    arithmetic = _window_arithmetic(n, add_class)
    mid = ceil(log(n))+6
    candidates = range(max(mid - 6, 1), min(mid + 6, n))
    return min(
        candidates,
        key=lambda l: (full_no_win(n, add_class) if l == 1 else _full_win(n, l, add_class, arithmetic))[axis],
    )

def full(n, add_class=None):
    """
    for n-bit full addition. n: qubits number.
    """
    add_class = _get_add_class(add_class)
    l = opt_full_win_size(n, add_class)
    return full_full(n, l, add_class)

def full_win_show(n, add_class=None):
    '''
    show the Toffoli depth for different window size.
    '''
    add_class = _get_add_class(add_class)
    mid = ceil(log(n))+6
    for l in range(max(mid-6,1), min(mid+10,n)):
        temp = full_full(n, l, add_class)
        print(f"window size: {l}: {temp}")

if __name__ == "__main__":
    add_class = AddClass(mode='depth')
    # for n in [25, 32, 64, 128, 192, 256, 521, 1024, 2048, 4096, 8192, 16384]:
    # for n in range(1258,1266):
    # for n in [2, 4, 8, 16, 32, 64, 128, 256]:
    for n in [256]:
        print('-'*20)
        # print(f"optimal multiplication size for {n}-bit full addition: {opt_mul_win_size(n, add_class)}")
        print(f'add({n}):\n  {add_class.add(n).format_summary()}')
        print(f'mul({n}):\n  {mul(n, add_class).format_summary()}')
        k = opt_mul_win_size(n, add_class)
        l = opt_full_win_size(n, add_class)
        print(f"Optimal windows for n={n}: k={k}, l={l}")
        print(f'point_add({n}):\n  {point_add(n, add_class).format_summary()}')
        print(f"full({n}):\n  {full_full(n, l, add_class).format_summary()}")
        widths = full_qubit_stages(n, k, l)
        print(f"Peak logical qubits: {max(widths.values())}")
        for stage, width in widths.items():
            print(f"  {stage}: {width}")
        # full_win_show(n, add_class)
        # print('-'*20)
        # mul_win_show(n)
    # full_win_show(256)
    # print(opt_full_win_size(256))
    # print(full_win(25,10))
