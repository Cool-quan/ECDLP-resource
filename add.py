# output format: Toffoli depth, Toffoli count, CNOT depth, CNOT count

from utils import *
from resource_metric import ADDITION, COMPARISON, DELTA, delta, metric
from math import ceil, floor
from numpy import log2 as log
import numpy as np


def _expected_comparison(comparison):
    # Appendix C.1 charges the comparator with probability one half.
    return comparison * 0.5


def fib(n):
    if n < 0:
        raise ValueError("n must be non-negative")
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a

class AddClass:
    def __init__(self, mode='depth'):
        self.mode = mode

    @staticmethod
    def w(n):
        '''
        return the number of 1s in binary representation of n
        '''
        return bin(n).count('1')

    @staticmethod
    def step1(n):
        '''
        for step 1 of carry computation.
        '''
        k = max(1, floor(log(n)))
        if k == 1:
            return metric([0, 0, 0, 0], ADDITION)
        full_blocks, remainder = divmod(n, k)
        t_depth = floor(log(k)) + floor(log(2*k/3))
        t_count = full_blocks * (4*k - 2*AddClass.w(k) - 2*floor(log(k)) - 2)
        if remainder > 1:
            t_count += 4*remainder - 2*AddClass.w(remainder) - 2*floor(log(remainder)) - 2
        c_depth = t_depth + 2
        c_count = t_count
        return metric([t_depth, t_count, c_depth, c_count], ADDITION)

    @staticmethod
    def step2(n):
        '''
        for step 2 of carry computation.
        '''
        k = max(1, floor(log(n))) # block size
        m = ceil(n / k) # number of blocks
        if m == 1:
            return metric([0, 0, 0, 0], ADDITION)
        t_depth = ceil(log(m))
        t_count = 0
        for i in range(m):
            t_count += 2*AddClass.w(i)
        c_depth = t_depth + 1
        c_count = t_count
        return metric([t_depth, t_count, c_depth, c_count], ADDITION)

    @staticmethod
    def step3(n):
        '''
        for step 3 of carry computation.
        '''
        k = max(1, floor(log(n)))
        if k == 1:
            return metric([0, 0, 0, 0], ADDITION)
        full_blocks, remainder = divmod(n, k)
        m = full_blocks + int(remainder > 0)
        t_depth = 2
        t_count = n-m
        # Prepare the merge GHZ states in three serial CNOT layers.
        c_depth = 3
        c_count = full_blocks * ghz(k - 1) + ghz(max(0, remainder - 1))
        return metric([t_depth, t_count, c_depth, c_count], ADDITION)

    @staticmethod
    def compute_carry(n):
        '''
        for carry computation. n: qubits number.
        '''
        return AddClass.step1(n) + AddClass.step2(n) + AddClass.step3(n)

        # better adder
        l = floor(log(n))
        size_opt = {32:56, 64:125, 128:272, 256:575, 512:1190, 1024:2437}
        if n in size_opt:
            size = size_opt[n]
            # print(size)
        elif AddClass.w(n) == 1:
            size = 3*n - fib(4+l)
        else:
            size = 3*n
        depth = l + 1
        count = 2 * size - n
        return metric([depth, count, depth, count], ADDITION)

    def cnot_ladder(self, n):
        """
        for n-bit CNOT ladder. n: qubits number.
        """
        if self.mode == 'depth':
            c_depth = 2
            c_count = 2*n
        elif self.mode == 'count':
            c_depth = floor(log(n)) + floor(log(2*n/3))
            c_count = 2*n - 2 - floor(log(n)) - floor(log(2*n/3))
        else:
            raise ValueError("Invalid mode. Please choose 'depth' or 'count'.")
        return delta([0, 0, c_depth, c_count])

    def add(self, n):
        """
        for n-bit add. n: qubits number.
        """
        if self.mode == 'depth':
            # The CNOT ladder and Carry include the high endpoint in the paper.
            carry = AddClass.compute_carry(n+1)
            cnot = AddClass.cnot_ladder(self, n+1)
            raw = 2*carry + 2*cnot + delta([0, 0, 4, 3*n+1])
            return raw.reclassify(ADDITION)
        elif self.mode == 'count':
            t_depth = n-1
            t_count = n-1
            c_depth = max(4*n-6, 1)
            c_count = max(6*n-9, 1)
            return metric([t_depth, t_count, c_depth, c_count], ADDITION)
        else:
            raise ValueError("Invalid mode. Please choose 'depth' or 'count'.")
        
    def const_add(self, n):
        """
        for n-bit const_add. n: qubits number.
        """
        if self.mode == 'depth':
            cnot = AddClass.cnot_ladder(self, n+1)
            raw = AddClass.add(self, n) + delta([0, 0, -2, -2*n]) - 2*cnot
            return raw.reclassify(ADDITION)
        elif self.mode == 'count':
            raw = AddClass.add(self, n) + delta([-1, -1, -1, -2])
            return raw.reclassify(ADDITION)
        else:
            raise ValueError("Invalid mode. Please choose 'depth' or 'count'.")

    def ctrl_add(self, n):
        """
        for n-bit ctrl_add. n: qubits number.
        """
        if self.mode == 'depth':
            raw = AddClass.add(self, n) + delta([2, n+1, 6, 4*n+1])
            return raw.reclassify(ADDITION)
        elif self.mode == 'count':
            t_depth = 2*n - 1
            t_count = 2*n - 1
            c_depth = max(5*n-7, 1)
            c_count = max(7*n-11, 1)
            return metric([t_depth, t_count, c_depth, c_count], ADDITION)
        else:
            raise ValueError("Invalid mode. Please choose 'depth' or 'count'.")

    def ctrl_const_add(self, n, *, carry_out=True):
        """Estimate controlled constant addition. The depth mode can omit the carry copy."""
        if self.mode == 'depth':
            cnot = AddClass.cnot_ladder(self, n+1)
            # Algorithm 3 has one controlled data layer and an optional
            # controlled copy of the leading carry (Appendix D.3.b).
            copy = int(carry_out)
            raw = AddClass.add(self, n) + delta([1 + copy, n + copy, 4, 2*n+1]) - 2*cnot
            return raw.reclassify(ADDITION)
        elif self.mode == 'count':
            raw = AddClass.ctrl_add(self, n) + delta([-1, -1, -1, -2])
            return raw.reclassify(ADDITION)
        else:
            raise ValueError("Invalid mode. Please choose 'depth' or 'count'.")

    def neg(self, n):
        """
        for n-bit neation. n: qubits number.
        """
        return AddClass.const_add(self, n)

    def ctrl_neg(self, n):
        """
        for n-bit ctrl_neg. n: qubits number.
        """
        return AddClass.ctrl_const_add(self, n)

    def com(self, n):
        '''
        for n-bit comparator.
        '''
        if self.mode == 'depth':
            levels = ceil(log(n-1))
            t_depth = levels + 1
            t_count = 2*n - levels - 3
            tree_depth = 2*(levels+1)
            tree_count = 2*(n-levels-2)
            cnot = AddClass.cnot_ladder(self, n+1)
            # The two outer XOR layers and the root copy each add one layer.
            raw = metric([t_depth, t_count, tree_depth+3, tree_count+2*n+1], COMPARISON) + 2*cnot
            return raw.reclassify(COMPARISON)
        elif self.mode == 'count':
            return AddClass.add(self, n).reclassify(COMPARISON)
        else:
            raise ValueError("Invalid mode. Please choose 'depth' or 'count'.")

    def ctrl_com(self, n):
        '''
        for n-bit ctrl_comparator.
        '''
        raw = AddClass.com(self, n) + delta([1, 1, -1, -1])
        return raw.reclassify(COMPARISON)

    def const_com(self, n):
        '''
        for n-bit const_comparator.
        '''
        if self.mode == 'depth':
            cnot = AddClass.cnot_ladder(self, n+1)
            raw = AddClass.com(self, n) + delta([0, 0, -2, -2*n]) - 2*cnot
            return raw.reclassify(COMPARISON)
        elif self.mode == 'count':
            raw = AddClass.com(self, n) + delta([-1, -1, -1, -2])
            return raw.reclassify(COMPARISON)
        else:
            raise ValueError("Invalid mode. Please choose 'depth' or 'count'.")

    def ctrl_const_com(self, n):
        '''
        for n-bit ctrl_const_comparator.
        '''
        if self.mode == 'depth':
            raw = AddClass.ctrl_com(self, n)
            cnot = AddClass.const_com(self, n) + delta([0, 0, -1, -1])
            raw = raw.with_cnot(cnot)
            return raw.reclassify(COMPARISON)
        elif self.mode == 'count':
            raw = AddClass.ctrl_com(self, n) + delta([-1, -1, -1, -2])
            return raw.reclassify(COMPARISON)
        else:
            raise ValueError("Invalid mode. Please choose 'depth' or 'count'.")

    def mod_add(self, n):
        """
        for n-bit modular addition.
        """
        # The extra register is the carry output of an n-bit addition (Fig. 13).
        adder = AddClass.add(self, n if self.mode == 'depth' else n+1)
        comparison = _expected_comparison(AddClass.com(self, n))
        return adder + AddClass.ctrl_const_add(self, n) + comparison

    def ctrl_mod_add(self, n):
        """
        for n-bit controlled modular addition.
        """
        adder = AddClass.ctrl_add(self, n if self.mode == 'depth' else n+1)
        comparison = _expected_comparison(AddClass.ctrl_com(self, n))
        return adder + AddClass.ctrl_const_add(self, n) + comparison

    def mod(self, n):
        """
        for n-bit modular reduction.
        """
        return AddClass.ctrl_const_add(self, n)

    def mod_doub(self, n, *, controlled=True):
        """Estimate modular doubling with or without an external control."""
        reduction = self.ctrl_const_add(n)
        if controlled:
            return ctrl_shift(n) + reduction + delta([1, 1, 0, 0])
        # The ordinary shift uses SWAPs. Output parity clears the reduction flag.
        return reduction + delta([0, 0, 1, 1])


def add(n, mode='depth'):
    return AddClass(mode=mode).add(n)


def const_add(n, mode='depth'):
    return AddClass(mode=mode).const_add(n)


def ctrl_add(n, mode='depth'):
    return AddClass(mode=mode).ctrl_add(n)


def ctrl_const_add(n, mode='depth'):
    return AddClass(mode=mode).ctrl_const_add(n)


def neg(n, mode='depth'):
    return AddClass(mode=mode).neg(n)


def ctrl_neg(n, mode='depth'):
    return AddClass(mode=mode).ctrl_neg(n)


def com(n, mode='depth'):
    return AddClass(mode=mode).com(n)


def ctrl_com(n, mode='depth'):
    return AddClass(mode=mode).ctrl_com(n)


def const_com(n, mode='depth'):
    return AddClass(mode=mode).const_com(n)


def ctrl_const_com(n, mode='depth'):
    return AddClass(mode=mode).ctrl_const_com(n)


def mod_add(n, mode='depth'):
    return AddClass(mode=mode).mod_add(n)


def ctrl_mod_add(n, mode='depth'):
    return AddClass(mode=mode).ctrl_mod_add(n)


def mod(n, mode='depth'):
    return AddClass(mode=mode).mod(n)


def mod_doub(n, mode='depth', *, controlled=True):
    return AddClass(mode=mode).mod_doub(n, controlled=controlled)

if __name__ == "__main__":
    # print(w(255))
    # print('mod_add(256):', mod_add(256))
    # print('mod(256):', mod(256))
    # print('mod_doub(256):', mod_doub(256))
    # print(ctrl_add(256))
    addClass = AddClass(mode='count')
    funcs = [
        addClass.step1, addClass.step2, addClass.step3,
        addClass.compute_carry, addClass.add, addClass.const_add, addClass.ctrl_add, addClass.ctrl_const_add,
        addClass.neg, addClass.ctrl_neg, addClass.com, addClass.ctrl_com,
        addClass.mod_add, addClass.ctrl_mod_add, addClass.mod, addClass.mod_doub
    ]
    n = 256
    for func in funcs:
        resource = func(n)
        print(f"{func.__name__}({n}):\n  {resource.format_summary()}")
