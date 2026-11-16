# output format: Toffoli depth, Toffoli count, CNOT depth, CNOT count

from math import ceil, log2

from add import AddClass
from resource_metric import delta
from utils import shift, table_unit, untable1


def _get_add_class(add_class):
    return add_class if add_class is not None else AddClass()


def big_mod(n, add_class=None):
    '''
    for modular reduction of a 2n-bit number. n: qubits number.
    '''
    add_class = _get_add_class(add_class)
    q_bit_number = 33
    q_one_number = 7
    res = (q_one_number+1) * add_class.add(n) + (q_one_number-1) * add_class.add(q_bit_number) + add_class.mod(n)
    return res


def mul_no_win(n, add_class=None):
    """Estimate n-bit modular multiplication without windowing."""
    add_class = _get_add_class(add_class)
    if n == 1:
        # Keep the single-bit estimate outside the historical window model.
        return 3 * add_class.add(1) + add_class.add(3) + big_mod(1, add_class)
    # Copy the parity bit before the controlled addition of p in each round.
    step = add_class.ctrl_add(n) + delta([0, 0, 1, 1]) + add_class.ctrl_const_add(n) + shift(n)
    return n * step + add_class.mod(n)


def mul_win(n, k, add_class=None):
    """Estimate n-bit modular multiplication using k-bit windows."""
    add_class = _get_add_class(add_class)
    num_full, num_remain = divmod(n, k)
    ctrl_add = add_class.ctrl_add(n)

    def window_cost(size):
        return (
            size * ctrl_add
            + delta([0, 0, 1, size])
            + add_class.add(n + size)
            + shift(n, size)
            + table_unit(size, n + size)
            + untable1(size)
        )

    result = num_full * window_cost(k) + add_class.mod(n)
    if num_remain:
        return result + window_cost(num_remain)

    return result


def mul_full(n, k, add_class=None):
    if k == 1:
        return mul_no_win(n, add_class)
    return mul_win(n, k, add_class)


def opt_mul_win_size(n, add_class=None):
    """Choose a window size by Toffoli depth or count, according to mode."""
    add_class = _get_add_class(add_class)
    axis = {"depth": 0, "count": 1}.get(add_class.mode)
    if axis is None:
        raise ValueError("Invalid mode. Mode should be 'depth' or 'count'.")
    if n == 1:
        return 1

    mid = max(ceil(log2(log2(n))), 1)
    candidates = range(max(mid - 5, 1), min(mid + 5, n))
    return min(candidates, key=lambda k: mul_full(n, k, add_class).only_value(axis))


def mul(n, add_class=None):
    """Estimate n-bit modular multiplication at the selected window size."""
    add_class = _get_add_class(add_class)
    return mul_full(n, opt_mul_win_size(n, add_class), add_class)


if __name__ == "__main__":
    n = 256
    add_class = AddClass(mode='depth')
    print(f"mul({n}):\n  {mul(n, add_class).format_summary()}")
