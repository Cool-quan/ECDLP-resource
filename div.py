# output format: Toffoli depth, Toffoli count, CNOT depth, CNOT count

from add import AddClass
from mul import mul
from resource_metric import ResourceMetric, axis_max, delta
from utils import and_gate, ctrl_shift, ctrl_swap


def _parallel(left, right):
    counts = left + right
    return ResourceMetric.from_axes([
        axis_max(left, right, 0), counts,
        axis_max(left, right, 2), counts,
    ])


def _uv_stages(n, add_class, *, first=False):
    swap = ctrl_swap(n)

    # In the first round the incoming enable is the known constant 1.
    flag_inputs = n if first else n + 1
    flag_cnots = 1 if first else 2
    flags = and_gate(flag_inputs) + delta([2, 2, flag_cnots, flag_cnots])
    comparison = _parallel(add_class.com(n), delta([0, 0, 2, 2]))
    # Compute the three-control predicate and XOR it into both a and l.
    branch_update = delta([2, 2, 2, 2])
    select = flags + comparison + branch_update + swap

    # A temporary logical AND combines f and !b for the controlled subtract.
    subtract = delta([1, 1, 0, 0]) + add_class.ctrl_add(n)
    clear_b_and_shift = _parallel(delta([0, 0, 2, 2]), ctrl_shift(n))
    return select, subtract + clear_b_and_shift + swap


def _rs_stages(n, add_class):
    swap = ctrl_swap(n)
    # XOR l' into a', condition the addition on f' and !a', then restore a'.
    add = swap + delta([1, 1, 2, 2]) + add_class.ctrl_add(n)
    double = add_class.mod_doub(n, controlled=False) + swap + delta([0, 0, 1, 1])
    return add, double


def mod_round(n, add_class=None):
    """Estimate a steady-state v32 round on current (u,v) and previous (r,s)."""
    add_class = add_class if add_class is not None else AddClass()
    uv_select, uv_subtract = _uv_stages(n, add_class)
    rs_add, rs_double = _rs_stages(n, add_class)
    # Copy and uncopy the incoming flag f' for the two lanes, then clear it
    # with the controlled XOR and CNOT shown at the end of the round.
    flag_overhead = delta([1, 1, 3, 3])
    # Each lane owns its arithmetic workspace and can advance independently.
    return _parallel(uv_select + uv_subtract, rs_add + rs_double) + flag_overhead


def inversion(n, add_class=None):
    """Estimate the fixed 2n iterations, pipelined over 2n+1 rounds."""
    add_class = add_class if add_class is not None else AddClass()
    fill = sum(_uv_stages(n, add_class, first=True))
    # Budget the r,s stages in the drain using the same primitives as an inner round.
    drain = sum(_rs_stages(n, add_class))
    return fill + (2 * n - 1) * mod_round(n, add_class) + drain + add_class.neg(n)


# The two directions now use the same circuit and resource model.
inv1 = inv2 = inversion
mod_round_doub = mod_round


def div(n, add_class=None):
    """Estimate compute, controlled copy, and uncompute for modular division."""
    add_class = add_class if add_class is not None else AddClass()
    # Retain the existing division copy convention for the before/after comparison.
    copy = delta([1, n, 0, 0])
    return 2 * (inversion(n, add_class) + mul(n, add_class)) + copy


if __name__ == "__main__":
    n = 256
    add_class = AddClass(mode='depth')
    for func in (mod_round, inversion, div):
        print(f"{func.__name__}({n}):\n  {func(n, add_class).format_summary()}")
