"""Logical width budgets using the v32 column layout and workspace reuse."""


def _adder_workspace(width):
    # Algorithm 3 uses Carry_(width+1) and an internal high endpoint w.
    return 3 * (width + 1) + 1


def _mod_add_workspace(n):
    # Constant p, adder workspace, and the high carry/reduction flag.
    return n + _adder_workspace(n) + 1


def add_qubits(n):
    return 2 * n + _adder_workspace(n) + 1


def mod_add_qubits(n):
    return 2 * n + _mod_add_workspace(n)


def _lookup_workspace(address_bits, output_bits):
    # One traversal stack and one GHZ copy per output after the first.
    return address_bits + output_bits - 1 if address_bits else 0


def _unlookup_workspace(address_bits, unary_bits):
    # The original lookup output has been measured and reset. The unary
    # fanout and phase-lookup fanout run at different times.
    unary_width = 1 << unary_bits
    phase_lookup = _lookup_workspace(address_bits - unary_bits, unary_width)
    return unary_width + max(unary_width // 2, phase_lookup)


def mul_qubits(n, k):
    """Budget the raw multiplier, including its result and history columns."""
    column_layout = 9 * n
    if k == 1:
        # The non-windowed path needs the high accumulator bit and a flag.
        return max(column_layout, 4 * n + 1 + _mod_add_workspace(n))

    k = min(k, n)
    # x, y, m, the widened accumulator, and a local copy of the table address.
    live = 3 * n + (n + k + 1) + k
    table_output = n + k
    return max(
        column_layout,
        live + table_output + _adder_workspace(n + k),
        live + table_output + _lookup_workspace(k, table_output),
        live + _unlookup_workspace(k, k - 1),
    )


def _division_stages(n, k):
    # u,v,r,s, l1,l2, y and lambda stay allocated during inversion. The
    # two lanes own separate workspace, with p and the high bit on r,s.
    # The comparator bank has room for its branch predicate. The next
    # predicate can reuse c, and the lower lane borrows idle p during add.
    flags = 6 + 1  # a, a', b, c, f, f', and the independent f' copy
    return {
        'division inversion': 8 * n + _adder_workspace(n) + _mod_add_workspace(n) + flags,
        # x, lambda and the two inversion-history columns survive the product.
        'division multiplication': mul_qubits(n, k) + 4 * n,
    }


def division_qubits(n, k):
    """Budget a controlled division, including the external control bit."""
    return max(_division_stages(n, k).values()) + 1


def point_add_qubits(n, k):
    # The division peak includes the same x, y, lambda and control registers.
    # Outside division, the product has only one extra destination column.
    return division_qubits(n, k)


def full_qubit_stages(n, k, l):
    """Return widths of stages sharing the lattice and its reusable columns."""
    if l == 1:
        return {'controlled point addition': point_add_qubits(n, k)}

    l = min(l, 2 * n)
    address_bits = l - 1
    # The signed window includes the sign bit. Its l qubits persist through
    # both divisions. In qpic/point-add-win, TD1 resets both outputs before
    # the first division and the second T1 runs after the second division.
    stages = _division_stages(n, k)
    stages.update({
        'point multiplication': mul_qubits(n, k) + n,
        'parallel point additions': 5 * n + 2 * _mod_add_workspace(n),
        'point lookup': 5 * n + _lookup_workspace(address_bits, 2 * n),
        'point lookup uncomputation': 3 * n + _unlookup_workspace(
            address_bits, address_bits - address_bits // 2,
        ),
    })
    return {stage: width + l for stage, width in stages.items()}


def full_qubits(n, k, l):
    return max(full_qubit_stages(n, k, l).values())
