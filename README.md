# Logical resource estimates for Shor’s elliptic-curve algorithm

This repository implements the finite-size cost model used in the [accompanying manuscript](paper/Resource%20analysis%20of%20Shor%E2%80%99s%20elliptic%20curve%20algorithm%20with%20an%20improved%20quantum%20adder.pdf). It estimates Toffoli and CNOT depth and count, along with peak logical qubit allocation.

## Quick start

Run these commands from the repository root. The estimator requires Python and NumPy. The plotting script also requires Matplotlib.

```sh
python -m pip install numpy
python full.py
python div.py
```

`full.py` prints the 256-bit estimates, selected windows, and qubit allocation by stage. `div.py` prints the 256-bit costs of one inversion round, inversion, and division. Both commands print to the terminal.

For the full 256-bit circuit in depth mode, the current model selects multiplication window `k = 4` and point-addition window `l = 14`.

| Toffoli depth | Toffoli count | CNOT depth | CNOT count | Peak logical qubits |
| ---: | ---: | ---: | ---: | ---: |
| 9,142,554.5 | 1,015,810,991.5 | 20,766,703.5 | 3,111,621,899.5 | 3,870 |

## Export tables and regenerate the plot

| Command from the repository root | Output |
| --- | --- |
| `python write_iterate.py` | `data/data.csv` and `data/data.json` for 25, 32, 192, 256, 384, and 521 bits in depth mode. |
| `python write_iterate.py --figure` | `data/figure1.csv` and `data/figure1.json` for every integer bit width from 2 through 2048. |
| `python output.py --mode depth` | `data/output_256_depth.csv` with gate metrics and category shares for 256-bit arithmetic subroutines and the full circuit. |
| `python output.py --mode count` | `data/output_256_count.csv` using the gate-count-oriented adder and window search. |
| `python plot.py` | `data/figure1.png`, rendered from `data/figure1.json` using Matplotlib. |

The export commands overwrite files at their output paths. The full figure sweep can take several minutes.

After generating `data/figure1.json`, install Matplotlib and render Fig. 1.

```sh
python -m pip install matplotlib
python plot.py
```

Use `python plot.py --help` for custom input and output paths.

[`write.py`](write.py) also exports formatted CSV, raw CSV, and LaTeX-formatted TSV files for the six table sizes to the current working directory.

## Counting conventions

- [`ResourceMetric`](resource_metric.py) stores metrics in the order **Toffoli depth, Toffoli count, CNOT depth, CNOT count**. The CNOT metrics cover algorithm-level CNOTs. Depth is tracked separately for each gate family, and measurement-conditioned uncomputation contributes its expected cost, which can produce fractional values.
- [`add.py`](add.py) provides the adder and comparator cost models. [`mul.py`](mul.py) and [`full.py`](full.py) search finite sets of integer window sizes against the selected Toffoli depth or count objective. The full-circuit estimator groups the `2n` exponent bits into outer windows and budgets any shorter final window separately.
- [`div.py`](div.py) models the inversion pipeline. The two lanes advance independently within each round, and flag cleanup follows both lanes. Their gate counts add, while their depths combine by taking the maximum for each gate family.
- [`qubits.py`](qubits.py) takes the maximum logical allocation across stages. It accounts for active window bits, widened multiplication registers, inversion flags and carry workspace, and reuse of lookup workspace. The resulting width is a logical column budget for the assumed schedule.

The full-circuit CSV includes the selected windows `k` and `l`. Its JSON counterpart stores arrays named `n`, `depth`, `number`, `c_depth`, `c_number`, `qubits`, and `ratios`.
