import csv
import json
from argparse import ArgumentParser
from pathlib import Path

from add import AddClass
from full import full_full, opt_full_win_size
from mul import opt_mul_win_size
from qubits import full_qubits
from resource_metric import ResourceMetric


TABLE_SIZES = (25, 32, 192, 256, 384, 521)
DATA_DIR = Path(__file__).resolve().parent / 'data'


def write_estimates(sizes, stem):
    add_class = AddClass(mode='depth')
    ratio_headers = ResourceMetric.ratio_headers()
    data = {key: [] for key in ('n', 'depth', 'number', 'c_depth', 'c_number', 'qubits', 'ratios')}
    DATA_DIR.mkdir(exist_ok=True)
    with (DATA_DIR / f'{stem}.csv').open('w', newline='', encoding='utf-8') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['n', 't_depth', 't_number', 'c_depth', 'c_number',
                         'qubit number', 'multiplication window k', 'point-addition window l',
                         *ratio_headers])
        for n in sizes:
            outer_window = opt_full_win_size(n, add_class)
            inner_window = opt_mul_win_size(n, add_class)
            resource = full_full(n, outer_window, add_class)
            qubits = full_qubits(n, inner_window, outer_window)
            ratios = resource.ratio_row()
            for field, value in zip(('n', 'depth', 'number', 'c_depth', 'c_number'), (n, *resource)):
                data[field].append(value)
            data['qubits'].append(qubits)
            data['ratios'].append(dict(zip(ratio_headers, ratios)))
            writer.writerow([n, *resource, qubits, inner_window, outer_window, *ratios])

    with (DATA_DIR / f'{stem}.json').open('w', encoding='utf-8') as json_file:
        json.dump(data, json_file, indent=4)
        json_file.write('\n')


if __name__ == '__main__':
    parser = ArgumentParser(description='Export the finite-size depth-optimized estimates.')
    parser.add_argument('--figure', action='store_true', help='Export every size from 2 to 2048 for Fig. 1.')
    if parser.parse_args().figure:
        write_estimates(range(2, 2049), 'figure1')
    else:
        write_estimates(TABLE_SIZES, 'data')
