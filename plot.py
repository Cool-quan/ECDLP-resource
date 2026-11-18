import json
from argparse import ArgumentParser
from pathlib import Path

from matplotlib.figure import Figure


DATA_DIR = Path(__file__).resolve().parent / 'data'
HIGHLIGHTED_SIZES = {25, 32, 64, 128, 192, 256, 384, 521}


def main():
    parser = ArgumentParser(description='Render the Toffoli resource estimates for Fig. 1.')
    parser.add_argument('--data', type=Path, default=DATA_DIR / 'figure1.json')
    parser.add_argument('--output', type=Path, default=DATA_DIR / 'figure1.png')
    args = parser.parse_args()

    with args.data.open(encoding='utf-8') as file:
        data = json.load(file)

    sizes = data['n']
    figure = Figure()
    depth_axis = figure.subplots()
    count_axis = depth_axis.twinx()

    depth_line, = depth_axis.plot(
        sizes, data['depth'], color='tab:blue', linewidth=1.5, label='Toffoli depth'
    )
    depth_axis.set_xlabel('input size (bit number)')
    depth_axis.set_ylabel('Toffoli depth', color='tab:blue')
    depth_axis.tick_params(axis='y', labelcolor='tab:blue')
    depth_axis.set_xscale('log')
    depth_axis.set_yscale('log')
    depth_axis.set_ylim(9e2, 2e10)
    depth_axis.grid(True, which='both', linestyle='--', linewidth=0.7, alpha=0.7)

    count_line, = count_axis.plot(
        sizes, data['number'], color='tab:red', linewidth=1.5, label='Toffoli number'
    )
    count_axis.set_ylabel('Toffoli number', color='tab:red')
    count_axis.tick_params(axis='y', labelcolor='tab:red')
    count_axis.set_yscale('log')

    for size, depth, count in zip(sizes, data['depth'], data['number']):
        if size in HIGHLIGHTED_SIZES:
            depth_axis.scatter(size, depth, color='tab:blue', edgecolor='black', zorder=5)
            depth_axis.text(size, depth * 0.8, str(size), fontsize=8,
                            color='tab:blue', ha='left', va='top')
            count_axis.scatter(size, count, color='tab:red', edgecolor='black', zorder=5)
            count_axis.text(size, count * 1.2, str(size), fontsize=8,
                            color='tab:red', ha='right', va='bottom')

    depth_axis.legend(handles=[depth_line, count_line], loc='upper left')
    figure.tight_layout()
    figure.savefig(args.output, dpi=300)


if __name__ == '__main__':
    main()
