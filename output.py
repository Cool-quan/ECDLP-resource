import csv
from argparse import ArgumentParser
from pathlib import Path

from add import AddClass
from div import div, inversion, mod_round
from full import full, point_add
from mul import mul
from resource_metric import METRIC_NAMES, ResourceMetric


def main(mode='depth'):
    n = 256
    addClass = AddClass(mode=mode)
    funcs_add = [
        addClass.add, addClass.const_add, addClass.ctrl_add, addClass.ctrl_const_add,
        addClass.neg, addClass.ctrl_neg, addClass.com, addClass.ctrl_com,
        addClass.mod_add, addClass.ctrl_mod_add, addClass.mod, addClass.mod_doub,
    ]
    funcs_others = [
        mul, mod_round, inversion, div,
        point_add, full
    ]
    data_dir = Path(__file__).resolve().parent / "data"
    data_dir.mkdir(exist_ok=True)
    with open(data_dir / f"output_{n}_{mode}.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['function', *METRIC_NAMES, *ResourceMetric.ratio_headers()])
        for func in funcs_add:
            resource = func(n)
            writer.writerow([func.__name__, *resource.csv_row()])
            print(f"{func.__name__}({n}):\n  {resource.format_summary()}")
        for func in funcs_others:
            resource = func(n, addClass)
            writer.writerow([func.__name__, *resource.csv_row()])
            print(f"{func.__name__}({n}):\n  {resource.format_summary()}")


if __name__ == "__main__":
    parser = ArgumentParser(description="Write the 256-bit ECDLP resource estimates.")
    parser.add_argument("--mode", choices=("count", "depth"), default="depth")
    main(parser.parse_args().mode)
