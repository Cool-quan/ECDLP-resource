import csv
import numpy as np
from math import log10, floor
from add import *
from mul import *
from div import *
from full import *
from qubits import add_qubits, mod_add_qubits, mul_qubits, division_qubits, point_add_qubits, full_qubits
from resource_metric import ResourceMetric

def scientific_format(n):
    if n < 10000:
        return n, 0
    else:
        k = floor(log10(n))           # 计算指数部分
        m = n / (10 ** k)             # 计算尾数部分
        m = round(m, 2)               # 保留两位小数
        return m, k
    
def format_number_str(n):
    if n < 10000:
        return str(n)
    else:
        m, k = scientific_format(n)
        return f"{m:.2f}e{k}"                # 使用科学计数法表示
    
def format_number_latex(n):
    if n < 10000:
        return f"${n}$"
    else:
        m, k = scientific_format(n)
        return f"${m:.2f} \\times 10^{{{k}}}$"   # 使用科学计数法表示

def format_ratio(n):
    return f"{100*n:.2f}%"


def write(n):
    with open(str(n)+'.csv', 'w', newline='') as csv1, open(str(n)+'_raw.csv', 'w', newline='') as csv2, open(str(n)+'.tsv', 'w', newline='') as tsv:
        writer1 = csv.writer(csv1)
        writer2 = csv.writer(csv2)
        writer3 = csv.writer(tsv, delimiter='\t')
        ratio_headers = ResourceMetric.ratio_headers()
        writer1.writerow(['function', 'qubit number', 'CNOT depth', 'CNOT count', 'Toffoli depth', 'Toffoli count', 'Toffoli volume', *ratio_headers])
        writer2.writerow(['function', 'qubit number', 'CNOT depth', 'CNOT count', 'Toffoli depth', 'Toffoli count', 'Toffoli volume', *ratio_headers])
        writer3.writerow(['function', 'qubit number', 'CNOT depth', 'CNOT count', 'Toffoli depth', 'Toffoli count', 'Toffoli volume', *ratio_headers])
        
        # Write the data
        def write_row(func, qubit_num, resource):
            volume = resource[0] * qubit_num
            ratios = resource.ratio_row()
            writer1.writerow([func, format_number_str(qubit_num), format_number_str(resource[2]), format_number_str(resource[3]), format_number_str(resource[0]), format_number_str(resource[1]), format_number_latex(volume), *[format_ratio(ratio) for ratio in ratios]])
            writer2.writerow([func, qubit_num, resource[2], resource[3], resource[0], resource[1], volume, *ratios])
            writer3.writerow([func, format_number_latex(qubit_num), format_number_latex(resource[2]), format_number_latex(resource[3]), format_number_latex(resource[0]), format_number_latex(resource[1]), format_number_latex(volume), *[format_ratio(ratio) for ratio in ratios]])

        k = opt_mul_win_size(n)
        l = opt_full_win_size(n)
        write_row('integer addition', add_qubits(n), add(n))
        write_row('modular addition', mod_add_qubits(n), mod_add(n))
        write_row('modular multiplication', mul_qubits(n, k), mul(n))
        write_row('modular division', division_qubits(n, k), div(n))
        write_row('point addition (without window)', point_add_qubits(n, k), point_add(n))
        write_row('full circuit', full_qubits(n, k, l), full_full(n, l))
            

if __name__ == "__main__":
    write(25)
    write(32)
    write(192)
    write(256)
    write(384)
    write(521)
    # print(mod_round(1.222))
