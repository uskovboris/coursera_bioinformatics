from fs_helpers import read_lines
from dna_lib import *

if __name__ == '__main__':
    dataset = read_lines(r"C:\src\coursera_bioinformatics\approximate_patterns_matchings_calculator_gbg.txt")
    pattern = dataset[0]
    text = dataset[1]
    d = int(dataset[2])

    result = approximate_pattern_positions(pattern, text, d)
    print(format_list_result(result))
