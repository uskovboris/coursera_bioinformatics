from fs_helpers import read_lines
from dna_lib import approximate_pattern_positions

if __name__ == '__main__':
    dataset = read_lines(r"C:\src\coursera_bioinformatics\approximate_patterns_count_dbg.txt")
    pattern = dataset[0]
    text = dataset[1]
    d = int(dataset[2])

    print(approximate_pattern_positions(text, pattern, d))