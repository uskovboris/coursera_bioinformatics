from fs_helpers import read_lines
from dna_lib import distance_between_pattern_and_strings

if __name__ == '__main__':
    dataset = read_lines(r"C:\Users\Boris\Downloads\dataset_5164_1 (1).txt")
    pattern = dataset[0]
    text = dataset[1]
    dna = text.split()

    print(distance_between_pattern_and_strings(pattern, dna))