import fs_helpers
import neighbors


def calculate_neighbours(path):
    raw_input = fs_helpers.read_lines(path)
    output = ' '.join(neighbors.neighbors(raw_input[0], int(raw_input[1])))
    print(output)


if __name__ == '__main__':
    calculate_neighbours(r"C:\Users\Boris\Downloads\dataset_3014_4 (1).txt")