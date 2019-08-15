import sys


def main():
    if len(sys.argv) != 2:
        show_usage()

    dataset_file = sys.argv[1]

    with open(dataset_file) as f:
        dataset = f.readlines()

    if len(dataset) != 1:
        print('Dataset should contains exactly 1 line with template strand')
        sys.exit(1)

    template_strand = dataset[0]
    complement_strand = ""

    for nucleotide in template_strand:
        if nucleotide == 'A':
            complement_strand = 'T' + complement_strand
        elif nucleotide == 'T':
            complement_strand = 'A' + complement_strand
        elif nucleotide == 'G':
            complement_strand = 'C' + complement_strand
        elif nucleotide == 'C':
            complement_strand = 'G' + complement_strand

    print(complement_strand)


if __name__ == '__main__':
    main()
