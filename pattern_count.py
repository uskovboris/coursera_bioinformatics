import sys

def show_usage():
	print("Usage:")
	print("python pattern_count.py <dataset_file>")

def count_pattern(chromosome, pattern):
	pattern_len = len(pattern)
	count = 0
	for i in range(len(chromosome)):
		if (chromosome[i:i+pattern_len] == pattern):
			count += 1
	return count

def main():

	if len(sys.argv) !=2:
		show_usage()

	dataset_file = sys.argv[1]

	with open(dataset_file) as f:
		dataset = f.readlines()

	if len(dataset) != 2:
		print('Dataset should contains 2 lines')
		sys.exit(1)

	chromosome = dataset[0]
	pattern = dataset[1]

	print("chromosome: %s" % chromosome)
	print("pattern: %s" % pattern)

	cnt = count_pattern(chromosome, pattern)
	print("Pattern inclusion count is %d" % cnt)

if __name__ == '__main__':
	main()
