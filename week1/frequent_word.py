import sys

def show_usage():
	print("Usage:")
	print("python pattern_count.py <dataset_file>")

def count_pattern(chromosome, pattern, pattern_len):
	count = 0
	for i in range(len(chromosome)):
		if (chromosome[i:i+pattern_len] == pattern):
			count += 1
	return count

def frequent_word(chromosome, pattern_len):
	
	frequent_words = []

	counts_dict = dict()

	i = 0
	max_pattern_count = 0
	for i in range(len(chromosome)-pattern_len):
		current_pattern = chromosome[i:i+pattern_len]
		counts_dict[current_pattern] = count_pattern(chromosome, current_pattern, pattern_len)
		print (current_pattern, counts_dict[current_pattern])
		if (max_pattern_count < counts_dict[current_pattern]):
			max_pattern_count = counts_dict[current_pattern]

	return list(pattern for pattern, cnt in counts_dict.items() if cnt == max_pattern_count)

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
	word_len = int(dataset[1])
	frequent_words = frequent_word(chromosome, word_len)

	print("chromosome: %s" % chromosome)
	print("frequent words:")
	for word in frequent_words:
		print(word)

if __name__ == '__main__':
	main()
