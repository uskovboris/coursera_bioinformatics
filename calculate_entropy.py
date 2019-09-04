import numpy
from dna_lib import matrix_entropy

entropy_val = matrix_entropy([
['T', 'C', 'G', 'G', 'G', 'G', 'G', 'T', 'T', 'T', 'T', 'T'],
['C', 'C', 'G', 'G', 'T', 'G', 'A', 'C', 'T', 'T', 'A', 'C'],
['A', 'C', 'G', 'G', 'G', 'G', 'A', 'T', 'T', 'T', 'T', 'C'],
['T', 'T', 'G', 'G', 'G', 'G', 'A', 'C', 'T', 'T', 'T', 'T'],
['A', 'A', 'G', 'G', 'G', 'G', 'A', 'C', 'T', 'T', 'C', 'C'],
['T', 'T', 'G', 'G', 'G', 'G', 'A', 'C', 'T', 'T', 'C', 'C'],
['T', 'C', 'G', 'G', 'G', 'G', 'A', 'T', 'T', 'C', 'A', 'T'],
['T', 'C', 'G', 'G', 'G', 'G', 'A', 'T', 'T', 'C', 'C', 'T'],
['T', 'A', 'G', 'G', 'G', 'G', 'A', 'A', 'C', 'T', 'A', 'C'],
['T', 'C', 'G', 'G', 'G', 'T', 'A', 'T', 'A', 'A', 'C', 'C']])

print(entropy_val)
print(numpy.around(entropy_val, decimals=4))
