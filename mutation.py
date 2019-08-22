import random
 
from terminaltables import AsciiTable
 
 
class Gen():
    def __init__(self, epoch, forward, reverse):
        self.epoch = epoch
        self.forward = forward
        self.reverse = reverse
 
    @property
    def total(self):
        return {
            'A': self.forward['A'] + self.reverse['A'],
            'T': self.forward['T'] + self.reverse['T'],
            'C': self.forward['C'] + self.reverse['C'],
            'G': self.forward['G'] + self.reverse['G'],
        }
 
    def mutate(self):
        mutate_forward = self.forward.copy()
        mutate_reverse = self.reverse.copy()
 
        # Some C'c get changed into T's on the forward half strand
        number_of_mutations = round(random.randint(3, 6) / 100.0 * mutate_forward['C'])
        mutate_forward['C'] -= number_of_mutations
        mutate_forward['T'] += number_of_mutations
 
        return Gen(
            self.epoch,
            mutate_forward,
            mutate_reverse
        )
 
    def simple_copy(self):
        def complementary(half_strand):
            return {
                'A': half_strand['T'],
                'T': half_strand['A'],
                'C': half_strand['G'],
                'G': half_strand['C'],
            }
 
        return Gen(
            self.epoch + 1,
            complementary(self.reverse),
            complementary(self.forward),
        )
 
    def copy(self):
        return self.mutate().simple_copy()
 
    def __str__(self):
        table_data = [
            ['', 'C', 'G', 'A', 'T'],
            ['total', str(self.total['C']), str(self.total['G']), str(self.total['A']), str(self.total['T'])],
            ['reverse half strand', str(self.reverse['C']), str(self.reverse['G']), str(self.reverse['A']),
             str(self.reverse['T'])],
            ['forward half strand', str(self.forward['C']), str(self.forward['G']), str(self.forward['A']),
             str(self.forward['T'])]
        ]
        table = AsciiTable(table_data)
        return "Generation #{}\n".format(self.epoch) + table.table
 
 
test = Gen(0, {'A': 10000, 'T': 10000, 'C': 10000, 'G': 10000}, {'A': 10000, 'T': 10000, 'C': 10000, 'G': 10000})
for i in range(0, 10):
    test = test.copy()
    print(test)