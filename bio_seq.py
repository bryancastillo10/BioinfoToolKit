from bio_struct import *
import random
class bio_seq:
    """ DNA sequence class. Default values: ATCG, DNA, None"""
    def __init__(self, seq="ATCG",seq_type="DNA",label='None'):
        """ Sequence initialization, validation"""
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided data doesn't seem to be a correct {self.seq_type} sequence."
    
    
    def __validate(self):
        """ Check the input sequences to verify if it is a valid DNA"""
        return set(Nucleotide).issuperset(self.seq)
    
    def get_seq_biotype(self):
        """Returns a sequence type """
        return self.seq_type

    def get_seq_info(self):
        """ Returns four strings. Full sequence information"""
        return f"[Label]: {self.label}\n [Sequence] {self.seq} \n [Biomolecule Type] {self.seq_type}\n [Sequence Length]{len(self.seq)}"

    def generate_rnd_seq(self,length=10,seq_type="DNA"):
        """ Generate a random DNA sequence, provided the length"""
        seq = ''.join([random.choice(Nucleotide)
            for x in range(length)])
        self.__init__(seq,seq_type,"Randomly generated sequence")