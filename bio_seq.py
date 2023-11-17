from bio_struct import *
class bio_seq:
    """ DNA sequence class. Default value: ATCG, DNA, No Label"""
    def __init__(self, seq="ATCG",seq_type="DNA",label='No Label'):
        """ Sequence initialization, validation"""
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.validate()
        assert self.is_valid, f"Provided data doesn't seem to be correct {self.seq}"
    
    
    def validate(self):
        """ Check the sequence to verify if it is a valid DNA string"""
        return set(Nucleotide).issuperset(self.seq)