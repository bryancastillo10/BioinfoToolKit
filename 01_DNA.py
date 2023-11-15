from DNATools import *
import random

rDNAstr = ''.join([random.choice(Nucleotide)
        for nuc in range(50)])

print(validateSeq(rDNAstr))
print(countNucFrequency(rDNAstr))