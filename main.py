from DNAToolkits import *
from utilities.colored import color
import random

randDNAstr = ''.join([random.choice(Nucleotide)
        for nuc in range(60)])

DNAStr = validateSeq(randDNAstr)
print(f'\n Sequence: {color(DNAStr)}\n')
print(f'[1]  Sequence Length: {len(DNAStr)}\n')
print(color(f'[2]  Nucleotide Frequency: {countNucFrequency(DNAStr)})\n'))
print(f'[3]  DNA to RNA Transcription: {color(transcription(DNAStr))}\n')

print(f"[4]  DNA String + Reverse Complement: \n5' {color(DNAStr)} 3'")
print(f"   {''.join(['|' for i in range(len(DNAStr))])}")
print(f"3' {color(reverse_complement(DNAStr))} 5'\n")