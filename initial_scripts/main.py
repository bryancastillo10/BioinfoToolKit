from DNAToolkits import *
from utilities.colored import color
import random

randDNAstr = "".join([random.choice(Nucleotide) for nuc in range(50)])

DNAStr = validateSeq(randDNAstr)
print(f"\n Sequence: {color(DNAStr)}\n")
print(f"[1]  Sequence Length: {len(DNAStr)}\n")
print(color(f"[2]  Nucleotide Frequency: {countNucFrequency(DNAStr)})\n"))
print(f"[3]  DNA to RNA Transcription: {color(transcription(DNAStr))}\n")

print(f"[4]  DNA String + Reverse Complement: \n5' {color(DNAStr)} 3'")
print(f"   {''.join(['|' for i in range(len(DNAStr))])}")
print(f"3' {color(reverse_complement(DNAStr)[::-1])} 5' [Complementary DNA]")
print(f"5' {color(reverse_complement(DNAStr))} 3' [Reverse Complementary]")

print(f"[5]  GC Content: {gc_content(DNAStr)}% \n")

print(f"[6]  GC Content in Subsection k=5: {gc_content_subseq(DNAStr, k=5)} \n")

print(f"[7]  Amino Acid Sequence from DNA: {translate_seq(DNAStr)} \n")

print(f"[8]  Codon Frequency (L): {codon_usage(DNAStr, 'L')} \n")

print(f"[9]  Reading Frame: ")
for frame in gen_reading_frames(DNAStr):
    print(frame)

print("\n [10] All Proteins in the 6 Open Reading Frames:")
for prot in all_proteins_from_orfs(NM_000207_3, 0, 0, True):
    print(f"{prot}")
