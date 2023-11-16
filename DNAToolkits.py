Nucleotide = ["G","T","A","C"]
import collections

def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotide:
            return False
    return tmpseq

def countNucFrequency(seq):
    # tmpFreqDict = {"A":0 ,"C":0,"T":0,"G":0}
    # for nuc in seq:
    #     tmpFreqDict[nuc] += 1
    # return tmpFreqDict
    return dict(collections.Counter(seq))