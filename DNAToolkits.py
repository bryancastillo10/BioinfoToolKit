from biomolecules import *
import collections

def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotide:
            return False
    return tmpseq

def countNucFrequency(seq):
    tmpFreqDict = {"A":0 ,"C":0,"T":0,"G":0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict
    # return dict(collections.Counter(seq))

def transcription(seq):
    # DNA -> RNA Transcription
    return seq.replace("T","U")

def reverse_complement(seq):
    return ''.join([DNA_ReverseComplement[nuc] for nuc in seq])[::-1]
    # mapping = str.maketrans('ATCG','TAGC')
    # return seq.translate(mapping)[::-1]

def gc_content(seq):
    """ GC content in a DNA/RNA Sequence"""
    return round(seq.count('C')+seq.count('G')/len(seq) * 100)

def gc_content_subseq(seq, k=20):
    """GC content in a DNA/RNA sub-sequence length k, k = 20 by default"""
    res = []
    for i in range(0, len(seq) - k + 1,k):
        subseq = seq[i:i+k]
        res.append(gc_content(subseq))
    return res

def translate_seq(seq, init_pos=0):
    """ Translates a DNA sequence into an amino acid sequence"""
    return [DNA_Codons[seq[pos:pos + 3]] for pos in range(init_pos, len(seq) - 2, 3)]

def codon_usage(seq,aminoacid):
    tmplist = []
    for i in range(0,len(seq) - 2, 3):
        if DNA_Codons[seq[i:i+3]] == aminoacid:
            tmplist.append(seq[i:i+3])
    
    freqDict = dict(collections.Counter(tmplist))
    totalWeight = sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq] = round(freqDict[seq]/totalWeight,2)
    return freqDict

    