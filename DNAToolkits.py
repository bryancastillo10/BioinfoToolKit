from biomolecules import *
import collections


def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotide:
            return False
    return tmpseq


def countNucFrequency(seq):
    tmpFreqDict = {"A": 0, "C": 0, "T": 0, "G": 0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict
    # return dict(collections.Counter(seq))


def transcription(seq):
    # DNA -> RNA Transcription
    return seq.replace("T", "U")


def reverse_complement(seq):
    return "".join([DNA_ReverseComplement[nuc] for nuc in seq])[::-1]
    # mapping = str.maketrans('ATCG','TAGC')
    # return seq.translate(mapping)[::-1]


def gc_content(seq):
    """GC content in a DNA/RNA Sequence"""
    return round(seq.count("C") + seq.count("G") / len(seq) * 100)


def gc_content_subseq(seq, k=20):
    """GC content in a DNA/RNA sub-sequence length k, k = 20 by default"""
    res = []
    for i in range(0, len(seq) - k + 1, k):
        subseq = seq[i : i + k]
        res.append(gc_content(subseq))
    return res


def translate_seq(seq, init_pos=0):
    """Translates a DNA sequence into an amino acid sequence"""
    return [DNA_Codons[seq[pos : pos + 3]] for pos in range(init_pos, len(seq) - 2, 3)]


def codon_usage(seq, aminoacid):
    tmplist = []
    for i in range(0, len(seq) - 2, 3):
        if DNA_Codons[seq[i : i + 3]] == aminoacid:
            tmplist.append(seq[i : i + 3])

    freqDict = dict(collections.Counter(tmplist))
    totalWeight = sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq] = round(freqDict[seq] / totalWeight, 2)
    return freqDict


def gen_reading_frames(seq):
    """Generate the 6 reading frames of a DNA sequence, including the reverse complement"""
    frames = []
    frames.append(translate_seq(seq, 0))
    frames.append(translate_seq(seq, 1))
    frames.append(translate_seq(seq, 2))
    frames.append(translate_seq(reverse_complement(seq), 0))
    frames.append(translate_seq(reverse_complement(seq), 1))
    frames.append(translate_seq(reverse_complement(seq), 2))
    return frames


def proteins_from_rf(aa_seq):
    """Returns a list of all possible proteins based on the amino acid sequence"""
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":  # Stop Codon
            for p in current_prot:
                proteins.append(p)
            current_prot = []
        else:
            if aa == "M":  # Start Codon
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins


# 1) Generate all Reading Frames
# 2) Extract all proteins
# 3) Return a sorted/unsorted list
def all_proteins_from_orfs(seq, startReadPos=0, endReadPos=0, ordered=False):
    """Compute all possible proteins for all open reading frames"""
    """ Protine Search DB: https://wwww.ncbi.nlm.nih.gov/nuccore/NM_001185097.2"""
    """ API can be used to pull protein info """
    if endReadPos > startReadPos:
        rfs = gen_reading_frames(seq[startRead:endRead])
    else:
        rfs = gen_reading_frames(seq)

    res = []
    for rf in rfs:
        prots = proteins_from_rf(rf)
        for p in prots:
            res.append(p)
    if ordered:
        return sorted(res, key=len, reverse=True)
    return res
