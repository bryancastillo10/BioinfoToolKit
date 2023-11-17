from bio_struct import *
from collections import Counter
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
    
    #DNA Toolkit Section
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
    
    def nucleotide_frequency(self):
        """ Counts the nucleotides in a given sequence. Returns a dictionary"""
        return dict(Counter(self.seq))

    def transcription(self):
        """ DNA -> RNA Transcription. Replacing Thymine with Uracil"""
        if self.seq_type == "DNA":
            return self.seq.replace("T","U")
        return "Not a DNA sequence"
    
    def reverse_complement(self):
        """ 
        Swapping adenine with thymine and guanine with cytosine.
        Reversing the nucleotide sequence to generate another string
        """
        mapping = str.maketrans('ATCG','TAGC')
        return self.seq.translate(mapping)[::-1]
    
    def gc_content(self):
        """ GC Content in a DNA/RNA sequence"""
        return round(((self.seq.count('C') + self.seq.count('G')) / len(self.seq)) * 100)

    def gc_content_subseq(self, k=20):
        """GC content in a DNA/RNA sub-sequence length k, k = 20 by default"""
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i : i + k]
            res.append(
                round(((subseq.count('C') + subseq.count('G')) / len(subseq)) * 100)

            )
        return res
    
    def translate_seq(self, init_pos=0):
        """Translates a DNA sequence into an amino acid sequence"""
        return [DNA_Codons[self.seq[pos : pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]

    def codon_usage(self, aminoacid):
        tmplist = []
        for i in range(0, len(self.seq) - 2, 3):
            if DNA_Codons[self.seq[i : i + 3]] == aminoacid:
                tmplist.append(self.seq[i : i + 3])

        freqDict = dict(Counter(tmplist))
        totalWeight = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq] / totalWeight, 2)
        return freqDict

    def gen_reading_frames(self):
        """Generate the six reading frames of a DNA sequence, including the reverse complement"""
        frames = []
        frames.append(self.translate_seq(0))
        frames.append(self.translate_seq(1))
        frames.append(self.translate_seq(2))
        tmp_seq = bio_seq(self.reverse_complement(), self.seq_type)
        frames.append(tmp_seq.translate_seq(0))
        frames.append(tmp_seq.translate_seq(1))
        frames.append(tmp_seq.translate_seq(2))
        del tmp_seq
        return frames

    def proteins_from_rf(self,aa_seq):
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
    
    def all_proteins_from_orfs(self, startReadPos=0, endReadPos=0, ordered=False):
        """Compute all possible proteins for all open reading frames"""
        """ Protein Search DB: https://wwww.ncbi.nlm.nih.gov/nuccore/NM_001185097.2"""
        """ API can be used to pull protein info """
        if endReadPos > startReadPos:
            tmp_seq = bio_seq(
                self.seq[startReadPos:endReadPos],self.seq_type)
            rfs = tmp_seq.gen_reading_frames()
        else:
            rfs = self.gen_reading_frames()

        res = []
        for rf in rfs:
            prots = self.proteins_from_rf(rf)
            for p in prots:
                res.append(p)
        if ordered:
            return sorted(res, key=len, reverse=True)
        return res
