from bio_seq import bio_seq

test_dna = bio_seq("ATCGGTTTG","DNA","TestLab")

print(test_dna.get_seq_info())
# print(test_dna.get_seq_biotype())

test_dna.generate_rnd_seq()
print(test_dna.get_seq_info())