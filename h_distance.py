dna_str_1 = "TTCGATCCATTG"
dna_str_2= "ATCAATCGATCG"

#Loop Approach
def h_d_loop(str_1,str_2):
    h_distance = 0
    for position in range(len(str_1)):
        if str_1[position] != str_2[position]:
            h_distance += 1
    return h_distance

#Set Approach
def h_d_set(str_1,str_2):
    nucleotide_set_1 = set([(x,y) for x, y in enumerate(str_1)])
    nucleotide_set_2 = set([(x,y) for x, y in enumerate(str_2)])
    # for x in range(len(nucleotide_set_1)):
    #     print(sorted(nucleotide_set_1)[x],sorted(nucleotide_set_2)[x])
    return len(nucleotide_set_1.difference(nucleotide_set_2))

#Zip Approach
def h_d_zip(str_1,str_2):
    return len([n1,n2] for n1, n2 in zip(str_1,str_2) if n1 != n2)


#Testing each functions
print("Loop Hamming Distance:", end='')
print(h_d_loop(dna_str_1,dna_str_2))

print("Set Hamming Distance:", end='')
print(h_d_set(dna_str_1,dna_str_2))

print("Zip Hamming Distance:", end='')
print(h_d_set(dna_str_1,dna_str_2))