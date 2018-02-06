#database of codon table? 
amino_acids = {
    'I': {'AUU', 'AUC', 'AUA'},
    'L': {'CUU', 'CUC', 'CUA', 'CUG', 'UUA', 'UUG'},
    'V': {'GUU', 'GUC', 'GUA', 'GUG'},
    'F': {'UUU', 'UUC'},
    'M': {'AUG'},
    'C': {'UGU', 'UGC'},
    'A': {'GCU', 'GCC', 'GCA', 'GCG'},
    'G': {'GGU', 'GGC', 'GGA', 'GGG'},
    'P': {'CCU', 'CCC', 'CCA', 'CCG'},
    'T': {'ACU', 'ACC', 'ACA', 'ACG'},
    'S': {'UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'},
    'Y': {'UAU','UAC'},
    'W': {'UGG'},
    'Q': {'CAA','CAG'},
    'N': {'AAU','AAC'},
    'H': {'CAU','CAC'},
    'E': {'GAA','GAG'},
    'D': {'GAU','GAC'},
    'K': {'AAA','AAG'},
    'R': {'CGU', 'CGC','CGA','CGG', 'AGA', 'AGG'},
    #'STOP': {'UAA', }
}

def find_stop_codon(a):
    if a == True:
        list = {template_strand_copy.find("TAA"), template_strand_copy.find("TAG"), template_strand_copy.find("TGA")}
    else:
        list = {complementary_strand_copy.find("TAA"), complementary_strand_copy.find("TAG"), complementary_strand_copy.find("TGA")}
    for x in list: 
        if x > -1:
            return x

amino_abbr = {
    'I': 'isoleucine',
    'L': 'leucine',
    'V': 'valine',
    'F': 'phenylalanine',
    'M': 'methionine',
    'C': 'cysteine',
    'A': 'alanine',
    'G': 'glycine',
    'P': 'proline',
    'T': 'threonine',
    'S': 'serine',
    'Y': 'tyrosine',
    'Q': 'tryptophan',
    'W': 'asparagine',
    'Q': 'glutamine',
    'N': 'asparagine',
    'H': 'histidine',
    'E': 'glutamic acid',
    'D': 'aspartic acid',
    'K': 'lysine',
    'R': 'arginine'
}

corresponding_nucleotides = {
    #purines
    'A': 'T', 
    'T': 'A',
    #pyradimines
    'G': 'C',
    'C': 'G'
}

#get the gene strand from user: test, TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAAATTTCATCCACT
#this ones better: TCAATGCGCGCTACCCGGTAAAGCTCTGGGCCCAAATTTCATCCACT
template_strand = raw_input('gene sequence: ')
complementary_strand = ""

#transcription -- forming complementary strand 
for nucleic_acid in template_strand:
    if (nucleic_acid not in list(corresponding_nucleotides.keys())):
        raise Exception(nucleic_acid + " is not a valid nucleotide")
    complementary_strand += corresponding_nucleotides[nucleic_acid]


#find largest open reading frame from dna sequence 
orf_length_template = 0
template_strand_copy = template_strand
template_strand_ORF = []

#find the first aug, then split off into 3 from there. -- FIXED VERSION OF FINDING AN ORF 
template_strand_ORF = [template_strand_copy[i:i+3] for i in range(template_strand_copy.find("ATG"), find_stop_codon(True), 3)] #replace the open reading frame mrna 

#find the first aug, then split off into 3 from there. -- FIXED VERSION OF FINDING AN ORF 
complementary_strand_copy = complementary_strand
complementary_strand_ORF = []
complementary_strand_ORF = [complementary_strand_copy[i:i+3] for i in range(complementary_strand_copy.find("ATG"), find_stop_codon(False), 3)] #replace the open reading frame mrna 

mrna = ""
#complementary_strand_ORF = str(complementary_strand_ORF
#go through array of complementary_strand_ORF and for each key, replace the thing w the corresponding nucelotide
'''
for nucleic_acid in template_strand_ORF:
    for nucelotide in nucleic_acid:
        if nucelotide == 'T':
            mrna += 'U'
        else:
            mrna += nucelotide

print("mrna" + str(mrna))
'''

#OPTIMIZED TEMPLATE -> MRNA 
mrna = template_strand_ORF
mrna = "".join(mrna)
mrna = mrna.replace('T', 'U')
#print("mrna" + str(mrna))

'''
else:
    mrna = ""
    #complementary_strand_ORF = str(complementary_strand_ORF)

    #go through array of complementary_strand_ORF and for each key, replace the thing w the corresponding nucelotide
    for nucleic_acid in complementary_strand_ORF: 
        for nucelotide in nucleic_acid:
            if nucelotide == 'A':
                mrna += 'U'
            else:
                mrna += corresponding_nucleotides[nucelotide]
'''

#previous method for finding the len of open reading frame, doesn't work because it won't find in increment of 3's, for ex it'll consider #attaaa having a stop, even tho its att and aaa, neither of which are stops but tta is.... this is bad 
# while template_strand_copy.find("ATG") > -1:
#     if (template_strand_copy.find("ATG") > -1 and template_strand_copy.find("TAA") > -1): #atg - start, taa - stop
#         diff = template_strand_copy.find("TAA") - template_strand_copy.find("ATG")
#         diff /= 3 #3 codons = one amino acid
#         if diff > orf_length_template: #orf_len in terms of # amino acids
#             orf_length_template = diff
#         template_strand_copy = template_strand_copy[template_strand_copy.find("ATG") + 3:]

#find largest open reading frame from dna sequence (complementary to template)
#|TCA,ATG,CGC,GCT,ACC,CGG,TAA,AGC,TCT,GGG,CCC,AAA,TTT,CAT,CCA,CT
#|AGT,TAC,GCG,CGA,TGG,GCC,ATT,TCG,AGA,CCC,GGG,TTT,AAA,GTA,GGT,GA
# orf_length_complementary = 0
# complementary_strand_copy = complementary_strand
# while complementary_strand_copy.find("ATG") > -1:
#     if (complementary_strand_copy.find("ATG") > -1 and complementary_strand_copy.find("TAA") > -1): #atg - start, taa - stop
#         diff = complementary_strand_copy.find("TAA") - complementary_strand_copy.find("ATG")
#         diff /= 3 #3 codons = one amino acid
#         if diff > orf_length_complementary: #orf_len in terms of # amino acids
#             orf_length_complementary = diff
#         complementary_strand_copy = complementary_strand_copy[complementary_strand_copy.find("ATG") + 3:]

# if orf_length_template > orf_length_complementary:
#     print("yea")
# else:
#     print("nah")

#transcription -- forming mrna from template strand

n = 3
mrna_codons = [mrna[i:i+n] for i in range(0, len(mrna), n)] #replace the open reading frame mrna 
protein_sequence = ""
#translation -- convert mrna strand to protein sequence 
#testing sequence: TCATAACGCGCTACCCGATGGAGCTCTGGGCCCAAATTTCATCCACT
print(mrna_codons)
for mrna_codon in mrna_codons:
    for amino_acid in amino_acids:
        for codon_option in amino_acids[amino_acid]:
            #first one -- nterminal, last one -- cterminal 
            if mrna_codon == codon_option:
                protein_sequence += (amino_acid + " - ")

print(protein_sequence)

if len(template_strand_ORF) > len(complementary_strand_ORF):
    print("Template strand is longer")
else:
    print("Complementary strand is longer")
