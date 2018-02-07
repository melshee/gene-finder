from random import *

#database of codon table
insertion = False;
amino_acids = {
    'AUU': "I",
    "AUC": "I",
    "AUA": "I",

    "CUU": 'L',
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "UUA": "L",
    "UUG": "L",
    
    'GUU': 'V',
    'GUC': 'V',
    'GUA': 'V',
    'GUG': 'V',

    'UUU': 'F',
    'UUC': 'F',

    'AUG': 'M',

    'UGU': 'C',
    'UGC': 'C',

    'GCU': 'A',
    'GCC': 'A',
    'GCA': 'A',
    'GCG': 'A',

    'GGU': 'G',
    'GGC': 'G',
    'GGA': 'G',
    'GGG': 'G',

    'CCU': 'P',
    'CCC': 'P',
    'CCA': 'P',
    'CCG': 'P',

    'ACU': 'T',
    'ACC': 'T',
    'ACA': 'T',
    'ACG': 'T',

    'UCU': 'S',
    'UCC': 'S',
    'UCA': 'S',
    'UCG': 'S',
    'AGU': 'S',
    'AGC': 'S',

    'UAU': 'Y',
    'UAC': 'Y',

    'UGG': 'W',

    'CAA': 'Q',
    'CAG': 'Q',

    'AAU': 'N',
    'AAC': 'N',

    'CAU': 'H',
    'CAC': 'H',
    
    'GAA': 'E',
    'GAG': 'E',

    'GAU': 'D',
    'GAC': 'D',
    
    'AAA': 'K',
    'AAG': 'K',

    'CGU': 'R',
    'CGC': 'R',
    'CGA': 'R',
    'CGG': 'R',
    'AGA': 'R',
    'AGG': 'R',
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
    'R': 'arginine',
}

corresponding_nucleotides = {
    #purines
    'A': 'T', 
    'T': 'A',
    #pyradimines
    'G': 'C',
    'C': 'G',
}

while (True):
    #get the gene strand from user: test, TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAAATTTCATCCACT
    #this ones better: TCAATGCGCGCTACCCGGTAAAGCTCTGGGCCCAAATTTCATCCACT
    if insertion == False:
        template_strand = raw_input('gene sequence: ')
    else: 
        index = randint(1, find_stop_codon(True))
        template_strand = template_strand[:index] + mutation.upper() + template_strand[index:]
        #check this
        insertion = False
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
    try:
        template_strand_ORF = [template_strand_copy[i:i+3] for i in range(template_strand_copy.find("ATG"), find_stop_codon(True), 3)] #replace the open reading frame mrna 
    except:
        template_strand_ORF = [template_strand_copy[i:i+3] for i in range(template_strand_copy.find("ATG"), len(template_strand_copy)- 1), 3 ] #replace the open reading frame mrna
    #find the first aug, then split off into 3 from there. -- FIXED VERSION OF FINDING AN ORF 
    complementary_strand_copy = complementary_strand
    complementary_strand_ORF = []
    complementary_strand_ORF = [complementary_strand_copy[i:i+3] for i in range(complementary_strand_copy.find("ATG"), find_stop_codon(False), 3)] #replace the open reading frame mrna 

    mrna = ""

    #OPTIMIZED TEMPLATE -> MRNA 
    mrna = template_strand_ORF
    mrna = "".join(mrna)
    mrna = mrna.replace('T', 'U')

    n = 3
    mrna_codons = [mrna[i:i+n] for i in range(0, len(mrna), n)] #replace the open reading frame mrna 
    protein_sequence = ""
    #translation -- convert mrna strand to protein sequence 
    #testing sequence: TCATAACGCGCTACCCGATGGAGCTCTGGGCCCAAATTTCATCCACT
    print(mrna_codons)
    for mrna_codon in mrna_codons:
        for codon_option in amino_acids:
                #first one -- nterminal, last one -- cterminal 
                if mrna_codon == codon_option:
                    protein_sequence += (amino_acids[codon_option] + " - ")

    print(protein_sequence)

    if len(template_strand_ORF) > len(complementary_strand_ORF):
        print("Template strand is longer")
    else:
        print("Complementary strand is longer")

    answer = raw_input('create mutation (press c) or input another gene sequence (press g): ')
    if (answer == 'c'):
        insertion = True
        mutation = raw_input("Enter a nucleotide (a, g, t, c) to randomly insert into DNA sequence")