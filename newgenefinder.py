from random import *

from dictionaries import amino_acids
from dictionaries import amino_abbr
from dictionaries import corresponding_nucleotides


def main():
    #gene_arr = ["","","","","",""] #array of size 6 to hold all 6 possible genes  
    #input strand =  TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAAATTTCATCCACT

    #get the gene strand from user: test, TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAAATTTCATCCACT
    #this ones better: TCAATGCGCGCTACCCGGTAAAGCTCTGGGCCCAAATTTCATCCACT


    #assume: given strand is the template strand (coding strand)
      #read dna strand (read all 3 ORFs)
          # gene_arr[0] = read(input_dna_strand)
          # gene_arr[1] = read(input_dna_strand - first index)
          # gene_arr[2] = read(input_dna_strand - second index)
    #assume: given strand is the complementary strand (negative of the coding strand)
      #negate(input_strand) #method to negate dna strand 
      #read dna strand (read all 3 ORFs)
          # gene_arr[3] = read(input_dna_strand)
          # gene_arr[4] = read(input_dna_strand - first index)
          # gene_arr[5] = read(input_dna_strand - second index)

    # longest_dna_strand = #the longest dna strand out of the 6 in gene_arr

    #longest_mrna_strand = toMRNA(longest_dna_strand)
    #longest_amino_acid = toAminoAcid(longest_mrna_strand)
    #print(longest_amino_acid)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# methods to implement:
# read(input_dna_strand) #returns length of LONGEST gene for that ORF (there might be multiple gene encodings on the SAME strand on the SAME ORF)
# toMRNA(longest_dna_strand) #converts input_dna_strand into an mRNA strand
# toAminoAcid(longest_mrna_strand) #converts input_mrna_strand into an amino acid sequence (one letter sequence)






# insertion = False;


# def find_stop_codon(a):
#     if a == True:
#         list = {template_strand_copy.find("TAA"), template_strand_copy.find("TAG"), template_strand_copy.find("TGA")}
#     else:
#         list = {complementary_strand_copy.find("TAA"), complementary_strand_copy.find("TAG"), complementary_strand_copy.find("TGA")}
#     for x in list: 
#         if x > -1:
#             return x

# while (True):
#     #get the gene strand from user: test, TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAAATTTCATCCACT
#     #this ones better: TCAATGCGCGCTACCCGGTAAAGCTCTGGGCCCAAATTTCATCCACT
#     if insertion == False:
#         template_strand = raw_input('gene sequence: ')
#     else: 
#         index = randint(1, find_stop_codon(True))
#         template_strand = template_strand[:index] + mutation.upper() + template_strand[index:]
#         #check this
#         insertion = False
#     complementary_strand = ""

#     #transcription -- forming complementary strand 
#     for nucleic_acid in template_strand:
#         if (nucleic_acid not in list(corresponding_nucleotides.keys())):
#             raise Exception(nucleic_acid + " is not a valid nucleotide")
#         complementary_strand += corresponding_nucleotides[nucleic_acid]

#     #find largest open reading frame from dna sequence 
#     orf_length_template = 0
#     template_strand_copy = template_strand
#     template_strand_ORF = []

#     #find the first aug, then split off into 3 from there. -- FIXED VERSION OF FINDING AN ORF 
#     try:
#         template_strand_ORF = [template_strand_copy[i:i+3] for i in range(template_strand_copy.find("ATG"), find_stop_codon(True), 3)] #replace the open reading frame mrna 
#     except:
#         template_strand_ORF = [template_strand_copy[i:i+3] for i in range(template_strand_copy.find("ATG"), len(template_strand_copy)- 1), 3 ] #replace the open reading frame mrna
#     #find the first aug, then split off into 3 from there. -- FIXED VERSION OF FINDING AN ORF 
#     complementary_strand_copy = complementary_strand
#     complementary_strand_ORF = []
#     complementary_strand_ORF = [complementary_strand_copy[i:i+3] for i in range(complementary_strand_copy.find("ATG"), find_stop_codon(False), 3)] #replace the open reading frame mrna 

#     mrna = ""

#     #OPTIMIZED TEMPLATE -> MRNA 
#     mrna = template_strand_ORF
#     mrna = "".join(mrna)
#     mrna = mrna.replace('T', 'U')

#     n = 3
#     mrna_codons = [mrna[i:i+n] for i in range(0, len(mrna), n)] #replace the open reading frame mrna 
#     protein_sequence = ""
#     #translation -- convert mrna strand to protein sequence 
#     #testing sequence: TCATAACGCGCTACCCGATGGAGCTCTGGGCCCAAATTTCATCCACT
#     print(mrna_codons)
#     for mrna_codon in mrna_codons:
#         for codon_option in amino_acids:
#                 #first one -- nterminal, last one -- cterminal 
#                 if mrna_codon == codon_option:
#                     protein_sequence += (amino_acids[codon_option] + " - ")

#     print(protein_sequence)

#     if len(template_strand_ORF) > len(complementary_strand_ORF):
#         print("Template strand is longer")
#     else:
#         print("Complementary strand is longer")

#     answer = raw_input('create mutation (press c) or input another gene sequence (press g): ')
#     if (answer == 'c'):
#         insertion = True
#         mutation = raw_input("Enter a nucleotide (a, g, t, c) to randomly insert into DNA sequence")