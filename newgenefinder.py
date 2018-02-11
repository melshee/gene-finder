from random import *
from dictionaries import amino_acids
from dictionaries import amino_abbr
from dictionaries import corresponding_nucleotides

#returns a list of indeces that substring (sub) is found at
def find_start_codons(dna_strand, sub="ATG"): 
  listindex = []
  i = dna_strand.find(sub)
  while i >= 0:
    listindex.append(i)
    i = dna_strand.find(sub, i + 1)
  return listindex

def clip_to_stop_codon(dna_strand):
  stop_codon_list = {dna_strand.find("TAA"), dna_strand.find("TAG"), dna_strand.find("TGA")}
  stop_codon_list.remove(-1)
  # stop_codon_list.remove(-1)
  first_stop_codon = min(stop_codon_list)
  first_stop_codon = first_stop_codon + 3 #to keep the 3 bases in the stop codon
  # print dna_strand[:first_stop_codon]
  return dna_strand[:first_stop_codon]

def read(input_dna):
  print input_dna
  start_codons_list = find_start_codons(input_dna)
  # print start_codons_list
  strands = [] #all strands from start codon to end of strand
  genes = [] #all strands from start codon to closest stop codon
  for start_codon_index in start_codons_list:
    strands.append(input_dna[start_codon_index:])
  print strands #list of all genes (each start with ATG)
  for strand in strands:
    gene = clip_to_stop_codon(strand)
    genes.append(gene)
  print genes
  longest = max(genes, key=len)
  print longest
  print len(longest)
  return longest



def main():
    print 'in main'

    gene_arr = ["","","","","",""] #array of size 6 to hold all 6 possible genes
    # print gene_arr  
    input_dna_strand =  "CCCATGCCCCCCCATGCCCCCCTGACCCCCATGCCCCTGA"
    # input_dna_strand =  "TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAAATTTCATCCACT"
    #get the gene strand from user: test, TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAAATTTCATCCACT
    #this ones better: TCAATGCGCGCTACCCGGTAAAGCTCTGGGCCCAAATTTCATCCACT

    #assummption 1: given strand is the template strand (coding strand)
      #read dna strand (read all 3 ORFs)
    orf1 = input_dna_strand
    gene_arr[0] = read(orf1)
    orf2 = input_dna_strand[1:]
    gene_arr[1] = read(orf2)
    orf3 = input_dna_strand[2:]
    gene_arr[2] = read(orf3)
    #assumption 2: given strand is the complementary strand (negative of the coding strand)
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

main()




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