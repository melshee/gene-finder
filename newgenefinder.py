from random import *
from dictionaries import amino_acids
from dictionaries import amino_abbr
from dictionaries import corresponding_nucleotides
import textwrap #to separate the strands into groups of 3
import sys #for terminating program 

#terminates program, displaying final results
def terminateprogram(results):
  final_results = ""
  final_results += "         DNA sequence ---------------------- | "
  final_results += "Amino Acid\n"
  roundnum = 1
  for dna_sequence in results:
    final_results += ("Round " + str(roundnum) + ": " + results[dna_sequence][0] + " | " + results[dna_sequence][1] + "\n")
    roundnum += 1
  print(final_results)

  # for dna_sequence in results:
  #   final_results += ("" + dna_sequence + " | " + results[dna_sequence] + "\n")
  # print(final_results)
  sys.exit()

#mutates the dna sequence 
def mutateGene(dna_sequence, num_nucleotides_changed):
  i = 0 
  while i <= num_nucleotides_changed:
    index = randint(1, len(dna_sequence) - 1)
    random_nucleotide = corresponding_nucleotides.keys()[randint(1, 3)]
    dna_sequence = dna_sequence[:index - 1] + random_nucleotide + dna_sequence[index:]
    i = i + 1
  return dna_sequence

#returns the MRNA strand of the longest ORF 
def toMRNA(longest_dna_strand): #will always read from 5' to 3', whether it be the complementary or template strand.
  mrna = "".join(longest_dna_strand)
  mrna = mrna.replace('T', 'U')
  return mrna

def toAminoAcid(longest_mrna_strand):
  amino_acid_sequence = ""
  mrna_codons = [longest_mrna_strand[i:i+3] for i in range(0, len(longest_mrna_strand), 3)] #replace the open reading frame mrna 
  for codon in mrna_codons:
    amino_acid_sequence += (amino_acids[codon] + ' ')
  return amino_acid_sequence

#returns a list of indeces that substring (sub) is found at
def find_start_codons(dna_strand, sub="ATG"):
  listindex = []
  codon_array = textwrap.wrap(dna_strand, 3)
  print codon_array
  i = codon_array.index("ATG") if "ATG" in codon_array else -1
  while i >= 0:
    listindex.append(i*3)
    codon_array = codon_array[i+1:]
    i = codon_array.index("ATG") if "ATG" in codon_array else -1
  # print listindex
  return listindex 

#returns the dna_strand minus all the base pairs following the first stop codon found
def clip_after_stop_codon(dna_strand):
  codon_array = textwrap.wrap(dna_strand, 3)
  # print codon_array
  #THIS LINE BELOW IS UGLY NEED HALP MAKING IT BETTER LOL
  stop_codon_list = [codon_array.index("TAA") if "TAA" in codon_array else -1, 
             codon_array.index("TAG") if "TAG" in codon_array else -1,
             codon_array.index("TGA") if "TGA" in codon_array else -1]
 #more ugliness I try to fix a better way somehow..later
  stop_codon_list.remove(-1) if -1 in stop_codon_list else None
  stop_codon_list.remove(-1) if -1 in stop_codon_list else None
  stop_codon_list.remove(-1) if -1 in stop_codon_list else None
  first_stop_codon = 0 #if no stop codon is found return an empty string (no gene possible/exists)
  if len(stop_codon_list) != 0:
    first_stop_codon = min(stop_codon_list) 
    first_stop_codon = first_stop_codon * 3
    first_stop_codon = first_stop_codon + 3  #to keep the 3 bases in the stop codon
  print dna_strand[:first_stop_codon]
  return dna_strand[:first_stop_codon]

def read(input_dna):
  longest = ""
  start_codons_list = find_start_codons(input_dna)
  if len(start_codons_list) != 0: #if ATG codon was found
    strands = [] #all strands from start codon to end of strand
    genes = [] #all strands from start codon to closest stop codon
    for start_codon_index in start_codons_list:
      strands.append(input_dna[start_codon_index:])
    # print strands #list of all genes (each start with ATG)
    for strand in strands:
      gene = clip_after_stop_codon(strand)
      genes.append(gene)

    longest = max(genes, key=len)
  return longest

def find_complement(input_dna):
  complementary_dna = ""
  for nucleic_acid in input_dna:
    if (nucleic_acid not in list(corresponding_nucleotides.keys())):
      raise Exception(nucleic_acid + " is not a valid nucleotide")
    complementary_dna += corresponding_nucleotides[nucleic_acid]
  # print "complementay dna = " + complementary_dna
  print "complementary (3' to 5') " + complementary_dna
  complementary_dna = complementary_dna[::-1]
  print "complementary (5' to 3') " + complementary_dna #reverse strand to read from 5' to 3'
  return complementary_dna

def main():
  rounds = 1
  input_dna_strand = "ATGAAACTATGATAAAAAATTACCCCCCCCCCTAA"
  results = {}
  while True:
    gene_arr = ["","","","","",""] #array of size 6 to hold all 6 possible genes
    print "\n---------------------- ROUND " + str(rounds) + " -------------------------"
    print "\noriginal strand = " + input_dna_strand
    # input_dna_strand = "ATGCCCCTAATGCTAAAAATTCAATAAAATAGAAATAA" #testing stop codon wit diff ORFs  
    # input_dna_strand =  "CCCATGCCCCCCCATGCCCCCCTGACCCCCATGCCCCTGA" #mel's ex on Sat
    # input_dna_strand =  "TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAAATTTCATCCACT"
    #get the gene strand from user: test, TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAAATTTCATCCACT
    #this ones better: TCAATGCGCGCTACCCGGTAAAGCTCTGGGCCCAAATTTCATCCACT

    #assummption 1: given strand is the template strand (coding strand). Read dna strand (read all 3 ORFs)
    print "~~~~~~ORF 1~~~~~~"
    orf1 = input_dna_strand
    gene_arr[0] = read(orf1)
    print "~~~~~~ORF 2~~~~~~"
    orf2 = input_dna_strand[1:]
    gene_arr[1] = read(orf2)
    print "~~~~~~ORF 3~~~~~~"
    orf3 = input_dna_strand[2:]
    gene_arr[2] = read(orf3)
    print
    
    #assumption 2: given strand is the complementary strand (negative of the coding strand). Read complement of dna strand (read all 3 ORFs)
    comp_input_dna_strand = find_complement(input_dna_strand) #method to negate dna strand 
    print "~~~~~~ORF 4~~~~~~"
    comp_orf1 = comp_input_dna_strand
    gene_arr[3] = read(comp_orf1)
    print "~~~~~~ORF 5~~~~~~"
    comp_orf2 = comp_input_dna_strand[1:]
    gene_arr[4] = read(comp_orf2)
    print "~~~~~~ORF 6~~~~~~"
    comp_orf3 = comp_input_dna_strand[2:]
    gene_arr[5] = read(comp_orf3)
    print
    print "open reading frames 1-6 (line below): "
    print gene_arr

    longest_orf = max(gene_arr, key=len)
    print("\nThe longest ORF is: " + longest_orf)
    longest_mrna_strand = toMRNA(longest_orf)
    longest_amino_acid = toAminoAcid(longest_mrna_strand)
    print("The longest ORF converted to protein is: " + longest_amino_acid)

    results[rounds] = [input_dna_strand, longest_amino_acid]

    answer = raw_input("Continue or terminate program? (c or t): ")
    if (answer == "t"):
      terminateprogram(results)

    answer = raw_input('Would you like to mutate the sequence? (y, n) ')
    if (answer == 'y'):
        num_nucleotides_changed = raw_input('How many nucleotides would you like to see changed? (1 - ' + str(len(input_dna_strand)) + '): ')
        print("og strand ------" + input_dna_strand)
        input_dna_strand = mutateGene(input_dna_strand, int(num_nucleotides_changed))
        print("mutated strand--" + input_dna_strand)
    
    rounds = rounds + 1
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