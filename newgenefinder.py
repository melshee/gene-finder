from random import * #used to generate random idices for mutations
from dictionaries import amino_acids # codons -> amino acids (abbrevations) database
from dictionaries import amino_abbr # abbreviations -> full amino acid name database
from dictionaries import corresponding_nucleotides # corresponding nucleotide (generating complementary strands) database
import textwrap # to separate the strands into groups of 3
import sys # to terminating program 

# Terminates program, displaying final results as a table
def terminateprogram(results):
  final_results = "\n"
  final_results += "         DNA sequence                                    | "
  final_results += "Number Mutations | Amino Acid Sequence \n"             
  roundnum = 1
  for dna_sequence in results:
    final_results += ("Round " + str(roundnum) + ": " + results[dna_sequence][0] + " |        " + str(results[dna_sequence][2]) + "         | " + results[dna_sequence][1] + "\n")
    roundnum += 1
  print(final_results)
  sys.exit()

# Mutation - Deletion
def mutateGeneDeletion(dna_sequence):
  index = randint(1, len(dna_sequence) - 1)
  temp_dna_sequence = list(dna_sequence)
  temp_dna_sequence[index] = ''
  dna_sequence = "".join(str(x) for x in temp_dna_sequence)
  return dna_sequence

# Mutation - Insertion
def mutateGeneInsertion(dna_sequence):
  index = randint(1, len(dna_sequence) - 1)
  random_nucleotide = corresponding_nucleotides.keys()[randint(1, 3)]
  temp_dna_sequence = list(dna_sequence)
  temp_dna_sequence[index] = random_nucleotide
  dna_sequence = "".join(str(x) for x in temp_dna_sequence)
  return dna_sequence

# Mutation - In place 
def mutateGeneInPlace(dna_sequence, num_nucleotides_changed):
  i = 0 
  while i <= num_nucleotides_changed:
    index = randint(1, len(dna_sequence) - 1)
    random_nucleotide = corresponding_nucleotides.keys()[randint(1, 3)]
    dna_sequence = dna_sequence[:index - 1] + random_nucleotide + dna_sequence[index:]
    i = i + 1
  return dna_sequence

# returns the MRNA strand of the longest ORF 
def toMRNA(longest_dna_strand): # will always read from 5' to 3', whether it be the complementary or template strand.
  mrna = "".join(longest_dna_strand)
  mrna = mrna.replace('T', 'U')
  return mrna

# returns amino acid sequence synthesized from an mRNA strand
def toAminoAcid(longest_mrna_strand):
  amino_acid_sequence = ""
  mrna_codons = [longest_mrna_strand[i:i+3] for i in range(0, len(longest_mrna_strand), 3)] #replace the open reading frame mrna 
  for codon in mrna_codons:
    amino_acid_sequence += (amino_acids[codon] + ' ')
  return amino_acid_sequence

# returns a list of indices that ATG (start) is found at
def find_start_codons(dna_strand):
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

# Returns the dna_strand minus all the base pairs following the first stop codon found
def clip_after_stop_codon(dna_strand):
  codon_array = textwrap.wrap(dna_strand, 3)
  # print codon_array
  # Finding all of the stop codons in the array 
  stop_codon_list = [codon_array.index("TAA") if "TAA" in codon_array else -1, 
             codon_array.index("TAG") if "TAG" in codon_array else -1,
             codon_array.index("TGA") if "TGA" in codon_array else -1]

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

# Returns longest open reading frame
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

# returns the complementary strand of a strand
def find_complement(input_dna):
  complementary_dna = ""
  for nucleic_acid in input_dna:
    if (nucleic_acid not in list(corresponding_nucleotides.keys())):
      raise Exception(nucleic_acid + " is not a valid nucleotide")
    complementary_dna += corresponding_nucleotides[nucleic_acid]
  # print "complementay dna = " + complementary_dna
  # print "complementary (3' to 5') " + complementary_dna
  complementary_dna = complementary_dna[::-1]
  print "complementary (5' to 3') " + complementary_dna #reverse strand to read from 5' to 3'
  return complementary_dna

# main runs the entire program.
def main():
  rounds = 1
  num_nucleotides_changed = 0
  input_dna_strand = "TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAAATTTCATCCACT"
  results = {}
  while True:
    gene_arr = ["","","","","",""] #array of size 6 to hold all 6 possible genes
    print "\n-------------------------------- ROUND " + str(rounds) + " --------------------------------"
    print "\noriginal strand = " + input_dna_strand
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
    print("The longest ORF converted to amino acids are: " + longest_amino_acid)

    results[rounds] = [input_dna_strand, longest_amino_acid, num_nucleotides_changed]

    answer = raw_input("Continue or terminate program? (c or t): ")
    if (answer == "t"):
      terminateprogram(results)

    answer = raw_input('Would you like to mutate the sequence? (y, n) ')
    if (answer == 'y'):
        mutation_type = raw_input("Type of mutation (i for insertion, d for deletion, p for in place): ")
    
        if mutation_type == 'p':
          num_nucleotides_changed = raw_input('How many nucleotides would you like to see changed? (1 - ' + str(len(input_dna_strand)) + '): ')
          print("original strand: " + input_dna_strand)
          input_dna_strand = mutateGeneInPlace(input_dna_strand, int(num_nucleotides_changed))
          print("mutated  strand: " + input_dna_strand)

        if mutation_type == 'i':
          print("original strand: " + input_dna_strand)
          print("Conducting insertion mutation....")
          input_dna_strand = mutateGeneInsertion(input_dna_strand)
          print("mutated  strand: " + input_dna_strand)
        
        if mutation_type == 'd':
          print("original strand: " + input_dna_strand)
          print("Conducting deletion mutation.....")
          input_dna_strand = mutateGeneDeletion(input_dna_strand)
          print("mutated  strand: " + input_dna_strand)
    
    rounds = rounds + 1
main()
