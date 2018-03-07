# gene finder 
program that generates an amino acid sequence from a DNA strand, in addition to exploring the effects of mutations on protein sequence formation.

**warning**
``` Input is provided in the program. Thus our program does not check for DNA sequence validation ```

## Running the program 
run python newgenefinder.py in the console.

### Input 
``` None. Input is provided in the program. ```

### Output 

``` 
-------------------------------- ROUND 1 --------------------------------

original strand = TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAAATTTCATCCACT
~~~~~~ORF 1~~~~~~
['TCA', 'ATG', 'TAA', 'CGC', 'GCT', 'ACC', 'CGG', 'AGC', 'TCT', 'GGG', 'CCC', 'AAA', 'TTT', 'CAT', 'CCA', 'CT']
ATGTAA
~~~~~~ORF 2~~~~~~
['CAA', 'TGT', 'AAC', 'GCG', 'CTA', 'CCC', 'GGA', 'GCT', 'CTG', 'GGC', 'CCA', 'AAT', 'TTC', 'ATC', 'CAC', 'T']
~~~~~~ORF 3~~~~~~
['AAT', 'GTA', 'ACG', 'CGC', 'TAC', 'CCG', 'GAG', 'CTC', 'TGG', 'GCC', 'CAA', 'ATT', 'TCA', 'TCC', 'ACT']

complementary (5' to 3') AGTGGATGAAATTTGGGCCCAGAGCTCCGGGTAGCGCGTTACATTGA
~~~~~~ORF 4~~~~~~
['AGT', 'GGA', 'TGA', 'AAT', 'TTG', 'GGC', 'CCA', 'GAG', 'CTC', 'CGG', 'GTA', 'GCG', 'CGT', 'TAC', 'ATT', 'GA']
~~~~~~ORF 5~~~~~~
['GTG', 'GAT', 'GAA', 'ATT', 'TGG', 'GCC', 'CAG', 'AGC', 'TCC', 'GGG', 'TAG', 'CGC', 'GTT', 'ACA', 'TTG', 'A']
~~~~~~ORF 6~~~~~~
['TGG', 'ATG', 'AAA', 'TTT', 'GGG', 'CCC', 'AGA', 'GCT', 'CCG', 'GGT', 'AGC', 'GCG', 'TTA', 'CAT', 'TGA']
ATGAAATTTGGGCCCAGAGCTCCGGGTAGCGCGTTACATTGA

open reading frames 1-6 (line below): 
['ATGTAA', '', '', '', '', 'ATGAAATTTGGGCCCAGAGCTCCGGGTAGCGCGTTACATTGA']

The longest ORF is: ATGAAATTTGGGCCCAGAGCTCCGGGTAGCGCGTTACATTGA
The longest ORF converted to amino acids are: M K F G P R A P G S A L H - STOP 
Continue or terminate program? (c or t): c
Would you like to mutate the sequence? (y, n) y
Type of mutation (i for insertion, d for deletion, p for in place): i

original strand: TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAAATTTCATCCACT
Conducting insertion mutation....
mutated  strand: TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAAATTGTCATCCACT

-------------------------------- ROUND 2 --------------------------------

original strand = TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAAATTGTCATCCACT
~~~~~~ORF 1~~~~~~
['TCA', 'ATG', 'TAA', 'CGC', 'GCT', 'ACC', 'CGG', 'AGC', 'TCT', 'GGG', 'CCC', 'AAA', 'TTG', 'TCA', 'TCC', 'ACT']
ATGTAA
~~~~~~ORF 2~~~~~~
['CAA', 'TGT', 'AAC', 'GCG', 'CTA', 'CCC', 'GGA', 'GCT', 'CTG', 'GGC', 'CCA', 'AAT', 'TGT', 'CAT', 'CCA', 'CT']
~~~~~~ORF 3~~~~~~
['AAT', 'GTA', 'ACG', 'CGC', 'TAC', 'CCG', 'GAG', 'CTC', 'TGG', 'GCC', 'CAA', 'ATT', 'GTC', 'ATC', 'CAC', 'T']

complementary (5' to 3') AGTGGATGACAATTTGGGCCCAGAGCTCCGGGTAGCGCGTTACATTGA
~~~~~~ORF 4~~~~~~
['AGT', 'GGA', 'TGA', 'CAA', 'TTT', 'GGG', 'CCC', 'AGA', 'GCT', 'CCG', 'GGT', 'AGC', 'GCG', 'TTA', 'CAT', 'TGA']
~~~~~~ORF 5~~~~~~
['GTG', 'GAT', 'GAC', 'AAT', 'TTG', 'GGC', 'CCA', 'GAG', 'CTC', 'CGG', 'GTA', 'GCG', 'CGT', 'TAC', 'ATT', 'GA']
~~~~~~ORF 6~~~~~~
['TGG', 'ATG', 'ACA', 'ATT', 'TGG', 'GCC', 'CAG', 'AGC', 'TCC', 'GGG', 'TAG', 'CGC', 'GTT', 'ACA', 'TTG', 'A']
ATGACAATTTGGGCCCAGAGCTCCGGGTAG

open reading frames 1-6 (line below): 
['ATGTAA', '', '', '', '', 'ATGACAATTTGGGCCCAGAGCTCCGGGTAG']

The longest ORF is: ATGACAATTTGGGCCCAGAGCTCCGGGTAG
The longest ORF converted to amino acids are: M T I W A Q S S G - STOP 
Continue or terminate program? (c or t): c
Would you like to mutate the sequence? (y, n) y
Type of mutation (i for insertion, d for deletion, p for in place): d
original strand: TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAAATTGTCATCCACT
Conducting deletion mutation.....
mutated  strand: TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAATTGTCATCCACT

-------------------------------- ROUND 3 --------------------------------

original strand = TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAATTGTCATCCACT
~~~~~~ORF 1~~~~~~
['TCA', 'ATG', 'TAA', 'CGC', 'GCT', 'ACC', 'CGG', 'AGC', 'TCT', 'GGG', 'CCC', 'AAT', 'TGT', 'CAT', 'CCA', 'CT']
ATGTAA
~~~~~~ORF 2~~~~~~
['CAA', 'TGT', 'AAC', 'GCG', 'CTA', 'CCC', 'GGA', 'GCT', 'CTG', 'GGC', 'CCA', 'ATT', 'GTC', 'ATC', 'CAC', 'T']
~~~~~~ORF 3~~~~~~
['AAT', 'GTA', 'ACG', 'CGC', 'TAC', 'CCG', 'GAG', 'CTC', 'TGG', 'GCC', 'CAA', 'TTG', 'TCA', 'TCC', 'ACT']

complementary (5' to 3') AGTGGATGACAATTGGGCCCAGAGCTCCGGGTAGCGCGTTACATTGA
~~~~~~ORF 4~~~~~~
['AGT', 'GGA', 'TGA', 'CAA', 'TTG', 'GGC', 'CCA', 'GAG', 'CTC', 'CGG', 'GTA', 'GCG', 'CGT', 'TAC', 'ATT', 'GA']
~~~~~~ORF 5~~~~~~
['GTG', 'GAT', 'GAC', 'AAT', 'TGG', 'GCC', 'CAG', 'AGC', 'TCC', 'GGG', 'TAG', 'CGC', 'GTT', 'ACA', 'TTG', 'A']
~~~~~~ORF 6~~~~~~
['TGG', 'ATG', 'ACA', 'ATT', 'GGG', 'CCC', 'AGA', 'GCT', 'CCG', 'GGT', 'AGC', 'GCG', 'TTA', 'CAT', 'TGA']
ATGACAATTGGGCCCAGAGCTCCGGGTAGCGCGTTACATTGA

open reading frames 1-6 (line below): 
['ATGTAA', '', '', '', '', 'ATGACAATTGGGCCCAGAGCTCCGGGTAGCGCGTTACATTGA']

The longest ORF is: ATGACAATTGGGCCCAGAGCTCCGGGTAGCGCGTTACATTGA
The longest ORF converted to amino acids are: M T I G P R A P G S A L H - STOP 
Continue or terminate program? (c or t): c
Would you like to mutate the sequence? (y, n) y
Type of mutation (i for insertion, d for deletion, p for in place): p
How many nucleotides would you like to see changed? (1 - 47): 4
original strand: TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAATTGTCATCCACT
mutated  strand: TCACTTTAACGCGCTACCGGGGGTTCTGGGCCCAATTGTCATCCACT

-------------------------------- ROUND 4 --------------------------------

original strand = TCACTTTAACGCGCTACCGGGGGTTCTGGGCCCAATTGTCATCCACT
~~~~~~ORF 1~~~~~~
['TCA', 'CTT', 'TAA', 'CGC', 'GCT', 'ACC', 'GGG', 'GGT', 'TCT', 'GGG', 'CCC', 'AAT', 'TGT', 'CAT', 'CCA', 'CT']
~~~~~~ORF 2~~~~~~
['CAC', 'TTT', 'AAC', 'GCG', 'CTA', 'CCG', 'GGG', 'GTT', 'CTG', 'GGC', 'CCA', 'ATT', 'GTC', 'ATC', 'CAC', 'T']
~~~~~~ORF 3~~~~~~
['ACT', 'TTA', 'ACG', 'CGC', 'TAC', 'CGG', 'GGG', 'TTC', 'TGG', 'GCC', 'CAA', 'TTG', 'TCA', 'TCC', 'ACT']

complementary (5' to 3') AGTGGATGACAATTGGGCCCAGAACCCCCGGTAGCGCGTTAAAGTGA
~~~~~~ORF 4~~~~~~
['AGT', 'GGA', 'TGA', 'CAA', 'TTG', 'GGC', 'CCA', 'GAA', 'CCC', 'CCG', 'GTA', 'GCG', 'CGT', 'TAA', 'AGT', 'GA']
~~~~~~ORF 5~~~~~~
['GTG', 'GAT', 'GAC', 'AAT', 'TGG', 'GCC', 'CAG', 'AAC', 'CCC', 'CGG', 'TAG', 'CGC', 'GTT', 'AAA', 'GTG', 'A']
~~~~~~ORF 6~~~~~~
['TGG', 'ATG', 'ACA', 'ATT', 'GGG', 'CCC', 'AGA', 'ACC', 'CCC', 'GGT', 'AGC', 'GCG', 'TTA', 'AAG', 'TGA']
ATGACAATTGGGCCCAGAACCCCCGGTAGCGCGTTAAAGTGA

open reading frames 1-6 (line below): 
['', '', '', '', '', 'ATGACAATTGGGCCCAGAACCCCCGGTAGCGCGTTAAAGTGA']

The longest ORF is: ATGACAATTGGGCCCAGAACCCCCGGTAGCGCGTTAAAGTGA
The longest ORF converted to amino acids are: M T I G P R T P G S A L K - STOP 
Continue or terminate program? (c or t): t

         DNA sequence                                    | Number Mutations | Amino Acid Sequence 
Round 1: TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAAATTTCATCCACT |        0         | M K F G P R A P G S A L H - STOP 
Round 2: TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAAATTGTCATCCACT |        0         | M T I W A Q S S G - STOP 
Round 3: TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAATTGTCATCCACT |        0         | M T I G P R A P G S A L H - STOP 
Round 4: TCACTTTAACGCGCTACCGGGGGTTCTGGGCCCAATTGTCATCCACT |        4         | M T I G P R T P G S A L K - STOP 

``` 

## Built Using 
* Python modules
