# gene finder 
program that generates an amino acid sequence from a DNA strand, in addition to exploring the effects of mutations on protein sequence formation.

**warning**
new

## Running the program 
run python genefinder.py in the console.

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
mutated  strand: TCAATGTAACGCGCTACCCGGAGCTCTGGGCCCAAATTTCATCCACT

-------------------------------- ROUND 2 --------------------------------

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
Continue or terminate program? (c or t): 

``` 

## Built Using 
* Python modules