#database of codon table

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