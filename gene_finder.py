"""
YOUR HEADER COMMENT HERE

@author: Jackie Zeng

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('G')
    'C'
    """
    if(nucleotide == 'A'):
        return 'T'
    elif(nucleotide == 'T'):
        return 'A'
    elif(nucleotide == 'G'):
        return 'C'
    return 'G'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("CTGATGATG")
    'CATCATCAG'
    """
    r_comp = ''
    for i in dna:
        comp = get_complement(i)
        r_comp = comp + r_comp
    return r_comp

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGTAG")
    'ATG'
    >>> rest_of_ORF("ATGGCC")
    'ATGGCC'

    """
    stop1 = 'TAG'
    stop2 = 'TAA'
    stop3 = 'TGA'

    for i in range(3,len(dna), 3):
        codon = dna[i:i+3]
        if(codon == stop1 or codon == stop2 or codon == stop3):
            return dna[0:i]
    return dna

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGTAG")
    ['ATG']
    >>> find_all_ORFs_oneframe("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG']
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG']

    """
    ORF_array = []
    ORF = ''

    start = "ATG"
    i = 0
    while (i < len(dna)-3):
        codon = dna[i:i+3]
        if (codon == start):
            ORF = rest_of_ORF(dna[i:])
            len_ORF = len(ORF)
            ORF_array.append(ORF)
            i = i + len_ORF - 3
        i +=3
    return ORF_array

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG']
    """
    all_ORF = []

    ORF1 = find_all_ORFs_oneframe(dna)
    ORF2 = find_all_ORFs_oneframe(dna[1:])
    ORF3 = find_all_ORFs_oneframe(dna[2:])

    all_ORF = ORF1 + ORF2 + ORF3

    return all_ORF

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    >>> find_all_ORFs_both_strands("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG', 'ATGCAT']
    """
    all_ORF_both_strand = []
    ORF1 = find_all_ORFs(dna)
    ORF2 = find_all_ORFs(get_reverse_complement(dna))

    all_ORF_both_strand = ORF1 + ORF2

    return all_ORF_both_strand

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("ATGCATGAATGTAG")
    'ATGCATGAATGTAG'
    """
    all_ORFs = find_all_ORFs_both_strands(dna)
    longest_ORF = max(all_ORFs, key=len)
    return longest_ORF

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF
    """
    i = num_trials
    ORFs_arr = []
    while (i > 0):
         shuffled_dna = shuffle_string(dna)
         ORF = longest_ORF(shuffled_dna)
         ORFs_arr.append(ORF)
         i = i-1
    overall_longest_ORF = max(ORFs_arr, key=len)
    return len(overall_longest_ORF)


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
        >>> coding_strand_to_AA("ATGTTACAGCGT")
        'MLQ'
    """
    aa_str = ''
    i = 0
    while (i < len(dna)-3):
        codon = str(dna[i: i+3])
        aa = aa_table[codon]
        aa_str += aa
        i += 3
    return aa_str

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    aa_arr = []
    threshold = longest_ORF_noncoding(dna, 1500)
    both_strands_ORFs = find_all_ORFs_both_strands(dna)
    for ORF in both_strands_ORFs:
        if (len(ORF) > threshold):
            aa_str = coding_strand_to_AA(ORF)
            aa_arr.append(aa_str)
    return sorted(aa_arr)


if __name__ == "__main__":
    from load import load_seq
    dna = load_seq("./data/X73525.fa")
    print(gene_finder(dna))

    # import doctest
    # doctest.testmod()
    # doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose=True)
