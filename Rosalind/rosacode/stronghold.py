from pathlib import Path
import numpy as np
from collections import (defaultdict,
                         Counter)
from itertools import product
"""
Answers to Rosalind problems: 
http://rosalind.info/problems/list-view/
"""
DATA = Path(__file__).parent.parent.joinpath('data')

def get_file_chunck(fname=None, char_max=1000):
    """
    Helper function to retrieve a portion of a string file
    from the start up to `char_max` number of characters.

    :param: fname: path of file
    :param: char_max, >=1: maximal number of characters
    """
    if fname is None:
        raise ValueError("get_file_chunck :: Missing file name")
    if char_max < 1:
        raise ValueError("get_file_chunck :: Expected char_max >= 1")

    with open(fname) as fh:
        chunk = fh.readline().strip()[:char_max]

    return chunk


"""DNA"""
def words_freq_from_file(fname=None, char_max=1000):
    line = get_file_chunck(fname=fname, char_max=char_max)
    cnt = Counter(line)
    print(cnt['A'], cnt['C'], cnt['G'], cnt['T'])
    return
    

"""RNA"""
def rna_from_dna_file(fname=None, char_max=1000):
    line = get_file_chunck(fname=fname, char_max=char_max)
    return line.replace('T', 'U')


"""REVC"""
def rev_complement_from_dna_file(fname=None, char_max=1000):
    line = get_file_chunck(fname=fname, char_max=char_max)
    d = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    print(''.join(d[nt] for nt in line[::-1]))
    return


"""FIB"""
def fib_pairs(months, p):
    """
    Return the total number of rabbit pairs that will be present after n months, 
    if we begin with 1 pair and in each generation, every pair of reproduction-age 
    rabbits produces a litter of p rabbit pairs (instead of only 1 pair).
    """
    a = 1
    b = 1
    if (months < 0) | (months > 40) | (p > 5):
        msg = "At least one parameter is invalid:\nmonths, "
        msg += "positive int ≤ 40;  p, # pairs in litter ≤ 5."
        return msg

    if months == 1: 
        return a 
    elif months == 2: 
        return b 
    else:
        for i in range(2, months): #pylint: disable=unused-variable
            a, b = a + b * p, a

        return a


def test_fib_pairs():
    months, pairs, ans = 5, 3, 19
    out = fib_pairs(months, pairs)
    
    ok =  out == ans
    print('test_fib_pairs(5,3) = {}?: {}'.format(ans, ok))


"""GC"""
def SimpleFastaParser(handle):
    """
    Iterate over Fasta records as string tuples.
    For each record a tuple of two strings is returned, the FASTA title
    line (without the leading '>' character), and the sequence (with any
    whitespace removed). The title line is not divided up into an
    identifier (the first word) and comment or description.
    >>> with open("Fasta/dups.fasta") as handle:
    ...     for values in SimpleFastaParser(handle):
    ...         print(values)
    ...
    ('alpha', 'ACGTA')
    ('beta', 'CGTC')
    ('gamma', 'CCGCC')
    ('alpha (again - this is a duplicate entry to test the indexing code)', 'ACGTA')
    ('delta', 'CGCGC')
    """
    # Skip any text before the first record (e.g. blank lines, comments),
    # save title if title line
    for line in handle:
        if line[0] == '>':
            title = line[1:].rstrip()
            break
        elif isinstance(line[0], int):
            # Same exception as for FASTQ files
            raise ValueError("Is this handle in binary mode not text mode?")
    else:
        # no line encountered - probably an empty file
        return

    # Main logic: remove trailing whitespace, and any internal spaces
    # (and any embedded \r which are possible in mangled files
    # when not opened in universal read lines mode)
    lines = []
    for line in handle:
        if line[0] == '>':
            yield title, ''.join(lines).replace(" ", "").replace("\r", "")
            # reset list
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())

    yield title, ''.join(lines).replace(" ", "").replace("\r", "")


"""GC"""
def gc_content(seq):
    cnt = Counter(seq)
    gc = 100 * (cnt['C'] + cnt['G']) / sum(cnt.values())
    return np.round(gc, 6)


"""GC, data/rosalind_gc.txt"""
def highest_gc_from_fasta_files(multifasta=None):
    if multifasta is None:
        raise ValueError("multifasta=None: Use data/rosalind_gc.txt?")
    gc = 0
    with open(multifasta) as fh:
        for seqs in SimpleFastaParser(fh):
            newgc = gc_content(seqs[1])
            if newgc > gc:
                gc = newgc
                gcid = seqs[0]
        
    print(F'{gcid}\n{gc}')
    return


"""HAMM"""
def dH(s, t):
    assert(len(s)==len(t))
    return sum([a!=b for a,b in zip(s, t)])


def test_dH():
    ans = dH('GAGCCTACTAACGGGAT', 'CATCGTAATGACGGCCT')
    print('ans == 7 ? {}'.format(ans == 7))


def two_seqs_from_file(fname=None, char_max=1000):
    if fname is None:
        raise ValueError("fname=None: Use data/rosalind_hamm.txt?")

    char_max += 1
    with open(fname) as fh:
        seqs = fh.readlines()
    # assume only 2 lines
    s = seqs[0].strip()[:char_max]
    t = seqs[1].strip()[:char_max]
    return s, t
#====================================

"""
Introduction to Mendelian Inheritance
"""
def trait_from_punnet_crossings(parent1_traits, parent2_traits, filter_trait):
    # Return a tuple of the count of filter_traits (str) in each of theses crosses:
    # crosses: 1 x 1, 1 x 2, 2 x 2 as per Punnett squares
    # parent[x]_traits: a string of letters representing a trait, e.g: 'Aa'.
    # filter_trait:  string of offspring trait for counting, e.g: 'AA', 'a', 'A'.
    # 
    # CAVEAT: Only meaningful for 2-letter traits (e.g. Dominant/Recessive)
    #         bc crossings :: cartesian product of 2 letters.
    # example:
    # parent1_traits='Aa', parent2_traits='aa' & filter_trait=='A'
    # Punnett squares:
    #   m x m:          m x n:        n x n:
    #    A   a           A   a        a   a
    # A AA  Aa        a  aA  aa     a aa  aa
    # a aA  aa        a  aA  aa     a aa  aa
    #   ------           ------       ------
    #  3 w/A             2 w/A        0 w/A
    #
    if len(parent1_traits) != len(parent2_traits):
        print("The strings denoting the parents' traits must have equal length.")
        return
    
    sep = ''
    if len(filter_trait) < len(parent1_traits):
        sep = ' '
    
    c11 = ['{}{}{}'.format(c[0], sep, c[1]) for c in product(parent1_traits, parent1_traits)]
    c12 = ['{}{}{}'.format(c[0], sep, c[1]) for c in product(parent1_traits, parent2_traits)]
    c22 = ['{}{}{}'.format(c[0], sep, c[1]) for c in product(parent2_traits, parent2_traits)]
    
    if sep == '':
        f11 = [(filter_trait in p)|(filter_trait[::-1] in p) for p in c11].count(True)/4
        f12 = [(filter_trait in p)|(filter_trait[::-1] in p) for p in c12].count(True)/4
        f22 = [(filter_trait in p)|(filter_trait[::-1] in p) for p in c22].count(True)/4
    else:
        f11 = [filter_trait in p for p in c11].count(True)/4
        f12 = [filter_trait in p for p in c12].count(True)/4
        f22 = [filter_trait in p for p in c22].count(True)/4
        
    return (f11, f12, f22)


def prob_dominant(k=2, m=2, n=2):
    """
    One trait two alleles, A and a;
    k # homozygous dominant   :: AA
    m # heterozygous          :: Aa
    n # homozygous recessive  :: aa
    Return: The probability that two randomly selected mating organisms will 
    produce an individual possessing a dominant allele (and thus displaying 
    the dominant phenotype). Assume that any two organisms can mate.
    """
    assert( (k>0) & (m>0) & (n>0))
    
    # Import comb (combination operation, 'N choose K') from scipy
    from scipy.special import comb
    
    pop = k + m + n 
    # Calculate the number of combos that could be made (valid or not):
    tot_combos = comb(pop, 2)
    
    # Calculate the number of combos that have a dominant allele:
    # A_k: dominant alleles if k pop selected:
    A_k = comb(k, 2) + k*m + k*n
    
    # A_mn: dominant alleles from other pops (use Punnett square):
    #   m x m:          m x n:
    #    A   a           A   a
    # A AA  Aa        a  aA  aa
    # a aA  aa        a  aA  aa
    #   ------           ------
    #  3 w/A             2 w/A
    #
    cross_coeffs = trait_from_punnet_crossings('Aa', 'aa', 'A')
    A_mn = cross_coeffs[0] * comb(m, 2) + cross_coeffs[1] * m * n

    return (A_k + A_mn)/tot_combos


def test_prob_dominant():
    k, m, n = 2, 2, 2
    ans = np.round(prob_dominant(k, m, n), 5)
    print(F'Given k:2, m:2, n2\nAns == 0.78333 ? {ans == 0.78333} ({ans})')
    
    
def prob_dominant_from_file(fname=None):
    if fname is None:
        raise ValueError("fname=None: Use data/rosalind_iprb.txt?")

    with open(fname) as fh:
        k, m, n = [int(n) for n in fh.readline().strip().split()]
    print('k, m, n:', k, m, n) 
    return prob_dominant(k, m, n)
#...............................................................

# AUG = 'M' = 'start'
rna_codons = {'UUU':'F', 'CUU':'L', 'AUU':'I', 'GUU':'V',
              'UUC':'F', 'CUC':'L', 'AUC':'I', 'GUC':'V',
              'UUA':'L', 'CUA':'L', 'AUA':'I', 'GUA':'V',
              'UUG':'L', 'CUG':'L', 'AUG':'M', 'GUG':'V',
              'UCU':'S', 'CCU':'P', 'ACU':'T', 'GCU':'A',
              'UCC':'S', 'CCC':'P', 'ACC':'T', 'GCC':'A',
              'UCA':'S', 'CCA':'P', 'ACA':'T', 'GCA':'A',
              'UCG':'S', 'CCG':'P', 'ACG':'T', 'GCG':'A',
              'UAU':'Y', 'CAU':'H', 'AAU':'N', 'GAU':'D',
              'UAC':'Y', 'CAC':'H', 'AAC':'N', 'GAC':'D',
              'UAA':'stop', 'CAA':'Q', 'AAA':'K', 'GAA':'E',
              'UAG':'stop', 'CAG':'Q', 'AAG':'K', 'GAG':'E',
              'UGU':'C', 'CGU':'R', 'AGU':'S', 'GGU':'G',
              'UGC':'C', 'CGC':'R', 'AGC':'S', 'GGC':'G',
              'UGA':'stop', 'CGA':'R', 'AGA':'R', 'GGA':'G',
              'UGG':'W', 'CGG':'R', 'AGG':'R', 'GGG':'G'}


def rna_from_file(fname=DATA.joinpath('rosalind_prot.txt')):
    if fname is None:
        raise ValueError("fname=None: Use data/rosalind_prot.txt?")

    MAX_LEN = 10_001

    with open(fname) as fh:
        seq = fh.readline().strip()[:MAX_LEN]

    return seq

    
def protein_seq_from_rna(rna):
    p = ''
    # assume no intervening seq
    for j in range(0, len(rna), 3):
        i = j
        aa = rna_codons[rna[i:j+3]]
        if aa == 'stop':
            break
        p += aa
    return p


def test_protein_seq_from_rna():
    rna = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
    p = protein_seq_from_rna(rna)
    test = 'MAMAPRTEINSTRING'
    ans = p == test
    print("protein_seq_from_rna() == '{}' ? {}".format(test, ans))
    
    return ans

#........................................................
def find_all(text, sub, ioffset=1, with_overlap=True):
    start = 0
    while True:
        start = text.find(sub, start)
        if start == -1: return
        
        yield start + ioffset
        
        if with_overlap:
            start += 1
        else:
            start += len(sub)


def get_motif_locs(dna, motif):
    
    if not len(motif) | (len(motif) > len(dna)):
        return 'Motif can neither be empty nor longer than the search string.'
    
    locs = [loc for loc in find_all(dna, motif)]
    return ' '.join(str(i) for i in find_all(dna, motif))


def test_get_motif_locs():
    dna = 'GATATATGCATATACTT'
    motif = 'ATAT'
    ok = '2 4 10'
    ans = get_motif_locs(dna, motif)
    
    check = ans == ok
    print("get_motif_locs() == '{}' ? {}".format(ans, check))
    
    return check

#...............................................................
def seqs_from_fasta_files(multifasta=DATA.joinpath('rosalind_xx.txt'),
                          split_seq=True,
                          max_len = 1000):
    max_len += 1
    seqs = []
    ids = []
    with open(multifasta) as fh:
        # SimpleFastaParser:: (id, seq)
        for item in SimpleFastaParser(fh):
            ids.append(item[0])
            if not split_seq:
                seqs.append(item[1][:max_len])
            else:
                seqs.append([b for b in item[1][:max_len]])

    return ids, seqs


def save_str_to_file(s, fname):
    with open(fname, 'w') as fh:
        fh.write(s)


def get_profile_printout(prof):
    """
    Return a string for each key (NT) in profile dict; 
    e.g. "A: 5 1 0 0 5 5 0 0" : 1st row
    """
    p = ''
    for k in prof.keys():
        p += '{}: {}\n'.format(k, ' '.join(str(v) for v in prof[k] ))
    return p
 

def get_consensus0(seqs):
    """
    seqs: an m sequences x n nucleotides array.
    """
    profile = defaultdict(list)
    consensus = ''
    nts = 'ACGT'
    
    for c in range(seqs.shape[1]):
        col = seqs[:, c].flatten()
        cnt = Counter(col)
        for bp in nts:
            profile[bp].append(cnt[bp])

    max_idx = np.argmax(np.array([v for v in profile.values()]), axis=0)
    consensus = ''.join(nts[i] for i in max_idx)
    
    return consensus, profile


def get_consensus(seqs):
    """
    seqs: an m sequences x n nucleotides array.
    """
    nts = 'ACGT'
    profile = defaultdict(list)
    seqs = np.array(seqs)
    for c in range(seqs.shape[1]):
        col = seqs[:, c].flatten()
        cnt = Counter(col)
        for base in nts:
            profile[base].append(cnt[base])

    max_idx = np.argmax(np.array([v for v in profile.values()]), axis=0)
    consensus = ''.join(nts[i] for i in max_idx)
    
    return consensus, profile


def test_get_consensus():
    """
    seqs: an m sequences x n nucleotides array.
    """
    test_consensus = "ATGCAACT"
    test_prof = {'A': [5, 1, 0, 0, 5, 5, 0, 0],
                 'C': [0, 0, 1, 4, 2, 0, 6, 1],
                 'G': [1, 1, 6, 3, 0, 1, 0, 0],
                 'T': [1, 5, 0, 0, 0, 1, 1, 6]}
    test_profile = get_profile_printout(test_prof)
    
    data = seqs_from_fasta_files(multifasta=DATA.joinpath('test_profile.fasta'))
    seqs = np.array(data[1])
    consensus, profile = get_consensus(seqs)

    ok1 = consensus == test_consensus
    ok2 = profile == test_prof

    print(F"test_get_consensus()\n  consensus: {ok1}; profile: {ok2}")

    
#...........................................................................
def fib_pairs_mortal(n, m, dbg=False):
    """
    Return: total # of pairs after n months have elapsed if all rabbits live
            for m months.
            Each pair of rabbits reaches maturity in one month and produces 
            a single pair of offspring (one male, one female) each subsequent
            month.
    n :: generations in months, n ≤ 100.
    m :: lifetime in months, m ≤ 20.
    dbg :: debug flag, to print intermediate result
    """
    if (n <= 0) | (n > 100):
        return "Incorrect number of generations, n: (1-100)"
    if (m <= 0) | (m > 20):
        return "Incorrect lifespan (in months), m: (1-20)"
    
    # Rabbits born in generation i are produced by rabbits born during 
    # the previous (m-1) generation, ie. births[i] = sum(births[i-m:i-1]). 
    # The sum of births during the last m generations is the current pop since
    # all rabbits older than m are dead.
    
    Pop = [1] + [0]*(m - 1)
    if dbg:
        pops = []
        pops.append([0, sum(Pop), Pop])
    
    for i in range(1, n):
        Pop = [sum(Pop[1:])] + Pop[:-1]
        if dbg:
            pops.append([i, sum(Pop), Pop])

    if dbg:
        print(pops)
    return sum(Pop)


def test_fib_pairs_mortal():
    tests = [(6, 3, 4),
             (94, 16, 19422747110843061063)]
    #for t in tests:
    #    n, m, a = t[0], t[1], t[2]
    for n, m, a in tests:
        out = fib_pairs_mortal(n, m)
        ok = out == a
        print(F'test_fib_pairs_mortal({m}, {n})\nAns = {a}?: {ok}: {out}')
        
#...........................................................................
def overlapped(seq1, seq2, k=3):
    """
    Seq1 and seq2 are different strings.
    """
    return seq2.startswith(seq1[-k:])


def print_overlap_edges(o_graph_gen):
    out = []
    for e in list(o_graph_gen):
        out.append(' '.join([e[0], e[1]]))
    return out
    
    
def overlap_graph(fast_file, k=3):
    """
    For a collection of strings and a positive integer 𝑘, the overlap graph
    is a directed graph 𝑂𝑘 in which each string is represented by a node, 
    and string 𝑠 is connected to string 𝑡 with a directed edge when there is 
    a length 𝑘 suffix of 𝑠 that matches a length 𝑘 prefix of 𝑡, as long as 𝑠≠𝑡.
    Return a list of tuples (s, t) = edges.
    """
    data = seqs_from_fasta_files(fast_file, split_seq=False, max_len=10_000)
    ids = data[0]
    seqs = np.array(data[1])
    n = seqs.shape[0]
    
    O = []
    for i in range(n):
        for j in range(n):
            ok = (i != j) & (seqs[i] != seqs[j]) 
            if ok & overlapped(seqs[i], seqs[j]):
                yield (ids[i], ids[j])


def test_overlap_graph():
    fname = multifasta=DATA.joinpath('test_overlap.fasta')
    ans = ['Rosalind_0498 Rosalind_2391',
           'Rosalind_0498 Rosalind_0442',
           'Rosalind_2391 Rosalind_2323']
    
    out = print_overlap_edges(overlap_graph(fname))
    ok = out == ans
    print('test_overlap_graph()\nAns = {}?: {}'.format(ans, ok))


def save_seqs_to_file(seqs, fname):
    with open(fname, 'w') as fw:
        for s in seqs:
            fw.write(s + '\n')
#...........................................................................
