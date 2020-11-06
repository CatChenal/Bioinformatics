import os
import numpy as np
from collections import (defaultdict,
                         Counter)
from itertools import product
"""
Answers to Rosalind problems: 
http://rosalind.info/problems/list-view/
"""
    
"""DNA"""
def words_freq_from_file(fname=None, char_max=1000):
    with open(fname) as fh:
        line = fh.readline()[:char_max]
            
    cnt = Counter(line)
    print(cnt['A'], cnt['C'], cnt['G'], cnt['T'])
    return
    

"""RNA"""
def rna_from_dna_file(fname=None, char_max=1000):
    with open(fname) as fh:
        t = fh.readline().strip()[:char_max]
    return t.replace('T', 'U')


"""REVC"""
def rev_complement_from_dna_file(fname=None, char_max=1000):
    with open(fname) as fh:
        dna = fh.readline().strip()[:char_max]
    d = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    print(''.join(d[nt] for nt in dna[::-1]))
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
        # no line/break encountered - probably an empty file
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

"""GC"""
def highest_gc_from_fasta_files(multifasta=os.path.join(os.curdir, 'rosalind_gc.txt')):
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

    mismatched = 0
    for a, b in list(zip(s, t)):
        if a != b:
            mismatched +=1
    return mismatched


def test_dH():
    ans = dH('GAGCCTACTAACGGGAT', 'CATCGTAATGACGGCCT')
    print('ans == 7 ? {}'.format(ans == 7))


def two_seqs_from_file(fname=os.path.join(os.curdir, 'rosalind_hamm.txt'),
                       char_max=1000):
    char_max += 1
    with open(fname) as fh:
        seqs = fh.readlines()
    # assume only 2 lines
    s = seqs[0].strip()[:char_max]
    t = seqs[1].strip()[:char_max]
    return s, t

"""
Introduction to Mendelian Inheritance

Figure 1. A Punnett square representing the possible outcomes 
of crossing a heterozygous organism (Yy) with a homozygous 
recessive organism (yy); here, the dominant allele Y corresponds 
to yellow pea pods, and the recessive allele y corresponds to 
green pea pods.
Modern laws of inheritance were first described by Gregor Mendel 
(an Augustinian Friar) in 1865. The contemporary hereditary model, 
called blending inheritance, stated that an organism must exhibit 
a blend of its parent's traits. This rule is obviously violated both 
empirically (consider the huge number of people who are taller than 
both their parents) and statistically (over time, blended traits would 
simply blend into the average, severely limiting variation).

Mendel, working with thousands of pea plants, believed that rather 
than viewing traits as continuous processes, they should instead be 
divided into discrete building blocks called factors. Furthermore, 
he proposed that every factor possesses distinct forms, called alleles.

In what has come to be known as his first law (also known as the law 
of segregation), Mendel stated that every organism possesses a pair 
of alleles for a given factor. If an individual's two alleles for a 
given factor are the same, then it is homozygous for the factor; if 
the alleles differ, then the individual is heterozygous. The first 
law concludes that for any factor, an organism randomly passes one 
of its two alleles to each offspring, so that an individual receives 
one allele from each parent.

Mendel also believed that any factor corresponds to only two possible 
alleles, the dominant and recessive alleles. An organism only needs to 
possess one copy of the dominant allele to display the trait represented 
by the dominant allele. In other words, the only way that an organism can 
display a trait encoded by a recessive allele is if the individual is 
homozygous recessive for that factor.

We may encode the dominant allele of a factor by a capital letter (e.g., A) 
and the recessive allele by a lower case letter (e.g., a). Because a 
heterozygous organism can possess a recessive allele without displaying 
the recessive form of the trait, we henceforth define an organism's 
genotype to be its precise genetic makeup and its phenotype as the 
physical manifestation of its underlying traits.

The different possibilities describing an individual's inheritance of 
two alleles from its parents can be represented by a Punnett square; 
see Figure 1 for an example.

Problem
Figure 2. The probability of any outcome (leaf) in a probability tree 
diagram is given by the product of probabilities from the start of the 
tree to the outcome. For example, the probability that X is blue and Y 
is blue is equal to (2/5)(1/4), or 1/10.
Probability is the mathematical study of randomly occurring phenomena. 
We will model such a phenomenon with a random variable, which is simply 
a variable that can take a number of different distinct outcomes depending 
on the result of an underlying random process.

For example, say that we have a bag containing 3 red balls and 2 blue balls. 
If we let X represent the random variable corresponding to the color of a 
drawn ball, then the probability of each of the two outcomes is given by 
Pr(X=red)=3/5 and Pr(X=blue)=2/5.

Random variables can be combined to yield new random variables. Returning 
to the ball example, let Y model the color of a second ball drawn from the 
bag (without replacing the first ball). The probability of Y being red depends 
on whether the first ball was red or blue. To represent all outcomes of X and 
Y, we therefore use a probability tree diagram. This branching diagram represents 
all possible individual probabilities for X and Y, with outcomes at the endpoints
 ("leaves") of the tree. The probability of any outcome is given by the product 
 of probabilities along the path from the beginning of the tree.

An event is simply a collection of outcomes. Because outcomes are distinct, the
 probability of an event can be written as the sum of the probabilities of its 
 constituent outcomes. For our colored ball example, let A be the event "Y is 
 blue." Pr(A) is equal to the sum of the probabilities of two different outcomes: 
 Pr(X=blue and Y=blue)+Pr(X=red and Y=blue), or 3/10+1/10=2/5 (see Figure 2 above).

Given: Three positive integers k, m, and n, representing a population containing 
k+m+n organisms:  
  k individuals are homozygous dominant for a factor:: AA, 
  m are heterozygous                                :: Aa 
  n are homozygous recessive                        :: aa

Return: The probability that two randomly selected mating organisms will 
produce an individual possessing a dominant allele (and thus displaying 
the dominant phenotype). Assume that any two organisms can mate.

Sample Dataset
2 2 2
Sample Output
0.78333
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
    k # homozygous dominant :: AA
    m # heterozygous        :: Aa
    n # homozyg. recessive  :: aa
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
    print('ans == 0.78333 ? {} ({})'.format(ans == 0.78333, ans))
    
    
def prob_dominant_from_file(fname=os.path.join(os.curdir, 'rosalind_iprb.txt')):
    with open(fname) as fh:
        k, m, n = [int(n) for n in fh.readline().strip().split()]
    print('k, m, n:', k, m, n) 
    return prob_dominant(k, m, n)
#...............................................................
"""
Problem
For a random variable X taking integer values between 1 and n, the expected value of X  
is E(X)=∑nk=1k×Pr(X=k). The expected value offers us a way of taking the long-term  
average of a random variable over a large number of trials.

As a motivating example, let X be the number on a six-sided die. Over a large number  
of rolls, we should expect to obtain an average of 3.5 on the die (even though it's  
not possible to roll a 3.5). The formula for expected value confirms that  
`E(X)=∑6k=1k×Pr(X=k)=3.5.`

More generally, a random variable for which every one of a number of equally spaced  
outcomes has the same probability is called a uniform random variable (in the die  
example, this "equal spacing" is equal to 1). We can generalize our die example to  
find that if X is a uniform random variable with minimum possible value a and  
maximum possible value b, then E(X)=a+b2. You may also wish to verify that for the  
dice example, if Y is the random variable associated with the outcome of a second  
die roll, then E(X+Y)=7.

Given:  
Six nonnegative integers, each of which does not exceed 20,000.  
The integers correspond to the number of couples in a population possessing each  
genotype pairing for a given factor. In order, the six given integers represent the  
number of couples having the following genotypes:
```
AA-AA
AA-Aa
AA-aa
Aa-Aa
Aa-aa
aa-aa
```
Return: The expected number of offspring displaying the dominant phenotype in  
the next generation, under the assumption that every couple has exactly two offspring.  

Sample Dataset  
```
1 0 0 1 0 1
```
Sample Output 
`3.5`
"""

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


def rna_from_file(fname=os.path.join(os.curdir, 'rosalind_prot.txt')):
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
def seqs_from_fasta_files(multifasta=os.path.join(os.curdir, 'rosalind_xx.txt'),
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


def get_profile_printout(prof):
    """
    Return a string for each key (NT) in profile dict; 
    e.g. "A: 5 1 0 0 5 5 0 0" : 1st row
    """
    p = ''
    for k in prof.keys():
        p += '{}: {}\n'.format(k, ' '.join(str(v) for v in prof[k] ))
    return p
 

def save_str_to_file(s, fname):
    with open(fname, 'w') as fh:
        fh.write(s)

    
def get_consensus(seqs):
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
    
    data = seqs_from_fasta_files(multifasta=os.path.join(os.curdir, 'test_profile.fasta'))
    seqs = np.array(data[1])
    consensus, profile = get_consensus(seqs)

    ok1 = consensus == test_consensus
    ok2 = profile == test_prof

    print("test_get_consensus()\n  consensus: {}; profile: {}".format(ok1, ok2))

    
#...........................................................................
def fib_pairs_mortal(n, m, dbg=False):
    """
    Return: total # of pairs after n months have elapsed if all rabbits live for m months.
            Each pair of rabbits reaches maturity in one month and produces a single 
            pair of offspring (one male, one female) each subsequent month.
    n :: generations in months
    m :: lifetime in months
    n ≤ 100 and m ≤ 20.
    # Rabbits born in generation i are produced by rabbits born during 
    # the previous (m-1) generation, ie. births[i] = sum(births[i-m:i-1]). 
    # The sum of births during the last m generations is the current pop since
    # all rabbits older than m are dead.
    """
    if (n <= 0) | (m > 100):
        return "Incorrect number of generations, n: (1-100)"
    if (m <= 0) | (m > 20):
        return "Incorrect lifespan (in months), m: (1-20)"
    
    i = 0
    
    pops = []
    Pop = [1] + [0]*(m - 1)
    pops.append([i, sum(Pop), Pop])
    
    for i in range(1, n):
        Pop = [sum(Pop[1:])] + Pop[:-1]
        pops.append([i, sum(Pop), Pop])

    if dbg:
        print(pops)
    return sum(Pop)


def test_fib_pairs_mortal():
    tests = [(6, 3, 4),
             (94, 16, 19422747110843061063)]
    for t in tests:
        n, m, a = t[0], t[1], t[2]
        out = fib_pairs_mortal(n, m)
        ok =  out == a
        print('test_fib_pairs_mortal({}, {}) = {}?: {}: {}'.format(n, m, a, ok, out))
        
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
    fname = multifasta=os.path.join(os.curdir, 'test_overlap.fasta')
    ans = ['Rosalind_0498 Rosalind_2391',
           'Rosalind_0498 Rosalind_0442',
           'Rosalind_2391 Rosalind_2323']
    
    out = print_overlap_edges(overlap_graph(fname))
    ok =  out == ans
    print('test_overlap_graph() = {}?: {}'.format(ans, ok))


def save_seqs_to_file(seqs, fname):
    with open(fname, 'w') as fw:
        for s in seqs:
            fw.write(s + '\n')
#...........................................................................
"""
Problem:: http://rosalind.info/problems/pmch/
A matching in a graph G is a collection of edges of G for which NO node belongs to more than 1 edge in the collection. See Figure 2 for examples of matchings. 
If G contains an even number of nodes (say 2n), then a matching on G is perfect if it contains n edges, which is clearly the maximum possible. An example of a graph containing a perfect matching is shown in Figure 3.

First, let Kn denote the complete graph on 2n labeled nodes, in which every node is connected to every other node with an edge, and let pn denote the total number of perfect matchings in Kn. For a given node x, there are 2n−1 ways to join x to the other nodes in the graph, after which point we must form a perfect matching on the remaining 2n−2 nodes. This reasoning provides us with the recurrence relation $p_n$=(2n−1)⋅$p_n−1$; using the fact that p1 is 1, this recurrence relation implies the closed equation $p_n$=(2n−1)(2n−3)(2n−5)⋯(3)(1).

Given an RNA string s=s1…sn, a bonding graph for s is formed as follows. 
* First, assign each symbol of s to a node, and arrange these nodes in order around a circle, connecting them with edges called adjacency edges. 
* Second, form all possible edges {A, U} and {C, G}, called basepair edges; we will represent basepair edges with dashed edges, as illustrated by the bonding graph in Figure 4.

Note that a matching contained in the basepair edges will represent one possibility for base pairing interactions in s, as shown in Figure 5. For such a matching to exist, s must have the same number of occurrences of 'A' as 'U' and the same number of occurrences of 'C' as 'G'.

Given: An RNA string s of length at most 80 bp having the same number of occurrences of 'A' as 'U' and the same number of occurrences of 'C' as 'G'.

Return: The total possible number of perfect matchings of basepair edges in the bonding graph of s.

Sample Dataset
>Rosalind_23
AGCUAGUCAU
Sample Output
12
"""

#.........................................................................
"""
Problem
For a random variable X taking integer values between 1 and n, the expected value of X is E(X)=∑nk=1k×Pr(X=k). The expected value offers us a way of taking the long-term average of a random variable over a large number of trials.

As a motivating example, let X be the number on a six-sided die. Over a large number of rolls, we should expect to obtain an average of 3.5 on the die (even though it's not possible to roll a 3.5). The formula for expected value confirms that E(X)=∑6k=1k×Pr(X=k)=3.5.

More generally, a random variable for which every one of a number of equally spaced outcomes has the same probability is called a uniform random variable (in the die example, this "equal spacing" is equal to 1). We can generalize our die example to find that if X is a uniform random variable with minimum possible value a and maximum possible value b, then E(X)=a+b2. You may also wish to verify that for the dice example, if Y is the random variable associated with the outcome of a second die roll, then E(X+Y)=7.

Given: Six nonnegative integers, each of which does not exceed 20,000. The integers correspond to the number of couples in a population possessing each genotype pairing for a given factor. In order, the six given integers represent the number of couples having the following genotypes:

AA-AA
AA-Aa
AA-aa
Aa-Aa
Aa-aa
aa-aa
Return: The expected number of offspring displaying the dominant phenotype in the next generation, under the assumption that every couple has exactly two offspring.

Sample Dataset
1 0 0 1 0 1
Sample Output
3.5
"""

#................................................................
"""
Problem: http://rosalind.info/problems/iprb/

Figure 2 (./images/balls_tree.png). The probability of any outcome (leaf) in a probability tree diagram is given by the product of probabilities from the start of the tree to the outcome. For example, the probability that X is blue and Y is blue is equal to (2/5)(1/4), or 1/10.
Probability is the mathematical study of randomly occurring phenomena. We will model such a phenomenon with a random variable, which is simply a variable that can take a number of different distinct outcomes depending on the result of an underlying random process.

For example, say that we have a bag containing 3 red balls and 2 blue balls. If we let X represent the random variable corresponding to the color of a drawn ball, then the probability of each of the two outcomes is given by Pr(X=red)=35 and Pr(X=blue)=25.

Random variables can be combined to yield new random variables. Returning to the ball example, let Y model the color of a second ball drawn from the bag (without replacing the first ball). The probability of Y being red depends on whether the first ball was red or blue. To represent all outcomes of X and Y, we therefore use a probability tree diagram. This branching diagram represents all possible individual probabilities for X and Y, with outcomes at the endpoints ("leaves") of the tree. The probability of any outcome is given by the product of probabilities along the path from the beginning of the tree; see Figure 2 for an illustrative example.

An event is simply a collection of outcomes. Because outcomes are distinct, the probability of an event can be written as the sum of the probabilities of its constituent outcomes. For our colored ball example, let A be the event "Y is blue." Pr(A) is equal to the sum of the probabilities of two different outcomes: Pr(X=blue and Y=blue)+Pr(X=red and Y=blue), or 310+110=25 (see Figure 2 above).

Given: Three positive integers k, m, and n, representing a population containing k+m+n organisms: k individuals are homozygous dominant for a factor, m are heterozygous, and n are homozygous recessive.

Return: The probability that two randomly selected mating organisms will produce an individual possessing a dominant allele (and thus displaying the dominant phenotype). Assume that any two organisms can mate.

Sample Dataset
2 2 2
Sample Output
0.78333

Hint
Consider simulating inheritance on a number of small test cases in order to check your solution.
"""

                