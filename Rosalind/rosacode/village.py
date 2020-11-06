"""
Problem
Given: Two positive integers a and b, each less than 1000.

Return: The integer corresponding to the square of the hypotenuse of the right triangle whose legs have lengths a and b.
"""
def sqr_hypo(a, b):
    a = int(a)
    b = int(b)
    assert ((a < 1000) & (b < 1000))

    return a**2 + b**2


def sqr_hypo_from_file():
    fname = os.path.join(os.curdir, 'rosalind_ini2.txt')
    with open(fname) as fh:
        data = [L.strip().split() for L in fh.readlines()]

    return [sqr_hypo(d[0], d[1]) for i, d in enumerate(data)]


"""
Problem
Given: A string s of length at most 200 letters and four integers a, b, c and d.

Return: The slice of this string from indices 
a through b and c through d (with space in between), inclusively.
In other words, we should include elements s[b] and s[d] in our slice.

Sample Dataset
HumptyDumptysatonawallHumptyDumptyhadagreatfallAlltheKingshorsesandalltheKingsmenCouldntputHumptyDumptyinhisplaceagain.
22 27 97 102
Sample Output
Humpty Dumpty
"""
import os


def str_slice(s, a, b, c, d):
    assert(len(s) < 201)
    a = int(a)
    b = int(b)
    c = int(c)
    d = int(d)

    return ' '.join(w for w in [s[a:b+1], s[c:d+1]])


def data_from_file3(fname):
    fname = os.path.join(os.curdir, 'rosalind_ini3.txt')
    with open(fname) as fh:
        data = fh.readlines()

    return data[0].strip(), data[1].split()

"""
Problem
Given: Two positive integers a and b (a<b<10000).
Return: The sum of all odd integers from a through b, inclusively.

Sample Dataset
100 200
Sample Output
7500
"""
def sum_odds(a, b):
    a = int(a)
    b = int(b)
    assert(all([a>0, b>0, b<10000, a<b]))

    return sum(i for i in range(a, b) if i % 2 == 1)

"""
Problem
Given: A file containing at most 1000 lines.
Return: A file containing all the even-numbered lines from the original file. Assume 1-based numbering of lines.
"""
def even_lines_from_file(fname=os.path.join(os.curdir, 'rosalind_ini5.txt')):
    i = 1
    line_max = 1001
    data = ''
    with open(fname) as fh:
        for L in fh.readlines():
            if (i < line_max) & (i % 2 == 0):
                data = data + '{}\n'.format(L)
            i += 1

    return data

"""
Problem
Given: A string s of length at most 10000 letters.

Return: The number of occurrences of each word in s, where words are separated by spaces. 
Words are case-sensitive, and the lines in the output can be in any order.

Sample Dataset
We tried list and we tried dicts also we tried Zen
Sample Output
and 1
We 1
tried 3
dicts 1
list 1
we 2
also 1
Zen 1
"""

def words_freq_from_file(fname=os.path.join(os.curdir, 'rosalind_ini6.txt')):
    from collections import defaultdict

    char_max = 10001
    d = defaultdict(int)
    with open(fname) as fh:
        for k in fh.readline()[:char_max].split():
            d[k] += 1
    
    for k, v in d.items():
        print(k,v)