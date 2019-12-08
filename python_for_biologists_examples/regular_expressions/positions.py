import re

dna = "CGATNCGGAACGATC"
m = re.search(r"[^ATGC]", dna)

if m:
    print("ambiguous base found!")
    print("at position " + str(m.start()))