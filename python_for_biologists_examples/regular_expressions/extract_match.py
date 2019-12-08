import re
dna = "CGATNCGGAACGATC"
m = re.search(r"[^ATGC]", dna)

# m is now a match object
if m:
    print("ambiguous base found!")
    ambig = m.group()
    print("the base is " + ambig)