import re

dna = "CGATGCTGAATTCGACTGC"
if re.search(r"GC[ATGC]GC", dna):
    print("found restriction site!")

dna = "ATCGCGAATTCAC"
if re.search(r"GG(A|T)CC", dna):
	print("restriction site found!")

dna = "ATCGCGAATTCAC"
if re.search(r"GAATTC", dna):
	print("restriction site found!")


