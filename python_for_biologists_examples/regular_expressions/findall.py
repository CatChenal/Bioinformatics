import re
dna = "CTGCATTATATCGTACGAAATTATACGCGCG!"

result = re.findall(r"[AT]{6,}", dna)

print(result)