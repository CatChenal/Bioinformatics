import re

accs = ["xkn59438", "yhdck2", "eihd39d9", "chdsye847", "hedle3455", "xjhd53e", "45da", "de37dp"]

print("contains 5:")
for acc in accs: 
    if re.search(r"5", acc): 
        print("\t" + acc) 

print("contains d or e:")
for acc in accs: 
    if re.search(r"(d|e)", acc): 
        print("\t" + acc) 

print("contains d and e in that order:")
for acc in accs: 
    if re.search(r"d.*e", acc): 
        print("\t" + acc)

print("contains d and e separated by a single letter:")
for acc in accs: 
    if re.search(r"(d.e)", acc): 
        print("\t" + acc) 

print("contains d and e in any order:")
for acc in accs: 
    if re.search(r"d.*e", acc) or re.search(r"e.*d", acc): 
        print("\t" + acc) 

print("starts with either x or y:")
for acc in accs: 
    if re.search(r"^(x|y)", acc): 
        print("\t" + acc) 

print("starts with either x or y and ends with e:")
for acc in accs: 
    if re.search(r"^(x|y).*e$", acc): 
        print("\t" + acc) 

print("contains three or more numbers in a row:")
for acc in accs: 
    if re.search(r"\d{3,}", acc): 
        print("\t" + acc) 

print("ends with d followed by either a, r or p:")
for acc in accs: 
    if re.search(r"d[arp]$", acc): 
        print("\t" + acc) 
