import re

dna = open("dna.txt").read().rstrip("\n") 
all_cuts = [0] 
for match in re.finditer(r"A[ATGC]TAAT", dna): 
    all_cuts.append(match.start() + 3) 
all_cuts.append(len(dna)) 
print(all_cuts) 

for i in range(1,len(all_cuts)): 
    this_cut_position = all_cuts[i] 
    previous_cut_position = all_cuts[i-1] 
    fragment_size = this_cut_position - previous_cut_position 
    print("one fragment size is "  + str(fragment_size)) 
