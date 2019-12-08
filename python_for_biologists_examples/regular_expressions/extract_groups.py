import re

scientific_name = "Homo sapiens"

m = re.search("(.+) (.+)", scientific_name)

if m:
    genus = m.group(1)
    species = m.group(2)
    print("genus is " + genus + ", species is " + species)