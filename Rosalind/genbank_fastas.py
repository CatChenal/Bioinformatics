from bs4 import BeautifulSoup
import requests
import csv


def getPage(url, fmt='lxml'):
    """
    Return the lxml (default) output of a request to url.
    """
	soup = BeautifulSoup(requests.get(url).text, "lxml")
	return soup.string


def save_seqs(seqs, fname):
    with open(fname, 'w') as fw:
        for i in seqs:
            fw.write(i + "\n\n")
            

#array of accession numbers
acc_values = [["JF778997", "andorraMadriu"], ["JF779064", "austriaAchenkirch"], 
              ["JF778998", "belgiumBryssel"], ["JF778877", "bulgariaVitosha"], 
              ["JF779199", "englandDorset"], ["JF779239", "finlandHelsinki"], 
              ["JF779172", "franceAzaysurcher"], ["DQ074380", "germanyBabenhausen"], 
              ["JF779092", "italyArborio"], ["JF778996", "kyrgysstanIssykkul"], 
              ["JF779024", "netherlandsUtrecht"], ["JF779217", "polandWarsaw"], 
              ["JF778987", "romaniaVoslobeni"], ["JF778983", "russiaBorisovka"], 
              ["JF778921", "swedenKrankesjon"], ["JF778997", "andorraMadiru"], 
              ["JF779064", "austriaAchenkirch"], ["JF778993", "usaCambridge"], 
              ["JF778892", "usaMaine"]]

           
#import csv file, read only
f = open("", "r")
reader = csv.reader(f)

#stores every row in the csv as an array
acc_info = [];
for row in reader:
        acc_info.append([row])

#stores the accession numbers and the unique identifier + city of each fire ant sample
acc_values = [];
for data in acc_info:
        acc_values.append([data[0], data[1] + "_" + data[2].split(",")[0]])

#stores fasta sequences from Genbank
fastas = [];
#part of the url. The accession number needs to be added to this
genbank_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="

for acc_num, identifier in acc_values:
	
	#adds the accession number to the url
	url = genbank_url + acc_num + "&rettype=fasta&retmode=text"

	#tries to connect and outputs and exception if not successful
	try:
	    	response = str(getPage(url))
		
		#the first line is the header
		header = response.splitlines()[0]
	
		#header is modified to fasta format
		newFasta = ">" + identifier + " " + header[1:] + "\n" + "".join(response.splitlines()[1:])
		
		#modified fasta sequence is stored
		fastas.append(str(newFasta))
		
	#throw an exception in case it breaks
	except requests.exceptions.Timeout:
	    print("The server didn't respond. Please, try again later.")

    

