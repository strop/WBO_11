from Bio import SeqIO
from Bio import Entrez
def fetch():
	#put your own email here
	Entrez.email="YOUR EMAIL PLEASE!"
	#search
	dict=Entrez.read(Entrez.esearch("protein", "hemoglobin AND mouse[Orgn]"))
	print dict["IdList"]
	# we have a list of IDs matching our
	# so we can fetch them in fasta
	seqs = []
	for x in dict["IdList"]:
	  print Entrez.efetch("protein",id=x, rettype="fasta").read()
	  seqs.append(SeqIO.read(Entrez.efetch("protein",id=x,rettype="fasta"),"fasta"))
	SeqIO.write(seqs, open("hemoglobin.prot", "w"), "fasta") 

#or rather parse them
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo
from Bio.Data import CodonTable

def align(): 
	matrix = MatrixInfo.blosum62
	for x in range(len(seqs)):
	  for y in range(x, len(seqs)):
	    a = pairwise2.align.globaldx(seqs[x].seq, seqs[y].seq, matrix)
	    print len(a)
	    a = a[0]
	    print pairwise2.format_alignment(*a)
