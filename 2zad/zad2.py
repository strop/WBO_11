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
	SeqIO.write(seqs, open("hemoglobin", "w"), "fasta") 

#or rather parse them
from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo

parser = SeqIO.parse("hemoglobin", "fasta")
seqs = [seq for seq in parser]
print len(seqs)
#seqs = seqs[1:3]
def align(): 
	matrix = MatrixInfo.blosum62
	for x in range(len(seqs)):
	  for y in range(x, len(seqs)):
	    a = pairwise2.align.globaldx(seqs[x].seq, seqs[y].seq, matrix)
	    print len(a)
	    a = a[0]
	    print pairwise2.format_alignment(*a)

# każde uliniowienie ma długość będąca wielokrotnościa 3.
# Każda przerwa (spacja) w uliniowieniu ma długość będącą wielokrotnością 3.
# Uliniowienia zawierające kodon stopu na pozycji innej niż ostatnia są niedozwolone
#
# subst_mat[x][y] -> kara za zmianę, która powoduje zamianę aminokwasu z x na y
# subst_mat[x][x] -> nagroda za niezmieniony kodon
# d -> kara za przerwę długości 3
# alpha -> kara za zmianę nukleotydu, która nie zmienia kodowanego aminokwasu
# zwraca jedno z optymalnych uliniowień i jego koszt
# Uwaga: prefix i suffix uliniawianych sekwencji nie muszą być wielokrotnościami 3
def align_coding(seq1, seq2, codon_table, subst_mat, d, alpha):
  print 2
