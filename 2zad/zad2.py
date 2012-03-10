# -*- coding: utf-8 -*-
from Bio import Entrez
def fetch():
	#put your own email here
	Entrez.email="YOUR EMAIL PLEASE!"
	#search
	dict=Entrez.read(Entrez.esearch("nucleotide", "hemoglobin AND mouse[Orgn]"))
	print dict["IdList"]
	# we have a list of IDs matching our
	# so we can fetch them in fasta
	seqs = []
	for x in dict["IdList"]:
	  print Entrez.efetch("nucleotide",id=x, rettype="fasta").read()
	  seqs.append(SeqIO.read(Entrez.efetch("nucleotide",id=x,rettype="fasta"),"fasta"))
	SeqIO.write(seqs, open("hemoglobin", "w"), "fasta") 

#or rather parse them
from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo
from Bio.Data import CodonTable

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

def_subst_mat = MatrixInfo.blosum62
def_codon_table = CodonTable.standard_dna_table.forward_table

class Direction:
  Left, Top, Diag, Start = range(4)

def codon_at(seq, pos):
  return seq[pos:pos+3]

def translate(codon, table):
  return table[codon]

def safe_get_cost(subst_mat, amino1, amino2):
  try:
    return subst_mat[amino1, amino2]
  except Exception:
    return 0

def calculateDiagonalScore(s1, s2, x, y, subst_mat, alpha, codon_table):
  c1 = codon_at(s1, x-3)
  c2 = codon_at(s2, y-3)
  if c1 == c2: # brak mutacji kodonu
    amino = translate(c1, codon_table)
    return safe_get_cost(subst_mat, amino, amino)
  else: # mutacja
    a1 = translate(c1, codon_table)
    a2 = translate(c2, codon_table)
    if a1 != a2: # mutacja zmieniająca aminokwas
      return safe_get_cost(subst_mat,a1,a2)
    else: # mutacja nie zmieniająca aminokwasu
      return alpha

def wasStop(scores, x, y):
  return scores[x][y][2]

def isStop(codon, codon_table): # zgodnie z implementacją w Biopythonie: kodonu stopu nie ma w tablicy translacji
  try:
    codon_table[codon]
    return False
  except Exception:
    return True

def print_scores(scores):
  for x in [ [x[0] for x in row ] for row in scores ]:
    print x

def backtrack(scores, bestScore, s1, s2):
  x = bestScore[0]
  y = bestScore[1]
  recovered1 = ""
  recovered2 = ""
  while x+y > 0 and (not scores[x][y][1] == Direction.Start):
    if scores[x][y][1] == Direction.Top:
      x = x-3
      recovered1 = "---" + recovered1
      recovered2 = codon_at(s2, x) + recovered2
    elif scores[x][y][1] == Direction.Left:
      y = y-3
      recovered2 = "---" + recovered2
      recovered1 = codon_at(s1, y) + recovered1
    elif scores[x][y][1] == Direction.Diag:
      y = y-3
      x = x-3
      recovered2 = codon_at(s2, x) + recovered2
      recovered1 = codon_at(s1, y) + recovered1
    else:
      raise Exception("Nieprawidlowy kierunek")
  return (recovered1, recovered2)

def align_coding(seq1, seq2, codon_table = def_codon_table, subst_mat = def_subst_mat, d = -0.5, alpha = -0.5):
  (s1, s2) = (seq1, seq2)
  bestScore = (0,0)
  scores = [ [(0, Direction.Start, False) for x in range(len(s1))] for x in range(len(s2)) ] # budujemy pustą macierz uliniowienia
  for x in range(3, len(s2)):
    for y in range(3, len(s1)):
      myResult = (0, Direction.Start, False)
      # z lewej:
      if wasStop(scores, x, y-3):
        horizontalResult = (d, Direction.Start, isStop(codon_at(s1, y-3), codon_table))
      else:
        horizontalResult = (d+scores[x][y-3][0], Direction.Left, isStop(codon_at(s1, y-3), codon_table))
      if horizontalResult[0] > myResult[0]:
        myResult = horizontalResult
      # z góry:
      if wasStop(scores, x-3, y):
        verticalResult = (d, Direction.Start, isStop(codon_at(s2, x-3), codon_table))
      else:
        verticalResult = (d+scores[x-3][y][0], Direction.Top, isStop(codon_at(s2, x-3), codon_table))
      if verticalResult[0] > myResult[0]:
        myResult = verticalResult
      # na skos:
      if wasStop(scores, x-3, y-3):
	if isStop(codon_at(s2, x-3), codon_table) or isStop(codon_at(s1, y-3), codon_table):
	  diagonalResult = (0, Direction.Start, False)
	else:
          diagonalResult = (calculateDiagonalScore(s2, s1, x, y, subst_mat, alpha, codon_table), Direction.Start, isStop(codon_at(s2, x-3), codon_table) or isStop(codon_at(s1, y-3), codon_table))
      else:
	if isStop(codon_at(s2, x-3), codon_table) or isStop(codon_at(s1, y-3), codon_table):
	  diagonalResult = (0, Direction.Start, False)
	else:
          diagonalResult = (calculateDiagonalScore(s2, s1, x, y, subst_mat, alpha, codon_table)+scores[x-3][y-3][0], Direction.Diag, isStop(codon_at(s2, x-3), codon_table) or isStop(codon_at(s1, y-3), codon_table))
      if diagonalResult[0] > myResult[0]:
        myResult = diagonalResult
      if myResult[0] > scores[bestScore[0]][bestScore[1]][0]:
        bestScore = (x,y)
      scores[x][y] = myResult
  print_scores(scores)
  a,b =  backtrack(scores, bestScore, s1, s2)
  print a
  print b
