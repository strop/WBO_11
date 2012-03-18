# -*- coding: utf-8 -*-
import random
import math
from Bio import Seq
from Bio import SeqIO
from Bio.Seq import Alphabet
from Bio.SubsMat import MatrixInfo
from Bio.Data import CodonTable

def_subst_mat = MatrixInfo.blosum62
def_codon_table = CodonTable.standard_dna_table.forward_table
def_protein_alphabet = Alphabet.IUPAC.IUPACProtein()
def_dna_alphabet = Alphabet.IUPAC.IUPACUnambiguousDNA()

""" Funkcja przeprowadza doświadczenie losowe z prawdopodobieństwem sukcesu mu """
def hit(rand, mu):
  return rand.random() <= mu;

""" Funkcja mutuje sekwencję seq na pozycji pos i zwraca wylosowany znak
    Znak po mutacji jest wybierany na podstawie alfabetu sekwencji, jeśli sekwencja go nie ma, to używamy domyślnego alfabetu białkowego """
def mutate(pos, seq):
  print "Mutating " + str(seq[pos]) + " into ", 
  try:
    other = [x for x in seq.alphabet.letters]
  except Exception:                                   # nasza sekwencja nie ma alfabetu, albo jest pusty
    other = [x for x in def_protein_alphabet.letters] # więc markujemy białko
  other.remove(seq[pos]) # tutaj są znaki, na które można zamienić mutowaną pozycję
  result = random.Random().choice(other)
  print str(result)
  return result

""" Funkcja generuje odwrotną tablicę translacji (to znaczy słownik, w którym każdemu aminokwasowi jest
    przypisana lista kodonów, które go kodują) """
def get_reverse_codon_table(forward = def_codon_table):
  reverse = {}
  for k in forward.keys():
    if not forward[k] in reverse:
      reverse[forward[k]] = []
    reverse[forward[k]].append(k)
  print reverse
  return reverse

""" Funkcja dokonuje odwrotnej translacji przypisując każdemu aminokwasowi losową kodującą go trójkę """
def reverse_translate(seq):
  reverse = get_reverse_codon_table()
  pchars = [x for x in seq] # tablica znaków w sekwencji
  # obcinamy po pierwszym stopie
  try:
    pchars = pchars[0:pchars.index('*')]
  except Exception: # nie było kodonu stopu
    pass
  r = random.Random()
  rand_mRNA = ''.join([r.choice(reverse[x]) for x in pchars ])
  return Seq.Seq(rand_mRNA, def_dna_alphabet)

""" Funkcja generuje losowy ciąg znaków podanego alfabetu o długości l """
def gen_random_dna(alphabet, l = 100):
  try:
    chars = [x for x in alphabet.letters]
  except Exception:
    chars = [x for x in def_dna_alphabet.letters]
  return Seq.Seq(''.join(random.Random().sample(chars*int(math.ceil(float(l)/float(len(chars)))), l)), alphabet)

""" seq - sekwencja białkowa
    n - liczba przebiegów
    mu - prawdopodobieństwo zajścia mutacji, dla każdej pozycji, w każdym przebiegu
    zwraca - sekwencję DNA, która jest odwrotnie przetłumaczoną zmutowaną sekwencją seq, z dopisanymi 100 losowymi nukleotydami na początku i końcu """
def mutated_DNA(seq, n = 1, mu = 0.5):
  rand = random.Random()
  s = seq.tomutable()
  for i in range(n):        # n razy
    for x in range(len(s)): # sprawdzamy na każdej pozycji
      if s[x] != '*':
        if hit(rand, mu):     # czy zaszła mutacja
          s[x] = mutate(x, s)
  s = reverse_translate(s)
  s.alphabet = def_dna_alphabet
  prefix = gen_random_dna(s.alphabet)
  suffix = gen_random_dna(s.alphabet)
  return prefix + s + suffix

def make_test_samples(seq, n = 1, mu = 0):
  s1 = mutated_DNA(seq, n, mu)
  s2 = mutated_DNA(seq, n, mu)
  return (s1,s2)

""" Każda przerwa (spacja) w uliniowieniu ma długość będącą wielokrotnością 3.
    Uliniowienia zawierające kodon stopu na pozycji innej niż ostatnia są niedozwolone
   
    subst_mat[x][y] -> kara za zmianę, która powoduje zamianę aminokwasu z x na y
    subst_mat[x][x] -> nagroda za niezmieniony kodon
    d -> kara za przerwę długości 3
    alpha -> kara za zmianę nukleotydu, która nie zmienia kodowanego aminokwasu
    zwraca jedno z optymalnych uliniowień i jego koszt
    Uwaga: prefix i suffix uliniawianych sekwencji nie muszą być wielokrotnościami 3 """

""" Enum kierunków poruszania się w macierzy """
class Direction:
  Left, Top, Diag, Start = range(4)

def codon_at(seq, pos):
#  print seq[pos:pos+3].tostring()
  return seq[pos:pos+3].tostring()

def translate(codon, table):
  return table[codon]

def safe_get_cost(subst_mat, amino1, amino2):
  try:
    return subst_mat[amino1, amino2]
  except Exception:
    return 0

""" Funkcja oblicza punktację, jaką ma komórka macierzy do której doszliśmy po skosie """
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

""" Sprawdza, czy dany kodon jest kodonem stopu według podanej tablicy """
def is_stop(codon, codon_table): # zgodnie z implementacją w Biopythonie: kodonu stopu nie ma w tablicy translacji
  try:
    codon_table[codon]
    return False
  except Exception:
    return True

""" Funkcja wypisuje wynikową macierz punktacji """
def print_scores(scores):
  for x in [ [x[0] for x in row ] for row in scores ]:
    print x

def backtrack(scores, bestScore, s1, s2):
  x = bestScore[0]
  y = bestScore[1]
  recovered1 = ""
  recovered2 = ""
  while x > 0 and y > 0 and (not scores[x][y][1] == Direction.Start) and scores[x][y][0] != 0:
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

""" Zgodnie z założeniem, że kodon stopu nie może wystepować wewnątrz uliniowienia, sprawdzamy, czy nasz
    aktualny ruch nie spowoduje dodania kodonu stopu do uliniowienia 
    Jest sporo warunków, ale sprawdzamy każdą możliwość zakradnięcia się kodonu stopu"""
def move_causes_stop_inclusion(s1, s2, y, x, direction, codon_table, scores):
  causes = False
  if direction == Direction.Diag:
    causes = causes or (is_stop(codon_at(s1, y-3), codon_table) or is_stop(codon_at(s2, x-3), codon_table))
    if scores[x-3][y-3][1] != Direction.Top:
      causes = causes or (is_stop(codon_at(s1, y-4), codon_table) or is_stop(codon_at(s1, y-5), codon_table))
    if scores[x-3][y-3][1] != Direction.Left:
      causes = causes or (is_stop(codon_at(s2, x-4), codon_table) or is_stop(codon_at(s2, x-5), codon_table))
  elif direction == Direction.Left:
    causes = causes or is_stop(codon_at(s1, y-3), codon_table)
    if scores[x][y-3][1] != Direction.Top:
      causes = causes or is_stop(codon_at(s1, y-4), codon_table) or is_stop(codon_at(s1, y-5), codon_table)
  elif direction == Direction.Top:
    causes = causes or is_stop(codon_at(s2, x-3), codon_table)
    if scores[x-3][y][1] != Direction.Left:
      causes = causes or is_stop(codon_at(s2, x-4), codon_table) or is_stop(codon_at(s2, x-5), codon_table)
  return causes

def align_coding(seq1, seq2, codon_table = def_codon_table, subst_mat = def_subst_mat, d = -1, alpha = 1):
  (s1, s2) = (seq1, seq2)
  bestScore = (0,0)
  scores = [ [(0, Direction.Start) for x in range(len(s1)+2)] for x in range(len(s2)+2) ] # budujemy pustą macierz uliniowienia
  for x in range(3, len(s2)):
    for y in range(3, len(s1)):
      # z lewej:
      horizontalScore = d
      if not move_causes_stop_inclusion(s1, s2, y, x, Direction.Left, codon_table, scores):
        horizontalScore = horizontalScore + scores[x][y-3][0]
      verticalScore = d
      if not move_causes_stop_inclusion(s1, s2, y, x, Direction.Top, codon_table, scores):
        verticalScore = verticalScore + scores[x-3][y][0]
      if not move_causes_stop_inclusion(s1, s2, y, x, Direction.Diag, codon_table, scores):
        diagonalScore = calculateDiagonalScore(s2, s1, x, y, subst_mat, alpha, codon_table) + scores[x-3][y-3][0]
      else:
        diagonalScore = 0
      finalScore = max(horizontalScore, verticalScore, diagonalScore, 0)
      direction = Direction.Diag
      if finalScore == horizontalScore:
        direction = Direction.Left
      elif finalScore == verticalScore:
        direction = Direction.Top
      scores[x][y] = (finalScore, direction)
      if finalScore > scores[bestScore[0]][bestScore[1]][0]:
        bestScore = (x,y)
  #print_scores(scores)
  a,b =  backtrack(scores, bestScore, s1, s2)
  print a
  print b
  return (a,b)
