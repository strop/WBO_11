# -*- coding: utf-8 -*-
from Bio import AlignIO
from Bio.SubsMat import MatrixInfo
from Bio.Phylo.BaseTree import Tree, Clade
from Bio import Phylo
import urllib2

def_format = "stockholm"
def_family = "PF05371"
def_subst_mat = MatrixInfo.blosum62

def print_matrix(m):
  for i in range(len(m)):
    print m[i]

def fetch_family(fam_id = def_family):
  url_template = "http://pfam.sanger.ac.uk/family/{0}/alignment/full/format?format={1}&alnType=full&order=t&case=l&gaps=dashes&download=0"
  url = url_template.format(fam_id, def_format)
  handle = urllib2.urlopen(url)
  return handle

def parse_family(handle):
  alignment = AlignIO.read(handle, def_format)
  return [x.seq for x in alignment]

""" oblicza punktację dla podanych dwóch aminokwasów (tak jakby miały być użyte na tej samej pozycji w uliniowieniu"""
def point_score(a1, a2, gap_penalty = -1, subst_mat = def_subst_mat):
  if a1 == '-' and a2 == '-':  # ignorujemy dwie spacje na raz
    pass
  elif a1 == '-' or a2 == '-': # kara za przerwę
    score = gap_penalty
  else:                                  # w wypadku, gdy są dwa aminokwasy punktujemy według macierzy substytucji
  try:
    score = subst_mat[seq1[c], seq2[c]]
  except KeyError: # czasem zdarza się... Uciekamy się do karania jak za przerwę
    score = gap_penalty
  return score
  
""" oblicza odległość dwóch sekwencji """
def calc_seq_distance(seq1, seq2, gap_penalty = -1, subst_mat = def_subst_mat):
  assert(len(seq1) == len(seq2))
  distance = 0
  for c in range(len(seq1)):               # porównujemy wszystkie pozycje po kolei:
  return distance;

def calc_distance_matrix(seqs, gap_penalty = -1, subst_mat = def_subst_mat):
  d = [[0 for x in seqs] for x in seqs]
  max_similarity = 0
  for i in range(len(seqs)):
    for j in range(i, len(seqs)):
      d[j][i] = calc_seq_distance(seqs[i], seqs[j])
      if d[j][i] > max_similarity:
        max_similarity = d[j][i]
  for i in range(len(seqs)):
    for j in range(i, len(seqs)):
      d[j][i] = max_similarity - d[j][i]
  print_matrix(d)
  return d

def remove_empty_columns(profile):
  mutables = [ s.to_mutable() for s in profile ]  # zamiana na mutable, żeby można było usuwac niektóre pozycje
  gap_chars = ['-', '_']
  for i in range(min([len(s) for s in profile])): # dla każdej pozycji sprawdzamy,
    if all((s[x] in gap_chars) for s in profile): # czy wszystkie sekwencje mają na tej pozycji przerwę
      for j in range(len(profile)):               # wtedy usuwamy kolumnę
        del profile[j][i]
  return profile

def back_track(scores, best_score_pos, profile1, profile2):
  return profile1 + profile2

""" True <=> lista w parametrze składa się z list o równej długości """
equal_lengths = lambda list_of_lists: all( len(one_list) == len(list_of_lists[0]) for one_list in list_of_lists )

class Direction:
  Left, Top, Diag = range(3)

def diagonal_score(scores, x, y, profile1, profile2, gap_penalty = -1, subst_mat = def_subst_mat):
  score = 0
  for s1 in profile1:
    for s2 in profile2:
      score = score + point_score(s1[x-1], s2[y-1], gap_penalty, subst_mat)
  return score

def horizontal_score(scores, x, y, profile1, profile2, gap_penalty = -1, subst_mat = def_subst_mat):
  return gap_penalty * len(profile1) * len(profile2) # każdy z każdym (opcja: pomijać dwie spacje?)

vertical_score = horizontal_score

""" profile - listy sekwencji. obliczamy punktację każdej pozycji na podstawie sumy punktów ze wszystkimi innymi sekwencjami """
def align_profiles(profile1, profile2, gap_penalty = -1, subst_mat = def_subst_mat):
  profile1 = remove_empty_columns(profile1)
  profile2 = remove_empty_columns(profile2)
  assert(equal_lengths(profile1))
  assert(equal_lengths(profile2))
  scores = [ [(0, Diag)] * len(profile2[0])+1 ] * (len(profile1[0])+1)
  best_score = (0,0)
  for y in range(1, len(profile2[0])+1):
    for x in range(1, len(profile1[0])+1):
      diag = (diagonal_score(scores, x, y, profile1, profile2, gap_penalty, subst_mat), Direction.Diag)
      horiz = (horizontal_score(scores, x, y, profile1, profile2, gap_penalty, subst_mat), Direction.Left)
      vert = (vertical_score(scores, x, y, profile1, profile2, gap_penalty, subst_mat), Direction.Top)
      best = max(diag, horiz, vert)

  
  return back_track(scores, best_score, profile1, profile2) # lista uliniowionych sekwencji

def init_tree():
  pass

def main(prot_family = def_family, subst_mat = def_subst_mat):
  initial_alignment = parse_family(fetch_family(prot_family))
  print_matrix(initial_alignment)
  distances = calc_distance_matrix(initial_alignment)
  max_distance = max([max(r) for r in distances])
  nodes = [ Clade(0, x) for x in range(len(initial_alignment) ]
  print max_distance
