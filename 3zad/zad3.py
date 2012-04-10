# -*- coding: utf-8 -*-
from Bio import AlignIO
import sys
from Bio.Seq import Seq, MutableSeq
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
  score = 0
  if a1 == '-' and a2 == '-':  # ignorujemy dwie spacje na raz
    pass
  elif a1 == '-' or a2 == '-': # kara za przerwę
    score = gap_penalty
  else:                                  # w wypadku, gdy są dwa aminokwasy punktujemy według macierzy substytucji
    try:
      score = subst_mat[a1, a2]
    except KeyError: # czasem zdarza się... Uciekamy się do karania jak za przerwę
      score = gap_penalty
  return score
  
""" oblicza odległość dwóch sekwencji """
def calc_seq_distance(seq1, seq2, gap_penalty = -1, subst_mat = def_subst_mat):
  assert(len(seq1) == len(seq2))
  distance = 0
  for c in range(len(seq1)):               # porównujemy wszystkie pozycje po kolei:
    distance = distance + point_score(seq1[c], seq2[c], gap_penalty, subst_mat)
  return distance;

def calc_distance_matrix(seqs, gap_penalty = -1, subst_mat = def_subst_mat):
  d = [[0 for x in seqs] for x in seqs]
  max_similarity = 0
  for i in range(len(seqs)):
    for j in range(i, len(seqs)):
      d[j][i] = calc_seq_distance(seqs[i], seqs[j])
      d[i][j] = d[j][i]
      if d[j][i] > max_similarity:
        max_similarity = d[j][i]
  for i in range(len(seqs)):
    for j in range(len(seqs)):
      d[j][i] = max_similarity - d[j][i]
  return d

def closest_neighbors(dist_matrix):
  d = dist_matrix
  r = len(d)
  min_pos = (1,0)
  min_val = (r-2)*d[1][0] - sum(d[1]) - sum(d[0])
  for i in range(1,r):
    for j in range(i+1, r):
      q = ((r-2)*d[i][j] - sum(d[i]) - sum(d[j])) 
      if q < min_val:
        min_val = q
	min_pos = (i, j)
  return min_pos

""" aktualizuje macierz odległości, aby zawierała prawidłowe informacje dotyczące węzła utworzonego przez połączenie sąsiadów
    odległości do nowego węzła znajdą sie w ostatnim wierszu wynikowej macierzy, a wiersze dotyczące sąsiadów zostaną z niej usunięte """
def update_dist_matrix(dist_matrix, neighbors):
  d = dist_matrix
  assert(equal_lengths(d))
  x = neighbors[0]
  y = neighbors[1]
  newdist = [0] * len(d)
  for i in range(len(d)): # liczenie odległości do nowego, połączonego węzła
    newdist[i] = 0.5 * (d[x][i] + d[y][i] - d[x][y])
  del d[max(x,y)]
  del d[min(x,y)]
  del newdist[max(x,y)]
  del newdist[min(x,y)]
  for i in range(len(d)):
    del d[i][max(x,y)]
    del d[i][min(x,y)]
    d[i].append(newdist[i])
  newdist.append(0)
  d.append(newdist) # dodajemy odległości do nowego węzła
  return d

def join_neighbors(dist_matrix, neighbors, nodes):
  d = dist_matrix
  x = neighbors[0]
  y = neighbors[1]
  # odległości łączonych węzłów do nowego wierzchołka
  x_dist = 0.5 * (d[x][y] + (sum(d[x]) - sum(d[y])) / (len(d) - 2))
  y_dist = d[x][y] - x_dist
  nodes[x].branch_length = x_dist
  nodes[y].branch_length = y_dist
  newnode = Clade(0, None, [nodes[x], nodes[y]])
  del nodes[max(x,y)]
  del nodes[min(x,y)]
  nodes.append(newnode)
  return (nodes, max(x_dist, y_dist)) # długości gałęzi przydatne do poszukiwania najdłuzszej

def remove_empty_columns(profile):
  mutables = [ s.tomutable() for s in profile ]  # zamiana na mutable, żeby można było usuwac niektóre pozycje
  gap_chars = ['-', '_']
  for i in reversed(range(min([len(s) for s in profile]))): # dla każdej pozycji sprawdzamy,
    if all((s[i] in gap_chars) for s in profile): # czy wszystkie sekwencje mają na tej pozycji przerwę
      for j in range(len(profile)):               # wtedy usuwamy kolumnę
        del mutables[j][i]
  return mutables

def back_track(scores, profile1, profile2):
  x = len(profile1[0])
  y = len(profile2[0])
  alignment1 = [ MutableSeq("") for z in profile1 ]
  alignment2 = [ MutableSeq("") for z in profile2 ]
  while x > 0 or y > 0:
    #print (x,y)
    direction = scores[x][y][1]
    if direction in [Direction.Diag, Direction.Left] :
      for i in range(len(profile1)): # dodajemy kolumnę z profile1
        alignment1[i].insert(0, profile1[i][x-1])
      x = x-1
    if direction in [Direction.Diag, Direction.Top] :
      for i in range(len(profile2)): # dodajemy kolumnę z profile2
        alignment2[i].insert(0, profile2[i][y-1])
      y = y-1
    if direction == Direction.Top:
      for s in alignment1: # dodajemy kolumnę przerw do alignment1
        s.insert(0, '-')
    if direction == Direction.Left:
      for s in alignment2: # dodajemy kolumnę przerw do alignment2
        s.insert(0, '-')
    #print_matrix(alignment2)
    #print "xOxOx"
  return (alignment1, alignment2)

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
  print "align_profiles1: "
  print_matrix(profile1)
  print "a_p2: "
  print_matrix(profile2)
  print ""
  assert(equal_lengths(profile1))
  assert(equal_lengths(profile2))
  # tablica punktacji zawiera elementy: (punkty, kierunek)
  scores = [ [(0, Direction.Left)] * (len(profile2[0])+1) ] * (len(profile1[0])+1)
  scores[0] = [(0, Direction.Top)] * (len(profile2[0])+1)
  for y in range(1, len(profile2[0])+1):
    for x in range(1, len(profile1[0])+1):
      diag = (diagonal_score(scores, x, y, profile1, profile2, gap_penalty, subst_mat), Direction.Diag)
      horiz = (horizontal_score(scores, x, y, profile1, profile2, gap_penalty, subst_mat), Direction.Left)
      vert = (vertical_score(scores, x, y, profile1, profile2, gap_penalty, subst_mat), Direction.Top)
      best = max(diag, horiz, vert)
      scores[x][y] = best
  #print_matrix(scores)
  (a1, a2) = back_track(scores, profile1, profile2) # lista uliniowionych sekwencji
  return a1 + a2

""" zwraca listę identyfikatorów liści drzewa Bio.Phylo """
def tree_2_leaves_ids(tree):
  return [int(x.name) for x in tree.get_terminals()]

""" główna funkcja - parametry tokod rodziny białek i macierz substytucji """
def main(prot_family = def_family, subst_mat = def_subst_mat):
  alignment = parse_family(fetch_family(prot_family))

  prev_max_branch_len = sys.maxint  
  while True:
    distances = calc_distance_matrix(alignment)
    max_branch_len = 0
    max_branch_parent = None # z tego wierzchołka wychodzi najdłuższa gałąź
    branch_len = 0
    nodes = [ Clade(0, str(x)) for x in range(len(alignment)) ]
    while len(nodes) > 2: # neighbor joining
      neighbors = closest_neighbors(distances)
      assert(neighbors[0] != neighbors[1])
      (nodes, branch_len) = join_neighbors(distances, neighbors, nodes)
      if branch_len > max_branch_len:
        max_branch_len = branch_len
        max_branch_parent = nodes[-1]
      distances = update_dist_matrix(distances, neighbors)
      assert(equal_lengths(distances))
      assert(len(distances) == len(nodes))
    # łączymy dwa ostatnie poddrzewa
    nodes[0].branch_length = distances[0][1] 
    nodes[1].clades.append(nodes[0])
    if distances[0][1] > max_branch_len:
      max_branch_len = distances[0][1]
      max_branch_parent = nodes[-1]
    if max_branch_len <= prev_max_branch_len:
      break
    prev_max_branch_len = max_branch_len
    # dzielimy drzewo na dwa względem najdłuższej ścieżki (wcześniej znalezionej)
    profile_1_root = None
    profile_2_root = None
    for i in range(len(max_branch_parent.clades)):
      if max_branch_parent.clades[i].branch_length == max_branch_len:
        profile_1_root = max_branch_parent.clades[i] # "korzeń" pierwszego poddrzewa
        del max_branch_parent.clades[i]
        profile_2_root = nodes[1]                    # "korzeń" drugiego poddrzewa
        break
    # separujemy profile
    profile1 = [ alignment[i] for i in tree_2_leaves_ids(profile_1_root) ]
    profile2 = [ alignment[i] for i in tree_2_leaves_ids(profile_2_root) ]
    # uliniawiamy
    alignment = align_profiles(profile1, profile2, -1, subst_mat)
  return alignment
