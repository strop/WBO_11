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

def calc_seq_distance(seq1, seq2, gap_penalty = -1, subst_mat = def_subst_mat):
  assert(len(seq1) == len(seq2))
  distance = 0
  for c in range(len(seq1)):
    if seq1[c] == '-' and seq2[c] == '-':
      pass
    elif seq1[c] == '-' or seq1[c] == '-':
      distance = distance + gap_penalty
    else:
      try:
        distance = distance + subst_mat[seq1[c], seq2[c]]
      except KeyError:
        distance = distance + gap_penalty
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

""" profile - listy sekwencji. obliczamy punktację każdej pozycji na podstawie sumy punktów ze wszystkimi innymi sekwencjami """
def align_profiles(profile1, profile2):
  pass

def init_tree():
  pass

def main(prot_family = def_family, subst_mat = def_subst_mat):
  initial_alignment = parse_family(fetch_family(prot_family))
  print_matrix(initial_alignment)
  distances = calc_distance_matrix(initial_alignment)
  max_distance = max([max(r) for r in distances])
  nodes = [ Clade(0, x) for x in range(len(initial_alignment) ]
  print max_distance

