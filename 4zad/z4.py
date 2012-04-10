# -*- coding: utf-8 -*-
from Bio import SeqIO
from Bio import Alphabet
import Bio.HMM.MarkovModel as MM
import pickle
from copy import deepcopy

class Dumpfiles:
  promotor_freq = "promotor_freq.dump"
  genome_freq = "genome_freq.dump"

''' Tworzy drzewko prefiksowe, które posłuży do zliczania podsłów sekwencji. (Użycie pythonowego substringa s[1:4] zużywa mnóstwo pamięci) '''
def prefix_tree(depth, initial_value = 0, symbols = ['A','T','G','C']):
  symbol_tree = {}
  for s in symbols:
    symbol_tree[s] = initial_value
  for d in range(depth-1):
    new_root = {}
    for s in symbols:
      new_root[s] = deepcopy(symbol_tree)
    symbol_tree = new_root
  return symbol_tree

''' Przywraca drzewko prefiksowe do 'normalnej' postaci - jednopoziomowego słownika <sekwencja> : <częstotliwość> '''
def fold_prefix_tree(tree, depth):
  for d in range(depth-1):
    tree2 = {}
    for k1 in tree.keys():
      for k2 in tree[k1].keys():
        tree2[k1+k2] = tree[k1][k2]
    tree = tree2
  return tree

symbol_lengths = [1,2,3,4]

def calc_symbol_freq(seqs, dump_filename):
  freqs = [0,0,0,0]
  for l in range(1,4):
    symbol_tree = prefix_tree(l)
    normalization = sum([len(s)-l for s in seqs])
    f = {}
    for s in seqs:
      for i in range(len(s) - l):
	curr_dict = symbol_tree
	for j in range(i, i+l-1):
	  curr_dict = curr_dict[s[j]]
        curr_dict[s[i+l-1]] = curr_dict[s[i+l-1]] + 1
    freqs[l] = symbol_tree
  dump_file = open(dump_filename, "w")
  pickle.dump(freqs, dump_file)
  return freqs

def calc_counts(seqs):
  counts = [0] * (len(symbol_lengths)+1)
  for i in symbol_lengths:
    counts[i] = sum([(len(s) - i + 1) for s in seqs])
  return counts

def normalize_freqs(freqs, count):
  for k in freqs.keys():
    freqs[k] = float(freqs[k]) / float(count)
  return freqs

def load_freqs(filename):
  return pickle.load(open(filename, 'r'))

''' Na podstawie danych plików z sekwencjami promotorowymi i danego genomu funkcja buduje model Markova'''
def main_build_Markov(promotor_filename = "promotor.fa", genome_filename = "genom.fa"):
  pass 

''' Na podstawie danego modelu Markova funkcja wyznacza prawdopodobne fragmenty promotorowe w danym genomie '''
def main_get_promotor_frags(markov_model, genome_filename = "genom.fa"):
  pass

def fold_and_normalize(tree, depth, count):
  tree = fold_prefix_tree(tree, depth)
  return normalize_freqs(tree, count)



promotor_sequences = [ x for x in SeqIO.parse("promotor.fa", "fasta")]
genome = [ x for x in SeqIO.parse("genom.fa", "fasta")]
#promotor_freqs = calc_symbol_freq(promotor_sequences, Dumpfiles.promotor_freq)
#genome_freqs = calc_symbol_freq(genome, Dumpfiles.genome_freq)

promotor_freqs = load_freqs(Dumpfiles.promotor_freq)
genome_freqs = load_freqs(Dumpfiles.genome_freq)
promotor_counts = calc_counts(promotor_sequences)
genome_counts = calc_counts(genome)

print promotor_counts

symbol_length = 2

promotor_freqs = fold_and_normalize(promotor_freqs[symbol_length], symbol_length, promotor_counts[symbol_length])

genome_freqs = fold_and_normalize(genome_freqs[symbol_length], symbol_length, genome_counts[symbol_length])

def build_markov(plain_freqs, promo_freqs):
  plain_promo_prob = 0.1
  promo_plain_prob = 0.1

  states = Alphabet.Alphabet()
  states.letters = ["PLAIN", "PROMO"]

  symbols = Alphabet.Alphabet()
  symbols.letters = plain_freqs.keys()

  b = MM.MarkovModelBuilder(states, symbols)
  #b.set_initial_probabilities({"PLAIN" : 0.5, "PROMO" : 0.5})
  b.allow_all_transitions()
  b.set_equal_probabilities()
  b.set_transition_score("PLAIN", "PROMO", plain_promo_prob)
  b.set_transition_score("PLAIN", "PLAIN", 1 - plain_promo_prob)
  b.set_transition_score("PROMO", "PLAIN", promo_plain_prob)
  b.set_transition_score("PROMO", "PROMO", 1 - promo_plain_prob)
  for (s, p) in plain_freqs.items():
    b.set_emission_score("PLAIN", s, p)
  for (s, p) in promo_freqs.items():
    b.set_emission_score("PROMO", s, p)
  return (b.get_markov_model(), states, symbols) 

(markov, states, symbols) = build_markov(genome_freqs, promotor_freqs)
genome_splitted = [0] * (len(genome[0])/2)
for i in range(0, len(genome[0])-1, 2):
  genome_splitted[i/2] = (genome[0][i] + genome[0][i+1])
genome[0].alphabet = symbols
markov.viterbi(genome_splitted, states)
