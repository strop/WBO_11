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

def fold_and_normalize(tree, depth, count):
  tree = fold_prefix_tree(tree, depth)
  return normalize_freqs(tree, count)

def build_markov(plain_freqs, promo_freqs):
  plain_promo_prob = 0.1
  promo_plain_prob = 0.1

  g_state = 'g'
  p_state='p'
  states = Alphabet.Alphabet()
  states.letters = [g_state, p_state]

  symbols = Alphabet.Alphabet()
  symbols.letters = plain_freqs.keys()

  b = MM.MarkovModelBuilder(states, symbols)
  #b.set_initial_probabilities({g_state : 0.5, p_state : 0.5})
  b.allow_all_transitions()
  b.set_equal_probabilities()
  b.set_transition_score(g_state, p_state, plain_promo_prob)
  b.set_transition_score(g_state, g_state, 1 - plain_promo_prob)
  b.set_transition_score(p_state, g_state, promo_plain_prob)
  b.set_transition_score(p_state, p_state, 1 - promo_plain_prob)
  for (s, p) in plain_freqs.items():
    b.set_emission_score(g_state, s, p)
  for (s, p) in promo_freqs.items():
    b.set_emission_score(p_state, s, p)
  return (b.get_markov_model(), states) 

''' Na podstawie danych plików z sekwencjami promotorowymi i danego genomu funkcja buduje model Markova'''
def main_build_Markov(promotor_filename = "promotor.fa", genome_filename = "genom.fa"):
  pass 

''' Na podstawie danego modelu Markova funkcja wyznacza prawdopodobne fragmenty promotorowe w danym genomie '''
def main_get_promotor_frags(markov_model, states, symbol_length = 2, genome_filename = "genom.fa"):
  genome = [ x for x in SeqIO.parse(genome_filename, "fasta")]
  genome = genome[0]
  prediction = markov_model.viterbi(genome, states)
  genome_splitted = [0] * (len(genome[0])/symbol_length)
  for i in range(0, len(genome)-1, symbol_length):
    genome_splitted[i/symbol_length] = (genome[i] + genome[i+1])
  viterbi = markov.viterbi(genome_splitted, states)
  pass

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

for k in promotor_freqs:
  assert(k in genome_freqs)
for k in genome_freqs:
  assert(k in genome_freqs)

print promotor_freqs

(markov, states) = build_markov(genome_freqs, promotor_freqs)

# Fragmenty odpowiedzialne za tłumaczenie wieloliterowych symboli emisji na symbole jednoznakowe
# okazały się niepotrzebne, bo Viterbi nie działał z powodu wieloznakowych nazw *stanów* a nie symboli

#def translate_symbols(freqs, translation):
#  assert(len(freqs) <= len(translation))
#  for k in freqs.keys():
#    freqs[translation[k]] = freqs[k]
#    del freqs[k]

#single_char_symbols = map(chr, range(ord('a'), ord('z')+1))
#single_char_trans = {}
#for (k, t) in zip(genome_freqs.keys(), single_char_symbols):
#  single_char_trans[k] = t

#translate_symbols(promotor_freqs, single_char_trans)
#translate_symbols(genome_freqs, single_char_trans)
#for i in range(len(genome_splitted)):
#  genome_splitted[i] = single_char_trans[genome_splitted[i]]
