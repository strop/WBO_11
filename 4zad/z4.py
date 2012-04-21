# -*- coding: utf-8 -*-
from Bio import SeqIO
from Bio import Alphabet
import Bio.HMM.MarkovModel as MM
import pickle
from copy import deepcopy

class Dumpfiles:
  promotor_freq = "promotor_freq.dump"
  genome_freq = "genome_freq.dump"
  viterbi_result = "viterbi.dump"

def prefix_tree(depth, initial_value = 0, symbols = ['A','T','G','C']):
  ''' Tworzy drzewko prefiksowe, które posłuży do zliczania podsłów sekwencji. (Użycie pythonowego substringa s[1:4] zużywa mnóstwo pamięci) '''
  symbol_tree = {}
  for s in symbols:
    symbol_tree[s] = initial_value
  for d in range(depth-1):
    new_root = {}
    for s in symbols:
      new_root[s] = deepcopy(symbol_tree)
    symbol_tree = new_root
  return symbol_tree

def fold_prefix_tree(tree, depth):
  ''' Przywraca drzewko prefiksowe do 'normalnej' postaci - jednopoziomowego słownika <sekwencja> : <częstotliwość> '''
  for d in range(depth-1):
    tree2 = {}
    for k1 in tree.keys():
      for k2 in tree[k1].keys():
        tree2[k1+k2] = tree[k1][k2]
    tree = tree2
  return tree

symbol_lengths = [1,2,3,4]

def calc_symbol_freq(seqs):
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

def dump_obj(obj, filename):
  dump_file = open(filename, "w")
  pickle.dump(obj, dump_file)

def load_obj(filename):
  return pickle.load(open(filename, 'r'))
  
def main_build_markov(promotor_filename = "promotor.fa", genome_filename = "genom.fa", symbol_length = 2, load_cached = False, save_cache = True):
  ''' Na podstawie plików z sekwencjami promotorowymi i genomem funkcja buduje model Markova'''
  promotor_sequences = [ x for x in SeqIO.parse("promotor.fa", "fasta")]
  genome = [ x for x in SeqIO.parse("genom.fa", "fasta")]
  if not load_cached:
    promotor_freqs = calc_symbol_freq(promotor_sequences)
    genome_freqs = calc_symbol_freq(genome)
    if save_cache:
      dump_obj(promotor_freqs, Dumpfiles.promotor_freq)
      dump_obj(genome_freqs, Dumpfiles.genome_freq)
  else:
    promotor_freqs = load_obj(Dumpfiles.promotor_freq)
    genome_freqs = load_obj(Dumpfiles.genome_freq)
  
  promotor_counts = calc_counts(promotor_sequences)
  genome_counts = calc_counts(genome)

  print promotor_counts

  promotor_freqs = fold_and_normalize(promotor_freqs[symbol_length], symbol_length, promotor_counts[symbol_length])
  genome_freqs = fold_and_normalize(genome_freqs[symbol_length], symbol_length, genome_counts[symbol_length])

  for k in promotor_freqs:
    assert(k in genome_freqs)
  for k in genome_freqs:
    assert(k in genome_freqs)

  print promotor_freqs
  (markov, states) = build_markov(genome_freqs, promotor_freqs)

  return (markov, states)

def main_get_promotor_frags(markov_model, states, symbol_length = 2, genome_filename = "genom.fa", load_cache = False, save_cache = True):
  ''' Na podstawie danego modelu Markova funkcja wyznacza prawdopodobne fragmenty promotorowe w danym genomie
      markov_model - bipythonowy obiekt modelu markowa
      states - alfabet stanów modelu markowa
      symbol_length - długość fragmentu kodu DNA, który uważamy za jeden "symbol" (1,2 lub 3)
      genome_filename - plik z genomem w formacie fasta, w którym chcemy wykryć fragmenty promotorowe
      load_cache - jeśli true, to wynik działania algorytmu viterbiego zostanie załadowany z wczesniej przygotowanego pliku (patrz save_cache)
      save_cache - jeśli true, to wynik działania algorytmu viterbiego zostanie zapisany do pliku do późniejszego wykorzystania
      
      zwraca - listę par (indeks_początku_odcinka_promotorowego, indeks_końca_odcinka_promotorowego) '''
  if load_cache:
    viterbi = load_obj(Dumpfiles.viterbi_result)
  else:
    genome = [ x for x in SeqIO.parse(genome_filename, "fasta")]
    genome = genome[0]
    # obiekt przekazany do viterbiego musi być sekwencją symboli, w tym wypadku symbole sa wieloznakowe, więc trzeba podzielić genom na kawałki odpowiedniej długości
    genome_splitted = [0] * ((len(genome)/symbol_length))
    for i in range(0, len(genome)-symbol_length+1, symbol_length):
      chunk = ''.join([genome[j] for j in range(i,i+symbol_length)])
      genome_splitted[i/symbol_length] = chunk
    # właściwe uruchomienie alg. viterbiego
    viterbi = markov_model.viterbi(genome_splitted, states)
  if save_cache:
    dump_obj(viterbi, Dumpfiles.viterbi_result)
  viterbi_string = viterbi[0]
  state = viterbi_string[0]
  promotor_positions = []
  start = 0
  # tworzenie listy par indeksów (początek,koniec) fragmentu promotorowego - liniowe przejście po sekwencji stanół zwróconej przez alg. viterbiego
  for i in range(1,len(viterbi_string)):
    if viterbi_string[i] != state:
      print "y"
      if state == 'p':
        promotor_positions.append((start,i))
      state = viterbi_string[i]
      start = i
  return promotor_positions
