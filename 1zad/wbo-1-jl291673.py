# WBO 11/12 - Zad 1
# Autor: Jarek Łukow (291673)

from Bio import SeqIO

# Uzycie:
# glowna funkcja to gen_probes(fileHandle, k, filtr)
# fileHandle powinien byc uchwytem do pliku (otrzymanym np. przez open("test.fa"))
# k okresla dlugosc wszystkich sekwencji w rozwiazaniu
# filtr to funkcja Seq -> Bool, ktora zwraca True wtw podana w parametrze sekwencja moze byc uzyta w rozwiazaniu
# Funkcja rzuca wyjatek w wypadku braku rozwiazan
#
# parametry k i filtr sa opcjonalne:
# jesli podano tylko k, funkcja zwroci sekwencje dlugosci dokladnie k
# jesli nie podano k, funkcja bedzie szukac rozwiazan dla kolejnych k zaczynajac od 1 i zwroci pierwsze znalezione rozwiazanie
# jesli podano filtr i k, to funkcja zwroci sekwencje dlugosci dokladnie k, ktore zostaly zaakceptowane przez filtr
# jesli podano filtr, ale nie podano k, funkcja zwroci rozwiazania wyszukane sposrod wszystkich podslow podanych sekwencji (to znaczy wszystkie sekwencje w rozwiazaniu beda akceptowane przez filtr, ale byc moze nie wszystkie beda tej samej dlugosci)

# funkcja pomocnicza generujaca podslowa dlugosci l
def gen_substrings(s, length):
  return [s[x:x+length] for x in range(0,len(s)-length+1)]

# funkcja pomocnicza generujaca wszystkie podslowa
def gen_all_substrings(s, dummyParam):
  return [ s[x:x+l] for x in range(len(s)) for l in range(1, len(s)+1) if x+l < len(s)+1 ]

# funkcja pomocnicza generujaca rozwiazania, mniej czula na brak parametrow
# generator podslow jest parametrem, bo czasem chcemy miec wszystkie podslowa, a czasem tylko te okreslonej dlugosci
def gen_probes_inner(sequences, k, filtr, substrings_generator):
  substrings = []
  uniques = []
  seq_count = len(sequences)
  for i in range(seq_count):
    uniques.append([x for x in substrings_generator(sequences[i].seq, k) if filtr(x)])
    for j in range(seq_count):
      if j != i:
        uniques[i] = [x for x in uniques[i] if sequences[j].seq.find(x) == -1] # eliminujemy fragmenty ktore sa tez podslowami w innej sekwencji
  probes = []
  for i in range(seq_count):
    if len(uniques[i]) == 0: # jakas sekwencja nie posiada unikalnego fragmentu dlugosci k
      return None            # wiec zwracamy porazke
    probes.append(uniques[i][0])
  return probes

def gen_probes(fileHandle, k = None, filtr = None):
  assert(type(fileHandle) == file)
  sequences = [ x for x in SeqIO.parse(fileHandle, "fasta")]
  defFiltr = lambda(x): True
  probes = None
  if k is None:
    if (filtr is None) or (type(filtr) != type(defFiltr)): # nie mamy ani k ani filtra, wiec szukamy rozwiazania dla najmniejszego k
      filtr = defFiltr
      # wyznaczamy dlugosc najkrotszej sekwencji, zeby wiedziec kiedy zakonczyc przeszukiwanie przy nieokreslonym k
      minLen = len(sequences[0])
      for i in range(len(sequences)):
        if len(sequences[i]) < minLen:
          minLen = len(sequences[i])
      k = 1
      # poszukujemy coraz dluzszych rozwiazan
      while (probes is None) and (k < minLen):
        probes = gen_probes_inner(sequences, k, filtr, gen_substrings)
        k += 1
    else: # nie mamy k, ale mamy filtr - rozwiazania moga miec rozna dlugosc, ale musza byc zaakceptowane przez filtr
      probes = gen_probes_inner(sequences, k, filtr, gen_all_substrings)
  else: # jesli okreslono k, to szukamy tylko sekwencji dlugosci k, jesli podano filtr, to dodatkowo je filtrujemy
    if (filtr is None) or (type(filtr) != type(defFiltr)):
      filtr = defFiltr
    probes = gen_probes_inner(sequences, k, filtr, gen_substrings)
  if probes is None:
    raise Exception("Nie istnieje rozwiazanie")
  for i in range(len(sequences)):
    sequences[i].seq = probes[i].reverse_complement()
  return sequences
