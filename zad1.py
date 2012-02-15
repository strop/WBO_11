from Bio import SeqIO

def gen_substrings(s, l):
  return [s[x:x+l] for x in range(0,len(s)-l+1)]

def gen_probes(fname = "test.fa", k = -1):
  substrings = []
  sequences = [ x for x in SeqIO.parse(fname, "fasta")]
  uniques = []
  seq_count = len(sequences)
  for x in sequences:
    substrings.append(set(gen_substrings(x.seq, k)))
  for i in range(seq_count):
    uniques.append(substrings[i])
    for j in [j for j in range(0,seq_count) if j != i]:
      uniques[i] = uniques[i] - substrings[j];
  for i in range(seq_count):
    print uniques[i].pop().reverse_complement()
