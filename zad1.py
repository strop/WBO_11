from Bio import SeqIO

def gen_substrings(s, l):
  return [s[x:x+l] for x in range(0,len(s)-l+1)]

def gen_probes_inner(sequences, k):
  substrings = []
  uniques = []
  seq_count = len(sequences)
  for x in sequences:
    substrings.append(gen_substrings(x.seq, k))
  for i in range(seq_count):
    uniques.append(substrings[i])
    for j in range(seq_count):
      if j != i:
        uniques[i] = [x for x in uniques[i] if sequences[j].seq.find(x) == -1]
  probes = []
  try:
    for i in range(seq_count):
      probes.append(uniques[i][0])
    return probes
  except IndexError:
    return None

def gen_probes(fname = "test.fa", k = -1):
  sequences = [ x for x in SeqIO.parse(fname, "fasta")]
  probes = None
  if k == -1:
    k = 1
    while probes is None:
      probes = gen_probes_inner(sequences, k)
      k += 1
  else:
    probes = gen_probes_inner(sequences, k)
  return probes
