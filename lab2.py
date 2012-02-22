from Bio import Entrez
from Bio import SeqIO
from Bio.Data import CodonTable
from random import Random

def get_amino_to_codon_dict():
  forward = CodonTable.unambiguous_dna_by_id[1].forward_table
  reverse = {}
  for k in forward.keys():
    if not forward[k] in reverse:
      reverse[forward[k]] = []
    reverse[forward[k]].append(k)
  print reverse
  return reverse

def get_random_mRNA_from_protein_sequence(pseq):
  reverse = get_amino_to_codon_dict()
  # generujemy tablice znakow
  pchars = [x for x in pseq]
  # obcinamy po pierwszym stopie
  pchars = pchars[0:pchars.index('*')]
  r = Random()
  rand_mRNA = ''.join([r.choice(reverse[x]) for x in pchars ])
  return rand_mRNA


Entrez.email = "A.N.Other@example.com"
handle = Entrez.efetch(db="nucleotide", rettype="fasta", id="6273291")
seq_record = SeqIO.read(handle, "fasta")
handle.close()
print "%s with %i features" % (seq_record.id, len(seq_record.features))
print seq_record

print get_random_mRNA_from_protein_sequence(seq_record.seq.transcribe().translate())
