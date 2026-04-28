# Script to extract ORF data from FASTA file (here: AtRTDv2_QUASI_19April2016.fa)
# - determines transcripts that have longest orfs among all transcripts of a gene; if several transcripts with longest orfs exist, 
#   one is arbitrarily choosen
# - creates plain text file containing entries consisting of <Transcript ID>, <ORF start position>, <ORF stop position>
# Created by Christoph Engel, Research Group Wachter, Department of Biology, Mainz University

# Example command line:
# python getORFs.py AtRTDv2_QUASI_19April2016.fa longest_orfs_per_gene.txt

import csv, re, sys
from Bio import SeqIO

usage = "Usage: " + sys.argv[0] + " <path/transcriptome file (fa file)>" + " <path/longest orf file (txt file)>"
# transcriptome file: contains transcript sequences, from which longest orfs have to be retrieved
# longest orf file: filename of file for export of longest orfs per gen

if len(sys.argv) != 3:
  print(len(sys.argv))
  print(usage)
  sys.exit()

# Subroutine according to https://www.youtube.com/watch?v=uTlVNl3zTP8 ("Finding the Longest Open Reading Frame (ORF) in an RNA Sequence" by Professor Hendrix, Oregon State University)
def find_longest_ORF4transcript(RNA):
   ORFs = []
   if 'AUG' in RNA:
     for startMatch in re.finditer('AUG', RNA):
       remaining = RNA[startMatch.start():]
       for stopMatch in re.finditer('UAA|UGA|UAG', remaining):
         substring = remaining[:stopMatch.end()]
         if len(substring) % 3 == 0:
           ORFs.append(substring)
           # first stopcodon alone relevant because translation ends here
           break
   ORFs.sort(key=len, reverse=True)
   # avoid error at return if list kept empty
   if len(ORFs) == 0:
     ORFs.append('')
   return ORFs[0]


filename = sys.argv[1]
# SeqIO usage according to https://www.youtube.com/watch?v=Y9POfa_PH-M ("Introducing the Biopython SeqIO Module: Reading FASTA files" by Professor Hendrix, Oregon State University)
sequences = SeqIO.parse(filename, 'fasta')
L = []

for record in sequences:
  L.append(record)


#1. Run:
#-------
#- find longest ORFs per gen
#- save these in Dictionary longest_orfs4genes

def save_current_data(current_gene, maxlen_orf, longest_orfs4genes):
  longest_orfs4genes[current_gene] = maxlen_orf

current_gene = ''
maxlen_orf = 0
longest_orfs4genes = {}

for record_ in L:
  RNA = str(record_.seq.transcribe())
  ORF = find_longest_ORF4transcript(RNA)
  if ORF != '':
    # when switching to a new gene...
    if record_.id[0:9] != current_gene:
      save_current_data(current_gene, maxlen_orf, longest_orfs4genes)
      # ... save its name in current_gene
      current_gene = record_.id[0:9]
      # ... reset maxlen
      maxlen_orf = 0
    if len(ORF) > maxlen_orf:
      maxlen_orf = len(ORF)


#2. Run:
#-------
#- build list records_ with transcript isoforms, which contain longest ORFs per accompanying gen

def save_possibly_current_record(current_gene, longest_orfs4genes, ORF, record_, records_):
  if current_gene != '':
    try:
      current_longest_orf_length4gen = longest_orfs4genes[current_gene]
    except:
      print('Laengster ORF zu Gen ', current_gene, ' existiert nicht.')
    else:
      if current_longest_orf_length4gen == len(ORF):
        records_.append(record_)

current_gene = ''
records_ = []

for record_ in L:
  RNA = str(record_.seq.transcribe())
  ORF = find_longest_ORF4transcript(RNA)
  if ORF != '':
    save_possibly_current_record(record_.id[0:9], longest_orfs4genes, ORF, record_, records_)


#Export list of transcript IDs with longest ORFs per gene along with start and end coordinate of ORF (1-based coordinates);
#if several transcripts of a gene have equally long longest ORFs, then export only one (here: the last scrutinized):
#--------------------------------------------------------------------------------------------------------------------------
FileHandle1 = open(sys.argv[2], 'w')

writer1 = csv.writer(FileHandle1, delimiter=',')
current_gene = ''
prev_element = object()
i = 0
for element in records_:
  if element.id[0:9] != current_gene:
    current_gene = element.id[0:9]
    # export data of sole or of last scrutinized transcript with longest orf per gene
    if i > 0:
      RNA = str(prev_element.seq.transcribe())
      ORF = find_longest_ORF4transcript(RNA)
      for match in re.finditer(ORF, RNA):
        row_ = []
        row_.append(prev_element.id)
        row_.append(match.start()+1)
        row_.append(match.end())
        writer1.writerow(row_)
      FileHandle1.flush()
  prev_element = element
  i += 1