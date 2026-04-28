# Script to extract exon data from GTF transcriptome file (here: AtRTDv2_QUASI_19April2016.gtf)
# - creates pkl file containing dictionary with transcript IDs as keys and lists of lists of exon and intron lengths as values
# Created by Christoph Engel, Research Group Wachter, Department of Biology, Mainz University

# Example command line:
# python getExondata.py AtRTDv2_QUASI_19April2016.gtf transcriptModels_AtRTDv2_QUASI_19April2016.pkl

import sys, re, pickle

usage = "Usage: " + sys.argv[0] + " <path/gtf file>" + " <path/exon file (pkl file)>"

if len(sys.argv) != 3:
  print(len(sys.argv))
  print(usage)
  sys.exit()

# Subroutine according to https://www.youtube.com/watch?v=4eyQrYGuUuk ("Extracting transcripts from a GTF file and genome FASTA file using basic python" by Professor Hendrix, Oregon State University)
def read_gtf_file(gtf_file):
# read data from gtf file; extract transcript-ID, exon start, exon end and strand; save extracted data as tuple and put tuple into dictionary
# return dictionary with exon data: {<transcript ID 1>:<exon list 1>, <transcript ID 2>:<exon list 2>, ...}
  data = {}
  GTF = open(gtf_file)
  lines_in_gtf = 0
  exons = 0
  exon = ()
  transcript_id = ''
  transcript_counter = 0
  for line in GTF:
    lines_in_gtf += 1
    chrom,source,seqtype,start,end,score,strand,frame,attr = line.strip().split('\t')
    if seqtype == 'exon':
      exons += 1
      if "transcript_id" in attr:
        match = re.search('transcript_id "(.*?)";',attr)
        transcript_id = match.group(1)
        if transcript_id not in data:
          data[transcript_id] = []
          transcript_counter += 1
        # gtf file coordinates: 1-based + inclusive - let it so
        exon = (chrom, int(start), int(end), strand)
        data[transcript_id].append(exon)
      else:
        print("error parsing GTF")
        exit()
  print(f"lines_in_gtf: {lines_in_gtf}, Exons: {exons}, Transcripts: {transcript_counter}")
  return data


def check_same_exon_start_for_txs_of_gene(transcriptModels, blacklist):
# transcriptModels: dictionary with transcript IDs as keys and lists of raw exon data (chromosome, start coordinate, end coordinate, strand) as values
# blacklist: contains genes for which transcripts with different start coordinates of the first exon exist - to be filled here
  gene_name = ''
  first_exon_start = -1
  for key, value in transcriptModels.items():
    if gene_name != key[0:9]:
      # first tx of new gene is reference
      if value[0][3] == "+":
        # plus strand transcript
        first_exon_start = value[0][1]
      else:
        # minus strand transcript
        first_exon_start = value[len(value)-1][2]
      gene_name = key[0:9]
    else:
      # another tx of same gene
      if value[0][3] == "+":
        # plus strand transcript
        if first_exon_start != value[0][1]:
          #print(f"Gene {key[0:9]}: Start coordinates of transcripts differ. ({first_exon_start} != {value[0][1]})(+)")
          if key[0:9] not in blacklist:
            blacklist.append(key[0:9])
      else:
        # minus strand transcript
        if first_exon_start != value[len(value)-1][2]:
          #print(f"Gene {key[0:9]}: Start coordinates of transcripts differ. ({first_exon_start} != {value[len(value)-1][2]})(-)")
          if key[0:9] not in blacklist:
            blacklist.append(key[0:9])


def get_edited_transcript_data(transcriptModels, blacklist):
# transcriptModels: dictionary with transcript IDs as keys and lists of raw exon data (chromosome, start coordinate, end coordinate, strand) as values
# blacklist: contains gene names for which transcripts with different start coordinates of the first exon exist
# from absolute exon coordinates determine exon and intron lengths
# rearrange exon and intron data for minus strand transcripts (invert sequence)
# drop transcripts which belong to genes in blacklist (31 genes concerned)
# return dictionary containing transcript IDs as keys and lists of lists of exon and intron lengths as values
  prev_transcript = ""
  curr_transcript = ""
  exon_length_list = []
  intron_length_list = []
  endcoord_predecessor_exon = 0
  transcriptModels1 = {}
  for key, value in transcriptModels.items():
    if key[0:9] in blacklist:
      continue
    else:
      # determine lengths of exons and introns
      for exon in value:
        exon_length = exon[2]-exon[1]+1
        exon_length_list.append(exon_length)
        if endcoord_predecessor_exon != 0:
          intron_length = exon[1]-endcoord_predecessor_exon-1
          intron_length_list.append(intron_length)
        endcoord_predecessor_exon = exon[2]
      # check first element of exon list, to which strand the transcript belongs
      if value[0][3] == "-":
        # minus strand: transcription in reverse direction -> reverse exon and intron length list
        exon_length_list.reverse()
        intron_length_list.reverse()
      exon_intron_lengths = [exon_length_list, intron_length_list]
      transcriptModels1[key] = exon_intron_lengths
      exon_length_list = []
      intron_length_list = []
      endcoord_predecessor_exon = 0
  return transcriptModels1


blacklist = []

transcriptModels = read_gtf_file(sys.argv[1])
check_same_exon_start_for_txs_of_gene(transcriptModels, blacklist)
transcriptModels1 = get_edited_transcript_data(transcriptModels, blacklist)
# save exon data in a pkl file:
with open(sys.argv[2], 'wb') as fp:
  pickle.dump(transcriptModels1, fp)