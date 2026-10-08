# Script to locate Alternative Splicing events in transcripts
# Created by Christoph Engel, Research Group Wachter, Department of Biology, Mainz University

# Example command line:
# python locateASEvents.py longest_orfs_per_gene.txt transcriptModels_AtRTDv2_QUASI_19April2016.pkl output_as_event_locations.txt

# longest orf file:
# - created by getORFs.py
# - content:
#   - lines consisting of transcript IDs, ORF start, ORF end
# - example content:
#   AT1G01010.1,130,1419
#   AT1G01020_P1,465,1202
#   AT1G01030_ID1,775,1851

# exon file:
# - created by getExondata.py
# - content after read in:
#   - dictionary containing as keys transcript IDs and as values lists with a list of exon lengths as the first element and a list of 
#     intron lengths as the second element
# - example content:
#   {'AT1G01010.1': [[283, 281, 120, 390, 153, 461], [82, 209, 100, 78, 112]], 'AT1G01020_P1': [[560, 48, 90, 46, 74, 86, 67, 76, 282], 
#   [106, 91, 248, 106, 112, 113, 151, 87]], 'AT1G01020_P2': [[560, 48, 90, 46, 74, 86, 294, 282], [106, 91, 248, 106, 112, 113, 87]], 
#   'AT1G01020_P3': [[560, 229, 46, 74, 86, 67, 76, 282], [106, 248, 106, 112, 113, 151, 87]], 'AT1G01020_P4': [[537, 229, 46, 74, 86, 67, 76, 282], 
#   [129, 248, 106, 112, 113, 151, 87]], 'AT1G01020_P5': [[537, 48, 90, 46, 74, 86, 67, 76, 282], [129, 91, 248, 106, 112, 113, 151, 87]], 
#   'AT1G01020_P6': [[537, 48, 90, 86, 294, 282], [129, 91, 586, 113, 87]], 'AT1G01030.1': [[380, 1525], [161]], 'AT1G01030_ID1': [[2066], []], 
#   'AT1G01030_P2': [[380, 750, 706], [161, 69]]}

import sys, pickle

usage = "Usage: python" + sys.argv[0] + " <path/longest orf file (txt file)>" + " <path/exon file (pkl file)>" + " <path/output file (txt file)>"
# longest orf file: output file of script getORFs.py
# exon file: output file of script getExondata.py

if len(sys.argv) != 4:
  print(len(sys.argv))
  print(usage)
  sys.exit()

log = False # set to True for log output

def read_longest_orf_file(orf_file):
# orf_file: expected: lines consisting of transcript IDs, ORF start, ORF end
# returns dictionary containing transcript IDs as keys and tuples of ORF boundaries as values
  transcript_data = {}
  orfFile = open(orf_file)
  for line in orfFile:
    # caution: file contains empty lines
    if line not in ['\n','\r\n']:
      try:
        tx_id,start,end = line.strip().split(",")
        transcript_data[tx_id] = (int(start),int(end))
        #data.append(transcript_data)
      except:
        print('read_orf_file1: Error in file reading')
  return transcript_data


def get_raw_exon_coordinates(ref_key, transcriptModels):
# ref_key: transcript id of reference transcript
# transcriptModels: dictionary containing transcript IDs as keys and lists of lists of exon and intron lengths as values
# Returns list of intron containing exon coordinates for handed over reference transcript (exon coordinates starting from 1)
  exon_lengths = transcriptModels[ref_key][0]
  intron_lengths = transcriptModels[ref_key][1]
  start_coord = 0
  end_coord = 0
  raw_exon_coordinates = []
  i = 0
  for exon_length in exon_lengths:
    raw_exon_coordinates.append((start_coord+1, start_coord+exon_length))
    start_coord += exon_length
    if i < len(intron_lengths):
      start_coord += intron_lengths[i]
    i += 1
  return raw_exon_coordinates


def check_for_cryptic_intron(raw_exon_coordinates_ref, raw_exon_coordinates_iso, i):
# raw_exon_coordinates_ref: list of intron containing exon coordinates for certain reference transcript
# raw_exon_coordinates_iso: list of intron containing exon coordinates for isoforms to certain reference transcript
# i = index of the exon list
# assumed as given: 1) end coordinates of currently examined exons of reference and isoform are unequal; 2) end coordinate of reference > end coordinate of isoform
# check whether end coordinate of current exon of reference coincides with end coordinate of following exon of isoform
# if yes, 
# return raw end and start coordinates of intron-preceding and intron-following exons, respectively, otherwise empty tuple
  global log
  ci_coordinates = ()
  try: 
    if len(raw_exon_coordinates_iso) <= (i+1):
      # following comparison not possible
      if log: print('Failure in check_for_cryptic_intron(): no following exon in isoform transcript existent')
      return ci_coordinates
    if raw_exon_coordinates_ref[i][1] == raw_exon_coordinates_iso[i+1][1]:
      if log: print('End of raw exon ', i, ' of reference = end of raw exon ', i+1, ' of isoform: ', raw_exon_coordinates_ref[i][1])
      ci_coordinates = (raw_exon_coordinates_iso[i][1], raw_exon_coordinates_iso[i+1][0])
      return ci_coordinates
    else:
      if log: print('no cryptic intron')
  except IndexError:
    print('IndexError in check_for_cryptic_intron()')
  except:
    print('Error in check_for_cryptic_intron()')
  return ci_coordinates


def check_for_intron_retention(raw_exon_coordinates_ref, raw_exon_coordinates_iso, i):
# raw_exon_coordinates_ref: list of intron containing exon coordinates for certain reference transcript
# raw_exon_coordinates_iso: list of intron containing exon coordinates for isoforms to certain reference transcript
# i: index of the exon list
# assumed as given: 1) end coordinates of currently examined exons of reference and isoform are unequal; 2) end coordinate of reference < end coordinate of isoform
# check whether end coordinate of following exon of reference coincides with end coordinate of current exon of isoform
# if yes, 
# return True, otherwise False
  global log
  try:
    if len(raw_exon_coordinates_ref) <= i+1:
      # no following exon in reference transcript existent, so that following comparison not possible
      if log: print('Failure in check_for_intron_retention(): no check for ir possible')
      return False
    if raw_exon_coordinates_ref[i+1][1] == raw_exon_coordinates_iso[i][1]:
      if log: print('End of raw exon ', i+1, ' of reference = end of raw exon ', i, ' of isoform: ', raw_exon_coordinates_ref[i+1][1])
      return True
    else:
      if log: print('no intron retention')
  except IndexError:
    print('IndexError in check_for_intron_retention()')
  except:
    print('Error in check_for_intron_retention()')
  return False


def check_for_alt_acceptor(raw_exon_coordinates_ref, raw_exon_coordinates_iso, i):
# raw_exon_coordinates_ref: list of intron containing exon coordinates for certain reference transcript
# raw_exon_coordinates_iso: list of intron containing exon coordinates for isoforms to certain reference transcript
# i: index of the exon list
# assumed as given: start coordinates of currently examined exons of reference and isoform are unequal
# check whether end coordinate of current exon of reference coincides with end coordinate of current exon of isoform
# returns start coordinate of modified exon (in relation to reference)
  global log
  raw_start_coordinate_modified_exon = ()
  try:
    if raw_exon_coordinates_ref[i][1] == raw_exon_coordinates_iso[i][1]:
      if log: print('End of raw exon ', i, ' of reference = end of raw exon ', i, ' of isoform: ', raw_exon_coordinates_ref[i][1])
      raw_start_coordinate_modified_exon = (raw_exon_coordinates_iso[i][0])
      return raw_start_coordinate_modified_exon
    else:
      if log: print('no alternative acceptor')
  except IndexError:
    print('IndexError in check_for_alt_acceptor()')
  except:
    print('Error in check_for_alt_acceptor()')
  return raw_start_coordinate_modified_exon


def check_for_alt_donor(raw_exon_coordinates_ref, raw_exon_coordinates_iso, i):
# raw_exon_coordinates_ref: list of intron containing exon coordinates for certain reference transcript
# raw_exon_coordinates_iso: list of intron containing exon coordinates for isoforms to certain reference transcript
# i: index of the exon list
# assumed as given: end coordinates of currently examined exons of reference and isoform are unequal
# check whether start coordinate of following exon of reference coincides with start coordinate of following exon of isoform
# returns end coordinate of modified exon (in relation to reference)
  global log
  raw_end_coordinate_modified_exon = ()
  try:
    if len(raw_exon_coordinates_ref) <= i+1 or len(raw_exon_coordinates_iso) <= i+1:
      # no following exon in reference or isoform transcript existent, so that following comparison not possible
      if log: print('Failure in check_for_alt_donor(): no check for altd possible')
      return raw_end_coordinate_modified_exon
    if raw_exon_coordinates_ref[i+1][0] == raw_exon_coordinates_iso[i+1][0]:
      if log: print('Start of raw exon of reference ', i+1, ' = start of raw exon of isoform ', i+1, ' : ', raw_exon_coordinates_ref[i+1][0])
      raw_end_coordinate_modified_exon = (raw_exon_coordinates_iso[i][1])
      return raw_end_coordinate_modified_exon
  except IndexError:
    print('IndexError in check_for_alt_donor()')
  except:
    print('Error in check_for_alt_donor()')
  return raw_end_coordinate_modified_exon


def check_for_exon_retention(raw_exon_coordinates_ref, raw_exon_coordinates_iso, i):
# raw_exon_coordinates_ref: list of intron containing exon coordinates for certain reference transcript
# raw_exon_coordinates_iso: list of intron containing exon coordinates for isoforms to certain reference transcript
# i: index of the exon list
# assumed as given: 1) start coordinates of currently examined exons of reference and isoform are unequal; 2) start coordinate of reference > start coordinate of isoform
# check whether start coordinate of current exon of reference coincides with start coordinate of the following exon of isoform; in this case 
# return True, otherwise False
  global log
  try:
    if len(raw_exon_coordinates_iso) <= i+1:
      # no following exon in isoform transcript existent, so that following comparison not possible
      if log: print('Failure in check_for_exon_retention(): no check for exonin possible')
      return False
    if raw_exon_coordinates_ref[i][0] == raw_exon_coordinates_iso[i+1][0]:
      if log: print('Start of raw exon of reference ', i, ' = start of raw exon of isoform ', i+1, ' : ', raw_exon_coordinates_ref[i][0])
      return True
  except IndexError:
    print('IndexError in check_for_exon_retention()')
  except:
    print('Error in check_for_exon_retention()')
  return False


def check_for_exon_skipping(raw_exon_coordinates_ref, raw_exon_coordinates_iso, i):
# raw_exon_coordinates_ref: list of intron containing exon coordinates for certain reference transcript
# raw_exon_coordinates_iso: list of intron containing exon coordinates for isoforms to certain reference transcript
# i: index of the exon list
# assumed as given: 1) start coordinates of currently examined exons of reference and isoform are unequal; 2) start coordinate of reference < start coordinate of isoform
# check whether start coordinate of the following exon of reference coincides with start coordinate of current exon of isoform
# if yes, 
# return start and end coordinates of skipped exon + end coordinate of previous exon + start coordinate of next exon, otherwise empty tuple
  global log
  exonskip_coordinates = ()
  try:
    if len(raw_exon_coordinates_ref) <= i+1:
      # no following exon in reference transcript existent, so that following comparison not possible
      if log: print('Failure in check_for_exon_skipping(): no following exon in reference transcript existent')
      return exonskip_coordinates
    if raw_exon_coordinates_ref[i+1][0] == raw_exon_coordinates_iso[i][0]:
      if log: print('Start of raw exon of reference ', i+1, ' = start of raw exon of isoform ', i, ' : ', raw_exon_coordinates_ref[i+1][0])
      exonskip_coordinates = (raw_exon_coordinates_ref[i][0], raw_exon_coordinates_ref[i][1], raw_exon_coordinates_ref[i-1][1], raw_exon_coordinates_ref[i+1][0])
      return exonskip_coordinates
  except IndexError:
    print('IndexError in check_for_exon_skipping()')
  except:
    print('Error in check_for_exon_skipping()')
  return exonskip_coordinates


def check_for_mutually_exclusive_exons(raw_exon_coordinates_ref, raw_exon_coordinates_iso, i):
# raw_exon_coordinates_ref: list of intron containing exon coordinates for certain reference transcript
# raw_exon_coordinates_iso: list of intron containing exon coordinates for isoforms to certain reference transcript
# i: index of the exon list
# assumed as given: 1) start coordinates of currently examined exons of reference and isoform are unequal
# check whether the current exon of isoform is located in front of the current exon of reference or behind it and 
# check whether start coordinates of the following exons of reference and isoform coincide
# if yes,
# return start and end coordinates of skipped exon + end coordinate of previous exon + start coordinate of next exon, otherwise empty tuple
  global log
  mxe_coordinates = ()

  try:
    if len(raw_exon_coordinates_ref) <= i+1 or len(raw_exon_coordinates_iso) <= i+1 or i == 0:
      # no following exon in reference transcript or isoform existent, so that following comparison not possible
      # if i = 0, no previous exon accessible
      if log: print('Failure in check_for_mutually_exclusive_exons(): no following exon in reference transcript existent')
      return mxe_coordinates
    if (raw_exon_coordinates_iso[i][1] < raw_exon_coordinates_ref[i][0] or raw_exon_coordinates_iso[i][0] > raw_exon_coordinates_ref[i][1]) and (raw_exon_coordinates_ref[i+1][0] == raw_exon_coordinates_iso[i+1][0]):
      if log: print('Start of raw exon of reference ', i+1, ' = start of raw exon of isoform ', i+1, ' : ', raw_exon_coordinates_ref[i+1][0])
      mxe_coordinates = (raw_exon_coordinates_ref[i][0], raw_exon_coordinates_ref[i][1], raw_exon_coordinates_ref[i-1][1], raw_exon_coordinates_ref[i+1][0])
      return mxe_coordinates
  except IndexError:
    print('IndexError in check_for_mutually_exclusive_exons()')
  except:
    print('Error in check_for_mutually_exclusive_exons()')
  return mxe_coordinates


def detect_as_event(transcriptModels, ref_key, iso_key, ir_list, alta_list, altd_list, er_list, ci_list, unresolved_list, es_list):
# transcriptModels: dictionary containing exon data for transcripts originating from gtf file (transcript ids as keys and a list of exon related tupels as value)
# ref_key: key of the reference transcript
# iso_key: key of the transcript isoform for comparison with the reference
# ir_list: list of IR events, containing tuples of transcript_id + number of retained intron - to be filled here
# alta_list: list of AltA events, containing tuples of transcript_id + number of affected exon - to be filled here
# altd_list: list of AltD events, containing tuples of transcript_id + number of affected exon - to be filled here
# er_list: list of ER (exon retention) events, containing tuples of transcript_id + number of affected exon - to be filled here
# ci_list: list of cryptic intron events, containing tuples of transcript_id + number of affected exon - to be filled here
# unresolved_list: list of unresolved AS events, containing tuples of transcript_id + number of affected exon - to be filled here
# es_list: list of ES events, containing tuples of transcript_id + number of affected exon - to be filled here
# compare exon boundaries of reference transcript and isoform in parallel; assign first deviation (from 5' site) to an AS event
# returns tuple containing first deviation coordinate and AS event code (ir = 1, alta = 2, altd = 3, er = 4, ci = 5, es = 6, mxe = 7) or 
#   tuple (0, 0) for unresolved AS event
  global log
  as_event = ()
  deviating_exon_coordinate = 0
  as_event_code = 0
  i = 0
  ci_t = ()
  altd_t = ()
  alta_t = ()
  es_t = ()
  mxe_t = ()
  raw_exon_coordinates_ref = get_raw_exon_coordinates(ref_key, transcriptModels)
  raw_exon_coordinates_iso = get_raw_exon_coordinates(iso_key, transcriptModels)
  if log: print(f"detect_as_event: raw_exon_coordinates_ref = {raw_exon_coordinates_ref}")
  if log: print(f"detect_as_event: raw_exon_coordinates_iso = {raw_exon_coordinates_iso}")
  for exon_ref in raw_exon_coordinates_ref:
    if log: print(f"detect_as_event: exon_ref = {exon_ref}")
    # for loop picks up exons of reference transcript
    try:
      if exon_ref[0] != raw_exon_coordinates_iso[i][0]:
        if log: print('Deviation in raw start coordinate of exon ', i, ': ', exon_ref[0], '(reference) != ', raw_exon_coordinates_iso[i][0], ' (isoform)')
        deviating_exon_coordinate = raw_exon_coordinates_iso[i][0]
        # mutually exclusive exon
        mxe_t = check_for_mutually_exclusive_exons(raw_exon_coordinates_ref, raw_exon_coordinates_iso, i)
        if mxe_t != ():
          if log: print('mutually_exclusive_exons')
          mxe_list.append((iso_key, mxe_t))
          as_event_code = 7
          break
        # alternative acceptor: isoform exon can be either shorter (starts later) or longer (starts earlier) than reference exon
        alta_t = check_for_alt_acceptor(raw_exon_coordinates_ref, raw_exon_coordinates_iso, i)
        if alta_t != ():
          if log: print('alternative acceptor')
          # save transcript id of isoform and coordinate of modified exon:
          alta_list.append((iso_key, alta_t))
          as_event_code = 2
          # treat only first coordinate deviation
          break
        # exon retention: reference exon starts later than (inserted) exon of isoform
        if exon_ref[0] > raw_exon_coordinates_iso[i][0]:
          if check_for_exon_retention(raw_exon_coordinates_ref, raw_exon_coordinates_iso, i) == True:
            if log: print('exon retention')
            # save transcript id of isoform and number of previous intron (where additional exon of isoform is inserted):
            er_list.append((iso_key, i-1))
            as_event_code = 4
            break
        # exon skipping: reference exon starts earlier than exon of isoform, which follows on not more present exon
        if exon_ref[0] < raw_exon_coordinates_iso[i][0]:
            es_t = check_for_exon_skipping(raw_exon_coordinates_ref, raw_exon_coordinates_iso, i)
            if es_t != ():
              if log: print('exon skipping')
              # save transcript id of isoform and coordinates of skipped exon of reference and coordinates of adjacent exons:
              es_list.append((iso_key, es_t))
              as_event_code = 6
              break
        if log: print('unresolved as event')
        unresolved_list.append((iso_key, i))
        # treat only first coordinate deviation
        break
      if exon_ref[1] != raw_exon_coordinates_iso[i][1]:
        if log: print('Deviation in raw end coordinate of raw exon ', i, ': ', exon_ref[1], '(reference) != ', raw_exon_coordinates_iso[i][1], ' (isoform)')
        deviating_exon_coordinate = raw_exon_coordinates_iso[i][1]
        # alternative donor: isoform exon can be either shorter (ends earlier) or longer (ends later) than reference exon
        altd_t = check_for_alt_donor(raw_exon_coordinates_ref, raw_exon_coordinates_iso, i)
        if altd_t != ():
          if log: print('alternative donor')
          # save transcript id of isoform and coordinate of modified exon:
          altd_list.append((iso_key, altd_t))
          as_event_code = 3
          break
        # cryptic intron: reference exon longer than intron-preceding exon of isoform
        if exon_ref[1] > raw_exon_coordinates_iso[i][1]:
          ci_t = check_for_cryptic_intron(raw_exon_coordinates_ref, raw_exon_coordinates_iso, i)
          if ci_t != ():
            if log: print('cryptic intron')
            # save transcript id of isoform and coordinates of adjacent exons:
            ci_list.append((iso_key, ci_t))
            as_event_code = 5
            break
        # intron retention: reference exon shorter than isoform exon (which contains former intron)
        if exon_ref[1] < raw_exon_coordinates_iso[i][1]:
          if check_for_intron_retention(raw_exon_coordinates_ref, raw_exon_coordinates_iso, i) == True:
            if log: print('intron retention')
            # save transcript id of isoform and number of retained intron:
            ir_list.append((iso_key, i))
            as_event_code = 1
            break
        if log: print('unresolved AS event')
        unresolved_list.append((iso_key, i))
        break
    except IndexError:
      print('IndexError in detect_as_event()')
      break
    except:
      print('Error other than IndexError in detect_as_event()')
      break  
    i += 1
  as_event = (deviating_exon_coordinate, as_event_code)
  return as_event


def get_intronless_exon_coordinates(transcript_id, transcriptModels):
# transcript_id: ID of transcript for which exon coordinates are to be calculated
# transcriptModels: dictionary containing transcript IDs as keys and lists of lists of exon and intron lengths as values
# returns dictionary consisting of transcript ID as key and list of exon coordinates as value
  intronless_exon_coordinates = []
  transcriptModels1 = {}
  previous_exon_length = 0
  for exon_length in transcriptModels[transcript_id][0]:
    exon_coordinate = (previous_exon_length+1, previous_exon_length+exon_length)
    intronless_exon_coordinates.append(exon_coordinate)
    previous_exon_length += exon_length
  transcriptModels1[transcript_id] = intronless_exon_coordinates
  return transcriptModels1


def get_intronless_translation_boundary_exons(ref_transcript_id, transcriptModels, orf_data):
# ref_transcript_id: ID of transcript for which translation boundaries are to be calculated
# transcriptModels: dictionary containing transcript IDs as keys and lists of lists of exon and intron lengths as values
# orf_data: dictionary containing transcript IDs as keys and tuples of ORF boundaries as values; ORF end after stop codon -> 
#   therefore 3 nucleotides to subtract
# Returns tuple containing handed over reference transcript id, number of start exon of translation (0 based), 
#   number of stop exon of translation; for each exon that is not found, -1 is returned
  global log
  transcriptModel = {}
  translation_data = ()
  start_exon_number = -1
  stop_exon_number = -1
  i = 0
  try:
    transcriptModel = get_intronless_exon_coordinates(ref_transcript_id, transcriptModels)
    if log: print(f"transcriptModel = {transcriptModel}")
    for exon in transcriptModel[ref_transcript_id]:
      if (exon[0] <= orf_data[ref_transcript_id][0]) and (orf_data[ref_transcript_id][0] <= exon[1]):
        start_exon_number = i
      if (exon[0] <= orf_data[ref_transcript_id][1]-3) and (orf_data[ref_transcript_id][1]-3 <= exon[1]):
        stop_exon_number = i
      i += 1
    translation_data = (ref_transcript_id, start_exon_number, stop_exon_number)
  except IndexError:
    print('IndexError in get_intronless_translation_boundary_exons()')
  except:
    print('Error in get_intronless_translation_boundary_exons()')
  return translation_data


def get_sum_intron_lengths(ref_transcript_id, transcriptModels, intronless_translation_boundary_exon):
# ref_transcript_id: ID of transcript for which translation coordinates are to be calculated
# transcriptModels: dictionary containing transcript IDs as keys and lists of lists of exon and intron lengths as values
# intronless_translation_boundary_exon: number of start or stop exon of translation (0 based)
  sum_intron_lengths = 0
  i = 0
  while i < intronless_translation_boundary_exon:
    sum_intron_lengths += transcriptModels[ref_transcript_id][1][i]
    i += 1
  return sum_intron_lengths


def get_raw_translation_coordinates(ref_transcript_id, transcriptModels, orf_data):
# ref_transcript_id: ID of transcript for which translation coordinates are to be calculated
# transcriptModels: dictionary containing transcript IDs as keys and lists of lists of exon and intron lengths as values
# orf_data: dictionary containing transcript IDs as keys and tuples of ORF boundaries as values
# Returns tuple containing raw start and end coordinate of translation for the transcript with handed over reference transcript ID
  intronless_translation_boundary_exons = ()
  raw_translation_start_coordinate = -1
  raw_translation_end_coordinate = -1
  intronless_translation_boundary_exons = get_intronless_translation_boundary_exons(ref_transcript_id, transcriptModels, orf_data)
  raw_translation_start_coordinate = orf_data[ref_transcript_id][0] + get_sum_intron_lengths(ref_transcript_id, transcriptModels, intronless_translation_boundary_exons[1])
  raw_translation_end_coordinate = orf_data[ref_transcript_id][1] - 3 + get_sum_intron_lengths(ref_transcript_id, transcriptModels, intronless_translation_boundary_exons[2])
  return (raw_translation_start_coordinate, raw_translation_end_coordinate)


def get_ci_localization(transcript_id, ci_list, raw_translation_coordinates):
# transcript_id: transcript id of isoform to be examined
# ci_list: list of cryptic intron events, containing tuples of transcript id + tuple of exon coordinates needed for localization
# raw_translation_coordinates: Tuple containing raw start and end coordinate of translation
# determine whether ci happened in 5'UTR, CDS, 3'UTR, 5'UTR+CDS, CDS+3'UTR, or 5'UTR+CDS+3'UTR
# Returns 0 for unsuccessful localization, 1 for 5'UTR, 2 for CDS, 3 for 3'UTR, 4 for 5'UTR+CDS, 5 for CDS+3'UTR, 6 for 5'UTR+CDS+3'UTR
  global log
  localization = 0
  raw_isoform_end_coordinate_exon_ahead_inserted_intron = 0
  raw_isoform_start_coordinate_exon_after_inserted_intron = 0

  # get raw start and end coordinate of adjacent exons:
  # - find corresponding list entry in ci_list: 
  for entry in ci_list:
    if entry[0] == transcript_id:
      raw_isoform_end_coordinate_exon_ahead_inserted_intron = entry[1][0]
      raw_isoform_start_coordinate_exon_after_inserted_intron = entry[1][1]
      break

  # comparisons
  if raw_translation_coordinates[1] <= raw_isoform_end_coordinate_exon_ahead_inserted_intron:
    # event in 3'UTR
    if log: print(f"* Cryptic intron in transcript isoform {transcript_id} happened in 3UTR")
    localization = 3
    return localization
  if raw_translation_coordinates[0] >= raw_isoform_start_coordinate_exon_after_inserted_intron:
    # event in 5'UTR
    if log: print(f"* Cryptic intron in transcript isoform {transcript_id} happened in 5UTR")
    localization = 1
    return localization
  if raw_translation_coordinates[0] <= raw_isoform_end_coordinate_exon_ahead_inserted_intron+1 and raw_translation_coordinates[1] >= raw_isoform_start_coordinate_exon_after_inserted_intron-1:
    # event in CDS
    if log: print(f"* Cryptic intron in transcript isoform {transcript_id} happened in CDS")
    localization = 2
    return localization
  if raw_translation_coordinates[0] > raw_isoform_end_coordinate_exon_ahead_inserted_intron+1 and raw_translation_coordinates[0] < raw_isoform_start_coordinate_exon_after_inserted_intron and raw_translation_coordinates[1] >= raw_isoform_start_coordinate_exon_after_inserted_intron-1:
    # event in 5UTR + CDS
    if log: print(f"* Cryptic intron in transcript isoform {transcript_id} happened in 5UTR and CDS")
    localization = 4
    return localization
  if raw_translation_coordinates[0] <= raw_isoform_end_coordinate_exon_ahead_inserted_intron+1 and raw_translation_coordinates[1] > raw_isoform_end_coordinate_exon_ahead_inserted_intron and raw_translation_coordinates[1] < raw_isoform_start_coordinate_exon_after_inserted_intron-1:
    # event in CDS + 3UTR
    if log: print(f"* Cryptic intron in transcript isoform {transcript_id} happened in CDS and 3UTR")
    localization = 5
    return localization
  if raw_translation_coordinates[0] > raw_isoform_end_coordinate_exon_ahead_inserted_intron+1 and raw_translation_coordinates[1] < raw_isoform_start_coordinate_exon_after_inserted_intron-1:
    # event in 5UTR + CDS + 3UTR
    if log: print(f"* Cryptic intron in transcript isoform {transcript_id} happened in 5UTR and CDS and 3UTR")
    localization = 6
    return localization

  return localization


def get_ir_er_localization(ref_transcript_id, transcript_id, ir_er_list, intronless_translation_boundary_exons):
# ref_transcript_id: transcript id of reference transcript
# transcript_id: transcript id of isoform to be examined
# ir_er_list: list of ir/er events, containing tuples of transcript_id + number of retained intron/exon
# intronless_translation_boundary_exons: tuple containing handed over reference transcript id, start exon of translation, stop exon of translation
# determine whether ir/er happened in 5'-UTR, CDS or 3'-UTR
# Returns 0 for unsuccessful localization, 1 for 5'UTR, 2 for CDS, 3 for 3'UTR
# logic:
# if translation start in or downstream to exon/intron succeeding retained intron/exon, then ir/er in 5'UTR;
# else if translation end upstream to exon/intron succeeding retained intron/exon, then ir/er in 3'UTR;
# else ir/er in CDS
  global log
  localization = 0
  # fish for number of retained intron
  intron_exon_number = -1
  for entry in ir_er_list:
    if entry[0] == transcript_id:
      intron_exon_number = entry[1]
      break
  if log: print('get_ir_er_localization(): intron_exon_number = ', intron_exon_number)
  # find number of exon of reference after retained intron
  succeeding_exon_intron_number = intron_exon_number+1
    
  if intronless_translation_boundary_exons[0] == ref_transcript_id and intron_exon_number != -1:
    if intronless_translation_boundary_exons[1] >= succeeding_exon_intron_number:
      if log: print(f"* IR/ER in transcript isoform {transcript_id} happened in 5UTR")
      localization = 1
    if intronless_translation_boundary_exons[2] < succeeding_exon_intron_number:
      if log: print(f"* IR/ER in transcript isoform {transcript_id} happened in 3UTR")
      localization = 3
    if intronless_translation_boundary_exons[1] < succeeding_exon_intron_number and intronless_translation_boundary_exons[2] >= succeeding_exon_intron_number:
      if log: print(f"* IR/ER in transcript isoform {transcript_id} happened in CDS")
      localization = 2

  return localization


def get_alta_localization(transcript_id, alta_list, raw_translation_coordinates):
# transcript_id: transcript id of isoform to be examined
# alta_list: list of alta events, containing tuples of transcript_id + number of modified exon
# raw_translation_coordinates: Tuple containing raw start and end coordinate of translation
# determine whether alta happened in 5'-UTR, CDS or 3'-UTR
# Returns 0 for unsuccessful localization, 1 for 5'UTR, 2 for CDS, 3 for 3'UTR
# logic:
# import raw start coordinate of modified (either elongated or shortened) exon
# if modified start coordinate at or downstream translation start (ts), event in CDS or in 3'UTR, otherwise (= upstream ts) event in 5'UTR
# in case of at or downstream ts: if modified start coordinate at or upstream translation end (te), event in CDS; 
#   else (= downstream te), event in 3'UTR
  global log
  localization = 0
  raw_start_coordinate_modified_exon = 0
  # get raw start coordinate of modified exon:
  for entry in alta_list:
    if entry[0] == transcript_id:
      raw_start_coordinate_modified_exon = entry[1]
      break
  # comparisons between event coordinate and translation coordinates:
  if raw_start_coordinate_modified_exon >= raw_translation_coordinates[0]:
    # event in CDS or 3UTR
    if raw_start_coordinate_modified_exon <= raw_translation_coordinates[1]:
      # event in CDS
      if log: print(f"* AltA in transcript isoform {transcript_id} happened in CDS - AS event at {raw_start_coordinate_modified_exon}")
      localization = 2
    else:
      if log: print(f"* AltA in transcript isoform {transcript_id} happened in 3UTR - AS event at {raw_start_coordinate_modified_exon}")
      localization = 3
  else:
    # event in 5UTR
    if log: print(f"* AltA in transcript isoform {transcript_id} happened in 5UTR - AS event at {raw_start_coordinate_modified_exon}")
    localization = 1

  return localization


def get_altd_localization(transcript_id, altd_list, raw_translation_coordinates):
# transcript_id: transcript id of isoform to be examined
# altd_list: list of altd events, containing Tupels of transcript_id + number of modified exon
# raw_translation_coordinates: Tupel containing raw start and end coordinate of translation
# determine whether altd happened in 5'-UTR, CDS or 3'-UTR
# Returns 0 for unsuccessful localization, 1 for 5'UTR, 2 for CDS, 3 for 3'UTR
# logic:
# import raw end coordinate of modified (either elongated or shortened) exon
# if modified end coordinate at or downstream translation start (ts), event in CDS or in 3'UTR, otherwise (= upstream ts) event in 5'UTR
# in case of at or downstream ts: if modified end coordinate at or upstream translation end (te), event in CDS; 
#   else (= downstream te), event in 3'UTR
  global log
  localization = 0
  raw_end_coordinate_modified_exon = 0
  # get raw end coordinate of modified exon:
  for entry in altd_list:
    if entry[0] == transcript_id:
      raw_end_coordinate_modified_exon = entry[1]
      break
  # comparisons between event coordinate and translation coordinates:
  if raw_end_coordinate_modified_exon >= raw_translation_coordinates[0]:
    # event in CDS or 3UTR
    if raw_end_coordinate_modified_exon <= raw_translation_coordinates[1]:
      # event in CDS
      if log: print(f"* AltD in transcript isoform {transcript_id} happened in CDS - AS event at {raw_end_coordinate_modified_exon}")
      localization = 2
    else:
      if log: print(f"* AltD in transcript isoform {transcript_id} happened in 3UTR - AS event at {raw_end_coordinate_modified_exon}")
      localization = 3
  else:
    # event in 5UTR
    if log: print(f"* AltD in transcript isoform {transcript_id} happened in 5UTR - AS event at {raw_end_coordinate_modified_exon}")
    localization = 1

  return localization


def get_exonskip_localization(transcript_id, es_list, raw_translation_coordinates):
# transcript_id: transcript id of isoform to be examined
# es_list_list: list of es_list events, containing tuples of transcript id + tuple of coordinates needed for localization
# raw_translation_coordinates: Tuple containing raw start and end coordinate of translation
# determine whether exon skipping happened in 5'UTR, CDS, 3'UTR, 5'UTR+CDS, CDS+3'UTR, or 5'UTR+CDS+3'UTR
# Returns 0 for unsuccessful localization, 1 for 5'UTR, 2 for CDS, 3 for 3'UTR, 4 for 5'UTR+CDS, 5 for CDS+3'UTR, 6 for 5'UTR+CDS+3'UTR
  global log
  localization = 0
  raw_end_coordinate_exon_ahead_skipped_one = 0
  raw_start_coordinate_exon_after_skipped_one = 0
  raw_start_coordinate_skipped_exon = 0
  raw_end_coordinate_skipped_exon = 0

  # get raw start and end coordinate of skipped exon:
  # - find corresponding list entry in es_list: 
  for entry in es_list:
    if entry[0] == transcript_id:
      raw_start_coordinate_skipped_exon = entry[1][0]
      raw_end_coordinate_skipped_exon = entry[1][1]
      raw_end_coordinate_exon_ahead_skipped_one = entry[1][2]
      raw_start_coordinate_exon_after_skipped_one = entry[1][3]
      break

  # comparisons
  if raw_translation_coordinates[1] <= raw_end_coordinate_exon_ahead_skipped_one:
    # event in 3'UTR
    if log: print(f"* Exon skip in transcript isoform {transcript_id} happened in 3UTR")
    localization = 3
    return localization
  if raw_translation_coordinates[0] >= raw_start_coordinate_exon_after_skipped_one:
    # event in 5'UTR
    if log: print(f"* Exon skip in transcript isoform {transcript_id} happened in 5UTR")
    localization = 1
    return localization
  if raw_translation_coordinates[0] <= raw_start_coordinate_skipped_exon and raw_translation_coordinates[1] >= raw_end_coordinate_skipped_exon:
    # event in CDS
    if log: print(f"* Exon skip in transcript isoform {transcript_id} happened in CDS")
    localization = 2
    return localization
  if raw_translation_coordinates[0] > raw_start_coordinate_skipped_exon and raw_translation_coordinates[0] < raw_start_coordinate_exon_after_skipped_one and raw_translation_coordinates[1] >= raw_end_coordinate_skipped_exon:
    # event in 5'UTR + CDS
    if log: print(f"* Exon skip in transcript isoform {transcript_id} happened in 5UTR and CDS")
    localization = 4
    return localization
  if raw_translation_coordinates[0] <= raw_start_coordinate_skipped_exon and raw_translation_coordinates[1] < raw_end_coordinate_skipped_exon and raw_translation_coordinates[1] > raw_end_coordinate_exon_ahead_skipped_one:
    # event in CDS + 3'UTR
    if log: print(f"* Exon skip in transcript isoform {transcript_id} happened in CDS and 3UTR")
    localization = 5
    return localization
  if raw_translation_coordinates[0] > raw_start_coordinate_skipped_exon and raw_translation_coordinates[1] < raw_end_coordinate_skipped_exon:
    # event in 5'UTR + CDS + 3'UTR
    if log: print(f"* Exon skip in transcript isoform {transcript_id} happened in 5UTR and CDS and 3UTR")
    localization = 6
    return localization

  return localization


def get_mxe_localization(transcript_id, mxe_list, raw_translation_coordinates):
# transcript_id: transcript id of isoform to be examined
# mxe_list_list: list of es_list events, containing tuples of transcript id + tuple of coordinates needed for localization
# raw_translation_coordinates: Tuple containing raw start and end coordinate of translation
# determine whether mxe happened in 5'UTR, CDS, 3'UTR
# Returns 0 for unsuccessful localization, 1 for 5'UTR, 2 for CDS, 3 for 3'UTR
  global log
  localization = 0
  raw_end_coordinate_exon_ahead_mutually_excluded_ones = 0
  raw_start_coordinate_exon_after_mutually_excluded_ones = 0
  raw_start_coordinate_skipped_exon = 0 # skipped in isoform
  raw_end_coordinate_skipped_exon = 0

  # get raw start and end coordinate of skipped exon:
  # - find corresponding list entry in es_list: 
  for entry in mxe_list:
    if entry[0] == transcript_id:
      #raw_start_coordinate_skipped_exon = entry[1][0]
      #raw_end_coordinate_skipped_exon = entry[1][1]
      raw_end_coordinate_exon_ahead_mutually_excluded_ones = entry[1][2]
      raw_start_coordinate_exon_after_mutually_excluded_ones = entry[1][3]
      break
  
  # comparisons
  if raw_translation_coordinates[1] <= raw_end_coordinate_exon_ahead_mutually_excluded_ones:
    # event in 3'UTR
    if log: print(f"* Mutually excluded exon in transcript isoform {transcript_id} happened in 3UTR")
    localization = 3
    return localization
  if raw_translation_coordinates[0] >= raw_start_coordinate_exon_after_mutually_excluded_ones:
    # event in 5'UTR
    if log: print(f"* Mutually excluded exon in transcript isoform {transcript_id} happened in 5UTR")
    localization = 1
    return localization
  if raw_translation_coordinates[0] <= raw_end_coordinate_exon_ahead_mutually_excluded_ones and raw_translation_coordinates[1] >= raw_start_coordinate_exon_after_mutually_excluded_ones:
    # event in CDS
    if log: print(f"* Mutually excluded exon in transcript isoform {transcript_id} happened in CDS")
    localization = 2
    return localization
  if raw_translation_coordinates[0] > raw_end_coordinate_exon_ahead_mutually_excluded_ones and raw_translation_coordinates[0] < raw_start_coordinate_exon_after_mutually_excluded_ones and raw_translation_coordinates[1] >= raw_start_coordinate_exon_after_mutually_excluded_ones:
    # event in 5'UTR + CDS
    if log: print(f"* Mutually excluded exon in transcript isoform {transcript_id} happened in 5UTR and CDS")
    localization = 4
    return localization
  if raw_translation_coordinates[0] <= raw_end_coordinate_exon_ahead_mutually_excluded_ones and raw_end_coordinate_exon_ahead_mutually_excluded_ones < raw_translation_coordinates[1] and raw_translation_coordinates[1] < raw_start_coordinate_exon_after_mutually_excluded_ones:
    # event in CDS + 3'UTR
    if log: print(f"* Mutually excluded exon in transcript isoform {transcript_id} happened in CDS and 3UTR")
    localization = 5
    return localization

  return localization


def detect_and_locate_as_events(ref_tx_for_gene, transcriptModels, orf_data, iso_txs_for_gene, ir_list, alta_list, altd_list, er_list, ci_list, es_list, mxe_list, unresolved_list, as_event_list, as_duplicates):
# ref_tx_for_gene: tuple containing ID of reference transcript for certain gene and associated list of lists of exon and intron lengths as value
# transcriptModels: dictionary containing transcript IDs as keys and lists of lists of exon and intron lengths as values
# orf_data: dictionary containing transcript IDs as keys and tuples of ORF boundaries as values
# iso_txs_for_gene: dictionary containing IDs of transcript isoforms for certain gene and associated lists of lists of exon and intron lengths as values
# ir_list: list of IR events, containing tuples of transcript_id + number of retained intron - to be filled by subroutine detect_as_event()
# alta_list: list of AltA events, containing tuples of transcript_id + number of affected exon - to be filled by subroutine detect_as_event()
# altd_list: list of AltD events, containing tuples of transcript_id + number of affected exon - to be filled by subroutine detect_as_event()
# er_list: list of ER events, containing tuples of transcript_id + number of affected exon - to be filled by subroutine detect_as_event()
# ci_list: list of CI events, containing tuples of transcript_id + number of affected exon - to be filled by subroutine detect_as_event()
# es_list: list of ES events, containing tuples of transcript_id + number of affected exon - to be filled by subroutine detect_as_event()
# mxe_list: list of MXE events, containing tuples of transcript_id + number of affected exon - to be filled by subroutine detect_as_event()
# unresolved_list: list of unresolved AS events, containing tuples of transcript_id + number of affected exon - to be filled by subroutine detect_as_event()
# as_event_list: serves to detect duplicate AS events - to be filled here
# as_duplicates: serves counting the different AS event duplicates - to be filled here
# returns False, if for reference transcript the determination of the exons that define the translation was not successful (means: concerned 
#   gene could not be examined), otherwise True
  global log
  checked_for_as_events = False
  if log: print(f"detect_and_locate_as_events()")
  if ref_tx_for_gene != ():
    intronless_translation_boundary_exons = get_intronless_translation_boundary_exons(ref_tx_for_gene[0], transcriptModels, orf_data)
    #print(f"intronless_translation_boundary_exons = {intronless_translation_boundary_exons}")
    if intronless_translation_boundary_exons[1] == -1 or intronless_translation_boundary_exons[2] == -1:
      # one or both boundary exons could not be determined
      if log: print(intronless_translation_boundary_exons[0], ': lacking exon data: no AS event investigation')
      return checked_for_as_events
    if log: print(f"{intronless_translation_boundary_exons[0]}: intronless translation coordinates: ({intronless_translation_boundary_exons[1]}, {intronless_translation_boundary_exons[2]})")
    raw_translation_coordinates = get_raw_translation_coordinates(ref_tx_for_gene[0], transcriptModels, orf_data)
    if log: print(f"{ref_tx_for_gene[0]}: raw translation coordinates: ({raw_translation_coordinates[0]}, {raw_translation_coordinates[1]})")
    # check transcript isoforms:
    for key1, value1 in iso_txs_for_gene.items():
      if log: print(f"{key1} is further isoform - comparison of exon boundaries required")
      as_event = detect_as_event(transcriptModels, ref_tx_for_gene[0], key1, ir_list, alta_list, altd_list, er_list, ci_list, unresolved_list, es_list)
      if log: print(f"as_event = {as_event}")
      # only count AS event, if it is new for the current gene (= deviating coordinate is different from previous deviating coordinates for a given as event type)
      # explains a discrepancy between a counted event and corresponding list length, which reflects all detected events
      if as_event[0] != 0:
        # resolved AS event
        # only process event if it's not the same event at the same position
        if as_event not in as_event_list:
          as_event_list.append(as_event)
          if log: print(f"raw exon coordinate {as_event[0]} of transcript {key1} deviates")
          # case differentiation according to AS event type:
          if as_event[1] == 1:
            # ir event
            localization = get_ir_er_localization(ref_tx_for_gene[0], key1, ir_list, intronless_translation_boundary_exons)
            if localization == 1:
              ir_events['5UTR'] += 1
            if localization == 2:
              ir_events['CDS'] += 1
            if localization == 3:
              ir_events['3UTR'] += 1
            if localization == 0:
              ir_events['not localized'] += 1

          if as_event[1] == 2:
            # alta event
            localization = get_alta_localization(key1, alta_list, raw_translation_coordinates)
            if localization == 1:
              alta_events['5UTR'] += 1
            if localization == 2:
              alta_events['CDS'] += 1
            if localization == 3:
              alta_events['3UTR'] += 1
            if localization == 0:
              alta_events['not localized'] += 1

          if as_event[1] == 3:
            # altd event
            localization = get_altd_localization(key1, altd_list, raw_translation_coordinates)
            if localization == 1:
              altd_events['5UTR'] += 1
            if localization == 2:
              altd_events['CDS'] += 1
            if localization == 3:
              altd_events['3UTR'] += 1
            if localization == 0:
              altd_events['not localized'] += 1

          if as_event[1] == 4:
            # er event
            localization = get_ir_er_localization(ref_tx_for_gene[0], key1, er_list, intronless_translation_boundary_exons)
            if localization == 1:
              er_events['5UTR'] += 1
            if localization == 2:
              er_events['CDS'] += 1
            if localization == 3:
              er_events['3UTR'] += 1
            if localization == 0:
              er_events['not localized'] += 1

          if as_event[1] == 5:
            # ci event
            localization = get_ci_localization(key1, ci_list, raw_translation_coordinates)
            if localization == 1:
              ci_events['5UTR'] += 1
            if localization == 2:
              ci_events['CDS'] += 1
            if localization == 3:
              ci_events['3UTR'] += 1
            if localization == 4:
              ci_events['5UTR+CDS'] += 1
            if localization == 5:
              ci_events['CDS+3UTR'] += 1
            if localization == 6:
              ci_events['5UTR+CDS+3UTR'] += 1
            if localization == 0:
              ci_events['not localized'] += 1

          if as_event[1] == 6:
            # es event
            localization = get_exonskip_localization(key1, es_list, raw_translation_coordinates)
            if localization == 1:
              es_events['5UTR'] += 1
            if localization == 2:
              es_events['CDS'] += 1
            if localization == 3:
              es_events['3UTR'] += 1
            if localization == 4:
              es_events['5UTR+CDS'] += 1
            if localization == 5:
              es_events['CDS+3UTR'] += 1
            if localization == 6:
              es_events['5UTR+CDS+3UTR'] += 1
            if localization == 0:
              es_events['not localized'] += 1

          if as_event[1] == 7:
            # mxe event
            localization = get_mxe_localization(key1, mxe_list, raw_translation_coordinates)
            if localization == 1:
              mxe_events['5UTR'] += 1
            if localization == 2:
              mxe_events['CDS'] += 1
            if localization == 3:
              mxe_events['3UTR'] += 1
            if localization == 4:
              mxe_events['5UTR+CDS'] += 1
            if localization == 5:
              mxe_events['CDS+3UTR'] += 1
            if localization == 0:
              mxe_events['not localized'] += 1
          as_event = (0, 0)
        else:
          # as_event in as_event_list - means: AS event already detected for the currently investigated gene
          if log: print('duplicate')
          if as_event[1] == 1:
            # ir event
            as_duplicates[0] += 1
          if as_event[1] == 2:
            # alta event
            as_duplicates[1] += 1
          if as_event[1] == 3:
            # altd event
            as_duplicates[2] += 1
          if as_event[1] == 4:
            # er event
            as_duplicates[3] += 1
          if as_event[1] == 5:
            # ci event
            as_duplicates[4] += 1
          if as_event[1] == 6:
            # es event
            as_duplicates[5] += 1
          if as_event[1] == 7:
            # mxe event
            as_duplicates[6] += 1
  checked_for_as_events = True
  return checked_for_as_events


# import reference transcript data
orf_file = sys.argv[1]
orf_data = read_longest_orf_file(orf_file)

# import exon data of all transcripts
with open(sys.argv[2], 'rb') as fp:
    transcriptModels = pickle.load(fp)
# transcriptModels: dictionary containing transcript IDs as keys and lists of lists of exon and intron lengths as values

ref_txs = []
for key in orf_data:
  ref_txs.append(key)
ref_tx = ''
not_investigated_genes = 0
ir_list = []
alta_list = []
altd_list = []
er_list = []
ci_list = []
es_list = []
mxe_list = []
unresolved_list = []
as_event = ()
# examine same AS events for a gene only once; therefore the following list:
as_event_list = []
# count duplicates:
as_duplicates = [0,0,0,0,0,0,0]
intronless_exons = []
intronless_translation_boundary_exons = ()
raw_translation_coordinates = ()
ir_events = {'5UTR':0, 'CDS':0, '3UTR':0, 'not localized':0}
alta_events = {'5UTR':0, 'CDS':0, '3UTR':0, 'not localized':0}
altd_events = {'5UTR':0, 'CDS':0, '3UTR':0, 'not localized':0}
er_events = {'5UTR':0, 'CDS':0, '3UTR':0, 'not localized':0}
ci_events = {'5UTR':0, 'CDS':0, '3UTR':0, '5UTR+CDS':0, 'CDS+3UTR':0, '5UTR+CDS+3UTR':0, 'not localized':0}
es_events = {'5UTR':0, 'CDS':0, '3UTR':0, '5UTR+CDS':0, 'CDS+3UTR':0, '5UTR+CDS+3UTR':0, 'not localized':0}
mxe_events = {'5UTR':0, 'CDS':0, '3UTR':0, '5UTR+CDS':0, 'CDS+3UTR':0,'not localized':0}
gene_name = ''
ref_tx_for_gene = ()
iso_txs_for_gene = {}
tx_counter = 0
examined_genes = 0

for key, value in transcriptModels.items():
  tx_counter += 1
  if gene_name != key[0:9] and gene_name != '':
    # new gene, but not first one
    # AS event analysis for previous gene:
    if detect_and_locate_as_events(ref_tx_for_gene, transcriptModels, orf_data, iso_txs_for_gene, ir_list, alta_list, altd_list, er_list, ci_list, es_list, mxe_list, unresolved_list, as_event_list, as_duplicates) == False:
      not_investigated_genes += 1
    # reset data structures for the new gene
    ref_tx_for_gene = ()
    iso_txs_for_gene = {}
    # set new gene name
    gene_name = key[0:9]
    examined_genes += 1
    # reset as_event_list, because each gene is checked separately for event replicates
    as_event_list = []
  elif gene_name == '':
    # first gene: only set gene name
    gene_name = key[0:9]
  # collect reference transcript and isoform data in data structures
  # when gene switch, do first AS event analysis for the previous gene with the filled data structures for that gene (see above)
  if key in ref_txs:
    # reference transcript 
    ref_tx_for_gene = (key, value)
    if log: print(f"{key}: reference transcript")
  else:
    # transcript isoform
    iso_txs_for_gene[key] = value
  if tx_counter == len(transcriptModels):
    # last transcript of last gene -> no more switch to a new gene will occur
    # therefore trigger here AS event analysis one more time
    if detect_and_locate_as_events(ref_tx_for_gene, transcriptModels, orf_data, iso_txs_for_gene, ir_list, alta_list, altd_list, er_list, ci_list, es_list, mxe_list, unresolved_list, as_event_list, as_duplicates) == False:
      not_investigated_genes += 1
    examined_genes += 1


with open(sys.argv[3],"w") as out_writer:
  out_writer.write(f"Event\t5'UTR\tCDS\t3'UTR\t5'UTR+CDS\tCDS+3'UTR\t5'UTR+CDS+3'UTR\tnot localized\tDuplicates\tSum")
  out_writer.write(f"\nIR\t{ir_events['5UTR']}\t{ir_events['CDS']}\t{ir_events['3UTR']}\t-\t-\t-\t{ir_events['not localized']}\t{as_duplicates[0]}\t{len(ir_list)}")
  out_writer.write(f"\nAltA\t{alta_events['5UTR']}\t{alta_events['CDS']}\t{alta_events['3UTR']}\t-\t-\t-\t{alta_events['not localized']}\t{as_duplicates[1]}\t{len(alta_list)}")
  out_writer.write(f"\nAltD\t{altd_events['5UTR']}\t{altd_events['CDS']}\t{altd_events['3UTR']}\t-\t-\t-\t{altd_events['not localized']}\t{as_duplicates[2]}\t{len(altd_list)}")
  out_writer.write(f"\nER\t{er_events['5UTR']}\t{er_events['CDS']}\t{er_events['3UTR']}\t-\t-\t-\t{er_events['not localized']}\t{as_duplicates[3]}\t{len(er_list)}")
  out_writer.write(f"\nES\t{es_events['5UTR']}\t{es_events['CDS']}\t{es_events['3UTR']}\t{es_events['5UTR+CDS']}\t{es_events['CDS+3UTR']}\t{es_events['5UTR+CDS+3UTR']}\t{es_events['not localized']}\t{as_duplicates[5]}\t{len(es_list)}")
  out_writer.write(f"\nCI\t{ci_events['5UTR']}\t{ci_events['CDS']}\t{ci_events['3UTR']}\t{ci_events['5UTR+CDS']}\t{ci_events['CDS+3UTR']}\t{ci_events['5UTR+CDS+3UTR']}\t{ci_events['not localized']}\t{as_duplicates[4]}\t{len(ci_list)}")
  out_writer.write(f"\nMXE\t{mxe_events['5UTR']}\t{mxe_events['CDS']}\t{mxe_events['3UTR']}\t{mxe_events['5UTR+CDS']}\t{mxe_events['CDS+3UTR']}\t-\t{mxe_events['not localized']}\t{as_duplicates[6]}\t{len(mxe_list)}")

  out_writer.write(f"\nSum\t{ir_events['5UTR'] + alta_events['5UTR'] + altd_events['5UTR'] + er_events['5UTR'] + ci_events['5UTR'] + es_events['5UTR'] + mxe_events['5UTR']}\t")
  out_writer.write(f"{ir_events['CDS'] + alta_events['CDS'] + altd_events['CDS'] + er_events['CDS'] + ci_events['CDS'] + es_events['CDS'] + mxe_events['CDS']}\t")
  out_writer.write(f"{ir_events['3UTR'] + alta_events['3UTR'] + altd_events['3UTR'] + er_events['3UTR'] + ci_events['3UTR'] + es_events['3UTR'] + mxe_events['3UTR']}\t")
  out_writer.write(f"{ci_events['5UTR+CDS'] + es_events['5UTR+CDS'] + mxe_events['5UTR+CDS']}\t")
  out_writer.write(f"{ci_events['CDS+3UTR'] + es_events['CDS+3UTR'] + mxe_events['CDS+3UTR']}\t")
  out_writer.write(f"{ci_events['5UTR+CDS+3UTR'] + es_events['5UTR+CDS+3UTR']}\t")
  out_writer.write(f"{ir_events['not localized'] + alta_events['not localized'] + altd_events['not localized'] + er_events['not localized'] + ci_events['not localized'] + es_events['not localized'] + mxe_events['not localized']}\t")
  out_writer.write(f"{as_duplicates[0] + as_duplicates[1] + as_duplicates[2] + as_duplicates[3] + as_duplicates[5] + as_duplicates[4] + as_duplicates[6]}\t")
  out_writer.write(f"-")

  out_writer.write(f"\n\nTranscripts with unresolved AS events:")
  out_writer.write(f"\nSum: {len(unresolved_list)}")
  #for einzelwert in unresolved_list:
  #  out_writer.write(f"\n{einzelwert}")