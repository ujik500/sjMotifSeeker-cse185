import argparse

# parse command line arguments

#https://docs.python.org/3/library/argparse.html

parser = argparse.ArgumentParser()
parser.add_argument("dir")
parser.add_argument("-m", required = True, type = int)
parser.add_argument("-s", required = True, type = int) 
parser.add_argument("-d", required = True, type = int)
parser.add_argument("-a", action = "store_true")

args = parser.parse_args()

filename = args.dir
match_score = args.m
mismatch_penalty = args.s
indel_penalty = args.d

'''
We would take annotated peaks.txt file and reference genome as input 
and output motifs for each transcription factor ranked by the most significant p-values
'''