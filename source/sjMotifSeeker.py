import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ucsc.api import Sequence
import logomaker

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("ref_genome", help = "path to reference genome")
    parser.add_argument("-c", required = True, type = str, help = "reference genome code for UCSC genome browser")
    parser.add_argument("-d", required = True, type = str, help = "space-delimited string of peak file paths")
    parser.add_argument("-k", required = False, type = int, default = 7, help = "size of k-mers to use in analysis of background vs. TF sample (default 7)")
    parser.add_argument("-q", required = False, type = bool, default = False, help = "set to True if using bed file instead of peaks.txt (default False)")


    args = parser.parse_args()
    peak_files = args.d.split()
    ref_genome = args.ref_genome
    ucsc_code = args.c
    k = args.k
    bed_status = args.q

    for peak_file in peak_files:
        print("Starting processing on", peak_file,"\n")
        motifs = analyze_peaks(peak_file, ref_genome, k, ucsc_code, bed_status)
        print("Generating motif logo...")
        for motif in motifs:
            generate_motif_image(motif, peak_file)

def analyze_peaks(peak_file, ref_genome, k, genome_code, bed_status):
    nucs = {"A": 0, "T": 1, "G": 2, "C": 3}
    f = open(peak_file, "r")
    lines = f.readlines()
    prev_hash = True
    k_mers = {}
    print("Generating kmers around TF binding sites...")
    count = 0
    for line in lines: # modify this line to help runtime, sacrificing some lines of peaks file
        if line.startswith("#"):
            continue
        if count == len(lines) // 4:
            print("25% complete")
        elif count == len(lines) // 2:
            print("50% complete")
        elif count == 3 * len(lines) // 4:
            print("75% complete")
        line = line.split("\t")
        if (prev_hash):
            if (bed_status):
                freqs = np.zeros((len(Sequence.get(genome=genome_code,chrom=line[0][3:],start=line[1],end=line[2]).dna), 4), dtype = int)
            else:
                freqs = np.zeros((len(Sequence.get(genome=genome_code,chrom=line[1],start=line[2],end=line[3]).dna), 4), dtype = int)
            prev_hash = False
        try:
            if (bed_status):
                seq = Sequence.get(genome=genome_code,chrom=line[0][3:],start=line[1],end=line[2]).dna.upper()[:len(freqs)]
            else:
                seq = Sequence.get(genome=genome_code,chrom=line[1],start=line[2],end=line[3]).dna.upper()
            for i in range(len(seq)):
                freqs[i][nucs[seq[i]]] += 1
                if (i <= len(seq) - k):
                    k_mer = seq[i:i+k]
                    if k_mer in k_mers:
                        k_mers[k_mer] += 1
                    else:
                        k_mers[k_mer] = 1
            count += 1
        except:
            continue
    tf_binding_kmers = dict(sorted(k_mers.items(), key=lambda item: item[1]))
    background_kmers = load_comparison(ref_genome, k)

    #normalization
    print("Normalizing...")
    size_data = 0
    for key in tf_binding_kmers:
        size_data += tf_binding_kmers[key]

    size_background = 0
    for key in background_kmers:
        size_background += background_kmers[key]

    norm_ratio = size_data/size_background

    for key in background_kmers:
        background_kmers[key] *= norm_ratio

    adjusted_kmers = {}
    for kmer in tf_binding_kmers:
        if kmer in background_kmers:
            adjusted_kmers[kmer] = round(tf_binding_kmers[kmer] - background_kmers[kmer], 3)
        else:
            adjusted_kmers[kmer] = 0
    adjusted_kmers = dict(sorted(adjusted_kmers.items(), key=lambda item: item[1], reverse = True))

    print("Top 25 Most Overrepresented " + str(k) + "-mers:")
    counter = 0
    for key in adjusted_kmers:
        print(key, adjusted_kmers[key])
        counter+= 1
        if counter >= 25:
            break

    print("Merging similar kmers into motifs...")
    motifs = merge_kmers(adjusted_kmers, k)
    for motif in motifs:
        for motif_info_block in motif:
            print("|".join([*motif_info_block[0]]), motif_info_block[3])

    return motifs

def load_comparison(ref_filename, k):
    print("Loading reference genome file...")
    f = open(ref_filename, "r")
    lines = f.readlines()
    print("Generating background kmers...")
    k_mer_dct = {}
    seq = ""
    for i in range(len(lines) // 50000):
        j = i * 49999
        if lines[j].startswith(">"):
            continue
        seq += lines[j].strip()

    for i in range(len(seq)):
        if (i <= len(seq) - k):
            k_mer = seq[i:i+k]
            if "N" in k_mer.upper():
                continue
            if k_mer in k_mer_dct:
                k_mer_dct[k_mer] += 1
            else:
                k_mer_dct[k_mer] = 1
    return dict(sorted(k_mer_dct.items(), key=lambda item: item[1]))

# call only on dictionaries sorted from highest to lowest score
def merge_kmers(kmer_dct, k):
    motifs = []
    forward_seen = []
    reverse_seen = []

    # FORWARD PASS 1
    pseudoindex = 1
    for key in kmer_dct:
        score = kmer_dct[key]
        # create root to build off of
        if pseudoindex == 1:
            motifs.append([[key, 0, k - 1, score]])
            forward_seen.append(key)
            pseudoindex += 1
            continue

        elif pseudoindex >= 50:
            break

        for motif in motifs:
            for seq_info in motif:
                seq_complete, start, end, seq_score = seq_info
                seq = seq_complete[start: end + 1]

                #CHECK CASE 1 (shifted 1 position)
                if (key[1:] == seq[:-1]):
                    if start == 0: # extend left
                        for info_block in motif:
                            info_block[0] = "-" + info_block[0]
                            info_block[1] += 1
                            info_block[2] += 1
                        start += 1
                        end += 1
                        seq_complete = "-" + seq_complete

                        constructed_key = ""
                        for i in range(len(seq_complete)):
                            if i < start - 1:
                                constructed_key += "-"
                            elif i < end:
                                constructed_key += key[i - start + 1]
                            else:                                
                                constructed_key += "-"
                        motif.append([constructed_key, start - 1, end - 1, score])

                    else: # no extension
                        constructed_key = ""
                        for i in range(len(seq_complete)):
                            if i < start - 1:
                                constructed_key += "-"
                            elif i < end:
                                constructed_key += key[i - start + 1]
                            else:                                
                                constructed_key += "-"
                        motif.append([constructed_key, start - 1, end - 1, score])
                    forward_seen.append(key)
                    break

                if (key[:-1] == seq[1:]):
                    if end == len(seq_complete) - 1: # extend right
                        for info_block in motif:
                            info_block[0] += "-"
                        seq_complete += "-"

                        constructed_key = ""
                        for i in range(len(seq_complete)):
                            if i < start + 1:
                                constructed_key += "-"
                            elif i < end + 2:
                                constructed_key += key[i - start - 1]
                            else:                                
                                constructed_key += "-"
                        motif.append([constructed_key, start + 1, end + 1, score])

                    else: # no extension
                        constructed_key = ""
                        for i in range(len(seq_complete)):
                            if i < start - 1:
                                constructed_key += "-"
                            elif i < end:
                                constructed_key += key[i - start - 1]
                            else:                                
                                constructed_key += "-"
                        motif.append([constructed_key, start + 1, end + 1, score])
                    forward_seen.append(key)
                    break

                #Check Case 2 (1 mismatch)
                if num_mismatch_same_len(key, seq) == 1: 
                    constructed_key = ""
                    for i in range(len(seq_complete)):
                        if i < start:
                            constructed_key += "-"
                        elif i <= end:
                            constructed_key += key[i - start]
                        else:                                
                            constructed_key += "-"
                    motif.append([constructed_key, start, end, score])
                    forward_seen.append(key)
                    break

        pseudoindex += 1

    # REVERSE PASS 1
    pseudoindex = 1
    for key in kmer_dct:
        score = kmer_dct[key]
        key = rev_com(key)
        if pseudoindex >= 50:
            break

        for motif in motifs:
            for seq_info in motif:
                seq_complete, start, end, seq_score = seq_info
                seq = seq_complete[start: end + 1]

                #CHECK CASE 1 (shifted 1 position)
                if (key[1:] == seq[:-1]):
                    if start == 0: # extend left
                        for info_block in motif:
                            info_block[0] = "-" + info_block[0]
                            info_block[1] += 1
                            info_block[2] += 1
                        start += 1
                        end += 1
                        seq_complete = "-" + seq_complete

                        constructed_key = ""
                        for i in range(len(seq_complete)):
                            if i < start - 1:
                                constructed_key += "-"
                            elif i < end:
                                constructed_key += key[i - start + 1]
                            else:                                
                                constructed_key += "-"
                        motif.append([constructed_key, start - 1, end - 1, score])

                    else: # no extension
                        constructed_key = ""
                        for i in range(len(seq_complete)):
                            if i < start - 1:
                                constructed_key += "-"
                            elif i < end:
                                constructed_key += key[i - start + 1]
                            else:                                
                                constructed_key += "-"
                        motif.append([constructed_key, start - 1, end - 1, score])
                    reverse_seen.append(key)
                    break

                if (key[:-1] == seq[1:]):
                    if end == len(seq_complete) - 1: # extend right
                        for info_block in motif:
                            info_block[0] += "-"
                        seq_complete += "-"

                        constructed_key = ""
                        for i in range(len(seq_complete)):
                            if i < start + 1:
                                constructed_key += "-"
                            elif i < end + 2:
                                constructed_key += key[i - start - 1]
                            else:                                
                                constructed_key += "-"
                        motif.append([constructed_key, start + 1, end + 1, score])

                    else: # no extension
                        constructed_key = ""
                        for i in range(len(seq_complete)):
                            if i < start - 1:
                                constructed_key += "-"
                            elif i < end:
                                constructed_key += key[i - start - 1]
                            else:                                
                                constructed_key += "-"
                        motif.append([constructed_key, start + 1, end + 1, score])
                    reverse_seen.append(key)
                    break

                #Check Case 2 (1 mismatch)
                if num_mismatch_same_len(key, seq) == 1: 
                    constructed_key = ""
                    for i in range(len(seq_complete)):
                        if i < start:
                            constructed_key += "-"
                        elif i <= end:
                            constructed_key += key[i - start]
                        else:                                
                            constructed_key += "-"
                    motif.append([constructed_key, start, end, score])
                    reverse_seen.append(key)
                    break

        pseudoindex += 1

    # FORWARD PASS 2
    pseudoindex = 1
    for key in kmer_dct:
        score = kmer_dct[key]
        if pseudoindex >= 50:
            break

        for motif in motifs:
            if key in forward_seen:
                continue
            for seq_info in motif:
                seq_complete, start, end, seq_score = seq_info
                seq = seq_complete[start: end + 1]

                #CHECK CASE 1 (shifted 1 position)
                if (key[1:] == seq[:-1]):
                    if start == 0: # extend left
                        for info_block in motif:
                            info_block[0] = "-" + info_block[0]
                            info_block[1] += 1
                            info_block[2] += 1
                        start += 1
                        end += 1
                        seq_complete = "-" + seq_complete

                        constructed_key = ""
                        for i in range(len(seq_complete)):
                            if i < start - 1:
                                constructed_key += "-"
                            elif i < end:
                                constructed_key += key[i]
                            else:                                
                                constructed_key += "-"
                        motif.append([constructed_key, start - 1, end - 1, score])

                    else: # no extension
                        constructed_key = ""
                        for i in range(len(seq_complete)):
                            if i < start - 1:
                                constructed_key += "-"
                            elif i < end:
                                constructed_key += key[i - start + 1]
                            else:                                
                                constructed_key += "-"
                        motif.append([constructed_key, start - 1, end - 1, score])
                    forward_seen.append(key)
                    break

                #Check Case 2 (1 mismatch)
                if num_mismatch_same_len(key, seq) == 1: 
                    constructed_key = ""
                    for i in range(len(seq_complete)):
                        if i < start:
                            constructed_key += "-"
                        elif i <= end:
                            constructed_key += key[i - start]
                        else:                                
                            constructed_key += "-"
                    motif.append([constructed_key, start, end, score])
                    forward_seen.append(key)
                    break
        pseudoindex += 1

    # REVERSE PASS 2
    pseudoindex = 1
    for key in kmer_dct:
        if pseudoindex >= 50:
            break

        for motif in motifs:
            if key in reverse_seen:
                continue
            for seq_info in motif:
                seq_complete, start, end, seq_score = seq_info
                seq = seq_complete[start: end + 1]

                #CHECK CASE 1 (shifted 1 position)
                if (key[1:] == seq[:-1]):
                    if start == 0: # extend left
                        for info_block in motif:
                            info_block[0] = "-" + info_block[0]
                            info_block[1] += 1
                            info_block[2] += 1
                        start += 1
                        end += 1
                        seq_complete = "-" + seq_complete

                        constructed_key = ""
                        for i in range(len(seq_complete)):
                            if i < start - 1:
                                constructed_key += "-"
                            elif i < end:
                                constructed_key += key[i]
                            else:                                
                                constructed_key += "-"
                        motif.append([constructed_key, start - 1, end - 1, score])

                    else: # no extension
                        constructed_key = ""
                        for i in range(len(seq_complete)):
                            if i < start - 1:
                                constructed_key += "-"
                            elif i < end:
                                constructed_key += key[i - start + 1]
                            else:                                
                                constructed_key += "-"
                        motif.append([constructed_key, start - 1, end - 1, score])
                    reverse_seen.append(key)
                    break

                #Check Case 2 (1 mismatch)
                if num_mismatch_same_len(key, seq) == 1: 
                    constructed_key = ""
                    for i in range(len(seq_complete)):
                        if i < start:
                            constructed_key += "-"
                        elif i <= end:
                            constructed_key += key[i - start]
                        else:                                
                            constructed_key += "-"
                    motif.append([constructed_key, start, end, score])
                    reverse_seen.append(key)
                    break
        pseudoindex += 1
    return motifs

def generate_motif_image(motif, peak_file):
    height = len(motif)
    cutoff = height * 0.8
    pfm = []
    blanks = []
    for motif_info_block in motif:
        seq = motif_info_block[0]
        if len(blanks) == 0:
            for nuc in seq:
                blanks.append(0)
                pfm.append([0,0,0,0])
        for i in range(len(seq)):
            if seq[i] == "-":
                blanks[i] += 1
            if seq[i] == "A":
                pfm[i][0] += 1
            if seq[i] == "C":
                pfm[i][1] += 1
            if seq[i] == "G":
                pfm[i][2] += 1
            if seq[i] == "T":
                pfm[i][3] += 1
    
    trimmed_pfm = []
    for i in range(len(blanks)):
        if blanks[i] < cutoff:
            trimmed_pfm.append(pfm[i])

    for row in trimmed_pfm:
        total = sum(row)
        for i in range(len(row)):
            row[i] = row[i] / total

    print(trimmed_pfm)
    trimmed_pfm = pd.DataFrame(trimmed_pfm, dtype = float)
    trimmed_pfm.columns = ["A", "C", "G", "T"]
    print(trimmed_pfm)
    ss_logo = logomaker.Logo(trimmed_pfm,
                        width=.4,
                        vpad=.05,
                        fade_probabilities=False,
                        stack_order='small_on_top',
                        color_scheme='classic')

    # style using Logo methods
    ss_logo.style_spines(spines=['left', 'right'], visible=False)

    # style using Axes methods
    ss_logo.ax.set_yticks([0, .5, 1])
    ss_logo.ax.axvline(2.5, color='k', linewidth=1, linestyle=':')
    ss_logo.ax.set_ylabel('probability')

    plt.savefig(peak_file + '.logo.png')
    plt.show()
    return None

def num_mismatch_same_len(seq1, seq2):
    mismatch = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            mismatch += 1
    return mismatch

def rev_com(seq):
    reverse = ""
    for char in seq:
        if char == "A":
            reverse = "T" + reverse
        elif char == "T":
            reverse = "A" + reverse
        elif char == "G":
            reverse = "C" + reverse
        else:
            reverse = "G" + reverse
    return reverse
    
if __name__ == "__main__":
    main()