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
    parser.add_argument("-k", required = False, type = int, default = 7, help = "size of k-mers to use in analysis of background vs. TF sample")
    parser.add_argument("-q", required = False, type = bool, default = False, help = "set to true if using bed file instead of peaks.txt")


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
    for line in lines[:2100]: # modify this line to help runtime, sacrificing some lines of peaks file
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
                freqs = np.zeros((len(Sequence.get(genome=genome_code,chrom=int(line[0][3:]),start=line[1],end=line[2]).dna), 4), dtype = int)
            else:
                freqs = np.zeros((len(Sequence.get(genome=genome_code,chrom=line[1],start=line[2],end=line[3]).dna), 4), dtype = int)
            prev_hash = False
        if (bed_status):
            seq = Sequence.get(genome=genome_code,chrom=int(line[0][3:]),start=line[1],end=line[2]).dna.upper()[:len(freqs)]
        else:
            seq = Sequence.get(genome=genome_code,chrom=line[1],start=line[2],end=line[3]).dna.upper()
        try:
            for i in range(len(seq)):
                freqs[i][nucs[seq[i]]] += 1
                if (i <= len(seq) - k):
                    k_mer = seq[i:i+k]
                    if k_mer in k_mers:
                        k_mers[k_mer] += 1
                    else:
                        k_mers[k_mer] = 1
            count += 1
        except ValueError:
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
    for i in range(len(lines) // 100000):
        j = i * 99999
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


# parse command line arguments

#https://docs.python.org/3/library/argparse.html


'''
We would take annotated peaks.txt file and reference genome as input 
and output motifs for each transcription factor ranked by the most significant p-values
'''


'''
(Taken from HOMER)
Usage: annotatePeaks.pl <peak file | tss> <genome version>  [additional options...]

        Available Genomes (required argument): (name,org,directory,default promoter set)
                        -- or --
                Custom: provide the path to genome FASTA files (directory or single file)
                If no genome is available, specify 'none'.
                If using FASTA file or none, may want to specify '-organism <...>'

        User defined annotation files (default is UCSC refGene annotation):
                annotatePeaks.pl accepts GTF (gene transfer formatted) files to annotate positions relative
                to custom annotations, such as those from de novo transcript discovery or Gencode.
                -gtf <gtf format file> (Use -gff and -gff3 if appropriate, but GTF is better)
                -gid (by default the GTF file is processed by transcript_id, use this option for gene_id)
                -ann <custom homer annotation file> (created by assignGenomeAnnotation, see website)

        Peak vs. tss/tts/rna mode (works with custom GTF file):
                If the first argument is "tss" (i.e. annotatePeaks.pl tss hg18 ...) then a TSS centric
                analysis will be carried out.  Tag counts and motifs will be found relative to the TSS.
                (no position file needed) ["tts" now works too - e.g. 3' end of gene]
                ["rna" specifies gene bodies, will automaticall set "-size given"]
                NOTE: The default TSS peak size is 4000 bp, i.e. +/- 2kb (change with -size option)
                -list <gene id list> (subset of genes to perform analysis [unigene, gene id, accession,
                         probe, etc.], default = all promoters)
                -cTSS <promoter position file i.e. peak file> (should be centered on TSS)

        Primary Annotation Options:
                -mask (Masked repeats, can also add 'r' to end of genome name)
                -m <motif file 1> [motif file 2] ... (list of motifs to find in peaks)
                        -mscore (reports the highest log-odds score within the peak)
                        -nmotifs (reports the number of motifs per peak)
                        -mdist (reports distance to closest motif)
                        -mfasta <filename> (reports sites in a fasta file - for building new motifs)
                        -fm <motif file 1> [motif file 2] (list of motifs to filter from above)
                        -rmrevopp <#> (only count sites found within <#> on both strands once, i.e. palindromic)
                        -matrix <prefix> (outputs a motif co-occurrence files:
                                prefix.count.matrix.txt - number of peaks with motif co-occurrence
                                prefix.ratio.matrix.txt - ratio of observed vs. expected  co-occurrence
                                prefix.logPvalue.matrix.txt - co-occurrence enrichment
                                prefix.stats.txt - table of pair-wise motif co-occurrence statistics
                                additional options:
                                -matrixMinDist <#> (minimum distance between motif pairs - to avoid overlap, default: 4)
                                -matrixMaxDist <#> (maximum distance between motif pairs)
                        -mbed <filename> (Output motif positions to a BED file to load at UCSC (or -mpeak))
                        -mlogic <filename> (will output stats on common motif orientations)
                -d <tag directory 1> [tag directory 2] ... (list of experiment directories to show
                        tag counts for) NOTE: -dfile <file> where file is a list of directories in first column
                -bedGraph <bedGraph file 1> [bedGraph file 2] ... (read coverage counts from bedGraph files)
                -wig <wiggle file 1> [wiggle file 2] ... (read coverage counts from wiggle files)
                -p <peak file> [peak file 2] ... (to find nearest peaks)
                        -pdist to report only distance (-pdist2 gives directional distance)
                        -pcount to report number of peaks within region
                -vcf <VCF file> (annotate peaks with genetic variation infomation, one col per individual)
                        -editDistance (Computes the # bp changes relative to reference)
                        -individuals <name1> [name2] ... (restrict analysis to these individuals)
                -gene <data file> ... (Adds additional data to result based on the closest gene.
                        This is useful for adding gene expression data.  The file must have a header,
                        and the first column must be a GeneID, Accession number, etc.  If the peak
                        cannot be mapped to data in the file then the entry will be left empty.
                -go <output directory> (perform GO analysis using genes near peaks)
                -genomeOntology <output directory> (perform genomeOntology analysis on peaks)
                        -gsize <#> (Genome size for genomeOntology analysis, default: 2e9)

        Annotation vs. Histogram mode:
                -hist <bin size in bp> (i.e 1, 2, 5, 10, 20, 50, 100 etc.)
                The -hist option can be used to generate histograms of position dependent features relative
                to the center of peaks.  This is primarily meant to be used with -d and -m options to map
                distribution of motifs and ChIP-Seq tags.  For ChIP-Seq peaks for a Transcription factor
                you might want to use the -center option (below) to center peaks on the known motif
                ** If using "-size given", histogram will be scaled to each region (i.e. 0-100%), with
                the -hist parameter being the number of bins to divide each region into.
                        Histogram Mode specific Options:
                        -nuc (calculated mononucleotide frequencies at each position,
                                Will report by default if extracting sequence for other purposes like motifs)
                        -di (calculated dinucleotide frequencies at each position)
                        -histNorm <#> (normalize the total tag count for each region to 1, where <#> is the
                                minimum tag total per region - use to avoid tag spikes from low coverage
                        -ghist (outputs profiles for each gene, for peak shape clustering)
                        -rm <#> (remove occurrences of same motif that occur within # bp)

        Peak Centering: (other options are ignored)
                -center <motif file> (This will re-center peaks on the specified motif, or remove peak
                        if there is no motif in the peak.  ONLY recentering will be performed, and all other
                        options will be ignored.  This will output a new peak file that can then be reanalyzed
                        to reveal fine-grain structure in peaks (It is advised to use -size < 200) with this
                        to keep peaks from moving too far (-mirror flips the position)
                -multi (returns genomic positions of all sites instead of just the closest to center)

        Genome comparisons (need genome & liftOver)
                -cmpGenome <genome1> [genome2] (Genomes to compare for sequence/motifs)
                -cmpLiftover <liftover1> [genome2] (Genomes to compare for sequence/motifs)

        Normalization options:
                -fpkm (normalize read counts to million reads or fragments per kilobase mapped)
                -raw (do not adjust the tag counts based on total tags sequenced, -noadj works too)
                -norm <#> (normalize tags to this tag count, default=1e7, 0=average tag count in all directories)
                -normLength <#> (Fragment length to normlize to for experiments with different lens, def: 100)
                -log (output tag counts as log2(x+1+rand) values - for scatter plots)
                -sqrt (output tag counts as sqrt(x+rand) values - for scatter plots)
                -ratio (process tag values as ratios - i.e. chip-seq, or mCpG/CpG)

        Advanced normalization options: (-rlog and -vst require R and DESeq2 to be installed)
                -rlog (quantile/variance normalization on reported genes using DESeq2 rlog funcition, needs R)
                -vst (quantile/variance normalization on reported genes using DESeq2 vst function, needs R)

        Advanced Options:
                -len <#> / -fragLength <#> (Fragment length, default=auto, might want to set to 1 for 5'RNA)
                -size <#> (Peak size[from center of peak], default=inferred from peak file)
                        -size #,# (i.e. -size -10,50 count tags from -10 bp to +50 bp from center)
                        -size "given" (count tags etc. using the actual regions - for variable length regions)
                -strand <+|-|both> (Count tags on specific strands relative to peak, default: both)
                -pc <#> (maximum number of tags to count per bp, default=0 [no maximum], -tbp <#> works too)
                -CpG (Calculate CpG/GC content)
                -nfr (report nuclesome free region scores instead of tag counts, also -nfrSize <#>)
                -norevopp (do not search for motifs on the opposite strand [works with -center too])
                -gwasCatalog <gwasCatalog file from UCSC> (list overlapping GWAS risk SNPs)
                -pdist (only report distance to nearest peak using -p, not peak name)
                -map <mapping file> (mapping between peak IDs and promoter IDs, overrides closest assignment)
                -noann, -nogene (skip genome annotation step, skip TSS annotation)
                -homer1/-homer2 (by default, the new version of homer [-homer2] is used for finding motifs)
                -cpu <#> (Number of processors to use when possible - only some parts utilize multiple cores)
                -noblanks (remove peaks/rows with missing data)
'''