# sjMotifSeeker
Jack Kissinger and Sujana Sreenivasan's SP24 CSE 185 Project

Utility: De Novo Motif Finder for Transcription Factors

Prerequisite Downloads: 
1. **numpy** for data processing: ```pip install numpy```
2. **UCSC Rest API** to query reference genomes: ```pip install ucsc-genomic-api```
3. **argparse** to parse command line arguments: ```pip install argparse```
4. You will need the GRCm38.fa reference genome file (can be found on DataHub at
```~/public/genomes/GRCm38.fa```)

Build Instructions:
To be able to run the sjMotifSeeker program, run ```pip install . ``` when your
working directory is the root of our repo.

Basic Usage:
```sjMotifSeeker <reference genome filename> -d <peak file name>```

Example Calls:
```sjMotifSeeker GRCm38.fa Oct4_peaks.txt```
```sjMotifSeeker GRCm38.fa Sox2_peaks.txt```

Detailed Usage:
sjMotifSeeker [-h] -d D ref_genome

positional arguments:
  ref_genome  path to reference genome

options:
  -h, --help  show this help message and exit
  -d D        space-delimited string of tag directories

Development Notes (v1.0):
1. Right now, the tool only works with the GRCm38 reference genome. We need to do more work to figure out translation between filenames and UCSC Genome Browser nicknames for the genomes.
2. So far, the k-mer size must be set in the code manually in the ```analyze_peaks()``` function. Only one k-mer size at a time. The top 25 k-mers will be listed in the output, as well as the entire dictionary of all k-mers and their scores for convenience (higher score = more common in TF regions than in background).
3. Only one peak file is supported per call.
4. For runtime purposes, only the first 100 peaks listed in the peaks.txt file are considered (can be changed in ```analyze_peaks()```)