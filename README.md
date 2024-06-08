# sjMotifSeeker
Jack Kissinger and Sujana Sreenivasan's SP24 CSE 185 Project

Utility: De Novo Motif Finder for Transcription Factors

Prerequisite Downloads: 
1. **numpy** for data processing: ```pip install numpy```
2. **UCSC Rest API** to query reference genomes: ```pip install ucsc-genomic-api```
3. **argparse** to parse command line arguments: ```pip install argparse```
4. **logomaker** to turn nucleotide frequencies into a spiffy motif logo ```pip install logomaker```
5. **pandas** for dataframe processing: ```pip install pandas```
6. **matplotlib** for displaying motif logo to the screen ```pip install matplotlib```
7. You will need the a reference genome file (mouse genome can be found on DataHub at
```~/public/genomes/GRCm38.fa```)

Build Instructions:
To be able to run the sjMotifSeeker program, run ```pip install . ``` when your
working directory is the root of our repo.

Basic Usage:
```sjMotifSeeker <reference genome filename> -d <peak file name>```

Example Calls:
```sjMotifSeeker GRCm38.fa -c "mm10" -d "Oct4_peaks.txt" -k 8```
```sjMotifSeeker GRCm38.fa -c "mm10" -d "Sox2_peaks.txt"```
```sjMotifSeeker hg19.fa -c "hg19" -d "some_file.bed" -q True```

Detailed Usage:
sjMotifSeeker ref_genome  [-h] -c C -d D [-k K] [-q Q]

positional arguments:
  ref_genome  path to reference genome

options:
  -h, --help  show this help message and exit
  -c C        reference genome code for UCSC genome browser (may not match filename!)
  -d D        space-delimited string of peak file paths
  -k K        size of k-mers to use in analysis of background vs. TF sample (default 7)
  -q Q        set to True if using bed file instead of peaks.txt (default False)

Input Files: tab-delimited text files, columns are 1-indexed
peak (.txt format): 
1. column 2 must be JUST the chromosome number (no "chr")
2. column 3 must be the start index of the peak
3. column 4 must be the start index of the peak

BED format:
1. column 1 must be the chromosome name in format "chrN", where N is the 
number/letter associated
with the chromosome
2. column 2 must be the start index of the peak
3. column 3 must be the start index of the peak

Columns beyond the ones listed are not required for sjMotifSeeker's function,
and can be omitted for memory-saving purposes.

Mixing and matching file formats within one call to sjMotifSeeker is not allowed.
To do this, perform two separate called with different args passed to -q.

Ouput Files:
Motif Logos will be saved as ```input_filename.logo.png```, and also displayed to
the screen before the program terminates. Only one can be viewed at a time.

Development Notes (v2.0):
Implemented -k and -q options, BED file functionality, motif logos, multiple
input files.

To find out the UCSC Genome Browser code for a particular reference genome, 
search [here](https://genome.ucsc.edu/cgi-bin/hgGateway). You will need it upon 
running the sjMotifSeeker command.
