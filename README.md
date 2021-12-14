# GPUDePiCT
# Language: CUDA
# Input: TXT
# Output: SCREEN
# Tested with: PluMA 1.0, CUDA 10

PluMA plugin for degenerate primer construction on the GPU (Cickovski et al, 2015)

The plugin accepts as input a tab-delimited file of keyword-value pairs:
inputfile: Sequence data (acceptable formats: .txt, .msf, .fasta, .sto)
rows: Number of sequences
seqlength: Length of each sequenes
AAorNuc: Amino acid (0) or nucleotide (1)
kmeans: Doing K-Means?  No (0) or Yes (1)
fuzzy: Doing Fuzzy C-Means?  If so provide fuzziness (0 for standard or K-Means)
tolerance: Tolerance if using Fuzzy C-Means.  0 for standard or K-Means
minlength: Minimum primer length.

Primers are output to the screen.  Supply 'none' for the output.
