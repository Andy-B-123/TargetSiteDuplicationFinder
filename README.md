# TargetSiteDuplicationFinder
_WorkInProgress!_ At the moment need to work on reducing false-positives. Please verify TSD calls.

A simple script to identify mobile element insertions characterised by Target Site Duplications (TSD) from a short-read data.

### Assumptions
This script is is focused on identifying mobile elements which cause Target Site Duplications, which are often not able to be identified correctly from other variant analysis tools (eg GatK, Smoove, LUMPY).
The following assumptions are made about the file:
 * A single sample mapped to a reference
 * Mapped using BWA-MEM2 (required for soft-clipping information from reads)
 * Sorted and index available
 * Average coverage ~20x[^1]
[^1]: Note that if you have less than this you can adjust the threshold, but much less than 10x is going to be hard to identify TSD reliably
### Requirements 
The parsing script takes advantage of the output from the _very_ fast cluster_identifier tool written for the [Scramble](https://github.com/GeneDx/scramble) program, which identifies ALU, LINE and SVA elements in humans. Please follow the installation instructions on that page to install cluster_identifier (you will need htslib available).

### Install 
The parsing script here uses a small number of non-standard python libraries (see Setup.py for list). Please install using pip with the following approach:
```
git clone https://github.com/Andy-B-123/TargetSiteDuplicationFinder.git
cd TargetSiteDuplicationFinder
pip install .
```

If installed correctly please try running the following:

```
python TSD_parser.cluster_wrapped.py --help
```

```python
usage: TSD_parser.cluster_wrapped.py [-h] [--input_bam INPUT_BAM] [--input_cluster_identifier_file INPUT_CLUSTER_IDENTIFIER_FILE] [--output_base OUTPUT_BASE] [--input_type INPUT_TYPE]

Process BAM files and identify potential Target Site Duplication events.

options:
  -h, --help            show this help message and exit
  --input_bam INPUT_BAM
                        Full path to input .BAM file to process. Assumed sorted and indexed.
  --output_base OUTPUT_BASE
                        Output base name. Default current working directory and "TSD_parse_out"
  --input_type INPUT_TYPE
                        One of SR, LR, AS for Short-Read, Long-Read or Assembly for the type of bam to process.
  --cluster_identifer_path
                        The path to the cluster_identifier executable if not availbale in current PATH. Assumed to be in current PATH.

```

### Running  
If the above is available, please try running with the provided example data:
```
cd TargetSiteDuplicationFinder
python3 TSD_parser.cluster_wrapped.py --input_bam ExampleData/SlimData.BothTargetRegions.sort.bam --output_base TSD_Check
```

### Output  
The output of the program is a slim .bed file with the coordinates of the start and end of the target site duplication event identified above, and the first coordinate from the cluster_identifier output in the 4th column. Sometimes there can be mild differences in the exact coordinates between the parsing script and region identified by cluster_identifier.
For the example data there are only two:
```
> cat TSD_Check.WindowSize_30.bed
scaffold_29     4822338 4822347 scaffold_29:4822338
scaffold_29     6632085 6632099 scaffold_29:6632085
```

This can be loaded into IGV as a .bed track or intersected with an annotation file using bedtools. 

For example, to look at intersections with coding sequences using bedtools from a gff file, first extract only the target annotations of interest (in our case, CDS):
```
grep -P "\tCDS\t" ${file.gff} > ${gff_file.onlyCDS.gff}
bedtools intersect -a TSD_Check.WindowSize_30.bed -b ${gff_file.onlyCDS.gff} -wb > intersections.out
```

### Follow-up

After identifying interesting TSD sites we can use the output from cluster_identifier for identification of potential repeat elements. Using the provided bash script "ConvertClusterIdentifierOutputToFasta.sh" you can convert the output into standard .fasta sequences. Note that this is for ALL candidate sites, and that the ones identified above are a subset of these. 

```
chmod a+x ConvertClusterIdentifierOutputToFasta.sh
./ConvertClusterIdentifierOutputToFasta.sh SlimData.BothTargetRegions.sort.cluster_identifier.out
```

This will produce two fasta files, containing the 'clipped' consensus sequence (eg those in the inesertion region) and the 'anchor' consensus sequence (eg those in the reference sequence).
```
$ head SlimData.BothTargetRegions.sort.cluster_identifier.anchoredConsensus.fasta

>scaffold_29:4821730_right_23_AnchoredReadSeq
gtttggccccgattaaatttaatttttaaatagattggggaggaaactagccagtaataattaattagaagttctgcagcttgtgaggtaaccgttttatgctaatgactttttatgataataaagctaaaatt
>scaffold_29:4821799_left_7_AnchoredReadSeq
tgctgcttttccacataattatgtatagctggagagctcattagacttaacatttatgtattagttattttattactctccatagctacatacatataaaacattgtaactcctagctgaatgtacaaatatac
>scaffold_29:4821901_left_25_AnchoredReadSeq
ttcaataatttttcatcttactgattacgttaaatctagtcgtagtcacactagtcactatgggttttgctaaaagtgcagttaaacatttctttaaatacttactgcagttcattccctgcagatcacttgtggtt
>scaffold_29:4822338_left_14_AnchoredReadSeq
atttgttccgtgcctcgcttgctcttccgctatctcaccgctggcaattactaaactttcgcattttcttaaatattgctctccgatgactgcaggactatctgcaattggaacaataataatacaaagccaagaga
>scaffold_29:4822347_right_10_AnchoredReadSeq
aaaaataaaaataccctttcatataaagagaataatgccatacacgctgtttctaatttgataaggtccatagctcccgtcgtagtcacacgttacaactgggttgccgaggaaatttgttcc


$ head SlimData.BothTargetRegions.sort.cluster_identifier.clipConsensus.fasta

>scaffold_29:4821730_right_23_ClipConsensus
gatttcgtggttgtcgagctgctgcttttccacataattatgtatagctggagagctcattagacttaacatttatgtattagttattttattactctccatagctacctac
>scaffold_29:4821799_left_7_ClipConsensus
aaccgttttatgctaatgactttttatgataataaagctaaaattgatttcgtggttgtcgagc
>scaffold_29:4821901_left_25_ClipConsensus
tagctggagagctcattagacttaacatttatgtattagttattttattactctccatagctacatacatataaaacattgtaactactagctgaatgtacaaatatactttgtaca
>scaffold_29:4822338_left_14_ClipConsensus
tcatctgcatgcgataagaactcaaaaaaaaagttgatctgctagacaatgatctgccacgttcacggacac
>scaffold_29:4822347_right_10_ClipConsensus
atgtcagtgaacttggcagatcattgtctagcagatcaactttttttttgagttcttatcgca
```

You can provide the 'clip' consensus flies to a program like [RepeatMasker](https://github.com/rmhubley/RepeatMasker) to identify repeats which can provide in preliminary identification or use blast+ or [mmseqs](https://github.com/soedinglab/MMseqs2) to search against known repeat databases. Alternatively, you can use these sequences as starting places for targeted assembly using tools such as [aTram](https://github.com/juliema/aTRAM).

The 'anchor' sequences can help identify the surrounding region to the target site insertion to assess for local structural similarities in these areas which may make them susceptible to insertion.
