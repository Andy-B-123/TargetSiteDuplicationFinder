# TargetSiteDuplicationFinder
A series of scripts to identify mobile element insertions characterised by Target Site Duplications from a short-read alignment file.

### Assumptions
This script is is focused on identifying mobile elements which cause Target Site Duplications, which are often not able to be identified correctly from other variant analysis tools.
The following assumptions are made about the file:
 * A single sample mapped to a reference
 * Mapped using BWA-MEM2 (required for soft-clipping from reads)
 * Sorted and index available
 * Average coverage ~20x

### Requirements 
The parsing script takes advantage of the output from the _very_ fast cluster_identifier tool written for the [Scramble](https://github.com/GeneDx/scramble) program, which identifies ALU, LINE and SVA elements in humans. Please follow the installation instructions on that page to install cluster_identifier (you will need htslib available).

### Install 
The parsing script here uses a small number of non-standard python libraries. Please install to user (if you do not have root permissions) for the following:
```
pip install tqdm --user
pip install icecream --user
pip install pandas --user
pip install pysam --user
pip install numpy --user
```

The python script does not require installation, please just clone the repository and run:
```
git clone https://github.com/Andy-B-123/TargetSiteDuplicationFinder.git
cd TargetSiteDuplicationFinder
python TSD_parser.py --help
```

```python
usage: TSD_parser.py [-h] [--input_bam INPUT_BAM] [--input_cluster_identifier_file INPUT_CLUSTER_IDENTIFIER_FILE] [--output_base OUTPUT_BASE] [--input_type INPUT_TYPE]

Process BAM files and identify potential Target Site Duplication events.

options:
  -h, --help            show this help message and exit
  --input_bam INPUT_BAM
                        Full path to input .BAM file to process. Assumed sorted and indexed.
  --input_cluster_identifier_file INPUT_CLUSTER_IDENTIFIER_FILE
                        Full path to the output of cluster_identifier run on your .BAM file. ASSUMES YOU HAVE RUN IT ON THE SAME BAM FILE AS PROVIDED ABOVE
  --output_base OUTPUT_BASE
                        Output base name. Default current working directory and "TSD_parse_out"
  --input_type INPUT_TYPE
                        One of SR, LR, AS for Short-Read, Long-Read or Assembly for the type of bam to process.

```

### Running  
If the above is available, please try running with the provided example data:
```
cd TargetSiteDuplicationFinder
python3 TSD_parser.py --input_bam ExampleData/SlimData.BothTargetRegions.sort.bam --input_cluster_identifier_file ExampleData/SlimData.BothTargetRegions.sort.clusterIdentifier.out --output_base TSD_Check
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

### Follow-up
After identifying interesting TSD sites we can use the output from cluster_identifier for identification of potential repeat elements (or other investigations).  
