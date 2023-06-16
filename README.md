# TargetSiteDuplicationFinder
A series of scripts to identify mobile element insertions characterised by Target Site Duplications from a short-read alignment file.


### Requiements and Install 
The parsing script takes advantage of the output from the _very_ fast cluster_identifier tool written for the [Scramble](https://github.com/GeneDx/scramble) program, which identifies ALU, LINE and SVA elements in humans. Please follow the installation instructions on that page to install cluster_identifier (you will need htslib available):


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
If the above is available, please try 

### Output  

