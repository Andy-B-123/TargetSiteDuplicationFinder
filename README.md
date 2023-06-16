# TargetSiteDuplicationFinder
A series of scripts to identify mobile element insertions characterised by Target Site Duplications from a short-read alignment file.


### Requiements and Install 
The parsing script takes advantage of the _very_ fast cluster_identifier tool written for the [Scramble](https://github.com/GeneDx/scramble) tool for mobile element insertion finder specifically for humans. Please follow the installation instructions on that page to install cluster_identifier:
```
$ cd cluster_identifier/src
$ make
```

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

```

### Running  
If the above is available, please try 

### Output  

