## Contents
- [Installation](#installation)
- [Running](#running)
  - [Main Modules](#main-modules)
  - [Using `run` submodule to evaluate genome completeness from genome assembly](#using-run-submodule-to-evaluate-genome-completeness-from-genome-assembly)
  - [Using `analysis` submodule to evaluate genome completeness from provided miniprot alignment](#using-analysis-submodule-to-evaluate-genome-completeness-from-provided-miniprot-alignment)
  - [Using `download` submodule to download lineage](#using-download-submodule-to-download-lineage)
  - [Using `list` submodule to list local or remote lineages](#using-list-submodule-to-list-local-or-remote-lineages)
  - [Using `run_miniprot` submodule to run miniprot](#using-run_miniprot-submodule-to-run-miniprot)

Installation
------------
Minibusco is developed on python3. If you only run the `download`, `list`, or `analysis` submodules, you only need to install python3 and pandas.
If you want to run the `run` or `run_miniprot` submodules, you need to install **miniprot**. 
Minubusco will search for miniprot in the directory where `minibusco.py` is located, the current execution directory, and system environment variables.
If you want to use the `--autolineage` mode, you need to install **sepp** and provide the path to sepp `--sepp_execute_path`.
- Prequisites:  
      [python3.*](https://www.python.org)  
      [miniprot](https://github.com/lh3/miniprot) (submodule: run, run_miniprot)  
      [sepp](https://github.com/smirarab/sepp) (autolineage mode)
- Dependencies:  
  [pandas](https://pandas.pydata.org/docs/getting_started/install.html#installing-from-pypi)  

## Running

### Main Modules:

```angular2html
download        Download specified BUSCO lineage
list            List local or remote BUSCO lineages
run_miniprot    Run miniprot alignment
analysis        Evaluate genome completeness from provided miniprot alignment
run             Run minibusco including miniprot alignment and completeness evaluation
```

### Using `run` submodule to evaluate genome completeness from genome assembly:

This will download the specified lineage (or automatically search for the best lineage with autolineage mode), align the
protein sequences in the lineage file to the genome sequence with miniprot, and parse the miniprot alignment result to
evaluate genome completeness.
#### Usage:
```angular2html
python minibusco.py run -a ASSEMBLY_PATH -o OUTPUT_DIR 
                        [-t THREADS] [-l LINEAGE] [--library_path LIBRARY_PATH] [--specified_contigs SPECIFIED_CONTIGS [SPECIFIED_CONTIGS ...]] 
                        [--miniprot_execute_path MINIPROT_EXECUTE_PATH] [--autolineage] [--sepp_execute_path SEPP_EXECUTE_PATH] 
                        [--min_diff MIN_DIFF] [--min_identity MIN_IDENTITY] [--min_length_percent MIN_LENGTH_PERCENT] 
                        [--min_complete MIN_COMPLETE] [--min_rise MIN_RISE][-h] 
```

#### Important parameters:
```angular2html
  -a, --assembly_path        Input genome file in FASTA format
  -o, --output_dir           The output folder
  -t, --threads              Number of threads to use
  -l, --lineage              Specify the name of the BUSCO lineage to be used. (e.g. eukaryota, primates, saccharomycetes etc.)
  --library_path             Folder path to download lineages or already downloaded lineages. 
                             If not specified, a folder named "mb_downloads" will be created on the current running path by default to store the downloaded lineage files.
  --miniprot_execute_path    Path to miniprot executable file. 
                             If not specified, minibusco will search for miniprot in the directory where minibusco.py is located, the current execution directory, and system environment variables.
  --autolineage              Automatically search for the best matching lineage without specifying lineage file.
  --sepp_execute_path        Path to sepp executable file. This is required if you want to use the autolineage mode.
  --specified_contigs        Specify the contigs to be evaluated, e.g. chr1 chr2 chr3. If not specified, all contigs will be evaluated.
```

#### Threshold parameters:
```angular2html
  --min_diff               The thresholds for the best matching and second best matching. default=0.2
  --min_identity           The identity threshold for valid mapping results. default=0.4
  --min_length_percent     The fraction of protein for valid mapping results. default=0.6
  --min_complete           The length threshold for complete gene. default=0.9
```

#### Example:
```angular2html
# with lineage specified
python minibusco.py run -a genome.fasta -o output_dir -l eukaryota -t 8

# autolineage mode
python minibusco.py run -a genome.fasta -o output_dir -t 8 --autolineage ---sepp_execute_path /path/to/run_sepp.py

# with custom specified already downloaded lineage folder
python minibusco.py run -a genome.fasta -o output_dir -l eukaryota -t 8 --library_path /path/to/lineages_folder

# specify contigs
python minibusco.py run -a genome.fasta -o output_dir -l eukaryota -t 8 --specified_contigs chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
```

### Using `analysis` submodule to evaluate genome completeness from provided miniprot alignment:
This will directly parse the provided miniprot alignment result to evaluate genome completeness. The execute command of miniprot should be like `miniprot -u -I --outs=0.95 --gff -t 8 ref-file protein.faa > output.gff`.
#### Usage:
```angular2html
python minibusco.py analysis -g GFF -o OUTPUT_DIR [--specified_contigs SPECIFIED_CONTIGS [SPECIFIED_CONTIGS ...]] 
                        [--min_diff MIN_DIFF] [--min_identity MIN_IDENTITY] [--min_length_percent MIN_LENGTH_PERCENT] 
                        [--min_complete MIN_COMPLETE] [--min_rise MIN_RISE] [-h]
```
#### Important parameters:
```angular2html
  -g, --gff                 Miniprot output gff file
  -o, --output_dir          Output analysis folder
  --specified_contigs       Specify the contigs to be evaluated, e.g. chr1 chr2 chr3. If not specified, all contigs will be evaluated.
```
Threshold parameters are same as `run` module.

#### Example:
```angular2html
# analysis with miniprot output gff file
python minibusco.py analysis -g miniprot.gff -o output_dir

# specify contigs
minibusco analysis -g miniprot.gff -o output_dir --specified_contigs chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
```

### Using `download` submodule to download lineage:
This will download the specified lineages and save to the specified folder.
#### Usage:
```angular2html
python minibusco.py download [-h] -l LINEAGE [LINEAGE ...] --library_path LIBRARY_PATH
```

#### Important parameters:
```angular2html
  -l, --lineage       Specify the names of the BUSCO lineages to be downloaded. (e.g. eukaryota, primates, saccharomycetes etc.)
  --library_path      The destination folder to store the downloaded lineage files.
```

#### Example:
```angular2html
python minibusco.py download -l primates brassicales -d /path/to/lineages_folder
```

### Using `run_miniprot` submodule to run miniprot alignment:
This will run miniprot alignment and output the gff file.
#### Usage:
```angular2html
python minibusco.py run_miniprot [-h] -a ASSEMBLY -p PROTEIN -o OUTDIR [-t THREADS] [--exec_path EXEC_PATH]
```

#### Important parameters:
```angular2html
  -a, --assembly        Input genome file in FASTA format
  -p, --protein         Input protein file
  -o, --outdir          Miniprot alignment output directory
  -t, --threads         Number of threads to use
  --exec_path           Path to miniprot executable file. 
                        If not specified, minibusco will search for miniprot in the directory where minibusco.py is located, the current execution directory, and system environment variables.
```

#### Example:
```angular2html
python minibusco.py run_miniprot -a genome.fasta -p protein.faa -o output_dir -t 8
```

### Using `list` submodule to show local or remote Busco lineages:
This will list the local or remote BUSCO lineages.
#### Usage:
```angular2html
python minibusco.py list [-h] [--remote] [--local] [--library_path LIBRARY_PATH]
```

#### Important parameters:
```angular2html
  --remote         List remote BUSCO lineages
  --local          List local BUSCO lineages
  --library_path   Folder path to stored lineages.
```

#### Example
```angular2html
# list local lineages
python minibusco.py list --local --library_path /path/to/lineages_folder

# list remote lineages
python minibusco.py list --remote
```

