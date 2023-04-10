Installation
------------

```
# create a conda environment and install dependencies
conda env create -f requirements.yml
# activate the environment
conda activate minibusco
# install minibusco
python setup.py install
```

## Running

### Main Modules:

```angular2html
run - Run minibusco including miniprot alignment and completeness evaluation
download - Download specified BUSCO lineage
analysis - Evaluate genome completeness from provided miniprot alignment
```

### Using minibusco run to evaluate genome completeness:

This will download the specified lineage (or automatically search for the best lineage with autolineage mode), align the
protein sequences in the lineage file to the genome sequence with miniprot, and parse the miniprot alignment result to
evaluate genome completeness.

#### Important parameters:
```angular2html
  -a, --assembly_path    Input genome file in FASTA format.
  -o, --output_dir       The output folder.
  -t, --threads          Number of threads to use
  -l, --lineage          Specify the name of the BUSCO lineage to be used. (e.g. eukaryota, primates, saccharomycetes etc.)
  --library_path         Folder path to download lineages or already downloaded lineages. 
                         If not specified, a folder named "downloads" will be created on the current running path by default to store the downloaded lineage files.
  --autolineage          Automatically search for the best matching lineage without specifying lineage file.
  --specified_contigs    Specify the contigs to be evaluated, e.g. chr1 chr2 chr3. If not specified, all contigs will be evaluated.
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
minibusco run -a genome.fasta -o output_dir -l eukaryota -t 8

# autolineage mode
minibusco run -a genome.fasta -o output_dir -t 8 --autolineage

# with custom specified already downloaded lineage folder
minibusco run -a genome.fasta -o output_dir -l eukaryota -t 8 --library_path /path/to/lineages_folder

# specified contigs
minibusco run -a genome.fasta -o output_dir -l eukaryota -t 8 --specified_contigs chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
```

### Using minibusco analysis to evaluate genome completeness:
This will directly parse the provided miniprot alignment result to evaluate genome completeness. The execute command of miniprot should be like `miniprot --gff ref-file protein.faa > output.gff`.
#### Important parameters:
```angular2html
  -g, --gff                 Miniprot output gff file
  -o, --output_dir          Output analysis folder
  -p,--path_to_proteins     Path to protein sequence file in lineage for count the number of total genes
  --specified_contigs       Specify the contigs to be evaluated, e.g. chr1 chr2 chr3. If not specified, all contigs will be evaluated.
```
Threshold parameters are same as `minibuso run` module.

#### Example:
```angular2html
# specified protein sequence file
minibusco analysis -g miniprot.gff -o output_dir -p /path/to/lineage_folder/lineage_name/refseq_db.faa.gz

# without specified protein sequence file
minibusco analysis -g miniprot.gff -o output_dir

# specified contigs
minibusco analysis -g miniprot.gff -o output_dir --specified_contigs chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
```

### Using minibusco download to download lineage:
This will download the specified lineage (or automatically search for the best lineage with autolineage mode) to the specified folder.
#### Important parameters:
```angular2html
  -l, --lineage          Specify the name of the BUSCO lineage to be downloaded. (e.g. eukaryota, primates, saccharomycetes etc.)
  -d, --destination      The destination folder to store the downloaded lineage files.
```

#### Example:
```angular2html
minibusco download -l primates -d /path/to/lineages_folder
```

