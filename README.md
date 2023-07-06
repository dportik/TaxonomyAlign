# TaxonomyAlign

Taxonomy-guided multiple sequence alignment for difficult genes, including 12S and 16S rRNA.


---------------

## Overview

Use to perform multiple sequence alignment for a highly variable gene and a large number of divergent sequences. 

**TaxonomyAlign** automatically divides the sequences into smaller subsets based on user-supplied taxonomic groupings, and builds sub-alignments for the groups using MAFFT (several algorithm options available). 

The defined groupings should be related to both taxonomy and phylogeny (e.g., genus, family, superfamily, etc). 

The resulting sub-alignments are then merged together with MAFFT, with or without an optional polishing step. If a grouping only contains a single sequence, that sequence will be added to the alignments during the merge step. 

The original use-case was for creating high-quality alignments from many 12S and 16S mtDNA sequences across a large phylogenetic scale. However, **TaxonomyAlign** should prove useful for other genes with many species and divergent sequences to align.


## Version

The current release of **TaxonomyAlign** is [v1.0](https://github.com/dportik/TaxonomyAlign/releases/).

## Installation

**TaxonomyAlign** is written in Python (compatible with 3.7) and functions as a command-line script (`TaxonomyAlign.py`). It can be downloaded and executed independently without the need to install **TaxonomyAlign** as a Python package or library. 


**TaxonomyAlign** requires MAFFT and BioPython (v1.77).
Installation of these requirements is fast and easy using `conda`. The [taxonomyalign-conda-env.yml](https://github.com/dportik/TaxonomyAlign/blob/main/taxonomyalign-conda-env.yml) file can be used to create the correct conda environment:

```
conda create -f taxonomyalign-conda-env.yml
```

The resulting conda environment can then be activated using:

```
conda activate taxonomyalign
```

You can then run **TaxonomyAlign** in this environment. 


## Citation

**TaxonomyAlign** publication coming soon.


## Instructions for Analyses

**TaxonomyAlign** provides an alternative to other alignment methods. It uses taxonomy-guided multiple sequence alignment for difficult genes. It was designed for large-scale 12S and 16S rRNA alignments for vertebrates (>3,000 sequences). However, it should prove useful for other genes with many species and divergent sequences to align. 

Example datasets are provided in the [example data folder](https://github.com/dportik/TaxonomyAlign/tree/main/example-data). There are two fasta files containing either 12S or 16S sequences, for all frogs. The map file (`Anura-Family-Map-File.txt`) can be used with either file to run the analysis.


### Quick Navigation for Instructions:

+ [**1. Create a taxonomy map file**](#OLTP)
+ [**2. Running the analysis**](#DGTP)
+ [**2. Outputs**](#OTP)


### 1. Create a taxonomy map file <a name="OLTP"></a>

This is a tab-delimited text file with two columns. The first column should contain the name of the taxon/sequence. The second column should contain the taxonomic group name that the taxon is a part of. This could be genus, family, or higher rank. For frogs, I used family level groupings. An example of this is shown below, and the full file can be found in the [example data folder](https://github.com/dportik/TaxonomyAlign/tree/main/example-data).

```
Semnodactylus_wealii	Hyperoliidae
Tachycnemis_seychellensis	Hyperoliidae
Leiopelma_archeyi	Leiopelmatidae
Leiopelma_hamiltoni	Leiopelmatidae
Leiopelma_hochstetteri	Leiopelmatidae
Leiopelma_pakeka	Leiopelmatidae
Adenomera_ajurauna	Leptodactylidae
Adenomera_andreae	Leptodactylidae
Adenomera_araucaria	Leptodactylidae
```

You can use any method to generate the map file, but it is important to ensure all sequences are assigned. Any sequence not appearing in this file will be excluded from the alignment. The map file can contain species that do not appear in your fasta file. 

There is a quick way to obtain all of the sequence names from your unaligned fasta file using the `--output_seqnames` flag: 

```
python TaxonomyAlign.py -f <fasta file> -m <map file> -o <output directory> --output_seqnames
```

The output from this option is a text file of the sequence names. This can be used to create the map file, but you will need to add the groupings.



### 2. Running the analysis <a name="DGTP"></a>

Once the map file has been created, the analysis can be run. The general usage is shown below:

```
python TaxonomyAlign.py -f <fasta file> -m <map file> -o <output directory> -a <alignment algorithm> -t <threads>
```

#### Argument Explanations:

##### `-f <path-to-file>`

> The full path to an unaligned fasta file that contains DNA sequences.

##### `-m <path-to-file>`

> The map file described in the above section.

##### `-o <path-to-directory>`

> The full path to an existing directory to write the output files.

##### `-a <choice>`

> The alignment algorithm to use. Choices include auto, FFT-NS-i, E-INS-i, L-INS-i, and G-INS-i. Defaults to auto.
auto: auto-select based on input sequences
FFT-NS-i: a progressive and general algorithm
E-INS-i: best for several conserved core regions and variable regions (e.g., 12S and 16S mtDNA)
L-INS-i: best for a single conserved core region with variable flanking regions (e.g., UCEs)
G-INS-i: best for full-length single core sequences

##### `-t <integer>`

> Specifies number of threads to use.


#### Example Use:

```
python TaxonomyAlign.py -f bin/Analysis/12S-Anura.fasta -m bin/Analysis/Anura-Family-Map-File.txt -o bin/Analysis/Output1/ -a E-INS-i -t 8
```

The above command will align the fasta file `12S-Anura.fasta` using the E-INS-i algorithm, creating sub-alignments from the map file `Anura-Family-Map-File.txt`. It will use 8 threads.

#### Additional arguments

##### `--polish`

> Performs additional distance measure refinement on final merged alignment, at cost of increasing run time (sometimes substantially).

##### `--verbose`

> Show mafft progress during final merge step instead of redirecting to log file. May be useful for monitoring the --polish option, which can be slow for larger datasets

##### `--output_seqnames`

> Write output file of sequence names and quit. See above section for usage.


### 3. Outputs <a name="OTP"></a>

There are three output directories created. Each contains files unique to particular steps of the analysis. An example of the output directory structure is shown below:

```
output_dir/
│
├── 1-Group-Fasta-Files/
│	├── Alytidae.multiples.fasta
│	├── Arthroleptidae.multiples.fasta
│	├── Ascaphidae.singleton.fasta
│	└── Batrachylidae.multiples.fasta
│
├── 2-Group-Alignments/
│	├── Alytidae.mafft-E-INS-i.fasta
│	├── Alytidae.mafft-E-INS-i.log
│	├── Arthroleptidae.mafft-E-INS-i.fasta
│	├── Arthroleptidae.mafft-E-INS-i.log
│	├── Ascaphidae.mafft-E-INS-i.fasta
│	├── Ascaphidae.mafft-E-INS-i.log
│	├── Batrachylidae.mafft-E-INS-i.fasta
│	└── Batrachylidae.mafft-E-INS-i.log
│
├── 3-Merged-Alignments/
│	├── Merged-alignment.fasta
│	├── Merge-Table.txt
│	├── Merged-alignment.log
│	└── Sub-Alignments.fasta
│
└── Taxonomy-Align_Aug-06-2022_12.16.14.log
```

`1-Group-Fasta-Files/`

This directory will contain fasta files for each grouping included in the map file. The `multiples` or `singleton` component of the name signify whether the file contains >= 2 sequences (`multiples`) or a single sequence (`singleton`). The `multiples` files are individually aligned in the next step. The `singleton` files only contain one sequence and they are aligned during the final sub-alignment merging step.


`2-Group-Alignments/`

This directory contains the aligned outputs of the group-specific fasta files. The associated log files show the MAFFT progress.


`3-Merged-Alignments/`

This directory contains the final outputs.

- `Merged-alignment.fasta`: This is the final alignment produced by TaxonomyAlign, and the main output file of interest.


- `Merge-Table.txt`: This MAFFT-specific file contains information about how the group merging was conducted. It is used in conjunction with the `Sub-Alignments.fasta` file.


- `Merged-alignment.log`: Contains all of the MAFFT alignment progress reports. If using `--verbose`, this will also appear on the screen. 


- `Sub-Alignments.fasta`: A concatenated sub-alignments fasta file that is specific to MAFFT. It is used in conjunction with the `Merge-Table.txt` file.


## License

GNU Lesser General Public License v3.0

## Contact and Bug Reporting

**TaxonomyAlign** is written and maintained by Daniel Portik. 

If you experience any problems running **TaxonomyAlign**, or have questions, please open an issue [here](https://github.com/dportik/TaxonomyAlign/issues).



























