RAPPAS
# Rapid Alignment-free Phylogenetic PLacement via Ancestral Sequences

## Description

RAPPAS (Rapid Alignment-free Phylogenetic PLacement via Ancestral Sequences) is a program dedicated to "Phylogenetic Placement" (PP) of metagenomic reads on a reference tree. As apposed to previous PP programs, RAPPAS uses a 2 step approach divided into a) the database build, and b) the placement itself.

The main advantage of RAPPAS is that it is alignment free, which means that after step (a) (the DB build) is performed, metagenomic reads can be directly placed on a referene tree _WHITHOUT_ aligning them to the reference alignment on which the tree was built (as required by other approaches).

The second advantage of RAPPAS is its algorithm based on ancestal k-mer matches, making its execution time linear with respect to the length of the placed sequences.

## Installation

### Prerequisites

- RAPPAS compilation requires a clean JDK 1.8 javac compiler installation. Java >=1.8 is a compulsory requirement as some operations are based on Lambda expressions.
- Apache Ant is used to facilitate the compilation.

We provide instructions for Debian-based Linux distributions. For compiling Java sources with Apache Ant on other operating systems, please perform analogous operations on your system.

Using OpenJDK 1.8:
```
#install packages
sudo apt-get update
sudo apt-get install openjdk-8-jdk
#update relevant symlinks to make v1.8 default
sudo update-java-alternatives --set java-1.8.0-openjdk-amd64

```
Using the proprietary Oracle JDK 1.8:
```
#install packages
sudo add-apt-repository ppa:webupd8team/java
sudo apt-get update
sudo apt-get install oracle-java8-installer
#update relevant symlinks to make v1.8 default
sudo apt-get install oracle-java8-set-default
```

Installation of Apache Ant:
```
sudo apt-get install ant
```

### Source download and compilation

```
#download git repository
git clone -b master https://gite.lirmm.fr/linard/RAPPAS.git
#compile
cd RAPPAS && ant -f build-cli.xml
```
The executable RAPPAS.jar can then be found in the ./dist directory.




## Usage

### Reference Dataset
First, one has to prepare a reference dataset designed to answer a biological question. Typically, in the context of metagenomics and taxonomic identifications, a marker gene (16S rRNA, cox1, rbcl...) is used to build a reference species tree. This species tree is the basis for phylogenetic placement of marker gene(s).
For RAPPAS, the reference dataset is composed of:
1. A reference alignment of all sequences of this marker gene
2. A phylogenetic tree inferred from this reference alignment

Such reference marker gene datasets can be found, for instance, via:
- "The All-Species Living Tree" Project (LTP, eukaryote rRNAs) :  <https://www.arb-silva.de/projects/living-tree/>,
- Greengenes (bacterial 16S) : <http://greengenes.secondgenome.com/>,
- The curated database of Eukref : <http://eukref.org/databases/>,
- Or built internally in the lab.

### RAPPAS database build 

__Basic command__

```
java -jar RAPPAS.jar -m b -s [nucl|prot] -b ARbinary -w workdir -r reference_alignment.fasta -t reference_tree.newick
```

where

option | expected value | description
--- | --- | ---
`-s (--states)` | "nucl" or "prot" | Set if we use a nucleotide or protein analysis.
`-b (--arbinary)` | a binary of PhyML or PAML | Set the path to the binary used for ancestral sequence reconstruction (see note below).
`-w (--workdir)` | a directory | Set the directory to save the database in.
`-r (--refalign)` | a file | The reference alignment, in fasta format.
`-t (--reftree)̀` | a file | The reference tree, in newick format.

__Note on PhyML and PAML binaries__:
Currently, the following programs are fully supported by RAPPAS for generating ancestral sequence posterior probabilities:
- PhyML : Fastest & strongly recommended but may require lots of RAM.
- PAML  : Slower,  but requires less memory.

You can use the latest versions provided on the authors' websites, but we recommand the _HACKED VERSIONS_ available in this git repository in the /depbin directory.
These are based on slightly modified sources of PhyML and PAML: no change in ML computations, but useless output is skipped, making the reconstruction process faster.

The reconstruction will result in the production of a directory structure and a database file in the given "workdir":

file/directory | description
--- | --- 
`*.union` | The RAPPAS database itself.
`workdir/extended_tree` | Temporary files used at DB construction, allowing the exploration of the phantom nodes.
`workdir/AR` | Temporary files used at DB construction, the raw output of PhyML or PAML.
`workdir/logs` | As the name says.


### Query placement

After building the RAPPAS DB, placement commands can be called numerous times on different query sequence datasets.
v1.00 of RAPPAS places 1,000,000 metagenomic of 150bp in ~40 minutes, using only a single core of a normal desktop PC.

```
java -jar RAPPAS.jar -m p -d database.union -q queries.fasta 
```

where

option | expected value | description
--- | --- | ---
`-d (--database)` | a file | the *.union file created at previous DB build step.
`-q (--queries)` | a file | The query reads, in fasta format.

The *.jplace describing the placements of all queries will be written in the ./workdir/logs directory.

__To know more about :__
- the [jplace format](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009).
- the [exploitation of phylogenetic placement results](https://matsen.github.io/pplacer/generated_rst/guppy.html#introduction) (OTU alpha diversity, Unifrac-like measures...).

### Other options

__Normal options:__

option | expected value | description
--- | --- | ---
`-k` | integer >=3 | The k-mer length used at DB build (default=8)


__Debug options:__

Avoid debug options if you are not involved in RAPPAS developpement !!!
Description coming soon...

## License

RAPPAS is available under the GPLv3 license.


