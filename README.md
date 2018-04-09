APPAS
# Rapid Alignment-free Phylogenetic PLacement via Ancestral Sequences

## Description

RAPPAS (Rapid Alignment-free Phylogenetic PLacement via Ancestral Sequences) is a software dedicated to "Phylogenetic Placement" (PP) of metagenomic reads on a reference tree. Compared to previous PP programs, RAPPAS uses a 2 step approach divided in a) database build, b) placement itself.

The main advantage of RAPPAS is that it is alignment free, which means after the first step (DB build) is performed, metagenomic reads can be directly placed on a referene tree WHITHOUT the step of aligning the reads to the reference alignment from which was built the tree (required by other approaches).

The second advantage of RAPPAS is its algorithm based on ancestal k-mer matches, which makes its execution time linear to the length of the placed sequences.

## Installation

### Prerequisite

- RAPPAS compilation requires a clean JDK 1.8 javac compiler installation. Java >=1.8 is a compulsory requirement as some operations are based on Lambda expressions.
- Apach Ant is used to facilitate the compilation.

Below are  instructions for a Debian based Unix distribution. For compiling Java sources with Apache Ant on other Operating Systems, please follow these instructions.

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

Installation of Apach Ant:
```
sudo apt-get install ant
```

### Source download and compilation

```
#download git repository
git clone -b master https://gite.lirmm.fr/linard/RAPPAS.git
#compile
ant -f build-cli.xml
```
The RAPPAS jar can then be found in the ./dist directory.




## Usage

### Reference Dataset
First, one has to prepare a reference dataset designed to answer its biological question. Typically, in the context of metagenomics and taxonomic identifications, a marker gene (16S rRNA, cox1, rbcl...) is used to build a reference species tree. This species tree is the basis for phylogenetic placement of marker gene(s).
For RAPPAS, the reference dataset is composed of:
1. A reference alignment of all sequences of this marker gene
2. A phylogenetic tree inferred from this reference alignment

Such reference marker gene datasets can be found, for instance, via:
- "The All-Species Living Tree" Project (LTP, eukaryote rRNAs) :  <https://www.arb-silva.de/projects/living-tree/>
- Greengenes (bacterial 16S) : <http://greengenes.secondgenome.com/>
- The curated database of Eukref : <http://eukref.org/databases/>
- Or built internally in the lab.

### RAPPAS Database build 

The basic command is 

```
java -jar RAPPAS.jar -m b -s [nucl|prot] -b ARbinary -w workdir -r reference_alignment.fasta -t reference_tree.newick
```

where

option | expected value | is used for
--- | --- | ---
`-s (--states)` | "nucl" or "prot" | Set if we plan a nucleotide or protein analysis.
`-b (--arbinary)` | an integer >=3 | Set the path to the AR binary used for ancestral sequence reconstruction (see note below).
`-w (--workdir)` | a directory | Set the directory in which will be saved the database.
`-r (--refalign)` | a file | The reference alignment, in fasta format.
`-t (--reftree)̀` | a file | The reference tree, in newick format.

Currently, PhyML or PAML are the binary called for ancestral sequeunce reconstruction and fully supported by RAPPAS.
You can use the versions used on the authors websites, but we recommand the following stable and hacked versions which are faster and require less disk space:
[link to hacked PhyML]() (faster, strongly recommended)
[link to hacked PAML]()  (require less RAM, but much slower)

--> The reconstruction will result to the production of a directory struture and a database file in the "workdir":

file/directory | description
--- | --- 
`*.union` | The RAPPAS database itself.
`workdir/extended_tree` | Temp files used at DB construction, allow to explore the phantom nodes.
`workdir/AR` | Temp files used at DB construction, the raw outputs of PhyML or PAML.
`workdir/logs` | As the name says.


### Queries placement

Once the database is built, a placement command can be called numerous times on different query sequence datasets. v1.00 of RAPPAS places around 1,000,000 metagenomic reads in 40 minutes using a single core of a normal desktop PC. Only the database file and the query sequences are used as parameters.


```
java -jar RAPPAS.jar -m p -d database.union -q queries.fasta 
```

where

option | expected value | is used for
--- | --- | ---
`-d (--database)` | a file | the *.union file created at previous DB build step.
`-q (--queries)` | a file | The query reads, in fasta format.

The *.jplace describing the placements will be in the workdir/logs directory.

To know more about the [jplace format](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031009).
To know more about the [exploitation of phylogenetic placement results](https://matsen.github.io/pplacer/generated_rst/guppy.html#introduction) (OTU alpha diversity, Unifrac-like measures...).


## License

RAPPAS is available under the GPLv3 license.


