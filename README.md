# SuperPang: non-redundant pangenome assemblies from individual genomes or bins

## Installation
Requires [graph-tool](https://graph-tool.skewed.de/) and [mOTUlizer](https://github.com/moritzbuck/mOTUlizer). Both can be easily installed from conda. I will make a proper bioconda package when the software becomes a bit more stable.

## Usage
`SuperPang.py --fasta <genome1.fasta> <genome2.fasta> <genomeN.fasta> --checkm <check_results> --output-dir <output_directory>`

**Arguments**

* *-f/--fasta*: Input fasta files with the sequences for each bin
* *-q/--checkm*: CheckM output for the bins. This can be the STDOUT of running checkm on all the fasta files passed in *--fasta*, or a tab-delimited file in the form `genome1 percent_completeness`. If empty, completeness will be estimated by [mOTUpan](https://www.biorxiv.org/content/10.1101/2021.06.25.449606v1) but this may lead to wrong estimations for very incomplete genomes.
* *-i/--identity_threshold*: Identity threshold (fraction) to initiate correction with minimap2. Default `0.9`
* *-m/--mismatch-size-threshold*: Maximum contiguous mismatch size that will be corrected. Default `100`
* *-g/--indel-size-threshold*: Maximum contiguous indel size that will be corrected. Default `100`
* *-r/--correction-repeats*: Maximum iterations for sequence correction. Default `1`
* *-k/--ksize*: Kmer-size. Default `101`
* *-l/--minlen*: Scaffold length cutoff. Default `0` (no cutoff)
* *-c/--mincov*: Scaffold coverage cutoff. Default `0` (no cutoff)
* *-a/--genome-assignment-threshold*. Fraction of shared kmers required to assign a contig to an input genome. Default `0` (a single shared kmer is enough)
* *-x/--default-completeness*: Default genome completeness to assume if a CheckM output is not provided with *--checkm*. Default `50`
* *-t/--threads*: Number of processors to use. Default `1`
* *-o/--output*: Output directory. Default `output`
* *--assume-complete*: Assume that the input genomes are complete (*--genome-assignment-threshold 1*, *--default-completeness 100*)
* *--minimap2-path*: Path to the minimap2 executable. Default `minimap2`
* *--keep-intermediate*: Keep intermediate files

**Output**
* `assembly.fasta`: contigs
* `nodes.fasta`: assembly nodes
* `core.fasta`: assembly nodes deemed to belong to the core genome of the species by [mOTUpan](https://www.biorxiv.org/content/10.1101/2021.06.25.449606v1)
* `graph.fastg`: assembly graph in a format compatible with [bandage](https://rrwick.github.io/Bandage/)
* `node2origins.tsv`: tab-separated file with the assembly nodes, and a comma-separated list of the input genome in which that node was deemed present

## About
*SuperPang* is developed by Fernando Puente-Sánchez (Sveriges lantsbruksuniversitet). Feel free to open an issue or reach out for support [fernando.puente.sanchez@slu.se](mailto:fernando.puente.sanchez@slu.se)).
