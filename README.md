# SuperPang: non-redundant pangenome assemblies from multiple genomes or bins

## Installation
Requires [graph-tool](https://graph-tool.skewed.de/), [mOTUlizer v0.2.4](https://github.com/moritzbuck/mOTUlizer), [minimap2](https://github.com/lh3/minimap2) and [mappy](https://pypi.org/project/mappy/). The easiest way to get it running is using conda.
```
# Install into a new conda environment
conda create -n SuperPang -c conda-forge -c bioconda -c fpusan superpang
# Check that it works for you!
conda activate SuperPang
test-SuperPang.py
```

## Usage
`SuperPang.py --fasta <genome1.fasta> <genome2.fasta> <genomeN.fasta> --checkm <check_results> --output-dir <output_directory>`

**Arguments**

* *-f/--fasta*: Input fasta files with the sequences for each bin/genome
* *-q/--checkm*: CheckM output for the bins. This can be the STDOUT of running checkm on all the fasta files passed in *--fasta*, or a tab-delimited file in the form `genome1 percent_completeness`. If empty, completeness will be estimated by [mOTUpan](https://www.biorxiv.org/content/10.1101/2021.06.25.449606v1) but this may lead to wrong estimations for very incomplete genomes.
* *-i/--identity_threshold*: Identity threshold (fraction) to initiate correction with minimap2. Default `0.9`.
* *-m/--mismatch-size-threshold*: Maximum contiguous mismatch size that will be corrected. Default `100`.
* *-g/--indel-size-threshold*: Maximum contiguous indel size that will be corrected. Default `100`.
* *-r/--correction-repeats*: Maximum iterations for sequence correction. Default `5`.
* *-n/--correction-repeats-min*: Minimum iterations for sequence correction. Default `5`.
* *-k/--ksize*: Kmer-size. Default `301`.
* *-l/--minlen*: Scaffold length cutoff. Default `0` (no cutoff).
* *-c/--mincov*: Scaffold coverage cutoff. Default `0` (no cutoff).
* *-b/--bubble-identity-threshold*: Minimum identity (matches / alignment length) required to remove a bubble in the sequence graph.
* *-a/--genome-assignment-threshold*. Fraction of shared kmers required to assign a contig to an input genome (0 means a shared kmer is enough). Default `0.5`.
* *-x/--default-completeness*: Default genome completeness to assume if a CheckM output is not provided with *--checkm*. Default `50`.
* *-t/--threads*: Number of processors to use. Default `1`.
* *-o/--output*: Output directory. Default `output`.
* *--assume-complete*: Assume that the input genomes are complete (*--genome-assignment-threshold 0.95*, *--default-completeness 95*).
* *--minimap2-path*: Path to the minimap2 executable. Default `minimap2`.
* *--keep-intermediate*: Keep intermediate files.

**Output**
* `assembly.fasta`: contigs.
* `assembly.info`: core/auxiliary and path information for each contig.
* `nodes.fasta`: assembly nodes.
* `core.fasta`: assembly nodes deemed to belong to the core genome of the species by [mOTUpan](https://www.biorxiv.org/content/10.1101/2021.06.25.449606v1).
* `auxiliary.fasta`: assembly nodes deemed to belong to the auxiliary genome of the species.
* `graph.fastg`: assembly graph in a format compatible with [bandage](https://rrwick.github.io/Bandage/).
* `node2origins.tsv`: tab-separated file with the assembly nodes, and a comma-separated list of the input genome in which that node was deemed present.
* `params.tsv`: parameters used in the run.

## About
*SuperPang* is developed by Fernando Puente-Sánchez (Sveriges lantsbruksuniversitet). Feel free to open an issue or reach out for support [fernando.puente.sanchez@slu.se](mailto:fernando.puente.sanchez@slu.se).
