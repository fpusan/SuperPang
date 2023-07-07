# SuperPang: non-redundant pangenome assemblies from multiple genomes or bins

**Check our paper:** Puente-Sánchez F, Hoetzinger M, Buck M and Bertilsson S. [Exploring intra-species diversity through non-redundant pangenome assemblies](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13826) _Molecular Ecology Resources_ (2023) DOI: 10.1111/1755-0998.13826

_... but note that performance is now better (3x less memory usage, 20% faster execution time) than when we first benchmarked Superpang._

## Installation
Requires [graph-tool](https://graph-tool.skewed.de/), [speedict](https://github.com/Congyuwang/RocksDict), [mOTUlizer v0.2.4](https://github.com/moritzbuck/mOTUlizer), [minimap2](https://github.com/lh3/minimap2) and [mappy](https://pypi.org/project/mappy/). The easiest way to get it running is using conda.
```
# Install into a new conda environment
conda create -n SuperPang -c conda-forge -c bioconda -c fpusan superpang
# Check that it works for you!
conda activate SuperPang
test-SuperPang.py
```

## Usage
`SuperPang.py --fasta <genome1.fasta> <genome2.fasta> <genomeN.fasta> --checkm <check_results> --output-dir <output_directory>`


## Input files and choice of parameters
- The input genomes can be genomes from isolates, MAGs (Metagenome-Assembled Genomes) or SAGs (Single-cell Assembled Genomes).
- The input genomes can have different qualities, for normal usage we recommend that you provide completeness estimates for each input genome through the `-q/--checkm` parameter.
- If you are certain that all your input genomes are complete, you can use the `--assume-complete` flag or manually tweak the `-a/--genome-assignment-threshold` and `-x/--default-completeness` parameters instead of providing a file with completeness estimates.
- The default parameter values in SuperPang assume that all of the input genomes come from the same species (ANI>=0.95). This can be controlled by changing the values of the `-i/--identity_threshold` and `-b/--bubble-identity-threshold` to the expected ANI. However SuperPang has currently only been tested in species-level clusters.


## Arguments
* *-f/--fasta*: Input fasta files with the sequences for each bin/genome, or a single file containing the path to one input fasta file per line.
* *-q/--checkm*: CheckM output for the bins. This can be the STDOUT of running checkm on all the fasta files passed in *--fasta*, or a tab-delimited file ended with a `.tsv` extension, in the form `genome1 percent_completeness`. Genome names should not contain the file extension (e.g. `.fna`). If empty, completeness will be estimated by [mOTUpan](https://www.biorxiv.org/content/10.1101/2021.06.25.449606v1) but this may lead to wrong estimations for very incomplete genomes.
* *-i/--identity_threshold*: Identity threshold (fraction) to initiate correction with minimap2. Values of 1 or higher will skip the correction step entirely. Default `0.95`.
* *-m/--mismatch-size-threshold*: Maximum contiguous mismatch size that will be corrected. Default `100`.
* *-g/--indel-size-threshold*: Maximum contiguous indel size that will be corrected. Default `100`.
* *-r/--correction-repeats*: Maximum iterations for sequence correction. Default `20`.
* *-n/--correction-repeats-min*: Minimum iterations for sequence correction. Default `5`.
* *-k/--ksize*: Kmer-size. Default `301`.
* *-l/--minlen*: Scaffold length cutoff. Default `0` (no cutoff).
* *-c/--mincov*: Scaffold coverage cutoff. Default `0` (no cutoff).
* *-b/--bubble-identity-threshold*: Minimum identity (matches / length) required to remove a bubble in the sequence graph. Default `0.95`.
* *-a/--genome-assignment-threshold*. Fraction of shared kmers required to assign a contig to an input genome (0 means a single shared kmer is enough). Default `0.5`.
* *-x/--default-completeness*: Default genome completeness to assume if a CheckM output is not provided with *--checkm*. Default `70`.
* *-t/--threads*: Number of processors to use. Default `1`.
* *-o/--output*: Output directory. Default `output`.
* *-d/--temp-dir*: Directory for temp files. Default `tmp`.
* *--assume-complete*: Assume that the input genomes are complete (*--genome-assignment-threshold 0.95*, *--default-completeness 99*).
* *--lowmem*: Use disk storages instead of memory when possible, reduces memory usage at the cost of execution time.
* *--minimap2-path*: Path to the minimap2 executable. Default `minimap2`.
* *--keep-intermediate*: Keep intermediate files.
* *--keep-temporary*: Keep temporary files.
* *--verbose-mOTUpan*: Print out mOTUpan logs.
* *--debug*: Run additional sanity checks (increases execution time).

## Output
* `assembly.fasta`: contigs.
* `assembly.info`: core/auxiliary and path information for each contig.
* `NBPs.fasta`: non-branching paths.
* `NBPs.core.fasta`: non-branching paths deemed to belong to the core genome of the species by [mOTUpan](https://www.biorxiv.org/content/10.1101/2021.06.25.449606v1).
* `NBPs.accessory.fasta`: non-branching paths deemed to belong to the accessory genome of the species.
* `graph.fastg`: assembly graph in a format compatible with [bandage](https://rrwick.github.io/Bandage/).
* `NBP2origins.tsv`: tab-separated file with the non-branching path IDs, a comma-separated list of the input sequences in which that NBP was deemed present, a comma-separated list of the input genomes in which that NBP was deemed present, and the number of input genomes in which that NBP was deemed present.
* `params.tsv`: parameters used in the run.

## About
*SuperPang* is developed by Fernando Puente-Sánchez (Sveriges lantbruksuniversitet). Feel free to open an issue or reach out for support [fernando.puente.sanchez@slu.se](mailto:fernando.puente.sanchez@slu.se).
