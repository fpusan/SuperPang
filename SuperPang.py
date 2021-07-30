#!/usr/bin/env python3

from lib.Assembler import Assembler
from lib.utils import read_fasta, write_fasta
from mOTUlizer.classes.mOTU import mOTU
from mOTUlizer.utils import parse_checkm

from uuid import uuid4
from collections import defaultdict
from subprocess import call
from argparse import ArgumentParser

from sys import argv
from os.path import dirname, realpath
path = dirname(realpath(argv[0]))

def main(args):

    uuid = uuid4().hex[:7]
    input_combined       = f'/tmp/{uuid}.combined.fasta'
    input_minimap2       = f'/tmp/{uuid}.pre.minimap2.fasta'
    input_minimap2_kept  = args.output.rsplit('.', 1)[0] + '.unassembled.corrected.fasta'
    outputPre_kept       = args.output.rsplit('.', 1)[0] + '.prelim.fasta'
    outputPre2origs_kept = args.output.rsplit('.', 1)[0] + '.prelim2origins.tsv'
    name2bin_kept        = args.output.rsplit('.', 1)[0] + '.orig2bin.tsv'
    node2origs_kept      = args.output.rsplit('.', 1)[0] + '.node2origins.tsv'

    ### Load sequences
    name2bin = {}
    seqDict = {}
    for f in args.fasta:
        bin_ = f.rsplit('/',1)[1].rsplit('.',1)[0]
        for name, seq in read_fasta(f).items():
            assert name not in seqDict
            seqDict[name] = seq
            name2bin[name] = bin_
    write_fasta(seqDict, input_combined)
    if args.keep_intermediate:
        with open(name2bin_kept, 'w') as outfile:
            for name, bin_ in name2bin.items():
                outfile.write(f'{name}\t{bin_}\n')
    
    ### Correct input sequences with minimap2
    if args.identity_threshold and args.identity_threshold < 1:
        call([path + '/' + 'run_minimap2.py', '-f', input_combined, '-o', input_minimap2, '-i', str(args.identity_threshold), '-m', str(args.mismatch_size_threshold),
              '-g', str(args.indel_size_threshold), '-r', str(args.correction_repeats), '-t', str(args.threads), '--minimap2-path', args.minimap2_path, '--silent'])
        if args.keep_intermediate:
            call(['cp', input_minimap2, input_minimap2_kept])
            outfiles = {bin_: open(f'{bin_}.corrected.fasta', 'w') for bin_ in set(name2bin.values())} ############## check out dir
            for name, seq in read_fasta(input_minimap2_kept).items():
                outfiles[name2bin[name]].write(f'>{name}\n{seq}\n')
            for of in outfiles.values():
                of.close()
    else:
        input_minimap2 = input_combined           

    ### Assemble
    contigs = Assembler(input_minimap2, args.ksize).run(args.minlen, args.mincov, args.genome_assignment_threshold)
    if args.keep_intermediate:
        prelim = {}
        with open(outputPre2origs_kept, 'w') as outfile:
            for id_, contig in contigs.items():
                origs = ','.join(contig.origins)
                prelimName = f'PRELIM_Sc{contig.scaffold}-{contig.i}_length_{len(contig.seq)}_cov_{round(contig.cov,2)};'
                outfile.write(f'{prelimName}\t{origs}\n')
                prelim[prelimName] = seq
        write_fasta(prelim, outputPre_kept)


    ### Tag contigs with mOTUpan
    id2tag = {}
    if len(set(name2bin.values())) == 1: # only one bin, treat everything as core (we should use a different tag, print a warning and don't output the *.core.fasta
        for id_ in contigs:
            id2tag[id_] = 'core'
    else:
        completeness = {bin_: args.default_completeness for bin_ in set(name2bin.values())} if not args.checkm else {bin_: vals['Completeness'] for bin_, vals in parse_checkm(args.checkm).items()}
        name2ids = defaultdict(set)
        for id_, contig in contigs.items():
            for ori in contig.origins:
                name2ids[ori].add(id_)
        featDict = defaultdict(set)
        for name, bin_ in name2bin.items():
            featDict[bin_].update(name2ids[name])
        
        motu = mOTU( name = "mOTUpan_core_prediction" , faas = {} , cog_dict = featDict, checkm_dict = completeness, max_it = 100, threads = args.threads, precluster = False, method = 'default')
        for id_ in motu.get_stats()['mOTUpan_core_prediction']['core']:
            id2tag[id_] = 'core'
        for id_ in motu.get_stats()['mOTUpan_core_prediction']['aux_genome']:
            id2tag[id_] = 'aux'
        for id_ in motu.get_stats()['mOTUpan_core_prediction']['singleton_cogs']:
            id2tag[id_] = 'singleton' # singletons are also aux, we overwrite them here
       
    
    ### Prepare assembly graph nodes
    assemblyNodes = {}
    assemblyNodesCore = {}
    nodeNames = {}
    for id_, contig in contigs.items():
        nodeNames[id_] = f'NODE_Sc{contig.scaffold}-{contig.i}-{id2tag[id_]}_length_{len(seq)}_cov_{round(contig.cov,2)}_tag_{id2tag[id_]};'
        if contig.tseq:
            assemblyNodes[nodeNames[id_]] = contig.tseq
            if id2tag[id_] == 'core':
               assemblyNodesCore[nodeNames[id_]] = contig.tseq

    ### Prepare assembly graph edges
    edgeNames = {}
    for id_, contig in contigs.items():
        edgeNames[id_] = f'EDGE_Sc{contig.scaffold}-{contig.i}-{id2tag[id_]}_length_{len(contig.seq)}_cov_{round(contig.cov,2)}_tag_{id2tag[id_]}'
    assemblyEdges = {}
    for id_, contig in contigs.items():
        succs = ','.join([edgeNames[succ] for succ in contig.successors])
        succs = f':{succs}' if succs else ''
        assemblyEdges[f'{edgeNames[id_]}{succs};'] = contig.seq

    ### Write final fasta output
    write_fasta(assemblyNodes, args.output)
    write_fasta(assemblyNodesCore, args.output.rsplit('.', 1)[0] + '.core.fasta')
    write_fasta(assemblyEdges, args.output.rsplit('.', 1)[0] + '.fastg')
    if args.keep_intermediate:
        with open(node2origs_kept, 'w') as outfile:
            for id_, contig in contigs.items():
                origs = ','.join(contig.origins)
                outfile.write(f'{nodeNames[id_]}\t{origs}\n')


def parse_args():
    parser = ArgumentParser(description='Create a consensus pangenome assembly from a set of bins from the same mOTU')
    parser.add_argument('-f', '--fasta', type = str, nargs='+',
                        help = 'Input fasta files with the sequences for each bin')
    parser.add_argument('-q', '--checkm', type = str,
                        help = 'CheckM output for the bins')
    parser.add_argument('-i', '--identity_threshold', type = float, default = 0.7,
                        help = 'Identity threshold (fraction) to initiate correction with minimap2')
    parser.add_argument('-m', '--mismatch-size-threshold', type = int, default = 100,
                        help = 'Maximum contiguous mismatch size that will be corrected')
    parser.add_argument('-g', '--indel-size-threshold', type = int, default = 100,
                        help = 'Maximum contiguous indel size that will be corrected')
    parser.add_argument('-r', '--correction-repeats', type = int, default = 1,
                        help = 'Maximum iterations for sequence correction')
    parser.add_argument('-k', '--ksize', type = int, default = 101,
                        help = 'Kmer size')
    parser.add_argument('-l', '--minlen', type = int, default = 0,
                        help = 'Scaffold length cutoff')
    parser.add_argument('-c', '--mincov', type = float, default = 0,
                        help = 'Scaffold coverage cutoff')
    parser.add_argument('-a', '--genome-assignment-threshold', default = 0, type = float,
                        help = 'Fraction of shared kmers required to assign a contig to an input genome')
    parser.add_argument('-x', '--default-completeness', type = float, default = 50,
                        help = 'Default genome completeness to assume if a CheckM output is not provided')
    parser.add_argument('-t', '--threads', type = int, default = 1,
                        help = 'Number of processors to use')
    parser.add_argument('-o', '--output', type = str, default = 'assembly.fasta',
                        help = 'Output name')
    parser.add_argument('--assume-complete', action='store_true',
                        help = 'Assume that the input genomes are complete (--genome-assignment-threshold 1 --default-completeness 100)')
    parser.add_argument('--minimap2-path', type = str, default = 'minimap2',
                        help = 'Path to the minimap2 executable')
    parser.add_argument('--keep-intermediate', action='store_true',
                        help = 'Keep intermediate files')
    args = parser.parse_args()
    if args.assume_complete:
        args.checkm = None
        args.genome_assignment_threshold = 1
        args.default_completeness = 100
    return args


if __name__ == '__main__':
    main(parse_args())

    
    
