#!/usr/bin/env python3

import sys
from os.path import dirname, realpath
from os import mkdir, remove
from shutil import rmtree, copy
path = dirname(realpath(__file__))
if path in sys.path:
    sys.path.remove(path)
sys.path.insert(0, realpath(path + '/../..'))

try:
    import graph_tool
except ModuleNotFoundError:
    print('\nCan\'t import the graph_tool python module! Make sure it\'s available in your environment\n')
    sys.exit(1)


from superpang.lib.Assembler import Assembler
from superpang.lib.utils import read_fasta, write_fasta, print_time
from mOTUlizer.classes.mOTU import mOTU
from mOTUlizer.utils import parse_checkm

from uuid import uuid4
from inspect import getfile
from hashlib import sha1
from collections import defaultdict
from subprocess import call, DEVNULL
from glob import glob
from argparse import ArgumentParser

libpath = dirname(realpath(getfile(Assembler)))

with open(libpath + '/../VERSION') as infile:
    VERSION = infile.read().strip()
CITATION = 'Puente-Sánchez F, Hoetzinger M, Buck M and Bertilsson S. Exploring intra-species diversity through non-redundant pangenome assemblies. Molecular Ecology Resources (2023) DOI: 10.1111/1755-0998.13826'

class ControlledExit(Exception):
    pass


def main(args, uuid):

    ### Welcome message
    print()
    print(f'This is SuperPang version {VERSION}')
    print()
    print('If using in publications or products please cite:')
    print()
    print(CITATION)
    print()
    if(args.debug):
        print('Running in debug mode. Execution speed will be slower...')
        print()
    print()


    ### File names
    temp_prefix          = f'{args.temp_dir}/{uuid}'
    corrected_dir        = args.output_dir + '/corrected_input'
    prelim_dir           = args.output_dir + '/preliminary'
    temp_saved_dir       = args.output_dir + '/temp'
    input_combined       = f'{temp_prefix}.combined.fasta'
    input_minimap2       = f'{temp_prefix}.pre.minimap2.fasta'
    diskdb               = f'{temp_prefix}.diskdb' if args.lowmem else None
    params               = args.output_dir + '/params.tsv'
    outputPre_kept       = prelim_dir      + '/prelim.fasta'
    outputPre2origs_kept = prelim_dir      + '/prelim2origins.tsv'
    name2bin_kept        = prelim_dir      + '/orig2bin.tsv'
    outputNodes          = args.output_dir + '/NBPs.fasta'
    outputCore           = args.output_dir + '/NBPs.core.fasta'
    outputAux            = args.output_dir + '/NBPs.accessory.fasta'
    outputNode2origs     = args.output_dir + '/NBP2origins.tsv'
    outputEdges          = args.output_dir + '/graph.fastg'
    outputName           = args.output_dir + '/assembly'


    ### Get sha1 of SuperPang scripts
    main_sha1 = sha1(open(path + '/' + 'SuperPang.py'     ).read().encode()).hexdigest()
    homo_sha1 = sha1(open(path + '/' + 'homogenize.py'    ).read().encode()).hexdigest()
    asse_sha1 = sha1(open(     getfile(Assembler)         ).read().encode()).hexdigest()
    cond_sha1 = sha1(open(path + '/' + 'condense-edges.py').read().encode()).hexdigest()
    

    ### Create output dirs
    try:
       mkdir(args.output_dir)
    except OSError as e:
        if e.errno != 17:
            raise
        elif args.force_overwrite:
            pass
        else:
            print(f'\nThe directory {args.output_dir} already contains results. Add --force-overwrite to ignore this message.\n')
            raise ControlledExit
    try:
        rmtree(corrected_dir)
    except:
        pass
    try:
        rmtree(prelim_dir)
    except:
        pass
    try:
        rmtree(temp_saved_dir)
    except:
        pass
    if args.keep_intermediate:
        mkdir(corrected_dir)
        mkdir(prelim_dir)
    if args.keep_temporary:
        mkdir(temp_saved_dir)


    ### Log params
    with open(params, 'w') as outfile:
        outfile.write(f'version\t{VERSION}\n'    )
        outfile.write(f'main_sha1\t{main_sha1}\n')
        outfile.write(f'homo_sha1\t{homo_sha1}\n')
        outfile.write(f'asse_sha1\t{asse_sha1}\n')
        outfile.write(f'cond_sha1\t{cond_sha1}\n')
        for arg in vars(args):
            outfile.write(f'{arg}\t{getattr(args, arg)}\n')
        outfile.write(f'uuid\t{uuid}\n')


    ### Load sequences
    name2bin = {}
    seqDict = {}
    genomes = open(args.fasta[0]).read().strip().split('\n') if len(args.fasta) == 1 else args.fasta
    for f in genomes:
        bin_ = f.rsplit('/',1)[1].rsplit('.',1)[0] if '/' in f else f.rsplit('.',1)[0]
        for name, seq in read_fasta(f).items():
            if name in seqDict:
                #raise Exception(f'Sequence "{name}" is duplicated in your input files')
                name = bin_ + '_' + name # assume the same name in two files is just a coincidence, and move on
            seqDict[name] = seq
            name2bin[name] = bin_
    write_fasta(seqDict, input_combined)
    if args.keep_intermediate:
        with open(name2bin_kept, 'w') as outfile:
            for name, bin_ in name2bin.items():
                outfile.write(f'{name}\t{bin_}\n')


    ### Ensure the checkm file is available and properly formatted, and parse it
    if not args.checkm:
        completeness = {bin_: args.default_completeness for bin_ in set(name2bin.values())}
    elif args.checkm.endswith('.tsv'):
        completeness = {bin_: float(completeness) for bin_, completeness in (line.strip().split('\t') for line in open(args.checkm))}
    else:
        completeness = {bin_: vals['Completeness'] for bin_, vals in parse_checkm(args.checkm).items()}
    missing_bins = set(name2bin.values()) - set(completeness)
    if missing_bins:
        print('\nThe following bins are missing from your CheckM/completeness file:\n')
        for bin_ in missing_bins:
            print(bin_)
        print('\nNote that SuperPang expect bin names in the CheckM/completeness file to NOT contain the file extension\n')
        raise ControlledExit()

    
    ### Correct input sequences with minimap2
    correct = args.identity_threshold and args.identity_threshold < 1
    if correct:
        # Check for minimap2
        try:
            call([args.minimap2_path, '-h'], stdout=DEVNULL, stderr=DEVNULL)
        except FileNotFoundError:
            print('\nCan\'t find minimap2. Make sure that it is present in your system and that the --minimap2-path is pointing to the executable\n')
            raise ControlledExit
        # Do stuff
        for i in range(1): # support for multiple calls to homogenize, but apparently it made no big difference
            if i:
                input_combined = input_minimap2
                input_minimap2 = f'{input_minimap2}.{i}'

            command = [sys.executable, path + '/' + 'homogenize.py', '-f', input_combined, '-o', input_minimap2, '-i', str(args.identity_threshold), '-m', str(args.mismatch_size_threshold),
                          '-g', str(args.indel_size_threshold), '-r', str(args.correction_repeats), '-n', str(args.correction_repeats_min), '-t', str(args.threads),
                          '--minimap2-path', args.minimap2_path, '--silent']
            if args.debug:
                command.append('--debug')
            ecode = call(command)

            if ecode:
                print('\nThere was an error running homogenize.py. Please open an issue\n')
                raise ControlledExit
                
        if args.keep_intermediate:
            outfiles = {bin_: open(f'{corrected_dir}/{bin_}.fasta', 'w') for bin_ in set(name2bin.values())}
            for name, seq in read_fasta(input_minimap2).items():
                outfiles[name2bin[name]].write(f'>{name}\n{seq}\n')
            for of in outfiles.values():
                of.close()
    else:
        input_minimap2 = input_combined          

    ### Assemble
    contigs = Assembler(input_minimap2, args.ksize, args.threads, diskdb, args.debug).run(args.minlen, args.mincov, args.bubble_identity_threshold, args.genome_assignment_threshold, args.threads, args.debug)
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
    print_time('Running mOTUpan')
    id2tag = {}
    if len(set(name2bin.values())) == 1: # only one bin, treat everything as core (we should use a different tag, print a warning and don't output the *.core.fasta
        for id_ in contigs:
            id2tag[id_] = 'core'
    else:
        # completeness dict for mOTUpan was defined above so that we can identify problems in the beginning
        name2ids = defaultdict(set)
        for id_, contig in contigs.items():
            for ori in contig.origins:
                name2ids[ori].add(id_)
        featDict = defaultdict(set)
        for name, bin_ in name2bin.items():
            featDict[bin_].update(name2ids[name])
        
        motu = mOTU( name = "mOTUpan_core_prediction" , faas = {} , cog_dict = featDict, checkm_dict = completeness, max_it = 100, threads = args.threads, precluster = False, method = 'default', quiet = not args.verbose_mOTUpan)
        if motu.get_stats()['mOTUpan_core_prediction']['core']:
            for id_ in motu.get_stats()['mOTUpan_core_prediction']['core']:
                id2tag[id_] = 'core'
            for id_ in motu.get_stats()['mOTUpan_core_prediction']['aux_genome']:
                id2tag[id_] = 'aux'
            for id_ in motu.get_stats()['mOTUpan_core_prediction']['singleton_cogs']:
                id2tag[id_] = 'singleton' # singletons are also aux, we overwrite them here
        else:
            print('mOTUpan was unable to predict the core genome for this dataset')
            for id_ in contigs:
                id2tag[id_] = 'noinfo'


    ### Header prefix
    header_prefix = '' if not args.header_prefix else f'{args.header_prefix}_'


    ### Prepare assembly graph nodes
    assemblyNodes = {}
    assemblyNodesCore = {}
    assemblyNodesAux = {}
    nodeNames = {}
    for id_, contig in contigs.items():
        nodeNames[id_] = f'{header_prefix}NODE_Sc{contig.scaffold}-{contig.i}-{id2tag[id_]}_length_{len(contig.tseq)}_cov_{round(contig.cov,2)}_tag_{id2tag[id_]};'
        if args.nice_headers:
            nodeNames[id_] = nodeNames[id_].replace('-', '_')
        if contig.tseq:
            assemblyNodes[nodeNames[id_]] = contig.tseq
            if id2tag[id_] == 'core':
                assemblyNodesCore[nodeNames[id_]] = contig.tseq
            elif id2tag[id_] == 'aux' or id2tag[id_] == 'singleton':
                assemblyNodesAux[nodeNames[id_]] = contig.tseq

    ### Prepare assembly graph edges
    edgeNames = {}
    for id_, contig in contigs.items():
        edgeNames[id_] = f'{header_prefix}EDGE_Sc{contig.scaffold}-{contig.i}-{id2tag[id_]}_length_{len(contig.seq)}_cov_{round(contig.cov,2)}_tag_{id2tag[id_]}'
        if args.nice_headers:
            edgeNames[id_] = edgeNames[id_].replace('-', '_')
    assemblyEdges = {}
    for id_, contig in contigs.items():
        succs = ','.join([edgeNames[succ] for succ in contig.successors])
        succs = f':{succs}' if succs else ''
        assemblyEdges[f'{edgeNames[id_]}{succs};'] = contig.seq

        
    ### Write final node, edge and core output
    write_fasta( assemblyNodes,     outputNodes )
    write_fasta( assemblyNodesCore, outputCore  )
    write_fasta( assemblyNodesAux,  outputAux   )
    write_fasta( assemblyEdges,     outputEdges )
    
    with open(outputNode2origs, 'w') as outfile:
        for id_, contig in contigs.items():
            origs = ','.join(contig.origins)
            bins = [name2bin[name] for name in contig.origins]
            lbins = len(set(bins))
            bins  = ','.join(bins)
            outfile.write(f'{nodeNames[id_]}\t{origs}\t{bins}\t{lbins}/{len(completeness)}\n')

    ### Condense edges
    print_time('Reconstructing contigs')
    command = [sys.executable, path + '/' + 'condense-edges.py', outputEdges, outputName, str(args.ksize)]
    if args.header_prefix:
        command.extend(['--header-prefix', args.header_prefix])
    ecode = call(command)
    if ecode:
        print('\nThere was an error running condense-edges.py. Please open an issue\n')
        raise ControlledExit

    if args.keep_temporary:
        for f in glob(f'{temp_prefix}*'):
            copy(f, temp_saved_dir)
    
    print_time('Finished')


def parse_args():
    parser = ArgumentParser(description='Create a consensus pangenome assembly from a set of bins/genomes from the same mOTU/species')
    parser.add_argument('-f', '--fasta', type = str, nargs='+', required = True,
                        help = 'Input fasta files with the sequences for each bin/genome, or a single file containing the path to one input fasta file per line')
    parser.add_argument('-q', '--checkm', type = str,
                        help = 'CheckM output for the bins, or *.tsv file with bin and completeness for each bin')
    parser.add_argument('-i', '--identity_threshold', type = float, default = 0.95,
                        help = 'Identity threshold (fraction) to initiate correction with minimap2')
    parser.add_argument('-m', '--mismatch-size-threshold', type = int, default = 100,
                        help = 'Maximum contiguous mismatch size that will be corrected')
    parser.add_argument('-g', '--indel-size-threshold', type = int, default = 100,
                        help = 'Maximum contiguous indel size that will be corrected')
    parser.add_argument('-r', '--correction-repeats', type = int, default = 20,
                        help = 'Maximum iterations for sequence correction')
    parser.add_argument('-n', '--correction-repeats-min', type = int, default = 5,
                        help = 'Minimum iterations for sequence correction')
    parser.add_argument('-k', '--ksize', type = int, default = 301,
                        help = 'Kmer size')
    parser.add_argument('-l', '--minlen', type = int, default = 0,
                        help = 'Scaffold length cutoff')
    parser.add_argument('-c', '--mincov', type = float, default = 0,
                        help = 'Scaffold coverage cutoff')
    parser.add_argument('-b', '--bubble-identity-threshold', type = float, default = 0.95,
                        help = 'Minimum identity (matches / alignment length) required to remove a bubble in the sequence graph')
    parser.add_argument('-a', '--genome-assignment-threshold', default = 0.5, type = float,
                        help = 'Fraction of shared kmers required to assign a contig to an input genome')
    parser.add_argument('-x', '--default-completeness', type = float, default = 70,
                        help = 'Default genome completeness to assume if a CheckM output is not provided')
    parser.add_argument('-t', '--threads', type = int, default = 1,
                        help = 'Number of processors to use')
    parser.add_argument('-o', '--output-dir', type = str, default = 'output',
                        help = 'Output directory')
    parser.add_argument('-d', '--temp-dir', type = str, default = '/tmp/',
                        help = 'Temp directory')
    parser.add_argument('-u', '--header-prefix', type = str,
                        help = 'Prefix to be added to output sequence names')
    parser.add_argument('--assume-complete', action='store_true',
                        help = 'Assume that the input genomes are complete (--genome-assignment-threshold 0.95 --default-completeness 99)')
    parser.add_argument('--lowmem', action='store_true',
                        help = 'Use disk storages instead of memory when possible, reduces memory usage at the cost of execution time')
    parser.add_argument('--minimap2-path', type = str, default = 'minimap2',
                        help = 'Path to the minimap2 executable')
    parser.add_argument('--keep-intermediate', action='store_true',
                        help = 'Keep intermediate files')
    parser.add_argument('--keep-temporary', action='store_true',
                        help = 'Keep temporary files')
    parser.add_argument('--verbose-mOTUpan', action='store_true',
                        help = 'Print out mOTUpan logs')
    parser.add_argument('--nice-headers', action='store_true',
                        help = 'Replace dashes with underscores in output sequence names')
    parser.add_argument('--force-overwrite', action='store_true',
                        help='Write results even if the output directory already exists')
    parser.add_argument('--debug', action='store_true',
                        help = 'Run additional sanity checks (increases execution time)')
    args = parser.parse_args()
    if args.assume_complete:
        args.checkm = None
        args.genome_assignment_threshold = 0.95
        args.default_completeness = 99
    if args.debug:
        args.verbose_mOTUpan = True
    return args


def cli():
    args        = parse_args()
    uuid        = uuid4().hex[:7]
    temp_prefix = f'{args.temp_dir}/{uuid}'
    try:
        main(args, uuid)
    except ControlledExit:
        sys.exit(1)
    finally:
        if not args.keep_temporary:
            call(['rm', '-r'] + glob(f'{temp_prefix}*'), stdout=DEVNULL, stderr=DEVNULL)

 
if __name__ == '__main__':
    cli()
