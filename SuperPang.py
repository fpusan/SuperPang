#!/usr/bin/env python3

from lib.Assembler import Assembler

from subprocess import call
from argparse import ArgumentParser
from uuid import uuid4

from sys import argv
from os.path import dirname, realpath
path = dirname(realpath(argv[0]))

def main(args):

    uuid = uuid4().hex[:7]
    input_pre = f'/tmp/{uuid}.pre'
    input_minimap2 = f'/tmp/{uuid}.pre.minimap2.fasta'
    output_ass = f'/tmp/{uuid}.ass'
    output_ass_fastg =f'/tmp/{uuid}.fastg'
    output_final_fastg = args.output.rsplit('.', 1)[0] + '.fastg'
    output_post = f'/tmp/{uuid}.post'
    input_minimap2_kept = args.output + '.pre_assembly.clustered'
    input_mhap_kept = args.output + '.pre_assembly.mhap'
    output_ass_kept = args.output + '.assembly.unclustered'

    if args.precluster and args.precluster < 1: ### remove precluster?
        call([args.cd_hit_est_path, '-i', args.fasta, '-o', input_pre, '-c', str(args.precluster), '-T', str(args.threads), '-M', '0'])
        if args.keep_intermediate:
            call(['cp', input_pre, input_pre_kept])
    else:
        input_pre = args.fasta

    if args.identity_threshold and args.identity_threshold < 1:
        call([path + '/' + 'run_minimap2.py', '-f', input_pre, '-o', input_minimap2, '-i', str(args.identity_threshold), '-m', str(args.mismatch_size_threshold),
              '-g', str(args.indel_size_threshold), '-r', str(args.correction_repeats), '-t', str(args.threads), '--minimap2-path', args.minimap2_path, '--silent'])
        if args.keep_intermediate:
            call(['cp', input_minimap2, input_minimap2_kept])
    else:
        input_minimap2 = input_pre           

    Assembler(input_minimap2, args.ksize, args.minlen, args.mincov, output_ass)

    if args.postcluster and args.postcluster < 1:
        call([args.cd_hit_est_path, '-i', output_ass, '-o', output_post, '-c', str(args.postcluster), '-T', str(args.threads), '-M', '0'])
        call(['mv', output_post, args.output])
        if args.keep_intermediate:
            call(['mv', output_ass, output_ass_kept])
    else:
        call(['mv', output_ass, args.output])
    call(['mv', output_ass_fastg, output_final_fastg])
    
    




def parse_args():
    parser = ArgumentParser(description='Create a consensus pangenome assembly from a set of mOTU sequences')
    parser.add_argument('-f', '--fasta', type = str, required = True,
                        help = 'Input fasta file with the sequences from the mOTU')
    parser.add_argument('-p', '--precluster', type = float, default = 0,
                        help = 'Pre-assembly clustering identity for CD-HIT')
    parser.add_argument('-i', '--identity_threshold', type = float, default = 0.8,
                        help = 'Identity threshold (fraction) to initiate correction with minimap2')
    parser.add_argument('-m', '--mismatch-size-threshold', type = int, default = 10,
                        help = 'Maximum contiguous mismatch size that will be corrected')
    parser.add_argument('-g', '--indel-size-threshold', type = int, default = 10,
                        help = 'Maximum contiguous indel size that will be corrected')
    parser.add_argument('-r', '--correction-repeats', type = int, default = 1,
                        help = 'Maximum iterations for sequence correction')
    parser.add_argument('-k', '--ksize', type = int, default = 101,
                        help = 'Kmer size')
    parser.add_argument('-l', '--minlen', type = int, default = 0,
                        help = 'Scaffold length cutoff')
    parser.add_argument('-c', '--mincov', type = float, default = 0,
                        help = 'Scaffold coverage cutoff')
    parser.add_argument('-d', '--postcluster', type = float, default = 0,
                        help = 'Post-assembly clustering identity for CD-HIT')
    parser.add_argument('-t', '--threads', type = int, default = 1,
                        help = 'Number of processors to use')
    parser.add_argument('-o', '--output', type = str, default = 'assembly.fasta',
                        help = 'Output name')
    parser.add_argument('--cd-hit-est-path', type = str, default = 'cd-hit-est',
                        help = 'Path to the cd-hit-est executable')
    parser.add_argument('--minimap2-path', type = str, default = 'minimap2',
                        help = 'Path to the minimap2 executable')
    parser.add_argument('--keep-intermediate', action='store_true',
                        help = 'Keep intermediate files')
    return parser.parse_args()


if __name__ == '__main__':
    main(parse_args())

    
    
