#!/usr/bin/env python3

import sys
from os.path import dirname, realpath
path = dirname(realpath(sys.argv[0]))
sys.path.remove(path)
sys.path.insert(0, realpath(path + '/../..'))

from SuperPang.lib.utils import read_fastq, write_fastq, fasta2fastq, fastq2fasta, reverse_complement, print_time

from argparse import ArgumentParser

from collections import defaultdict
from subprocess import call, DEVNULL
import re


def main(args):

    assert args.fasta != args.output
    prefix = args.output.rsplit('.', 1)[0].split('/')[-1]
    current1 = '/tmp/' + prefix + '.current1.fastq'
    current2 = '/tmp/' + prefix + '.current2.fastq'
    fasta2fastq(args.fasta, current1)
    corrected = {}
    correction_tries = defaultdict(int)
    prev_tries = sum(correction_tries.values())
    rounds = 0

    print_time('Correcting input sequences with minimap2')
    print_time('round\tcorrectedFull\tcorrectedPartial\tmismatches\tindels\terrors')
    while True:
        rounds += 1
        mLen, iLen, oLen, corrected_thisround, nErrors = iterate(args, rounds, prefix, current1, current2, corrected, correction_tries,
                                                                 identity_threshold = args.identity_threshold,
                                                                 mismatch_size_threshold = args.mismatch_size_threshold,
                                                                 indel_size_threshold = args.indel_size_threshold)
        for b in corrected_thisround:
            if b not in corrected:
                corrected[b] = set()
            corrected[b].update(corrected_thisround[b])
        call(['mv', current2, current1])

        current_tries = sum(correction_tries.values())
        corrected_partial = current_tries - prev_tries

        print_time('\t'.join(map(str, [rounds, len(corrected_thisround), corrected_partial, iLen-mLen, oLen-iLen, nErrors])))

        if not corrected_thisround and not corrected_partial and rounds >= args.correction_repeats_min:
            # In some cases we report no corrected_thisround because fullCorrection = False for all sequences
            # We actually had corrected parts of them, and we will properly correct them in later rounds
            # correction_repeats_min ensures that we will try enough in order to get to fullCorrection = True and progress from there
            # in some cases we also reach correction_repeats_min without doing any full correction. we check with corrected_partial
            break
        if rounds == args.correction_repeats:
            break
        prev_tries = current_tries

    fastq2fasta(current1, args.output)


def iterate(args, rounds, prefix, infastq, outfastq, corrected_previously, correction_tries, identity_threshold, mismatch_size_threshold, indel_size_threshold):

    seqOrderOriginal = []
    seqs = read_fastq(infastq)
    seqOrderOriginal = list(seqs.keys()) # trust that dict remembers insertion order
    sortedSeqs = list(sorted(seqs, key = lambda a: len(seqs[a]), reverse = True))
    lexi = list(sorted((f'{i}' for i in range(len(sortedSeqs))), reverse = True))
    ori2minimap = {name: lexi[i] for i, name in enumerate(sortedSeqs)} # with "--dual=no" minimap2 will skip query-target pairs wherein the query name
    minimap2ori = {mm: name for name, mm in ori2minimap.items()}       #  is lexicographically greater than the target name.

    mLen = 0
    iLen = 0
    oLen = 0

    correctedSeqs = {}
    corrected_thisround = {}
    nErrors = 0

    tfile   = '/tmp/' + prefix + '.t.fastq'
    mmifile = '/tmp/' + prefix + '.t.mmi'
    qfile   = '/tmp/' + prefix + '.q.fastq'
    paf     = '/tmp/' + prefix + '.temp.paf'

    querybuf = set()
    targetbuf= []
    fLen = sum(len(s) for s in seqs.values())

    for i, target0 in enumerate(sortedSeqs):

        if target0 in corrected_thisround or target0 in correctedSeqs:
            continue

        for query0 in sortedSeqs[i+1:]:
            if query0 in corrected_previously and target0 in corrected_previously[query0]:
                continue
            if query0 in correctedSeqs or query0 in corrected_thisround:
                continue
            querybuf.add(query0)

        if not querybuf:
            continue
        
        
        targetbuf.append(target0)

        tLen = sum(len(seqs[target]) for target in targetbuf)
        qLen = sum(len(seqs[query]) for query in querybuf)

        if tLen * qLen < 200 * fLen * len(seqs) and i+1 < len(sortedSeqs):
            continue

        queries = [n for n in sortedSeqs if n in querybuf]
        targets = [n for n in targetbuf]
        querybuf  = set()
        targetbuf = []

        print_time(f'{rounds}\t{i+1}/{len(sortedSeqs)}', end = '\r')

        write_fastq({ori2minimap[target0]: seqs[target0] for target0 in targets}, tfile)
        write_fastq({ori2minimap[query0]:  seqs[query0]  for query0  in queries}, qfile)
        ecode = call([args.minimap2_path, '-Hk19', '-w5', '-e0', '-m100', '--rmq=yes', '--dual=no', '-DP', '--no-long-join', '-U50,500', '-g10k', '-s200',
                      tfile, qfile, '-t', str(args.threads), '-c', '--eqx', '-L'],
                      stdout=open(paf, 'w'), stderr = DEVNULL if args.silent else None)
        if ecode:
            print('\nError running minimap2\n')
            sys.exit(1)

        overlaps = {target0: defaultdict(list) for target0 in targets}
        bestCorr = {}

        for line in open(paf):
            # positions are zero indexed, start pos is closed, end pos is open
            query, queryLen, queryStart, queryEnd, isRC, target, targetLen, targetStart, targetEnd, matches, alignLen, qual, *_, cigar = line.strip().split('\t')
            query, target = minimap2ori[query], minimap2ori[target]
            if query == target:
                continue
            if query in corrected_thisround or query in correctedSeqs:
                continue
            if query in corrected_previously and target in corrected_previously[query]:
                continue

            queryLen, queryStart, queryEnd, targetLen, targetStart, targetEnd, matches, alignLen, qual = map(int, [queryLen, queryStart, queryEnd, targetLen, targetStart, targetEnd, matches, alignLen, qual])
                
            isRC = isRC == '-'
            cigar = cigar.replace('cg:Z:','')
            cigar = re.split('(\d+)',cigar)[1:]
            cigar = [(int(cigar[i]), cigar[i+1]) for i in range(0, len(cigar), 2)]
            idlen = sum([L for L, op in cigar if op in {'=', 'X'}])

            if len(cigar) == 1 and cigar[0][1] == '=':
                continue # perfect match, ignore

            if sum([L for L, op in cigar if op == '=']) / idlen < identity_threshold:
                continue # too many mismatches, ignore (indels are not counted here)

##            if sum([L for L, op in cigar if op == '=']) / alignLen < identity_threshold:
##                continue # too mamy mismatches and indels

            assert targetLen >= queryLen

            overlaps[target][query].append( ( queryStart, queryEnd, targetStart, targetEnd, cigar, isRC, matches, idlen, alignLen) )
        
        # query is the sequence we want to correct, target is the template
        # before we were using target0 and query0, but the code should work even if we kept using query and target as variable names
        for target in targets:
            if target in corrected_thisround or target in correctedSeqs:
                continue
            for query, ols in overlaps[target].items():
                # Check whether we can make this correction
                if query in correctedSeqs or query in corrected_thisround:
                    continue
                if query in corrected_previously:
                    assert target not in corrected_previously[query]

                #print(ori2minimap[target], ori2minimap[query], len(ols)) ####################
                # Ensure that the different regions to correct do not overlap
                fullCorrection = True
                ols = sorted(ols, key = lambda x: x[0]) # sort by query start
                lastEnd = -1
                goodOls = []
                for queryStart, queryEnd, targetStart, targetEnd, cigar, isRC, matches, idlen, alignLen in ols:
                    if queryStart > lastEnd:
                        goodOls.append( (queryStart, queryEnd, targetStart, targetEnd, cigar, isRC, matches, idlen, alignLen) )
                        lastEnd = queryEnd
                    else:
                        fullCorrection = False

                # Double-check
                lastEnd = -1
                for queryStart, queryEnd, targetStart, targetEnd, cigar, isRC, matches, idlen, alignLen in goodOls:
                    assert queryStart > lastEnd
                    lastEnd = queryEnd
        
                    mLen += matches
                    iLen += idlen
                    oLen += alignLen

                # Correct
                try:
                    correctedSeq = []
                    lastEnd = 0
                    for  queryStart, queryEnd, targetStart, targetEnd, cigar, isRC, matches, idlen, alignLen in goodOls:
                        correctedSeq.append(seqs[query][lastEnd:queryStart])
                        correctedSeq.append( correct_query(seqs[query], queryStart, queryEnd, seqs[target], targetStart, targetEnd,
                                                           cigar = cigar, isRC = isRC, mismatch_size_threshold = mismatch_size_threshold, indel_size_threshold = indel_size_threshold)
                                           )
                        lastEnd = queryEnd
                    correctedSeq.append( seqs[query][lastEnd:-1] )

                    correctedSeqs[query] = ''.join(correctedSeq)
                    if fullCorrection or correction_tries[(query, target)] == 5:
                        # we'll try to achieve full correction in 5 rounds, otherwise we'll consider it as corrected so we can use it as a template
                        corrected_thisround[query] = {target}
                    else:
                        correction_tries[(query, target)] += 1

                except AssertionError: # we failed some assertion in correct_query. I checked a case and it was actually because the alignment/CIGAR string generated by minimap2 was wrong
                    nErrors += 1        # so we just ignore this correction completely


    for header in seqs:
        if header not in correctedSeqs:
            correctedSeqs[header] = seqs[header]

    with open(outfastq, 'w') as outfile:
        for header in seqOrderOriginal:
            seq = correctedSeqs[header]
            quals = 'I'*len(seq)
            outfile.write(f'@{header}\n{seq}\n+\n{quals}\n')

    print_time(' '*40, end = '\r')
    
    return mLen, iLen, oLen, corrected_thisround, nErrors



def switch_cigar(cigar): # unused
    newCigar = []
    for L, op in cigar:
        if op == 'D':
            op = 'I'
        elif op == 'I':
            op = 'D'
        newCigar.append( (L, op) )
    return newCigar


def correct_query(query, queryStart, queryEnd, target, targetStart, targetEnd, cigar, isRC = False, mismatch_size_threshold = 10, indel_size_threshold = 100):
    ops = {'=', 'X', 'I', 'D'}
    query = query if not isRC else reverse_complement(query)
    queryPos = queryStart if not isRC else len(query) - queryEnd
    targetPos = targetStart
    correction = []
        
    for i, (L, op) in enumerate(cigar):
        assert op in ops
        if op == '=' or op == 'X':
            queryFrag = query[queryPos:queryPos+L]
            newFrag = target[targetPos:targetPos+L]
            if op == '=':
                assert newFrag == queryFrag
##                try:
##                    assert newFrag == queryFrag
##                except:
##                    print('ERROR!')
##                    print(i, L, op, isRC)
##                    print('>NEWFRAG')
##                    print(newFrag)
##                    print('>QUERY')
##                    print(query[queryPos:queryPos+L])
##                    print(cigar)
##                    print('qstart', queryStart)
##                    print(query)
##                    print('tstart', targetStart)
##                    print(target)
##                    raise
            if op == '=' or L <= mismatch_size_threshold:
                correction.append(newFrag) # add matches, correct mismatches
            else:
                correction.append(queryFrag)
            queryPos += L
            targetPos += L
        elif op == 'I':
            if L > indel_size_threshold:
                correction.append(query[queryPos:queryPos+L])
            queryPos += L
                
        elif op == 'D':
            if L <= indel_size_threshold:
                correction.append(target[targetPos:targetPos+L])
            targetPos += L

    correction = ''.join(correction)

##    assert correction == target[targetStart:targetEnd] # should be true if we are removing all indels no matter how big they are and queryStart and targetStart equal 0

    if isRC:
        correction = reverse_complement(correction)

    return correction



def parse_args():
    parser = ArgumentParser(description='Correct a set of sequences using mhap')
    parser.add_argument('-f', '--fasta', type = str, required = True,
                        help = 'Input fasta')
    parser.add_argument('-i', '--identity_threshold', type = float, default = 0.7,
                        help = 'Identity threshold (fraction) to initiate correction with minimap2')
    parser.add_argument('-m', '--mismatch-size-threshold', type = int, default = 100,
                        help = 'Maximum contiguous mismatch size that will be corrected')
    parser.add_argument('-g', '--indel-size-threshold', type = int, default = 100,
                        help = 'Maximum contiguous indel size that will be corrected')
    parser.add_argument('-r', '--correction-repeats', type = int, default = 1,
                        help = 'Maximum iterations for sequence correction')
    parser.add_argument('-n', '--correction_repeats_min', type = int, default = 1,
                        help = 'Minimum iterations for sequence correction')
    parser.add_argument('-o', '--output', type = str, default = 'corrected.fasta',
                        help = 'Output name')
    parser.add_argument('-t', '--threads', type = int, default = 1,
                        help = 'Number of processors to use')
    parser.add_argument('--minimap2-path', type = str, default = 'minimap2',
                        help = 'Path to the minimap2 executable')
    parser.add_argument('--silent', action = 'store_true',
                        help = 'Ignore stderr from minimus2')

    return parser.parse_args()

                
if __name__ == '__main__':
    main(parse_args())
