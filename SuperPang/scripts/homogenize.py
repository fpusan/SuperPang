#!/usr/bin/env python3

import sys
from os.path import dirname, realpath
path = dirname(realpath(sys.argv[0]))
sys.path.remove(path)
sys.path.append(realpath(path + '/../..'))

from SuperPang.lib.utils import read_fastq, fasta2fastq, fastq2fasta, fasta2fastq, reverse_complement, print_time

from argparse import ArgumentParser

import numpy as np
from subprocess import call, DEVNULL
import re


def main(args):

    assert args.fasta != args.output
    prefix = args.output.rsplit('.', 1)[0].split('/')[-1]
    current1 = '/tmp/' + prefix + '.current1.fastq'
    paf = '/tmp/' + prefix + '.current.paf'
    current2 = '/tmp/' + prefix + '.current2.fastq'
    fasta2fastq(args.fasta, current1)
    corrected = {}
    rounds = 0

    print_time('Correcting input sequences with minimap2')
    print_time('round\t#corrected\tmLen\toLen\tsim\terrors')
    while True:
        rounds += 1
        #call([args.minimap2_path, '-x', 'asm20', '-f', '100', current1, current1, '-t', str(args.threads), '-c', '--eqx', '-L'], stdout=open(paf, 'w'), stderr = DEVNULL if args.silent else None)
        call([args.minimap2_path, '-x', 'ava-pb', '-f', '100', current1, current1, '-t', str(args.threads), '-c', '--eqx', '-L'], stdout=open(paf, 'w'), stderr = DEVNULL if args.silent else None)
        mLen, oLen, corrected_thisround, nErrors = iterate(current1, paf, current2, corrected, identity_threshold = args.identity_threshold, mismatch_size_threshold = args.mismatch_size_threshold, indel_size_threshold = args.indel_size_threshold)
        print_time('\t'.join(map(str, [rounds, len(corrected_thisround), mLen, oLen, round(mLen/oLen, 4), nErrors])))
        for b in corrected_thisround:
            if b not in corrected:
                corrected[b] = set()
            corrected[b].update(corrected_thisround[b])
        call(['mv', current2, current1])
        if not corrected_thisround and rounds >= args.correction_repeats_min:
            # In some cases we report no corrected_thisround because fullCorrection = False for all sequences
            # We actually had corrected parts of them, and we will properly correct them in later rounds
            # correction_repeats_min ensures that we will try enough in order to get to fullCorrection = True and progress from there
            break
        if rounds == args.correction_repeats:
            break
    fastq2fasta(current1, args.output)
            
        
        
    
def iterate(infastq, paf, outfastq, corrected_previously, identity_threshold, mismatch_size_threshold, indel_size_threshold):

    corrected_thisround = {}
    
    seqOrderOriginal = []
    seqs = read_fastq(infastq)
    minimusNameConversion = {}
        
    seqOrderOriginal = list(seqs.keys()) # trust that dict remembers insertion order

    for name in seqOrderOriginal: # since minimus2 removes everything after the first space
        newName = name.split(' ')[0]
        assert newName not in minimusNameConversion # make sure that names are still different. If not, kill ourselves and blame society
        minimusNameConversion[newName] = name

            
    sortedSeqs = list(sorted(seqs, key = lambda a: len(seqs[a]), reverse = True))

    mLen = 0
    oLen = 0
    overlaps = {}
    for line in open(paf):
        # positions are zero indexed, start pos is closed, end pos is open
        query, queryLen, queryStart, queryEnd, isRC, target, targetLen, targetStart, targetEnd, matches, alignLen, qual, *_, cigar = line.strip().split('\t')
        if query == target:
            continue
        query, target = minimusNameConversion[query], minimusNameConversion[target]
        queryLen, queryStart, queryEnd, targetLen, targetStart, targetEnd, matches, alignLen, qual = map(int, [queryLen, queryStart, queryEnd, targetLen, targetStart, targetEnd, matches, alignLen, qual])
        isRC = isRC == '-'
        cigar = cigar.replace('cg:Z:','')
        cigar = re.split('(\d+)',cigar)[1:]
        cigar = [(int(cigar[i]), cigar[i+1]) for i in range(0, len(cigar), 2)]
        mLen += matches
        oLen += alignLen

        if len(cigar) == 1 and cigar[0][1] == '=':
            continue # perfect match, ignore

        if sum([L for L, op in cigar if op == '=']) / sum([L for L, op in cigar if op in {'=', 'X'}]) < identity_threshold:
            continue # too many mismatches, ignore (indels are not counted here)

        if targetLen >= queryLen:
            if target not in overlaps:
                overlaps[target] = {}
            if query not in overlaps[target]:
                overlaps[target][query] = list()
            overlaps[target][query].append( ( queryStart, queryEnd, targetStart, targetEnd, cigar, isRC) )
        elif queryLen > targetLen:
            if query not in overlaps:
                overlaps[query] = {}
            if target not in overlaps[query]:
                overlaps[query][target] = list()
            overlaps[query][target].append( ( targetStart, targetEnd, queryStart, queryEnd, switch_cigar(cigar) if not isRC else switch_cigar(cigar)[::-1], isRC) )
        else:
            pass # do something different if lenghts are equal????
    
    correctedSeqs = {}
    nErrors = 0
    # query is the sequence we want to correct, target is the template
    for target in sortedSeqs: # go from the longest target to the shortest
        if target not in overlaps: # if it doesn't overlap, go on
            continue
        if target in correctedSeqs: # if this target was corrected by a longer sequence, go on
            continue
        for query, ols in overlaps[target].items():
            # Check whether we can make this correction
            if query in correctedSeqs: # if this query was already corrected by a longer sequence in this round, go on
                continue
            if query in corrected_previously:
                if target in corrected_previously[query]: # if this query was already corrected by this target in a previous round, go on
                    continue
##                if a in corrected_previously:
##                    if corrected_previously[query] & corrected_previously[target]: # if both were corrected by the same target in a previous round, go on
##                        continue
            assert query not in correctedSeqs and query not in corrected_thisround

            # Ensure that the different regions to correct do not overlap
            fullCorrection = True
            ols = sorted(ols, key = lambda x: x[0]) # sort by query start
            lastEnd = -1
            goodOls = []
            for queryStart, queryEnd, targetStart, targetEnd, cigar, isRC in ols:
                if queryStart > lastEnd:
                    goodOls.append( (queryStart, queryEnd, targetStart, targetEnd, cigar, isRC) )
                    lastEnd = queryEnd
                else:
                    fullCorrection = False
                    
            # Double-check
            lastEnd = -1
            for queryStart, queryEnd, targetStart, targetEnd, cigar, isRC in goodOls:
                assert queryStart > lastEnd
                lastEnd = queryEnd
            # Correct
            try:
                correctedSeq = []
                lastEnd = 0
                for  queryStart, queryEnd, targetStart, targetEnd, cigar, isRC in goodOls:
                    correctedSeq.append(seqs[query][lastEnd:queryStart])
                    correctedSeq.append( correct_query(seqs[query], queryStart, queryEnd, seqs[target], targetStart, targetEnd,
                                                       cigar = cigar, isRC = isRC, mismatch_size_threshold = mismatch_size_threshold, indel_size_threshold =indel_size_threshold)
                                         )
                    lastEnd = queryEnd
                correctedSeq.append( seqs[query][lastEnd:-1] )
                
                correctedSeqs[query] = ''.join(correctedSeq)
                if fullCorrection:
                    corrected_thisround[query] = {target}
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

    return mLen, oLen, corrected_thisround, nErrors


def switch_cigar(cigar):
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
