#!/usr/bin/env python3

from lib.utils import read_fasta, write_fasta
from collections import defaultdict
from sys import argv

infasta, outfasta, ksize = argv[1], argv[2], int(argv[3])
outinfo = outfasta.rsplit('.',1)[0] + '.info.tsv'

seqs = read_fasta(infasta)

predecessors = defaultdict(set)
successors = defaultdict(set)

name2id  = {}
name2covLen = {}
name2longName = {}


for longName, seq in seqs.items():
    name = longName.strip(';')
    name, *succs = name.split(':')
    if succs:
        succs = succs[0].split(',')
    fields = name.split('_')
    assert name not in name2id
    id_ = fields[1]
    name2id[name] = id_
    name2covLen[name] = (float(fields[5]), len(seq), id_) # add id_ to break ties deterministically
    name2longName[name] = longName
    for succ in succs:
        predecessors[succ].add(name)
        successors[name].add(succ)

realPredecessors = predecessors
realSuccessors = successors
addedEdges = set()
addedNodes = set()

newSeqs = {}

from itertools import combinations
for v1, v2 in combinations(name2covLen.values(), 2):
    if v1 == v2:
        print(v1, v2)

with open(outinfo, 'w') as outfile:
    outfile.write('id\tregions\tpath\ttrim\n')
    while True:
        origins = sorted({n for n in name2id if not predecessors[n] and n not in addedNodes}, key = lambda n: name2covLen[n], reverse = True)
        if not origins:
            if len(addedNodes) == len(name2id):
                break
            else: # We have cycles
                # Remove short loops (this must be improved to cover more cases, we should probably rewrite the whole thing using graph_tool)
                missing = set(name2id) - addedNodes
                badEdges = []
                for n in missing:
                    for pred in predecessors[n]:
                        if n in predecessors[pred] and len(successors[pred]) == 1:
                            badEdges.append( (n, pred) )
                for source, target in badEdges:
                    predecessors[target] = {pred for pred in predecessors[target] if pred != source}
                    successors[source]   = {succ for succ in successors[source]   if succ != target}
                origins = sorted({n for n in name2id if not predecessors[n] and n not in addedNodes}, key = lambda n: name2covLen[n], reverse = True)
                
        for n in origins:
            path = [n]
            while True:
                edges = [(path[-1], succ) for succ in sorted(successors[path[-1]], key = lambda n: name2covLen[n], reverse = True) if succ not in addedNodes]
                edges = [e for e in edges if e not in addedEdges] # this allows us to visit the same node twice in this path as long as we don't do it from the same predecessor (i.e. this is some kind of cycle)
                if not edges:
                    path = path[:-1] if path[-1] in addedNodes else path
                    assert path
                    # Build seq
                    trim_left  = any(pred in addedNodes for pred in realPredecessors[path[0]])
                    trim_right = any(succ in addedNodes for succ in realSuccessors[path[-1]])
                    if len(path) == 1:
                        seq = seqs[name2longName[path[0]]]
                        if trim_left:
                            seq = seq[ksize:]
                        if trim_right:
                            seq = seq[:-ksize]
                    else:
                        seq = ''.join(seqs[name2longName[n]][:-ksize] for n in path[:-1])
                        if trim_right:
                            seq += seqs[name2longName[path[-1]]][:-ksize]
                        else:
                            seq += seqs[name2longName[path[-1]]]
                        if trim_left:
                            seq = seq[ksize:]

                    # Go for it
                    if seq:
                        cid = f'SUPERPANG_{len(newSeqs)}_length={len(seq)}'
                        pathsStr = '[{}]'.format(','.join([name2id[n] for n in path]))
                        tagDict0 = defaultdict(set) # zero-indexed
                        pos = 0
                        for i, n in enumerate(path):
                            tag = name2id[n].split('-')[2]
                            s = seqs[name2longName[n]] # we don't trim for now
                            assert tag in ('core', 'aux', 'singleton')
                            if i > 0:
                                pos -= ksize
                            for _ in range(len(s)):
                                tagDict0[pos].add(tag)
                                pos += 1
                        if trim_left:
                            tagDict0 = {pos: tags for pos, tags in tagDict0.items() if pos >= ksize}
                        tagDict = {} # one-indexed
                        for i, pos in enumerate(sorted(tagDict0)):
                            tags = tagDict0[pos]
                            if 'core' in tags:
                                tag = 'core'
                            elif 'aux' in tags:
                                tag = 'aux'
                            else:
                                tag = 'singleton'
                            tagDict[i+1] = tag
                        if trim_right:
                            tagDict = {pos: tag for pos, tag in tagDict.items() if pos <= len(tagDict) - ksize}
                        assert len(tagDict) == len(seq)
                        last = None
                        lastStart = 0
                        tagsStr = []
                        for pos in sorted(tagDict):
                            tag = tagDict[pos]
                            if tag != last:
                                if last:
                                    tagsStr.append(f'({lastStart},{pos-1},{last})')
                                lastStart = pos
                                last = tag
                        tagsStr.append(f'({lastStart},{len(tagDict)},{last})')
                        tagsStr = '[{}]'.format(','.join(tagsStr))
                        if not trim_left and not trim_right:
                            trimStr = 'none'
                        elif trim_left and trim_right:
                            trimStr = 'both'
                        elif trim_left:
                            trimStr = 'left'
                        elif trim_right:
                            trimStr = 'right'
                        outfile.write(f'{cid}\t{tagsStr}\t{pathsStr}\t{trimStr}\n')
                        newSeqs[cid] = seq
                    addedNodes.update(path) # we can visit the same node twice in the same path, but not in different paths
                    break
                edge = edges[0]
                succ = edge[1]
    ##            addedNodes.add(succ) # uncomment to disallow visiting the same node twice in the same path
                addedEdges.add(edge)
                path.append(succ)

        newPreds = defaultdict(set)
        newSuccs = defaultdict(set)
        for n in predecessors:
            if n in addedNodes:
                continue
            newPreds[n] = {pred for pred in predecessors[n] if pred not in addedNodes}
            newSuccs[n] = {succ for succ in successors[n]   if succ not in addedNodes}
        predecessors = newPreds
        successors   = newSuccs

assert len(addedNodes) == len(name2id)

write_fasta(newSeqs, outfasta)

    
