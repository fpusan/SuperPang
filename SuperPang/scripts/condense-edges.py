#!/usr/bin/env python3

import sys
from os.path import dirname, realpath
path = dirname(realpath(sys.argv[0]))
sys.path.remove(path)
sys.path.insert(0, realpath(path + '/../..'))

from SuperPang.lib.utils import read_fasta, write_fasta

from collections import defaultdict

import graph_tool as gt
from graph_tool.all import Graph

def main():

    infasta, outname, ksize, minlen = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(200)
##    outinfo = outfasta.rsplit('.',1)[0] + '.info.tsv'

    name2id  = {}
    name2covLen = {}
    name2longName = {}
    name2vertex = {}
    vertex2name = {}

    seqs = read_fasta(infasta)
    edges = set()
    nodes = set()

    GS = Graph()

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
        if name not in name2vertex:
            v = GS.add_vertex()
            name2vertex[name] = v
            vertex2name[v] = name
        else:
            v = name2vertex[name]
        for succ in succs:
            if succ not in name2vertex:
                sv = GS.add_vertex()
                name2vertex[succ] = sv
                vertex2name[sv] = succ
            else:
                sv = name2vertex[succ]
            edges.add( (v, sv) )
            nodes.update( (v, sv) )

    GS.add_edge_list(edges)
    GS.vp.vfilt = GS.new_vertex_property('bool')
    GS.ep.efilt = GS.new_edge_property('bool')

    predecessors, successors = G2dicts(GS, name2vertex, vertex2name)

    realPredecessors = predecessors
    realSuccessors = successors
    addedEdges = set()
    addedNodes = set()
    
    newSeqs = {}
    newCore = {}

##    sources = {n for n in name2vertex if not predecessors[n]}
##    sinks = {n for n in name2vertex if not successors[n]}
##    for source in sources:
##        best = (None, 0)
##        for sink in sinks:
##            if source == sink:
##                continue
##            p = gt.topology.shortest_path(GS, name2vertex[source], name2vertex[sink])[0]
##            l = sum([len(seqs[name2longName[vertex2name[int(v)]]]) for v in p])
##            if l > best[1]:
##                best = (sink, l)
##        print(source, best[0], best[1])

    with open(outname + '.info.tsv', 'w') as outfile:
        outfile.write('id\tregions\tpath\ttrim\n')
        while True:
            origins = sorted({n for n in name2id if not predecessors[n] and n not in addedNodes}, key = lambda n: name2covLen[n], reverse = True)
            if not origins:
                if len(addedNodes) == len(name2id):
                    break
                else: # We have cycles
                    # Remove short loops
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
                    if not origins:
                        # Properly search for cycles
                        cycles = sorted(gt.topology.all_circuits(GS), key = lambda c: len(c))
                        origins = [sorted({vertex2name[c[0]] for c in cycles}, key = lambda n: name2covLen[n], reverse = True)[0]] # break cycles one at a time
                    
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
                                assert tag in ('core', 'aux', 'singleton', 'noinfo')
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
                            coreAdded = 0
                            for pos in sorted(tagDict):
                                tag = tagDict[pos]
                                if tag != last:
                                    if last:
                                        tagsStr.append(f'({lastStart},{pos-1},{last})')
                                        if last == 'core':
                                            coreSeq = seq[lastStart-1:pos-1]
                                            newCore[f'SUPERPANG_{len(newSeqs)}-core-{coreAdded}_length={len(coreSeq)}'] = coreSeq
                                            coreAdded += 1
                                    lastStart = pos
                                    last = tag
                            if last == 'core':
                                coreSeq = seq[lastStart-1:pos-1]
                                newCore[f'SUPERPANG_{len(newSeqs)}-core-{coreAdded}_length={len(coreSeq)}'] = coreSeq
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
                            if len(seq) >= minlen:
                                outfile.write(f'{cid}\t{tagsStr}\t{pathsStr}\t{trimStr}\n')
                                newSeqs[cid] = seq
                        addedNodes.update(path) # we can visit the same node twice in the same path, but not in different paths
                        break
                    edge = edges[0]
                    succ = edge[1]
                    if False:
                        addedNodes.add(succ) # uncomment to disallow visiting the same node twice in the same path
                    addedEdges.add(edge)
                    path.append(succ)

            

            newPreds = defaultdict(set)
            newSuccs = defaultdict(set)

            addedNodesV = {name2vertex[n] for n in addedNodes}
            addedEdgesV = {(name2vertex[s], name2vertex[t]) for s, t in addedEdges}
            set_edge_filter(GS, addedEdgesV, inverted = True)
            set_vertex_filter(GS, addedNodesV, inverted = True)

            predecessors, successors = G2dicts(GS, name2vertex, vertex2name)
            

    assert len(addedNodes) == len(name2id)

    write_fasta(newSeqs, outname + '.fasta')
##    write_fasta(newCore, outname + '.core.fasta') # this is somehow not working anymore (has very low completeness), not sure why but we were not using it anyways


def G2dicts(GS, name2vertex, vertex2name):
    predecessors = defaultdict(set)
    successors   = defaultdict(set)
    for n in name2vertex:
        v = name2vertex[n]
        for pv in GS.get_in_neighbors(v):
            if pv != v: # avoid self loops
                predecessors[n].add(vertex2name[pv])
        for sv in GS.get_out_neighbors(v):
            if sv != v:
                successors[n].add(vertex2name[sv])
    return predecessors, successors
        



def set_vertex_filter(GS, vSS, inverted = False):
    for vS in GS.vertices():
        GS.vp.vfilt[vS] = vS in vSS
    GS.set_vertex_filter(GS.vp.vfilt, inverted = inverted)

def set_edge_filter(GS, edges, inverted = False):
    for e in GS.edges():
        GS.ep.efilt[e] = (e.source(), e.target()) in edges
    GS.set_edge_filter(GS.ep.efilt, inverted = inverted)


if __name__ == '__main__':
    main()
