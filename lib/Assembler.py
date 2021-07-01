from collections import defaultdict
from itertools import combinations

import numpy as np

from lib.utils import read_fasta, write_fasta, reverse_complement, print_time

import graph_tool as gt
from graph_tool.all import Graph

### RIGHT NOW WE REMOVE Ns FROM THE SEQUENCES IN read_fasta. What do they mean and what should we do with them?

### see if we can put GS.clear_filters() inside the set filters method, instead of having a lot of them throughout the code

class Assembler:
##    @profile
    def __init__(self, fasta, ksize, minlen, mincov, out_fasta = '/home/fer/test.fasta'):

        self.ksize = ksize
        confRadius = ksize + 1000
        out_fastg = out_fasta.rsplit('.',1)[0] + '.fastg'        

        seqDict = read_fasta(fasta)
       
        print_time(f'Creating DBG from {len(seqDict)} sequences')
        
        self.edgeAbund = defaultdict(int)
        self.includedSeqs = 0
        self.kmer2vertex = {}
        self.seqPaths = {} # name: pathSet
        self.addedEdges = set()
        totalSize = sum(len(seq) for seq in seqDict.values())
        elapsedSize = 0
        lastSize = 0
        rcKmers = {}

        assembly_nodes = {}
        assembly_edges = {}

        ### Get kmers from sequences and create DBG
        G = Graph()
        G.set_fast_edge_removal(True)
        for name, seq in seqDict.items():
            elapsedSize += len(seq)
            if len(seq) <= self.ksize:
                continue
            else:
                self.includedSeqs += 1
                rcSeq = reverse_complement(seq)
                kmers = self.seq2kmers(seq)
                revkmers = self.seq2kmers(rcSeq)
                for k, rk in zip(kmers, revkmers[::-1]):
                    if k not in rcKmers:
                        rcKmers[k] = rk
                        rcKmers[rk] = k
                for i, k in enumerate(kmers):
                    if k not in self.kmer2vertex:
                        self.kmer2vertex[k] = int(G.add_vertex())
                        self.kmer2vertex[rcKmers[k]] = self.kmer2vertex[k]

                self.seqPaths[seq] = set(np.array([self.kmer2vertex[kmer] for kmer in kmers], dtype=np.uint32))

                for i in range(len(kmers)-1):
                    k1, k2 = kmers[i],kmers[i+1]
                    v1, v2 = self.kmer2vertex[k1], self.kmer2vertex[k2]
                    edge = (v1, v2)
                    if edge not in self.addedEdges: # Igual preferimos añadir el edge varias veces para contar el coverage?
                        G.add_edge(v1, v2) # O igual preferimos dejarlo undirected?
                        G.add_edge(v2, v1)
                        self.addedEdges.add( (v1, v2) )
                        self.addedEdges.add( (v2, v1) )
                if elapsedSize - lastSize > totalSize/100:
                    print_time(f'\t{round(100*elapsedSize/totalSize, 2)}% bases added, {G.num_vertices()} vertices, {G.num_edges()} edges         ', end = '\r')
                    lastSize = elapsedSize
                    
        G.shrink_to_fit()
        del(seqDict)
        print_time(f'\t100% bases added, {G.num_vertices()} vertices, {G.num_edges()} edges         ')
        
        ### Each kmer maps to one vertex, each vertex maps to two kmers
        self.vertex2kmers = defaultdict(list)
        for k, v in self.kmer2vertex.items():
            self.vertex2kmers[v].append(k)
        for v, kmers in self.vertex2kmers.items():
            assert len(kmers) == 2
        
        ###Link each vertex to the original sequence in which it appears
        vertex2names = defaultdict(set)
        for name, path in self.seqPaths.items():
            for v in path:
                vertex2names[v].add(name)

        ### Identify connected components
        print_time('Identifying connected components')
        vertex2comp = {v: c for v, c in enumerate(gt.topology.label_components(G, directed = False)[0])}
        comp2vertices = defaultdict(set)
        currentScaffold = 0
        for v, c in vertex2comp.items():
            comp2vertices[c].add(v)

        ### Process connected components
        for c, vs in comp2vertices.items():
            print_time(f'Working on comp {c+1}/{len(comp2vertices)}')
            ### Reconstruct non-branching paths
            print_time('\tCollecting non-branching paths')
            paths = []
            path_limits = defaultdict(list)
            inits = set()
            for v in vs:
                if len(set(G.get_all_neighbors(v)) - {v}) != 2: # -{v} to avoid self loops being init points
                    inits.add(v)

            # Break cycle if required
            isCycle = False
            ignoreEdge = set()
            if not inits:
                #is this a cycle?
                if min([len(set(G.get_all_neighbors(v))) for v in vs]) == 2: # no sources or sinks
                    for v in vs:
                        neighbors = G.get_all_neighbors(v)
                        if len(set(neighbors)) == 2: # find the first extender edge, and break there
                            inits.add(v)
                            inits.add(neighbors[0])
                            isCycle = True
                            ignoreEdge = {v, neighbors[0]}
                            break
                   
            for j, ini in enumerate(inits):
                for n in set(G.get_all_neighbors(ini)):
                    p = [ini, n]
                    last = ini
                    while True:
                        s = set(G.get_all_neighbors(p[-1])) - {last, p[-1]} #p[-1] to avoid stopping at self loops
                        if len(s) != 1:
                            break
                        last = p[-1]
                        s = s.pop()
                        if isCycle: # Break cycle at the pre-defined point
                            if last in ignoreEdge and s in ignoreEdge:
                                break
                        p.append(s)
                    p = tuple(p)
                    paths.append(p)

            # Corner case
            if len(paths) == 2: #### case in which only a starting/ending vertex is common to two sequences,
                                #### so both are joined by that vertex but that vertex seems an extender (1 in 1 out)
                                #### If there was really one path we are cutting it, but we'll join it later
                p = paths[0] # cut it in half
                if len(p)%2:
                    cut = int((len(p)+1)/2)
                    inits.add(p[cut-1])
                    paths = [p[:cut], p[:cut][::-1], p[-cut:], p[-cut:][::-1]]

            # Check that there are no len==1 paths
            for p in paths:
                assert len(p) > 1

            ### Build a graph of connected non-branching paths, and split into fwd and rev component
            print_time(f'\tBuilding sequence graph out of {len(paths)} paths')

            # Translate paths into sequences
            path2seq = {p: self.reconstruct_sequence(p) for p in paths}

            # Collect prefixes-suffixes for locating potential overlaps
            start2paths = defaultdict(set) # first k nucleotides
            sStart = {}
            sEnd = {}
            for path, seq in path2seq.items():
                start, end = seq[:ksize], seq[-ksize:]
                start2paths[start].add(path)
                sStart[path] = start
                sEnd[path] = end

            # Build sequence graph
            GS = Graph()
            path2vertexGS = {}
            for p in paths:
                assert p not in path2vertexGS
                path2vertexGS[p] = int(GS.add_vertex())
            for p1 in paths:
                confSeqs = vertex2names[p1[-1]]
                for p2 in start2paths[sEnd[p1]]:
                    cp = set(p1+p2[1:])
                    for name in confSeqs:
                        if cp.issubset(self.seqPaths[name]) or self.seqPaths[name].issubset(cp):
                            GS.add_edge(path2vertexGS[p1], path2vertexGS[p2])
                            break
                           
            vertex2path = {v: p for p,v in path2vertexGS.items()}
            
            ### Extract forward component
            print_time('\tExtracting forward component')
            comp2vertexGS = defaultdict(set) ## ahora uso v para todos los vértices. debería usar v para vértices de bruijn, y otra terminología para vertices secuencia. YA HE CAMBIADO PARTE
            comps =  gt.topology.label_components(GS, directed = False)[0]
            GS.vp.comps = comps
            GS.vp.vfilt = GS.new_vertex_property('bool')
            GS.ep.efilt = GS.new_edge_property('bool')

            for vS in GS.vertices():
                comp2vertexGS[GS.vp.comps[vS]].add(int(vS))

            # Try to find separated components in the seq graph that contain the same vertices in the kmer graph (thus are RC of each other)
            rcComps = {}
            for cS1, cS2 in combinations(comp2vertexGS,2):
                v1 = {v for vS in comp2vertexGS[cS1] for v in vertex2path[vS]}
                v2 = {v for vS in comp2vertexGS[cS2] for v in vertex2path[vS]}
                if v1 == v2:
                    rcComps[cS1] = cS2
                    rcComps[cS2] = cS1

            # If we didn't find a RC for all components this means that in some case fwd and rev are joined in the same component
            noRC = {cS for cS in comp2vertexGS if cS not in rcComps} 
            for cS in noRC: # for each of those evil components...
                psets = defaultdict(list)
                for vS in comp2vertexGS[cS]:
                    pset = frozenset(vertex2path[vS])
                    psets[pset].append(vS)
                assert len(psets) == len(comp2vertexGS[cS]) / 2 # each sequence path should have a reverse complement equivalent (same vertices in the kmer graph, reverse order)
                vS2rc = {}
                for vS1, vS2 in psets.values():
                    vS2rc[vS1] = vS2
                    vS2rc[vS2] = vS1
                self.set_vertex_filter(GS, comp2vertexGS[cS])
                assert {int(vS) for vS in GS.vertices()} == comp2vertexGS[cS] # did we alter vS indices when filtering?
                assert(len(set(gt.topology.label_components(GS, directed = False)[0]))) == 1
            
                # Separate fwd and rev
                vSs1 = set()
                vSs2 = set()
                added = set()

                for e in gt.search.bfs_iterator(GS): # is BFS the best choice? I need all fwd vertices to be visited before all rev vertices,
                    for vS in (e.source(), e.target()): # but it's hard to be sure without a DAG/ topoSort
                        if vS not in added:             # breaking cycles to get a DAG was too expensive
                            rc = vS2rc[vS]
                            added.add(vS)
                            added.add(rc)
                            vSs1.add(vS)
                            vSs2.add(rc)                   

                # Check that fwd and rev have the same vertices
                v1 = {v for vS in vSs1 for v in vertex2path[vS]}
                v2 = {v for vS in vSs2 for v in vertex2path[vS]}
                assert v1 == v2
                GS.clear_filters()

                # Identify subcomponents in vS1
                self.set_vertex_filter(GS, vSs1)
                subcomps1 = defaultdict(set)
                for n,c in zip(GS.vertices(), gt.topology.label_components(GS, directed = False)[0]):
                    subcomps1[c].add(n)
                GS.clear_filters()
                # Identify subcomponents in vS2
                self.set_vertex_filter(GS, vSs2)
                subcomps2 = defaultdict(set)
                for n,c in zip(GS.vertices(), gt.topology.label_components(GS, directed = False)[0]):
                    subcomps2[c].add(n)
                GS.clear_filters()
                # Identify RC components
                first = True
                assert len(subcomps1) == len(subcomps2)
                for scS1, svS1 in subcomps1.items():
                    v1 = {v for vS in svS1 for v in vertex2path[vS]}
                    for scS2, svS2 in subcomps2.items():
                        v2 = {v for vS in svS2 for v in vertex2path[vS]}
                        if v1 == v2:
                            i1 = cS if first else len(comp2vertexGS) # recycle the previous component index for cS1
                            comp2vertexGS[i1] = svS1
                            first = False
                            i2 = len(comp2vertexGS)
                            comp2vertexGS[i2] = svS2
                            rcComps[i1] = i2
                            rcComps[i2] = i1
                            break
                    else: # assert that all our subcomponents have RC components
                        assert False
                
                GS.clear_filters()

            assert len(rcComps) == len(comp2vertexGS) # all our components have a rc component, splitting was successful


            # Treat the component with more ATG as fwd and report it
            ATGs = defaultdict(int)
            for cS, vSs in comp2vertexGS.items():
                for vS in vSs:
                    ATGs[cS] += path2seq[vertex2path[vS]].count('ATG')
            added = set()

            compS2paths = defaultdict(set)
            
            for cS1, cS2 in rcComps.items():
                if cS1 in added:
                    continue
                fwdComp = cS1 if ATGs[cS1] >= ATGs[cS2] else cS2
                fwdPaths = [vertex2path[vS] for vS in comp2vertexGS[fwdComp]]
                for p in fwdPaths:
                    compS2paths[currentScaffold].add(p)
                
                added.add(cS1)
                added.add(cS2)
                currentScaffold += 1

##            # Check that the  non-init vertices appear only once
##            vcov = defaultdict(int)
##            for p in paths:
##                for v in p:
##                    vcov[v] += 1
##            for v, cov in vcov.items(): # THIS ASSERTION FAILS IN SOME CASES IN WHICH THE FWD AND REV VERSION OF THE KMER STILL ARE IN THE SAME SEQUENCE
##                if v not in inits:      # CD-HIT SHOULD TAKE CARE OF THE PROBLEM AFTER THIS
##                    assert cov == 1
            
            # Check that we included all the input vertices
            assert vs == {v for p in paths for v in p}

            ### Process each scaffold in the sequence graph
            for scaffold, paths in compS2paths.items():
                # Test that the scaffolds are really connected
                GS.clear_filters() # NEVER FORGET!!!
                self.set_vertex_filter(GS, {path2vertexGS[p] for p in paths})
                msg = f'\tScaffold {scaffold}, {GS.num_vertices()} vertices, {GS.num_edges()} edges'
                print_time(msg, end = '\r')
                assert len(set(gt.topology.label_components(GS, directed = False)[0])) == 1

                # Compute scaffold length and skip if too short
                scaffoldLen = 0
                for p in paths:
                    seq = path2seq[p]
                    if GS.get_out_neighbors(path2vertexGS[p]).shape[0]:
                        seq = seq[:-ksize] # trim the overlap unless this is a terminal node
                    scaffoldLen += len(seq)
                msg += f', {scaffoldLen} bases'
                if minlen and scaffoldLen < minlen:
                    msg += ' (< minlen, IGNORED)'
                    print_time(msg)
                    continue

                # Compute scaffold coverage
                scaffoldCov = np.mean([len(vertex2names[v]) for p in paths for v in p])
                msg += f', cov {round(scaffoldCov, 2)}'
                if mincov and scaffoldCov < mincov:
                    msg += ' (< mincov, IGNORED)'
                    print_time(msg)
                    continue
                
                print_time(msg, end='\r')
                
                # Identify and flatten bubbles
                bubble_paths = set()
                path_limits = defaultdict(list)
                for p in paths:
                    path_limits[(p[0], p[-1])].append(p)    
                m1 = 0
                for (start, end), ps in path_limits.items():
                    minpathlen = min(len(p) for p in ps)
                    if len(ps) > 1: # more than one path with similar start and end
                        if minpathlen <= 2*ksize-1: # fast way if mismatches are separated enough
                            bubble_paths.update(ps[1:])
                            m1 += len(ps)-1
                paths = [p for p in paths if p not in bubble_paths]
                msg += f', removed {m1} bubbles'
                print_time(msg, end = '\r')
                GS.clear_filters()
                self.set_vertex_filter(GS, {path2vertexGS[p] for p in paths})

                

                # Join contiguous unbranching paths (can happen after solving the fwd component and popping bubbles
                predecessors = defaultdict(list)
                successors = defaultdict(list)
                for p1 in paths:
                    for vS2 in GS.get_out_neighbors(path2vertexGS[p1]):
                        if vertex2path[vS2] not in bubble_paths:
                            successors[p1].append(vertex2path[vS2])
                    for vS2 in GS.get_in_neighbors(path2vertexGS[p1]):
                        if vertex2path[vS2] not in bubble_paths:
                            predecessors[p1].append(vertex2path[vS2])

                
                extenders = {p for p in paths if len(successors[p]) <= 1 and (successors[p] or predecessors[p])}
                extenders = extenders | {p for p in paths if len(predecessors[p]) == 1} ### added later, consider removing this line if we start failing assertions
                joinStarters =  {p for p in extenders if len(predecessors[p]) != 1 or len(successors[predecessors[p][0]]) > 1}

                if not joinStarters: # Is this a cycle
                    if min(len(ps) for ps in predecessors.values()) == 1: # This does not cover the case in which there is a cycle but also a linear component
                        joinStarters = {list(extenders)[0]}
                
                joinExtenders = {p for p in extenders if p not in joinStarters}
                assert extenders == joinStarters | joinExtenders

                paths = [p for p in paths if p not in joinStarters | joinExtenders] # this also removes some joinStarters that can't be extended (no successors), but we'll add them back in the next bit
                
                visited = set()
                for start in joinStarters:
                    new_path = start
                    end = start
                    succs = successors[start]
                    succ = succs[0] if succs else None
                    while succ in joinExtenders:
                        visited.add(succ)
                        new_path += succ[1:]
                        end = succ
                        succs = successors[succ]
                        succ = succs[0] if succs else None
                    paths.append(new_path)
                    sStart[new_path] = sStart[start]
                    sEnd[new_path] = sEnd[end]
                    path2seq[new_path] = self.reconstruct_sequence(new_path)
                    predecessors[new_path] = predecessors[start]
                    for pred in predecessors[new_path]:
                        successors[pred].append(new_path)
                    successors[new_path] = successors[end]
                    for succ in successors[new_path]:
                        predecessors[succ].append(new_path)

                assert visited == joinExtenders

                # Update scaffold graph
                pathSet = set(paths)
                for p, preds in predecessors.items():
                    predecessors[p] = [pred for pred in preds if pred in pathSet]
                predecessors = defaultdict(list, [(p, preds) for p, preds in predecessors.items() if p in pathSet])
                for p, succs in successors.items():
                    successors[p]   = [succ for succ in succs if succ in pathSet]
                successors   = defaultdict(list, [(p, succs) for p, succs in successors.items()   if p in pathSet])

                    
                # Check that paths with one predecessor and successor happen only around branching points
                for p in paths:
                    if len(predecessors[p])==1 and len(successors[p])==1:
                        pred = predecessors[p][0]
                        succ = successors[p][0]
                        assert len(successors[pred]) > 1 or len(predecessors[pred]) > 1
                        assert len(successors[succ]) > 1 or len(predecessors[succ]) > 1

                msg += f', found {len(paths)} paths'

                # Get path coverages
                path2cov = {p: np.mean([len(vertex2names[v]) for v in p]) for p in paths}
                
                # Write assembly graph edges
                path2edgeName = {}
                for i, p in enumerate(paths):
                    path2edgeName[p] = f'EDGE_Sc{scaffold}-{i}_length_{len(path2seq[p])}_cov_{path2cov[p]}'
                
                for p in paths:
                    succs = ','.join(path2edgeName[s] for s in successors[p])
                    succs = f':{succs}' if succs else ''
                    assembly_edges[f'{path2edgeName[p]}{succs};'] = path2seq[p]
                # Write contigs
                for i, p  in enumerate(paths):
                    cov = 1
                    seq = path2seq[p]
                    if successors[p]: # trim the overlap unless this is a terminal node
                        seq = seq[:-ksize]
                    assembly_nodes[f'NODE_Sc{scaffold}-{i}_length_{len(seq)}_cov_{path2cov[p]}'] = seq

                print_time(msg)

            
                
        ### Output
        write_fasta(assembly_nodes, out_fasta)
        write_fasta(assembly_edges, out_fastg)


    @staticmethod
    def set_vertex_filter(GS, vSS):
        for vS in GS.vertices():
            GS.vp.vfilt[vS] = vS in vSS
        GS.set_vertex_filter(GS.vp.vfilt)


    def seq2kmers(self, seq):
        return [seq[i:i+self.ksize] for i in range(len(seq)-self.ksize+1)]


    def reconstruct_sequence(self, path):
        kmers = []
        # Each vertex corresponds to two kmers (fwd and rev)
        # Which kmers give us a meaningful sequence?
        # We are reading this in a given direction (from path[0] to path[-1])
        # For the starting kmer (in position 0), pick the one that overlaps with any of the two kmers in position 1
        if any(k[:-1]==self.vertex2kmers[path[0]][0][1:] for k in self.vertex2kmers[path[1]]):
            kmers.append(self.vertex2kmers[path[0]][0])
        if any(k[:-1]==self.vertex2kmers[path[0]][1][1:] for k in self.vertex2kmers[path[1]]):
            kmers.append(self.vertex2kmers[path[0]][1])
        assert len(kmers) == 1
        # Then just keep advancing in the path, for each vertex pick the kmer that overlaps with the previously-picked kmer
        for v in path[1:]:
            ks = self.vertex2kmers[v]
            toadd = []
            for k in ks:
                if k[:-1] == kmers[-1][1:]: #and k in self.kmerFollow[kmers[-1]]:
##                    assert k in self.kmerFollow[kmers[-1]] # the check for self.kmerFollow[kmers[-1]] was commented as it was always true if k[:-1] == kmers[1][1:]
                    toadd.append(k)
            assert len(toadd) <= 2
            if not toadd:
                break
            kmers.extend(toadd)
        return self.seq_from_kmers(kmers)
        
    @staticmethod
    def seq_from_kmers(kmers):
        fullseq = [kmers[0]]
        extra = [kmer[-1] for kmer in kmers[1:]]
        fullseq.extend(extra)
        return ''.join(fullseq)

        

