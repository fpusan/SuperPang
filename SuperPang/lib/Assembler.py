from collections import defaultdict, namedtuple
from itertools import combinations

import numpy as np
from scipy.sparse import lil_matrix, csr_matrix

from SuperPang.lib.utils import read_fasta, reverse_complement, print_time

import graph_tool as gt
from graph_tool.all import Graph

from mappy import Aligner


### see if we can put GS.clear_filters() inside the set filters method, instead of having a lot of them throughout the code

class Assembler:
##    @profile
    def __init__(self, fasta, ksize):

        self.ksize = ksize

        seqDict = read_fasta(fasta, Ns = 'split')

        print_time(f'Creating DBG from {len(seqDict)} sequences')
        
        self.edgeAbund = defaultdict(int)
        self.includedSeqs = 0
        self.hash2vertex = {}
        self.vertex2source = {}
        self.seqPaths = {} # name: pathSet
        totalSize = sum(len(seq) for seq in seqDict.values())
        elapsedSize = 0
        lastSize = 0
        
        ### Get kmers from sequences and create DBG
        edges = set()
        for name, seq in seqDict.items():
            elapsedSize += len(seq)
            if len(seq) <= self.ksize:
                continue
            else:
                self.includedSeqs += 1
                rcSeq = reverse_complement(seq)
                hashes = [hash(k) for k in self.seq2kmers(seq)]
                revhashes = [hash(rk) for rk in self.seq2kmers(rcSeq)]
                rcDict = {}
                for i, (h, rh) in enumerate(zip(hashes, revhashes[::-1])):
                    rcDict[h] = rh
                    rcDict[rh] = h
                    if h not in self.hash2vertex and rh not in self.hash2vertex:
                        v = len(self.vertex2source)
                        self.hash2vertex[h] = v
                        self.vertex2source[v] = (name, i)
                    else:
                        if h in self.hash2vertex:
                            v = self.hash2vertex[h]
                        else:
                            v = self.hash2vertex[rh]
                sp = set()
                for h in hashes:
                    if h in self.hash2vertex:
                        sp.add(self.hash2vertex[h])
                    else:
                        sp.add(self.hash2vertex[rcDict[h]])
                self.seqPaths[name] = sp
                    

                for i in range(len(hashes)-1):
                    h1, h2 = hashes[i],hashes[i+1]
                    if h1 in self.hash2vertex:
                        v1 = self.hash2vertex[h1]
                    else:
                        v1 = self.hash2vertex[rcDict[h1]]
                    if h2 in self.hash2vertex:
                        v2 = self.hash2vertex[h2]
                    else:
                        v2 = self.hash2vertex[rcDict[h2]]
                    edges.add( (v1, v2) )
                    edges.add( (v2, v1) )
                if elapsedSize - lastSize > totalSize/100:
                    print_time(f'\t{round(100*elapsedSize/totalSize, 2)}% bases added, {len(self.vertex2source)} vertices, {len(edges)} edges         ', end = '\r')
                    lastSize = elapsedSize

##        import sys; print('\n'); sys.exit()

        self.G = Graph()
        self.G.add_edge_list(edges)
        self.seqDict = seqDict
        del(edges)
        print_time(f'\t100% bases added, {self.G.num_vertices()} vertices, {self.G.num_edges()} edges         ')
        

##    @profile
    def run(self, minlen, mincov, bubble_identity_threshold, genome_assignment_threshold):
        ### Start assembly
        Contig = namedtuple('Contig', ['scaffold', 'i', 'seq', 'tseq', 'cov', 'origins', 'successors'])
        contigs = {}
        addedPaths = set()
        
        ### Identify connected components
        print_time('Identifying connected components')
        vertex2comp = {v: c for v, c in enumerate(gt.topology.label_components(self.G, directed = False)[0])}
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
                if len(set(self.G.get_all_neighbors(v)) - {v}) != 2: # -{v} to avoid self loops being init points
                    inits.add(v)

            # Break cycle if required
            isCycle = False
            ignoreEdge = set()
            if not inits:
                #is this a cycle?
                if min([len(set(self.G.get_all_neighbors(v))) for v in vs]) == 2: # no sources or sinks
                    for v in vs:
                        neighbors = self.G.get_all_neighbors(v)
                        if len(set(neighbors)) == 2: # find the first extender edge, and break there
                            inits.add(v)
                            inits.add(neighbors[0])
                            isCycle = True
                            ignoreEdge = {v, neighbors[0]}
                            break
                   
            for j, ini in enumerate(inits):
                for n in set(self.G.get_all_neighbors(ini)):
                    p = [ini, n]
                    last = ini
                    while True:
                        s = set(self.G.get_all_neighbors(p[-1])) - {last, p[-1]} #p[-1] to avoid stopping at self loops
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

            # Remove paths with terminally duplicated vertices
            # Those broke our code once. we had three vSs with paths (0,0,1,2) (0,1,2) and (2,1,0) under the same pset (see below)
            paths = [p for p in paths if p[0] != p[1] and p[-1] != p[-2]]

            # Translate paths into sequences
            path2seq = {p: self.reconstruct_sequence(p) for p in paths}

            # Collect prefixes-suffixes for locating potential overlaps
            start2paths = defaultdict(set) # first k nucleotides
            sStart = {}
            sEnd = {}
            for path, seq in path2seq.items():
                start, end = seq[:self.ksize], seq[-self.ksize:]
                start2paths[start].add(path)
                sStart[path] = start
                sEnd[path] = end

            # Build sequence graph
            GS = Graph()
            path2vertexGS = {}
            good = 0
            lessGood = 0

            for p in paths:
                assert p not in path2vertexGS
                path2vertexGS[p] = int(GS.add_vertex())
            for p1 in paths:
                for p2 in start2paths[sEnd[p1]]:
                    GS.add_edge(path2vertexGS[p1], path2vertexGS[p2])
                           
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
                for vSs in psets.values():
                    assert len(vSs) == 2
                                
                assert len(psets) == len(comp2vertexGS[cS]) / 2 # each sequence path should have a reverse complement equivalent (same vertices in the kmer graph, reverse order)
                vS2rc = {}
                for vS1, vS2 in psets.values():
                    vS2rc[vS1] = vS2
                    vS2rc[vS2] = vS1
                self.set_vertex_filter(GS, comp2vertexGS[cS])
                assert {int(vS) for vS in GS.vertices()} == comp2vertexGS[cS] # did we alter vS indices when filtering?
                assert(len(set(gt.topology.label_components(GS, directed = False)[0]))) == 1 # check that this is only one component

                # Map the input sequences to nodes in the sequence graph (in the original orientation of the input sequence)
                seqPathsGS = defaultdict(set)
                addedNodes = set()
                for pset, (vS1, vS2) in psets.items():
                    assert vertex2path[vS1] == vertex2path[vS2][::-1]
                    breakOut = False
                    for name, path in self.seqPaths.items():
                        if pset.issubset(path):
                            for vS in (vS1, vS2):
                                if path2seq[vertex2path[vS]] in self.seqDict[name]: # check the string to make sure that we only add nodes from the same orientation
                                    addedNodes.update( (vS1, vS2) )
                                    seqPathsGS[name].add(vS)
                                    break

                # Add nodes in the sequence graph that weren't fully contained in an original sequence
                for vS in comp2vertexGS[cS] - addedNodes:
                    fakeName = hash( (len(seqPathsGS), vS) )
                    seqPathsGS[fakeName] = {vS} | set(GS.get_all_neighbors(vS)) # add its neighbors so that it'll overlap with the real sequences in at least one vS
                
                seqPathsGS_rev = {name: {vS2rc[vS] for vS in vSs} for name, vSs in seqPathsGS.items()}

                x = {vS for vSs in seqPathsGS.values() for vS in vSs}
                y = {vS2rc[vS] for vS in x}
                
                # Separate fwd and rev
                vSs1 = set()
                refSeqNames = sorted(seqPathsGS, key = lambda name: seqPathsGS[name]) # sort by decreasing number of nodes in the GS

                addedNames = {refSeqNames[0]}
                vSs1.update(seqPathsGS[refSeqNames[0]]) # start from the longest seq

                while True:
                    remaining = [n for n in refSeqNames if n not in addedNames] # keep them sorted
                    if not remaining:
                        break
                    best_name = ''
                    best_ol = -1 # -1 so that even if there are no overlaps (there is more than one subcomponent, and the current one is full) we add something
                    for name in remaining: # keep starting from the longest seq, so that if there are no overlaps we add the first longest one
                        fwd = seqPathsGS[name]
                        rev = seqPathsGS_rev[name]
                        fwd_ol = len(fwd & vSs1) / len(fwd) # find the best overlap (looking in both orientations) to the existing fwd component
                        rev_ol = len(rev & vSs1) / len(rev)
                        if fwd_ol > best_ol and fwd_ol >= rev_ol:
                            best_ol, best_name, best_vSs = fwd_ol, name, fwd
                        elif rev_ol > best_ol and rev_ol > fwd_ol:
                            best_ol, best_name, best_vSs = rev_ol, name, rev
                    assert best_name
                    addedNames.add(best_name)
                    best_vSs = {vS for vS in best_vSs if vS2rc[vS] not in vSs1} # don't add reverse complements at this stage
                                                                                # some are really there and needed to keep everything into a single component
                                                                                # some will be duplications (e.g. because one of the input genomes has an inverted region)
                    vSs1.update(best_vSs)

                vSs2 = {vS2rc[vS] for vS in vSs1}
                v1 = {v for vS in vSs1 for v in vertex2path[vS]}
                v2 = {v for vS in vSs2 for v in vertex2path[vS]}

                assert (vSs1 | vSs2) == comp2vertexGS[cS]

                # Identify subcomponents in vS1
                GS.clear_filters()
                self.set_vertex_filter(GS, vSs1)
                subcomps1 = defaultdict(set)
                for vS,c in zip(GS.vertices(), gt.topology.label_components(GS, directed = False)[0]):
                    subcomps1[c].add(vS)
                GS.clear_filters()
                # Identify subcomponents in vS2
                self.set_vertex_filter(GS, vSs2)
                subcomps2 = defaultdict(set)
                for vS,c in zip(GS.vertices(), gt.topology.label_components(GS, directed = False)[0]):
                    subcomps2[c].add(vS)
                GS.clear_filters()
                # Identify RC components
                assert len(subcomps1) == len(subcomps2)
                
                subcomp2vertex1 = {scS1: {v for vS in svS1 for v in vertex2path[vS]} for scS1, svS1 in subcomps1.items()} # so we don't have to compute this inside the loop
                subcomp2vertex2 = {scS2: {v for vS in svS2 for v in vertex2path[vS]} for scS2, svS2 in subcomps2.items()}


                # Try to merge into one fwd and one rev subcomp
                # In some cases the fwd and rev sequences can be actually be present in the same genome
                # So we identify the largest fwd subcomponent, throw in the rest of components in fwd and rev, and hope that they become connected when everything is together
                largestSubcomp = list(sorted(subcomps1, key = lambda scS: len(subcomps1[scS]), reverse = True))[0]
                combined = subcomps1[largestSubcomp]
                for scS in subcomps1:
                    if scS != largestSubcomp:
                        fwd = subcomps1[scS]
                        rev = {vS2rc[vS] for vS in fwd}
                        combined = combined | fwd | rev
                
                subcomps1_corrected = defaultdict(set)
                subcomps2_corrected = defaultdict(set)
                self.set_vertex_filter(GS, combined)
                for vS, c  in zip(GS.vertices(), gt.topology.label_components(GS, directed = False)[0]):
                    subcomps1_corrected[c].add(vS)
                    subcomps2_corrected[c].add(vS2rc[vS])


                # The lines above may make a rev subcomponent may end up being contained in a fwd component
                # E.g. let us have initially three components A, B, C (fwd) and A', B', C' (rev)
                # We correct and it turns out that B' got inside A so by the code above we end up with AB', B, C (fwd) and A'B, B', C' (rev)
                # So we want to remove B from fwd and B' from rev, so that we don't end up having duplicate paths
                bad1 = set()
                bad2 = set()
                for scS1, svS1 in subcomps1_corrected.items():
                    for scS2, svS2 in subcomps2_corrected.items():
                        if   svS2.issubset(svS1) and not svS1.issubset(svS2):
                            bad2.add(scS2)
                        elif svS1.issubset(svS2) and not svS2.issubset(svS1):
                            bad1.add(scS1)
                subcomps1_corrected = {c: vSs for c, vSs in subcomps1_corrected.items() if c not in bad1}
                subcomps2_corrected = {c: vSs for c, vSs in subcomps2_corrected.items() if c not in bad2}

                # Keep going
                subcomps1 = subcomps1_corrected
                subcomps2 = subcomps2_corrected

                subcomp2vertex1 = {scS1: {v for vS in svS1 for v in vertex2path[vS]} for scS1, svS1 in subcomps1.items()} # so we don't have to compute this inside the loop
                subcomp2vertex2 = {scS2: {v for vS in svS2 for v in vertex2path[vS]} for scS2, svS2 in subcomps2.items()}   

                # Assertions
                v1 = set.union(*list(subcomp2vertex1.values()))
                v2 = set.union(*list(subcomp2vertex2.values()))
                assert v1 == v2

                for svS1, svS2 in combinations(subcomps1.values(), 2):
                    assert not svS1 & svS2 # if they are not RC they should not share any paths
                for svS1, svS2 in combinations(subcomps2.values(), 2):
                    assert not svS1 & svS2 # if they are not RC they should not share any paths

                
                for scS1, svS1 in subcomps1.items():
                    for cS2, vSs in comp2vertexGS.items(): # they should not share any paths with other components in the sequence graph either
                        if cS2 != cS:
                            assert not vSs & svS1
                for scS2, svS2 in subcomps1.items():
                    for cS2, vSs in comp2vertexGS.items(): # they should not share any paths with other components in the sequence graph either
                        if cS2 != cS:
                            assert not vSs & svS2

                GS.clear_filters()                

                first = True
                for scS1, svS1 in subcomps1.items():
                    v1 = subcomp2vertex1[scS1]
                    for scS2, svS2 in subcomps2.items():
                        v2 = subcomp2vertex2[scS2]
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

            # Check that we included all the input vertices
            assert vs == {v for p in paths for v in p}
            # Check that no paths appear in two different scaffolds
            for cs1, cs2 in combinations(compS2paths, 2):
                assert not compS2paths[cs1] & compS2paths[cs2]

            ### Process each scaffold in the sequence graph
            for scaffold, paths in compS2paths.items():
                # Test that the scaffolds are really connected
                GS.clear_filters() # NEVER FORGET!!!
                self.set_vertex_filter(GS, {path2vertexGS[p] for p in paths})
                msg = f'\tScaffold {scaffold}, {GS.num_vertices()} seq vertices, {GS.num_edges()} seq edges'
                print_time(msg, end = '\r')
                assert len(set(gt.topology.label_components(GS, directed = False)[0])) == 1

                # Compute scaffold length and skip if too short
                scaffoldLen = 0
                for p in paths:
                    seq = path2seq[p]
                    if GS.get_out_neighbors(path2vertexGS[p]).shape[0]:
                        seq = seq[:-self.ksize] # trim the overlap unless this is a terminal node
                    scaffoldLen += len(seq)
                msg += f', {scaffoldLen} bases'
                if minlen and scaffoldLen < minlen:
                    msg += ' (< minlen, IGNORED)'
                    print_time(msg)
                    continue

                # Compute scaffold coverage
                scaffoldCov = np.mean([len(self.vertex2origins(v)) for p in paths for v in p])
                msg += f', cov {round(scaffoldCov, 2)}'
                if mincov and scaffoldCov < mincov:
                    msg += ' (< mincov, IGNORED)'
                    print_time(msg)
                    continue
                
                print_time(msg, end='\r')
                
                # Identify and flatten bubbles
                bubble_paths = set()
                bubble_vertex2origins = defaultdict(set) # I don't want to update self.vertex2origins since it's created in __init__ and this is run. Stupid, I know...
                if bubble_identity_threshold and bubble_identity_threshold < 1:
                    path_limits = defaultdict(list)
                    for p in paths:
                        path_limits[(p[0], p[-1])].append(p)    
                    m1 = 0
                    for (start, end), ps in path_limits.items():
                        ps = sorted(ps, key = lambda x: len(x), reverse = True)
                        minpathlen = len(ps[0])
                        if len(ps) > 1: # more than one path with similar start and end
                            p1 = ps[0] # the longest one
                            s1 = path2seq[p1]
                            ref = Aligner(seq = s1, n_threads = 1)
                            newOris = set()
                            for p2 in ps[1:]:                       # note that I am considering the 1st and last kmer here, which will always be shared!
                                al = list(ref.map(path2seq[p2]))[0] # this results in higher identities for short bubbles, which is not the worst
                                if al.mlen/al.blen >= bubble_identity_threshold: 
                                    bubble_paths.add(p2)  
                                    m1 += 1
                                    newOris.update( {ori for v2 in p2 for ori in self.vertex2origins(v2)} )
                            for v1 in p1:
                                bubble_vertex2origins[v1].update(newOris) # count vertices in p1 as appearing in all of the origins in the bubble paths that we just discarded
                                    
                    paths = [p for p in paths if p not in bubble_paths]
                    msg += f', removed {m1} bubbles'
                    print_time(msg, end = '\r')
                    GS.clear_filters()
                    self.set_vertex_filter(GS, {path2vertexGS[p] for p in paths})

                

                # Join contiguous unbranching paths (can happen after solving the fwd component and popping bubbles
                predecessors = defaultdict(set)
                successors = defaultdict(set)
                for p1 in paths:
                    for vS2 in GS.get_out_neighbors(path2vertexGS[p1]):
                        if vertex2path[vS2] not in bubble_paths:
                            successors[p1].add(vertex2path[vS2])
                    for vS2 in GS.get_in_neighbors(path2vertexGS[p1]):
                        if vertex2path[vS2] not in bubble_paths:
                            predecessors[p1].add(vertex2path[vS2])

                
                extenders = {p for p in paths if len(successors[p]) <= 1 and (successors[p] or predecessors[p])}
                extenders = extenders | {p for p in paths if len(predecessors[p]) == 1} ### added later, consider removing this line if we start failing assertions
                joinStarters =  {p for p in extenders if len(predecessors[p]) != 1 or len(successors[list(predecessors[p])[0]]) > 1}

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
                    succ = list(succs)[0] if succs else None
                    extended = False
                    while succ in joinExtenders:
                        extended = True
                        visited.add(succ)
                        new_path += succ[1:]
                        end = succ
                        succs = successors[succ]
                        succ = list(succs)[0] if succs else None
                    paths.append(new_path)
                    if extended:
                        sStart[new_path] = sStart[start]
                        sEnd[new_path] = sEnd[end]
                        path2seq[new_path] = self.reconstruct_sequence(new_path)
                        predecessors[new_path] = predecessors[start]
                        for pred in predecessors[new_path]:
                            successors[pred].add(new_path)
                        successors[new_path] = successors[end]
                        for succ in successors[new_path]:
                            predecessors[succ].add(new_path)
                        

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

                # By default we'll trim overlaps at 3', but we make exceptions
                # The idea is to remove also overlaps between to sequences sharing the same predecessor
                trim_left  = set()
                keep_right = set()
                for p in paths:
                    if len(successors[p]) > 1:
                        for succ in successors[p]:
                            if len(succ) > 1:
                                trim_left.add(succ)
                        keep_right.add(p)
                for p in paths:
                    if not successors[p] or len(p) < 2 or (len(p)<=2 and p in trim_left):
                        keep_right.add(p)
                    for succ in successors[p]:
                        if succ in trim_left:
                            keep_right.add(p)

                trimmedPaths = []
                for p in paths:
                    tp = tuple(p) # copy
                    if p not in keep_right:
                        tp = tp[:-1]
                    if p in trim_left:
                        tp = tp[1:]
                    assert tp
                    trimmedPaths.append(tp)

                vs = {v for p in paths for v in p}
                vst = {v for p in trimmedPaths for v in p}
                assert vs == vst
                
                trimmedSeqs = []
                for p in paths:
                    tseq = path2seq[p]
                    if p not in keep_right:
                        tseq = tseq[:-self.ksize]
                    if p in trim_left:
                        tseq = tseq[self.ksize:]
                    trimmedSeqs.append(tseq)
                

                # Update global results
                    
                path2id = {p: (scaffold, i) for i, p in enumerate(paths)}

                for p, tp, tseq in zip(paths, trimmedPaths, trimmedSeqs):
                    assert p not in addedPaths
                    id_ = path2id[p]
                    # Don't take into account overlaps with other paths for mean cov calculation and origin reporting
                    vertices = tp
                    ori2cov = defaultdict(int)
                    for v in vertices:
                        for ori in self.vertex2origins(v) | bubble_vertex2origins[v]: # note how I'm also getting origins not really associated to this path
                            ori2cov[ori] += 1                                         # but to bubbles that were originally in this path before we popped them
                    ori2cov   = {ori: prev/len(vertices) for ori, prev in ori2cov.items()}
                    # Even for complete genomes, we can have non-branching paths that belong to the core but don't have full coverage in some sequences
                    # This can happen if there are missing bases at the beginning or end of some of the input sequences (contigs in the original assemblies)
                    # So there is no branching, but some are missing
                    goodOris  = {ori for ori, cov in ori2cov.items() if cov >= genome_assignment_threshold}
                    goodOris.add( list(sorted(ori2cov, key = lambda ori: ori2cov[ori], reverse = True))[0] ) # add at least the origin with the highest cov
                    goodOris = {ori.split('_Nsplit_')[0] for ori in goodOris}
                    # Go on
                    contigs[id_] = Contig(scaffold   = id_[0],
                                          i          = id_[1],
                                          seq        = path2seq[p],
                                          tseq       = tseq,
                                          cov        = np.mean([len(self.vertex2origins(v) | bubble_vertex2origins[v]) for v in vertices]),
                                          origins    = goodOris,
                                          successors = [path2id[succ] for succ in successors[p]])
                    addedPaths.add(p)
                print_time(msg)

                    
        ### Finished!
        return contigs


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
        k0x = self.vertex2kmer(path[0])
        k1x = self.vertex2kmer(path[1]) 
        for k1 in (k1x, reverse_complement(k1x)):
            for k0 in (k0x, reverse_complement(k0x)):
                if k1[:-1] in k0:
                    kmers.append(k0)

        # Then just keep advancing in the path, for each vertex pick the kmer that overlaps with the previously-picked kmer
        for v in path[1:]:
            k = self.vertex2kmer(v)
            ks = (k, reverse_complement(k))
            toadd = []
            for k in ks:
                if k[:-1] == kmers[-1][1:]:
                    toadd.append(k)
            assert len(toadd) <= 2
            if not toadd:
                break
            kmers.extend(toadd)
        return self.seq_from_kmers(kmers)


    def vertex2origins(self, v):
        return {name for name, vs in self.seqPaths.items() if v in vs}


    def vertex2kmer(self, v):
        name, i = self.vertex2source[v]
        return self.seqDict[name][i:i+self.ksize]
            

        
    @staticmethod
    def seq_from_kmers(kmers):
        fullseq = [kmers[0]]
        extra = [kmer[-1] for kmer in kmers[1:]]
        fullseq.extend(extra)
        return ''.join(fullseq)

        

