from collections import defaultdict, namedtuple
from itertools import combinations
from multiprocessing import Pool
from multiprocessing.managers import SharedMemoryManager
from functools import partial

import numpy as np

from superpang.lib.utils import read_fasta, write_fasta, print_time
from superpang.lib.cutils import reverse_complement
from superpang.lib.vtools import compress_vertices, vertex_overlap, isInSeqPath
from superpang.lib.Compressor import Compressor

import graph_tool as gt
from graph_tool.all import Graph

from mappy import Aligner
from speedict import Rdict, Options, SliceTransform, PlainTableFactoryOptions

### see if we can put NBPG.clear_filters() inside the set filters method, instead of having a lot of them throughout the code

class Assembler:

    ### MULTIPROCESS METHODS
    # In order to save memory minimize serialization/deserialization overheads when multiprocessing we will pass the
    #   input variables to our target function as class attributes (a.k.a. gloryfied globals).
    #
    # We set them before opening the process pool: this way they get inherited when forking the original process,
    #   and linux's Copy On Write makes it so memory is actually not copied unless it is modified.
    #
    # The output of our function still needs to be serialized by the workers and deserialized by the main process,
    #   so this will not work if outputs are large and frequent
    #
    # For that we can use shared memory (see __init__), particularly for output amenable to be stored in
    #   an array (ints, bytes...)

    multiprocessing_globals = tuple()

    @classmethod
    def set_multiprocessing_globals(cls, *args):
        # Multiprocessing.Pool.map comes with a significant overhead in serializing the arguments and sending them to the workers,
        #    which beats the purpose of multiprocessing.
        #
        # What we'll do instead is store the required variables as a class attribute before starting the pool. Those will be copied
        #    when forking the process, which is much faster.
        cls.multiprocessing_globals = args


    @classmethod
    def clear_multiprocessing_globals(cls):
        cls.multiprocessing_globals = tuple()


    @classmethod
    def get_multiprocessing_globals(cls):
        return cls.multiprocessing_globals if len(cls.multiprocessing_globals) > 1 else cls.multiprocessing_globals[0]


    @classmethod
    def multimap(cls, fun, threads, iterable, *args, single_thread_threshold = 24, imap = False, pool = None):
        args = [iterable] + list(args)
        cls.set_multiprocessing_globals(*args)
        if threads == 1 or len(iterable) < single_thread_threshold:
            res = list(map(fun, range(len(iterable))))
        else:
            with Pool(threads) as pool:
                res = pool.map(fun, range(len(iterable)))
        cls.clear_multiprocessing_globals()
        return res
    ###############################



##    @profile
    def __init__(self, fasta, ksize, threads, diskdb = None, debug = False):
        """
        Build a De-Bruijn Graph from an input fasta file
        """

        self.ksize        = ksize
        self.includedSeqs = 0
        self.seqPaths     = {} # name: compressedVertices
        self.seqLimits    = set()
        MAX_KMER_CUTOFF   = 1000000000

        ### Read input sequences
        print_time('Reading sequences')
        self.seqDict = read_fasta(fasta, ambigs = 'as_Ns', Ns = 'split')
        nNs = len({name.split('_Nsplit_')[0] for name in self.seqDict})
        self.ref2name = list(self.seqDict)
        max_kmers = sum(len(seq) - self.ksize + 1 for seq in self.seqDict.values() if len(seq) > self.ksize)
        max_kmers = max_kmers if max_kmers < MAX_KMER_CUTOFF else MAX_KMER_CUTOFF

        print_time('Allocating resources')
        ### Kmer hashing related parameters
        self.compressor = Compressor(self.ksize)
        clength = self.ksize // 4 + self.ksize % 4                # the length of a hash
        maxseqlength = max(len(s) for s in self.seqDict.values()) # length of the largest input sequence
        maxnkmers = maxseqlength-self.ksize+1                     # number of kmers in the largest input sequence

        ### Set up kmer storage
        self.vertex2coords = np.empty(shape = (max_kmers, 2), dtype = np.uint32)
        if diskdb:
            opts = Options()
            opts.increase_parallelism(threads)
            opts.set_max_background_jobs(threads)
            opts.set_prefix_extractor(SliceTransform.create_max_len_prefix(8))
            opts.set_plain_table_factory(PlainTableFactoryOptions())
            hash2vertex = Rdict(path = diskdb, options = opts)
        else:
            hash2vertex = {}

        ### Set up edge storage
        maxint    = np.iinfo(np.uint32).max
        edges     = np.empty(shape = (max_kmers, 2), dtype = np.uint32)
        edges[:]  = maxint

        ### Go for the eyes, boo!
        print_time(f'Creating DBG from {len(self.seqDict)} sequences ({nNs} originals)')
        includedSeqs   = 0
        self.seqLimits = set()
        totalSize      = sum(len(seq) for seq in self.seqDict.values())
        elapsedSize    = 0
        lastSize       = 0
        nVertices      = np.uint32(0)
        one            = np.uint32(1) # This way increasing the nVertices counter may be slightly faster
        nEdges         = 0 # actual edges added so far (includes duplicated edges)
        nEdges_ndEst   = 0 # estimate of non-duplicated edges added

        with SharedMemoryManager() as smm:
            # seqMem and rcSeqMem are two shared mem buffers that will hold the fwd and rev sequences
            #    they are pre-allocated using the size of the largest sequence, for each sequence we will
            #    store/recover the correct amount of bytes according to its length
            # hashMem and revhashMem are two continuous shared mem buffers that will hold all the hashes concatenated
            #    they are pre-allocated using the number of kmers in the largest sequence multiplied by the hash size
            #    to store/recover a particular hash we can just use the right slice of the mem buffer, since hash size is constant

            seqMem     = smm.SharedMemory(size=maxseqlength)
            rcSeqMem   = smm.SharedMemory(size=maxseqlength)
            hashMem    = smm.SharedMemory(size=clength*maxnkmers)
            self.set_multiprocessing_globals(seqMem, rcSeqMem, hashMem, self.ksize, clength, self.compressor)

            with Pool(threads) as pool:
 
                for ref, (name, seq) in enumerate(self.seqDict.items()):
                    if len(seq) <= self.ksize:
                        elapsedSize += len(seq)
                        continue
                    else:
                        includedSeqs += 1
                        rcSeq = reverse_complement(seq)
                        seq, rcSeq = seq.encode(), rcSeq.encode()
                        nkmers = len(seq)-self.ksize+1
                        elapsedSize += len(seq) - nkmers

                        seqMem.buf  [:len(seq)  ] = seq
                        rcSeqMem.buf[:len(rcSeq)] = rcSeq

                        seq2hashes_ = partial(self.seq2hashes, seqlength=len(seq)) # we can't pass seqlength as a global when forking
                                                                                   #  since it changes for every sequence
                        pool.map(seq2hashes_, range(nkmers))

                        # Retrieve the kmer hashes from the memory buffers
                        # Note the use of a generator expression so that all the hashes are not residing in memory at the same time
                        #  ... so we can't close the pool until we finish going through such generator expression!
                        hashes    = (hashMem.buf[i*clength:(i+1)*clength].tobytes() for i in range(nkmers))
                        # Populate the DBG
                        sp = np.empty(nkmers, dtype = np.uint32)
                        prev_v, prev_new, past_first = -1, False, False # store the previous vertex here so we can build (prev, v) edge tuples
                        
                        for idx, h in enumerate(hashes):
                            if h in hash2vertex:
                                newVertex = False
                                v = hash2vertex[h]
                            else:
                                newVertex = True
                                v = nVertices
                                hash2vertex[h] = v
                                self.vertex2coords[v] = (ref, idx)
                                nVertices += one

                            sp[idx] = v
                            
                            if past_first:
                                edges[nEdges] = prev_v, v
                                nEdges += 1
                                if (newVertex or prev_new): # Just so we have an ongoing estimate of the unique edges added so far
                                    nEdges_ndEst += 1       # could remove this for extra efficiency

                            prev_v, prev_new, past_first = v, newVertex, True
                            if idx == 0:
                                self.seqLimits.add(v)
                            if idx == nkmers - 1:
                                self.seqLimits.add(v)

                            elapsedSize += 1
                            
                            if elapsedSize - lastSize > totalSize/100:
                                print_time(f'\t{round(100*elapsedSize/totalSize, 2)}% bases processed, {includedSeqs} sequences, {nVertices} vertices, {nEdges_ndEst} edges         ', end = '\r')
                                lastSize = elapsedSize

                        
                        self.seqPaths[name] = compress_vertices(sp)
                        if max_kmers == MAX_KMER_CUTOFF and nEdges / max_kmers > 0.75: # Pack the edge array by removing duplicates if we are size-capped and close to filling it
                            edges  = edges[:nEdges]
                            edges.sort(axis=1)                     # Sort each row so duplicated edges in opposite directions (e.g. 2,3 and 3,2) become identical (2,3 and 2,3)
                            edges  = np.unique(edges, axis=0)      # Deduplicate
                            nEdges = nEdges_ndEst = edges.shape[0] # Adjust to the number of deduplicated edges
                            ext    = np.empty(shape = (max_kmers - nEdges, 2), dtype = np.uint32)              
                            ext[:] = maxint
                            edges  = np.vstack((edges, ext))       # Add empty space up to max_kmers  
                        if debug:
                            assert idx+1 == nkmers
                            assert edges.shape[0] == max_kmers

        ### Populate DBG and cleanup
        del seqMem, rcSeqMem, hashMem
        if diskdb:
            hash2vertex.close()
            Rdict.destroy(diskdb, opts)
        del hash2vertex
        self.clear_multiprocessing_globals()

        if debug:
            for edge in edges[:nEdges]:
                if maxint in edge:
                    assert False

        self.vertex2coords = self.vertex2coords[:nVertices]
        edges = edges[:nEdges]
        edges.sort(axis=1)               # Sort each row so duplicated edges in opposite directions (e.g. 2,3 and 3,2) become identical (2,3 and 2,3)
        edges = np.unique(edges, axis=0) # Deduplicate
        
        self.DBG = Graph(directed = False)
        self.DBG.add_edge_list(edges)
        del edges
        self.DBG.vp.comp  = self.DBG.new_vertex_property('int16_t')
        self.DBG.ep.efilt = self.DBG.new_edge_property('bool')
        print_time(f'\t100% bases processed, {includedSeqs} sequences, {self.DBG.num_vertices()} vertices, {self.DBG.num_edges()} edges         ')


##    @profile
    def run(self, minlen, mincov, bubble_identity_threshold, genome_assignment_threshold, max_threads, debug = False):
        """
        Create contigs from a De-Bruijn Graph
        """

        Contig = namedtuple('Contig', ['scaffold', 'i', 'seq', 'tseq', 'cov', 'origins', 'successors'])
        contigs = {}
        addedPaths = set()

        ### Identify connected components
        print_time('Identifying connected components')
        self.DBG.vp.comp = gt.topology.label_components(self.DBG, directed = False)[0]
        comps, sizes = np.unique(self.DBG.vp.comp.get_array(), return_counts=True) # Process from larger to smaller
        comps = [comp for comp, size in sorted(zip(comps, sizes), key = lambda x: x[1], reverse = True)]
        nComps = len(comps)

        currentScaffold = 0

        ### Process connected components
        for ic, c in enumerate(comps):
            if False: # _debug CODE: Skip unwanted components
                TGT_COMP = 1000
                if ic + 1 != TGT_COMP:
                    continue
                print(f'Comp{ic+1}, c={c}')
            vs = np.where(self.DBG.vp.comp.get_array() == np.int16(c))[0]
            print_time(f'Working on comp {ic+1}/{nComps}, {len(vs)} vertices')

            # The overhead of opening a multiprocessing pool increases with resident memory, even if we do nothing inside it
            #  so I don't think this is due to COW being triggered, just something the way fork and/or python multiprocessing pool works
            #  This becomes a problem when I have lots of small components in an otherwise large project, as I pay a large overhead several
            #  times per component. So instead we will use one thread for small-ish components.
            threads = max_threads if len(vs) > 150000 else 1

            # Get the seqPaths that appear in this DBG component
            #  just check that the first vertex in the seqPath belongs to this component
            seqPathsThisComp = {name: sp for name, sp in self.seqPaths.items() if sp[0][0] in vs}
            ### Reconstruct non-branching paths
            print_time('\tCollecting non-branching paths')
            NBPs = []
            inits = set()
            badEdges = set()
            inits = {v for v, isInit in zip(vs, self.multimap(self.is_NBP_init, threads, vs, self.DBG, self.seqLimits,
                                                                     self.vertex2coords, self.ref2name, self.seqDict, self.ksize)) if isInit}
            # Break cycle if required
            isCycle = False
            ignoreEdge = set()
            if not inits:
                #is this a cycle?
                if min([len(set(self.DBG.get_out_neighbors(v))) for v in vs]) == 2: # no sources or sinks
                    for v in vs:
                        neighbors = self.DBG.get_out_neighbors(v)
                        if len(set(neighbors)) == 2: # find the first extender edge, and break there
                            inits.add(v)
                            inits.add(neighbors[0])
                            isCycle = True
                            ignoreEdge = {v, neighbors[0]}
                            badEdges.add(self.DBG.edge(v, neighbors[0]))
                            badEdges.add(self.DBG.edge(neighbors[0], v))
                            break

            # Filter out bad edges
            if badEdges:
                for e in self.DBG.edges():
                    self.DBG.ep.efilt[e] = e not in badEdges
                self.DBG.set_edge_filter(self.DBG.ep.efilt)
                    
            for j, ini in enumerate(inits):
                for n in set(sorted(self.DBG.get_out_neighbors(ini))):
                    NBP = [ini, n]
                    last = ini
                    while True:
                        if NBP[-1] in inits:
                            break
                        s = [x for x in self.DBG.get_out_neighbors(NBP[-1]) if x != last and x != NBP[-1]][0] #NBP[-1] to avoid stopping at self loops
                        last = NBP[-1]
                        NBP.append(s) # before we casted to np.uint32, but it doesn't seem to make a big difference memory wise and takes a bit.
                        #NBP.append(np.uint32(s))
                    NBP = tuple(NBP)
                    NBPs.append(NBP)

            if isCycle: # make them circular by adding the first vertex also in the end so there's an overlap
                assert len(NBPs) == 2
                NBPs = [NBP+(NBP[0],) for NBP in NBPs]

            # Check that there are no len==1 NBPs
            for NBP in NBPs:
                assert len(NBP) > 1

            ### Build a graph of connected non-branching paths, and split into fwd and rev component
            print_time(f'\tBuilding sequence graph out of {len(NBPs)} non-branching paths')

            # Remove NBPs with terminally duplicated vertices
            # Those broke our code once. we had three nvs with NBPs (0,0,1,2) (0,1,2) and (2,1,0) under the same pset (see below)
            NBPs = [NBP for NBP in NBPs if NBP[0] != NBP[1] and NBP[-1] != NBP[-2]]

            # Translate NBPs into sequences
            NBP2seq = dict(zip(NBPs, self.multimap(self.reconstruct_sequence, threads, NBPs, self.vertex2coords, self.ref2name, self.seqDict, self.ksize)))
            # Sort based on the sequence for extra reproducibility (e.g. when the order of the input sequences is different)
            # Right now (20 Jul 23) it seems this actually made no difference except for the order in which bubble representatives are selected
            #  and I modified the bubble sorting code so sorting here makes no difference at all as far as I could test
            # But it is cheap to do anyways so I'll keep it
            NBP2seq = dict(sorted(NBP2seq.items(), key = lambda x: (len(x[1]), x[1]), reverse = True)) # Assuming that python dicts preserve insertion ordering, which they should nowadays
            NBPs = list(NBP2seq)

            # Collect prefixes-suffixes for locating potential overlaps
            start2NBPs = defaultdict(set) # first k nucleotides
            sStart = {}
            sEnd = {}
            for NBP, seq in NBP2seq.items():
                start, end = seq[:self.ksize], seq[-self.ksize:]
                start2NBPs[start].add(NBP)
                sStart[NBP] = start
                sEnd[NBP] = end

            # Build graph of non-branching paths
            NBPG = Graph()
            NBP2vertex = {}

            for NBP in NBPs:
                assert NBP not in NBP2vertex
                NBP2vertex[NBP] = int(NBPG.add_vertex())
            for NBP1 in NBPs:
                for NBP2 in start2NBPs[sEnd[NBP1]]:
                    NBPG.add_edge(NBP2vertex[NBP1], NBP2vertex[NBP2])
                           
            vertex2NBP = {v: NBP for NBP,v in NBP2vertex.items()}
                
            ### Write component's input sequences, DBG and NBPG
            if False:
                oris = {ori for oris in self.multimap(self.get_vertex2origins, threads, list(vs), seqPathsThisComp, single_thread_threshold = 10000) for ori in oris}
                with open(f'comp_{c+1}.fasta', 'w') as outfile:
                    for ori in oris:
                        outfile.write(f'>{ori}\n{self.seqDict[ori]}\n')
                self.DBG.vp.vfilt = self.DBG.new_vertex_property('bool')
                self.set_vertex_filter(self.DBG, set(vs))
                DBG2 = Graph(self.DBG, prune = True)
                print(f'DBG: {DBG2.num_vertices()} vertices, {DBG2.num_edges()} edges')
                DBG2.save(f'comp_{c+1}.DBG.graphml')
                print(f'NBPG: {NBPG.num_vertices()} vertices, {NBPG.num_edges()} edges')
                NBPG.save(f'comp_{c+1}.NBPG.graphml')
            del vs


            ### Extract forward component
            def get_psets(nvs):
                """Get a dict with a frozenset of vertices as keys, and the two complementary non-branching paths sharing those vertices as values"""
                psets = defaultdict(list)
                for nv in nvs:
                    pset = frozenset(vertex2NBP[nv]) # Tried to use hash(frozenset(...)) as keys before
                    psets[pset].append(nv)           #  but led to hash collisions and the assertion below failing
                for nvs_ in psets.values():
                    assert len(nvs_) == 2
                return psets


            print_time('\tExtracting forward component')
            comp2nvs = defaultdict(set)
            comps =  gt.topology.label_components(NBPG, directed = False)[0]
            NBPG.vp.comps = comps
            NBPG.vp.vfilt = NBPG.new_vertex_property('bool')
            NBPG.ep.efilt = NBPG.new_edge_property('bool')

            for nv in NBPG.vertices():
                comp2nvs[NBPG.vp.comps[nv]].add(int(nv))

            # Try to find separated components in the seq graph that contain the same vertices in the kmer graph (thus are RC of each other)
            rcComps = {}
            for nc1, nc2 in combinations(comp2nvs,2):
                v1 = {v for nv in comp2nvs[nc1] for v in vertex2NBP[nv]}
                v2 = {v for nv in comp2nvs[nc2] for v in vertex2NBP[nv]}
                if v1 == v2:
                    rcComps[nc1] = nc2
                    rcComps[nc2] = nc1
                del v1, v2

            # Count how many times each Non-Branching Path is in the same orientation as the input sequences
            nv2rightOrientation = defaultdict(int)
            added = set()
            for nc1, nc2 in rcComps.items(): # For now we only do it for the components with a clear RC
                if nc1 in added:
                    continue
                # Find reverse complement Non-Branching Paths
                psets = get_psets(comp2nvs[nc1] | comp2nvs[nc2])
                assert len(psets) == len(comp2nvs[nc1]) == len(comp2nvs[nc2]) # each sequence path should have a reverse complement equivalent (same vertices in the kmer graph, reverse order)
                # Count
                names = list(seqPathsThisComp.keys())
                stt = threads if len(psets) > threads else 1000000000
                seqPaths_NBPs = dict(zip(names, self.multimap(self.get_seqPaths_NBPs, threads, names, seqPathsThisComp, psets, NBP2seq, vertex2NBP, self.seqDict)))
                seqPaths_NBPs = {name: spGSs for name, spGSs in seqPaths_NBPs.items() if spGSs}
                for spGS in seqPaths_NBPs.values():
                    for nv in spGS:
                        nv2rightOrientation[nv] += 1
                added.add(nc1)
                added.add(nc2)


            # If we didn't find a RC for all components this means that in some case fwd and rev are joined in the same component
            noRC = {nc for nc in comp2nvs if nc not in rcComps}
            for nc in noRC: # for each of those evil components...

                # Find reverse complement Non-Branching Paths
                psets = get_psets(comp2nvs[nc])             
                assert len(psets) == len(comp2nvs[nc]) / 2 # each sequence path should have a reverse complement equivalent (same vertices in the kmer graph, reverse order)
                nv2rc = {}
                for nv1, nv2 in psets.values():
                    nv2rc[nv1] = nv2
                    nv2rc[nv2] = nv1
                self.set_vertex_filter(NBPG, comp2nvs[nc])

                if debug:
                    assert {int(nv) for nv in NBPG.vertices()} == comp2nvs[nc] # did we alter nv indices when filtering?
                    assert(len(set(gt.topology.label_components(NBPG, directed = False)[0]))) == 1 # check that this is only one component

                    for pset, (nv1, nv2) in psets.items():
                        assert vertex2NBP[nv1] == vertex2NBP[nv2][::-1]

                # Map the input sequences to nodes in the sequence graph (in the original orientation of the input sequence)
                names = list(seqPathsThisComp.keys())
                seqPaths_NBPs = dict(zip(names, self.multimap(self.get_seqPaths_NBPs, threads, names, seqPathsThisComp, psets, NBP2seq, vertex2NBP, self.seqDict)))
                seqPaths_NBPs = {name: spGSs for name, spGSs in seqPaths_NBPs.items() if spGSs}

                # Get the times each node in the sequence graph mapped to an input sequence in the original orientation
                addedNodes = set()
                for spGS in seqPaths_NBPs.values():
                    for nv in spGS:
                        nv2rightOrientation[nv] += 1
                        addedNodes.update( (nv, nv2rc[nv]) )
                    
                # Add nodes in the sequence graph that weren't fully contained in an original sequence
                for nv in comp2nvs[nc] - addedNodes:
                    fakeName = hash( (len(seqPaths_NBPs), nv) )
                    seqPaths_NBPs[fakeName] = {nv} | set(NBPG.get_all_neighbors(nv)) # add its neighbors so that it'll overlap with the real sequences in at least one nv

                seqPaths_NBPs_rev = {name: {nv2rc[nv] for nv in nvs} for name, nvs in seqPaths_NBPs.items()}
                
                # Separate fwd and rev
                nvs1 = set()
                refSeqNames = sorted(seqPaths_NBPs, key = lambda name: len(seqPaths_NBPs[name]), reverse = True) # sort by decreasing number of nodes in the NBPG

                addedNames = {refSeqNames[0]}
                nvs1.update(seqPaths_NBPs[refSeqNames[0]]) # start from the longest seq

                while True:
                    remaining = [n for n in refSeqNames if n not in addedNames] # keep them sorted
                    if not remaining:
                        break
                    best_name = ''
                    best_ol = -1 # -1 so that even if there are no overlaps (there is more than one subcomponent, and the current one is full) we add something
                    for name in remaining: # keep starting from the longest seq, so that if there are no overlaps we add the first longest one
                        fwd = seqPaths_NBPs[name]
                        rev = seqPaths_NBPs_rev[name]
                        fwd_ol = len(fwd & nvs1) / len(fwd) # find the best overlap (looking in both orientations) to the existing fwd component
                        rev_ol = len(rev & nvs1) / len(rev)
                        if fwd_ol > best_ol and fwd_ol >= rev_ol:
                            best_ol, best_name, best_nvs = fwd_ol, name, fwd
                        elif rev_ol > best_ol and rev_ol > fwd_ol:
                            best_ol, best_name, best_nvs = rev_ol, name, rev
                    assert best_name
                    addedNames.add(best_name)
                    best_nvs = {nv for nv in best_nvs if nv2rc[nv] not in nvs1} # don't add reverse complements at this stage
                                                                                # some are really there and needed to keep everything into a single component
                                                                                # some will be duplications (e.g. because one of the input genomes has an inverted region)
                    nvs1.update(best_nvs)

                nvs2 = {nv2rc[nv] for nv in nvs1}

                assert (nvs1 | nvs2) == comp2nvs[nc]

                # Identify subcomponents in nvs1
                NBPG.clear_filters()
                self.set_vertex_filter(NBPG, nvs1)
                subcomps1 = defaultdict(set)
                for nv,snc1 in zip(NBPG.vertices(), gt.topology.label_components(NBPG, directed = False)[0]):
                    subcomps1[snc1].add(nv)
                NBPG.clear_filters()
                # Identify subcomponents in nvs2
                self.set_vertex_filter(NBPG, nvs2)
                subcomps2 = defaultdict(set)
                for nv,snc2 in zip(NBPG.vertices(), gt.topology.label_components(NBPG, directed = False)[0]):
                    subcomps2[snc2].add(nv)
                NBPG.clear_filters()
                # Identify RC components
                assert len(subcomps1) == len(subcomps2)

                
                # Try to merge into one fwd and one rev subcomp
                # In some cases the fwd and rev sequences can be actually be present in the same genome
                # So we identify the largest fwd subcomponent, throw in the rest of components in fwd and rev, and hope that they become connected when everything is together
                largestSubcomp = list(sorted(subcomps1, key = lambda snc: len(subcomps1[snc]), reverse = True))[0]
                combined = subcomps1[largestSubcomp]
                for snc in subcomps1:
                    if snc != largestSubcomp:
                        fwd = subcomps1[snc]
                        rev = {nv2rc[nv] for nv in fwd}
                        combined = combined | fwd | rev

                # We recalculate the components after adding the rev subcomponents
                subcomps1_corrected = defaultdict(set)
                subcomps2_corrected = defaultdict(set)
                self.set_vertex_filter(NBPG, combined)
                for nv, snc  in zip(NBPG.vertices(), gt.topology.label_components(NBPG, directed = False)[0]):
                    subcomps1_corrected[snc].add(nv)
                    subcomps2_corrected[snc].add(nv2rc[nv])


                # The lines above may make a rev subcomponent may end up being contained in a fwd component
                # E.g. let us have initially three components A, B, C (fwd) and A', B', C' (rev)
                # We correct and it turns out that B' got inside A so by the code above we end up with AB', B, C (fwd) and A'B, B', C' (rev)
                # So we want to remove B from fwd and B' from rev, so that we don't end up having duplicate paths
                bad1 = set()
                bad2 = set()
                for snc1, snv1 in subcomps1_corrected.items():
                    for snc2, snv2 in subcomps2_corrected.items():
                        if   snv2.issubset(snv1) and not snv1.issubset(snv2):
                            bad2.add(snc2)
                        elif snv1.issubset(snv2) and not snv2.issubset(snv1):
                            bad1.add(snc1)
                subcomps1_corrected = {snc: nvs for snc, nvs in subcomps1_corrected.items() if snc not in bad1}
                subcomps2_corrected = {snc: nvs for snc, nvs in subcomps2_corrected.items() if snc not in bad2}

                # Keep going
                subcomps1 = subcomps1_corrected
                subcomps2 = subcomps2_corrected


                # Assertions
                if debug:
                    subcomp2vertex1 = {snc1: {v for nv in snv1 for v in vertex2NBP[nv]} for snc1, snv1 in subcomps1.items()}
                    subcomp2vertex2 = {snc2: {v for nv in snv2 for v in vertex2NBP[nv]} for snc2, snv2 in subcomps2.items()}
                    v1 = set.union(*list(subcomp2vertex1.values()))
                    v2 = set.union(*list(subcomp2vertex2.values()))
                    assert v1 == v2
                    del subcomp2vertex1, subcomp2vertex2

                    for snv1, snv2 in combinations(subcomps1.values(), 2):
                        assert not snv1 & snv2 # if they are not RC they should not share any NBPs or have the same vertices
                        assert not {v for nv in snv1 for v in vertex2NBP[nv]} == {v for nv in snv2 for v in vertex2NBP[nv]}
                    for snv1, snv2 in combinations(subcomps2.values(), 2):
                        assert not snv1 & snv2 # if they are not RC they should not share any NBPs or have the same vertices
                        assert not {v for nv in snv1 for v in vertex2NBP[nv]} == {v for nv in snv2 for v in vertex2NBP[nv]}
                    for snc1, snv1 in subcomps1.items():
                        for nc2, nvs in comp2nvs.items(): # they should not share any paths with other components in the sequence graph either
                            if nc2 != nc:
                                assert not nvs & snv1
                    for snc2, snv2 in subcomps1.items():
                        for nc2, nvs in comp2nvs.items(): # they should not share any paths with other components in the sequence graph either
                            if nc2 != nc:
                                assert not nvs & snv2

                NBPG.clear_filters()                

                first = True
                for snv1 in subcomps1.values():
                    v1 = sorted([v for nv in snv1 for v in vertex2NBP[nv]]) # list instead of set to reduce memory usage
                    i1 = nc if first else len(comp2nvs) # recycle the previous component index for nc1
                    comp2nvs[i1] = snv1
                    first = False
                    for snv2 in subcomps2.values():                         #  but ofc it sucks
                        v2 = sorted([v for nv in snv2 for v in vertex2NBP[nv]])
                        if v1 == v2:
                            i2 = len(comp2nvs)
                            comp2nvs[i2] = snv2
                            rcComps[i1] = i2
                            rcComps[i2] = i1
                            break
                    else: # assert that all our subcomponents have RC components
                        assert False

                del v1, v2
                NBPG.clear_filters()

            assert len(rcComps) == len(comp2nvs) # all our components have a rc component, splitting was successful

            # For each pair of RC components, we select the component with the most Non-Branching paths in the same orientation as the input sequences
            compS2paths = {}
            added = set()
            for nc1, nc2 in rcComps.items():
                if nc1 in added:
                    continue
                nvsRight1 = sum(nv2rightOrientation[nv] for nv in comp2nvs[nc1])
                nvsRight2 = sum(nv2rightOrientation[nv] for nv in comp2nvs[nc2])
                if nvsRight1 > nvsRight2:
                    fwd = nc1
                elif nvsRight1 < nvsRight2:
                    fwd = nc2
                else: # nvsRight1 == nvsRight2; break the tie with ATGs
                    atg1 = sum(NBP2seq[vertex2NBP[nv]].count('ATG') for nv in comp2nvs[nc1])
                    atg2 = sum(NBP2seq[vertex2NBP[nv]].count('ATG') for nv in comp2nvs[nc2])
                    if atg1 > atg2:
                        fwd = nc1
                    elif atg1 < atg2:
                        fwd = nc2
                    else: # break the tie howewer
                        ord1 = sum(ord(b) for nv in comp2nvs[nc1] for b in NBP2seq[vertex2NBP[nv]])
                        ord2 = sum(ord(b) for nv in comp2nvs[nc2] for b in NBP2seq[vertex2NBP[nv]])
                        fwd = nc1 if ord1 >= ord2 else nc2
                compS2paths[currentScaffold] =  [vertex2NBP[nv] for nv in comp2nvs[fwd]]
                added.add(nc1)
                added.add(nc2)
                currentScaffold += 1

            if debug:
                # Check that we included all the input vertices
                vs = np.where(self.DBG.vp.comp.get_array() == np.int16(c))[0]
                assert set(vs) == {v for NBP in NBPs for v in NBP}
                # Check that no paths appear in two different scaffolds
                for cs1, cs2 in combinations(compS2paths, 2):
                    assert not set(compS2paths[cs1]) & set(compS2paths[cs2])


            ### Process each scaffold in the sequence graph
            # First order by length for extra reproducibility
            compS2length = defaultdict(int)
            for scaffold, NBPs in compS2paths.items():
                for NBP in NBPs:
                    seq = NBP2seq[NBP]
                    if NBPG.get_out_neighbors(NBP2vertex[NBP]).shape[0]:
                        seq = seq[:-self.ksize] # trim the overlap unless this is a terminal node
                    compS2length[scaffold] += len(seq)
            compS2length = dict( sorted(compS2length.items(), key = lambda x: x[1], reverse = True) )
            compS2paths = {scaffold: compS2paths[scaffold] for scaffold in compS2length}

            # Now go
            for scaffold, NBPs in compS2paths.items():

                # Test that the scaffolds are really connected
                NBPG.clear_filters() # NEVER FORGET!!!
                self.set_vertex_filter(NBPG, {NBP2vertex[NBP] for NBP in NBPs})
                msg = f'\tScaffold {scaffold}, {NBPG.num_vertices()} NBP vertices, {NBPG.num_edges()} NBP edges'
                print_time(msg, end = '\r')

                if debug:
                    assert len(set(gt.topology.label_components(NBPG, directed = False)[0])) == 1

                # Skip if too short
                scaffoldLen = compS2length[scaffold]
                msg += f', {scaffoldLen} bases'
                if minlen and scaffoldLen < minlen:
                    msg += ' (< minlen, IGNORED)'
                    print_time(msg)
                    continue

                if mincov:
                    # Compute scaffold coverage
                    scaffoldCov = np.mean(self.multimap(self.get_vertex2originslen, threads, [v for NBP in NBPs for v in NBP], seqPathsThisComp, single_thread_threshold = 10000))
                    msg += f', cov {round(scaffoldCov, 2)}'
                    if mincov and scaffoldCov < mincov:
                        msg += ' (< mincov, IGNORED)'
                        print_time(msg)
                        continue
                
                print_time(msg, end='\r')
                
                # Identify and flatten bubbles
                sortBy = 'length'         # can be 'length' or 'cov'.
                bubble_paths_all = list() # length had slighly higher assembly completeness when I tested it, and is faster
                bubble_paths     = set()
                bubble_vertex2origins = defaultdict(list) # Keep track of the new ones here,
                                                          #  the originals are calculated on the fly to save memory
                                                          #  use list instead of set to save memory
                if bubble_identity_threshold and bubble_identity_threshold < 1:
                    path_limits = defaultdict(list)
                    for NBP in NBPs:
                        path_limits[(NBP[0], NBP[-1])].append(NBP)
                    for (start, end), nbps in path_limits.items():
                        if len(nbps) > 1: # more than one path with similar start and end
                            bubble_paths_all.extend(nbps)
                    if sortBy == 'cov':
                        path2cov = dict(zip(bubble_paths_all, self.multimap(self.get_avgPathCov, threads, bubble_paths_all, seqPathsThisComp)))
                    m1 = 0
                    for (start, end), nbps in path_limits.items():
                        
                        if len(nbps) > 1: # more than one path with similar start and end
                            if sortBy == 'length': # We will also add the seq itself as a second ordering factor to break ties
                                nbps = sorted(nbps, key = lambda x: (len(x),      NBP2seq[x]), reverse = True) 
                            elif sortBy == 'cov':
                                nbps = sorted(nbps, key = lambda x: (path2cov[x], NBP2seq[x]), reverse = True)
                            else:
                                assert False
                            nbp1 = nbps[0] # the longest / most covered one
                            s1 = NBP2seq[nbp1][self.ksize:-self.ksize] # remove the first and last kmer
                            ref = Aligner(seq = s1, preset = 'sr', n_threads = 1)
                            newOris = set()
                            for nbp2 in nbps[1:]:
                                s2 = NBP2seq[nbp2][self.ksize:-self.ksize]
                                al = list(ref.map(s2)) 
                                if not al:
                                    continue
                                mlen = sum(a.mlen for a in al)
                                if mlen/len(s2) >= bubble_identity_threshold: # note that s2 might be smaller than s1 if sorting by cov, so I should correct this
                                    bubble_paths.add(nbp2)  
                                    m1 += 1
                                    newOris.update( {ori for v2 in nbp2 for ori in self.vertex2origins(v2, seqPathsThisComp)} )
                            for v1 in nbp1:
                                bubble_vertex2origins[v1].extend(newOris) # count vertices in p1 as appearing in all of the origins in the bubble paths that we just discarded
                           
                    NBPs = [NBP for NBP in NBPs if NBP not in bubble_paths]
                    msg += f', removed {m1} bubbles'
                    print_time(msg, end = '\r')
                    NBPG.clear_filters()
                    self.set_vertex_filter(NBPG, {NBP2vertex[NBP] for NBP in NBPs})

                

                # Join contiguous non-branching paths (can happen after solving the fwd component and popping bubbles
                predecessors = defaultdict(list)
                successors = defaultdict(list)
                for p1 in NBPs:
                    for nv2 in NBPG.get_out_neighbors(NBP2vertex[p1]):
                        if vertex2NBP[nv2] not in bubble_paths:
                            successors[p1].append(vertex2NBP[nv2])
                    for nv2 in NBPG.get_in_neighbors(NBP2vertex[p1]):
                        if vertex2NBP[nv2] not in bubble_paths:
                            predecessors[p1].append(vertex2NBP[nv2])
                
                extenders = {NBP for NBP in NBPs if len(successors[NBP]) <= 1 and (successors[NBP] or predecessors[NBP])}
                extenders = extenders | {NBP for NBP in NBPs if len(predecessors[NBP]) == 1} ### added later, consider removing this line if we start failing assertions
                joinStarters =  {NBP for NBP in extenders if len(predecessors[NBP]) != 1 or len(successors[list(predecessors[NBP])[0]]) > 1}

                if not joinStarters: # Is this a cycle
                    if min(len(ps) for ps in predecessors.values()) == 1: # This does not cover the case in which there is a cycle but also a linear component
                        joinStarters = {list(extenders)[0]}
                
                joinExtenders = {NBP for NBP in extenders if NBP not in joinStarters}
                joinSE = joinStarters | joinExtenders
                assert extenders == joinSE

                NBPs = [NBP for NBP in NBPs if NBP not in joinSE] # this also removes some joinStarters that can't be extended (no successors), but we'll add them back in the next bit


                newPaths = []
                
                visited = set()
                for start in joinStarters:
                    new_path = start
                    end = start
                    succs = successors[start]
                    succ = succs[0] if succs else None
                    extended = False
                    while succ in joinExtenders:
                        extended = True
                        visited.add(succ)
                        new_path += succ[1:]
                        end = succ
                        succs = successors[succ]
                        succ = succs[0] if succs else None
                    NBPs.append(new_path)
                    if extended:
                        sStart[new_path] = sStart[start]
                        sEnd[new_path] = sEnd[end]
                        newPaths.append(new_path)
                        predecessors[new_path] = predecessors[start]
                        for pred in predecessors[new_path]:
                            successors[pred].append(new_path)
                        successors[new_path] = successors[end]
                        for succ in successors[new_path]:
                            predecessors[succ].append(new_path)

                NBP2seq.update(dict(zip(newPaths, self.multimap(self.reconstruct_sequence, threads, newPaths, self.vertex2coords, self.ref2name, self.seqDict, self.ksize))))

                assert visited == joinExtenders

                # Update scaffold graph
                pathSet = set(NBPs)
                for NBP, preds in predecessors.items():
                    predecessors[NBP] = [pred for pred in preds if pred in pathSet]
                predecessors = defaultdict(list, [(NBP, preds) for NBP, preds in predecessors.items() if NBP in pathSet])
                for NBP, succs in successors.items():
                    successors[NBP]   = [succ for succ in succs if succ in pathSet]
                successors   = defaultdict(list, [(NBP, succs) for NBP, succs in successors.items()   if NBP in pathSet])

                    
                # Check that paths with one predecessor and successor happen only around branching points or in cycles
                if debug:
                    if isCycle:
                        assert len(NBPs) == 1
                    else:
                        for NBP in NBPs:
                            if len(predecessors[NBP])==1 and len(successors[NBP])==1:
                                pred = predecessors[NBP][0]
                                succ = successors[NBP][0]
                                assert len(successors[pred]) > 1 or len(predecessors[pred]) > 1
                                assert len(successors[succ]) > 1 or len(predecessors[succ]) > 1

                msg += f', processing {len(NBPs)} non-branching paths'
                print_time(msg, end = '\r')
                

                # By default we'll trim overlaps at 3', but we make exceptions
                # The idea is to remove also overlaps between two sequences sharing the same predecessor
                trim_left  = set()
                keep_right = set()
                for NBP in NBPs:
                    if len(successors[NBP]) > 1:
                        for succ in successors[NBP]:
                            if len(succ) > 1:
                                trim_left.add(succ)
                        keep_right.add(NBP)
                for NBP in NBPs:
                    if not successors[NBP] or len(NBP) < 2 or (len(NBP)<=2 and NBP in trim_left):
                        keep_right.add(NBP)
                    for succ in successors[NBP]:
                        if succ in trim_left:
                            keep_right.add(NBP)

                trimmedPaths = []
                for NBP in NBPs:
                    tp = tuple(NBP) # copy
                    if NBP not in keep_right:
                        tp = tp[:-1]
                    if NBP in trim_left:
                        tp = tp[1:]
                    assert tp
                    trimmedPaths.append(tp)

                if debug:
                    vs = {v for NBP in NBPs for v in NBP}
                    vst = {v for NBP in trimmedPaths for v in NBP}
                    assert vs == vst
                
                trimmedSeqs = []
                for NBP in NBPs:
                    tseq = NBP2seq[NBP]
                    if NBP not in keep_right:
                        tseq = tseq[:-self.ksize]
                    if NBP in trim_left:
                        tseq = tseq[self.ksize:]
                    trimmedSeqs.append(tseq)

                # Update global results
                    
                ori2covs = self.multimap(self.get_ori2cov, threads, trimmedPaths, seqPathsThisComp, bubble_vertex2origins)
                avgcovs = [sum(ori2cov.values()) for ori2cov in ori2covs]
                # Using trimmedPaths/trimmedSeqs: don't take into account overlaps with other paths for mean cov calculation and origin reporting

                # Sort based on trimmed seq length and cov
                idx = [i for i, tseq in sorted([(i, tseq) for i, tseq in enumerate(trimmedSeqs)], key = lambda x: (len(x[1]), round(avgcovs[x[0]], 2), x[1]), reverse = True)]
                NBPs = [NBPs[i] for i in idx]
                ori2covs = [ori2covs[i] for i in idx]
                avgcovs = [avgcovs[i] for i in idx]
                trimmedSeqs = [trimmedSeqs[i] for i in idx]
                del idx

                path2id  = {NBP: (scaffold, i) for i, NBP in enumerate(NBPs)}

                for NBP, ori2cov, avgcov, tseq in zip(NBPs, ori2covs, avgcovs, trimmedSeqs):
                    if debug:
                        assert NBP not in addedPaths
                    id_ = path2id[NBP]
                    # Even for complete genomes, we can have non-branching paths that belong to the core but don't have full coverage in some sequences
                    # This can happen if there are missing bases at the beginning or end of some of the input sequences (contigs in the original assemblies)
                    # So there is no branching, but some are missing
                    goodOris = {ori for ori, cov in ori2cov.items() if cov >= genome_assignment_threshold}
                    goodOris.add( list(sorted(ori2cov, key = lambda ori: ori2cov[ori], reverse = True))[0] ) # add at least the origin with the highest cov
                    goodOris = {ori.split('_Nsplit_')[0] for ori in goodOris}
                    # Go on
                    contigs[id_] = Contig(scaffold   = id_[0],
                                          i          = id_[1],
                                          seq        = NBP2seq[NBP],
                                          tseq       = tseq,
                                          cov        = avgcov,
                                          origins    = goodOris,
                                          successors = [path2id[succ] for succ in successors[NBP]])
                    
                    if debug:
                        addedPaths.add(NBP)

                msg += ', done'
                print_time(msg)

        ### Finished!
        return contigs


    ################################ AUX METHODS

    @staticmethod
    def set_vertex_filter(NBPG, nvs):
        """
        Set a vertex filter in a graph_tool Graph
        """
        for nv in NBPG.vertices():
            NBPG.vp.vfilt[nv] = nv in nvs
        NBPG.set_vertex_filter(NBPG.vp.vfilt)


    @classmethod
    def seq2hashes(cls, i, seqlength):
        """
        Obtain sequence hashes for the kmers of a forward and a reverse sequence
        Stores results in a shared memory buffer
        """
        # hashMem and revhashMem are two continuous shared mem buffers that will hold all the hashes concatenated
        # for each index, we just read/write from/to the proper slices, since hash size is constant
        seqMem, rcSeqMem, hashMem, ksize, clength, compressor = cls.get_multiprocessing_globals()
        seq, rcSeq = seqMem.buf[:seqlength], rcSeqMem.buf[:seqlength]
        ri = len(rcSeq) - i - ksize
        h  = compressor.compress(seq[i:i+ksize])
        rh = compressor.compress(rcSeq[ri:ri+ksize])
        h = h if h >= rh else rh
        hashMem.buf[i*clength:(i+1)*clength] = h


    @classmethod
    def is_NBP_init(cls, i):
        """
        Check if a given vertex of a DBG will be a non-branching path init
        This will happen if
        1) It has more than 2 outgoing edges (meaning that it is just an extender)
            (note that the edges in our DBG are bidirectional)
        2) The node has 2 outgoing edges but it is a fake extender
            - This happened for a example in which two input sequences diferred only in their
              last kmer
            - So original sequences were kA->kB, kC->kB (where kA,kB,kC are kmers)
            - But due to how we build our DBG, the DBG ends up being A=B=C (where = denotes a
               bidirectional link)
            - But a sequence can not really be reconstructed for the paths A->B->C or C->B->A
            - So we mark B as an init even if it looks like an extender in the DBG
        """
        vs, DBG, seqLimits, vertex2coords, ref2name, seqDict, ksize = cls.get_multiprocessing_globals()
        v = vs[i]
        succs = set(DBG.get_out_neighbors(v)) - {v} # - {v} to avoid self loops being init points
        isInit = False
        if len(succs) != 2: 
            isInit = True
        elif v in seqLimits:
            # We are reconstructing the sequence for the edges B->A and B->C
            # Where A,B,C are vertices in the DBG such that A=B=C ("=" denotes bidirectional link)
            #  so A,C are successors of B
            # Normally when reconstructing B->A and B->C, the valid kmer for B would be in forward
            #  orientation in one case, and in reverse complement orientation in the other,
            #  so len(kmers) == 2
            # However if the node is a fake extender, the valid kmer for B will be the same in
            #  B->A and B->C
            kmers = set()
            k1x = cls.vertex2kmer(v, vertex2coords, ref2name, seqDict, ksize)
            for s in succs:
                for k1 in (k1x, reverse_complement(k1x)):
                    ext = cls.extendKmer(k1, s, vertex2coords, ref2name, seqDict, ksize)
                    if ext:
                        kmers.add(k1)
            if len(kmers) == 1:
                isInit = True
        return isInit


    @classmethod
    def reconstruct_sequence(cls, i):
        """
        Reconstruct a nucleotide sequence for a De-Bruijn Graph path
        """
        NBPs, vertex2coords, ref2name, seqDict, ksize = cls.get_multiprocessing_globals()
        path = NBPs[i]
        kmers = []
        # Each vertex corresponds to two kmers (fwd and rev)
        # Which kmers give us a meaningful sequence?
        # We are reading this in a given direction (from path[0] to path[-1])
        # For the starting kmer (in position 0), pick the one that overlaps with any of the two kmers in position 1
        kmers = cls.edge2kmer(path[0], path[1],  vertex2coords, ref2name, seqDict, ksize)
        
        # Then just keep advancing in the path, for each vertex pick the kmer that overlaps with the previously-picked kmer
        for v in path[2:]:
            kmers.append(cls.extendKmer(kmers[-1], v, vertex2coords, ref2name, seqDict, ksize, force = True))
            
        assert len(kmers) == len(path) # We can remove this assertion, since edge2kmer and extendKmer(force=True)
                                       #  have internal assertions
        
        return cls.seq_from_kmers(kmers)


    @classmethod
    def vertex2kmer(cls, v, vertex2coords, ref2name, seqDict, ksize):
        """
        Return the kmer corresponding to a given vertex in the De-Bruijn Graph
        """
        ref, idx = vertex2coords[v]
        name = ref2name[ref]
        return seqDict[name][ idx : (idx + ksize) ]


    @classmethod
    def edge2kmer(cls, v1, v2,  vertex2coords, ref2name, seqDict, ksize):
        """
        Return the kmers corresponding to a given edge in the De-Bruijn Graph
        """
        # Each vertex corresponds to two kmers (fwd and rev)
        # Which kmers give us a meaningful sequence?
        # We are reading this in a given direction (from path[0] to path[-1])
        # For the starting kmer (in position 0), pick the one that overlaps with any of the two kmers in position 1
        kmers = []
        k1x = cls.vertex2kmer(v1, vertex2coords, ref2name, seqDict, ksize)
        for k1 in (k1x, reverse_complement(k1x)):
            k2 = cls.extendKmer(k1, v2,  vertex2coords, ref2name, seqDict, ksize)
            if k2:
                kmers.extend([k1, k2])
        
        assert len(kmers) == 2
        return kmers


    @classmethod
    def extendKmer(cls, k1, v2, vertex2coords, ref2name, seqDict, ksize, force = False):
        """
        Return the kmer that extends a given first kmer k1, where k1 is a 
         kmer coming from a vertex v1, and v1,v2 is an edge in the De-Bruijn Graph
        """
        kmers = []
        k2x = cls.vertex2kmer(v2, vertex2coords, ref2name, seqDict, ksize)
        for k2 in (k2x, reverse_complement(k2x)):
            if k2[:-1] == k1[1:]:
                kmers.append(k2)
        if force:
            assert len(kmers) == 1
        else:
            assert len(kmers) <= 1 # can be 0 if k1 is in the wrong orientation
        return kmers[0] if kmers else None


    @staticmethod
    def seq_from_kmers(kmers):
        """
        Reconstruct a nucleotide sequence from a list of kmers
        """
        fullseq = [kmers[0]]
        extra = [kmer[-1] for kmer in kmers[1:]]
        fullseq.extend(extra)
        return ''.join(fullseq)

    
    @classmethod
##    @profile
    def get_seqPaths_NBPs(cls, i):
        """
        Return the vertices in the Sequence Graph that are contained in a given input sequence in its original orientation
        """
        names, seqPaths, psets, NBP2seq, vertex2NBP, seqDict = cls.get_multiprocessing_globals()
        name = names[i]
        seqPath = seqPaths[name]     
        
        spGS = set()
        for pset, (nv1, nv2) in psets.items():
            NBP = vertex2NBP[nv1]
            cNBP = compress_vertices(np.array(NBP, dtype=np.uint32))
            refLen = len(NBP) if NBP[0] != NBP[-1] else len(NBP) - 1 # A NBP will only have duplicated vertices if it is circular
            if vertex_overlap(cNBP, seqPaths[name]) == refLen: # the pset is a subset of the seqPath.
                for nv in (nv1, nv2):
                    if NBP2seq[vertex2NBP[nv]] in seqDict[name]: # check the string to make sure that we only add nodes from the same orientation
                        spGS.add(nv)
        return spGS



    @classmethod
    def get_avgPathCov(cls, i):
        """
        For a given path in the De-Bruijn Graph, calculate the average coverage of its vertices in the input sequences
        """
        NBPs, seqPaths = cls.get_multiprocessing_globals()
        path = NBPs[i]
        return np.mean([len(cls.vertex2origins(v, seqPaths)) for v in path])



    @classmethod
    def get_ori2cov(cls, i):
        """
        For a given path in the De-Bruijn Graph, calculate the fraction of its vertices that are contained in each input sequence
        """
        NBPs, seqPaths, bubble_vertex2origins = cls.get_multiprocessing_globals()
        path = NBPs[i]
        
        ori2cov = defaultdict(int)
        for v in path:
            for ori in cls.vertex2origins(v, seqPaths) | set(bubble_vertex2origins[v]): # note how I'm also getting origins not really associated to this path
                ori2cov[ori] += 1                                                       # but to bubbles that were originally in this path before we popped them
        return {ori: cov/len(path) for ori, cov in ori2cov.items()}


    @classmethod
    def get_vertex2origins(cls, i):
        vertices, seqPaths = cls.get_multiprocessing_globals()
        return cls.vertex2origins(vertices[i], seqPaths)
        

    @classmethod
    def vertex2origins(cls, v, seqPaths):
        """
        Return the names of the input sequences that contain a given vertex in the De-Bruijn Graph
        """
        return {name for name in seqPaths if isInSeqPath(v, seqPaths[name])}



    @classmethod
    def get_vertex2originslen(cls, i):
        """
        Return the number of input sequences that contain a given vertex in the De-Bruijn Graph
        """
        vertices, seqPaths = cls.get_multiprocessing_globals()
        return len(cls.vertex2origins(vertices[i], seqPaths))





        

