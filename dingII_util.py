from ilp_util_adj import *
import networkx as nx
from argparse import ArgumentParser, FileType

HEAD = 'h'
TAIL = 't'
TELO = 'o'



def read_genomes(args):
    allgnms = readGenomes(args.unimog.readlines())
    if args.pair:
        gnms = [g for g in allgnms if g[0] in args.pair]
    else:
        indices = [0,1]
        if args.pairnumber:
            indices = args.pairnumber
        gnms = []
        for i in indices:
            if i < 0 or i >= len(allgnms):
                raise Exception('Pair index %i greater than number of genomes provided.'%i)
            gnms.append(allgnms[i])
    if len(gnms) > 2:
        raise Exception('Multiple genomes of that name found!')
    elif len(gnms) < 2:
        raise Exception('Some genome names could not be found in the provided file!')
    return gnms

TELOGENE = ''

def raw_extremities(chrm, idgen):
    exts = []
    if chrm[0] == CHR_LINEAR:
        t1id = idgen.get_new()
        exts.append((t1id, TELOGENE , TELO))
    for o,g in chrm[1]:
        xs = [(idgen.get_new(),g, e) for e in [TAIL, HEAD]]
        if o == ORIENT_NEGATIVE:
            xs = xs[::-1]
        exts.extend(xs)
    if chrm[0] == CHR_LINEAR:
        t2id = idgen.get_new()
        exts.append((t2id, TELOGENE, TELO))
    return exts
        

SELFEDGE = 'self'
ADJACENCY = 'adj'
EXTREMITYEDGE = 'ext'

GENOME1 = '1'
GENOME2 = '2'

def genome_graph(gnm, idgen, name='A', pad = 0, return_circulars=True):
    gg = nx.MultiGraph()
    #for each gene store all adjacencies of all occurrences
    gene_index = {}
    gene_index[HEAD] = {}
    gene_index[TAIL] = {}
    gene_index[TELO] = {}
    circs = []
    extremity_genome = []
    for chrm in gnm[1]:
        exts = raw_extremities(chrm, idgen)
        extremity_genome.append(exts)
        #store self edges for circular chromosomes
        selfs = []
        for eid, g, xt in exts:
            gg.add_node(eid, genome=name, gene=g, extremity=xt)
            insertl(gene_index[xt], g, eid)
        #whether the first edge is a self edge
        if chrm[0] == CHR_CIRCULAR:
            laste = exts[-1][0]
        else:
            laste = exts[0][0]
            exts = exts[1:]
        selfedge = False
        for thise, _,_ in exts:
            if selfedge:
                gg.add_edge(laste, thise, etype=SELFEDGE)
                thisselfedgeid = [k for k, data in gg[laste][thise].items() if data['etype']==SELFEDGE]
                if len(thisselfedgeid) != 1:
                    raise Exception('Wrong number of selfedges (%i) between %i and %i.'%(len(thisselfedgeid),laste, thise))
                selfs.append((laste, thise, thisselfedgeid[0]))
            else:
                gg.add_edge(laste, thise, etype=ADJACENCY)
            selfedge = not selfedge
            laste = thise
        if chrm[0] == CHR_CIRCULAR:
            circs.append(selfs)
    for i in range(pad):
        t1 = idgen.get_new()
        t2 = idgen.get_new()
        gg.add_node(t1, genome=name, gene=TELOGENE, extremity=TELO)
        gg.add_node(t2, genome=name, gene=TELOGENE, extremity=TELO)
        gg.add_edge(t1,t2, etype=ADJACENCY)
        insertl(gene_index[TELO], TELOGENE, t1)
        insertl(gene_index[TELO], TELOGENE, t2)
    return gg, gene_index, circs, extremity_genome
    
def relational_diagram(gg1, gg2, gi1, gi2):
    rd = nx.union(gg1, gg2)
    for x in [HEAD, TAIL, TELO]:
        for g, occs1 in gi1[x].items():
            occs2 = getl(gi2[x], g)
            for id1 in occs1:
                for id2 in occs2:
                    rd.add_edge(id1,id2, etype=EXTREMITYEDGE)
    return rd
    
def full_relational_diagram(genomes, LOG):
    gen = Simple_Id_Generator()
    lchrs1 = len([c[0] for c in genomes[0][1] if c[0] == CHR_LINEAR])
    lchrs2 = len([c[0] for c in genomes[1][1] if c[0] == CHR_LINEAR])
    diff = lchrs1 - lchrs2
    pad1 = 0
    pad2 = 0
    if diff > 0:
        pad2 = diff
        LOG.info('Padding genome %s with %i empty chromosomes.'%(genomes[1][0], pad2))
    if diff < 0:
        pad1 = -diff
        LOG.info('Padding genome %s with %i empty chromosomes.'%(genomes[0][0], pad1))
    gg1, gi1, circs, eg1 = genome_graph(genomes[0], gen, name=GENOME1, pad=pad1)
    gg2, gi2, circs2, eg2 = genome_graph(genomes[1], gen, name=GENOME2, pad=pad2)
    circs.extend(circs2)
    rd = relational_diagram(gg1, gg2, gi1, gi2)
    return rd, gi1, gi2, circs, [eg1, eg2]

def add_unimog_parsing_groups(parser):
    parser.add_argument('unimog', type=FileType('r'), help='The genomes provided in UniMOG format.')
    pairs = parser.add_mutually_exclusive_group()
    pairs.add_argument('-p', '--pair', type=str, nargs=2, help='Give the two names of the genomes you want to compare, as specified in the unimog header.')
    pairs.add_argument('-pn', '--pairnumber', type=int, nargs=2, help='Chose the two genomes via their position in the file (starting at 0). Default: 0,1')

def read_gurobi(fil):
    ''' Read a gurobi solution file and return maxfvalue, {variableclass : {vid:value}}
    '''
    obj = 'N/A'
    vrs = {}
    for line in fil:
        if line.startswith('#'):
            s = (line.split('# Objective value ='))
            if len(s) > 1:
                obj = int(float(s[1].strip()))
            continue
        s = line.split(' ')
        vname = s[0][0]
        vid = s[0][2:]
        if vname not in vrs:
            vrs[vname] = {}
        vrs[vname][vid] = int(round(float(s[1].strip())))
    return obj, vrs
