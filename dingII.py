#!/usr/bin/python3
from argparse import ArgumentParser, FileType, ArgumentTypeError
import logging
from ilp_util_adj import *
import networkx as nx
import math
import sys
from dingII_util import *

def matching_range(v):
    iv = float(v)
    if iv < 0.0 or iv > 1.0:
        raise ArgumentTypeError("%s is an invalid range bound. Please only provide values in [0,1]." % v)
    return iv

def mm_bound(a,b):
    m = min(a,b)
    return (m,m)
    
def em_bound(a,b):
    if a == 0 or b == 0:
        return (0,0)
    return (1,1)

def im_bound(a,b):
    if a == 0 or b == 0:
        return (0,0)
    m = min(a,b)
    return (1,m)

class RangedBound:
    def r_bound(self, a, b):
        m = min(a,b)
        return (int(math.ceil(self.lower*m)), int(math.ceil(self.upper*m)))



def create_model(amults, bmults, boundf):
    bounds = {}
    for gene, mult in amults.items():
        bounds[gene] = boundf(mult, dgetc(bmults, gene))
    for gene, mult in bmults.items():
        bounds[gene] = boundf(mult, dgetc(amults, gene))
    return bounds
        

def create_matching_model(args, gnms):
    amults = get_multiplicities(gnms[0])
    bmults = get_multiplicities(gnms[1])
    needs_supplement = True
    model = {}
    if args.custom:
        model = read_matching_model(args.custom)
        enforce_model_integrity(model, amults, bmults)
        needs_supplement = False
        if is_incomplete(model, amults, bmults):
            LOG.warning('The matching model provided via custom file is incomplete. It will be supplemented by your chosen model.')
            needs_supplement = True
    if needs_supplement:    
        if args.exemplary:
            boundf = em_bound
        elif args.intermediate:
            boundf = im_bound
        elif args.range:
            rb = RangedBound()
            bs = sorted(args.range)
            rb.lower = bs[0]
            rb.upper = bs[1]
            boundf = rb.r_bound
        else:
            if not args.maximal:
                LOG.info('No matching model chosen. Default is maximal matching.')
            boundf = mm_bound
        model2 = create_model(amults, bmults, boundf)
        #supplement model so far
        for gene in model2:
            if not gene in model:
                model[gene] = model2[gene]
    return model
    
def canon_disp_tup(tp):
    tp_ = tp[0:2]
    return '%i_%i_%i'%(min(tp_),max(tp_),tp[2])

def get_multiplicities(gnm):
    mlt = {}
    for _, chrm in gnm[1]:
        for _, gene in chrm:
            dadd(mlt,gene,1)
    return mlt


def write_matching_model(model, fl):
    fl.write('Gene\tLower\tUpper\n')
    for gene, bnds in model.items():
        l, u = bnds
        fl.write('%s\t%i\t%i\n'%(gene,l,u))

def read_matching_model(fl):
    lines = fl.readlines()
    #throw away header
    lines = lines[1:]
    model = {}
    for line in lines:
        entries = line.strip().split('\t')
        model[entries[0]] = (int(entries[1]),int(entries[2]))
    return model

def enforce_model_integrity(model, amults, bmults):
    for gene, bnds in model.items():
        l, u = bnds
        m = min(dgetc(amults, gene),dgetc(bmults, gene))
        if l < 0:
            LOG.warning('Inconsistend bound %i will be reset to 0.'%l)
            l= 0
        if u < 0:
            LOG.warning('Inconsistend bound %i will be reset to 0.'%l)
            l= 0
        if l > u:
            LOG.warning('Inconsistent bounds: upper below lower (lower %i, upper %i) for gene %s. Both bounds will be set to %i.'%(l,u,gene, u))
            l = u
        if l > m:
            LOG.warning('Too high bound %i for gene %s (smallest family has only %i occurrences) will be trimmed to %i.'%(l,gene,m,m))
            l = m
        if u > m:
            LOG.warning('Too high bound %i for gene %s (smallest family has only %i occurrences) will be trimmed to %i.'%(u,gene,m,m))
            u = m
        if l == 0 and m > 0:
            LOG.info('Lower bound found for gene %s is 0. If this affects many genes and is not counteracted by deletion penalties this may lead to very large indel stretches.'%gene)
        model[gene] = (l,u)

def is_incomplete(model, amults, bmults):
    for gene in amults:
        if not gene in model:
            return True
    for gene in bmults:
        if not gene in model:
            return True
    return False

def get_first_extremity(orient):
    if orient == ORIENT_POSITIVE:
        return TAIL
    else:
        return HEAD

def get_last_extremity(orient):
    if orient == ORIENT_POSITIVE:
        return HEAD
    else:
        return TAIL



def siblings(rd, gi1, gi2):
    sibs = {}
    for g, occs1 in gi1[HEAD].items():
        occs2 = getl(gi2[HEAD],g)
        for hid1 in occs1:
            for hid2 in occs2:
                edgeid1 = list(rd[hid1][hid2].keys())
                if len(edgeid1) != 1:
                    raise Exception('Wrong number of extremity edges (%i) between vertices %i %i.'%(len(edgeid1),hid1,hid2))
                sister = (hid1, hid2,edgeid1[0])
                tid1 = maybe_adjacent_selfedge_raw(rd,hid1)
                if len(tid1) != 1:
                    raise Exception('Vertex %i connected to %i self edges!'%(hid1, len(tid1)))
                else:
                    tid1 = tid1[0][1]
                tid2 = maybe_adjacent_selfedge_raw(rd,hid2)
                if len(tid2) != 1:
                    raise Exception('Vertex %i connected to %i self edges!'%(hid2, len(tid2)))
                else:
                    tid2 = tid2[0][1]
                edgeid2 = list(rd[tid1][tid2].keys())
                if len(edgeid2) != 1:
                    raise Exception('Wrong number of extremity edges (%i) between vertices %i %i.'%(len(edgeid2),tid1,tid2))
                brother = (tid1, tid2,edgeid2[0])
                insertl(sibs, g, (brother, sister))
    return sibs
                

def add_constraint(ilp, constraint_type, constraint_string):
    cid = len(getl(ilp, constraint_type))
    insertl(ilp, constraint_type, (cid, constraint_string))
    

def maybe_adjacent_selfedge(rd, u):
    return [canon_disp_tup(t) for t in maybe_adjacent_selfedge_raw(rd, u)]
    
def maybe_adjacent_selfedge_raw(rd,u):
    return [(u,v,k) for v, edges in rd[u].items() for k, data in edges.items() if data['etype']==SELFEDGE]

def edge_constraints(rd, ilp,max_indels=False):
    '''
    Add constraints C01, 04, 05, 07, 08, 09, 10 by iterating through edges once.
    '''
    for u,v, k, data in rd.edges(data=True, keys=True):
        e = canon_disp_tup((u,v,k))
        etype= data['etype']
        insertl(ilp, 'binary', 'x_%s'%e)
        insertl(ilp, 'binary', 't_%s'%e)
        if etype == ADJACENCY:
            #c01
            cs = 'x_%s = 1'%e
            add_constraint(ilp, 'c01', cs)
        if etype == SELFEDGE:
            #c05
            for i in [u,v]:
                if not max_indels:
                    cs = 'y_%i + %i x_%s <= %i'%(i,i,e,i)
                    add_constraint(ilp, 'c05', cs)
                else:
                    cs = 'x_%s - d_%i <= 0'%(e,i)
            #c07
            if rd.nodes[u]['genome'] == GENOME1:
                #if (u,v) is selfedge genome is the same for u,v
                for v_ in [u,v]:
                    cs = 'r_%i + x_%s <= 1'%(v_,e)
                    add_constraint(ilp, 'c07', cs)
            else:
                #case for GENOME2
                for v_ in [u,v]:
                    cs = 'r_%i - x_%s >= 0'%(v_,e)
                    add_constraint(ilp, 'c07', cs)
        if etype != SELFEDGE or max_indels:
            #c04 DONE: this is only relevant for non-self edges
            for i,j in [(u,v), (v,u)]:
                cs = 'y_%i - y_%i + %i x_%s <= %i'%(i,j,i,e, i)
                add_constraint(ilp, 'c04', cs)
        #c08
        for u_,v_ in [(u,v),(v,u)]:
            cs = 't_%s - r_%i + r_%i - x_%s >= -1'%(e, v_,u_,e)
            add_constraint(ilp, 'c08', cs)
        if etype == ADJACENCY and rd.nodes[u]['genome'] == GENOME1:
            #c09
            d1 = maybe_adjacent_selfedge(rd, u)
            d2 = maybe_adjacent_selfedge(rd, v)
            d1.extend(d2)
            d1 = list(set(d1))
            if len(d1) > 0:
                cs = ' + '.join(['x_%s'%d for d in d1 ]) + ' - t_%s >= 0'%(e)
                add_constraint(ilp, 'c09', cs)
        else:
            #c10
            cs = 't_%s = 0'%e
            add_constraint(ilp, 'c10', cs)
        insertl(ilp, 'obj', (0.5, 't_%s'%e))
        if etype == EXTREMITYEDGE:
            insertl(ilp, 'obj', (0.5, 'x_%s'%e))

def sibs_constraints(sibs, ilp, matching_model):
    '''
    Implement Constraints c03, c12.
    '''
    for fam in sibs:
        for d_, e_ in sibs[fam]:
            d = canon_disp_tup(d_)
            e = canon_disp_tup(e_)
            cs = 'x_%s - x_%s = 0'%(d,e)
            add_constraint(ilp, 'c03', cs)
        headsum = ' + '.join(['x_%s'%canon_disp_tup(d) for d,_ in sibs[fam]])
        lower, upper = matching_model[fam]
        ch1 = '%s >= %i'%(headsum, lower)
        ch2 = '%s <= %i'%(headsum, upper)
        add_constraint(ilp,'c12',ch1)
        add_constraint(ilp,'c12',ch2)
        
                
def vertex_constraints(rd, ilp):
    '''
    Implement constraints c02, c06.
    '''
    for v in rd.nodes():
        insertl(ilp, 'binary', 'z_%i'%v)
        insertl(ilp, 'binary', 'r_%i'%v)
        insertl(ilp, 'general', ('y_%i'%v, 0, v))
        #c02
        adjacent_edges = [canon_disp_tup((v,u,k)) for u, edges in rd[v].items() for k in edges]
        smd = ' + '.join(['x_%s'%e for e in adjacent_edges])
        add_constraint(ilp, 'c02', '%s = 2'%smd)
        #c06
        cs = '%i z_%i - y_%i <= 0'%(v,v,v)
        add_constraint(ilp, 'c06',cs)
        insertl(ilp, 'obj', (-1, 'z_%i'%v))

def singleton_constraints(circs, ilp):
    for i, c in enumerate(circs):
        k = len(c)
        smd = ' + '.join(['x_%s'%canon_disp_tup(e) for e in c])
        add_constraint(ilp, 'c11', '%s - s_%i <= %i'%(smd, i, k-1))
        insertl(ilp, 'obj', (1, 's_%i'%i))
        insertl(ilp, 'binary', 's_%i'%i)

def maximize_indels_constraints(rd,ilp):
    for u,v, k, data in rd.edges(data=True, keys=True):
        e = canon_disp_tup((u,v,k))
        etype= data['etype']
        for i,j in [(u,v),(v,u)]:
            cns = 'd_%i - d_%i + x_%s <= 1'%(i,j,e)
            add_constraint(ilp,'c13id',cns)
        if etype==SELFEDGE:
            for i in [u,v]:
                cns = 'x_%s - d_%i <= 0'%(e,i)
                add_constraint(ilp,'c05id',cns)
    for v in rd.nodes():
        insertl(ilp,'binary','d_%i'%v)
        insertl(ilp,'binary','p_%i'%v)
        add_constraint(ilp,'c15id','z_%i + d_%i - p_%i <= 1'%(v,v,v))
        insertl(ilp,'obj',(1,'p_%i'%v))
    num_nodes = len(list(rd.nodes()))
    insertl(ilp,'general',('max_indels',0,num_nodes))
    tsum = ' - '.join(['0.5 t_%s'%canon_disp_tup((u,v,k)) for u,v,k,data in rd.edges(data=True,keys=True)])
    psum = ' - '.join(['p_%i'%v for v in rd.nodes()])
    add_constraint(ilp, 'c17id','max_indels - %s - %s = 0'%(tsum,psum))
    ilp['obj1'] = [(-1,'max_indels')]

def maximize_dcj_constraints(rd,ilp):
    num_nodes = len(list(rd.nodes()))
    for v in rd.nodes():
        insertl(ilp,'binary','w_%i'%v)
        insertl(ilp,'binary','d_%i'%v)
        #insertl(ilp,'general', ('a_%i'%v,0,num_nodes+1))
        #insertl(ilp,'general', ('b_%i'%v,0,num_nodes+1))
        cns = 'a_%i - b_%i - %i w_%i - %i d_%i - %i y_%i <= -1'%(v,v,num_nodes,v, num_nodes,v,num_nodes,v)
        add_constraint(ilp,'c13dcj',cns)
        cns = 'b_%i - a_%i - %i w_%i + %i d_%i - %i y_%i <= %i'%(v,v,num_nodes,v, num_nodes,v,num_nodes,v,num_nodes -1)
        add_constraint(ilp,'c13dcj',cns)
    for u,v, k, data in rd.edges(data=True, keys=True):
        e = canon_disp_tup((u,v,k))
        etype= data['etype']
        if etype==ADJACENCY:
            if rd.nodes[u]['genome'] == GENOME1:
                for i,j in [(u,v),(v,u)]:
                    cns = 'a_%i - a_%i - %i t_%s <= 0'%(i,j,num_nodes,e)
                    add_constraint(ilp,'c14dcj',cns)
            else:
                cns = 'a_%i - a_%i = 0'%(u,v)
                add_constraint(ilp,'c14dcj',cns)
        else:
            for i,j in [(u,v),(v,u)]:
                cns = 'b_%i - b_%i + %i x_%s <= %i'%(i,j,num_nodes,e,num_nodes)
                add_constraint(ilp,'c15dcj',cns)
        cns = 'd_%i + d_%i + x_%s <= 2'%(u,v,e)
        add_constraint(ilp,'c16dcj',cns)
        cns = 'd_%i + d_%i - x_%s >= 0'%(u,v,e)
        add_constraint(ilp,'c16dcj',cns)
    insertl(ilp,'binary','n')
    insertl(ilp,'general',('min_indels',0,num_nodes))
    tsum = ' - '.join(['t_%s'%canon_disp_tup((u,v,k)) for u,v,k,data in rd.edges(data=True,keys=True)])
    cns = '%i n - %s >= 0'%(num_nodes,tsum)
    add_constraint(ilp,'c17dcj',cns)
    wsum = ' - '.join(['w_%i'%u for u in rd.nodes()])
    cns = 'min_indels - %s - 2 n = 0'%wsum
    add_constraint(ilp,'c18dcj',cns)
    ilp['obj1']=[(1,'min_indels')]


def print_constraints(ilp, ctype, file=sys.stdout):
    for num, cs in getl(ilp, ctype):
        print(' %s.%i: %s'%(ctype,num, cs),file=file)

def print_all_constraints(ilp, file=sys.stdout,max_indels=False,max_rearrangements=False):
    for ctype in ['c%02d'%i for i in range(1,13)]:
        print_constraints(ilp, ctype, file=file)
    if max_indels:
        for ctype in ['c05id','c13id','c15id','c17id']:
            print_constraints(ilp,ctype,file=file)
    if max_rearrangements:
        for ctype in ['c13dcj','c14dcj','c15dcj','c16dcj','c17dcj','c18dcj']:
            print_constraints(ilp,ctype,file=file)

def print_objective(ilp, file=sys.stdout):
    print('Minimize ',file=file)
    obj_strpos = ' + '.join(['%f %s'%(k, var) for (k, var) in ilp['obj'] if k > 0])
    obj_strneg =' - '.join(['%f %s'%(-k, var) for (k, var) in ilp['obj'] if k < 0])
    if obj_strneg == '':
        obj_str=obj_strpos
    else:
        obj_str = ' - '.join([obj_strpos,obj_strneg])
    print(' obj: %s'%obj_str, file=file)

def print_with_secondary_objective(ilp,file=sys.stdout):
    print('Minimize multi-objectives',file=file)
    obj_strpos = ' + '.join(['%f %s'%(k, var) for (k, var) in ilp['obj'] if k > 0])
    obj_strneg=' - '.join(['%f %s'%(-k, var) for (k, var) in ilp['obj'] if k < 0])
    if obj_strneg == '':
        obj_str=obj_strpos
    else:
        obj_str = ' - '.join([obj_strpos,obj_strneg])
    obj_strbpos = ' + '.join(['%f %s'%(k, var) for (k, var) in ilp['obj1'] if k > 0])
    obj_strbneg =' - '.join(['%f %s'%(-k, var) for (k, var) in ilp['obj1'] if k < 0])
    if obj_strbneg == '':
        obj_strb=obj_strbpos
    else:
        obj_strb = ' - '.join([obj_strbpos,obj_strbneg])
    print(' obj0: Priority=2 Weight=1 AbsTol=0 RelTol=0',file=file)
    print('  '+obj_str,file=file)
    print(' obj1: Priority=1 Weight=1 AbsTol=0 RelTol=0',file=file)
    print('  '+obj_strb,file=file)

def print_domains(ilp, file=sys.stdout):
    for var, lower, upper in ilp['general']:
        print(' %i <= %s <= %i'%(lower, var, upper), file=file)

def print_generals(ilp, file=sys.stdout):
    for var, _ , _ in ilp['general']:
        print(' %s'%var, file=file)

def print_binaries(ilp, file=sys.stdout):
    for var in ilp['binary']:
        print(' %s'%var, file=file)

def print_ilp(ilp, file=sys.stdout,max_indels=False,max_rearrangements=False):
    if not (max_indels or max_rearrangements):
        print_objective(ilp, file=file)
    else:
        print_with_secondary_objective(ilp,file=file)
    print('Subject To', file=file)
    print_all_constraints(ilp, file=file,max_indels=max_indels,max_rearrangements=max_rearrangements)
    print('Bounds', file=file)
    print_domains(ilp, file=file)
    print('General', file=file)
    print_generals(ilp, file=file)
    print('Binary', file=file)
    print_binaries(ilp, file=file)
    print('End', file=file)

LOGLEVELS = {'debug':logging.DEBUG,'info':logging.INFO,'warning':logging.WARNING,'error':logging.ERROR,'critical':logging.CRITICAL}    

def main():
    parser = ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-mm','--maximal', action='store_true', help='Set matching model to maximal matching.')
    group.add_argument('-em','--exemplary', action='store_true', help='Set matching model to exemplary matching.')
    group.add_argument('-im','--intermediate', action='store_true', help='Set matching model to intermediate matching.')
    group.add_argument('-r', '--range', type=matching_range, nargs=2, help='Provide upper and lower percentiles to be matched per marker in range [0,1]. Actual discrete bounds will always be rounded up.')
    parser.add_argument('-c','--custom', type=FileType('r'), action='store', help='Provide a custom matching file.')
    add_unimog_parsing_groups(parser)
    writewhat = parser.add_mutually_exclusive_group(required=True)
    writewhat.add_argument('--writemodel', type=FileType('w'), help='Write the matching model to a file in order to customize it.')
    writewhat.add_argument('--writeilp', type=FileType('w'), help='Write the resulting ILP to the specified file.')
    writewhat.add_argument('--toheuristicsol',type=FileType('w'), help='Convert a user specified matching to a heuristic starting solution.')
    parser.add_argument('--heuristicmatching',type=FileType('r'), help='User specified matching to convert to heuristic solution. Only to be used with --heuristicsol.')
    maximizewhat = parser.add_mutually_exclusive_group()
    maximizewhat.add_argument('--maximizeindels',action='store_true',help='**Warning: Highly experimental, do not use if you are not a ding dev - gurobi only** Under all co-optimal solutions prefer those that maximize the number of indel operations and minimize the number of DCJs.')
    maximizewhat.add_argument('--maximizerearrangements',action='store_true',help='**Warning: Highly experimental, do not use if you are not a ding0 dev - gurobi only** Under all co-optimal solutions prefer those that minimize the number of indel operations and maximize the number of DCJs.')
    parser.add_argument('--log-level',choices=LOGLEVELS.keys(),default='warning')
    args = parser.parse_args()
    st = logging.StreamHandler(sys.stderr)
    LOG.setLevel(LOGLEVELS[args.log_level])
    st.setLevel(LOGLEVELS[args.log_level])
    st.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(st)
    if args.range:
        if min(args.range) == 0.0:
            LOG.warning('Lowest percentile to be matched is 0. If this is not counteracted by deletion penalties, almost only whole chromosome deletions will be modeled!')
    genomes = read_genomes(args)
    model = create_matching_model(args, genomes)
    if args.writemodel:
        write_matching_model(model, args.writemodel)
        sys.exit(0)
    rd, gi1, gi2, circs, exts = full_relational_diagram(genomes, LOG)
    LOG.debug('Graph:')
    for c in rd.nodes(data=True):
        LOG.debug(c)
    for c in rd.edges(data=True):
        LOG.debug(c)
    sibs = siblings(rd, gi1, gi2)
    ilp = {}
    edge_constraints(rd, ilp,max_indels=args.maximizeindels)
    sibs_constraints(sibs, ilp, model)
    vertex_constraints(rd, ilp)
    singleton_constraints(circs, ilp)
    if args.maximizeindels:
        maximize_indels_constraints(rd,ilp)
    if args.maximizerearrangements:
        maximize_dcj_constraints(rd,ilp)
    print_ilp(ilp, file=args.writeilp,max_indels=args.maximizeindels,max_rearrangements=args.maximizerearrangements)
    
LOG = logging.getLogger(__name__)





main()


