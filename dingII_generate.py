import math
import sys
from dingII.dingII_util import *
import logging

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
            logging.warning('The matching model provided via custom file is incomplete. It will be supplemented by your chosen model.')
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
                logging.info('No matching model chosen. Default is maximal matching.')
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
            logging.warning('Inconsistend bound %i will be reset to 0.'%l)
            l= 0
        if u < 0:
            logging.warning('Inconsistend bound %i will be reset to 0.'%l)
            l= 0
        if l > u:
            logging.warning('Inconsistent bounds: upper below lower (lower %i, upper %i) for gene %s. Both bounds will be set to %i.'%(l,u,gene, u))
            l = u
        if l > m:
            logging.warning('Too high bound %i for gene %s (smallest family has only %i occurrences) will be trimmed to %i.'%(l,gene,m,m))
            l = m
        if u > m:
            logging.warning('Too high bound %i for gene %s (smallest family has only %i occurrences) will be trimmed to %i.'%(u,gene,m,m))
            u = m
        if l == 0 and m > 0:
            logging.info('Lower bound found for gene %s is 0. If this affects many genes and is not counteracted by deletion penalties this may lead to very large indel stretches.'%gene)
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

def edge_constraints(rd, ilp):
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
                cs = 'y_%i + %i x_%s <= %i'%(i,i,e,i)
                add_constraint(ilp, 'c05', cs)
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
        else:
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


def print_constraints(ilp, ctype, file=sys.stdout):
    for num, cs in getl(ilp, ctype):
        print(' %s.%i: %s'%(ctype,num, cs),file=file)

def print_all_constraints(ilp, file=sys.stdout):
    for ctype in ['c%02d'%i for i in range(1,13)]:
        print_constraints(ilp, ctype, file=file)

def print_objective(ilp, file=sys.stdout):
    obj_str = ' + '.join(['%f %s'%(k, var) for (k, var) in ilp['obj'] if k > 0])
    obj_str+=' - '
    obj_str+=' - '.join(['%f %s'%(-k, var) for (k, var) in ilp['obj'] if k < 0])
    print(' obj: %s'%obj_str, file=file)

def print_domains(ilp, file=sys.stdout):
    for var, lower, upper in ilp['general']:
        print(' %i <= %s <= %i'%(lower, var, upper), file=file)

def print_generals(ilp, file=sys.stdout):
    for var, _ , _ in ilp['general']:
        print(' %s'%var, file=file)

def print_binaries(ilp, file=sys.stdout):
    for var in ilp['binary']:
        print(' %s'%var, file=file)

def print_ilp(ilp, file=sys.stdout):
    print('Minimize', file=file)
    print_objective(ilp, file=file)
    print('Subject To', file=file)
    print_all_constraints(ilp, file=file)
    print('Bounds', file=file)
    print_domains(ilp, file=file)
    print('General', file=file)
    print_generals(ilp, file=file)
    print('Binary', file=file)
    print_binaries(ilp, file=file)
    print('End', file=file)


def main(args):
    if args.range:
        if min(args.range) == 0.0:
            logging.warning('Lowest percentile to be matched is 0. If this is not counteracted by deletion penalties, almost only whole chromosome deletions will be modeled!')
    genomes = read_genomes(args)
    model = create_matching_model(args, genomes)
    if args.writemodel:
        write_matching_model(model, args.writemodel)
        sys.exit(0)
    rd, gi1, gi2, circs, exts = full_relational_diagram(genomes, logging)
    logging.debug('Graph:')
    for c in rd.nodes(data=True):
        logging.debug(c)
    for c in rd.edges(data=True):
        logging.debug(c)
    sibs = siblings(rd, gi1, gi2)
    ilp = {}
    edge_constraints(rd, ilp)
    sibs_constraints(sibs, ilp, model)
    vertex_constraints(rd, ilp)
    singleton_constraints(circs, ilp)
    print_ilp(ilp, file=args.writeilp)

