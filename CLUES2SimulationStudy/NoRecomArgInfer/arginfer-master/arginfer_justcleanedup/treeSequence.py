"""
Tree sequence conversion to ARGnode
"""
import random
import numpy as np
import bintrees
import collections
import msprime
import math
import argbook

class TreeSeq(object):
    """
    convert tree sequences (ts_full) to Node data structure, and vice versa
    """
    def __init__(self, ts_full):
        self.ts_full = ts_full
        self.tables= ts_full.tables
        self.edges_dict= collections.defaultdict(list)
        for edge in self.tables.edges:
            self.edges_dict[edge.parent].append(edge)
        self.n = ts_full.sample_size
        self.seq_length = ts_full.sequence_length
        self.mutation_map = self.map_mutation()#[[SNP]] where index is ts_node
        # print("selfmutations", self.mutation_map)
        self.arg = argbook.ARG()
        self.arg.nextname = self.n
        self.parent_nodes = collections.defaultdict(list)
        for k in range(self.n):
            node = self.arg.alloc_node(k, 0)#index, time,
            samples = bintrees.AVLTree()
            samples.__setitem__(k, k)
            s = self.arg.alloc_segment(0, math.ceil(self.seq_length), node, samples)
            node.first_segment = s
            self.add_mutation(node)
            node = self.arg.add(node)
            x = self.arg.alloc_segment(0, math.ceil(self.seq_length), node, samples)
            self.parent_nodes[k] = x

    def add_mutation(self, node):
        '''put the mutations on the node'''
        if self.mutation_map[node.index]:
            node.snps.update({k:k for k in self.mutation_map[node.index]})

    def ts_to_argnode(self):
        while self.edges_dict:
            parent = next(iter(self.edges_dict))
            time = self.tables.nodes[parent].time
            if self.tables.nodes[parent].flags == msprime.NODE_IS_RE_EVENT:
                pass
            else: # CA
                child0 = self.edges_dict[parent][0].child
                child1 = self.edges_dict[parent][-1].child
                assert child0 != child1
                self.commom_ancestor_event(time, child0, child1, parent)
                del  self.edges_dict[parent]

    def argnode_to_ts(self):
        '''TODO: get ts_full from argnode'''

    def commom_ancestor_event(self, time, child0, child1, p):
        x = self.parent_nodes[child0]
        y = self.parent_nodes[child1]
        x = x.get_first_segment()
        y = y.get_first_segment()
        assert x is not None
        assert y is not None
        index =  self.arg.new_name()
        node = self.arg.alloc_node(index, time, x.node, y.node)
        if self.mutation_map[node.index]:
            self.add_mutation(node)
            self.arg.push_all_mutations_down(node)
        # self.C[node.index] = node.index
        self.arg.coal[node.index] = node.index
        x.node.left_parent = node
        x.node.right_parent = node
        y.node.left_parent = node
        y.node.right_parent = node
        z = None
        defrag_required = False
        while x is not None or y is not None:
            alpha = None
            if x is None or y is None:
                if x is not None:
                    alpha = x
                    x = None
                    assert alpha.left < alpha.right
                if y is not None:
                    alpha = y
                    y = None
                    assert alpha.left < alpha.right
            else:
                if y.left < x.left:
                    beta = x
                    x = y
                    y = beta
                if x.right <= y.left:
                    alpha = x
                    x = x.next
                    alpha.next = None
                    assert alpha.left < alpha.right
                elif x.left != y.left:
                    alpha = self.arg.alloc_segment(x.left, y.left, node, x.samples)
                    x.left = y.left
                    assert alpha.left < alpha.right
                else:
                    left = x.left
                    r_max = min(x.right, y.right)
                    right = r_max
                    alpha = self.arg.alloc_segment(left, right, node, x.union_samples(y))
                    assert alpha.left < alpha.right
                    if alpha.is_mrca(self.n):
                        alpha = None
                    if x.right == right:
                        x = x.next
                    else:
                        x.left = right
                    if y.right == right:
                        y = y.next
                    else:
                        y.left = right
            if alpha is not None:
                if z is None:
                    self.parent_nodes[p] = alpha
                if z is not None:
                    defrag_required |= z.right == alpha.left
                    z.next = alpha
                alpha.prev = z
                z = alpha
        assert node is not None
        if z is not None:
            z = z.get_first_segment()
            node.first_segment = z
            self.arg.store_node(z, node)
        else:
            self.arg.add(node)
            self.arg.roots[node.index] = node.index

    def map_mutation(self):
        '''returns a list of mutation positions  where the index
        is the node on which the mutation took place
        '''
        mutation_map = [[] for _ in range(self.ts_full.num_nodes)]
        position = self.ts_full.tables.sites.position
        site = self.ts_full.tables.mutations.site
        node = self.ts_full.tables.mutations.node
        for mutation_id in range(self.ts_full.num_mutations):
            site_position = position[site[mutation_id]]
            mutation_map[node[mutation_id]].append(math.ceil(site_position))
        return mutation_map

def get_ts_genotype(ts_full):
    return ts_full.genotype_matrix()

def get_arg_genotype(ts_full):
    '''
    get the genotypes from ts and implements a discretising function,
     which  rounds upwards to the nearest int, and store in a dict where keys
    are SNP positions and values are the sequences for which
    the allele is derived. ANDing if the results if 2 or more variants end up at the same
    integer position.
    '''
    simple_ts = ts_full.simplify()
    new_data = {}
    derived = 1
    genotypes = None
    position = 0
    for v in simple_ts.variants():
        if int(math.ceil(v.position)) != position:
            genotypes = v.genotypes
            position = int(math.ceil(v.position))
        else:
            raise NameError("Not two SNPs can have identical genomic position")
            # genotypes = np.logical_and(genotypes, v.genotypes)
        b = bintrees.AVLTree()
        b.update({k: k for k in list(np.where(genotypes == derived)[0])})
        new_data[math.ceil(position)] = b
    return new_data
