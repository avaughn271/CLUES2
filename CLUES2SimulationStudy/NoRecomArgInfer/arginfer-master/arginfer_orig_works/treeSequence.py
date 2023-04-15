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

    def find_break(self, z,break_point):
        while z.prev is not None:
            z = z.prev
        while z is not None:
            if z.contains(break_point):
                return z
            z = z.next

    def add_mutation(self, node):
        '''put the mutations on the node'''
        if self.mutation_map[node.index]:
            node.snps.update({k:k for k in self.mutation_map[node.index]})

    def ts_to_argnode(self):
        while self.edges_dict:
            parent = next(iter(self.edges_dict))
            time = self.tables.nodes[parent].time
            if self.tables.nodes[parent].flags == msprime.NODE_IS_RE_EVENT:
                # print("R"*100)
                child = self.edges_dict[parent][0].child
                parent2 = parent + 1#self.edges_dict[parent + 1] # the other rec parent
                assert self.tables.nodes[parent2].time == time
                # find breakpoint
                if self.edges_dict[parent][-1].right == self.edges_dict[parent2][0].left:
                    l_break = math.ceil(self.edges_dict[parent][-1].right)
                    r_break = None
                else:
                    r_break = math.ceil(self.edges_dict[parent2][0].left)
                    l_break = math.ceil(self.edges_dict[parent][-1].right)
                self.recombination_event(time, l_break, r_break,
                                         parent, parent2, child)
                del self.edges_dict[parent]
                del self.edges_dict[parent2]
            else: # CA
                child0 = self.edges_dict[parent][0].child
                child1 = self.edges_dict[parent][-1].child
                assert child0 != child1
                self.commom_ancestor_event(time, child0, child1, parent)
                del  self.edges_dict[parent]

    def argnode_to_ts(self):
        '''TODO: get ts_full from argnode'''

    def recombination_event(self, time, l_break,
                            r_break, p1, p2, child):
        s = self.parent_nodes[child]
        assert s is not None
        if r_break == None:#  ancestral REC
            y = self.find_break(s, l_break)
            y.node.breakpoint = l_break
            if y.prev is not None:
                assert y.left <= l_break
            else:
                assert y.left < l_break
            z = self.arg.alloc_segment(l_break, y.right, y.node,
                                       y.samples, None, y.next)
            assert l_break < y.right
            if y.next is not None:
                y.next.prev = z
            y.next = None
            y.right = l_break
            assert  y.left < l_break
            lhs_tail = y
        else: # non ancestral REC
            y = self.find_break(s, r_break)
            x = y.prev
            break_point = random.choice(range(x.right, y.left + 1))
            y.node.breakpoint = break_point
            assert x.right <= break_point <= y.left
            x.next = None
            y.prev = None
            z = y
            lhs_tail = x
        assert z is not None
        assert lhs_tail is not None
        self.parent_nodes[p1] = lhs_tail
        self.parent_nodes[p2] = z
        # if z.node.left_child is not None and z.node.left_child.index == z.node.right_child.index:
        #     # If the child node of the recombination is itself a recombination node,
        #     # then recombination being in progress renders the child node non-removable
        #     del self.R[z.node.index]
        # TODO: Implement a way to check whether a recombination is removable, and
        # delete those recombination nodes which are not. Will have to be done once the
        # ARG has been constructed by checking which ones will violate snps

        node = self.arg.alloc_node(self.arg.new_name(), time,
                                   lhs_tail.node, lhs_tail.node)
        lhs_tail.node.left_parent = node
        if self.mutation_map[node.index]:
            self.add_mutation(node)
            self.arg.push_all_mutations_down(node)
        self.arg.store_node(lhs_tail, node)
        # node = self.arg.add(node1)
        # self.R[node.index] = node.index
        self.arg.rec[node.index] = node.index
        node = self.arg.alloc_node(self.arg.new_name(), time,  z.node, z.node)
        if self.mutation_map[node.index]:
            self.add_mutation(node)
            self.arg.push_all_mutations_down(node)
        z.node.right_parent = node
        self.arg.store_node(z, node)
        # self.R[node.index] = node.index
        self.arg.rec[node.index] = node.index

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
        if defrag_required:
            z.defrag_segment_chain()
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

    def print_state(self):
        print("node", "time", "left", "right", "l_chi", "r_chi", "l_par", "r_par",
              "l_bp", "snps", "fir_seg_sam",
              sep="\t")
        for j in self.arg.nodes:
            node = self.arg.nodes[j]
            if node.left_parent is not None or node.left_child is not None:
                s = node.first_segment
                if s is None:
                    print(j, "%.5f" % node.time, "root", "root",
                              node.left_child.index,
                              node.right_child.index, "Root", "Root",
                              node.breakpoint,
                              node.snps ,None, sep="\t")

                while s is not None:
                    l = s.left
                    r = s.right
                    if node.left_child is None:
                        print(j, "%.5f" % node.time, l,r, "Leaf", "Leaf",
                              node.left_parent.index,node.right_parent.index,
                              node.breakpoint,
                              node.snps,  s.samples,  sep="\t")#
                    elif  node.left_parent is None:
                        print(j, "%.5f" % node.time, l, r,
                              node.left_child.index,
                              node.right_child.index, "Root", "Root",
                              node.breakpoint,
                              node.snps ,s.samples, sep="\t")
                    else:
                        # print(node.time, node.index, int(l), int(r), node.left_child.index, node.right_child.index)
                        print( j, "%.5f" % node.time, l, r,
                             node.left_child.index, node.right_child.index,
                              node.left_parent.index, node.right_parent.index,
                              node.breakpoint,
                              node.snps, s.samples, sep="\t")
                    s = s.next

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

# def get_arg_genotype(ts_full):
#     '''
#     get the genotyped from ts and store in a dict where keys
#     are SNP positions and values are the sequences for which
#     the allele is derived.
#     '''
#     derived = 1
#     ts_data = get_ts_genotype(ts_full)
#     # print(ts_data)
#     snp_positions = ts_full.tables.sites.position
#     new_data = {}#collections.defaultdict(list)
#     for ind in range(len(ts_data)):
#         if math.ceil(snp_positions[ind]) in new_data:
#             raise NameError("Not two SNPs can have identical genomic position")
#         else:
#             b = bintrees.AVLTree()
#             b.update({k: k for k in list(np.where(ts_data[ind] == derived)[0])})
#             new_data[math.ceil(snp_positions[ind])] = b
#                 #(np.where(ts_data[ind] == derived)[0])#list(np.where(ts_data[ind] == derived)[0]) # derived
#     return new_data

def test_run():
    recombination_rate=1e-8
    Ne= 5000
    sample_size = 5
    length = 2e5
    ts_full = msprime.simulate(sample_size = sample_size, Ne = Ne, length = length, mutation_rate = 1e-8,
                        recombination_rate = recombination_rate, random_seed = 20, record_full_arg = True)

    # print(ts_full.tables.edges)
    tsarg= TreeSeq(ts_full)
    tsarg.ts_to_argnode()
    tsarg.print_state()

if __name__ == "__main__":
    test_run()




