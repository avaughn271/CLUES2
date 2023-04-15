'''
build an arg from the given data
use the built arg as the initial ARG for the MCMC
IDEA:
1. choose a time for the next event from exp(rate)
    rate = CwR rate for time
2. choose an event "CA"or "REC" proportional to their rate
3. If a "REC":
     choose a lineage proportional to its links
        and split from a random breakpoint
4. else: randomly choose a lineage (A), among all
        nodes compatible with A, choose one with
        highest overlapping and merge them. If there is none,
        reserve the time and go back to 2.
'''
import random
from operator import itemgetter
import argbook
from argbook import *

class Initial(object):

    def __init__(self, data, sample_size, seq_length,
                 Ne, mutation_rate, recombination_rate):
        self.data = data # must be in our  format {snp:AVLTree(sequences)}
        self.n = sample_size
        self.seq_length = seq_length
        self.Ne = Ne
        self.mu = mutation_rate
        self.r = recombination_rate
        self.number_of_lineages = sample_size
        self.number_of_links = self.number_of_lineages * (self.seq_length - 1)
        self.arg = argbook.ARG()
        self.arg.nextname = self.n
        self.P = collections.defaultdict(list)
        self.remaining_SNPs = bintrees.AVLTree()
        self.remaining_SNPs.update({k: k for k in self.data.keys()})# all snps for now
        for k in range(self.n):
            node = self.arg.alloc_node(k, 0)#index, time,
            samples = bintrees.AVLTree()
            samples.__setitem__(k, k)
            s = self.arg.alloc_segment(0, math.ceil(self.seq_length), node, samples)
            node.first_segment = s
            node = self.arg.add(node)
            x = self.arg.alloc_segment(0, math.ceil(self.seq_length), node, samples)
            self.P[k] = x
        # put the singletone mutations on the nodes
        for k in self.data.keys():
            if len(self.data[k]) == 1:
                self.arg.nodes[self.data[k].min_key()].snps.__setitem__(k, k)
                self.remaining_SNPs.discard(k)
        self.t = 0
        self.next_event = True

    def generate_time(self):
        c=1
        self.number_of_lineages = len(self.P)
        # print("self.number_of_lineages", self.number_of_lineages)
        coal_rate = (self.number_of_lineages * (self.number_of_lineages - 1)
                    / (4*self.Ne))
        rec_rate = (self.number_of_links/c * (self.r))
        # print("self.-------------------number_of_links", self.number_of_links)
        tot_rate = coal_rate + rec_rate
        if self.next_event:
            self.t += random.expovariate(tot_rate)
        if random.random()<= coal_rate/tot_rate:
            return "coal"
        else:
            return "rec"

    def get_node_variants(self, u):
        '''return the SNPs on a node that havent mutated yet
        :param u is  a node'''
        node_variants = bintrees.AVLTree()
        seg = u.first_segment
        while seg is not None:
            seg_snps = self.remaining_SNPs[seg.left:seg.right]
            if seg_snps:
                node_variants = node_variants.union(seg_snps)
            seg = seg.next
        return node_variants

    def compatibility_check(self,S1,  S2, snp, new_muts):
        '''
        if the coalescence of two nodes with samples S1 and S2
        is compatible for  snp .
        S1: is the node1 samples for this segment
        S2: is the node2 samples for this segment
        snp: the focal SNP
        new_muts: an AVLTree on new mutations on the
            parent of node1 and node2.
        '''
        ret = True
        D = self.data[snp]
        # symmetric difference between S1  and D
        A = S1.union(S2)
        symA_D = A.difference(D)
        if len(symA_D) == 0:# subset or equal
            if len(A) == len(D): #put the mutation on this node
                new_muts.__setitem__(snp, snp)
        elif len(symA_D) == len(A): # distinct
            pass
        else:#
            symD_A = D.difference(A)
            if len(symD_A) > 0: # incompatible
                ret = False
        return ret, new_muts

    def two_seqs_are_consistent(self, u, v):
        '''check if two sequences are consistent
         according to infinite sites model
         1. find the variants of nodes u and v from the
            remaining_SNPs (those that havent mutated yet)
        2. for the intersect sites check compatibility
        :return True/False consistent, and the new mutations
            if they are consistent.
         '''
        consistent = True
        new_muts=bintrees.AVLTree()
        if self.remaining_SNPs:
            #1. get the variants of both nodes
            u_variants = self.get_node_variants(u.node)
            v_variants = self.get_node_variants(v.node)
            #2. get the intersects
            intersect_variants = u_variants.intersection(v_variants)
            # now check if these sites are consistent

            for snp in intersect_variants:
                S1 = u.node.x_segment(snp).samples
                S2 = v.node.x_segment(snp).samples
                consistent, new_muts = self.compatibility_check(S1,
                                                    S2, snp, new_muts)
                if not consistent:
                    break
        return consistent, new_muts

    def commom_ancestor_event(self, x, y):
        assert x is not None
        assert y is not None
        index =  self.arg.new_name()
        node = self.arg.alloc_node(index, self.t, x.node, y.node)
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
                    alpha = self.arg.alloc_segment(x.left, y.left,
                                                   node, x.samples)
                    x.left = y.left
                    assert alpha.left < alpha.right
                else:
                    left = x.left
                    r_max = min(x.right, y.right)
                    right = r_max
                    alpha = self.arg.alloc_segment(left, right,
                                                   node, x.union_samples(y))
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
                    self.P[node.index] = alpha
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
        return node

    def merge_event(self):
        ''''''
        u1_index = random.sample(self.P.keys(), 1)[0]
        u1_head = self.P.pop(u1_index)
        u1_tail = u1_head.node.get_tail()
        # now find another seq in P that u1 can join to--> is compatible
        #TODO: parallele
        all_consistents = []#u2_index, new_muts, overlapping_degree
        for ind in self.P.keys():
            consistent, new_muts = self.two_seqs_are_consistent(u1_head, self.P[ind])
            if consistent:
                u2_index = ind
                u2_head = self.P[ind]
                u2_tail = u2_head.node.get_tail()
                overlapping = min(u1_tail.right, u2_tail.right)-\
                    max(u1_head.left, u2_head.left)
                all_consistents.append([u2_index, new_muts, overlapping])
                if overlapping >=  u1_tail.right- u1_head.left:
                    break
        if all_consistents:
            #take the one with highest overalpping
            all_consistents = sorted(all_consistents, key=itemgetter(2), reverse=True)
            # u1_index and u2_head can merge
            chosen = all_consistents[0]
            #----- added this if: 28 Sep 2020--> only overlapping can coal
            if chosen[2]>0:
                u2_index = chosen[0]; new_muts = chosen[1]
                u2_head = self.P.pop(u2_index)
                node = self.commom_ancestor_event(u1_head, u2_head)
                node.snps = new_muts
                #remove the mutated one from SNPs
                self.remaining_SNPs = self.remaining_SNPs.difference(new_muts)
                #---------
                self.number_of_links -= self.arg.__getitem__(u1_index).num_links()
                self.number_of_links -= self.arg.__getitem__(u2_index).num_links()
                if node.first_segment != None:
                    self.number_of_links += node.num_links()
                self.next_event = True
            else:
                self.P[u1_index] = u1_head
                self.next_event = False
        else:
            #there is no node consistent with u1, then no coal
            # put u1 back to P
            self.P[u1_index] = u1_head
            self.next_event = False

    def build(self):
        count =0
        while self.P:
            event_type = self.generate_time()
            if event_type == "coal":
                self.merge_event()
            else:
                count+=1
                print("problem18")
        assert self.remaining_SNPs.is_empty
        assert self.number_of_links == 0