''' This module is responsible for ARG classes'''
import collections
import math
from sortedcontainers import SortedSet
from sortedcontainers import SortedList
import bintrees
import pickle
import numpy as np
import pandas as pd

class Segment(object):
    """
    A class representing a single segment. Each segment has a left and right, denoting
    the loci over which it spans, a node giving the node to which it belongs in an ARG,
     a prev and next, giving the previous and next segments in the chain,
     and samples representing the samples underneath the segment.
    """
    def __init__(self):
        self.left = None
        self.right = None
        self.node = None
        self.prev = None
        self.next = None
        self.samples = bintrees.AVLTree()

    def __str__(self):
        s = "({}:{}-{}->{}: prev={} next={})".format(
            self.left, self.right, self.node.index, self.samples,repr(self.prev),
            repr(self.next))
        return s

    def __lt__(self, other):
        return ((self.left, self.right)
                < (other.left, other.right))

    def copy(self):
        ''':return a copy of this segment'''
        s = Segment()
        s.left = self.left
        s.right = self.right
        s.node = self.node
        s.samples = self.samples
        s.next = self.next
        s.prev = self.prev
        return s

    def contains(self, x):
        return x >= self.left and x < self.right

    def is_mrca(self, sample_size):#V3
        assert len(self.samples) <= sample_size
        return len(self.samples) == sample_size

    def equal_samples(self, other):
        return sorted(self.samples) == sorted(other.samples)

    def union_samples(self, other):
        return self.samples.union(other.samples)

    def equal(self, other):
        if self.node is None and other.node is None:
            return self.left == other.left and \
                   self.right == other.right and self.node == other.node and \
                   sorted(self.samples) == sorted(other.samples)
        elif self.node is not None and other.node is not None:
            return self.left == other.left and \
                   self.right == other.right and self.node.index == other.node.index and \
                   sorted(self.samples) == sorted(other.samples)
        else:
            return False

    def get_first_segment(self):
        '''get the fist segment in a chain of segments'''

        seg = self
        while seg.prev is not None:
            seg = seg.prev
        return seg

    def get_variants(self, data):
        '''get all snps of data that are in  self
        TODO: find an efficient way
        '''
        seg_variants = bintrees.AVLTree()
        for item in data.keys():
            if self.contains(item):
                seg_variants[item] = item
        return seg_variants

    def add_snps(self, node , data):
        '''add snps to node that are in this segment'''
        for snp in data.keys():
            if self.contains(snp) and \
                    sorted(self.samples) == sorted(data[snp]):
                node.snps[snp] = snp

    def get_seg_variants(self, data):
        '''TODO: implement efficiently'''
        assert self is not None
        seg_variants = bintrees.AVLTree()
        for item in data.keys():
            if self.contains(item):
                seg_variants[item] = item
        return seg_variants

    def get_intersect(self, start, end):
        '''
        return the intersection of self and [start, end)
        '''
        ret = []
        if self.left <= start and self.right > start:
            if self.right<= end:
                ret = [start, self.right]
            else:
                ret = [start, end]
        elif self.left > start and self.left < end:
            if self.right<=end:
                ret = [self.left, self.right]
            else:
                ret= [self.left, end]
        return ret

class Node(object):
    """
    A class representing a single node. Each node has a left and right child,
    and a left and right parent. If a node arises from a recombination, then
    left_child == right_child, while if it ends in a coalescence, then
    left_parent == right_parent. Each node also has a time at which it
    appears in the tree, and a list of segments of ancestral material, and a list of
    snps. The snps represent mutations arising on the branch of which the node is a
    child. The fact that snps are stored on a single branch is for computational
    convenience; the MCMC algorithm marginalises over all contiguous branches which
    subtend the same leaves at the site of the snp, and hence the snp could just as well
    be stored on any one of them.
    """

    def __init__(self, index):
        self.left_child = None
        self.right_child = None
        self.left_parent = None
        self.right_parent = None
        self.first_segment = None
        self.snps = bintrees.AVLTree()
        self.time = None
        self.breakpoint = None
        self.index = index

    def copy(self):
        '''a copy of the node'''
        cpy_node = Node(self.index)
        s = self.first_segment
        if s is not None:
            x = s.copy()
            cpy_node.first_segment = x
            x.node = cpy_node
            x = x.next
            while x is not None:
                s = x.copy()
                s.prev.next = s
                x.node = cpy_node
                x = x.next
        else:
            cpy_node.first_segment = None
        cpy_node.time = self.time
        cpy_node.breakpoint = self.breakpoint
        cpy_node.snps = self.snps.copy()
        return cpy_node

    def contains(self, x):
        seg = self.first_segment
        while seg is not None:
            if seg.contains(x):
                return True
            seg = seg.next
        return False

    def x_segment(self, x):
        '''return the segment containing x
        given that we know self contains x'''
        seg = self.first_segment
        while seg is not None:
            if seg.contains(x):
                return seg
            seg = seg.next
        raise ValueError("x is not in node")

    def num_links(self):
        seg = self.first_segment
        left = seg.left
        while seg.next is not None:
            seg = seg.next
        return seg.right - left -1

    def is_leaf(self):
        return self.left_child == None

    def is_root(self):
        return  self.left_parent == None

    def equal(self, other):
        if self is not None and other is not None:
            if self.time != other.time or\
                    sorted(self.snps) != sorted(other.snps) or\
                    self.index != other.index:
                return False
            else:
                seg = self.first_segment
                sego = other.first_segment
                while seg is not None and sego is not None:
                    if not seg.equal(sego):
                        return False
                    seg = seg.next
                    sego = sego.next
                if seg is None and sego is None:
                    return True
                else:
                    return False
        else:
            raise ValueError("one or both nodes are NONE")

    def arg_node_age(self):
        '''the arg branch length of a node '''
        if self.left_parent is not None:
            return self.left_parent.time - self.time
        else:
            return 0

    def upward_path(self, x):
        '''for position x check if we can move upward.
        this is used in finding the branch length at
        position x in a tree'''
        if self.left_parent is None:
            block = True
            return self, block
        elif self.left_parent.index is not self.right_parent.index:
            if self.left_parent.contains(x) + self.right_parent.contains(x) != 1:
                print("in upward_path x is", x, "left_aprent", self.left_parent.index,
                      "right parent",self.right_parent.index, "node", self.index)
            assert self.left_parent.contains(x) + self.right_parent.contains(x) == 1
            block = False
            if self.left_parent.contains(x):
                return self.left_parent, block
            else:
                return self.right_parent, block
        else:#CA
            sib = self.sibling()
            #--- after spr before clean up, sib might be NAM
            if sib.first_segment != None and sib.contains(x):
                block = True
                return self.left_parent, block
            else:
                block = False
                return self.left_parent, block

    def tree_node_age(self, x, return_parent_time= False):
        '''
        the tree branch length of
        node self, at position x
        :param x the site
        :param return_parent_time: if we only want to
            report parent time ---> in the case of alelle age
         '''
        node = self
        child_time = node.time
        block = False
        while not block:
            node, block = node.upward_path(x)
        assert node.time - child_time > 0
        if not return_parent_time:
            return node.time - child_time
        else:
            return node.time

    def sibling(self):
        '''
        Find and return the sibling node of u
        where u is a child of a CA
        '''
        assert self.left_parent is not None
        assert self.left_parent.index == self.right_parent.index
        p = self.left_parent
        v = p.left_child
        if v.index == self.index:
            v = p.right_child
        assert v.index is not self.index
        return v

    def push_snp_down(self, x):
        # Push the snp at position x down one branch from node to one of its children
        # provided only one is ancestral at x.
        if self.left_child is None:
            block = True
            return self, block
        elif self.left_child is not self.right_child:
            if self.left_child.contains(x) and self.right_child.contains(x):
                block = True
                return self, block
            elif self.left_child.contains(x):
                self.left_child.snps.__setitem__(x, x)
                self.snps.discard(x)
                block = False
                return self.left_child, block
            else:
                self.right_child.snps.__setitem__(x, x)
                self.snps.discard(x)
                block = False
                return self.right_child, block
        else:# rec
            self.left_child.snps.__setitem__(x, x)
            self.snps.discard(x)
            block = False
            return self.left_child, block

    def get_tail(self):
        seg = self.first_segment
        while seg.next is not None:
            seg = seg.next
        return seg

    def get_variants(self, data):
        '''get all snps in data lay in the self segments
        TODO: an efficient way
        it is not efficient ot loop over data SNPs for each segment
        '''
        node_variants = bintrees.AVLTree()
        seg = self.first_segment
        while seg is not None:
            for item in data.keys():
                if seg.contains(item):
                    node_variants[item] = item
            seg = seg.next
        return node_variants

    def update_child(self, oldchild, newchild):
        '''update self child from oldchild to newchild'''
        if self.left_child != None:
            if self.left_child.index == oldchild.index:
                self.left_child = newchild
        if self.right_child != None:
            if self.right_child.index == oldchild.index:
                self.right_child = newchild

    def reconnect(self, child):# BUG7
        '''from child--> self--> parent: TO child ---> parent '''
        leftparent = self.left_parent
        rightparent = self.right_parent
        child.left_parent = leftparent
        child.right_parent = rightparent
        child.breakpoint = self.breakpoint
        leftparent.update_child(self, child)
        rightparent.update_child(self, child)

class ARG(object):
    '''
    Ancestral Recombination Graph
    '''
    def __init__(self):
        self.nodes = {}
        self.roots = bintrees.AVLTree()# root indexes
        self.rec = bintrees.AVLTree() # arg rec parents nodes
        self.coal = bintrees.AVLTree() # arg CA parent node
        self.num_ancestral_recomb = 0
        self.num_nonancestral_recomb = 0
        self.branch_length = 0
        self.nextname = 1 # next node index
        self.available_names = SortedSet()

    def __iter__(self):
        '''iterate over nodes in the arg'''
        return list(self.nodes)

    def __len__(self):
        '''number of nodes'''
        return len(self.nodes)

    def __getitem__(self, index):
        '''returns node by key: item'''
        return self.nodes[index]

    def __setitem__(self, index, node):
        '''adds a node to the ARG'''
        node.index = index
        self.add(node)

    def __contains__(self, index):
        '''if ARG contains node key '''
        return index in self.nodes

    def copy(self):
        '''return a copy of the ARG'''
        arg = ARG()
        for node in self.nodes.values():
            arg.nodes[node.index] = node.copy()
        # connect nodes
        for node in self.nodes.values():
            node2 = arg.__getitem__(node.index)
            if node.left_child != None:
                node2.left_child = arg.__getitem__(node.left_child.index)
                node2.right_child = arg.__getitem__(node.right_child.index)
            if node.left_parent != None:
                node2.left_parent = arg.__getitem__(node.left_parent.index)
                node2.right_parent = arg.__getitem__(node.right_parent.index)
        arg.roots = self.roots.copy()# root indexes
        arg.rec = self.rec.copy()# arg rec parents nodes
        arg.coal = self.coal.copy() # arg CA parent node
        arg.num_ancestral_recomb = self.num_ancestral_recomb
        arg.num_nonancestral_recomb = self.num_nonancestral_recomb
        arg.branch_length = self.branch_length
        arg.nextname = self.nextname # next node index
        arg.available_names = self.available_names.copy()
        return arg

    def equal(self, other):
        '''if self is equal with other (structural equality)
        TODO : complete this'''
        if self.__len__() != other.__len__():
            return False
        else:
            for node in self.nodes.values():
                if node.index not in other:
                    return False
                if not node.equal(other[node.index]):
                    return False
            return True

    def leaves(self, node=None):
        """
        Iterates over the leaves of the ARG.
        """
        if node is None:
            for node in self.nodes.values():
                if node.left_child == None:
                    yield node
        else:
            for node in self.preorder(node):
                if node.left_child == None:
                    yield node

    def preorder(self, node=None):
        """
        Iterates through nodes in preorder traversal.
        """
        visit = set()
        if node is None:
            node = self.__getitem__(self.roots.max_key())
        queue = [node]
        for node in queue:
            if node in visit:
                continue
            yield node
            visit.add(node)
            if node.left_child != None:
                queue.append(node.left_child)
                if node.left_child.index != node.right_child.index:
                    queue.append(node.right_child)

    def get_times(self):
        '''return a sorted set of the ARG node.time'''
        times = SortedSet()
        for node in self.nodes.values():
            times.add(node.time)
        return times

    #==========================
    # node manipulation
    def alloc_segment(self, left = None, right = None, node = None,
                      samples = bintrees.AVLTree(), prev = None, next = None):
        """
        alloc a new segment
        """
        s = Segment()
        s.left = left
        s.right = right
        s.node = node
        s.samples = samples
        s.next = next
        s.prev = prev
        return s

    def alloc_node(self, index = None, time = None,
                    left_child = None, right_child = None):
        """
        alloc a new Node
        """
        node = Node(index)
        node.time = time
        node.first_segment = None
        node.left_child = left_child
        node.right_child = right_child
        node.left_parent = None
        node.right_parent = None
        node.breakpoint = None
        node.snps = bintrees.AVLTree()
        return node

    def store_node(self, segment, node):
        '''store node with segments: segment'''
        x = segment
        if x is not None:
            while x.prev is not None:
                x = x.prev
            s = self.alloc_segment(x.left, x.right, node, x.samples.copy())
            node.first_segment = s
            x.node = node
            x = x.next
            while x is not None:
                s = self.alloc_segment(x.left, x.right, node, x.samples.copy(), s)
                s.prev.next = s
                x.node = node
                x = x.next
        else:#
            node.first_segment = None
        self.nodes[node.index] = node

    def copy_node_segments(self, node):
        '''
        copy the segments of a node,
        in CA event or Rec events, we need to copy the first node
        in order to make changes on them
        '''
        x = node.first_segment
        if x is None:
            return None
        else:
            assert x.prev is None
            s = self.alloc_segment(x.left, x.right, node, x.samples.copy())
            x.node = node
            x = x.next
            while x is not None:
                s = self.alloc_segment(x.left, x.right, node, x.samples.copy(), s)
                s.prev.next = s
                x.node = node
                x = x.next
            return s

    def get_available_names(self):
        '''get free names from 0 to max(nodes)'''
        self.available_names = SortedSet()
        current_names = SortedSet(self.__iter__())
        counter = 0
        prev = current_names[0]
        while counter < len(current_names):
            if current_names[counter] != prev + 1:
                self.available_names.update(range(prev+1, current_names[counter]))
            prev = current_names[counter]
            counter += 1

    def new_name(self):
        '''returns a new name for a node'''
        if self.available_names:
            name = self.available_names.pop()
        else:
            name = self.nextname
            self.nextname += 1
        return name

    def add(self, node):
        ''' add a ready node to the ARG:
        '''
        self.nodes[node.index] = node
        return node

    def rename(self, oldindex, newindex):
        '''renames a node in the ARG'''
        node = self.nodes[oldindex]
        node.index = newindex
        del self.nodes[oldindex]
        self.nodes[newindex] = node

    def total_branch_length(self):
        '''the ARG total branch length'''
        total_material = 0
        for node in self.nodes.values():
            if node.left_parent is not None:
                age = node.left_parent.time - node.time
                seg = node.first_segment
                while seg is not None:
                    total_material += ((seg.right - seg.left)* age)
                    seg = seg.next
        return total_material
    #=======================
    #spr related

    def detach(self, node, sib):
        '''
        Detaches a specified coalescence node from the rest of the ARG
        '''
        # print("Detach()",node.index, "sib", sib.index, "p",node.left_parent.index)
        assert node.left_parent.index == node.right_parent.index
        parent = node.left_parent
        sib.left_parent = parent.left_parent
        sib.right_parent = parent.right_parent
        sib.breakpoint = parent.breakpoint
        grandparent = parent.left_parent
        if grandparent is not None:
            grandparent.update_child(parent, sib)
            grandparent = parent.right_parent
            grandparent.update_child(parent, sib)

    def reattach(self, u, v, t, new_names):
        # Reattaches node u above node v at time t, new_names is a avltree of all
        #new nodes.index in a new ARG in mcmc
        assert t > v.time
        # assert v.left_parent == None or t < v.left_parent.time
        if u.left_parent is None:# new_name
            new_name = self.new_name()
            new_names[new_name] = new_name
            # self.coal[new_name] = new_name # add the new CA parent to the ARG.coal
            parent = self.add(self.alloc_node(new_name))
            parent.left_child = u
            u.left_parent = parent
            u.right_parent = parent
        else:
            assert u.left_parent.index == u.right_parent.index
            parent = u.left_parent
        parent.time = t
        parent.breakpoint = v.breakpoint
        v.breakpoint = None
        parent.left_parent = v.left_parent
        grandparent = v.left_parent
        if grandparent is not None:
            grandparent.update_child(v, parent)
        parent.right_parent = v.right_parent
        grandparent = v.right_parent
        if grandparent is not None:
            grandparent.update_child(v, parent)
        v.left_parent = parent
        v.right_parent = parent
        if parent.left_child.index == u.index:
            parent.right_child = v
        else:
            parent.left_child = v
        return new_names

    def push_mutation_down(self, node, x):
        '''
        for a given node push the mutation (at x) as down as possible
        normally mutations automatically should stay at their
        lowest possible position. This might be useful for initial ARG
        '''
        block = False
        while not block:
            node, block = node.push_snp_down(x)

    def push_all_mutations_down(self, node):
        '''push down all mutations on node as low as possible'''
        snp_keys = [k for k in node.snps]
        for x in snp_keys:
            self.push_mutation_down(node, x)

    #========== probabilites
    def log_likelihood(self, mutation_rate, data):
        '''
        log_likelihood of mutations on a given ARG up to a normalising constant
         that depends on the pattern of observed mutations, but not on the ARG
         or the mutation rate.
         Note after spr and berfore clean up we might have NAM lineages,
         this method covers take this into account.
         :param m : is number of snps
         '''
        snp_nodes = [] # nodes with len(snps) > 0
        total_material = 0
        number_of_mutations = 0
        #get total matereial and nodes with snps
        for node in self.nodes.values():
            if node.first_segment != None:
                assert node.left_parent != None
                age = node.left_parent.time - node.time
                seg = node.first_segment
                assert seg.prev == None
                while seg is not None:
                    total_material += ((seg.right - seg.left)* age)
                    seg = seg.next
                if node.snps:
                    number_of_mutations += len(node.snps)
                    snp_nodes.append(node)
        self.branch_length = total_material
        # print("number_of_mutations", number_of_mutations, "m", len(data))
        assert number_of_mutations == len(data) # num of snps
        if mutation_rate == 0:
            if number_of_mutations == 0:
                ret = 0
            else:
                ret = -float("inf")
        else:
            ret = number_of_mutations * math.log(total_material * mutation_rate) -\
                (total_material * mutation_rate)
        # now calc prob of having this particular mutation pattern
        for node in snp_nodes:
            # num_branches = collections.defaultdict(list)
            for x in node.snps:
                potential_branch_length = node.tree_node_age(x)
                ret += math.log(potential_branch_length / total_material)
        return ret

    def log_prior(self, sample_size, sequence_length, recombination_rate, Ne,
                  NAM = True, new_roots = False , kuhner = False):
        '''
        probability of the ARG under coalescen with recombination
        this is after a move and before clean up. then there might be some
         extra NAM lineages, we ignore them.
         :param NAM: no-ancestral material node. If NAm node is allowed. note after spr and
            before clean up step there might be some NAM in the ARG which is ok. But after clean up
            or on the initial ARG there should not be any.
         '''
        # order nodes by time
        #TODO: find an efficient way to order nodes
        ordered_nodes = [v for k, v in sorted(self.nodes.items(),
                                     key = lambda item: item[1].time)]
        number_of_lineages = sample_size
        number_of_links = number_of_lineages * (sequence_length - 1)
        number_of_nodes = self.__len__()
        counter = sample_size
        time  = 0
        ret = 0
        rec_count = 0
        coal_count = 0
        roots = bintrees.AVLTree()
        new_coal = bintrees.AVLTree()
        if kuhner:
            self.rec.clear()
        self.num_ancestral_recomb = 0
        self.num_nonancestral_recomb = 0
        while counter < number_of_nodes:
            node = ordered_nodes[counter]
            assert node.time >= time # make sure it is ordered]
            rate = (number_of_lineages * (number_of_lineages - 1)
                    / (4*Ne)) + (number_of_links * (recombination_rate))
            # ret -= rate * (node.time - time)
            if node.left_child.index == node.right_child.index: #rec
                assert node.left_child.first_segment != None
                assert node.left_child.left_parent.first_segment != None
                assert node.left_child.right_parent.first_segment != None
                ret -= rate * (node.time - time)
                gap = node.left_child.num_links()-\
                      (node.left_child.left_parent.num_links() +
                       node.left_child.right_parent.num_links())
                ret += math.log(recombination_rate)
                assert gap >= 1
                if gap == 1:
                    self.num_ancestral_recomb += 1
                else:
                    self.num_nonancestral_recomb += 1
                number_of_links -= gap
                number_of_lineages += 1
                if kuhner:# add rec
                    self.rec[node.index] = node.index
                    self.rec[ordered_nodes[counter+1].index] = ordered_nodes[counter+1].index
                counter += 2
                time = node.time
                rec_count += 1
            elif node.left_child.first_segment != None and\
                        node.right_child.first_segment != None:
                ret -= rate * (node.time - time)
                ret -=  math.log(2*Ne)
                if node.first_segment == None:
                    node_numlink = 0
                    number_of_lineages -= 2
                    counter += 1
                    if new_roots:
                        roots[node.index] = node.index
                else:
                    node_numlink = node.num_links()
                    number_of_lineages -= 1
                    counter += 1
                lchild_numlink = node.left_child.num_links()
                rchild_numlink = node.right_child.num_links()
                number_of_links -= (lchild_numlink + rchild_numlink) - node_numlink
                time = node.time
                coal_count += 1
                if new_roots:
                    new_coal[node.index] = node.index
            else:
                counter += 1
            if not NAM:
                assert node.left_child.first_segment != None
                assert node.right_child.first_segment != None
        if new_roots:
            return ret, roots, new_coal
        else:
            return ret

    def dump(self, path = ' ', file_name = 'arg.arg'):
        output = path + "/" + file_name
        with open(output, "wb") as file:
            pickle.dump(self, file)

    def load(self, path = ' '):
        with open(path, "rb") as file:
            return pickle.load(file)

    def print_state(self):
        print("self.arg.coal", self.coal)
        print("self.arg.rec", self.rec)
        print("self.arg.roots", self.roots)
        print("node", "time", "left", "right", "l_chi", "r_chi", "l_par", "r_par",
              "l_bp", "snps", "fir_seg_sam",
              sep="\t")
        for j in self.nodes:
            node = self.__getitem__(j)
            if node.left_parent is not None or node.left_child is not None:
                s = node.first_segment
                if s is None:
                    print(j, "%.5f" % node.time, "root", "root",
                              node.left_child.index,
                              node.right_child.index,
                              node.left_parent,node.right_parent,
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
                        print( j, "%.5f" % node.time, l, r,
                             node.left_child.index, node.right_child.index,
                              node.left_parent.index, node.right_parent.index,
                              node.breakpoint,
                              node.snps, s.samples, sep="\t")
                    s = s.next