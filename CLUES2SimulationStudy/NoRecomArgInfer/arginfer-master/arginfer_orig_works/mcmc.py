''' This module is responsible for the mcmc '''
import treeSequence
import copy
from tqdm import tqdm
import sys
import os
import shutil
from initialARG import *
from plots import *
# for recursion issue
# print("current recursion limit is ", sys.getrecursionlimit())
sys.setrecursionlimit(500000)

class TransProb(object):
    '''transition probability calculation'''
    def __init__(self):
        self.log_prob_forward = 0
        self.log_prob_reverse = 0

    def spr_choose_detach(self, numCoals, forward = True):
        '''
        choose a coal parent randomly and
        then choose one child with half prob
        '''
        if forward:
            self.log_prob_forward += math.log((1/numCoals) *(1/2))
        else: #reverse prob
            self.log_prob_reverse += math.log((1/numCoals) *(1/2))

    def spr_choose_reattach(self, numReattaches, forward = True):
        '''choose a reattach node among all possible nodes: numReattaches'''
        if forward:
            self.log_prob_forward += math.log(1/numReattaches)
        else: #reverse
            self.log_prob_reverse += math.log(1/numReattaches)

    def spr_reattach_time(self, new_time,lower_bound = 0, upper_bound = 0,
                          reattach_root = True, forward = True, lambd = 10000):
        '''if reattach_root: new time  is lower_bound + exp(1)
        else: new_time is  lower_bound  + U(lower_bound, upper_bound)
        '''
        if forward:
            if reattach_root: #expo
                self.log_prob_forward += math.log(lambd) - \
                                         (lambd * (new_time- lower_bound))
            else: #uniform
                self.log_prob_forward += math.log(1/(upper_bound -lower_bound))
                # self.log_prob_forward +=  math.log(lambd) - (lambd * (new_time - lower_bound))-\
                #            math.log(1- math.exp(- lambd *(upper_bound-lower_bound)))
        else:
            if reattach_root:
                self.log_prob_reverse += math.log(lambd)- (lambd * (new_time- lower_bound))
            else:
                self.log_prob_reverse += math.log(1/(upper_bound -lower_bound))
                # self.log_prob_reverse += math.log(lambd) - (lambd * (new_time - lower_bound))-\
                #            math.log(1- math.exp(- lambd *(upper_bound-lower_bound)))

    def spr_recomb_simulate(self, l_break, r_break, forward = True):
        '''simulate the recomb breakpoint if a window is given
        there are four scenarios on the recombination:
        1. ancestral to non-ancestral: no forward, yes reverse
        2. ancestral to ancestral: no forward, no reverse
        3. non ancestral to ancestral: yes forward, no reverse
        4. non ancestral to non ancestral: yes forward, yes reverse
        If a rec parent is NAM or rec child is NAM: No transition +
            no chnage of breakpoints
        :param l_break: left_breakpoint
        :param r_break: right_breakpoint
        :param forward: Forward transition if True, else, reverse.
        '''
        if forward:
            self.log_prob_forward += math.log(1/(r_break - l_break +1))
        else:
            self.log_prob_reverse += math.log(1/(r_break - l_break +1))

    def rem_choose_remParent(self, numRecPs, forward = True):
        '''choose a rec parent to remove'''
        if forward:
            self.log_prob_forward += math.log(1/numRecPs)
        else:
            self.log_prob_reverse += math.log(1/numRecPs)

    def add_choose_node(self, numnodes, forward = True):
        ''' prob of choosing a node to which add a rec'''
        if forward:
            self.log_prob_forward += math.log(1/numnodes)
        else:
            self.log_prob_reverse += math.log(1/numnodes)

    def add_choose_breakpoint(self, numlinks, forward = True):
        '''choose a genomic position as recomb breakpoint'''
        assert numlinks > 0
        if forward:
            self.log_prob_forward += math.log(1/numlinks)
        else:
            self.log_prob_reverse += math.log(1/numlinks)

    def add_choose_node_to_float(self, forward = True):
        '''randomly choose one of new parents to
        follow the child path, the other flow'''
        if forward:
            self.log_prob_forward += math.log(1/2)
        else:
            self.log_prob_reverse += math.log(1/2)

    def kuhner_num_nodes(self, num_nodes, forward =  True):
        '''resiprocal of the number of nodes in the ARG excluding the roots'''
        if forward:
            self.log_prob_forward = math.log(1/num_nodes)
        else:
            self.log_prob_reverse = math.log(1/num_nodes)

class MCMC(object):

    def __init__(self, ts_full = None, sample_size = 5, Ne =5000, seq_length= 3e5, mutation_rate=1e-8,
                 recombination_rate=1e-8,
                 input_data_path = '',
                 haplotype_data_name = None,
                 ancAllele_data_name='',
                 snpPos_data_name='',
                 outpath = os.getcwd()+"/output", verbose=False):
        '''
        :param ts_full: TS with full_record=True, if given, simulation data else: real data in
        :param sample_size:
        :param Ne:
        :param seq_length:
        :param mutation_rate:
        :param recombination_rate:
        :param input_data_path: the imput path to the data
        :param haplotype_data_name: a \t separater '' delimiter txt file
            :[T G C T...
              C G T A], nrow= num seqs, ncol= num od SNPs
        :param ancAllele_data_name: a '' delimiter txt file of ancestral allele for each site [G, T, A,...] (1*m)
        :param snpPos_data_name: a '' delimiter txt file of SNP positions on chromosome [1, 2000, ...] (1*m)
        :param outpath:
        :param verbose:
        '''
        self.ts_full = ts_full
        self.data = {} #a dict, key: snp_position- values: seqs with derived allele
        self.input_data_path = input_data_path
        self.haplotype_data_name = haplotype_data_name
        self.ancAllele_data_name = ancAllele_data_name
        self.snpPos_data_name = snpPos_data_name
        self.arg = ARG()
        self.mu = mutation_rate #mu
        self.r = recombination_rate
        self.n = sample_size
        self.Ne = Ne
        self.seq_length = seq_length
        self.m = 0 # number of snps
        self.log_lk = 0
        self.log_prior = 0
        self.outpath = outpath
        self.verbose = verbose
        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath)
        else:# if exists, first delete it and then create the directory
            shutil.rmtree(self.outpath)
            os.mkdir(self.outpath)
        self.transition_prob = TransProb()
        self.floatings = bintrees.AVLTree()#key: time, value: node index, d in the document
        self.floatings_to_ckeck = bintrees.AVLTree()#key index, value: index; this is for checking new roots
        self.new_names = bintrees.AVLTree()# index:index, of the new names (new parent indexes)
        self.lambd = 1/(self.Ne) # lambd in expovariate
        self.NAM_recParent = bintrees.AVLTree() # rec parent with seg = None
        self.NAM_coalParent = bintrees.AVLTree() # coal parent with seg = None
        self.coal_to_cleanup = bintrees.AVLTree()# ca parent with atleast a child.seg = None
        self.accept = False
        #------- To start get msprime output as initial
        self.get_initial_arg()
        self.arg.nextname = max(self.arg.nodes) + 1
        self.summary = pd.DataFrame(columns=('likelihood', 'prior', "posterior",
                                             'ancestral recomb', 'non ancestral recomb',
                                                'branch length',"mu", "r", "Ne", 'setup'))
        #---- kuhner
        self.floats = bintrees.AVLTree()#floatings: index:index
        self.partial_floatings = collections.defaultdict(list)# index:[a, b, num]
        self.need_to_visit = bintrees.AVLTree()#index:index, all the nodes we need to visit
        self.higher_times = collections.defaultdict(set)#time: (index)
        self.active_nodes = bintrees.AVLTree()#node.index
        self.active_links = 0
        self.original_interval = collections.defaultdict(list)# node.index:[left, right]
        #--- test
        self.detail_acceptance = collections.defaultdict(list)# transiton:[total, accepted]
        #---------- current ARG
        self.original_ARG = copy.deepcopy(self.arg) #  self.arg.copy()
        # for sd in update parameters
        self.default_mutation_rate = mutation_rate
        self.default_recombination_rate = recombination_rate
        self.default_Ne = Ne

    def read_convert_data(self):
        '''
        convert the imput data for the required data
        :param: general_path to the data
        :param haplotype_name: the haplotype name, important, this file should be a \t sep txt file
            with no header and n*m of the nuceotides
        :param ancestral_allele_name: the ancestral name. the file in this
            path should be a txt file with '' delimiter
        :param snp_pos_name: ANP chromosome posiitions. a txt file 1*m,  '' delimiter.
        '''
        haplotype_df = pd.read_csv(self.input_data_path + '/'+
                                   self.haplotype_data_name, sep="\t",header=None)
        if self.verbose:
            print("######succesfuly read haplotype file ")
        haplotype_df.columns=list(range(haplotype_df.shape[1]))
        ancestral_allele_np= np.loadtxt(self.input_data_path+'/'+
                                        self.ancAllele_data_name, dtype='str')
        if self.verbose:
            print("######succesfuly read ancestral_allele file")
        snp_pos_np= np.loadtxt(self.input_data_path+'/'+self.snpPos_data_name)
        if self.verbose:
            print("######succesfuly read SNP positions file ")
        new_data = {}
        for ind  in range(haplotype_df.shape[1]):
            derived_df= haplotype_df[ind]!= ancestral_allele_np[ind]
            seq_carry_derived = derived_df.index[derived_df].tolist()
            position = int(math.ceil(snp_pos_np[ind]))
            b = bintrees.AVLTree()
            b.update({k: k for k in seq_carry_derived})
            new_data[position] = b
        for key in new_data.keys():
            assert 0<len(new_data[key])<self.n
        self.data = new_data

    def get_initial_arg(self):
        '''
        TODO: build an ARG for the given data.
        '''
        if self.haplotype_data_name != None:
            #real data
            self.read_convert_data()
        else: # 'test' is for test in tests/
            ts_full= self.ts_full
            tsarg = treeSequence.TreeSeq(ts_full)
            tsarg.ts_to_argnode()
            self.data = treeSequence.get_arg_genotype(ts_full)
            #----- true values
            self.arg = tsarg.arg
            if  self.verbose:
                print("true number of rec", len(self.arg.rec)/2)
            self.log_lk = self.arg.log_likelihood(self.mu, self.data)
            self.log_prior = self.arg.log_prior(self.n, self.seq_length,
                                                self.r, self.Ne, False)
            np.save(self.outpath+"/true_values.npy", [self.log_lk, self.log_prior,
                                                  self.log_lk + self.log_prior,
                                                 self.arg.branch_length,
                                                  self.arg.num_ancestral_recomb,
                                                 self.arg.num_nonancestral_recomb,
                                                  self.mu, self.r, self.Ne])

        #-----initial
        init= Initial(self.data, self.n, self.seq_length,
                     self.Ne, self.mu, self.r)
        init.build()
        self.arg = init.arg
        if self.verbose:
            print("initial num rec:", len(self.arg.rec)/2)
        self.m = len(self.data)
        self.log_lk = self.arg.log_likelihood(self.mu, self.data)
        self.log_prior = self.arg.log_prior(self.n, self.seq_length,
                                            self.r, self.Ne, False)

    def truncated_expo(self, a, b, lambd):
        '''
        generate a random number from trucated exponential with rate lambd  in (a, b)
        '''
        assert b > a
        u= random.random()
        trunc_number = -(1/lambd) * (math.log(math.exp(-lambd*a) -\
                                              (u*(math.exp(-lambd*a)- math.exp(-lambd*b)))))
        if not (trunc_number < b) and not (a<trunc_number):
            return self.truncated_expo(a = a, b = b, lambd = lambd)
        else:
            return trunc_number

    def Metropolis_Hastings(self, new_log_lk, new_log_prior,
                            trans_prob = True, kuhner = False):
        '''if trans_prob: the ratio includes likelihood,
            prior and transiton probabilities,
            Otherwaise: only likelihood. This is
            for cases where the transition is based on the prior
        '''
        self.accept = False
        if trans_prob:
            ratio = new_log_lk + new_log_prior + \
                    self.transition_prob.log_prob_reverse - \
                    (self.log_lk + self.log_prior +
                     self.transition_prob.log_prob_forward)
        else:
            ratio = new_log_lk - self.log_lk
            if kuhner:
                ratio += (self.transition_prob.log_prob_reverse -\
                          self.transition_prob.log_prob_forward)
        if  self.verbose:
            print("forward_prob:", self.transition_prob.log_prob_forward)
            print("reverse_prob:", self.transition_prob.log_prob_reverse)
            print("ratio:", ratio)
            print("new_log_lk", new_log_lk, "new_log_prior", new_log_prior)
            print("old.log_lk", self.log_lk,"old.log_prior", self.log_prior)
        if math.log(random.random()) <= ratio: # accept
            self.log_lk = new_log_lk
            self.log_prior = new_log_prior
            self.accept = True

    def find_break_seg(self, z, break_point):
        '''z is the chain of segments
        return: the segment that includes breakpoint
            or the fist segment after breakpoint
            if None, then  second parent is empty
            '''
        while z.prev is not None:
            z = z.prev
        while z is not None:
            if z.contains(break_point):
                return z
            elif z.prev is not None and z.prev.right <= break_point and\
                    break_point < z.left:
                return z
            z = z.next
        return None

    def real_parent(self, node):
        ''' to be used in spr_reattachment_nodes
        return the original parent of a node
        if exist: return it otherwise None
        also this assumes that all rec parents
         are original nodes
        '''
        original_parent = None
        while node.left_parent != None:
            if self.new_names.__contains__(node.index):
                node = node.left_parent
            else:
                original_parent = node
                break
        return original_parent

    def spr_reattachment_nodes(self, detach_time, forward = True):
        '''
        return all possible reattachment nodes
        for detach
        return those exist at or after detach_time  and their
        segment != None,
        since we already update ARG, then if a segment is None,
        we should reattach to them
        :param forward: if false, for clean up step, when we check the new roots to calc the
            reverse prob, we need to also include those original nodes that will be floating
        '''
        reattach_nodes = bintrees.AVLTree()
        if forward:
            for node in self.arg.nodes.values():
                if node.time > detach_time:
                    if node.first_segment != None: # make sure it is not a NAM
                        reattach_nodes[node.index] = node.index
                elif node.left_parent != None and node.left_parent.time > detach_time and\
                        node.first_segment != None:
                    reattach_nodes[node.index] = node.index
                elif not self.floatings.is_empty(): #BUG2 in notes
                    reattach_nodes.update({v:v for v in self.floatings.values()})
        else:
            for node in self.arg.nodes.values():
                if node.time > detach_time:
                    if node.first_segment != None:
                        reattach_nodes[node.index] = node.index
                    elif not self.new_names.__contains__(node.index):
                        # that is a original node that will float in reverse
                        reattach_nodes[node.index] =  node.index
                elif node.left_parent != None:
                    if node.left_parent.time > detach_time:
                        if node.first_segment != None:
                            reattach_nodes[node.index] = node.index
                        elif not self.new_names.__contains__(node.index):
                            reattach_nodes[node.index] = node.index
                    else:# BUG5 NOTES
                        #its parent might not be its original
                        # (a new name with lower time than its original)
                        original_parent = self.real_parent(node.left_parent)
                        if original_parent != None and original_parent.time > detach_time:
                            if node.first_segment != None:
                                reattach_nodes[node.index] = node.index
                            elif not self.new_names.__contains__(node.index):
                                reattach_nodes[node.index] = node.index
        return reattach_nodes

    def incompatibility_check(self,node,  S1, S2, s, detach_snps, completed_snps):
        '''
        if the coalescence of child 1 and child 2 compatible for this snp.
        All are AVLTrees()
        S1: samples for first child segment
        S2: samples for the second child segment
        s:  the focal SNP
        node: the parent node
        '''
        ret = True
        D = self.data[s]
        # symmetric difference between S1 union S2  and D
        A = S1.union(S2)
        symA_D = A.difference(D)
        if len(symA_D) == 0:# subset or equal
            if len(A) == len(D): #put the mutation on this node
                node.snps.__setitem__(s, s)
                # delete s from S_F
                detach_snps.discard(s)
                # add to completed_snps
                completed_snps[s] = s
        elif len(symA_D) == len(A): # distinct
            pass
        else:#
            symD_A = D.difference(A)
            if len(symD_A) > 0: # incompatible
                ret = False
        return ret, detach_snps, completed_snps

    def update_ancestral_material(self, node_index, nodes_to_update,
                                  nodesToUpdateTimes, rem = None):
        '''
        update the materials of the parent(s) of node
        do not change any parent/child node, only update materials
        '''
        s = nodes_to_update.pop(node_index)
        node = self.arg.__getitem__(node_index)
        if rem != None and node.index == rem.index:#BUG7
            s = self.arg.copy_node_segments(node)
        if node.left_parent is None:
            raise ValueError("The root node shouldn't be added "
                             "to node_to_check at the first place")
        elif node.left_parent.index == node.right_parent.index:
            # common ancestor event
            node_sib = node.sibling()
            if nodes_to_update.__contains__(node_sib.index):
                sib_segs = nodes_to_update.pop(node_sib.index)
                if rem != None and node_sib.index == rem.index:
                    sib_segs = self.arg.copy_node_segments(node_sib)
            else:# copy them to remain intact
                sib_segs = self.arg.copy_node_segments(node_sib)
            if s is not None and sib_segs is not None:
                x = s.get_first_segment()
                y = sib_segs.get_first_segment()
                assert x is not None
                assert y is not None
                z = None
                defrag_required = False
                while x is not None or y is not None:
                    alpha = None
                    if x is None or y is None:
                        if x is not None:
                            alpha = x
                            x = None
                        if y is not None:
                            alpha = y
                            y = None
                    else:
                        if y.left < x.left:
                            beta = x
                            x = y
                            y = beta
                        if x.right <= y.left:
                            alpha = x
                            x = x.next
                            alpha.next = None
                        elif x.left != y.left:
                            alpha = self.arg.alloc_segment(x.left, y.left,
                                                    node.left_parent, x.samples)
                            x.left = y.left
                        else:
                            left = x.left
                            r_max = min(x.right, y.right)
                            right = r_max
                            alpha = self.arg.alloc_segment(left, right,
                                                    node.left_parent, x.union_samples(y))
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
                        if z is not None:
                            defrag_required |= z.right == alpha.left
                            z.next = alpha
                        alpha.prev = z
                        z = alpha
                if defrag_required:
                    z.defrag_segment_chain()
                if node.left_parent.left_parent is not None:# if not Root
                        nodes_to_update[node.left_parent.index] = z
                        nodesToUpdateTimes[node.left_parent.left_parent.time].add(node.left_parent.index)
                if z is not None:
                    z = z.get_first_segment()
                    #--- this is where we should add to floatings
                    if node.left_parent.left_parent is None:
                        # this is floating
                        self.floatings[node.left_parent.time] = node.left_parent.index
                        self.floatings_to_ckeck[node.left_parent.index] = node.left_parent.index
                    node.left_parent.first_segment = z
                    self.arg.store_node(z, node.left_parent)
                    self.NAM_coalParent.discard(node.left_parent.index)
                else:
                    node.left_parent.first_segment = None
                    self.NAM_coalParent[node.left_parent.index] = node.left_parent.index
                self.coal_to_cleanup.discard(node.left_parent.index)
            elif s is None and sib_segs is None:
                node.left_parent.first_segment = None
                if node.left_parent.left_parent is not None:
                    nodes_to_update[node.left_parent.index] = None
                    nodesToUpdateTimes[node.left_parent.left_parent.time].add(node.left_parent.index)
                self.NAM_coalParent[node.left_parent.index] = node.left_parent.index
                self.coal_to_cleanup[node.left_parent.index] = node.left_parent.index
            else: # s is  None or  sib_seg is None
                if sib_segs is None:
                    z = s.get_first_segment()
                else:# s is None
                    z = sib_segs.get_first_segment()
                if node.left_parent.left_parent is not None:
                    nodes_to_update[node.left_parent.index] = z
                    nodesToUpdateTimes[node.left_parent.left_parent.time].add(node.left_parent.index)
                else:
                    self.floatings[node.left_parent.time] = node.left_parent.index
                    self.floatings_to_ckeck[node.left_parent.index] = node.left_parent.index
                node.left_parent.first_segment = z
                self.arg.store_node(z, node.left_parent)
                self.NAM_coalParent.discard(node.left_parent.index)
                self.coal_to_cleanup[node.left_parent.index] = node.left_parent.index
        else:
            if s is None: # both parents are None
                z = None
                lhs_tail = None
            else:
                y = self.find_break_seg(s, node.breakpoint)
                if y is not None: # parent2 is not empty
                    x = y.prev
                    if y.left < node.breakpoint < y.right:#y.contains(node.breakpoint):# new is ancestral
                        # no forward + no reverse
                        z = self.arg.alloc_segment(node.breakpoint, y.right,
                                                   node, y.samples, None, y.next)
                        assert node.breakpoint < y.right
                        if y.next is not None:
                            y.next.prev = z
                        y.next = None
                        y.right = node.breakpoint
                        assert y.left < node.breakpoint
                        lhs_tail = y
                    elif x is not None:
                        # no forward+ yes reverse
                        # assert x.right is not y.left
                        assert x.right <= node.breakpoint <= y.left
                        x.next = None
                        y.prev = None
                        z = y
                        lhs_tail = x
                    else: # first parent is empty
                        # no update of breakpoint + no transition
                        z = y
                        lhs_tail = None
                else: # second parent is empty
                    # dont change the breakpoint no transition
                    z = None
                    lhs_tail = s
            if node.left_parent.left_parent is not None:# BUG1 in NOTES
                nodes_to_update[node.left_parent.index] = lhs_tail
                nodesToUpdateTimes[node.left_parent.left_parent.time].add(node.left_parent.index)
            else:
                #must be in self.floatings
                assert node.left_parent.index in [ind for ind in self.floatings.values()]
            if node.right_parent.left_parent is not None:# BUG1 in NOTES
                nodes_to_update[node.right_parent.index] = z
                nodesToUpdateTimes[node.right_parent.left_parent.time].add(node.right_parent.index)
            else:
                #must be in self.floatings
                assert node.right_parent.index in [ind for ind in self.floatings.values()]
            # TODO: Implement a way to check whether a recombination is removable, and
            # delete those recombination nodes which are not. Will have to be done once the
            # ARG has been constructed by checking which ones will violate snps
            # node = self.arg.alloc_node(self.arg.new_name(), time, lhs_tail.node, lhs_tail.node)
            # lhs_tail.node.left_parent = node
            self.arg.store_node(lhs_tail, node.left_parent)
            self.arg.store_node(z, node.right_parent)
            # if NAM put them in NAM_recParent
            if lhs_tail is None:
                self.NAM_recParent[node.left_parent.index] = node.left_parent.index
            else:
                self.NAM_recParent.discard(node.left_parent.index)
            if z is None:
                self.NAM_recParent[node.right_parent.index] = node.right_parent.index
            else:
                self.NAM_recParent.discard(node.right_parent.index)
        return nodes_to_update, nodesToUpdateTimes

    def get_detach_SF(self, detach, sub_interval=[]):
        '''get snps from detach
         that we need to check for incompatibility
          1. for a seg find all snps within.
          2. for each snp check if the mut has happened or will happen in future
          3. if in future add to detach_snps
          :param sub_interval: only check those in [start, end) interval
                    this is currenctly using in adjust_breakpoint algorithm
          '''
        detach_snps = bintrees.AVLTree()
        seg = detach.first_segment
        while seg is not None:
            seg_snps =[]
            if sub_interval:# adjust breakpoint
                intersect = seg.get_intersect(sub_interval[0], sub_interval[1])
                if intersect:
                    seg_snps = [key for key in self.data.keys() if intersect[0] <= key < intersect[1]]
            else:
                seg_snps = [key for key in self.data.keys() if seg.left <= key < seg.right]
            for item in seg_snps:
                D = self.data[item]
                A = seg.samples
                symA_D = A.difference(D)
                if len(symA_D) == 0: # A is equal or subset
                    if len(A) == len(D):
                        pass # its on detach or already occured
                        # assert detach.snps.__contains__(item)
                    else:# A is subset+ mut in future
                        detach_snps[item] = item
                elif len(symA_D) == len(A):# distinct
                    # A is all ancestral allel + might in future
                    # if self.n - len(D) < len(A): # if not all ancestral are here
                    detach_snps[item] = item
                elif len(D.difference(A)) > 0:# len(A.symD) >0 and len(D.symA)>0
                    raise ValueError("The original ARG was incompatible")
                else:# len(D.symmetric_difference(A)) == 0#mutation already happened
                    pass
            seg = seg.next
        return detach_snps

    def update_all_ancestral_material(self, node, rem = None):
        '''update the ancestral materials of all the ancestors of node
        :param node: a list of nodes [node]
        forward the node entry if of length 2. the first is child and
         the second is the remParentSib. The later segments may change during update,
         so we cant copy its segments at first as it might not be the true one.'''
        nodes_to_update = {}
        their_times = collections.defaultdict(set) #key:time, value (set of indexes)
        for n in node:
            if n.left_parent is not None:
                #BUG7
                nodes_to_update[n.index] = self.arg.copy_node_segments(n)
                their_times[n.left_parent.time].add(n.index)
        if rem != None:# only one case remove forward node =[child, remParentSib]
            nodes_to_update[rem.index] = None
        current_time = min(their_times)
        while nodes_to_update:
            min_time = min(their_times)
            if min_time < current_time:
                print("current_time", current_time)
                print("min_time", min_time)
                print("nodes to update", nodes_to_update)
                print("theit time", their_times)
                self.print_state()
            assert min_time >= current_time
            next_ind = their_times.pop(min(their_times))
            while next_ind:
                next_index = next_ind.pop()
                if nodes_to_update.__contains__(next_index):# BUG4_2
                    if self.arg.__contains__(next_index):# might pruned already in backtrack
                        nodes_to_update, nodes_times = self.update_ancestral_material(next_index,
                                                        nodes_to_update, their_times, rem)
                    else:
                        del nodes_to_update[next_index]
            current_time = min_time

    def find_sc_original_parent(self, node):
        '''find second child of original_parent the original parent'''
        node = node.left_parent # this is original parent
        sc_original_parent = None; valid = True
        while node.left_parent is not None:
            if node.left_parent.index != node.right_parent.index:
                if node.left_parent.first_segment == None or\
                    node.right_parent.first_segment == None:
                    valid = False
                    break
                assert not self.new_names.__contains__(node.left_parent.index)
                sc_original_parent = node.left_parent
                break
            elif self.new_names.__contains__(node.left_parent.index):
                node = node.left_parent
            else:
                sc_original_parent = node.left_parent
                break
        return sc_original_parent, valid

    def find_original_child(self, node):
        '''find the original child for node, the child in the original ARG
        '''
        original_child = None
        while node.left_child != None:
            leftch_ind =node.left_child.index
            rightch_ind = node.right_child.index
            if leftch_ind == rightch_ind:
                if not self.floatings_to_ckeck.__contains__(leftch_ind):
                    original_child = node
                    break
                else:
                    break
            else:
                if not self.floatings_to_ckeck.__contains__(leftch_ind) and \
                    not self.floatings_to_ckeck.__contains__(rightch_ind):
                    original_child = node
                    break
                elif self.floatings_to_ckeck.__contains__(leftch_ind) and \
                    self.floatings_to_ckeck.__contains__(rightch_ind):
                    break
                elif self.floatings_to_ckeck.__contains__(leftch_ind):
                    node = node.right_child
                else:
                    node = node.left_child
        return original_child

    def generate_new_time(self, lower_bound = 0, upper_bound=1, root= True):
        '''generate a new time'''
        if root:
            return lower_bound + random.expovariate(self.lambd)
        else:
            #from uniform
            return  random.uniform(lower_bound, upper_bound)
            # from truncate exponential
            # return self.truncated_expo(lower_bound, upper_bound, self.lambd)

    def spr_reattach_floatings(self, detach, sib, old_merger_time):
        '''reattach all the floatings including detach'''
        while self.floatings:
            min_time = self.floatings.min_key()
            node = self.arg.nodes[self.floatings[min_time]]
            self.floatings.discard(min_time)
            # check if it is still floating
            still_floating = False
            if node.index == detach.index:
                still_floating = True
            elif node.first_segment is not None and node.left_parent is None:
                still_floating = True
            if still_floating:
                # first check if the move is valid ---> do this in cleanup
                # if  node.left_child.first_segment is None or \
                #     node.right_child.first_segment is None:
                #     # this move is not valid
                #     pass
                # --- reattach
                # 1. find potential reattachment
                all_reattachment_nodes = self.spr_reattachment_nodes(min_time)
                if node.left_parent is not None: # only detach
                    assert node.index == detach.index
                    all_reattachment_nodes.discard(node.left_parent.index)
                    #C is not added to all_reattach_nodes: because t_c<t_F
                    # and sib.left_parent = None ---after detach if P is root
                    all_reattachment_nodes[sib.index] = sib.index
                all_reattachment_nodes.discard(node.index)
                if all_reattachment_nodes.is_empty(): # BUG2 in notes, the grand root
                    print("node", node.index)
                    print("self.floatings", self.floatings)
                    raise ValueError("all_reattachment_nodes is empty."
                                     " There must be atleast one node to reattach.")
                reattach = self.arg.__getitem__(random.choice(list(all_reattachment_nodes)))
                if  self.verbose:
                    print("node", node.index,"with time", node.time, "rejoins to ",
                      reattach.index, "with time", reattach.time)

                #---trans_prob for choose reattach
                self.transition_prob.spr_choose_reattach(len(all_reattachment_nodes))
                max_time = max(min_time, reattach.time)
                if node.index == detach.index and reattach.index == sib.index:
                    self.floatings.discard(old_merger_time) # this is for removing sib
                if reattach.left_parent is None:
                    new_merger_time = self.generate_new_time(max_time)
                    if  self.verbose:
                        print("new time is from an exponential distribution"
                              " + reattach.time:", reattach.time)
                        print("new_merger_time", new_merger_time)
                    self.transition_prob.spr_reattach_time(new_merger_time, max_time, 0,
                                                           True, True, self.lambd)
                else:
                    new_merger_time = self.generate_new_time(max_time, reattach.left_parent.time, False)
                    if  self.verbose:
                        print("new_merger_time", new_merger_time, "reattach.left_parent.time",
                          reattach.left_parent.time)
                    if max_time >= reattach.left_parent.time:
                        self.print_state()
                        print("all_reattachment_nodes\n", all_reattachment_nodes)
                        print("node_time", node.time, "reattach_time", reattach.time,
                              "reattach_left_parent_time", reattach.left_parent.time)
                        print("node", node.index, "reattach", reattach.index,
                              "reattach.left_parent.index", reattach.left_parent.index)
                        print("new_merger_time", new_merger_time)
                    self.transition_prob.spr_reattach_time(new_merger_time, max_time,
                                                   reattach.left_parent.time, False, True, self.lambd)
                #-- reattach
                self.new_names = self.arg.reattach(node, reattach, new_merger_time, self.new_names)
                #---- update
                self.update_all_ancestral_material([node])
                #---
                self.floatings.discard(reattach.time) #if any

    def clean_up(self, coal_to_cleanup):
        '''clean up the Accepted ARG. the ancestral material and snps
        has already set up. This method is only for cleaning the ARG from NAM lineages
        NAM is No ancestral Material nodes that are not a root.
        order the nodes by time and then clean up'''
        def reconnect(child, node):# BUG7
            '''from child--> node--> parent: TO child ---> parent '''
            leftparent = node.left_parent
            rightparent = node.right_parent
            child.left_parent = leftparent
            child.right_parent = rightparent
            child.breakpoint = node.breakpoint
            leftparent.update_child(node, child)
            rightparent.update_child(node, child)

        while coal_to_cleanup:
            node = self.arg.__getitem__(coal_to_cleanup.pop(coal_to_cleanup.min_key()))
            if node.left_child == None and node.right_child == None:
                if node.left_parent is not None:
                    assert node.left_parent.index == node.right_parent.index
                    node.left_parent.update_child(node, None)
            elif node.left_child != None and node.right_child is None:
                if node.left_parent is not None:
                    reconnect(node.left_child, node)
            elif node.right_child != None and node.left_child is None:
                if node.left_parent is not None:
                    reconnect(node.right_child, node)
            else: # both not None
                if node.left_child.first_segment == None and\
                        node.right_child.first_segment == None:
                    if node.left_parent is not None:
                        assert node.first_segment is  None
                        assert node.left_parent.index == node.right_parent.index
                        node.left_child.left_parent = None
                        node.left_child.right_parent = None
                        node.right_child.left_parent = None
                        node.right_child.right_parent = None
                        node.left_parent.update_child(node, None)
                elif node.left_child.first_segment != None and\
                        node.right_child.first_segment is None:
                    assert node.left_parent is not None
                    reconnect(node.left_child, node)
                elif node.right_child.first_segment != None and\
                        node.left_child.first_segment is None:
                    assert node.left_parent is not None
                    reconnect(node.right_child, node)
                else: # non None
                    raise ValueError("both children have seg, so "
                                     "this shouldn't be in coal_to_clean")
            del self.arg.nodes[node.index]

    def spr_validity_check(self, node,  clean_nodes,
                           detach_snps,  completed_snps, reverse_done):
        '''after rejoining all the floating, it is time to check if
        the changes cancels any recombination.
        Also, if there is a new root (node.segment == None and node.left_parent!=None),
        calculate the reverse prob for them.
        In addition, whether there is a incompatibility.
        This method is responsible for checking the above mentioned for node.
        :param node: the node we need to check its validability and revers prob (if applicable)
        :param detach: the detach node in spr
        :param clean_nodes: a dict of k: time, v: node.indexes the nodes we need to check for validity
        '''
        valid = True
        if node.first_segment != None:# not a new root
            assert node.left_parent != None
            # check incompatibility if both children have segments
            if node.left_child.index != node.right_child.index:
                # first delete all the detach_snps in node.snps
                node.snps = node.snps.difference(detach_snps)
                node.snps = node.snps.difference(completed_snps)
                #find all SNP on left_child
                lc_variants = node.left_child.get_variants(self.data)
                #find all SNP on right child
                rc_variants = node.right_child.get_variants(self.data)
                # find the intersection of them
                intersect_variants = detach_snps.intersection(lc_variants, rc_variants)
                for snp in intersect_variants:
                    S1 = node.left_child.x_segment(snp).samples
                    S2 = node.right_child.x_segment(snp).samples
                    valid, detach_snps, completed_snps = \
                        self.incompatibility_check(node, S1, S2, snp, detach_snps,
                                                   completed_snps)
                    if not valid:
                        break
            else:
                assert len(node.snps) == 0
            clean_nodes[node.left_parent.time].add(node.left_parent.index)
            clean_nodes[node.right_parent.time].add(node.right_parent.index)
            #add the parents to the clean_node
        elif node.left_parent != None:
            if node.left_parent.index != node.right_parent.index:
                valid = False
            elif node.left_child.index == node.right_child.index:
                valid = False
            else:
                # might be new root
                # is it really new root?
                # find origina_child: it might be a parent of two floating lins, then NAM in reverse.
                original_child = self.find_original_child(node)
                if original_child != None and\
                        not reverse_done.__contains__(node.index):# NOTE AFTER BUG5:
                    if not self.original_ARG.__contains__(original_child.index):
                        if self.original_ARG.__contains__(original_child.left_child.index):
                            original_child = original_child.left_child
                        elif self.original_ARG.__contains__(original_child.right_child.index):
                            original_child = original_child.right_child
                        else:
                            valid = False
                            # if one or both of original_child is new_name, then
                            # let say in original p1= 4, c1= 1. F=2 joins c1, new parent p2= 3.
                            # now parent of c1 is p2. if F2 rejoins F, new parent is p3. now (p3, c1) are
                            #sibling with parent p2. If p2 become root, none of its child are
                            # in floating, then we guess p2 is in original while it is not.
                            # if F3 joins c1 at a time before p1 with new paren p4,
                            # then p3 and p4 both are new name and child of p2. if p2 is root, again
                            # the function says it is in original but it is not, while this means c1
                            # can be floating from time of p2. finding this is difficult so I ignore.
                    original_parent = None
                    if valid:
                        if original_child.time == self.original_ARG.__getitem__(original_child.index).time:
                            original_parent = self.original_ARG.__getitem__(original_child.index).left_parent
                        else:# the original_child is detach_parent and its time is different
                            valid= False
                    if original_parent is not None:# original_child was not a root in G_j
                        if original_parent.left_child.index == original_parent.right_child.index:
                            # it is a rec parent
                            valid = False
                        else:
                            second_child = original_parent.left_child
                            if second_child.index == original_child.index:
                                second_child = original_parent.right_child
                            all_reattachment_nodes = self.spr_reattachment_nodes(node.time, False)
                            all_reattachment_nodes.discard(original_parent.index)
                            all_reattachment_nodes.discard(node.index)
                            # ---- reverse of choosin g a lineage to rejoin
                            self.transition_prob.spr_choose_reattach(len(all_reattachment_nodes), False)
                            #----- reverse prob time
                            # if node.index != detach.index: # already done for detach
                            original_grandparent = original_parent.left_parent
                            # assert node.index != second_child.index
                            if original_grandparent is None:
                                self.transition_prob.spr_reattach_time(original_parent.time,
                                        max(node.time, second_child.time), 0, True, False, self.lambd)
                            else:
                                if original_grandparent.time <= max(node.time, second_child.time):
                                    self.original_ARG.print_state()
                                    print("original_grandparent.time:", original_grandparent.time)
                                    print("node.time", node.time, "second_child.time", second_child.time)
                                    print("original_grandparent", original_grandparent.index)
                                    print("original_parent:", original_parent.index, "its time",
                                          original_parent.index)
                                    print("second_child", second_child.index, "original_child",
                                          original_child.index)
                                    print("node:", node.index)
                                self.transition_prob.spr_reattach_time(original_parent.time,
                                        max(node.time, second_child.time), original_grandparent.time
                                                                        ,False, False, self.lambd)
                            reverse_done[second_child.index] = second_child.index
                clean_nodes[node.left_parent.time].add(node.left_parent.index)
        return valid, clean_nodes, detach_snps, completed_snps, reverse_done

    def all_validity_check(self, clean_nodes, detach_snps):
        '''do spr_validity_check()for all the needed nodes
        '''
        valid = True # move is valid
        completed_snps = bintrees.AVLTree() # those of detach_snps that completed already
        reverse_done = bintrees.AVLTree()# reverse prob in done for them
        while valid and clean_nodes:
            # get the min_time one
            nodes = clean_nodes.pop(min(clean_nodes))
            assert 0 < len(nodes) <= 2
            if len(nodes) == 2:# two rec parents
                nodes = [self.arg.__getitem__(nodes.pop()), self.arg.__getitem__(nodes.pop())]
                assert nodes[0].left_child.index == nodes[0].right_child.index
                assert nodes[1].left_child.index == nodes[1].right_child.index
                assert nodes[0].left_child.index == nodes[1].left_child.index
                if nodes[0].first_segment is None or nodes[1].first_segment is None:
                    valid = False # cancels a rec
                    break
            else:
                assert len(nodes) == 1
                nodes = [self.arg.__getitem__(nodes.pop())]
            while nodes:
                node = nodes.pop(0)
                valid, clean_nodes, detach_snps, completed_snps, reverse_done = \
                    self.spr_validity_check(node, clean_nodes, detach_snps,
                                             completed_snps, reverse_done)
                if not valid:
                    break
        return valid

    #============
    # remove recombination

    def empty_containers(self):
        '''empty all the containers'''
        self.floatings.clear()
        self.NAM_recParent.clear()
        self.NAM_coalParent.clear()
        self.coal_to_cleanup.clear()
        self.new_names.clear()
        self.floatings_to_ckeck.clear()
        self.transition_prob.log_prob_forward = 0
        self.transition_prob.log_prob_reverse = 0
        self.arg.nextname = max(self.arg.nodes) + 1
        self.arg.get_available_names()
        self.original_ARG = ARG()

    def detach_otherParent(self, remParent, otherParent, child):
        '''child ---rec---(remPrent, otherParent), now that we've removed
        remParent, rejoin child to otherParent.left_parent and detach otherParent
        a. If remParent and otherParent coalesce back: parent of both rem and other is equal
            then we should rejoin child to  otherParent.left_parent.left_parent (must exist)
        b. Otherwise:
            rejoin child to otherParent.left_parent
        '''
        if remParent.left_parent.index == otherParent.left_parent.index:
            invisible = True # rec is invisible
            assert remParent.left_parent.left_parent is not None
            othergrandparent = otherParent.left_parent
            child.left_parent = othergrandparent.left_parent
            child.right_parent = othergrandparent.right_parent
            child.breakpoint = othergrandparent.breakpoint
            parent = othergrandparent.left_parent
            parent.update_child(othergrandparent, child)
            parent = othergrandparent.right_parent
            parent.update_child(othergrandparent, child)
        else:
            invisible = False# rec is not invisible
            child.left_parent = otherParent.left_parent
            child.right_parent = otherParent.right_parent
            child.breakpoint = otherParent.breakpoint
            parent = otherParent.left_parent
            parent.update_child(otherParent, child)
            parent = otherParent.right_parent
            parent.update_child(otherParent, child)
        return invisible

    #=============
    #  ADD recombination

    def check_root_parents(self, new_roots):
        '''after mh acceptance, make sure all the new roots dont have any parent
        This is not taken care of in clean_up(), because coal_to_clean does
        not contain nodes without seg where both children have segments (roots)
        '''
        for ind in new_roots:
            root = self.arg.__getitem__(ind)
            assert root.first_segment is None
            if root.left_parent is not None:
                assert not self.arg.nodes.__contains__(root.left_parent.index)
                root.left_parent = None
                root.right_parent = None
            else:
                root.right_parent = None

    def add_choose_child(self):
        '''choose a node to put a recombination on it'''
        ind = random.choice(list(self.arg.nodes.keys()))
        if not self.arg.roots.__contains__(ind):
            return self.arg.__getitem__(ind)
        else:
            return self.add_choose_child()

    def split_node(self, child, k, t):
        '''split a node (child) to two parent node from k at time t
        and add the parents to the arg
        '''
        s = self.arg.copy_node_segments(child) #child segments
        y = self.find_break_seg(s, k)
        if y is None:
            print("in split node", "child is", child.index, "k:", k,
                  "child_head", child.first_segment.left,
                  "child tail right", child.get_tail().right, "y", y)
        assert y is not None
        x = y.prev
        if y.left < k:
            assert k < y.right
            z = self.arg.alloc_segment(k, y.right, y.node,
                                       y.samples, None, y.next)
            if y.next is not None:
                y.next.prev = z
            y.next = None
            y.right = k
            lhs_tail = y
        else:
            assert x is not None
            x.next = None
            y.prev = None
            z = y
            lhs_tail = x
        leftparent = self.arg.alloc_node(self.arg.new_name(),
                                         t, lhs_tail.node, lhs_tail.node)
        lhs_tail.node.left_parent = leftparent
        self.arg.store_node(lhs_tail, leftparent)
        rightparent = self.arg.alloc_node(self.arg.new_name(),
                                          t, z.node, z.node)
        z.node.right_parent = rightparent
        self.arg.store_node(z, rightparent)
        #--- update breakpoint
        child.breakpoint = k
        #--- add to rec
        self.arg.rec[leftparent.index] = leftparent.index
        self.arg.rec[rightparent.index] = rightparent.index
        return leftparent, rightparent

    def add_recombination(self):
        '''
        Transition number 3
        add a recombination to the ARG
        1. randomly choose a node excluding roots
        2. randomly choose a time on the node
        3. randomly choose a breakpoint on the node and split the node to two
        4. randomly choose a parent to follow the path, the other to float
        5. update ancestral material
        6. do spr on the floating node to reattach it to the ARG
            a) randomly choose the potential node to rejoin the detach
            b) rejoin the floating and update ancestral material
            c) rejoin all the floatings
            d) check the validity and compatibility
            e) m-h
        forward transition: also reattach detachPrent, choose time to reattach
        '''
        # self.original_ARG = copy.deepcopy(self.arg)#copy.deepcopy(self.arg)
        self.arg.dump(path = self.outpath, file_name = 'arg.arg')
        self.original_ARG =  self.arg.load(path = self.outpath+'/arg.arg')
        assert self.arg.__len__() == (len(self.arg.coal) + len(self.arg.rec) + self.n)
        child = self.add_choose_child()
        assert child.first_segment != None and child.left_parent != None
        self.transition_prob.add_choose_node(len(self.arg.nodes) - len(self.arg.roots))#1
        head = child.first_segment
        tail = child.get_tail()
        if tail.right - head.left <= 1:# no link
            if  self.verbose:
                print("the node has less than 2 links")
            valid = False
        else:
            #breakpoint and time
            self.transition_prob.add_choose_breakpoint(tail.right - head.left - 1)#2
            break_point = random.choice(range(head.left + 1, tail.right))
            new_rec_time = self.generate_new_time(child.time, child.left_parent.time, False)
            self.transition_prob.spr_reattach_time(new_rec_time, child.time,
                                                   child.left_parent.time, False, True, self.lambd)#3
            oldleftparent = child.left_parent
            oldrightparent = child.right_parent
            oldbreakpoint = child.breakpoint
            newleftparent, newrightparent = self.split_node(child, break_point, new_rec_time)
            #--- choose one to follow child path
            followParent = random.choice([newleftparent, newrightparent])
            if followParent.index == newleftparent.index:
                detachParent = newrightparent
            else:
                detachParent = newleftparent
            self.transition_prob.add_choose_node_to_float()#4
            if  self.verbose:
                print("child", child.index, "followparent", followParent.index,
                      "detach", detachParent.index, "old_leftpanre",
                      oldleftparent.index, "oldrightparent", oldrightparent.index)
            #----
            followParent.breakpoint = oldbreakpoint
            followParent.left_parent = oldleftparent
            followParent.right_parent = oldrightparent
            oldleftparent.update_child(child, followParent)
            oldrightparent.update_child(child, followParent)
            #now update ancestral material
            if  self.verbose:
                print("child.left_parent", child.left_parent.index)
                print("childrp", child.right_parent.index)
            self.update_all_ancestral_material([followParent])
            #--- reattach detachParent
            self.floatings[detachParent.time] = detachParent.index
            self.floatings_to_ckeck[detachParent.index] = detachParent.index
            self.spr_reattach_floatings(child, child, child.time) # fake *args
            if self.NAM_recParent: # rec is canceled, reject
                if  self.verbose:
                    print("not valid due to removing rec")
                valid = False
            else:
                #validability, compatibility
                detach_snps = self.get_detach_SF(detachParent)
                clean_nodes = collections.defaultdict(set) #key: time, value: nodes
                clean_nodes[child.left_parent.time].add(child.left_parent.index)
                clean_nodes[child.right_parent.time].add(child.right_parent.index)
                valid = self.all_validity_check(clean_nodes, detach_snps)
            if valid:
                new_log_lk = self.arg.log_likelihood(self.mu, self.data)
                new_log_prior, new_roots, new_coals = self.arg.log_prior(self.n,
                                                self.seq_length, self.r, self.Ne, True, True)
                #--- reverse prob--
                self.transition_prob.rem_choose_remParent(len(self.original_ARG.rec)+2, False)
                self.Metropolis_Hastings(new_log_lk, new_log_prior)
                if  self.verbose:
                    print("mh_accept ", self.accept)
                if self.accept:
                    #update coal, and rec
                    self.clean_up(self.coal_to_cleanup)
                    self.check_root_parents(new_roots)
                    self.arg.roots = new_roots # update roots
                    self.arg.coal = new_coals
                    assert len(self.arg.rec) == len(self.original_ARG.rec)+2
                else:
                    self.arg= self.original_ARG
            else:
                self.arg = self.original_ARG
        self.empty_containers()

    #=========
    # adjust times move

    def adjust_times(self, calc_prior = True):
        '''
        Transition number 4
        modify the node times according to CwR
        also calculate the prior in place if calc_prior = true
        '''
        # self.original_ARG = copy.deepcopy(self.arg)
        ordered_nodes = [v for k, v in sorted(self.arg.nodes.items(),
                                     key = lambda item: item[1].time)]
        number_of_lineages = self.n
        number_of_links = number_of_lineages * (self.seq_length - 1)
        number_of_nodes = self.arg.__len__()
        counter = self.n
        prev_t  = 0
        log_prior = 0
        original_t = [0 for i in range(number_of_nodes)] # for reverting the move
        while counter < number_of_nodes:
            node = ordered_nodes[counter]
            rate = (number_of_lineages * (number_of_lineages - 1)
                    / (4*self.Ne)) + (number_of_links * (self.r))
            t = prev_t + random.expovariate(rate)
            # ret -= rate * (node.time - time)
            if node.left_child.index == node.right_child.index: #rec
                original_t[counter] = node.time
                original_t[counter +1] = node.time
                node.time = t
                ordered_nodes[counter+1].time = t
                gap = node.left_child.num_links()-\
                          (node.left_child.left_parent.num_links() +
                           node.left_child.right_parent.num_links())
                if calc_prior:
                    log_prior -= (rate * (t - prev_t))
                    log_prior += math.log(self.r)
                number_of_links -= gap
                number_of_lineages += 1
                counter += 2
            else: #CA
                original_t[counter] = node.time
                node.time = t
                if calc_prior:
                    log_prior -= (rate * (t - prev_t))
                    log_prior -=  math.log(2*self.Ne)
                if node.first_segment == None:
                    assert node.left_parent == None
                    node_numlink = 0
                    number_of_lineages -= 2
                    counter += 1
                else:
                    node_numlink = node.num_links()
                    number_of_lineages -= 1
                    counter += 1
                lchild_numlink = node.left_child.num_links()
                rchild_numlink = node.right_child.num_links()
                number_of_links -= (lchild_numlink + rchild_numlink) - node_numlink
            prev_t = t
        # m-h, without prior because it get cancels out with transition probs
        old_branch_length = self.arg.branch_length
        new_log_lk = self.arg.log_likelihood(self.mu, self.data)
        self.Metropolis_Hastings(new_log_lk, log_prior, trans_prob = False)
        if  self.verbose:
            print("mh_accept ", self.accept)
        if not self.accept:
            # self.arg = self.original_ARG
            self.arg.branch_length = old_branch_length
            self.revert_adjust_times(ordered_nodes, original_t)
        self.empty_containers()
        ordered_nodes=[]
        original_t=[]

    def revert_adjust_times(self, ordered_nodes, original_t):
        '''revert the proposed ARG by
        adjust_times move to its orgiginal
        '''
        counter = self.n
        number_of_nodes = len(ordered_nodes)
        while counter < number_of_nodes:
            node = ordered_nodes[counter]
            node.time = original_t[counter]
            if node.left_child.index == node.right_child.index: #rec
                ordered_nodes[counter + 1].time = original_t[counter + 1]
                counter += 2
            else:
                counter += 1

    #==============
    #transition 5: adjust recombination position

    def adjust_breakpoint(self):
        '''transition number 5
        change the  breakpoint of
        an existing recombination
        1. randomly choose a recombination event
        2. simulate a new breakpoint for the rec
        3. update ancestral material
        4. if floating: reattach
        5. compatibility/ validity check
        transition would only be for floatings
        '''
        # self.original_ARG = copy.deepcopy(self.arg)#copy.deepcopy(self.arg)
        self.arg.dump(path = self.outpath, file_name = 'arg.arg')
        self.original_ARG =  self.arg.load(path = self.outpath+'/arg.arg')
        assert not self.arg.rec.is_empty()
        recparent = self.arg.__getitem__(random.choice(list(self.arg.rec.keys())))
        assert recparent.left_child.index == recparent.right_child.index
        child = recparent.left_child
        old_breakpoint = child.breakpoint
        if  self.verbose:
            print("child is ", child.index)
        # TODO complete this
        leftparent = child.left_parent
        rightparent  = child.right_parent
        assert leftparent.index != rightparent.index
        # simulate a new breakpoint
        child_head = child.first_segment
        child_tail = child.get_tail()
        new_breakpoint = random.choice(range(child_head.left + 1, child_tail.right))
        if  self.verbose:
            print("old_breakpoint is",  child.breakpoint)
            print("new_breakpoint:", new_breakpoint)
        y = self.find_break_seg(child_head, child.breakpoint)
        if y.prev is not None  and self.verbose:
            print("x.right", y.prev.right)
        if  self.verbose:
            print("y left", y.left, "y.right", y.right)
        if new_breakpoint == old_breakpoint or\
                (not y.contains(old_breakpoint) and\
                 y.prev.right <= old_breakpoint<= y.left and\
                 y.prev.right <= new_breakpoint<= y.left):
            if  self.verbose:
                print("new_breakpoint is still non ancestral and at the same interval")
            # no prob changes, the breakpoint is still in the same non ancestral int
            child.breakpoint = new_breakpoint
            self.accept = True
        else:
            assert old_breakpoint != new_breakpoint
            start = old_breakpoint
            end = new_breakpoint
            if new_breakpoint < old_breakpoint:#
                start = new_breakpoint
                end = old_breakpoint
            # update ancestral material
            child.breakpoint = new_breakpoint
            self.update_all_ancestral_material([child])
            #reattach flaotings if any
            self.spr_reattach_floatings(child, child, child.time)#fake *args
            # is there any canceled rec?
            if self.NAM_recParent:
                if  self.verbose:
                    print("not valid due to removing rec")
                valid = False
            else: # check validity
                # get the affected snps from [start, end) interval
                child_snps = self.get_detach_SF(child, [start, end])
                clean_nodes = collections.defaultdict(set) #key: time, value: nodes
                clean_nodes[child.left_parent.time].add(child.left_parent.index)
                clean_nodes[child.right_parent.time].add(child.right_parent.index)
                valid = self.all_validity_check(clean_nodes, child_snps)
            if valid:
                new_log_lk = self.arg.log_likelihood(self.mu, self.data)
                new_log_prior, new_roots, new_coals = self.arg.log_prior(self.n,
                                                self.seq_length, self.r, self.Ne, True, True)
                #--- now mh
                self.Metropolis_Hastings(new_log_lk, new_log_prior)
                if  self.verbose:
                    print("mh_accept ", self.accept)
                if self.accept:
                    #update coal, and rec
                    self.clean_up(self.coal_to_cleanup)
                    self.check_root_parents(new_roots)
                    self.arg.roots = new_roots # update roots
                    self.arg.coal = new_coals
                else:#revert
                    self.arg = self.original_ARG
            else:
                self.arg = self.original_ARG
        self.empty_containers()

    #===================
    #transition 6 : Kuhner move

    def find_active_nodes(self, t):
        '''nodes higher than t and also active nodes immediately after t'''
        for node in self.arg.nodes.values():
            if node.time >t:
                self.higher_times[node.time].add(node.index)
            elif node.left_parent != None and node.left_parent.time > t:
                self.active_nodes[node.index] = node.index

    def update_prune_parent(self, prune, node, deleted_nodes, parents):
        '''if prune is a child of a CA event:
        easily detach it and remove its parent.
        if sibling is root, depending on its time add to float
        If prune is a child of a rec, delete both parents,
        and continue deleting them till reach a CA, then detach
        '''
        if node.left_parent is None:
            # is possible only if node is  parent of a rec,
            # so it should have already remove in else option
            assert not self.arg.__contains__(node.index)
            assert not self.floats.__contains__(node.index)
            assert not self.partial_floatings.__contains__(node.index)
            assert not self.active_nodes.__contains__(node.index)
        elif node.left_parent.index == node.right_parent.index:
            sib = node.sibling()
            if deleted_nodes.__contains__(sib.index):
                if node.left_parent.left_parent != None:
                    parents[node.left_parent.index] = node.left_parent
                    if parents.__contains__(sib.index):
                        del parents[sib.index]
                        self.need_to_visit.discard(sib.index)
            else:
                self.arg.detach(node, sib)
                self.need_to_visit[sib.index] = sib.index
                if sib.left_parent is None:
                    assert sib.left_parent == None
                    if sib.time <= prune.time:
                        assert sib.first_segment != None
                        self.floats[sib.index] = sib.index
                        self.active_links += sib.num_links()
                        self.active_nodes.discard(sib.index)
                        if self.partial_floatings.__contains__(sib.index):#B
                            ch = self.partial_floatings.pop(sib.index)
                            self.active_links -= ch[2]
                    if self.floats.__contains__(node.left_parent.index):#B
                        self.floats.discard(node.left_parent.index)
                        self.active_links -= node.left_parent.num_links()
                    self.need_to_visit.discard(node.left_parent.index)
            deleted_nodes[node.left_parent.index] = node.left_parent.index
            self.need_to_visit.discard(node.left_parent.index)
            if self.arg.__contains__(node.left_parent.index):
                del self.arg.nodes[node.left_parent.index]
                del self.higher_times[node.left_parent.time]
            else:
                assert not self.higher_times.__contains__(node.left_parent.time)
        else:
            parents[node.left_parent.index] = node.left_parent
            parents[node.right_parent.index] = node.right_parent
            deleted_nodes[node.left_parent.index] = node.left_parent.index
            deleted_nodes[node.right_parent.index] = node.right_parent.index
            if self.floats.__contains__(node.left_parent.index):
                self.floats.discard(node.left_parent.index)
                self.active_links -= node.left_parent.num_links()
            self.need_to_visit.discard(node.left_parent.index)
            if self.floats.__contains__(node.left_parent.index):
                self.floats.discard(node.right_parent.index)
                self.active_links -= node.right_parent.num_links()
            self.need_to_visit.discard(node.right_parent.index)
            del self.higher_times[node.left_parent.time]
            del self.arg.nodes[node.left_parent.index]
            del self.arg.nodes[node.right_parent.index]
        return deleted_nodes, parents

    def update_prune_parents(self, prune):
        parents = {prune.index: prune}
        deleted_nodes = bintrees.AVLTree()
        while parents:
            node = parents.pop(min(parents))
            deleted_nodes, parents = self.update_prune_parent(prune, node,
                                                              deleted_nodes, parents)
        prune.left_parent = None
        prune.right_parent = None
        prune.breakpoint = None

    def get_active_links(self):
        self.active_links = 0
        for ind in self.floats:
            num_link = self.arg.__getitem__(ind).num_links()
            self.active_links += num_link
        for ind in self.partial_floatings:
            self.active_links += self.partial_floatings[ind][2]

    def new_event_time(self, lower_time = 0, upper_time = 1, passed_gmrca= False):
        '''
        simulate the time of a new event
        given a time interval (lower, upper), if the new time is in between,
        there is a new event. Three types of events:
            1. coal between two floatings
            2. coal between a floating and one from the rest actve nodes
            3. a rec on a floating lineage
        if passed_gmrca is True: we already passed the original GMRCA,
            so there is no upper_time and any time is acceptable. and the
            events are 1. a coal between two floatings or a rec on a floating :)
        :return:
        '''
        assert len(self.floats) != 0 or self.active_links != 0
        assert self.floats.is_disjoint(self.active_nodes)
        coalrate_bothF = (len(self.floats) * (len(self.floats) - 1))/(4*self.Ne)
        # coalrate_bothF = (len(self.active_nodes) * (len(self.active_nodes)-1)/2)/(2*self.Ne)
        coalrate_1F1rest = (len(self.floats) * len(self.active_nodes))/(2*self.Ne)
        recrate = self.r * self.active_links
        totrate = coalrate_bothF + coalrate_1F1rest + recrate
        new_time = lower_time + random.expovariate(totrate)
        if not passed_gmrca:
            if new_time  >= upper_time:
                # no new event in the time interval
                return False
            else:
                if random.random() < (recrate/totrate):
                    return ["REC", new_time]
                elif random.random() < (coalrate_bothF/(coalrate_1F1rest+coalrate_bothF)):
                    return ["CABF", new_time]
                else:
                    return ["CA1F", new_time]
        else:
            assert coalrate_1F1rest == 0
            assert len(self.floats) > 1
            if random.random() < (recrate/totrate):
                return ["REC", new_time]
            else:
                return ["CABF", new_time]

    def general_incompatibility_check(self,node,  S, s):
        '''
        whether the coalescence of child 1 and child 2 compatible for this snp.
        All are AVLTrees()
        this is applicable to Kuhner and initial
        S: is the parent samples for this segment
        s:  the focal SNP
        node: the parent node
        '''
        valid = True
        D = self.data[s]
        # symmetric difference between S1  and D
        A = S
        symA_D = A.difference(D)
        if len(symA_D) == 0:# subset or equal
            if len(A) == len(D): #put the mutation on this node
                node.snps.__setitem__(s, s)
                # delete s from S_F
                # detach_snps.discard(s)
                # # add to completed_snps
                # completed_snps[s] = s
        elif len(symA_D) == len(A): # distinct
            pass
        else:#
            symD_A = D.difference(A)
            if len(symD_A) > 0: # incompatible
                valid = False
        return valid

    def merge(self, leftchild, rightchild, parent = None):
        '''CA event between two lineages,
        also check compatibility and put mut on node
        :param parent: if None, the parent node is a new node in ARG
            else: the parent already exists in the ARG and we need to update
             the segments and the snps
        '''
        x = self.arg.copy_node_segments(leftchild)
        y = self.arg.copy_node_segments(rightchild)
        x = x.get_first_segment()
        y = y.get_first_segment()
        assert x is not None
        assert y is not None
        index =  self.arg.new_name()
        if parent == None:
            node = self.arg.alloc_node(index, time, x.node, y.node)
        else:
            node = parent
            node.snps.clear()
        self.arg.coal[node.index] = node.index
        x.node.left_parent = node
        x.node.right_parent = node
        y.node.left_parent = node
        y.node.right_parent = node
        z = None
        defrag_required = False
        valid = True
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
                    else:# check compatibility, add snps
                        seg_snps = alpha.get_seg_variants(self.data)
                        for snp in seg_snps:#intersect_variants:
                            valid = self.general_incompatibility_check(node,  alpha.samples, snp)
                            if not valid:# break for
                                break
                        if not valid:# break while
                            break
                    if x.right == right:
                        x = x.next
                    else:
                        x.left = right
                    if y.right == right:
                        y = y.next
                    else:
                        y.left = right
            if alpha is not None:
                # if z is None:
                #     self.parent_nodes[p] = alpha
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
            if parent == None:
                self.arg.store_node(z, node)
            else:# already exist in ARG, update segments
                node.first_segment = z
                while z is not None:
                    z.node = node
                    z = z.next
        else:
            if parent == None:
                self.arg.add(node)
            else:
                node.first_segment = None
            self.arg.roots[node.index] = node.index
        return valid, node

    def make_float(self, child, oldbreakpoint, oldleftparent, oldrightparent,
                       floatParent, followparent, old_left, old_right, ch):
        '''to be used in new_recombination'''
        followparent.breakpoint = oldbreakpoint
        followparent.left_parent = oldleftparent
        followparent.right_parent = oldrightparent
        oldleftparent.update_child(child, followparent)
        oldrightparent.update_child(child, followparent)
        # assert (old_right <= floatParent.first_segment.left) or \
        #        (floatParent.get_tail().right <= old_left)
        #update partial floating
        rphead = followparent.first_segment
        rptail = followparent.get_tail()
        self.original_interval[followparent.index] = [old_left, old_right]
        self.add_to_partial_floating(followparent, old_left,
                old_right, rphead.left, rptail.right)
        self.floats[floatParent.index] = floatParent.index
        # active links
        self.active_links -= ch[2]
        self.active_links += floatParent.num_links()
        # diff = ch[2] - floatParent.num_links()
        # assert diff>=0
        # self.active_links -= diff
        self.need_to_visit[floatParent.index] = floatParent.index
        self.need_to_visit[followparent.index] = followparent.index
        self.need_to_visit.discard(child.index)
        self.active_nodes.discard(child.index)

    def depricated_new_recombination(self, t):
        '''
        The new event is a recomb at t
        1. choose a lineage to put the rec on from
            self.floats and self.partial_floatings,
            proportional to their num of links
        2. if a floating is chosen, randomly choose a breakpoint
                it and split to two and put both parents in self.floats
            Else: it is from a partial floating,
                choose a breakpoint randomly from the new sites (a', a), (b, b')
                split the lineage to two. The parent with all new sites will
                float and the other follows the child path.
        '''
        valid = True
        partial_links = [self.partial_floatings[item][2] for item in self.partial_floatings]
        float_keys = list(self.floats.keys())
        float_links  = [self.arg.__getitem__(item).num_links() for item in float_keys]
        assert sum(partial_links) + sum(float_links) == self.active_links
        type = "f" #from floating
        if self.partial_floatings:
            type = random.choices(["pf", "f"], [sum(partial_links), sum(float_links)])[0]
        if type == "f": #from floating
            child_ind = random.choices(float_keys, float_links)[0]
            child = self.arg.__getitem__(child_ind)
            # choose a breakpoint on child
            head = child.first_segment
            tail = child.get_tail()
            if tail.right - head.left>1:
                break_point = random.choice(range(head.left + 1,
                                                  tail.right))
                leftparent, rightparent = self.split_node(child, break_point, t)
                assert leftparent.get_tail().right <= break_point
                assert break_point <= rightparent.first_segment.left
                self.floats.discard(child_ind)
                self.floats[leftparent.index] = leftparent.index
                self.floats[rightparent.index] = rightparent.index
                self.need_to_visit.discard(child_ind)
                self.need_to_visit[leftparent.index] = leftparent.index
                self.need_to_visit[rightparent.index] = rightparent.index
                self.active_links -= (child.num_links() -\
                                      (leftparent.num_links() + rightparent.num_links()))
            else:
                valid = False
        else: # rec on one of the partials
            child_ind = random.choices(list(self.partial_floatings.keys()), partial_links)[0]
            ch = self.partial_floatings.pop(child_ind)
            if ch[2]>1:
                bp = random.choice(range(ch[2]))
                child = self.arg.__getitem__(child_ind)
                head = child.first_segment
                tail = child.get_tail()
                a = ch[0]
                b = ch[1]
                oldleftparent = child.left_parent
                oldrightparent = child.right_parent
                oldbreakpoint = child.breakpoint
                assert self.original_interval.__contains__(child.index)
                old_interval = self.original_interval.pop(child.index)
                old_left = old_interval[0]
                old_right = old_interval[1]
                if a != None and b != None:
                    assert a == old_left and b == old_right
                    assert ch[2] == (a - head.left) + (tail.right - b)
                    if bp < a - head.left:
                        # the left parent is floating
                        break_point = head.left + bp + 1
                        assert break_point <= old_left or old_right <= break_point
                        assert head.left < break_point <= a
                        leftparent, rightparent = self.split_node(child, break_point, t)
                        assert leftparent.get_tail().right <= break_point
                        assert break_point <= rightparent.first_segment.left
                        self.make_float(child, oldbreakpoint, oldleftparent, oldrightparent,
                                            leftparent, rightparent, old_left, old_right, ch)
                    else:
                        #right parent is floating
                        break_point = bp - (a- head.left) +b
                        assert b<= break_point < tail.right
                        assert break_point <= old_left or old_right <= break_point
                        leftparent, rightparent = self.split_node(child, break_point, t)
                        assert leftparent.get_tail().right <= break_point
                        assert break_point <= rightparent.first_segment.left
                        self.make_float(child, oldbreakpoint, oldleftparent, oldrightparent,
                                                rightparent, leftparent, old_left, old_right, ch)
                elif a!= None:
                    break_point = head.left + bp + 1
                    assert tail.right - head.left -1 == ch[2] or old_left- head.left == ch[2]
                    assert break_point <= old_left or old_right <= break_point
                    # assert head.left < break_point <= a
                    leftparent, rightparent = self.split_node(child, break_point, t)
                    assert leftparent.get_tail().right <= break_point
                    assert break_point <= rightparent.first_segment.left
                    # choose which one to float
                    if head.left< old_left and old_left< tail.right: # case 3 in notebook
                        self.make_float(child, oldbreakpoint, oldleftparent, oldrightparent,
                           leftparent, rightparent, old_left, old_right, ch)
                    elif (head.left < old_left and tail.right<= old_left) or \
                            (old_right <= head.left and old_right < tail.right):# case 5, 8, 4, 7
                        if random.random() < 0.5:# leftparent is floating
                            self.make_float(child, oldbreakpoint, oldleftparent, oldrightparent,
                                leftparent, rightparent, old_left, old_right, ch)
                        else:#right parent is floating
                            self.make_float(child, oldbreakpoint, oldleftparent, oldrightparent,
                                rightparent, leftparent, old_left, old_right, ch)
                    else:
                        raise ValueError("the new interval violates the vaild cases")
                elif b!= None:
                    #right parent is floating
                    assert tail.right - old_right == ch[2]
                    break_point = bp + b
                    assert break_point <= old_left or old_right <= break_point
                    assert b<= break_point < tail.right
                    leftparent, rightparent = self.split_node(child, break_point, t)
                    assert leftparent.get_tail().right <= break_point
                    assert break_point <= rightparent.first_segment.left
                    self.make_float(child, oldbreakpoint, oldleftparent, oldrightparent,
                                rightparent, leftparent, old_left, old_right, ch)
                else:
                    raise ValueError("both a and b are None")
            else:
                valid = False
        return valid

    def new_recombination(self, t):
        '''
        New recombination in Kuhner move
        1. choose a lineage to put the rec on from
            self.floats and self.partial_floatings,
            proportional to their num of links
        2. if a floating is chosen, randomly choose a breakpoint
                it and split to two and put both parents in self.floats
            Else: it is from a partial floating,
                choose a breakpoint randomly from the new sites (a', a), (b, b')
                split the lineage to two. The parent with all new sites will
                float and the other follows the child path.

         10 sep 2020: in partially_floating the recomb parent is chosen
            randomly to float, reghardless if the float parent's sites
            are all new sites or a part of it are new sites.
        The new event is a recomb at t
        '''
        valid = True
        partial_links = [self.partial_floatings[item][2] for item in self.partial_floatings]
        partial_keys = list(self.partial_floatings.keys())
        float_keys = list(self.floats.keys())
        float_links  = [self.arg.__getitem__(item).num_links() for item in float_keys]
        partial_links.extend(float_links)
        partial_keys.extend(float_keys)
        assert sum(partial_links) == self.active_links
        child_ind = random.choices(partial_keys, partial_links)[0]
        if self.floats.__contains__(child_ind):
            child = self.arg.__getitem__(child_ind)
            # choose a breakpoint on child
            head = child.first_segment
            tail = child.get_tail()
            if tail.right - head.left>1:
                break_point = random.choice(range(head.left + 1,
                                                  tail.right))
                leftparent, rightparent = self.split_node(child, break_point, t)
                assert leftparent.get_tail().right <= break_point
                assert break_point <= rightparent.first_segment.left
                self.floats.discard(child_ind)
                self.floats[leftparent.index] = leftparent.index
                self.floats[rightparent.index] = rightparent.index
                self.need_to_visit.discard(child_ind)
                self.need_to_visit[leftparent.index] = leftparent.index
                self.need_to_visit[rightparent.index] = rightparent.index
                self.active_links -= (child.num_links() - leftparent.num_links())
                self.active_links += rightparent.num_links()
            else:
                valid = False
        else: # rec on one of the partials
            assert self.partial_floatings.__contains__(child_ind)
            ch = self.partial_floatings.pop(child_ind)
            if ch[2]>1:
                bp = random.choice(range(ch[2]))
                child = self.arg.__getitem__(child_ind)
                head = child.first_segment
                tail = child.get_tail()
                a = ch[0]
                b = ch[1]
                oldleftparent = child.left_parent
                oldrightparent = child.right_parent
                oldbreakpoint = child.breakpoint
                assert self.original_interval.__contains__(child.index)
                old_interval = self.original_interval.pop(child.index)
                old_left = old_interval[0]
                old_right = old_interval[1]
                if a != None and b != None:
                    assert a == old_left and b == old_right
                    assert ch[2] == (a - head.left) + (tail.right - b)
                    if bp < a - head.left:
                        # the left parent is floating
                        break_point = head.left + bp + 1
                        assert break_point <= old_left or old_right <= break_point
                        assert head.left < break_point <= a
                        leftparent, rightparent = self.split_node(child, break_point, t)
                        assert leftparent.get_tail().right <= break_point
                        assert break_point <= rightparent.first_segment.left
                        # 10 sep2020: choose float parent randomly
                        if random.random() < 0.5:# leftparent is floating
                            self.make_float(child, oldbreakpoint, oldleftparent, oldrightparent,
                                                leftparent, rightparent, old_left, old_right, ch)
                        else:#rightparent is floating
                            self.make_float(child, oldbreakpoint, oldleftparent, oldrightparent,
                                                 rightparent,leftparent, old_left, old_right, ch)
                    else:
                        #right parent is floating
                        break_point = bp - (a- head.left) +b
                        assert b<= break_point < tail.right
                        assert break_point <= old_left or old_right <= break_point
                        leftparent, rightparent = self.split_node(child, break_point, t)
                        assert leftparent.get_tail().right <= break_point
                        assert break_point <= rightparent.first_segment.left
                        # 10 sep2020: choose float parent randomly
                        if random.random() < 0.5:# leftparent is floating
                            self.make_float(child, oldbreakpoint, oldleftparent, oldrightparent,
                                                leftparent, rightparent, old_left, old_right, ch)
                        else:
                            self.make_float(child, oldbreakpoint, oldleftparent, oldrightparent,
                                                 rightparent,leftparent, old_left, old_right, ch)
                elif a!= None:
                    break_point = head.left + bp + 1
                    assert tail.right - head.left -1 == ch[2] or old_left- head.left == ch[2]
                    assert break_point <= old_left or old_right <= break_point
                    # assert head.left < break_point <= a
                    leftparent, rightparent = self.split_node(child, break_point, t)
                    assert leftparent.get_tail().right <= break_point
                    assert break_point <= rightparent.first_segment.left
                    # choose which one to float
                    if (head.left< old_left and old_left< tail.right):# case 3 in notebook
                        # 10 sep2020: choose float parent randomly
                        if random.random() < 0.5:# leftparent is floating
                            self.make_float(child, oldbreakpoint, oldleftparent, oldrightparent,
                                                leftparent, rightparent, old_left, old_right, ch)
                        else:
                            self.make_float(child, oldbreakpoint, oldleftparent, oldrightparent,
                                                 rightparent,leftparent, old_left, old_right, ch)
                    elif (head.left < old_left and tail.right<= old_left) or \
                            (old_right <= head.left and old_right < tail.right):# case 5, 8, 4, 7
                        if random.random() < 0.5:# leftparent is floating
                            self.make_float(child, oldbreakpoint, oldleftparent, oldrightparent,
                                leftparent, rightparent, old_left, old_right, ch)
                        else:#right parent is floating
                            self.make_float(child, oldbreakpoint, oldleftparent, oldrightparent,
                                rightparent, leftparent, old_left, old_right, ch)
                    else:
                        raise ValueError("the new interval violates the vaild cases")
                elif b!= None:
                    #right parent is floating
                    assert tail.right - old_right == ch[2]
                    break_point = bp + b
                    assert break_point <= old_left or old_right <= break_point
                    assert b<= break_point < tail.right
                    leftparent, rightparent = self.split_node(child, break_point, t)
                    assert leftparent.get_tail().right <= break_point
                    assert break_point <= rightparent.first_segment.left
                    # 10 sep2020: choose float parent randomly
                    if random.random() < 0.5:# leftparent is floating
                        self.make_float(child, oldbreakpoint, oldleftparent, oldrightparent,
                                            leftparent, rightparent, old_left, old_right, ch)
                    else:
                        self.make_float(child, oldbreakpoint, oldleftparent, oldrightparent,
                                             rightparent,leftparent, old_left, old_right, ch)
                else:
                    raise ValueError("both a and b are None")
            else:
                valid = False
        return valid

    def coal_bothF(self, t):
        '''
        coalescence of two floating lineages at t
        1. choose two lineages in self.floats
        2. update materials and check compatibility
            also put new mut on parent node if any
        '''
        valid = True
        assert len(self.floats) > 1
        inds = random.sample(list(self.floats), 2)
        leftchild = self.arg.__getitem__(inds[0])
        rightchild = self.arg.__getitem__(inds[1])
        assert leftchild.breakpoint == None
        assert rightchild.breakpoint == None
        valid, parent = self.merge(leftchild, rightchild)
        if valid:
            parent.time = t
            if parent.first_segment != None:
                self.floats[parent.index] = parent.index
                self.need_to_visit[parent.index] = parent.index
                self.active_links -= (leftchild.num_links() + rightchild.num_links())
                self.active_links += parent.num_links()
            else:
                self.active_links -= (leftchild.num_links() + rightchild.num_links())
            self.floats.discard(inds[0])
            self.floats.discard(inds[1])
            self.need_to_visit.discard(inds[0])
            self.need_to_visit.discard(inds[1])
        return valid

    def coal_1F_1active(self, t):
        '''
        coalescence between one floating node
        and one active node
        1. choose one from active, one from flaots
        '''
        valid = True
        active_ind = random.choice(list(self.active_nodes))
        float_ind = random.choice(list(self.floats))
        if  self.verbose:
            print("CA1F", "floats:", float_ind, "reattaches to", active_ind)
        activenode = self.arg.__getitem__(active_ind)
        floatnode = self.arg.__getitem__(float_ind)
        #--- start
        if self.original_interval.__contains__(activenode.index):
            old_interval = self.original_interval.pop(activenode.index)
            old_left = old_interval[0]
            old_right = old_interval[1]
        else:
            old_left = activenode.first_segment.left
            old_right = activenode.get_tail().right
            assert self.original_ARG.__getitem__(activenode.index).first_segment.left == old_left
            assert self.original_ARG.__getitem__(activenode.index).get_tail().right == old_right
        float_head = floatnode.first_segment
        float_tail = floatnode.get_tail()
        oldleftparent = activenode.left_parent
        oldrightparent = activenode.right_parent
        oldbreakpoint = activenode.breakpoint
        activenode.breakpoint = None
        floatnode.breakpoint = None
        valid, parent = self.merge(activenode, floatnode)
        parent.time = t
        if valid:
            parent.breakpoint = oldbreakpoint
            parent.left_parent = oldleftparent
            parent.right_parent = oldrightparent
            oldleftparent.update_child(activenode, parent)
            oldrightparent.update_child(activenode, parent)
            if parent.first_segment != None:
                parent_head = parent.first_segment
                parent_tail = parent.get_tail()
                if parent.left_parent == None: #future floating
                    self.floats[parent.index] = parent.index
                    self.active_links += parent_tail.right - parent_head.left-1
                else:
                    self.original_interval[parent.index] = [old_left, old_right]
                    self.add_to_partial_floating(parent, old_left,
                                old_right, parent_head.left, parent_tail.right)
                self.need_to_visit[parent.index] = parent.index
            elif parent.left_parent != None:
                # delete its parent--> same way we did for prune
                # self.update_prune_parents(parent)
                self.set_nodechild_to_None(parent)
                self.need_to_visit.discard(parent.index)
            # else:# root
            self.active_links -= float_tail.right - float_head.left -1
            if self.partial_floatings.__contains__(active_ind):
                assert self.partial_floatings[active_ind][2] > 0
                ch = self.partial_floatings.pop(active_ind)
                self.active_links -= ch[2]
            self.active_nodes.discard(active_ind)
            self.floats.discard(float_ind)
            self.need_to_visit.discard(active_ind)
            self.need_to_visit.discard(float_ind)
        return valid

    def apply_new_event(self, newevent_type, new_time):
        valid = True
        if newevent_type == "REC":
            valid = self.new_recombination(new_time)
        elif newevent_type == "CABF":
            valid = self.coal_bothF(new_time)
        else:#CA1F
            valid = self.coal_1F_1active(new_time)
        return valid

    def split_node_kuhner(self, child, left_parent, right_parent):
        '''
        update the materials of leftparent and rightparent
        this is for an exsting recombination event in the ARG
        example: in kuhner when there is no new event.
        '''
        s = self.arg.copy_node_segments(child) #child segments
        y = self.find_break_seg(s, child.breakpoint)
        if y != None:
            x = y.prev
            if y.left < child.breakpoint < y.right:
                z = self.arg.alloc_segment(child.breakpoint, y.right, y.node,
                                           y.samples, None, y.next)
                if y.next is not None:
                    y.next.prev = z
                y.next = None
                y.right = child.breakpoint
                lhs_tail = y
            elif x is not None:
                assert x.right <= child.breakpoint <= y.left
                x.next = None
                y.prev = None
                z = y
                lhs_tail = x
            else: # first parent is empty
                z = y
                lhs_tail = None
        else: # second parent is empty
            z = None
            lhs_tail = s
        #=======
        # the parents already in ARG just update segments
        if lhs_tail != None:
            seg = lhs_tail.get_first_segment()
            left_parent.first_segment = seg
            while seg is not None:
                seg.node = left_parent
                seg = seg.next
        else:
            left_parent.first_segment = None
        if z != None:
            seg = z.get_first_segment()
            right_parent.first_segment = z
            while seg is not None:
                seg.node = right_parent
                seg = seg.next
        else:
            right_parent.first_segment = None
        return left_parent, right_parent

    def add_to_partial_floating(self,node, old_left,
                                old_right, new_left, new_right):
        '''check and add the node to partial floating'''
        assert node.first_segment is not None
        if node.left_parent is None:
            # add to floats
            assert not self.floats.__contains__(node.index)
            self.floats[node.index] = node.index
            self.active_links += node.num_links()
        else:
            a = None
            b = None
            new_links = 0
            if (new_left < old_left) and (new_right > old_left):# cases 1, 3
                if new_right <= old_right:# case 1
                    a = old_left
                    new_links += old_left - new_left
                else:# 3
                    a = old_left
                    b = old_right
                    new_links += (old_left - new_left) + (new_right - old_right)
            elif (new_left < old_left)  and (new_right <= old_left):# 5, 8
                new_links += new_right - new_left - 1
                if new_links > 0:
                    a = new_left
            elif (new_left >= old_left) and (new_left< old_right):# 2, 6
                if new_right > old_right: # 2
                    b = old_right
                    new_links += new_right - old_right
            elif (new_left >= old_left) and (new_left >= old_right):# 4, 7
                new_links += new_right - new_left -1
                if new_links > 0:
                    a = new_left
            if (a != None) or (b != None):
                self.partial_floatings[node.index] = [a, b, new_links]
                self.active_links += new_links
            self.active_nodes[node.index] = node.index

    def check_material(self, leftchild, rightchild,
                       leftparent, rightparent):
        '''
         for a alredy existing event, update the ancestral material
         this is for kuhner move, the cases with no event between a time interval
        '''
        valid = True
        if leftparent.index != rightparent.index:# rec
            assert leftchild.index == rightchild.index
            child = leftchild
            oldleftparent_left = leftparent.first_segment.left
            oldleftparent_right = leftparent.get_tail().right
            oldrightparent_left = rightparent.first_segment.left
            oldrightparent_right = rightparent.get_tail().right
            assert oldleftparent_left ==\
                   self.original_ARG.__getitem__(leftparent.index).first_segment.left
            assert oldleftparent_right ==\
                   self.original_ARG.__getitem__(leftparent.index).get_tail().right
            assert oldrightparent_left ==\
                   self.original_ARG.__getitem__(rightparent.index).first_segment.left
            assert oldrightparent_right ==\
                   self.original_ARG.__getitem__(rightparent.index).get_tail().right
            leftparent, rightparent = self.split_node_kuhner(child, leftparent, rightparent)
            if self.partial_floatings.__contains__(child.index):
                    self.active_links -= self.partial_floatings.pop(child.index)[2]
            if leftparent.first_segment != None and\
                    rightparent.first_segment != None:
                leftparent_head = leftparent.first_segment
                leftparent_tail = leftparent.get_tail()
                rightparent_head = rightparent.first_segment
                rightparent_tail = rightparent.get_tail()
                assert leftparent_tail.right <= rightparent_head.left
                self.add_to_partial_floating(leftparent, oldleftparent_left,
                                oldleftparent_right,leftparent_head.left, leftparent_tail.right)
                self.original_interval[leftparent.index] = [oldleftparent_left, oldleftparent_right]
                self.add_to_partial_floating(rightparent, oldrightparent_left,
                                oldrightparent_right,rightparent_head.left, rightparent_tail.right)
                self.original_interval[rightparent.index] = [oldrightparent_left, oldrightparent_right]
                assert leftparent.left_parent != None
                assert rightparent.left_parent != None
                self.active_nodes[leftparent.index] = leftparent.index
                self.active_nodes[rightparent.index] = rightparent.index
                # elif leftparent.left_parent != None and rightparent.left_parent != None:
                #     self.active_nodes[leftparent.index] = leftparent.index
                #     self.active_nodes[rightparent.index] = rightparent.index
                # elif leftparent.left_parent != None:# right parent was future floating
                #     assert not self.floats.__contains__(rightparent.index)
                #     self.floats[rightparent.index] = rightparent.index
                #     self.active_links += rightparent.num_links()
                #     self.active_nodes[leftparent.index] = leftparent.index
                # elif rightparent.left_parent != None:
                #     assert not self.floats.__contains__(leftparent.index)
                #     self.floats[leftparent.index] = leftparent.index
                #     self.active_links += leftparent.num_links()
                #     self.active_nodes[rightparent.index] = rightparent.index
                # else:# both are None
                #     self.floats[leftparent.index] = leftparent.index
                #     self.active_links += leftparent.num_links()
                #     self.floats[rightparent.index] = rightparent.index
                #     self.active_links += rightparent.num_links()
                self.need_to_visit[leftparent.index] = leftparent.index
                self.need_to_visit[rightparent.index] = rightparent.index
                self.need_to_visit.discard(child.index)
                self.active_nodes.discard(child.index)
            elif leftparent.first_segment != None and rightparent.first_segment == None:
                assert leftparent.left_parent != None
                leftparent.reconnect(child)
                leftparent_head = leftparent.first_segment
                leftparent_tail = leftparent.get_tail()
                assert leftparent_head.left == child.first_segment.left
                assert leftparent_tail.right == child.get_tail().right
                self.add_to_partial_floating(child, oldleftparent_left,
                                oldleftparent_right,leftparent_head.left, leftparent_tail.right)
                self.original_interval[child.index] = [oldleftparent_left, oldleftparent_right]
                # child.left_parent = leftparent.left_parent
                # child.right_parent = leftparent.right_parent
                # child.breakpoint = leftparent.breakpoint
                # parent = leftparent.left_parent
                # if parent is None:
                #     #child will be floating
                #     self.floats[child.index] = child.index
                #     self.active_links += child.num_links()
                #     self.active_nodes.discard(child.index)
                #     if self.partial_floatings.__contains__(child.index):
                #         ch = self.partial_floatings.pop(child.index)
                #         self.active_links -= ch[2]
                #     assert self.need_to_visit.__contains__(child.index)
                # else:
                #     parent.update_child(leftparent, child)
                #     parent = leftparent.right_parent
                #     parent.update_child(leftparent, child)
                    # child stays active and if in partial-->stays same
                self.set_nodechild_to_None(rightparent)
                del self.arg.nodes[leftparent.index]
                del self.arg.nodes[rightparent.index]
                # self.update_prune_parents(rightparent)
                #parents might alredy added to need_to_visit
                self.need_to_visit.discard(rightparent.index)
                self.need_to_visit.discard(leftparent.index)
            elif rightparent.first_segment != None and leftparent.first_segment == None:
                assert rightparent.left_parent != None
                rightparent.reconnect(child)
                rightparent_head = rightparent.first_segment
                rightparent_tail = rightparent.get_tail()
                assert rightparent_head.left == child.first_segment.left
                assert rightparent_tail.right == child.get_tail().right
                self.add_to_partial_floating(child, oldrightparent_left,
                                oldrightparent_right,rightparent_head.left, rightparent_tail.right)
                self.original_interval[child.index] = [oldrightparent_left, oldrightparent_right]
                self.set_nodechild_to_None(leftparent)
                del self.arg.nodes[leftparent.index]
                del self.arg.nodes[rightparent.index]
                # self.update_prune_parents(rightparent)
                #parents might alredy added to need_to_visit
                self.need_to_visit.discard(rightparent.index)
                self.need_to_visit.discard(leftparent.index)
                # child.left_parent = rightparent.left_parent
                # child.right_parent = rightparent.right_parent
                # child.breakpoint = rightparent.breakpoint
                # parent = rightparent.left_parent
                # if parent is None:
                #     #child will be floating
                #     self.floats[child.index] = child.index
                #     self.active_links += child.num_links()
                #     self.active_nodes.discard(child.index)
                #     if self.partial_floatings.__contains__(child.index):
                #         ch = self.partial_floatings.pop(child.index)
                #         self.active_links -= ch[2]
                #     assert self.need_to_visit.__contains__(child.index)
                # else:
                #     parent.update_child(rightparent, child)
                #     parent = rightparent.right_parent
                #     parent.update_child(rightparent, child)
                #     # child stays active
                # del self.arg.nodes[leftparent.index]
                # del self.arg.nodes[rightparent.index]
                # self.update_prune_parents(leftparent)
                # self.need_to_visit.discard(leftparent.index)
                # self.need_to_visit.discard(rightparent.index)
            else: # both None
                raise ValueError("at least one parent must have segment")
        else: #coal
            assert leftchild.index != rightchild.index
            parent = leftparent
            if parent.first_segment != None:
                oldparent_left = parent.first_segment.left
                oldparent_right = parent.get_tail().right
                assert oldparent_left ==\
                        self.original_ARG.__getitem__(parent.index).first_segment.left
                assert oldparent_right ==\
                       self.original_ARG.__getitem__(parent.index).get_tail().right
            else:
                oldparent_left = None
                oldparent_right = None
            valid, parent = self.merge(leftchild, rightchild, parent)
            if valid:
                if parent.first_segment!= None:
                    parent_head = parent.first_segment
                    parent_tail = parent.get_tail()
                    if parent.left_parent == None: #future floating
                        self.floats[parent.index] = parent.index
                        self.active_links += parent_tail.right - parent_head.left-1
                    else:
                        assert oldparent_left != None
                        assert oldparent_right != None
                        self.add_to_partial_floating(parent, oldparent_left, oldparent_right,
                                                     parent_head.left, parent_tail.right)
                        self.original_interval[parent.index] = [oldparent_left, oldparent_right]
                    self.need_to_visit[parent.index] = parent.index
                elif parent.left_parent != None:
                    #this is new root
                    # self.update_prune_parents(parent)
                    self.set_nodechild_to_None(parent)
                    self.need_to_visit.discard(parent.index)
                    # del self.arg.nodes[parent.index]
                else:# parent is root
                    self.need_to_visit.discard(parent.index) # might have been added for future float
                #----
                if self.partial_floatings.__contains__(leftchild.index):
                    self.active_links -= self.partial_floatings[leftchild.index][2]
                    del self.partial_floatings[leftchild.index]
                if self.partial_floatings.__contains__(rightchild.index):
                    self.active_links -= self.partial_floatings[rightchild.index][2]
                    del self.partial_floatings[rightchild.index]
                self.active_nodes.discard(leftchild.index)
                self.active_nodes.discard(rightchild.index)
                self.need_to_visit.discard(leftchild.index)
                self.need_to_visit.discard(rightchild.index)
        return valid

    def set_nodechild_to_None(self, child):
        '''
        For NAM or pruned nodes:
        set the parent.child of node_child to None.
        Also, set the child.parent to None'''
        if child.left_parent != None:
            if child.left_parent.left_child != None:
                if child.left_parent.left_child.index == child.index:
                    child.left_parent.left_child = None
            if child.left_parent.right_child != None:
                if child.left_parent.right_child.index == child.index:
                    child.left_parent.right_child = None
            if child.right_parent.left_child != None:
                if child.right_parent.left_child.index == child.index:
                    child.right_parent.left_child = None
            if child.right_parent.right_child != None:
                if child.right_parent.right_child.index == child.index:
                    child.right_parent.right_child = None
            self.need_to_visit[child.left_parent.index] = child.left_parent.index
            self.need_to_visit[child.right_parent.index] = child.right_parent.index
            child.left_parent = None
            child.right_parent = None
            child.breakpoint = None

    def kuhner(self):
        '''
        1.randomly choose a lineage other than the root,
            a.put it in need_check, floats
        2. find all the time and nodes greater than prune.time
        3. find mutations need to be checked in future
        4. if prune is a child of CA:
            a) detach prune
            b)remove P
            c) put C in need_check,
            d) if P root, put C in future floating
            else: prune parents, and all the further parents
        '''
        #---- forward transition
        self.arg.dump(path = self.outpath, file_name = 'arg.arg')
        self.original_ARG =  self.arg.load(path = self.outpath+'/arg.arg')
        # self.original_ARG = self.arg.copy()
        # self.original_ARG = copy.deepcopy(self.arg)#copy.deepcopy(self.arg)

        self.transition_prob.kuhner_num_nodes(self.arg.__len__() - len(self.arg.roots))
        valid = True
        prune = self.add_choose_child()
        if  self.verbose:
            print("prune is ", prune.index)
        assert prune.left_parent != None
        self.floats[prune.index] = prune.index
        self.need_to_visit[prune.index] = prune.index
        self.need_to_visit[prune.left_parent.index] = prune.left_parent.index
        self.find_active_nodes(prune.time)# 2
        # self.update_prune_parents(prune)# 4
        self.set_nodechild_to_None(prune)
        self.get_active_links()# active links
        self.active_nodes = self.active_nodes.difference(self.floats)
        lower_time = prune.time
        while self.need_to_visit:
            if self.higher_times:# time interval
                upper_time = min(self.higher_times)
                parent_indexes = self.higher_times[upper_time]
                # check for future floating
                if self.active_links > 0 or len(self.floats) > 0:
                    new_event = self.new_event_time(lower_time, upper_time)
                else:
                    new_event = False
                if new_event:
                    new_time = new_event[1]; newevent_type = new_event[0]
                    if self.verbose:
                        print("new event type", newevent_type, "at time", new_time)
                    valid = self.apply_new_event(newevent_type, new_time)
                    lower_time = new_time
                    if not valid:
                        break
                else: # no new event
                    del self.higher_times[upper_time]
                    assert 0 < len(parent_indexes) <= 2
                    if len(parent_indexes) == 2:#rec
                        parent1 = self.arg.__getitem__(parent_indexes.pop())
                        parent2 = self.arg.__getitem__(parent_indexes.pop())
                        assert parent1.index != parent2.index
                        assert parent1.left_parent != None
                        assert parent2.left_parent != None
                        if parent1.left_child != None:#
                            assert parent1.left_child.index == parent1.right_child.index
                            child = parent1.left_child
                            leftparent = child.left_parent
                            rightparent = child.right_parent
                            if self.need_to_visit.__contains__(child.index):
                                valid = self.check_material(child, child, leftparent, rightparent)
                                if not valid:
                                    break
                            else:
                                self.active_nodes.discard(child.index)
                                assert leftparent.first_segment != None
                                assert rightparent.first_segment != None
                                assert not self.partial_floatings.__contains__(child.index)
                                assert not self.original_interval.__contains__(child.index)
                                self.active_nodes[leftparent.index] = leftparent.index
                                self.active_nodes[rightparent.index] = rightparent.index
                        else: # child is NAM or pruned
                            self.set_nodechild_to_None(parent1)
                            self.set_nodechild_to_None(parent2)
                            self.need_to_visit.discard(parent1.index)
                            self.need_to_visit.discard(parent2.index)
                            del self.arg.nodes[parent1.index]
                            del self.arg.nodes[parent2.index]
                    else: # coal
                        parent = self.arg.__getitem__(parent_indexes.pop())
                        assert len(parent_indexes) == 0
                        leftchild = parent.left_child
                        rightchild = parent.right_child
                        if leftchild != None and rightchild != None:
                            assert leftchild.index != rightchild.index
                            assert leftchild.first_segment != None
                            assert rightchild.first_segment != None
                            if self.need_to_visit.__contains__(leftchild.index) or\
                                self.need_to_visit.__contains__(rightchild.index):
                                valid = self.check_material(leftchild, rightchild, parent, parent)
                                if not valid:
                                    break
                            else:
                                assert not self.partial_floatings.__contains__(leftchild.index)
                                assert not self.partial_floatings.__contains__(rightchild.index)
                                assert not self.original_interval.__contains__(leftchild.index)
                                assert not self.original_interval.__contains__(rightchild.index)
                                self.active_nodes.discard(leftchild.index)
                                self.active_nodes.discard(rightchild.index)
                                if parent.left_parent != None:
                                    assert parent.first_segment != None
                                    self.active_nodes[parent.index] = parent.index
                                elif parent.first_segment != None:
                                    #it is a future floating
                                    assert self.need_to_visit.__contains__(parent.index)
                                    self.floats[parent.index] = parent.index
                                    self.active_links += parent.num_links()
                                else:# might have been added to check future floating
                                    self.need_to_visit.discard(parent.index)
                        elif leftchild == None and rightchild != None:
                            assert rightchild.first_segment != None
                            if self.partial_floatings.__contains__(rightchild.index):#B
                                self.active_links -= self.partial_floatings.pop(rightchild.index)[2]
                            if parent.left_parent == None:#
                                self.floats[rightchild.index] = rightchild.index
                                self.active_links += rightchild.num_links()
                                self.active_nodes.discard(rightchild.index)
                                self.need_to_visit[rightchild.index] = rightchild.index
                            else:
                                self.original_interval[rightchild.index] = [parent.first_segment.left,
                                                                            parent.get_tail().right]
                                self.add_to_partial_floating(rightchild,parent.first_segment.left,
                                        parent.get_tail().right, rightchild.first_segment.left,
                                                             rightchild.get_tail().right)
                                parent.reconnect(rightchild)
                                self.need_to_visit[rightchild.index] = rightchild.index
                                assert self.active_nodes.__contains__(rightchild.index)
                            self.need_to_visit.discard(parent.index)
                            del self.arg.nodes[parent.index]
                        elif leftchild != None and rightchild == None:
                            assert leftchild.first_segment != None
                            if self.partial_floatings.__contains__(leftchild.index):#B
                                self.active_links -= self.partial_floatings.pop(leftchild.index)[2]
                            if parent.left_parent == None:#
                                self.floats[leftchild.index] = leftchild.index
                                self.active_links += leftchild.num_links()
                                self.active_nodes.discard(leftchild.index)
                                self.need_to_visit[leftchild.index] = leftchild.index
                            else:
                                self.original_interval[leftchild.index] = [parent.first_segment.left,
                                                                           parent.get_tail().right]
                                self.add_to_partial_floating(leftchild, parent.first_segment.left,
                                        parent.get_tail().right, leftchild.first_segment.left,
                                                             leftchild.get_tail().right)
                                parent.reconnect(leftchild)
                                self.need_to_visit[leftchild.index] = leftchild.index
                                assert self.active_nodes.__contains__(leftchild.index)
                            self.need_to_visit.discard(parent.index)
                            del self.arg.nodes[parent.index]
                        else: # both children are None
                            if parent.left_parent == None:#
                                self.need_to_visit.discard(parent.index)
                                del self.arg.nodes[parent.index]
                            else:
                                self.set_nodechild_to_None(parent)
                                self.need_to_visit.discard(parent.index)
                                del self.arg.nodes[parent.index]
                            # check for paritial floating
                    lower_time = upper_time
            else: #passed gmrca
                if self.verbose:
                    print("PASSED GMRCA")
                assert len(self.floats) > 1
                assert self.active_nodes.is_empty()
                new_event = self.new_event_time(lower_time, passed_gmrca=True)
                new_time = new_event[1]; newevent_type = new_event[0]
                if self.verbose:
                    print("new event type", newevent_type, "at time", new_time)
                valid = self.apply_new_event(newevent_type, new_time)
                lower_time = new_time
                if not valid:
                    break
        if  self.verbose:
            print("valid is", valid)
        if valid:
            assert self.floats.is_empty()
            assert self.active_links == 0
            # assert self.active_nodes.is_empty()
            assert len(self.partial_floatings) == 0
            assert self.need_to_visit.is_empty()
            #----
            new_log_lk = self.arg.log_likelihood(self.mu, self.data)
            new_log_prior, new_roots, new_coals = self.arg.log_prior(self.n,
                                            self.seq_length, self.r, self.Ne,
                                            False, True, True)
            #--- reverse transition
            self.transition_prob.kuhner_num_nodes(self.arg.__len__() - len(new_roots), False)
            self.Metropolis_Hastings(new_log_lk, new_log_prior, trans_prob = False, kuhner =True)
            if self.accept:
                if  self.verbose:
                    print("Kuhner MH accepted")
                self.arg.roots = new_roots
                self.arg.coal = new_coals
            else:
                self.arg = self.original_ARG
        else:
            self.arg = self.original_ARG
        self.floats.clear()
        self.need_to_visit.clear()
        self.active_nodes.clear()
        self.active_links = 0
        self.partial_floatings = collections.defaultdict(list)
        self.higher_times = collections.defaultdict(set)#time: (index)
        self.original_interval = collections.defaultdict(list)
        self.empty_containers()
        self.original_ARG = ARG()
        # gc.collect()

    #============= transition 7: update parameters
    def random_normal(self,  mean, sd):
        '''non negative random number from normal distribution'''
        v = np.random.normal(mean, sd, 1)[0]
        if v < 0:
            v = -v
        return v

    def normal_log_pdf(self, mean, sd, x):
        var = float(sd)**2
        denom = (2*math.pi*var)**0.5
        num = math.exp((-(float(x)-float(mean))**2)/(2*var))
        return math.log(num/denom)

    def update_parameters(self, mu= False, r = False, Ne= False):
        sd_mu = self.default_mutation_rate/4
        sd_r = self.default_recombination_rate/4
        sd_N  = 5
        N0 = self.default_Ne
        # propose mu
        new_mu = self.mu
        new_r= self.r
        new_N = self.Ne
        if mu:
            new_mu = self.random_normal(self.mu, sd_mu)
        if r:
            new_r = self.random_normal(self.r, sd_r)
        if Ne:
            new_N = self.random_normal(self.Ne, sd_N)
        new_log_lk = self.arg.log_likelihood(new_mu, self.data)
        new_log_prior = self.arg.log_prior(self.n,
                               self.seq_length, new_r, new_N, False)
        #------- MH
        self.accept = False
        ratio = new_log_lk + new_log_prior + \
                self.normal_log_pdf(new_r, sd_r, self.r)+\
                self.normal_log_pdf(new_mu, sd_mu, self.mu)+\
                self.normal_log_pdf(new_N, sd_N, self.Ne)+\
                self.transition_prob.log_prob_reverse - \
                (self.log_lk + self.log_prior +
                 self.transition_prob.log_prob_forward +\
                 self.normal_log_pdf(self.r, sd_r, new_r)+\
                self.normal_log_pdf(self.mu, sd_mu, new_mu)+\
                self.normal_log_pdf(self.Ne, sd_N, new_N))
        if  self.verbose:
            print("forward_prob:", self.transition_prob.log_prob_forward)
            print("reverse_prob:", self.transition_prob.log_prob_reverse)
            print("ratio:", ratio)
            print("new_log_lk", new_log_lk, "new_log_prior", new_log_prior)
            print("old.log_lk", self.log_lk,"old.log_prior", self.log_prior)

        if math.log(random.random()) <= ratio: # accept
            self.log_lk = new_log_lk
            self.log_prior = new_log_prior
            self.accept = True
            self.mu = new_mu
            self.r = new_r
            self.Ne = new_N

    def write_summary(self, row = [], write = False):
        '''
        write the ARG summary to a pd.df
        [lk, prior, numAncRec, numNonancRec, branchLength, setup]
        for setup: iteration, thining, burnin, n, seq_length, m,
            Ne, theta (mu), rho (r), acceptance rate, time took]
        '''
        if not write:
            self.summary.loc[0 if math.isnan(self.summary.index.max())\
                else self.summary.index.max() + 1] = row
        else:
            self.summary.to_hdf(self.outpath + "/summary.h5", key = "df")

    def run_transition(self, w):
        '''
        choose a transition move proportional to the given weights (w)
        Note that order is important:
        orders: "spr", "remove", "add", "at"(adjust_times), "ab" (breakpoint),
            kuhner, update__mu, update_rho, update_Ne
        NOTE: if rem and add weights are not equal, we must include
                that in the transition probabilities
        '''
        if len(self.arg.rec) == 0:
            w[4] = 0# adjust breakpoint
        ind = random.choices([i for i in range(len(w))], w)[0]
        if ind == 0:
            print("PROBLME")
        elif ind == 1:
            if len(self.arg.rec) != 0:
                print("PROBLEM22222")
        elif ind == 2:
            self.add_recombination()
        elif ind == 3:
            self.adjust_times()
        elif ind == 4:
            self.adjust_breakpoint()
        elif ind == 5:
            self.kuhner()
        elif ind == 6:
            # print("---------update_parameters")
            self.update_parameters(True, True, False)
        # elif ind == 7:
        #     print("---------update_rho")
        #     self.update_parameters(r=True)
        # elif ind == 8:
        #     print("---------update_Ne")
        #     self.update_parameters(Ne=True)
        self.detail_acceptance[ind][0] += 1
        if self.accept:
            self.detail_acceptance[ind][1] += 1

    def run(self, iteration = 20, thin= 1, burn= 0,
            verify = False):
        it = 0
        accepted = 0
        t0 = time.time()
        #---test #1 for no division by zero
        self.detail_acceptance = {i:[1, 0] for i in range(9)}
        #----test
        for it in tqdm(range(iteration), ncols=100, ascii=False):#
        # while it < iteration:
            self.run_transition(w = [0, 1, 1, 5, 1, 5, 0])#[1, 1, 1, 5, 1, 4, 0][0, 1, 1, 6, 1, 6, 0]
            if self.accept:
                accepted += 1
            if it > burn and it % thin == 0:
                self.write_summary([self.log_lk, self.log_prior,
                                    self.log_lk + self.log_prior,
                                self.arg.num_ancestral_recomb,
                                self.arg.num_nonancestral_recomb,
                                self.arg.branch_length,
                                    self.mu, self.r, self.Ne, -1])
                #dump arg
                self.print_state(self.outpath  + "/" + "output"+str(int(it))+".txt")
            if verify:
                self.arg.verify()
            it += 1
            self.accept = False
            #self.write_summary(write = True)
        # print("detail acceptance", self.detail_acceptance)

    def print_state(self, filename): ###WHAT WE ACTUALLY PRINT
        f = open(filename, "a")
        f.writelines(["node" + "\t" +  "time" + "\t"  "par"+"\n"])
        for j in self.arg.nodes:
            node = self.arg.__getitem__(j)
            if node.left_parent is not None or node.left_child is not None:
                s = node.first_segment
                if s is None:
                    f.writelines([ str(j)  + "\t" +  ("%.5f" % node.time) + "\t"  + "-1"+"\n"])
                else:
                    f.writelines([ str(j)  + "\t" +  ("%.5f" % node.time) + "\t" +str(node.left_parent.index) +"\n"])
        f.close()

def infer_sim(
        ts_full,
        sample_size,
        iteration = 100,
        thin = 20,
        burn = 0,
        Ne =5000,
        seq_length= 6e4,
        mutation_rate=1e-8,
        recombination_rate=1e-8,
        outpath = "/output",
        plot= True,
        verbose=False,
        verify= False
):
    """

    Takes `msprime`  tree sequence with `record_full_arg= True` and
    - converts tree_sequence to `Augmented Tree Sequence (ATS)` format.
    - calculates true likelihood/ prior/ branch length / the number of ancestral and non-ancestral
    recombinations and returns them at  `outpath+"/true_values.npy`
    - runs the mcmc and returns samples of ARGs from their posterior.

    :param ts_full: `msprime` tree sequence with `record_full_arg= True`.
    :param int sample_size: The number of sampled genomes.
    :param int iteration: The number of MCMC iterations. Must be `>20`. The default = `100`
    :param int thin: This specifies how often to write ARG samples to file. By default, the ARG is written every `10` iterations after `burn-in`.
    :param int burn: This specifies how many ARGs to discard as `burn-in`. Default is `0`.
    :param float Ne: The effective (diploid) population size for the population. This defaults to `5000` if not specified.
    :param float seq_length: The length of the sequences in bases.This defaults to `6e4` if not specified.
    :param float mutation_rate: The rate of mutation per base per generation. This defaults to `1e-8` if not specified.
    :param float recombination_rate: The rate of recombination per base per generation. This defaults to `1e-8` if not specified.
    :param outpath: The path to store the outputs. This defaults to `./output` if not specified.
    :param bool plot: plots the trace plots if `True`
    :param bool verbose: verbose.
    :param bool verify: for debugging purposes. Checks if the ARG is a valid ATS.
    :return: None. inferred ARGs are stored in  `outpath`.
    """
    if ts_full !=None:#else real data
        try:
            ts_full = msprime.load(ts_full.name) #trees is a fh
        except AttributeError:
            ts_full = msprime.load(ts_full)
    else:
        IOError("ts_full is required")
    # random.seed(args.random_seed)
    # np.random.seed(args.random_seed+1)
    mcmc = MCMC(ts_full= ts_full,
                sample_size= sample_size,
                Ne=Ne,
                seq_length= seq_length,
                mutation_rate= mutation_rate,
                recombination_rate=recombination_rate,
                outpath= outpath,
                verbose= verbose)
    mcmc.run(iteration = iteration, thin= thin, burn= burn,
            verify = verify)
    if plot:
        p= Trace(outpath)
        p.arginfer_trace()
    if verbose:
        mcmc.print_state()

if __name__ == "__main__":
    mcmc = MCMC()
    mcmc.run(200, 1, 0, verify=False)
