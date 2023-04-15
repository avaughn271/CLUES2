''' This module is responsible for the mcmc '''
import treeSequence
import copy
from tqdm import tqdm
import sys
import os
import shutil
from initialARG import *
sys.setrecursionlimit(500000)

class TransProb(object):
    '''transition probability calculation'''
    def __init__(self):
        self.log_prob_forward = 0
        self.log_prob_reverse = 0

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

    def add_choose_child(self):
        '''choose a node to put a recombination on it'''
        ind = random.choice(list(self.arg.nodes.keys()))
        if not self.arg.roots.__contains__(ind):
            return self.arg.__getitem__(ind)
        else:
            return self.add_choose_child()

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

    def new_recombination(self, t):
        print("HUGEPROBLEM")

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
            print("ttttttttttt")
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
                        print("big old problem")
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
            print("PROBLEM2")
        elif ind == 1:
            print("problem4")
        elif ind == 2:
            print("problem6")
        elif ind == 3:
            self.adjust_times()
        elif ind == 4:
            print("problem3")
        elif ind == 5:
            self.kuhner()
        elif ind == 6:
            print("PROBLEM")
        self.detail_acceptance[ind][0] += 1
        if self.accept:
            self.detail_acceptance[ind][1] += 1

    def run(self, iteration = 20, thin= 1, burn= 0):
        it = 0
        accepted = 0
        t0 = time.time()
        #---test #1 for no division by zero
        self.detail_acceptance = {i:[1, 0] for i in range(9)}
        #----test
        for it in tqdm(range(iteration), ncols=100, ascii=False):#
        # while it < iteration:
            self.run_transition(w = [0, 0, 0, 5, 0, 5, 0])
            if self.accept:
                accepted += 1
            if it > burn and it % thin == 0:
                #dump arg
                self.arg.dump(self.outpath,file_name="arg"+str(int(it))+".arg")
            it += 1
            self.accept = False
        if iteration > 18:
            self.summary.setup[0:18] = [iteration, thin, burn, self.n, self.seq_length,
                                       self.m, self.Ne, self.mu, self.r,
                                        round(accepted/iteration, 2), round((time.time() - t0), 2),
                                        round(self.detail_acceptance[0][1]/self.detail_acceptance[0][0], 2),
                                        round(self.detail_acceptance[1][1]/self.detail_acceptance[1][0], 2),
                                        round(self.detail_acceptance[2][1]/self.detail_acceptance[2][0], 2),
                                        round(self.detail_acceptance[3][1]/self.detail_acceptance[3][0], 2),
                                        round(self.detail_acceptance[4][1]/self.detail_acceptance[4][0], 2),
                                        round(self.detail_acceptance[5][1]/self.detail_acceptance[5][0], 2),
                                        round(self.detail_acceptance[6][1]/self.detail_acceptance[6][0], 2)]
        # print("detail acceptance", self.detail_acceptance)

    def print_state(self):
        print("self.arg.coal", self.arg.coal)
        print("self.arg.rec", self.arg.rec)
        print("self.arg.roots", self.arg.roots)
        print("self.floatins", self.floatings)
        print("node", "time", "left", "right", "l_chi", "r_chi", "l_par", "r_par",
              "l_bp", "snps", "fir_seg_sam",
              sep="\t")
        for j in self.arg.nodes:
            node = self.arg.__getitem__(j)
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

if __name__ == "__main__":
    mcmc = MCMC()
    mcmc.run(200, 1, 0)
