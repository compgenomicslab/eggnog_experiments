
from ete3 import PhyloTree, NCBITaxa, TreeStyle, TextFace, NodeStyle
from ete3 import SeqGroup
import sys
from collections import Counter, OrderedDict, defaultdict
import numpy as np
import os
import string
import random
import json
import pandas as pd
from ete3.treeview import random_color
from ete3 import ProfileFace ,ImgFace, RectFace
from ete3.treeview.faces import add_face_to_node
import statistics 
import re
import time
#from layout_pfam import layout
import pickle 

start_time = time.time()

ncbi = NCBITaxa('/data/jhc/cold/eggnog6/build/00_level_clades/ref_trees_damian_taxonomy/etetoolkit/taxa.sqlite')

class OrderedCounter(Counter, OrderedDict):
     'Counter that remembers the order elements are first seen'
     def __repr__(self):
         return '%s(%r)' % (self.__class__.__name__,
                            OrderedDict(self))
     def __reduce__(self):
         return self.__class__, (OrderedDict(self),)


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def most_frequent(List): 
    return max(set(List), key = List.count) 

def get_lca(n):
    lineages = OrderedCounter()
    nleaves = 0
    for n in n.get_leaves():        
        nleaves += 1
        lineages.update(n.lineage)

    lca = [l for l, count in lineages.items() if count == nleaves]
    if not lca: 
        lca = ['Unk']
    return lca[-1]


def load_node_scores(n):
    leaf_targets = []
    for l in n.leaves_in:
        leaf_targets.append(l.split('.')[0])
    n.nspcs = len(set(leaf_targets))   
    
    
    dups_per_sp = Counter(leaf_targets)
    n.inparalogs_rate = np.median(list(dups_per_sp.values()))    
    
    n.nseqs = len(n.leaves_in)    
    n.score1 = (n.nspcs/SPTOTAL) 
    n.score2 = n.nspcs/n.nseqs if n.nseqs else 0.0


def is_leaf_og(node):
    # process all leaves in the tree containing TARGET_LCA in their lineage
    # and calculates Last Common Ancestor. Other leaves not matching TARGET_LCA
    # are ignored for score computation. This allows, for instance, finding euk-OGs 
    # in a tree where there are also bacteria


    # if not hasattr(node, 'lca'):
        # load_node_scores(node)
    if "rank" in node.features:
        rank =  node.rank
    else:                
        rank = 'not rank'   

    if node.is_leaf():
        # If we reached a leave, just return it as an end point. We cannot do anything        
        return True
    else:     
        # If the node is internal, let's evaluate if it could be considered an OG,
        # or should be split
                
        if rank in ['species', 'strain', 'varietas', 'genus', 'subspecies']: 
            return True
        else:           
            return False

#FUNCIONES PARA CALCULAR LAS PERDIDAS
def get_dup_score(n, n2content):
    
    if len(n.sp_in_ch1) >1 and len(n.sp_in_ch2) >1:
        sp1 = n.sp_in_ch1
        sp2 = n.sp_in_ch2

        dup_score = len(sp1 & sp2) / min(len(sp1), len(sp2))
    else: 
        dup_score = 0
        sp1 = None
        sp2 = None
    return dup_score, sp1, sp2

def process_tree(node):

    #n2content = node.get_cached_content()   
    dup_score, sp1, sp2 = get_dup_score(node, node.sp_in)
    
    if dup_score > 0.0:    
        losses1 = count_losses(expected_sp=sp1|sp2, found_sp=sp1)
        losses2= count_losses(expected_sp=sp1|sp2, found_sp=sp2)

        return dup_score, losses1, losses2
    else:
        return None, None, None
        #print("  Dup score %0.6f, losses branchA: %s, losses branchB %s" %(dup_score, losses1, losses2))


def count_losses(expected_sp, found_sp):
    def is_leaf_2(_n):       
        # Returns true if node is an actual leaf, or a branch without target species
        if not _n.children: 
            return True
        elif not (refn2sp[node] & found_sp):
            return True
        else: 
            return False

    if len(expected_sp) == 1: 
        return 0
    root = reftree.get_common_ancestor(expected_sp)
    losses = 0
    n2losses = {}
    for post, node in reftree.iter_prepostorder():
        if post: 
            n2losses[node] = sum([n2losses[ch] for ch in node.children])
        else: 
            if not node.children:
                if node.name not in found_sp:
                    n2losses[node] = 1
                else: 
                    n2losses[node] = 0
    miss_sp = set()
    for node in root.traverse(is_leaf_fn=is_leaf_2): 
        if not (refn2sp[node] & found_sp):
            losses += 1
    return losses


#FUNCIONES PARA ENCONTRAR OUTLIERS

def proces_node(n):
    #count_lin: dict para contar por cada nivel de profundidad para cada taxid el nº de seqs (ej, para Bilateria que seria nivel 6, cuantas seqs hay en ese nodo {6:{33213:x}})
    count_lin = defaultdict(lambda: defaultdict(int))

    #count_lin_mem: dict para guardar por cada nivel de profundidad para cada taxid los seqs.name 
    count_lin_mem = defaultdict(lambda: defaultdict(list))

    #para el nodo guarda todas las seqs, independientemente del linaje que tengan
    sp_in_node = set()
    sp_per_level = defaultdict(set)

    leaves = CONTENT[n]
    leaves_names  = []
    for l in leaves:
        leaves_names.append(l.name)
        #las sp solo las tengo en cuenta la primera vez que aparecen (me da igual que haya varias seqs q la misma sp)
        if l.taxid not in sp_in_node: #and l.taxid not in n.sp_out_up:
            sp_in_node.add(l.taxid)
            for index, tax in enumerate(l.lineage):
                count_lin[index][tax] += 1
                count_lin_mem[index][tax].append(l.name)
                sp_per_level[tax].add(l.taxid)

    sp2remove = set()

    if n.up and properties[n.up.name] and properties[n.up.name]['sp_out']:
        for sp in properties[n.up.name]['sp_out']:
            if sp in sp_in_node:
                sp2remove.add(sp)

    for level, linages in count_lin.items():
            
        #si para un mismo nivel hay mas de un linaje, ej nivel 0 son todas root pero en nivel 1 hay euk y bact
        if len(linages.keys())>1:
                
            for tax, num in linages.items():
                #num = numero de sp para ese nivel taxonomico en el nodo
                #sp_in_node = todas las sp que hay en el nodo
                #global_linages = numero de sp para ese nivel taxonomico en todo el arbol
                #SPTOTAL = numero de sp en total del arbol
                #levels_eggnog = numero de sp para ese nivel taxonomico en eggnog

                #% linaje en el nodo
                per_Node = num / len(sp_in_node)

                #% linaje  el arbol
                per_Tree = num / global_linages[level][tax]
                per_Tree_total = num /SPTOTAL
                
                #% linaje en eggnog
                per_Egg = num / levels_eggnog[str(tax)]
                per_Egg_total = levels_eggnog[str(tax)] / len(reftree)
                
                
                if per_Egg <0.05:
                    if per_Node < 0.01 and per_Tree_total< 0.01:
                        for seq in count_lin_mem[level][tax]:
                            sp2remove.add(str(seq.split('.')[0]))
                           
    
    return sp2remove            
    



print('-load taxonomy')
levels_eggnog = {}
with open('/home/plaza/projects/eggnog6/pfamA_families/eggnog_experiments/data/levels2numSp.json') as levels:
    levels_eggnog = json.load(levels)


print('load reftree')
reftree = PhyloTree('/home/plaza/projects/eggnog6/pfamA_families/eggnog_experiments/data/totalNCBITree.nw')
refn2sp = reftree.get_cached_content(store_attr='name', container_type=set)
print(len(reftree))


print('READ TREE, SEQ FILE, REF TREE')
t = PhyloTree(sys.argv[1])
name_tree = os.path.basename(sys.argv[1])
fasta = SeqGroup(sys.argv[2])
path_out = sys.argv[3]


print('START PROCESSING TREE')
t.resolve_polytomy()
if  len(sys.argv)>4  and sys.argv[4] == "midpoint":
    print('MIDPOINT ROOTING')
    t.set_outgroup(t.get_midpoint_outgroup())

t.set_species_naming_function(lambda x: x.split('.')[0])
ncbi.annotate_tree(t,  taxid_attr="species")

ref_leaves_names = []
for l in reftree.get_leaves():
    ref_leaves_names.append(l.name)

get_ref_ancestor = ncbi.get_topology(ref_leaves_names)


#count how many species are in each taxonomic level
taxa = {}
sp_set = set()
counter_taxono = defaultdict(int)
global_linages = defaultdict(lambda: defaultdict(int))

total_mems_in_tree = set()


miss_lineages = []
for leaf in t:
    total_mems_in_tree.add(leaf.name)

    if len(leaf.lineage) == 0:
        miss_lineages.append(leaf.name)
    
    if leaf.taxid not in sp_set:
        sp_set.add(leaf.taxid)
        taxa[leaf] = set(leaf.lineage)
        for index, term in enumerate(leaf.lineage):
            global_linages[index][term] += 1
            counter_taxono[term] +=1

print('MISS LING: ', len(miss_lineages))
SPTOTAL = len(sp_set)

# What is the LCA of the whole tree? (just for information)
e_tree_lca = get_lca(t)


t.get_descendant_evol_events() 

print('DETECT OUTLIERS AND DUPS SCORE FUNCTION')
dup_score_d = defaultdict(list)
seq2info = defaultdict(list)
pnames = set()
ogs = set()
pfams = set()
properties = defaultdict(dict)

for n in t.traverse("preorder"):
    if not n.is_leaf():
        n.name = id_generator()

# Let's cache the list of leaves under each internal node. 
# This is a global variable.
CONTENT = t.get_cached_content()

for n in t.traverse("preorder"):
    if n.is_leaf():
        sci_name = ncbi.get_taxid_translator([n.taxid])[n.taxid]
        properties[n.name]['sci_name'] = sci_name

    if not n.is_leaf():
        sp_out = []
        leaves_out = set()
        sp_in = set()
        leaves_in = set()
        

        if n.is_root():
            all_mems = []
            for l in n.get_leaves():
                all_mems.append(l.name)

            properties[n.name]['is_root'] = 'True'
            properties[n.name]['all_mems'] = all_mems

        #only detect outliers if node is a duplication node
        if getattr(n, 'evoltype', None) == 'D':
            sp_out = proces_node(n)
            if len(sp_out) > 0:
                all_leafs = n.get_leaves()

                for l in all_leafs:
                    if str(l.taxid) in sp_out:
                        leaves_out.add(l.name)   
                    else:
                        sp_in.add(str(l.taxid))
                        leaves_in.add(l.name)

            else:
                for l in  n.get_leaves():
                    sp_in.add(str(l.taxid))
                    leaves_in.add(l.name)
        else:
            for l in  n.get_leaves():
                sp_in.add(str(l.taxid))
                leaves_in.add(l.name)
    

        all_spcs = set()
        ch1 = n.children[0]
        ch1_name = n.children[0].name
        sp_ch1 = set()
        leaves_ch1 = set()
        ch2 = n.children[1]
        ch2_name = n.children[1].name
        sp_ch2 = set()
        leaves_ch2 = set()

        for l in ch1:
            if str(l.taxid) not in sp_out:
                all_spcs.add(str(l.taxid))
                sp_ch1.add(str(l.taxid))
                leaves_ch1.add(l.name)
        for l in ch2:
            if str(l.taxid) not in sp_out:
                all_spcs.add(str(l.taxid))
                sp_ch2.add(str(l.taxid))
                leaves_ch2.add(l.name)

        overlaped_spces = set(sp_ch1 & sp_ch2)
        if len(overlaped_spces)>0:
            so_score = float(len(overlaped_spces) / len(all_spcs))
        else:
            so_score = 0.0
        

        properties[n.name]['bootstrap'] = n.support



        lca_node = ncbi.get_topology(list(sp_in)).taxid
        rank = ncbi.get_rank([lca_node])[lca_node]
        lin_lca = ncbi.get_lineage(lca_node)
        
        if n.up:
            lca_node_up = n.up.lca_node
            lin_lca_node_up = ncbi.get_lineage(lca_node_up)
        else:
            lin_lca_node_up = []

    
        if len(lin_lca) < len(lin_lca_node_up):
            n.add_feature('lca_node', lca_node_up)
            properties[n.name]['lca_node'] = lca_node_up
            lca_node_name = ncbi.get_taxid_translator([lca_node_up])[lca_node_up]
            
        else:
            n.add_feature('lca_node', lca_node)
            properties[n.name]['lca_node'] = lca_node
            lca_node_name = ncbi.get_taxid_translator([lca_node])[lca_node]
        
        #n.add_feature('lca_node_name', lca_node_name)
        properties[n.name]['lca_node_name'] = lca_node_name
        
        #n.add_feature('rank', rank)
        properties[n.name]['rank'] = rank
        
        n.add_feature('sp_in', sp_in)
        properties[n.name]['sp_in'] = list(sp_in)
        
        n.add_feature('sp_in_ch1', sp_ch1)
        properties[n.name]['sp_in_ch1'] = list(sp_ch1)
        
        n.add_feature('sp_in_ch2', sp_ch2)
        properties[n.name]['sp_in_ch2'] = list(sp_ch2)

        #n.add_feature('ch1_name', ch1_name)
        properties[n.name]['ch1_name'] = ch1_name

        #n.add_feature('ch2_name', ch2_name)
        properties[n.name]['ch2_name'] = ch2_name

        n.add_feature('leaves_in', leaves_in)
        properties[n.name]['leaves_in'] = list(leaves_in)

        #n.add_feature('leaves_ch1', leaves_ch1)
        properties[n.name]['leaves_ch1'] = list(leaves_ch1)

        #n.add_feature('leaves_ch2', leaves_ch2)
        properties[n.name]['leaves_ch2'] = list(leaves_ch2)

        #n.add_feature('sp_out', sp_out)
        properties[n.name]['sp_out'] = list(sp_out)

        #n.add_feature('leaves_out', leaves_out)
        properties[n.name]['leaves_out'] = list(leaves_out)

        #n.add_feature('so_score', so_score)
        if so_score > 0.0:
            n.add_feature('evoltype_2', 'D')
            properties[n.name]['evoltype_2'] = 'D'
        else:
            n.add_feature('evoltype_2', 'S')
            properties[n.name]['evoltype_2'] = 'S'

        properties[n.name]['so_score'] = so_score

        load_node_scores(n)
        #n.add_feature('score1', n.score1)
        properties[n.name]['score1'] = n.score1

        #n.add_feature('score2', n.score2)
        properties[n.name]['score2'] = n.score2

        #n.add_feature('inparalogs_rate', n.inparalogs_rate) 
        properties[n.name]['inparalogs_rate'] = n.inparalogs_rate

        dup_score, loss1, loss2 = process_tree(n)
        if dup_score != None:
            properties[n.name]['dup_score'] = [dup_score, loss1, loss2]
            #dup_score_d[n.name] = [dup_score, loss1, loss2]
        else:
            properties[n.name]['dup_score'] = []

        
       


og_dict = defaultdict(dict)
tax_dups_og = set()
explore_dups = set()
total_mems_in_ogs = set()



print('ITER DUPS')
for node in t.traverse("preorder"):
    if node.is_root():
        
        og_dict[node.name]['mems'] = properties[node.name]['all_mems']
        og_dict[node.name]['lca'] = properties[node.name]['lca_node']
        og_dict[node.name]['dup_lineage'] = list(ncbi.get_lineage(properties[node.name]['lca_node']))
        og_dict[node.name]['so_score'] = properties[node.name]['so_score']
        og_dict[node.name]['is_root'] = 'True'
        og_dict[node.name]['dup_node_name'] = node.name
                        
   
    if not node.is_leaf() and properties[node.name]['evoltype_2'] == 'D' and 0.2 < float(properties[node.name]['so_score']) \
    and len(properties[node.name]['leaves_in'])>2 and len(properties[node.name]['sp_in'])> 2 :
    #and float( len(properties[node.name]['sp_in']) / counter_taxono[(node.lca_node)] ) > 0.6:
        
        dups_under_node = []
        for n in node.search_nodes(evoltype_2='D'):
            if n.name!= node.name:
                dups_under_node.append(n)
        
        
        #Hay mas dups bajo ese nodo
        if len(dups_under_node) > 0:
                
            lca_target = node.lca_node

            dups_under_ch1 = node.children[0].search_nodes(lca_node=lca_target)
            dups_under_ch2 = node.children[1].search_nodes(lca_node=lca_target)
            save_dups_ch1 = defaultdict()
            save_dups_ch2 = defaultdict()
            
            if len(dups_under_ch1)> 0:
                #dups por debajo del nodo principal anotadas al mismo nivel
                for n_ in  dups_under_ch1:
                    if float(properties[n_.name]['so_score'])> 0.2 and lca_target == properties[n_.name]['lca_node'] \
                    and len(properties[n_.name]['leaves_ch1'])>2 and len(properties[n_.name]['sp_in_ch1'])> 2 :
                    #and float( len(properties[n_.name]['sp_in']) / counter_taxono[(lca_target)] ) > 0.6:

                        root2node = node.get_distance(node,n_, topology_only=True)
                        save_dups_ch1[n_.name] = root2node
            
            if len(dups_under_ch2)> 0:
                for n_ in  dups_under_ch2:
                    if float(properties[n_.name]['so_score'])> 0.2  and lca_target == properties[n_.name]['lca_node'] \
                    and len(properties[n_.name]['leaves_ch2'])>2 and len(properties[n_.name]['sp_in_ch2'])> 2 :
                    #and float( len(properties[n_.name]['sp_in']) / counter_taxono[(lca_target)] ) > 0.6:

                        root2node = node.get_distance(node,n_, topology_only=True)
                        save_dups_ch2[n_.name] = root2node
                           
            # sort_dups_1 = {k: v for k, v in sorted(save_dups_ch1.items(), key=lambda item: item[1] ,reverse = True)}
            # sort_dups_2 = {k: v for k, v in sorted(save_dups_ch2.items(), key=lambda item: item[1] ,reverse = True)}


            if len(save_dups_ch1) == 0:  
                #si las dups que hay por debajo del hijo1 no cumplen los requisitos, entonces hijo1 es OG
                #mems_in_og = set()
                og_ch1_mems = set(properties[node.name]['leaves_ch1'])
                og_name_ch1 = properties[node.name]['ch1_name']
                sp_ch1 = properties[node.name]['sp_in_ch1']
                
                #if  len(og_ch1_mems.intersection(mems_in_og)) == 0 :
                if len(sp_ch1) > 2 and len(og_ch1_mems) > 2:
                    og_dict[og_name_ch1]['mems'] = list((og_ch1_mems))
                    og_dict[og_name_ch1]['lca'] = properties[node.name]['lca_node']
                    og_dict[og_name_ch1]['dup_lineage'] = list(ncbi.get_lineage(properties[node.name]['lca_node']))
                    og_dict[og_name_ch1]['so_score'] = properties[node.name]['so_score']
                    og_dict[og_name_ch1]['dup_node_name'] = node.name
                    properties[og_name_ch1]['node_is_og'] = 'True'
                    
                    #mems_in_og.update(og_ch1_mems)
                    tax_dups_og.add(node.lca_node)
                    total_mems_in_ogs.update(og_ch1_mems)
                
            
                
            if len(save_dups_ch2) == 0 :  
                #si las dups que hay por debajo del hijo2 no cumplen los requisitos, entonces hijo2 es OG
                #mems_in_og = set()
                og_ch2_mems = set(properties[node.name]['leaves_ch2'])
                og_name_ch2 = properties[node.name]['ch2_name']
                sp_ch2 = properties[node.name]['sp_in_ch2']
                
                #if  len(og_ch2_mems.intersection(mems_in_og)) == 0 :
                if len(sp_ch2) > 2 and len(og_ch2_mems) > 2:   
                    og_dict[og_name_ch2]['mems'] = list((og_ch2_mems))
                    og_dict[og_name_ch2]['lca'] = properties[node.name]['lca_node']
                    og_dict[og_name_ch2]['dup_lineage'] = list(ncbi.get_lineage(properties[node.name]['lca_node']))
                    og_dict[og_name_ch2]['so_score'] = properties[node.name]['so_score']
                    og_dict[og_name_ch2]['dup_node_name'] = node.name
                    properties[og_name_ch2]['node_is_og'] = 'True'
                    
                    #mems_in_og.update(og_ch2_mems)
                    tax_dups_og.add(node.lca_node)
                    total_mems_in_ogs.update(og_ch2_mems)
                
                
        #No hay mas dups   
        elif len(dups_under_node) == 0:

            #mems_in_og = set()

            og_ch1_mems = set(properties[node.name]['leaves_ch1'])
            og_name_ch1 = properties[node.name]['ch1_name']
            sp_ch1 = properties[node.name]['sp_in_ch1']
            og_ch2_mems = set(properties[node.name]['leaves_ch2'])
            og_name_ch2 = properties[node.name]['ch2_name']
            sp_ch2 = properties[node.name]['sp_in_ch2']
            
            
            og_dict[og_name_ch1]['mems'] = list((og_ch1_mems))
            og_dict[og_name_ch1]['lca'] = properties[node.name]['lca_node']
            og_dict[og_name_ch1]['dup_lineage'] = list(ncbi.get_lineage(properties[node.name]['lca_node']))
            og_dict[og_name_ch1]['so_score'] = properties[node.name]['so_score']
            og_dict[og_name_ch1]['dup_node_name'] = node.name
            properties[og_name_ch1]['node_is_og'] = 'True'
            
            #mems_in_og.update(og_ch1_mems)
            tax_dups_og.add(node.lca_node)
            total_mems_in_ogs.update(og_ch1_mems)
            
            og_dict[og_name_ch2]['mems'] = list((og_ch2_mems))
            og_dict[og_name_ch2]['lca'] = properties[node.name]['lca_node']
            og_dict[og_name_ch2]['dup_lineage'] = list(ncbi.get_lineage(properties[node.name]['lca_node']))
            og_dict[og_name_ch2]['so_score'] = properties[node.name]['so_score']
            og_dict[og_name_ch2]['dup_node_name'] = node.name
            properties[og_name_ch2]['node_is_og'] = 'True'
            
            #mems_in_og.update(og_ch2_mems)
            tax_dups_og.add(node.lca_node)
            total_mems_in_ogs.update(og_ch2_mems)
                



print('ogs tree')
for name_og, info in og_dict.items():
    mems_og = set(info['mems'])
    og_dict[name_og]['anc_og'] = []
    node_intersection = []
    max_dist = t.get_distance(t,name_og, topology_only=True)
    #save all og with shared members
    for n, i in og_dict.items():
        og2compare = set(i['mems'])
        if mems_og.intersection(og2compare):
            node_intersection.append(n)

    #get distance from node with shared memers to root
    save_dups = defaultdict()
    for node in node_intersection:
        root2node = t.get_distance(t,node, topology_only=True)
        if root2node < max_dist:
            save_dups[node] = root2node
    sort_dups = {k: v for k, v in sorted(save_dups.items(), key=lambda item: item[1] ,reverse = True)}
    prev_og = defaultdict(dict)
    if len(sort_dups) > 0:
        
        for og_anc, dist in sort_dups.items():
            lca = properties[og_anc]['lca_node']
            prev_og[og_anc]['dist'] = dist
            prev_og[og_anc]['lca'] = lca
        og_dict[name_og]['anc_og'] = prev_og
    else:
        
        root_node = t.get_tree_root().name
        root_lca = properties[root_node]['lca_node']
        prev_og[root_node]['dist'] = 0.0
        prev_og[root_node]['lca'] = root_lca
        og_dict[root_node]['anc_og'] = prev_og
        



# for name_og, info in og_dict.items():
    # print(name_og, og_dict[name_og]['anc_og'])
    

                

for leaf in t.iter_leaves(is_leaf_fn=is_leaf_og):
    properties[leaf.name]['draw_descendants'] = 'False'    
    

diff = (total_mems_in_tree.difference(total_mems_in_ogs))

name_fam = name_tree.split('.')[0]
not_og_fasta = path_out+name_fam+'_not_og.faa'

with open(not_og_fasta, 'w') as  f_out:
    for name_seq in diff:
        aa = fasta.get_seq(name_seq)
        f_out.write('>'+name_seq+'\n'+aa+'\n')

for name_og, info in og_dict.items():
    lca = str(info['lca'])
    with open(path_out+name_og+'_'+lca+'.faa', 'w') as f_out:
        
        for m in info["mems"]:
            aa = fasta.get_seq(m)
            f_out.write('>'+m+'\n'+aa+'\n')


post_tree = path_out+'post_'+name_fam+'.nw'
t.write(format=1, outfile=post_tree)



# for name, info in properties.items():
    # if 'node_is_og' in properties.keys():
        # properties[name]['']

with open(path_out+'properties.json', 'w') as fp:
    json.dump(properties, fp,  indent=1)



# not_og_mems_out = path_out+name_tree+'_not_og_seqs.json'
# with open(not_og_mems_out, 'w') as f:
    # json.dump(seqs_out_og, f)


og_mems_out = path_out+name_tree+'_mems.json'
with open(og_mems_out, 'w') as f:
    json.dump(og_dict, f, indent=2)

print("--- %s min---" % ((time.time() - start_time)/60))