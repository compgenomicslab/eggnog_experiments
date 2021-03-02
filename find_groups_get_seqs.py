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


class OrderedCounter(Counter, OrderedDict):
     'Counter that remembers the order elements are first seen'
     def __repr__(self):
         return '%s(%r)' % (self.__class__.__name__,
                            OrderedDict(self))
     def __reduce__(self):
         return self.__class__, (OrderedDict(self),)


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def get_lca(n, lin_filter):
    lineages = OrderedCounter()
    
    nleaves = 0
    #for n in CONTENT[n]:
    for n in n.get_leaves():        
        if lin_filter in n.lineage: 
            nleaves += 1
            lineages.update(n.lineage)

    lca = [l for l, count in lineages.items() if count == nleaves]
    if not lca: 
        lca = ['Unk']
        

    return lca[-1]


def most_frequent(List): 
    return max(set(List), key = List.count) 

def load_node_scores(n):
    leaf_targets = [n.taxid for n in CONTENT[n] if TARGET_LCA in taxa[n]]
    n.nspcs = len(set(leaf_targets))   
    
    n.lca = get_lca(n, TARGET_LCA)
    
    dups_per_sp = Counter(leaf_targets)
    n.inparalogs_rate = np.median(list(dups_per_sp.values()))    
    
    n.nseqs = len(leaf_targets)    
    n.score1 = (n.nspcs/SPTOTAL) 
    n.score2 = n.nspcs/n.nseqs if n.nseqs else 0.0

def is_leaf_og(node):
    # process all leaves in the tree containing TARGET_LCA in their lineage
    # and calculates Last Common Ancestor. Other leaves not matching TARGET_LCA
    # are ignored for score computation. This allows, for instance, finding euk-OGs 
    # in a tree where there are also bacteria
    
    node_ =  get_taxo_groups2(node)

    if not hasattr(node_, 'lca'):
        load_node_scores(node_)
    
        
    if node_.is_leaf():
        # If we reached a leave, just return it as an end point. We cannot do anything        
        return True
    else:     
        # If the node is internal, let's evaluate if it could be considered an OG,
        # or should be split
        
        ch1 = node_.children[0]
        ch2 = node_.children[1]
        load_node_scores(ch1)
        load_node_scores(ch2)
                
        if TARGET_LCA not in [ch1.lca, ch2.lca]: 
            return True
        else:           
            return False


# def get_taxo_groups(node):

    # count_lin = defaultdict(lambda: defaultdict(int))
    # count_lin_mem = defaultdict(lambda: defaultdict(list))
    # count_lin_perc = defaultdict(lambda: defaultdict(float))

    # parent_node = node.up
    # lca_parent_node = get_lca(parent_node, TARGET_LCA)
    # lin_parent_node = ncbi.get_lineage(lca_parent_node)
    # filter_content = defaultdict(list)

    # for l in CONTENT[parent_node]:
        # if TARGET_LCA in l.lineage:
            # for index, (term) in enumerate(l.lineage):
                # count_lin[index][term] += 1
                # count_lin_mem[index][term].append(l.name)
                # filter_content[parent_node.name].append(l.name)
    
    # for index, term in count_lin.items():
        # if index == len(lin_parent_node):
            # for taxid, num in term.items():
                # leng = len(filter_content[parent_node.name])
                # p = float(num/leng)
                # p_taxid = float(num/counter_taxono[taxid])
                # print(p_taxid, counter_taxono[taxid], num )
                # count_lin_perc[index][taxid] = p
                # if p <= 0.001 and p_taxid < 0.01 :
                    # print(leaf.name, p, taxid)
                    # return False

    # return True

def get_taxo_groups2(node):
    #count_lin: dict para contar por cada nivel de profundidad para cada taxid el nº de seqs (ej, para Bilateria que seria nivel 6, cuantas seqs hay en ese nodo {6:{33213:x}})
    count_lin = defaultdict(lambda: defaultdict(int))
    #count_lin_mem: dict para guardar por cada nivel de profundidad para cada taxid los seqs.name 
    count_lin_mem = defaultdict(lambda: defaultdict(list))
    
    #count_lin_perc = defaultdict(lambda: defaultdict(float))

    lca_node = get_lca(node, TARGET_LCA)
    if lca_node != 'Unk':
        lin_node = ncbi.get_lineage(lca_node)
    else:
        lin_node = ['Unk']
    
    #para el nodo guarda todas las seqs, independientemente del linaje que tengan
    filter_content = defaultdict(list)

    for l in CONTENT[node]:
        if TARGET_LCA in l.lineage:
            for index, (term) in enumerate(l.lineage):
                count_lin[index][term] += 1
                count_lin_mem[index][term].append(l.name)
                filter_content[node.name].append(l.name)
    
    #para nada nivel de profundidad, para todos los grupos(taxid) en ese nivel
    for index, term in count_lin.items():
        
        #miro solo los taxids que estan un nivel por debajo del lca del nodo, es decir si el lca del nodo es bilateria solo miro protostomos, deuterostomos
        if index == len(lin_node):
            
            #para cada taxid del nivel posterior al lca del nodo calculo el porcentaje de seqs que hay de ese taxid en ese nodo
            for taxid, num in term.items():
                
                #calculo el porcentaje de seqs que hay de ese taxid (num) en ese nodo (leng)
                leng = len(filter_content[node.name])
                p = float(num/leng)
                
                #calculo el porcenje de seqs que hay de ese taxid en comparacion con todas las seqs que hay de ese taxid en todo el grupo 
                #no es lo mismo 1 protostomo si en total hay 3, que 1 protostomo si en total hay 1000
                p_taxid = float(num/counter_taxono[taxid])
                
                #count_lin_perc[index][taxid] = p

                #elimino las seqs que pertenezcan a taxid muy muy muy minoritarios
                if p <= 0.001 and p_taxid < 0.01 :
                    for le in count_lin_mem[index][taxid]:
                        det = node.search_nodes(name=le)
                        for d in det:
                            d.detach()


    remain_leaves =  []
    for l in node.get_leaves():
        remain_leaves.append(l.name)
    node.prune(remain_leaves)
    
    node_clean = node
    return node_clean

    

ncbi = NCBITaxa('/data/jhc/cold/eggnog6/build/00_level_clades/ref_trees_damian_taxonomy/etetoolkit/taxa.sqlite')
t = PhyloTree(sys.argv[1])
t.resolve_polytomy()
if sys.argv[5] == "midpoint":
    t.set_outgroup(t.get_midpoint_outgroup())

t.set_species_naming_function(lambda x: x.split('.')[0])
t.annotate_ncbi_taxa()

# Let's cache the list of leaves under each internal node. 
# This is a global variable.
CONTENT = t.get_cached_content()

total_seqs = []
for n in CONTENT[t]:
    total_seqs.append(n.name)
events = t.get_descendant_evol_events()

# Which is the target level to create OGs? This is global variable used by several funcions
TARGET_LCA = int(sys.argv[2])
lin_target = ncbi.get_lineage(TARGET_LCA)
print (TARGET_LCA)
t2n=ncbi.get_taxid_translator([TARGET_LCA])
print(t2n[TARGET_LCA])
print(lin_target)
 
        
for node in t.traverse("preorder"):
    if not node.is_leaf():
        node.name = id_generator()

    

with open('/home/plaza/projects/eggnog6/pfamA_families/eggnog6_pfamParse/domains_sorted.json') as d: 
    data = json.load(d) 

data2 = {}
for seq, dom in data.items():
    dom = dom.split(',')
    d = set(dom)
    data2[seq] = str(list(d))

data_names = {}
with open(sys.argv[6]) as names:
    data_names = json.load(names)

     

sos_dict = {}
for ev in events:
    in_seq = ev.in_seqs
    out_seq = ev.out_seqs
    all_seq = in_seq
    all_seq.update(out_seq)
    ancestor = t.get_common_ancestor(all_seq)
    sos_dict[ancestor.name] = str(ev.sos)


taxa = {}
counter_taxono = defaultdict(int)
for leaf in t:
    taxa[leaf] = set(leaf.lineage)
    for taxid in leaf.lineage:
        counter_taxono[taxid] +=1






# What is the LCA of the whole tree? (just for information)
e_tree_lca = get_lca(t, 'root')

# number of species in the tree mathching the requested LCA. For instance, 
# if Eukaryota is requested, all bacteria are ignored for the computations. 
# Note that SPTOTAL is a global variable used by several functions like get_node_score. 
SPTOTAL = len(set([n.taxid for n, lin in taxa.items() if TARGET_LCA in lin]))
print(SPTOTAL, len(taxa), e_tree_lca)



# finds the nodes that can be collpased as OGs given the TARGET_LCA. 
# The decision is made by the is_leaf_of() function, which calculates several scores per node
scores = []

og_doms_final = defaultdict(list)
count_tax_levels = defaultdict(dict)



for leaf in t.iter_leaves(is_leaf_fn=is_leaf_og):

    scores.append([leaf.score1, leaf.score2, len(CONTENT[leaf]), leaf.nspcs, leaf.lca, leaf])
    leaf.img_style['draw_descendants'] = False
    og_doms_prefilter = defaultdict(list)
    level_len = defaultdict(list)
    d_list = []

    if leaf.lca != 'Unk':
        tax2name = ncbi.get_taxid_translator([leaf.lca])
        name_lca = tax2name[leaf.lca]
    else:
        name_lca = leaf.lca

    if leaf.name in sos_dict.keys():
        so_node = float(sos_dict[leaf.name])
    else:
        so_node = 0.0

    if leaf.up.name in sos_dict.keys():
        so_parent = float(sos_dict[leaf.up.name])
    else:
        so_parent = 0.0

    sp_nodes = []
    if so_parent > 0.3:
        for l in CONTENT[leaf]:
            sp_ = l.name.split('.')[0]
            if TARGET_LCA in ncbi.get_lineage(sp_) :
                og_doms_prefilter[leaf.name].append(l.name)
                level_len[leaf.name].append(l.name)
                sp_nodes.append(sp_)
    else:
        leaf.img_style["bgcolor"] = "Salmon"
    
    # if len(sp_nodes) >0:
        # top = ncbi.get_topology(sp_nodes)
        #print(top)
    #-------------------------------------------------------
    # if so_node <= 0.3:
        # for l in CONTENT[leaf]:
            # sp_ = l.name.split('.')[0]
            # if TARGET_LCA in ncbi.get_lineage(sp_):
                # og_doms_prefilter[leaf.name].append(l.name)
                # level_len[leaf.name].append(l.name)
    # elif so_node >= 0.3 and leaf.lca != 'Unk':
        # leaf.img_style["bgcolor"] = "Salmon"
        # for l in CONTENT[leaf]:
            # sp_ = l.name.split('.')[0]
            # if TARGET_LCA in ncbi.get_lineage(sp_):
                # level_len[leaf.name].append(l.name)
    #---------------------------------------------------------        
    for node, leafs in og_doms_prefilter.items():
        for l in leafs:
            d = data2[l]
            d_list.append(d)
    count_doms = (Counter(d_list))
    common = count_doms.most_common(1)
    gnamesl = list()
    
    if len(d_list) != 0:
        for node, leafs in og_doms_prefilter.items():
            for l in leafs:
                if common[0][0] == data2[l]:
                    og_doms_final[leaf.name].append(l)
                    if l in data_names.keys():
                        gnamesl.append(data_names[l])
            if len(og_doms_final[leaf.name]) > 3:
                leaf.img_style["bgcolor"] = "LightGreen"
            else:
                leaf.img_style["bgcolor"] = "Orange"
            

    count_names = Counter(gnamesl)
    common_name = count_names.most_common(1)

    leaf.img_style['size'] = len(og_doms_final[leaf.name])
    leaf.add_face(TextFace(("LCA:%s; DOMS:%s; GENE:%s ;LEN:%s; LEVEL_LEN:%s; TOTAL_LEN:%s" % (name_lca, common, common_name, len(og_doms_final[leaf.name]), len(level_len[leaf.name]),len(CONTENT[leaf])))), column=0, position="branch-right")



scores.sort(reverse=True, key=lambda x: x[:-1])
covered_seqs = 0
for s in scores[:10]:    
    covered_seqs += s[2]
    print(s[:5], "%0.2f" %covered_seqs, "%0.2f" %(covered_seqs/len(taxa)), sep='\t')
   
def layout(node):
    if not node.is_leaf():
        node.add_face(TextFace(node.name), column=0, position = "branch-right")
        node.add_face(TextFace(sos_dict[node.name]), column=0, position = "branch-right")

        sp_l = []
        for l in CONTENT[node]:
            sp_l.append(l.name.split('.')[0])
        t_anc = ncbi.get_topology(sp_l)
        name_anc = t_anc.sci_name

        node.add_face(TextFace(name_anc), column=0, position = "branch-right")   
    

    if getattr(node, 'evoltype', None) == 'S':
        node.img_style["fgcolor"] = 'red'
        
    elif getattr(node, 'evoltype', None) == 'D':
        node.img_style["fgcolor"] = 'green'
    
        

ts = TreeStyle()
ts.layout_fn = [layout]
ts.show_leaf_name = False

name_tree = os.path.basename(sys.argv[1])
path_out = sys.argv[4]  

outfile =path_out+name_tree+'_'+str(TARGET_LCA)+'_'+'doms_14'+'.pdf'
print(len(t))
t.render( outfile, w=350, units="mm", tree_style=ts)

seqs = SeqGroup(sys.argv[3])
seqs_dict = {}

for n, (name,seq,_) in enumerate(seqs):
    seqs_dict[name] = seq

seqs_in_og = []
# for name_node,seqs in og_doms_final.items():
    # if len(seqs) >3:
        # out_name = name_node+'_'+name_tree+'_'+str(TARGET_LCA)
        # with open(path_out+out_name+'.faa', 'a') as f_out:
            # for s in seqs:
                # seqs_in_og.append(s)
                # aa = seqs_dict[s]
                # f_out.write('>'+s+'\n'+aa+'\n')

# with open(path_out+'not_og-'+name_tree+'.faa', 'a') as f_out:
    # for s in total_seqs:
        # if s not in seqs_in_og:
            # sp_ = s.split('.')[0]
            # if TARGET_LCA in ncbi.get_lineage(sp_):
                # aa = seqs_dict[s]
                # f_out.write('>'+s+'\n'+aa+'\n')
