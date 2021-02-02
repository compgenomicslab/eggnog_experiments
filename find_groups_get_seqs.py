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
    for n in CONTENT[n]:        
        if lin_filter in n.named_lineage: 
            nleaves += 1
            lineages.update(n.named_lineage)

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
    if not hasattr(node, 'lca'):
        load_node_scores(node)
    
        
    if node.is_leaf():
        # If we reached a leave, just return it as an end point. We cannot do anything        
        return True
    else:     
        # If the node is internal, let's evaluate if it could be considered an OG,
        # or should be split                        
        ch1 = node.children[0]
        ch2 = node.children[1]
        load_node_scores(ch1)
        load_node_scores(ch2)
                
        if TARGET_LCA not in [ch1.lca, ch2.lca]: 
            return True
        else:
            return False


ncbi = NCBITaxa('/data/jhc/cold/eggnog6/build/00_level_clades/ref_trees_damian_taxonomy/etetoolkit/taxa.sqlite')
t = PhyloTree(sys.argv[1])

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
TARGET_LCA = sys.argv[2]
 
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
     

sos_dict = {}
for ev in events:
    in_seq = ev.in_seqs
    out_seq = ev.out_seqs
    all_seq = in_seq
    all_seq.update(out_seq)
    ancestor = t.get_common_ancestor(all_seq)
    sos_dict[ancestor.name] = str(ev.sos)


taxa = {}
for leaf in t:
    taxa[leaf] = set(leaf.named_lineage)


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
og_doms = defaultdict(list)

for leaf in t.iter_leaves(is_leaf_fn=is_leaf_og):
    scores.append([leaf.score1, leaf.score2, len(CONTENT[leaf]), leaf.nspcs, leaf.lca, leaf])
    leaf.img_style['draw_descendants'] = False
    leaf.img_style['size'] = leaf.nspcs

    d_list = []
    for l in CONTENT[leaf]:
        d = data2[l.name]
        d_list.append(d)

    count_doms = (Counter(d_list))
    common = count_doms.most_common(1)
    
    for l in CONTENT[leaf]:
        if common[0][0] == data2[l.name]:
            og_doms[leaf.name].append(l.name)


    leaf.add_face(TextFace(("LCA:%s; DOMS:%s; LEN:%s;" % (leaf.lca, common, len(CONTENT[leaf]), ))), column=0, position="branch-right")

    
scores.sort(reverse=True, key=lambda x: x[:-1])
covered_seqs = 0
for s in scores[:10]:    
    covered_seqs += s[2]
    print(s[:5], "%0.2f" %covered_seqs, "%0.2f" %(covered_seqs/len(taxa)), sep='\t')
   
def layout(node):
    if not node.is_leaf():
        node.add_face(TextFace(node.name), column=0, position = "branch-right")
    

    if getattr(node, 'evoltype', None) == 'S':
        node.img_style["fgcolor"] = 'red'
        
    elif getattr(node, 'evoltype', None) == 'D':
        node.img_style["fgcolor"] = 'green'
    
        

ts = TreeStyle()
ts.layout_fn = [layout]
ts.show_leaf_name = False

name_tree = os.path.basename(sys.argv[1])
path_out = sys.argv[4]  

outfile =path_out+name_tree+'_'+TARGET_LCA+'_'+'doms'+'.pdf'
print(len(t))
t.render( outfile, w=350, units="mm", tree_style=ts)

seqs = SeqGroup(sys.argv[3])
seqs_dict = {}

for n, (name,seq,_) in enumerate(seqs):
    seqs_dict[name] = seq

seqs_in_og = []
for name_node,seqs in og_doms.items():
    if len(seqs) >3:
        out_name = name_node+'_'+name_tree+'_'+TARGET_LCA
        with open(path_out+out_name+'.faa', 'a') as f_out:
            for s in seqs:
                seqs_in_og.append(s)
                aa = seqs_dict[s]
                f_out.write('>'+s+'\n'+aa+'\n')

with open(path_out+'not_og-'+name_tree+'.faa', 'a') as f_out:
    for s in total_seqs:
        if s not in seqs_in_og:
            aa = seqs_dict[s]
            f_out.write('>'+s+'\n'+aa+'\n')



