from ete3 import PhyloTree, NCBITaxa, TreeStyle, TextFace, NodeStyle
from ete3 import SeqGroup
import sys
from collections import Counter, OrderedDict
import numpy as np
import os
import string
import random
import json
import pandas as pd
from ete3.treeview import random_color
from ete3 import ProfileFace ,ImgFace, RectFace
from ete3.treeview.faces import add_face_to_node


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
            if float(parser[node.name][13]) < SO_FILTER or getattr(node, 'evoltype', None) == 'S':
                return True
            else:
                return False


def parse_leaves(l, sp_ch1, euk, euk_seq, bact, bact_seq, archa, archa_seq):
    sp = l.name.split('.')[0]
    sp_ch1.add(sp)
            
    lin = ncbi.get_lineage(sp)
    if 2759 in lin:
        euk.add(sp)
        euk_seq.append(l.name)
    elif 2 in lin:
        bact.add(sp)
        bact_seq.append(l.name)
    elif 2157 in lin:
        archa.add(sp)
        archa_seq.append(l.name)
    try:
        arc = data[l.name]
    except:
        n = re.sub('_', ':',l.name)
        arc = data[n]
    return arc


def parse_node(node):
    
    total_leaves = len(node)
    ch1 = node.children[0]
    ch2 = node.children[1]
    
    
    childs = node.get_children()
    child_names = []
    for c in childs:
        child_names.append(c.name)
    try:
        parent_node = node.up
        parent_name = parent_node.name
        
    except:
        parent_name = 'root'
        
    n_sister = node.get_sisters()
    if len(n_sister) == 0:
        sister_name = 'root'
    elif n_sister[0].is_leaf():
        sister_name = 'leaf'
    else:
        sister_name = n_sister[0].name
    
    leaves_under_ch1 = str(len(ch1))
    leaves_under_ch2 = str(len(ch2))

    euk = set()
    euk_seq = []
    bact = set()
    bact_seq = []
    archa = set()
    archa_seq = []
    
    
    sp_ch1 = set()
    arc_ch1 = []
    leaves_1= ch1.get_leaves()
    
    for l in leaves_1:
        arc_ = parse_leaves(l, sp_ch1, euk, euk_seq, bact, bact_seq, archa, archa_seq)
        arc_ch1.append(arc_)

    count_arc1 = dict(Counter(arc_ch1))

    for k,v in count_arc1.items():
        v = int(v)/int(leaves_under_ch1)
        v = round(v, 3)
        count_arc1[k] = str(v)
    
    sp_ch2 = set()
    arc_ch2 = []
    leaves_2= ch2.get_leaves()
    for l in leaves_2:
        arc_ = parse_leaves(l, sp_ch2, euk, euk_seq, bact, bact_seq, archa, archa_seq)
        arc_ch2.append(arc_)

    count_arc2 = dict(Counter(arc_ch2))
    for k,v in count_arc2.items():
        v = int(v)/int(leaves_under_ch2)
        v = round(v, 3)
        count_arc2[k] = str(v)

    #print(count_arc1)
    #print(count_arc2)

    dups_ch1 = 0
    for n in ch1.traverse("preorder"):
        if getattr(n, 'evoltype', None) == 'D':
            dups_ch1+=1
    dups_ch2 = 0
    for n in ch2.traverse("preorder"):
        if getattr(n, 'evoltype', None) == 'D':
            dups_ch2+=1

    total_dups = dups_ch1+dups_ch2

    t1 = ncbi.get_topology(list(sp_ch1))
    t2 = ncbi.get_topology(list(sp_ch2))
    t1_sci_name = t1.sci_name
    t2_sci_name = t2.sci_name

    shared_sp = set(sp_ch1&sp_ch2)
    len_shared_sp = str(len(shared_sp))
    total_sp = set()
    total_sp.update(sp_ch1)
    total_sp.update(sp_ch2)

    euk_len = 'euk_'+str(len(euk))+'_'+str(len(euk_seq))
    bact_len = 'bact_'+str(len(bact))+'_'+str(len(bact_seq))
    arc_len = 'arc_'+str(len(archa))+'_'+str(len(archa_seq))
    
    tall = ncbi.get_topology(list(total_sp))
    tall_sci_name = tall.sci_name
    
    shared_arc = list(set(arc_ch1).intersection(arc_ch2))
    len_shared_arc = str(len(shared_arc))
    total_arc = set()
    total_arc.update(set(arc_ch1))
    total_arc.update(set(arc_ch2))
    
    #nw_t = str(node.write(format=1))
    #json_out[node.name] = nw_t
    results = [node.name, parent_name, sister_name, child_names[0], child_names[1] , str(len(node)), leaves_under_ch1, leaves_under_ch2, str(dups_ch1), str(dups_ch2), str(len(sp_ch1)), str(len(sp_ch2)), len_shared_sp, sos_dict[node.name], str(len(set(arc_ch1))), str(len(set(arc_ch2))), len_shared_arc, tall_sci_name, t1_sci_name, t2_sci_name, euk_len, bact_len, arc_len, (count_arc1), (count_arc2), (total_arc)]

    return results 


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
SO_FILTER = float(sys.argv[3])
#result_tab = open(sys.argv[4], 'w')

with open('/home/plaza/projects/eggnog6/pfamA_families/eggnog6_pfamParse/domains_sorted.json') as json_file:
    data = json.load(json_file)


with open(sys.argv[4]) as emap:
    emapper_annot = json.load(emap)

 
for node in t.traverse("preorder"):
    if not node.is_leaf():
        node.name = id_generator()

sos_dict = {}
for ev in events:
    in_seq = ev.in_seqs
    out_seq = ev.out_seqs
    all_seq = in_seq
    all_seq.update(out_seq)
    ancestor = t.get_common_ancestor(all_seq)
    sos_dict[ancestor.name] = str(ev.sos)


parser = {}
for node in t.traverse("preorder"):
    if not node.is_leaf():
        r = parse_node(node)
        parser[node.name] = r


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


#RATIO_SEQ_SP = float(sys.argv[3])

# finds the nodes that can be collpased as OGs given the TARGET_LCA. 
# The decision is made by the is_leaf_of() function, which calculates several scores per node
scores = []
results_ = []
go_results = {}

for leaf in t.iter_leaves(is_leaf_fn=is_leaf_og):
    scores.append([leaf.score1, leaf.score2, len(CONTENT[leaf]), leaf.nspcs, leaf.lca, leaf])
    leaf.img_style['draw_descendants'] = False
    if leaf.is_leaf() == False:

        r = parser[leaf.name]

    go_terms = []
    for s in CONTENT[leaf]:
        try:
            go_ = emapper_annot[s.name]
            for g in go_:
                go_terms.append(g)
        except:
            print(s.name)
            go_terms.append('-')
        
    go_to_count = (go for go in go_terms if go[:1].isupper())
    c = Counter(go_to_count)
    go_freq = (c.most_common(5))
    go_results[leaf.name] = []
    for g in go_freq:
        go_results[leaf.name].append(g[0])

    seqs_no_go = c['-']
    

    leaf.add_face(TextFace("LCA:%s; NSPCS:%03s; NSEQS:%05d; INP_RATE:%0.2f; SCORE1:%03f; SCORE2:%03f; DOMS:%03d; SO:%s; total_GO:%s; seqs_no_GO:%s; GO_freq:%s" \
        %(leaf.lca, leaf.nspcs, leaf.nseqs, leaf.inparalogs_rate, leaf.score1, leaf.score2, len(r[25]), r[13], len(c), seqs_no_go, go_freq)), column=0, position = "branch-right")
        
    leaf.img_style['size'] = leaf.nspcs
    
scores.sort(reverse=True, key=lambda x: x[:-1])
covered_seqs = 0
for s in scores[:10]:    
    covered_seqs += s[2]
    print(s[:5], "%0.2f" %covered_seqs, "%0.2f" %(covered_seqs/len(taxa)), sep='\t')
    #print(s[-1].get_ascii(attributes=["named_lineage"]))
    

print(go_results)

go_order = set()
node2go = {}
for k, val in go_results.items():
    for v in val:
        go_order.add(v)

go_order = list(go_order)


for k, val in go_results.items():
    node2go[k] = []
    for g in go_order:
        if g in val:
            node2go[k].append('1')
        else:
            node2go[k].append('0')
print(node2go)


colors = random_color(num=len(go_order), l=0.5, s=0.5)
go2color = {go:colors[i] for i, go in enumerate(go_order)}

img_face = {}
# Create faces based on external images
for colidx, go in enumerate(go_order): 
    color = go2color[go]
    img_face[colidx] = RectFace(10, 10, fgcolor=color, bgcolor=color)
img_face['white'] = RectFace(10, 10, fgcolor="white", bgcolor="white")


# Draw resulting grouping
def layout(node):
    if not node.is_leaf():
        node.add_face(TextFace(node.lca+'_'+str(len(node))), column=0, position = "branch-right")
    
    
    if node.name in node2go:
        for colidx, go_present in enumerate(node2go[node.name]):
            if go_present=='1': 
                color = go2color[go_order[colidx]]
            else: 
                color = "white"
            f = RectFace(10, 10, fgcolor=color, bgcolor=color)
            add_face_to_node(f, node, column=colidx, position="aligned")

    if getattr(node, 'evoltype', None) == 'S':
        node.img_style["fgcolor"] = 'red'
        
    elif getattr(node, 'evoltype', None) == 'D':
        node.img_style["fgcolor"] = 'green'
        #node.img_style["draw_descendants"] == False

    
        

ts = TreeStyle()
ts.layout_fn = [layout]
ts.show_leaf_name = False

for colidx, name in enumerate(go_order): 
    nameF = TextFace(name)
    nameF.rotation = 90
    ts.aligned_header.add_face(nameF, column=colidx)


name = os.path.basename(sys.argv[1])  

outfile ='new2_'+name+'_'+TARGET_LCA+'_'+str(SO_FILTER)+'.pdf'
print(len(t))
t.render( outfile, w=350, units="mm", tree_style=ts)