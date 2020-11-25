from ete3 import PhyloTree, NCBITaxa, TreeStyle
import sys
from collections import Counter, OrderedDict
import numpy as np

class OrderedCounter(Counter, OrderedDict):
     'Counter that remembers the order elements are first seen'
     def __repr__(self):
         return '%s(%r)' % (self.__class__.__name__,
                            OrderedDict(self))
     def __reduce__(self):
         return self.__class__, (OrderedDict(self),)


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


ncbi = NCBITaxa()
t = PhyloTree(sys.argv[1])

# If tree needs to be rooted
t.set_outgroup(t.get_midpoint_outgroup())

t.set_species_naming_function(lambda x: x.split('.')[0])
t.annotate_ncbi_taxa()

# Let's cache the list of leaves under each internal node. 
# This is a global variable.
CONTENT = t.get_cached_content()

# Which is the target level to create OGs? This is global variable used by several funcions
TARGET_LCA = sys.argv[2]

# convert the lineage of each leaf into a set, for faster look ups
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
for leaf in t.iter_leaves(is_leaf_fn=is_leaf_og):
    scores.append([leaf.score1, leaf.score2, len(CONTENT[leaf]), leaf.nspcs, leaf.lca, leaf])
    leaf.img_style['draw_descendants'] = False
    leaf.name = "LCA:%s NSPCS:%03s NSEQS:%05d INP_RATE:%0.2f" %(leaf.lca, leaf.nspcs, leaf.nseqs, leaf.inparalogs_rate)
    leaf.img_style['size'] = leaf.nspcs


# Print the scores inferred for each node
scores.sort(reverse=True, key=lambda x: x[:-1])
covered_seqs = 0
for s in scores[:10]:    
    covered_seqs += s[2]
    print(s[:5], "%0.2f" %covered_seqs, "%0.2f" %(covered_seqs/len(taxa)), sep='\t')
    #print(s[-1].get_ascii(attributes=["named_lineage"]))
    

# Draw resulting grouping
def layout(node):
    if node.img_style["draw_descendants"] == False:
        node.img_style["size"] = node.nspcs

ts = TreeStyle()
ts.layout_fn = [layout]
t.show(tree_style=ts)