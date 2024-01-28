import dendropy
from Utils import *
import myCanvas as myCanvas
import importlib
importlib.reload(myCanvas)

class myTree:
    rt = None  # Reference Tree
    tc = None  # Treelist
    rt_canvas = None
    rt_view_support = False

    def __init__(self,treefile = None,type="newick"):
        self.read_rt(treefile = treefile, type = type)

    def read_rt(self,treefile,type):
        self.rt = dendropy.Tree.get(path=treefile, schema=type)

    def reference_tree(self,view_support=False):
        height = get_leaf_node_amount(self.rt) * Y_INTERVAL + Y_INTERVAL

        if not self.rt_canvas or self.rt_canvas.view_support != view_support:
            self.rt_canvas = myCanvas.rtCanvas(self.rt,width = CANVAS_MAX_WIDTH,height = height,view_support = view_support)

        return self.rt_canvas

    def set_outgroup(self,outgroup_taxon):
        # ["Uronema sp", "Monomastix opisthostigma", "Pyramimonas parkeae", "Nephroselmis pyriformis"]
        mrca = self.rt.mrca(taxon_labels=outgroup_taxon)
        self.rt.reroot_at_edge(mrca.edge, update_bipartitions=False)

def print_tree(tree):
    print(tree.as_ascii_plot())

def get_node_list(tree):
    return tree.leaf_node_iter()

def get_leaf_node_amount(tree):
    return sum(1 for node in tree.leaf_node_iter())

# mytree = myTree(treefile = "../../Data/69species/astral.FAA.trim50genes.final.tre")
# print(get_leaf_node_amount(mytree.rt))