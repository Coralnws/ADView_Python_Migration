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

    def add_tree_collection(self,treefile=None,type="newick"):
        self.tc = dendropy.TreeList.get_from_path(treefile, schema=type)


def print_tree(tree):
    print(tree.as_ascii_plot())

def get_node_list(tree):
    return tree.leaf_node_iter()

def get_leaf_node_amount(tree):
    return sum(1 for node in tree.leaf_node_iter())

'''
mytree = myTree(treefile = "Data/69species/astral.FAA.trim50genes.final.tre")
mytree.add_tree_collection(treefile = "Data/69species/tc.tre")
for tree in mytree.tc:
    print(tree)
'''

"""
def filter_node_selected(node, x, y):
    if hasattr(node, 'x1'):
        if node.x1 <= x <= node.x2 and node.y1 <= y <= node.y2:
            print(node.taxon)
            return True
        return False


def add_coordinate(node):
    if node.is_leaf() and node.taxon.label == "Uronema sp":
        node.x1 = 10
        node.x2 = 20
        node.y1 = 10
        node.y2 = 20

mytree = myTree(treefile = "../../Data/69species/astral.FAA.trim50genes.final.tre")
print(get_leaf_node_amount(mytree.rt))
mytree.rt.find_node(lambda node: add_coordinate(node))
node_selected = mytree.rt.find_node(lambda node: filter_node_selected(node, 13, 13))
print(node_selected)
"""